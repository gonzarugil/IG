#version 420
#define MAX_ITERATIONS 2000


in vec3 entryPoint;

out vec4 outColor;

uniform sampler3D dataTex;
uniform sampler2D exitPointTex;
uniform ivec2    screenSize;

//TODO esto no es eficiente, pasar VMatrix y VMatrixInv y vmNMatrix por uniform
//matriz para el punto
uniform mat4 model;
uniform mat4 view;
uniform mat4 proy;


float isoValue = 0.4;
float stepSize = 0.25;
vec4 color = vec4(1,1,1,1);
vec4 backgroundColor = vec4(0, 0, 1, 0);
//float  isovalue = 0.5; 
vec3	  volumeDim = vec3(256,256,64); //Esto es texSize
float alphaStepFactor = 1.0f / stepSize;

vec3 getVoxelID(in vec3 v, in vec3 p){ //pasar incDir como parámetro?¿
	vec3 s = sign(v);
	vec3 pointAux = p + s * 0.01f;
	vec3 c = ceil(pointAux) - vec3(0.5f);
	vec3 f = floor(pointAux) + vec3(0.5f);
	vec3 pos = max(s, vec3(0));
	vec3 neg = -min(s, vec3(0));
	return pos * f + neg * c;
}

void computeDensity(in vec3 voxelId, out float d0, out float d1, out float d2,
	out float d3, out float d4, out float d5, out float d6, out float d7)
{
	
	d0 = texelFetch(dataTex, ivec3(voxelId + vec3(-0.5f)), 0).r;
	d1 = texelFetch(dataTex, ivec3(voxelId + vec3(-0.5f, -0.5f, 0.5f)), 0).r;
	d2 = texelFetch(dataTex, ivec3(voxelId + vec3(0.5f, -0.5f, 0.5f)), 0).r;
	d3 = texelFetch(dataTex, ivec3(voxelId + vec3(0.5f, -0.5f, -0.5f)), 0).r;
	d4 = texelFetch(dataTex, ivec3(voxelId + vec3(-0.5f, 0.5f, -0.5f)), 0).r;
	d5 = texelFetch(dataTex, ivec3(voxelId + vec3(-0.5f, 0.5f, 0.5f)), 0).r;
	d6 = texelFetch(dataTex, ivec3(voxelId + vec3(0.5f)), 0).r;
	d7 = texelFetch(dataTex, ivec3(voxelId + vec3(0.5f, 0.5f, -0.5f)), 0).r;
}

void main()
{
	mat4 vmMatrix = view*model;
	mat4 vmMatrixInverse = inverse(vmMatrix);
	mat4 vmNMatrix = transpose(inverse(vmMatrix));

	vec2 screenCoords = gl_FragCoord.st/screenSize;
    vec3 exitPoint = texture(exitPointTex, screenCoords).xyz;
	    
	if (entryPoint.z <= exitPoint.z)
	{
		outColor = backgroundColor;
	}
	else{
		//Direccion del rayo. Calculamos el rayo un cubo en el que el 
		vec3 volumeDimInv = 1.0f / volumeDim;
		vec3 dir = exitPoint - entryPoint;
		//vector escalado en coordenadas de la textura
		dir = volumeDim*(vmMatrixInverse * vec4(dir, 0)).xyz;
		vec3 pto = volumeDim*(vmMatrixInverse * vec4(entryPoint, 1)).xyz;

		//Direcciones de la recta

		vec3 incDir = sign(dir);
		
		//variables sobre las que iteramos (t de la recta)
		float t = 0; //la t del entrypoint es 0 (la recta se calcula usando este punto)
		float tfin = 1; //la ultima t será 1 ya que el vector está calculado como exitpoint-entrypoint
		int i = 0;
		vec3 voxelID;
		vec3 currentPoint = pto; //primer punto de corte
		vec4 colorAcum = vec4(0);
		while (t < tfin && i < MAX_ITERATIONS){
			//obtenemos el voxelID del punto
			voxelID = getVoxelID(dir, currentPoint);
			//AQUI SE HACEN LOS CALCULOS DE DENSIDADES Y ETC
			//colorAcum += vec4(vec3(0.003), 0);
			float d0, d1, d2, d3, d4, d5, d6, d7;
			computeDensity(voxelID, d0, d1, d2, d3, d4, d5, d6, d7);
			//colorAcum += 0.1 * vec4(vec3(d6), 0);
			//calculamos siguiente entrypoint y su t correspondiente

			//primero calculamos que planos se cortan
			vec3 faces = incDir * 0.5f + voxelID;

			/*//t del plano x
			float tx = (faces.x - pto.x) / dir.x;
			//t del plano y
			float ty = (faces.y - pto.y) / dir.y;
			//t del plano z
			float tx = (faces.z - pto.z) / dir.z;*/
			vec3 tcortes = (faces - pto) / dir;
			//obtenemos la menor t
			float tnext = min(tcortes.x,min( tcortes.y, tcortes.z));
			//calculamos el entrypoint correspondiente
			vec3 pointnext = pto + dir * tnext;

			//actualizamos las variables del while
			t = tnext;
			currentPoint = pointnext;
		}
		outColor = colorAcum;


		//float len = length(dir);
		//vec3 deltaDir = (volumeDimInv)*(stepSize / len)*dir;



		//Puntos de la texturas
		//vec3 voxelCoord = (vmMatrixInverse * vec4(entryPoint, 1)).xyz;

	}
	/*
	//Color y longitud acumulados en el bucle    
	vec4 colorAcum = vec4(0.0); 
    float lengthAcum = 0.0;
	float lastDensity = 0.0;

	for(int i = 0; i < MAX_ITERATIONS; i++)
    {
    	float intensity =  texture(dataTex, voxelCoord).x;
		vec3 texelCoord = voxelCoord*volumeDim;


		//Esto es para asegurarnos de que no cogemos la "capa exterior" del cubo (Que daria error de acceso a la textura)
		if (!(any(lessThan (texelCoord, vec3(0.5f))) || any(greaterThan (texelCoord, volumeDim-vec3(0.5f)))))
		{
			
			if ((isoValue - lastDensity)*(isoValue - intensity) < 0.0)
			{
				colorAcum = color;
				break;
			}
		
		}


    	voxelCoord += deltaDir;
    	lengthAcum += stepSize;
    	
		if (lengthAcum >= len )
    	{	
    		break;  	
    	}	
    	

		lastDensity = intensity;
    }//For
	
	outColor = colorAcum;
    }*/
}
