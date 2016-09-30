#version 420
#define MAX_ITERATIONS 3000


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


uniform float isoValue;

//propiedades de la luz

vec3 Ia = vec3(0.6);
vec3 Id = vec3(0.5);
vec3 Is = vec3(1.0);
vec3 lpos = vec3(0,0,0);


//propiedades del material (HARDCODEADOS)
vec3 Ka = vec3(1, 0.5, 0.0);
vec3 Kd = vec3(1, 1, 1);
vec3 Ks = vec3(1.0);
float alpha = 500.0;


//float isoValue = 0.35;
float stepSize = 0.25;
vec4 color = vec4(1,1,1,1);
vec4 backgroundColor = vec4(0, 0, 1, 0);
//float  isovalue = 0.5; 
uniform ivec3 texSize;
vec3 volumeDim = texSize; //vec3(256,256,64); //Esto es texSize

float alphaStepFactor = 1.0f / stepSize;

vec3 getVoxelID(in vec3 v, in vec3 p){ //pasar incDir como parámetro?¿
	vec3 s = sign(v);
	vec3 pointAux = p + s * 0.01f;
	vec3 c = floor(pointAux) + vec3(0.5f);
	/*vec3 f = floor(pointAux) + vec3(0.5f);
	vec3 pos = max(s, vec3(0));
	vec3 neg = -min(s, vec3(0));
	return pos * f + neg * c;*/
	return c;

}

//---------------Solver ecuaciones cuadráticas---------------

int solveQuad(in vec3 coff, out vec3 res){
	res.x = (-coff.y + sqrt(pow(coff.y, 2) - 4 * coff.x*coff.z)) / 2 * coff.x;
	res.y = (-coff.y - sqrt(pow(coff.y, 2) - 4 * coff.x*coff.z)) / 2 * coff.x;
	return 2;
}
//---------------Solver ecuaciones cúbicas-------------------

float cbrt(float x){
	float EPSILON = 0.00000001f;
	if (abs(x) < EPSILON) return 0.0;
	if (x> 0.0) return pow(x, 0.33333333333333333333f);
	return -pow(-x, 0.33333333333333333333f);
}

float computeX(float As, float Ch, float Dh, float delta){
	float T0 = sign(Dh) * abs(As) * sqrt(-delta);
	float T1 = -Dh + T0;
	float p = cbrt(T1 * 0.5f);
	// TODO: Poner epsilon?
	float q = (T1 == T0) ? -p : -Ch / (p);
	return (Ch <= 0) ? p + q : -Dh / (p*p + q*q + Ch);
}

int solveCubic(in vec4 coff, out vec3 res){
	res = vec3(0.0);
	int nsol = 0;
	coff.yz *= 0.33333333333333;
	vec3 sigma = vec3(coff.x*coff.z - coff.y*coff.y,
		coff.x*coff.w - coff.y*coff.z,
		coff.y*coff.w - coff.z*coff.z);
	float delta = 4 * sigma.x *sigma.z - sigma.y*sigma.y;
	if (delta <= 0.0f){
		if (coff.y*coff.y*coff.y*coff.w > coff.x * coff.z*coff.z*coff.z){
			//Algoritmo A
			float x = computeX(coff.x, sigma.x, -2.0f*coff.y*sigma.x + coff.x*sigma.y, delta);
			//TODO: Controlar el caso A=0 (Exacto)
			res.x = (x - coff.y) / coff.x;
			nsol += 1;
		}
		else{
			//ALgoritmo D
			float x = computeX(coff.w, sigma.z, -coff.w * sigma.y + 2.0f * coff.z * sigma.z, delta);
			res.x = (-coff.w / (x + coff.z));
			nsol += 1;
		}
	}
	else{
		float Ch = sigma.x;
		float Dh = -2.0f * coff.y * sigma.x + coff.x * sigma.y;
		float theta = 0.3333333333333333f * abs(atan(-Dh, coff.x*sqrt(delta)));
		float sqCh2 = 2.0f * sqrt(-Ch);
		float cosTheta = cos(theta);
		float sinTheta = sin(theta);
		float x1 = sqCh2 * cosTheta;
		float x3 = sqCh2 * (-0.5f * cosTheta - 0.86602540378f * sinTheta);
		float xL = (x1 + x3 > 2.0f * coff.y) ? x1 : x3;
		if (coff.x == 0){
			res.x = (-coff.z + sqrt(pow(coff.z, 2) - 4 * coff.y*coff.w)) / 2 * coff.y;
			nsol += 1;
			res.y = (-coff.z - sqrt(pow(coff.z, 2) - 4 * coff.y*coff.w)) / 2 * coff.y;
			nsol += 1;
		}
		else{
			res.x = (xL - coff.y) / coff.x;
			nsol += 1;
			//TODO: Controlar el caso de A = 0 (exacto) Poner Infinito
			//Sol2 xs
			Ch = sigma.z;
			Dh = -coff.w * sigma.y + 2.0f * coff.z * sigma.z;
			theta = 0.33333333333333333f * abs(atan(-Dh, coff.w*sqrt(delta)));
			sqCh2 = 2.0f *sqrt(-Ch);
			cosTheta = cos(theta);
			sinTheta = sin(theta);
			x1 = sqCh2 * cosTheta;
			x3 = sqCh2 * (-0.5f * cosTheta - 0.86602540378f * sinTheta);
			float xS = (x1 + x3 < 2.0f * coff.z) ? x1 : x3;
			res.y = -coff.w / (xS + coff.z);
			nsol += 1;
			//The only time this will go wrong is when F^2=EG, which is when the quadratic has a double root. In
			//this case — when the largest and smallest roots of the cubic are equal — it means that we must have a
			//triple-root cubic. But in fact, this will only really cause problems if the triple root is at zero or infinity. We
			//will deal with this another time.
			//https://courses.cs.washington.edu/courses/cse590b/13au/lecture_notes/solvecubic_p5.pdf
			float E = (xS + coff.z)*coff.x;
			float F = -(xL - coff.y)*(xS + coff.z) - coff.x * (-coff.w);
			float G = -coff.w * (xL - coff.y);
			res.z = (xL - coff.y) / coff.x;
			nsol += 1;
		}
	}
	return nsol;
}


//-----------------------------------------------------------

vec4 shade(in vec3 norm, in vec3 colorpos){
	
	vec3 c = vec3(0.0);
	c = Ia * Ka;

	vec3 L = normalize(lpos - colorpos);
	vec3 diffuse = Id * Kd * dot(L, norm);
	c += clamp(diffuse, 0.0, 1.0);

	vec3 V = normalize(-colorpos);
	vec3 R = normalize(reflect(-L, norm));
	float factor = max(dot(R, V), 0.001);
	vec3 specular = Is * Ks * pow(factor, alpha);
	c += clamp(specular, 0.0, 1.0);

	return vec4(c, 0.0f);
	//return vec4(norm, 1);
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

bool checkIsoValue(in float d0, in float d1, in float d2, in float d3, in float d4,
	in float d5, in float d6, in float d7, in float isovalue)
{

	float v1sig = d0 - isoValue;
	if ((v1sig*(d1 - isovalue)) < 0.0) return true;
	if ((v1sig*(d2 - isovalue)) < 0.0) return true;
	if ((v1sig*(d3 - isovalue)) < 0.0) return true;
	if ((v1sig*(d4 - isovalue)) < 0.0) return true;
	if ((v1sig*(d5 - isovalue)) < 0.0) return true;
	if ((v1sig*(d6 - isovalue)) < 0.0) return true;
	if ((v1sig*(d7 - isovalue)) < 0.0) return true;
	return false;

}
/*
// Returns true when isoValue is inside the voxel
bool CheckIsoValue(in float d000, in float d001, in float d101, in float d100, in float d010,
	in float d011, in float d111, in float d110, in float isoValue)
{
	// Comparar todos los valores del voxel con el isoValue -> signos -1, 0, 1
	const float signs[] = { sign(d000 - isoValue), sign(d001 - isoValue), sign(d101 - isoValue), sign(d100 - isoValue),
		sign(d010 - isoValue), sign(d011 - isoValue), sign(d111 - isoValue), sign(d110 - isoValue), };

	// Identificar signos positivos y negativos (valores mayores y menores que isoValue)
	bool foundPositive = false;
	bool foundNegative = false;

	for (int i = 0; i < 8; i++)
	{
		if (signs[i] < 0.0)
			foundNegative = true;
		else
			foundPositive = true;
	}

	// Si hay positivo y negativo, hay una isosuperficie en el voxel
	bool changeInVoxel = foundPositive && foundNegative;

	return changeInVoxel
}
*/

bool checkSolution(in vec3 solutions, in float t,in float tnext, in int nsol, out float s){
	bvec3 validsolution = bvec3(false);
	validsolution.x = solutions.x > t && solutions.x <= tnext;
	if (nsol > 1) validsolution.y = solutions.y > t && solutions.y <= tnext;
	if (nsol > 2) validsolution.z = solutions.z > t && solutions.z <= tnext;
	vec3 valids = solutions * vec3(validsolution);
	vec3 invalid = tnext * vec3(not(validsolution));
	vec3 sum = valids + invalid;
	float temp = min(min(sum.x, sum.y), sum.z);
	s = temp;
	return s != tnext;
	

	
	
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
		volumeDim -= vec3(1.0f);
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
		while (t <= tfin && i < MAX_ITERATIONS){
			i++;
			//obtenemos el voxelID del punto
			voxelID = getVoxelID(dir, currentPoint);
			/*
			colorAcum = vec4((voxelID - vec3(0.5f))*volumeDimInv, 1f);
			break;*/

			//calculamos siguiente entrypoint y su t correspondiente

			//primero calculamos que planos se cortan
			vec3 faces = incDir * 0.5f + voxelID;

			vec3 tcortes = (faces - pto) / dir;
			//obtenemos la menor t
			float tnext = min(tcortes.x, min(tcortes.y, tcortes.z));
			//calculamos el entrypoint correspondiente
			vec3 pointnext = pto + dir * tnext;


			//AQUI SE HACEN LOS CALCULOS DE DENSIDADES Y ETC
			float d0, d1, d2, d3, d4, d5, d6, d7;
			computeDensity(voxelID, d0, d1, d2, d3, d4, d5, d6, d7);
			if (checkIsoValue(d0, d1, d2, d3, d4, d5, d6, d7, isoValue)){


				float A = d0;
				float X = d3 - d0;
				float Y = d4 - d0;
				float Z = d1 - d0;
				float XY = d7 + d0 - d3 - d4;
				float ZX = d2 + d0 - d4 - d1;
				float YZ = d5 + d0 - d4 - d3;
				float XYZ = d6 + d3 + d4 + d1 - d0 - d5 - d2 - d7;

				vec3 ptotrans = currentPoint - (voxelID - vec3(0.5f)); //translacion del voxel al origen
				vec4 ec;
				ec.x = XYZ*dir.x*dir.y*dir.z;
				ec.y = XYZ*(ptotrans.z*dir.y*dir.x + dir.z*ptotrans.y*dir.x + dir.z*dir.y*ptotrans.x) + XY*dir.y*dir.x + ZX*dir.z*dir.x + YZ*dir.z*dir.y;
				ec.z = XY*(ptotrans.y*dir.x + dir.y*ptotrans.x) + YZ*(ptotrans.z*dir.y + dir.z*ptotrans.y) + XYZ*(dir.z*ptotrans.y*ptotrans.x + ptotrans.z*ptotrans.y*dir.x + ptotrans.z*dir.y*ptotrans.x) + ZX*(ptotrans.z*dir.x + dir.z*ptotrans.x) + Y*dir.y + Z*dir.z + X*dir.x;
				ec.w = A - isoValue + XY*ptotrans.y*ptotrans.x + ZX*ptotrans.z*ptotrans.x + YZ*ptotrans.z*ptotrans.y + XYZ*ptotrans.z*ptotrans.y*ptotrans.x + X*ptotrans.x + Z*ptotrans.z + Y*ptotrans.y;

				//colorAcum += 0.1*vec4(ec.w);
				//Ecuación cúbica en funcion de T At^3 + Bt^2 + Ct + D -isoValue = 0 
				//Para encontrar los puntos de interseccion con el isovalor
				vec3 sol;
				int nsol;
				//Hay que controlar el caso en el que A=0 (Ecuacion de segundo grado)
				if (ec.x == 0){
					nsol = solveQuad(vec3(ec.y, ec.z, ec.w), sol);
				}
				else {
					nsol = solveCubic(ec, sol);
				}
				

				//Esto saca los t de corte con la isosuperficie
				//Hay que comprobar cuales caen dentro del voxel y dentro de estos, cual es el menor
				// tmin = t, tmax = tnext

				//comprobamos las soluciones obtenidas
				float s;
				bool validsolution;
				validsolution = checkSolution(sol,t,tnext,nsol,s);
			

				if (validsolution)
				{
					colorAcum = 2*vec4(isoValue,0,0,0);
					//calculo del punto
					vec3 ptosol = ptotrans + dir*s;
					//deshacemos la translacion de voxelID
					ptosol = ptosol + voxelID - vec3(0.5f);
					ptosol = ptosol / volumeDim;
					ptosol = (vmMatrix * vec4(ptosol, 1)).xyz;

					//Calculo de normales

					vec3 norm = vec3(0.0f);
					norm.x = X + XY * ptotrans.y + ZX * ptotrans.z + XYZ * ptotrans.y * ptotrans.z;
					norm.y = Y + XY * ptotrans.x + YZ * ptotrans.z + XYZ * ptotrans.x * ptotrans.z;
					norm.z = Z + ZX * ptotrans.x + YZ * ptotrans.y + XYZ * ptotrans.x * ptotrans.y;
					norm = norm * -1;

					norm = norm / volumeDim;// solo escalado porque es un vector y no es homogeneo
					norm = (vmNMatrix * vec4(norm, 0)).xyz; //coordenadas de la camara
					norm = normalize(norm);

					colorAcum = vec4(voxelID/volumeDim, 1);//shade(norm, ptotrans);
					break;
					
				}
				
			


			}


			

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
