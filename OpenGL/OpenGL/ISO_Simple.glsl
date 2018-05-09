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


uniform float isoValue;
float stepSize = 0.25;
vec4 color = vec4(1,1,1,1);
vec4 backgroundColor = vec4(0, 0, 1, 0);
vec3	  volumeDim = vec3(256,256,64);
float alphaStepFactor = 1.0f / stepSize; 

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
	}else{
	//Direccion del rayo. Calculamos el rayo un cubo en el que el 
	vec3 volumeDimInv = 1.0f / volumeDim;
	vec3 dir = exitPoint - entryPoint;
	dir = volumeDim*(vmMatrixInverse * vec4 (dir,0)).xyz;
	float len = length(dir);
	vec3 deltaDir = (volumeDimInv)*(stepSize/len)*dir;

	//Puntos de la texturas
	vec3 voxelCoord = (vmMatrixInverse * vec4 (entryPoint,1)).xyz;

	//Color y longitud acumulados en el bucle    
	vec4 colorAcum = vec4(0.0); 
    float lengthAcum = 0.0;
	float lastDensity = 0.0;

	for(int i = 0; i < MAX_ITERATIONS; i++)
    {
    	float intensity =  texture(dataTex, voxelCoord).x;
		vec3 texelCoord = voxelCoord*volumeDim;

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
    }
}
