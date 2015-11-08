#version 420
#define MAX_ITERATIONS 4000

in vec3 entryPoint;

uniform sampler3D dataTex;
uniform sampler2D exitPointTex;

uniform vec2      screenSize;
uniform float     stepSize;

//matriz para el punto
uniform mat4 rot;
uniform mat4 view;
uniform mat4 proy;

out vec4 outColor;



void main()
{
	mat4 vmMatrix = view*rot;
	mat4 vmMatrixInverse = inverse(vmMatrix);

	vec4 color = vec4(0,0,0,0);
	vec2 screenCoords = gl_FragCoord.st/screenSize;
    vec3 exitPoint = texture(exitPointTex, screenCoords).xyz;
    
	if (entryPoint.z == exitPoint.z)
	{
		discard;
	}else{
		vec3 v = exitPoint - entryPoint;
		v = (vmMatrixInverse * vec4 (v,0)).xyz;
        vec3 pto   = (vmMatrixInverse * vec4 (entryPoint,1)).xyz;
        vec3 ptofin   = (vmMatrixInverse * vec4 (exitPoint,1)).xyz;
		//pto *= texSize;
		//v *= texSize;


        float tfin = (ptofin.z-pto.z)/v.z;
        float t = 0;
		//aqui va el while
		while (t<tfin){
			vec3 uvnew = pto + v *t;
			uvnew.z *= -1.0;
			color = color + texture(dataTex,uvnew);
			t = t + 0.01;
		}
	}
	outColor = color;
	//outColor = vec4(1.0,0.0,0.0,0.0);
}
