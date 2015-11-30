
#version 420

uniform mat4 proy;
uniform mat4 view;
uniform mat4 model;

in vec3 inVertex;  
	
out vec3 entryPoint;



void main()
{
	mat4 vmMatrix = view*model;
	mat4 pvmMatrix = proy * vmMatrix;
    entryPoint = (vmMatrix * vec4(inVertex,1.0f)).xyz;
	gl_Position = pvmMatrix*vec4(inVertex,1.0f); 
}