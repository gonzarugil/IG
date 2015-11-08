#version 330
// for raycasting

in vec3 inVertex;

uniform mat4 view;
uniform mat4 rot;
uniform mat4 proy;

out vec3 vertex;


void main()
{
	mat4 vmMatrix = view*rot;
	mat4 pvmMatrix = proy * vmMatrix;
    vertex = (vmMatrix*vec4(inVertex, 1.0)).xyz;
    gl_Position = pvmMatrix * vec4(inVertex, 1.0);
}