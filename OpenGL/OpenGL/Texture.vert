#version 330 core

in vec2 vertexUV;
in vec3 inVertex;

out vec2 UV;
out vec3 fragVert;

uniform mat4 proy;
uniform mat4 rot;
uniform mat4 view;




void main()
{
	UV=vertexUV;

	mat4 model = view*rot;

	fragVert = (model * vec4(inVertex,1)).xyz;


    gl_Position =  proy*model*vec4(inVertex,1);
}