#version 330 core

in vec3 inColor;
in vec2 vertexUV;
in vec3 inVertex;
in vec3 vertexNormal;

out vec3 vcolor;
out vec2 UV;
out vec3 fragNormal;
out vec3 fragVert;

uniform mat4 proy;
uniform mat4 rot;
uniform mat4 view;




void main()
{
	UV=vertexUV;
	vcolor = inColor;

	mat4 model = view*rot;

	fragNormal = inverse(transpose(mat3(model)))*vertexNormal;
	fragVert = (model * vec4(inVertex,1)).xyz;


    gl_Position =  proy*vec4(fragVert,1);
}