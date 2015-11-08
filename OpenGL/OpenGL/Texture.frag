#version 330 core

in vec2 UV;
in vec3 fragVert;

uniform sampler2D myTextureSampler;

out vec4 outColor;

void main(){
	//color plano
	//outColor = vec4(1,0,0,0); 

	//color normal
	outColor = vec4(texture(myTextureSampler, UV).rgb,0);
	
}