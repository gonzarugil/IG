#version 330 core

in vec3 vcolor;
in vec2 UV;
in vec3 fragNormal;
in vec3 fragVert;

uniform vec3 lightposition;
uniform vec3 lightintensity;
uniform float lightambientCoefficient;
uniform float lightattenuation;
uniform float materialShininess;
uniform vec3 materialSpecularColor;

out vec4 outColor;

uniform bool transparency;
uniform sampler3D densityTex;

uniform mat4 rot;
uniform mat4 view;

void main()
{
	
	vec4 color = vec4(0,0,0,0);
	float z = 0;

	//Raycasting Ortogonal
	
	/*
	while(z<=1){
		color +=  texture(densityTex,vec3(UV,z));
		z=z+0.01;
	}*/
	
	
	//Raycasting en perspectiva desde la camara  (No se si esta bien)
	/*mat4 invModel = inverse (view*rot);
	vec3 pto = ((invModel*vec4(fragVert,1)).xyz+vec3(2.0,2.0,-2.0))*0.25;
	vec3 v = (invModel*vec4(fragVert,0)).xyz/4.0;
	*/
	mat4 model = view*rot;
	mat4 invModel = inverse(model);
	vec3 pto = (invModel*vec4(fragVert,1)).xyz;
	vec3 v = (invModel*vec4(fragVert,0)).xyz;

	float tfin = (-1-pto.z)/v.z;
	float t = 0;
	while (t<tfin)
	{
		vec3 uvnew = pto + v *t;
		uvnew.z *= -1.0;
		color = color + texture(densityTex,uvnew);
		t = t + 0.01;
	}
	
	outColor = color.rgbr;
	//outColor = vec4 (pto.zzz+vec3(1.5),0);
    //outColor = vec4(diffuse * lightintensity *surfaceColor.rgb , surfaceColor.a);
	//outColor = vec4(brightness*lightintensity, surfaceColor.a);
}
