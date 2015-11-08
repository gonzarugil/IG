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
uniform sampler2D myTextureSampler;

void main()
{

	vec3 cameraPosition;
	cameraPosition.x = 0;
	cameraPosition.y = 0;
	cameraPosition.z = 0;


	//calculate the vector from this pixels surface to the light source
	vec3 surfacePos = fragVert;
    vec3 surfaceToLight = normalize(lightposition - surfacePos);
	//Color con textura
	vec4 surfaceColor = vec4(texture(myTextureSampler, UV).rgb,transparency);
	//Color solo luces
	//vec4 surfaceColor = vec4(1.0,1.0,1.0,transparency);
	//calculate the cosine of the angle of incidence
    float diffuseCoefficient = max(0.0, dot(fragNormal, surfaceToLight));
	//vec3 diffuse = diffuseCoefficient * surfaceColor.rgb * lightintensity;
	vec3 diffuse = diffuseCoefficient * lightintensity;


	//ambient component
	vec3 ambient = lightambientCoefficient * surfaceColor.rgb * lightintensity;

	//specular component 
	//vec3 incidenceVector = -surfaceToLight; //a unit vector
	//vec3 reflectionVector = reflect(incidenceVector, fragNormal); //also a unit vector
	//vec3 surfaceToCamera = normalize(cameraPosition - surfacePos); //also a unit vector
	//float cosAngle = max(0.0, dot(surfaceToCamera, reflectionVector));
	//float specularCoefficient = pow(cosAngle, materialShininess);
	//vec3 specularComponent = specularCoefficient * materialSpecularColor * lightintensity;

	//condensed specular
	vec3 surfaceToCamera = normalize(cameraPosition - surfacePos); //also a unit vector
	float specularCoefficient = 0.0;
	if(diffuseCoefficient > 0.0)
		specularCoefficient = pow(max(0.0, dot(surfaceToCamera, reflect(-surfaceToLight, fragNormal))), materialShininess);
	vec3 specular = specularCoefficient * materialSpecularColor * lightintensity;

	//attenuation
	float distanceToLight = length(lightposition - surfacePos);
	float attenuation = 1.0 / (1.0 + lightattenuation * pow(distanceToLight, 2));

	//linear color (color before gamma correction)
    vec3 linearColor = ambient + attenuation*(diffuse + specular);
	
	//final color (after gamma correction)
    vec3 gamma = vec3(1.0/2.2);
    outColor = vec4(pow(linearColor, gamma), surfaceColor.a);
	//outColor = vec4(1.0,0,0,0);
    //outColor = vec4(diffuse * lightintensity *surfaceColor.rgb , surfaceColor.a);
	//outColor = vec4(brightness*lightintensity, surfaceColor.a);
}

