#pragma once

class CLight
{
public:
	CLight(void);
	~CLight(void);

	

	void setColour(const float *ambient, const float *diffuse, const float *specular);
	void getColour(float *ambient, float *diffuse, float *specular);

	void setColour(const float *ambient, const float *diffuse, const float *specular, const float *att);
	void getColour(float *ambient, float *diffuse, float *specular, float *att);

	void setPosition(const float *pos, const float angle = 180.0f, const float *direction = 0L);
	void getPosition(float *pos, float &angle, float *direction);


protected:
	//El alumno podrá definir nuevos parametros de configuración.
	float ambient[4];
	float diffuse[4];
	float specular[4];
	float att[3];     //Atenuación(cte,lineal,cuadrática)

	float direction[4];
	float pos[4];
	float angle; //180 es una luz puntual. 
};


