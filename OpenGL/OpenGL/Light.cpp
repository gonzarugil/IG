#include "Light.h"




CLight::CLight(void)
{
	ambient[0] = 0.2f;
	ambient[1] = 0.2f;
	ambient[2] = 0.2f;
	ambient[3] = 1.0f;

	diffuse[0] = 0.5f;
	diffuse[1] = 0.5f;
	diffuse[2] = 0.5f;
	diffuse[3] = 1.0f;

	specular[0] = 1.0f;
	specular[1] = 1.0f;
	specular[2] = 1.0f;
	specular[3] = 1.0f;

	att[0] = 1.0f;
	att[1] = 0.0f;
	att[2] = 0.0f;

	pos[0] = 0.0f;
	pos[1] = 0.0f;
	pos[2] = 0.0f;
	pos[3] = 1.0f;

	angle = 180.0f;

	direction[0] = 0.0f;
	direction[1] = 0.0f;
	direction[2] = -1.0f;
	direction[3] = 0.0f;
}

CLight::~CLight(void)
{
}




void CLight::setColour(const float *ambient, const float *diffuse, const float *specular)
{
	this->ambient[0] = ambient[0];
	this->ambient[1] = ambient[1];
	this->ambient[2] = ambient[2];
	this->ambient[3] = ambient[3];

	this->diffuse[0] = diffuse[0];
	this->diffuse[1] = diffuse[1];
	this->diffuse[2] = diffuse[2];
	this->diffuse[3] = diffuse[3];

	this->specular[0] = specular[0];
	this->specular[1] = specular[1];
	this->specular[2] = specular[2];
	this->specular[3] = specular[3];
}

void CLight::setColour(const float *ambient, const float *diffuse, const float *specular, const float *att)
{
	this->ambient[0] = ambient[0];
	this->ambient[1] = ambient[1];
	this->ambient[2] = ambient[2];
	this->ambient[3] = ambient[3];

	this->diffuse[0] = diffuse[0];
	this->diffuse[1] = diffuse[1];
	this->diffuse[2] = diffuse[2];
	this->diffuse[3] = diffuse[3];

	this->specular[0] = specular[0];
	this->specular[1] = specular[1];
	this->specular[2] = specular[2];
	this->specular[3] = specular[3];

	this->att[0] = att[0];
	this->att[1] = att[1];
	this->att[2] = att[2];
}

void CLight::getColour(float *ambient, float *diffuse, float *specular)
{
	ambient[0] = this->ambient[0];
	ambient[1] = this->ambient[1];
	ambient[2] = this->ambient[2];
	ambient[3] = this->ambient[3];

	diffuse[0] = this->diffuse[0];
	diffuse[1] = this->diffuse[1];
	diffuse[2] = this->diffuse[2];
	diffuse[3] = this->diffuse[3];

	specular[0] = this->specular[0];
	specular[1] = this->specular[1];
	specular[2] = this->specular[2];
	specular[3] = this->specular[3];
}

void CLight::getColour(float *ambient, float *diffuse, float *specular, float *att)
{
	ambient[0] = this->ambient[0];
	ambient[1] = this->ambient[1];
	ambient[2] = this->ambient[2];
	ambient[3] = this->ambient[3];

	diffuse[0] = this->diffuse[0];
	diffuse[1] = this->diffuse[1];
	diffuse[2] = this->diffuse[2];
	diffuse[3] = this->diffuse[3];

	specular[0] = this->specular[0];
	specular[1] = this->specular[1];
	specular[2] = this->specular[2];
	specular[3] = this->specular[3];

	att[0] = this->att[0];
	att[1] = this->att[1];
	att[2] = this->att[2];

}

void CLight::setPosition(const float *pos, const float angle, const float *direction)
{
	this->pos[0] = pos[0];
	this->pos[1] = pos[1];
	this->pos[2] = pos[2];
	this->pos[3] = pos[3];

	if (angle >= 0.0f && angle <= 90.0f)
		this->angle = angle;

	if (direction)
	{
		this->direction[0] = direction[0];
		this->direction[1] = direction[1];
		this->direction[2] = direction[2];
		this->direction[3] = direction[3];
	}
}

void CLight::getPosition(float *pos, float &angle, float *direction)
{
	pos[0] = this->pos[0];
	pos[1] = this->pos[1];
	pos[2] = this->pos[2];
	pos[3] = this->pos[3];

	angle = this->angle;

	direction[0] = this->direction[0];
	direction[1] = this->direction[1];
	direction[2] = this->direction[2];
	direction[3] = this->direction[3];
}
