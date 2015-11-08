#include "FBO.h"
#include <iostream>
#include "Dependencies\glew\glew.h"
#include "Dependencies\freeglut\freeglut.h"






CFBO::CFBO(void)
{
	//Bloque 1
	glGenTextures(1, &colorTex);
	glBindTexture(GL_TEXTURE_2D, colorTex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	//NULL means reserve texture memory, but texels are undefined
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 256, 256, 0, GL_BGRA, GL_UNSIGNED_BYTE, NULL);
	//-------------------------
	//Bloque 2
	resize(250, 250); //Se explicará más adelante

	//Bloque 3
	glGenFramebuffers(1, &id);
	glBindFramebuffer(GL_FRAMEBUFFER, id);

	//Bloque 4
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
		GL_TEXTURE_2D, colorTex, 0);
	
	//Bloque 5
	const GLenum buffs[1] =
	{GL_COLOR_ATTACHMENT0};
	glDrawBuffers(1, buffs);

	//Bloque 6
	if (GL_FRAMEBUFFER_COMPLETE !=
		glCheckFramebufferStatus(GL_FRAMEBUFFER))
	{
		std::cerr << "Error configurando el FBO" << std::endl;
		exit(-1);
	}

	//Bloque 7
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

CFBO::~CFBO(void)
{
	glDeleteFramebuffers(1, &id);
	glDeleteTextures(1, &colorTex);
}


void CFBO::resize(unsigned int w, unsigned int h)
{
	//Bloque 1
	glBindTexture(GL_TEXTURE_2D, colorTex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, w, h, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
		GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
		GL_NEAREST);

	glBindTexture(GL_TEXTURE_2D, 0);

}


void CFBO::activate()   {
	glBindFramebuffer(GL_FRAMEBUFFER, id);
}
void CFBO::deactivate() {
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

unsigned int CFBO::getColorTexId()		{ return colorTex; }

