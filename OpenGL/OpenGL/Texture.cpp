#include "Dependencies\glew\glew.h"
#include "Dependencies\freeglut\freeglut.h"
#include "Dependencies\freeimage\FreeImage.h"
#include <iostream>
#include <fstream>
#include <windows.h>
#include "Texture.h"

void CTexture::buildTexture()
{
	if (size[0] != 0 && size[1] != 0 && !_isInit)
	{
		_isInit = true;
		glGenTextures(1, &id);
		glBindTexture(GL_TEXTURE_2D, id);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, size[0], size[1], 0,
			GL_RGBA, GL_UNSIGNED_BYTE, (GLvoid*)map);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
			GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
			GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
			GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,
			GL_REPEAT);
	}
}



void CTexture::activate()
{
	if (_isInit)
		glBindTexture(GL_TEXTURE_2D, id);
	else
		glBindTexture(GL_TEXTURE_2D, 0);
}


void CTexture::freeTexture()
{
	if (_isInit)
	{
		glDeleteTextures(1, &id);
		_isInit = false;
	}
}


CTexture::CTexture(void)
{
	size[0] = 0;
	size[1] = 0;

	map = NULL;
	_isInit = false;
}

CTexture::~CTexture(void)
{
	freeTexture();
	if (map != NULL) delete[]   map;
}



bool CTexture::isInit()
{
	return _isInit;
}


void CTexture::loadFile(const char *fileName)
{
	FreeImage_Initialise(TRUE);

	FREE_IMAGE_FORMAT format = FreeImage_GetFileType(fileName, 0);
	if (format == FIF_UNKNOWN)
		format = FreeImage_GetFIFFromFilename(fileName);
	if ((format == FIF_UNKNOWN) || !FreeImage_FIFSupportsReading(format))
		return;

	FIBITMAP* img = FreeImage_Load(format, fileName);
	if (img == NULL)
		return;

	FIBITMAP* tempImg = img;
	img = FreeImage_ConvertTo32Bits(img);
	FreeImage_Unload(tempImg);

	size[0] = FreeImage_GetWidth(img);
	size[1] = FreeImage_GetHeight(img);

	//BGRA a RGBA
	map = new char[4 * size[0] * size[1]];
	char *buff = (char*)FreeImage_GetBits(img);

	for (unsigned int j = 0; j<size[0] * size[1]; j++){
		map[j * 4 + 0] = buff[j * 4 + 2];
		map[j * 4 + 1] = buff[j * 4 + 1];
		map[j * 4 + 2] = buff[j * 4 + 0];
		map[j * 4 + 3] = buff[j * 4 + 3];
	}

	FreeImage_Unload(img);

	FreeImage_DeInitialise();
}

CTexture &CTexture::operator =(const CTexture &m)
//CTexture::CTexture(const CTexture &m)
{
	size[0] = m.size[0];
	size[1] = m.size[1];

	freeTexture();

	if (map != NULL) {
		delete[]   map;
		map = NULL;
	}

	if (m.map != NULL) {
		map = new char[m.size[0] * m.size[1] * 4];
		memcpy(map, m.map, m.size[0] * m.size[1] * 4);
	}

	return *this;
}
