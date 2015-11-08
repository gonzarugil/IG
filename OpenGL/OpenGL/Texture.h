#pragma once


class CTexture
{
public:
	CTexture(void);
	~CTexture(void);

	void buildTexture();
	bool isInit();

	void freeTexture();
	void loadFile(const char *fileName);

	CTexture &operator =(const CTexture &m);

	void activate();



protected:
	unsigned int size[2];
	char *map;
	unsigned int id;
	bool _isInit;
};


