#pragma once
class CFBO
{
public:
	CFBO(unsigned int w,unsigned int h);
	~CFBO(void);

	void resize(unsigned int w, unsigned int h);
	void activate();
	void deactivate();

	unsigned int getColorTexId();


protected:
	unsigned int id;

	unsigned int colorTex;
	
};

