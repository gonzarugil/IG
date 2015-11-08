#pragma once
class CFBO
{
public:
	CFBO(void);
	~CFBO(void);

	void resize(unsigned int w, unsigned int h);
	void activate();
	void deactivate();

	unsigned int getColorTexId();


protected:
	unsigned int id;

	unsigned int colorTex;
	
};

