#define GLM_FORCE_RADIANS
#define _CRT_SECURE_NO_WARNINGS
#include "Dependencies\glew\glew.h"
#include "Dependencies\freeglut\freeglut.h"
#include "Dependencies\freeimage\FreeImage.h"
#include <iostream>
#include <fstream>
#include <windows.h>
#include "Texture.h"
#include "FBO.h"
#include "Dependencies\glm\glm\glm.hpp"
#include "Dependencies\glm\glm\gtc\type_ptr.hpp"
#include "Dependencies\glm\glm\gtc\matrix_transform.hpp"
#include <math.h>
#include "preprocessor.h"
#include "ddsbase.h"




CTexture texture;
CFBO *fbo;

//Performance
unsigned int fpsCurrent = 0;
unsigned int fpsCount = 0;

//Camara
// angle of rotation for the camera direction
float angle = 0.0;
// actual vector representing the camera's direction
float lx = 0.0f, lz = -1.0f;
// XZ position of the camera
float x = 0.0f, z = 5.0f;

//estructura que define las luces
struct Light {
	glm::vec3 position;
	glm::vec3 intensities; //a.k.a. the color of the light
	float attenuation; 
	float ambientCoefficient; 
};

struct Material {
	const float shininess = 80.0;
	const glm::vec3 specularColor = glm::vec3(1.0f, 1.0f, 1.0f);
};

Light gLight;
Material gMaterial;

GLfloat vertexData[] = {
	//  X     Y     Z       U     V          Normal
	// bottom
	0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
	1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, -1.0f, 0.0f,
	0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f,

	1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, -1.0f, 0.0f,
	1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f, -1.0f, 0.0f,
	0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f,

	// top
	0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,
	0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f,
	1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
	
	1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
	0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f,
	1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,
	

	// front
	0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f,
	0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
	1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,

	1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
	0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
	1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,

	// back
	
	
	
	1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,
	0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, -1.0f,
	0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
	
	
	1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,
	1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, -1.0f,
	0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, -1.0f,
	

	// left
	
	
	
	0.0f, 1.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 1.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f,

	
	
	0.0f, 1.0f, 1.0f, 1.0f, 1.0f, -1.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 1.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f,

	// right
	1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
	1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f,
	1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,

	1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
	1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,
	1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f
};


//lista de indices de los triangulos para las caras

const int boxNTriangleIndex = 36;

const unsigned int BoxTriangleIndex[] = {
	//Cara bottom
	0, 1, 2, 3, 4, 5,
	//Cara top
	6, 7, 8, 9, 10, 11,
	//Cara front
	12, 13, 14, 15, 16, 17,
	//Cara back
	18, 19, 20, 21, 22, 23,
	//Cara left
	24, 25, 26, 27, 28, 29,
	//Cara right
	30, 31, 32, 33, 34, 35,
};

//Datos para la textura 3D
GLfloat profcoord = 0;
GLint coordtex = -1;
GLuint texId[1];
GLfloat angulo = 274.0f;

GLint texLoc = -1;
//GLint isoLoc = -1;
//GLint sizeLoc = -1;

GLuint w3d = 512;
GLuint h3d= 512;
GLuint d3d= 134;

GLfloat iso = 0.20f;

//aqui va el archivo pvm que cargaremos (lo de arriba es el tamaño)
//const char *fileName = "MRI-Head.pvm";
const char *fileName = "VisMale.pvm";
//const char *fileName = "Baby.pvm";
//const char *fileName = "Orange.pvm";
//const char *fileName = "Porsche.pvm";
//const char *fileName = "Lobster.pvm";

//VAO
GLuint vao;
//GLuint fboVao;

//VBO
GLuint vbo[2];
#define vbuffer  vbo[0] //Vertex Data
#define tribuffer vbo[1] //Indices de los triangulos

//coordenadas de la pantalla
GLsizei scrw;
GLsizei scrh;

//Shaders
GLuint vshader = 0;
GLuint fshader = 0;
GLuint v2shader = 0;
GLuint f2shader = 0;
GLuint program = 0;
GLuint program2 = 0;

//Variables Uniform 
GLint uProy = -1;
GLint uView = -1;
GLint uModel = -1;
GLint u2Proy = -1;
GLint u2View = -1;
GLint u2Model = -1;
GLint scrSize = -1;
GLint stpSize = -1;
GLint textureSize = -1;
GLint isoValue = -1;

GLint uTex = -1;



//Atributos
GLint inVertex = -1;
GLint in2Vertex = -1;


//Matrices
glm::mat4 proy(1.0f);

glm::mat4 view ( 1.0f, 0.0f, 0.0f, 0.0f,
0.0f, 1.0f, 0.0f, 0.0f,
0.0f, 0.0f, 1.0f, 0.0f,
0.0f, 0.0f, -10.0f, 1.0f );

glm::mat4 model(1.0f);
glm::mat4 model2(1.0f);



glm::mat4 modelstatic(1.0f);





unsigned char test[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
						 0,0,0,0,255,0,0,0,0,
						 0, 0, 0, 0, 0, 0, 0, 0, 0};

void loadTexture(unsigned char **texels, unsigned int &w, unsigned int &h, unsigned int &d, const char *fileName)
{
	unsigned int ww, hh, dd;
	/**/
	delete[](*texels);

	
	*texels = readPVMvolume((char *)fileName, &ww, &hh, &dd);
	
	w = ww;
	h = hh;
	d = dd;

	model2[1].y = (float)h / (float)w;
	model2[2].z = (float)d / (float)w;

	/**//**
	*texels = test;

	w = 3;
	h = 3;
	d = 3;/**/

	return;

}

void loadTextureRaw(unsigned char *texels, unsigned int w, unsigned int h, unsigned int d, const char *fileName)
{
	unsigned int size = w * h * d;
	memset(texels, 0, size);

	std::fstream file(fileName, std::ios_base::in | std::ios::binary);
	if (!file.is_open())
		return;

	//Comprobacion de tamaño
	file.seekg(0, std::ios::end);
	if (size != file.tellg())
		return;

	//Vuelves a poner el cursor al principio
	file.seekg(0, std::ios::beg);
	file.read((char *)texels, (int)size);
	if (!file.good())
	{
		file.close();
		return;
	}
	file.close();
}

void countFPS(int value) {
	char title[120];
	fpsCount = fpsCurrent;
	fpsCurrent = 0;

	sprintf(title, "IsoSurface Render. FPS: %d.", fpsCount);
	glutSetWindowTitle(title);
	glutTimerFunc(1000, countFPS, 1);
}

//Calcula la matriz de proyección
void buildProjectionMatrix(float fov, float ratio, float nearPlane, float farPlane)
{
	float f = 1.0f / tan(fov * (3.141599f / 360.0f));

	proy[0][0] = f / ratio;
	proy[1][1] = f;
	proy[2][2] = (farPlane + nearPlane) / (nearPlane - farPlane);
	proy[3][2] = (2.0f * farPlane * nearPlane) / (nearPlane - farPlane);
	proy[2][3] = -1.0f;
	proy[3][3] = 0.0f;
}
//funcion que Carga los ficheros de los shaders y los crea
GLuint loadShader(const char *fileName, GLenum type)
{
	//Se carga el fichero
	std::ifstream file;
	file.open(fileName, std::ios::in);
	if (!file) return 0;

	//Se calcula la longitud del fichero
	file.seekg(0, std::ios::end);
	unsigned int fileLen = file.tellg();
	file.seekg(std::ios::beg);

	//Se lee el fichero
	char *source = new char[fileLen + 1];

	int i = 0;
	while (file.good())
	{
		source[i] = file.get();
		if (!file.eof()) i++;
		else fileLen = i;
	}
	source[fileLen] = '\0';
	file.close();

	//////////////////////////////////////////////
	//Creación y compilación del Shader
	GLuint shader;
	shader = glCreateShader(type);
	glShaderSource(shader, 1, (const GLchar **)&source, (const GLint *)&fileLen);
	glCompileShader(shader);
	delete source;

	//Comprobamos que se compilo bien
	GLint compiled;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
	if (!compiled)
	{
		//Calculamos una cadena de error
		GLint logLen;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLen);

		char *logString = new char[logLen];
		glGetShaderInfoLog(shader, logLen, NULL, logString);
		std::cout << "Error: " << logString << std::endl;
		delete logString;

		glDeleteShader(shader);
		return 0;
	}

	return shader;
}
//inicializa los shaders y los linkea (vshader ->vertices)(fshader->fragmentos)
void shaderInit()
{
	//Creamos los shaders
	const char vname[] = "vExitpoints.glsl";
	const char fname[] = "fExitpoints.glsl";

	const char v2name[] = "vvoxelshader.glsl";
	const char f2name[] = "ISO_DATA_S2OUTFrag.glsl";
	

	//Compilación
	
	vshader = loadShader(vname, GL_VERTEX_SHADER);
	fshader = loadShader(fname, GL_FRAGMENT_SHADER);
	v2shader = loadShader(v2name, GL_VERTEX_SHADER);
	f2shader = loadShader(f2name, GL_FRAGMENT_SHADER);
	

	//Link
	program = glCreateProgram();
	program2 = glCreateProgram();

	

	glAttachShader(program, vshader);
	glAttachShader(program, fshader);

	glAttachShader(program2, v2shader);
	glAttachShader(program2, f2shader);
	

	glBindAttribLocation(program, 0, "inVertex");
	glLinkProgram(program);

	glBindAttribLocation(program2, 0, "inVertex");
	glLinkProgram(program2);
	

	int linked2;
	glGetProgramiv(program2, GL_LINK_STATUS, &linked2);
	int linked;
	glGetProgramiv(program, GL_LINK_STATUS, &linked);
	if (!linked || !linked2)
	{
		//Calculamos una cadena de error
		GLint logLen;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLen);

		char *logString = new char[logLen];
		glGetProgramInfoLog(program, logLen, NULL, logString);
		std::cout << "Error: " << logString << std::endl;
		delete logString;

		glDeleteProgram(program);
		program = 0;
		return;
	}
	//Variable uniform
	uProy = glGetUniformLocation(program, "proy");
	uView = glGetUniformLocation(program, "view");
	uModel = glGetUniformLocation(program, "model");

	u2Proy = glGetUniformLocation(program2, "proy");
	u2View = glGetUniformLocation(program2, "view");
	u2Model = glGetUniformLocation(program2, "model");
	uTex = glGetUniformLocation(program2, "exitPointTex");
	texLoc = glGetUniformLocation(program2, "dataTex");
	scrSize = glGetUniformLocation(program2, "screenSize");
	stpSize = glGetUniformLocation(program2, "stepsize");
	textureSize = glGetUniformLocation(program2, "texSize");
	isoValue = glGetUniformLocation(program2, "isoValue");

	
	//Atributos
	inVertex = glGetAttribLocation(program, "inVertex");
	in2Vertex = glGetAttribLocation(program2, "inVertex");

	
	
}
//libera memoria al finalizar la ejecucion
void destroy()
{
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glDeleteBuffers(2, vbo);
	

	glBindVertexArray(0);
	glDeleteVertexArrays(1, &vao);
	

	glDetachShader(program, vshader);
	glDetachShader(program, fshader);
	
	glDeleteShader(vshader);
	glDeleteShader(fshader);
	
	glDeleteProgram(program);
	
}


//Inicializa y configura la escena
void sceneInit(){


	glEnable(GL_DEPTH_TEST);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	buildProjectionMatrix(45.0f, 4.0f / 3.0f, 0.1f, 50.0f);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_CULL_FACE);

	//glDisable(GL_CULL_FACE);

	glGenVertexArrays(1, &vao);
	//glGenVertexArrays(1, &fboVao);
	glBindVertexArray(vao);

	//Buffers
	glGenBuffers(2, vbo);
	
	glBindBuffer(GL_ARRAY_BUFFER, vbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexData), vertexData, GL_STATIC_DRAW);
	

	glVertexAttribPointer(inVertex, 3, GL_FLOAT, GL_FALSE, sizeof(float)*8, 0);
	glEnableVertexAttribArray(inVertex);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, tribuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(BoxTriangleIndex), BoxTriangleIndex, GL_STATIC_DRAW);

	glBindVertexArray(0);


	unsigned char *texels = new unsigned char[w3d*h3d*d3d];
	loadTexture(&texels, w3d, h3d, d3d, fileName);
	//glEnable(GL_TEXTURE_3D);
	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, texId);
	glBindTexture(GL_TEXTURE_3D, texId[0]);
	//GL_MAX_3D_TEXTURE_SIZE

	GLint max3dtexsize;
	glGetIntegerv(GL_MAX_3D_TEXTURE_SIZE, &max3dtexsize);
	std::cout << "GL_MAX_3D_TEXTURE_SIZE: " << max3dtexsize << std::endl;

	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);


	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage3D(GL_TEXTURE_3D, 0, GL_R8, w3d, h3d, d3d, 0, GL_RED,
		GL_UNSIGNED_BYTE, texels);

	//TODO: Volver a poner
	delete[]texels;

}

void keyPressed(unsigned char key, int x, int y) {
	float fraction = 0.1f;
	float incV = 0.05f;
	float inciso = 0.01f;
	float angulo = 0.01f;

	switch (key) {
	case 'r':
		view = glm::rotate(view, angulo, glm::vec3(1, 0, 0));
		break;
	case 'f':
		view = glm::rotate(view, -angulo, glm::vec3(1, 0, 0));
		break;
	case 't':
		view = glm::rotate(view, angulo, glm::vec3(0, 0, 1));
		break;
	case 'g':
		view = glm::rotate(view, -angulo, glm::vec3(0, 0, 1));
		break;
	case 'q':
		view = glm::rotate(view, angulo, glm::vec3(0, 1, 0));
		break;
	case 'e':
		view = glm::rotate(view, -angulo, glm::vec3(0, 1, 0));
		break;
	case 'z':
		iso -= inciso;
		break;
	case 'x':
		iso += inciso;
		break;
	case 'w':
		view[3].y += incV;
		break;
	case 'a':
		view[3].x += incV;
		break;
	case 's':
		view[3].y -= incV;
		break;
	case 'd':
		view[3].x -= incV;
		break;
	case '+':
		view[3].z += incV;
		break;
	case '-':
		view[3].z -= incV;
		break;
	}
}
void processSpecialKeys(int key, int xx, int yy) {

	float fraction = 0.1f;
	float incV = 0.05f;

	switch (key) {
	case GLUT_KEY_ALT_L:
		iso = iso + 0.05;
		break;
	case GLUT_KEY_ALT_R:
		iso = iso - 0.05;
		break;
	case GLUT_KEY_UP:
		view[3].y += incV;
		std::cout << "a";
		//iso = iso + 0.01;
		break;
	case GLUT_KEY_DOWN:
		view[3].y -= incV;
		std::cout << "a";
		//iso = iso - 0.01;
		break;
	case GLUT_KEY_LEFT:
		view[3].x += incV;
		std::cout << "a";
		break;
	case GLUT_KEY_RIGHT:
		view[3].x -= incV;
		break;
	case GLUT_KEY_PAGE_DOWN:
		view[3].z += incV;
		std::cout << "a";
		break;
	case GLUT_KEY_PAGE_UP:
		view[3].z -= incV;
		break;
	}
}

//DisplayFunc
void renderScene(void){
	//Cargar primer programa (primera pasada)

	fbo->activate();
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glDepthFunc(GL_GREATER);
	glFrontFace(GL_CW);
	glUseProgram(program);
	
	glUniformMatrix4fv(uModel, 1, GL_FALSE, glm::value_ptr(model));
	glUniformMatrix4fv(uView, 1, GL_FALSE, glm::value_ptr(view));
	glUniformMatrix4fv(uProy, 1, GL_FALSE, glm::value_ptr(proy));
	

	//Pintado del buffer    
	glBindVertexArray(vao);
	glDrawElements(GL_TRIANGLES, boxNTriangleIndex, GL_UNSIGNED_INT, (void*)0);
	
	fbo->deactivate();
	/**/
	//Cargar Segundo programa (Segunda pasada)
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(program2);
	glFrontFace(GL_CCW);

	//Variables
	glUniformMatrix4fv(u2Model, 1, GL_FALSE, glm::value_ptr(model));
	glUniformMatrix4fv(u2View, 1, GL_FALSE, glm::value_ptr(view));
	glUniformMatrix4fv(u2Proy, 1, GL_FALSE, glm::value_ptr(proy));
	glUniform2i(scrSize,scrw,scrh);
	glUniform1f(stpSize, 0.01f);
	glUniform3i(textureSize, w3d, h3d, d3d);
	glUniform1f(isoValue, iso);

	

	glActiveTexture(GL_TEXTURE0);
	glUniform1i(uTex, 0);
	glBindTexture(GL_TEXTURE_2D, fbo->getColorTexId());

	glActiveTexture(GL_TEXTURE0+1);
	glBindTexture(GL_TEXTURE_3D, texId[0]);
	glUniform1i(texLoc, 1);

	//Pintado del cubo
	glBindVertexArray(vao);
	glDrawElements(GL_TRIANGLES, boxNTriangleIndex, GL_UNSIGNED_INT, (void*)0);/**/
	glLoadIdentity();
	gluLookAt(x, 1.0f, z,
		x + lx, 1.0f, z + lz,
		0.0f, 1.0f, 0.0f);


	glUseProgram(NULL);
	glutSwapBuffers();

	fpsCurrent++;
}
//ReshapeFunc
void funcionDeReescalado(GLsizei w, GLsizei h)
{
	glViewport(0, 0, w, h);
	buildProjectionMatrix(45.0f, 4.0f / 3.0f, 0.1f, 50.0f);
	//fbo->resize(w,h); cambio esto porlo siguiente para no reescribir en el fbo sino crear uno nuevo
	fbo->~CFBO();
	fbo = new CFBO(w, h);
	scrw = w;
	scrh = h;
}
//IdleFunc
void funcionIdle(){
	
	angulo = 0.0f;
	angulo = (angulo<3.141599f*2.0f) ? angulo + 0.003f : 0.0f;

	
	 model = glm::rotate(model, angulo, glm::vec3(1, 0, 0));
	 model = glm::rotate(model, angulo, glm::vec3(0, 1, 0));
	 
	 model = glm::mat4(1);
	 model = glm::scale(model, glm::vec3(3, 3, 3));

	 model = model * model2;
	 model = glm::translate(model, glm::vec3(-0.333f));




	//Sleep(5);
	glutPostRedisplay();
}


int main(int argc, char **argv)
{
	
	glutInit(&argc, argv);
	glutInitContextVersion(4, 3);
	glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
	//glutInitContextProfile(GLUT_CORE_PROFILE);
	glutInitContextProfile(GLUT_COMPATIBILITY_PROFILE);

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(512, 512);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("glsl3");
	glutReshapeFunc(funcionDeReescalado);
	glutDisplayFunc(renderScene);
	glutIdleFunc(funcionIdle);
	glutTimerFunc(1000, countFPS, 1);
	glutKeyboardFunc(keyPressed);
	glutSpecialFunc(processSpecialKeys);

	//Esperamos a que se cree el contexto de OpenGL
	glewExperimental = GL_TRUE;
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		printf("Error: %s\n", glewGetErrorString(err));
	}
	const GLubyte *oglVersion = glGetString(GL_VERSION);
	//CFBO fboaux;
	fbo = new CFBO(512,512);// fboaux;

	printf("This system supports OpenGL Version %s.\n", oglVersion);
	shaderInit();
	sceneInit();

	glutMainLoop();
	destroy();
	delete fbo;
	return 0;

}