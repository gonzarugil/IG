// (c) by Stefan Roettger

#ifndef DDSBASE_H
#define DDSBASE_H

#include <stdio.h>
#include <stdlib.h>


#define ERRORMSG() errormsg(__FILE__,__LINE__)

inline void errormsg(char *file,int line)
{
	fprintf(stderr,"fatal error in <%s> at line %d!\n",file,line);
	exit(EXIT_FAILURE);
}

#define BOOLINT char
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif



void writeDDSfile(char *filename,unsigned char *data,unsigned int bytes,unsigned int skip=0,unsigned int strip=0,int nofree=0);
unsigned char *readDDSfile(char *filename,unsigned int *bytes);

void writeRAWfile(char *filename,unsigned char *data,unsigned int bytes,int nofree=0);
unsigned char *readRAWfile(char *filename,unsigned int *bytes);

void writePNMimage(char *filename,unsigned char *image,unsigned int width,unsigned int height,unsigned int components,int dds=0);
unsigned char *readPNMimage(char *filename,unsigned int *width,unsigned int *height,unsigned int *components);

void writePVMvolume(char *filename,unsigned char *volume,
                    unsigned int width,unsigned int height,unsigned int depth,unsigned int components=1,
                    float scalex=1.0f,float scaley=1.0f,float scalez=1.0f,
                    unsigned char *description=NULL,
                    unsigned char *courtesy=NULL,
                    unsigned char *parameter=NULL,
                    unsigned char *comment=NULL);

unsigned char *readPVMvolume(char *filename,
                             unsigned int *width,unsigned int *height,unsigned int *depth,unsigned int *components=NULL,
                             float *scalex=NULL,float *scaley=NULL,float *scalez=NULL,
                             unsigned char **description=NULL,
                             unsigned char **courtesy=NULL,
                             unsigned char **parameter=NULL,
                             unsigned char **comment=NULL);

unsigned int checksum(unsigned char *data,unsigned int bytes);

void swapbytes(unsigned char *data,unsigned int bytes);
void convbytes(unsigned char *data,unsigned int bytes);
void convfloat(unsigned char *data,unsigned int bytes);

unsigned char *quantize(unsigned char *volume,
                        unsigned int width,unsigned int height,unsigned int depth,
                        BOOLINT nofree=FALSE,
                        BOOLINT linear=FALSE,
                        BOOLINT verbose=FALSE);

#endif
