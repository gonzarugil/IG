#version 440
#define MAX_ITERATIONS 4000

in vec3 entryPoint;

uniform sampler3D dataTex;
uniform sampler2D exitPointTex;
uniform vec3 texSize;

uniform vec2      screenSize;
uniform float     stepSize;
float isoValue  = 0.40;



//matriz para el punto
uniform mat4 rot;
uniform mat4 view;
uniform mat4 proy;

out vec4 outColor;

float cbrt(float x){
	float EPSILON = 0.00000001f;
	if (abs(x) < EPSILON) return 0.0;
	if (x> 0.0) return pow(x,0.33333333333333333333f);
	return -pow(-x,0.33333333333333333333f);
}

float computeX (float As,float Ch,float Dh,float delta){
	float T0 = sign(Dh) * abs (As) * sqrt(-delta);
	float T1 = -Dh + T0;
	float p = cbrt (T1 * 0.5f);
	// TODO: Poner epsilon?
	float q = (T1==T0)? -p : -Ch/(p);
	return (Ch<=0)? p+q: -Dh/(p*p+q*q+Ch);
}

int solveCubic(in vec4 coff, out vec3 res){
	int nsol = 0;
	coff.yz *= 0.33333333333333;
	vec3 sigma = vec3(coff.x*coff.z - coff.y*coff.y,
					coff.x*coff.w - coff.y*coff.z,
					coff.y*coff.w - coff.z*coff.z);
	float delta = 4 * sigma.x *sigma.z - sigma.y*sigma.y;
	if (delta <= 0.0f){
		if (coff.y*coff.y*coff.y*coff.w > coff.x * coff.z*coff.z*coff.z){
			//Algoritmo A
			float x = computeX (coff.x,sigma.x,-2.0f*coff.y*sigma.x + coff.x*sigma.y,delta);
			//TODO: Controlar el caso A=0 (Exacto)
			res.x = (x - coff.y)/coff.x;
			nsol += 1;
		}else{
			//ALgoritmo D
			float x = computeX (coff.w,sigma.z,-coff.w * sigma.y + 2.0f * coff.z * sigma.z,delta);
			res.x = (-coff.w/(x + coff.z));
			nsol +=1;
		}
	}else{
		float Ch = sigma.x;
		float Dh = -2.0f * coff.y * sigma.x + coff.x * sigma.y;
		float theta = 0.3333333333333333f * abs(atan(-Dh,coff.x*sqrt(delta)));
		float sqCh2 = 2.0f * sqrt (-Ch);
		float cosTheta = cos(theta);
		float sinTheta = sin(theta);
		float x1 = sqCh2 * cosTheta;
		float x3 = sqCh2 * (-0.5f * cosTheta - 0.86602540378f * sinTheta);
		float xL = (x1+x3 > 2.0f * coff.y)?x1:x3;
		res.x = (xL - coff.y)/coff.x;
		nsol += 1;
		//TODO: Controlar el caso de A = 0 (exacto) Poner Infinito
		//Sol2 xs
		Ch = sigma.z;
		Dh = -coff.w * sigma.y + 2.0f * coff.z * sigma.z;
		theta = 0.33333333333333333f * abs(atan(-Dh,coff.w*sqrt(delta)));
		sqCh2= 2.0f *sqrt(-Ch);
		cosTheta = cos(theta);
		sinTheta = sin(theta);
		x1 = sqCh2 * cosTheta;
		x3 = sqCh2 * (-0.5f * cosTheta - 0.86602540378f * sinTheta);
		float xS = (x1+x3 < 2.0f * coff.z)? x1:x3;
		res.y = -coff.w / (xS+ coff.z);
		nsol +=1;
		//The only time this will go wrong is when F^2=EG, which is when the quadratic has a double root. In
		//this case — when the largest and smallest roots of the cubic are equal — it means that we must have a
		//triple-root cubic. But in fact, this will only really cause problems if the triple root is at zero or infinity. We
		//will deal with this another time.
		//https://courses.cs.washington.edu/courses/cse590b/13au/lecture_notes/solvecubic_p5.pdf
		float E = (xS + coff.z)*coff.x;
		float F = -(xL - coff.y)*(xS + coff.z) - coff.x * (-coff.w);
		float G = -coff.w * (xL - coff.y);
		res.z = (xL - coff.y)/coff.x;
		nsol +=1;
	}
	return nsol;
}

bool checkIsoValue(in float isoValue, in vec3 voxelID){
	vec3 vertex[8];
	vertex[0] = voxelID;
	vertex[1] = vec3(voxelID.x,voxelID.y,voxelID.z-1);
	vertex[2] = vec3(voxelID.x,voxelID.y+1,voxelID.z); 
	vertex[3] = vec3(voxelID.x,voxelID.y+1,voxelID.z-1); 
	vertex[4] = vec3(voxelID.x+1,voxelID.y,voxelID.z); 
	vertex[5] = vec3(voxelID.x+1,voxelID.y,voxelID.z-1); 
	vertex[6] = vec3(voxelID.x+1,voxelID.y+1,voxelID.z); 
	vertex[7] = vec3(voxelID.x+1,voxelID.y+1,voxelID.z-1);
	
	
	float v1sig = texture(dataTex,vertex[0]*vec3(1,1,-1)/texSize).r-isoValue;
	for (int i=1; i < 8 ;i++){
		float vIsig = texture(dataTex,vertex[i]*vec3(1,1,-1)/texSize).r-isoValue;
		if (v1sig*vIsig <= 0.0){
			return true;
		}
	}
	return false;
}


void main()
{
	mat4 vmMatrix = view*rot;
	mat4 vmMatrixInverse = inverse(vmMatrix);

	vec4 color = vec4(0,0,0,0);
	vec2 screenCoords = gl_FragCoord.st/screenSize;
    vec3 exitPoint = texture(exitPointTex, screenCoords).xyz;
    
	if (entryPoint.z <= exitPoint.z)
	{
		discard;
	}else{
		
		vec3 pto   = (vmMatrixInverse * vec4 (entryPoint,1)).xyz;
        vec3 ptofin   = (vmMatrixInverse * vec4 (exitPoint,1)).xyz;
		vec3 v = ptofin -pto;

		pto *= texSize;
		v *= texSize;
		ptofin *= texSize;
		vec3 incDir = sign(v);
		vec3 voxelId = ceil(pto);
		vec3 voxelCenter;
		vec3 texCoords;

        float tfin = (ptofin.z-pto.z)/(v.z + 0.0000000001f);
        float t = 0; // t del entrypoint es 0 por la forma en la que se crea la recta---

		//aqui va el while
		while (t<tfin){
			
			//calculo del siguiente, esto debe hacerse antes ya que se necesita el valor de tnext
			float tnext;
			vec3 newvoxelId;
			vec3 sig = voxelId + incDir;
			vec3 tsig = (sig - pto)/v;
			if (tsig.x<tsig.y && tsig.x<tsig.z){
				newvoxelId = voxelId + vec3(1,0,0)*incDir;
				tnext = tsig.x;
			}
			else if (tsig.y<tsig.x && tsig.y<tsig.z){
				newvoxelId = voxelId + vec3(0,1,0)*incDir;
				tnext = tsig.y;
			}else {
				newvoxelId = voxelId + vec3(0,0,1)*incDir;
				tnext = tsig.z;			
			}

		    //Calculo del centro del voxel actual
			voxelCenter = voxelId + vec3(0.5);

			if (checkIsoValue(isoValue,voxelId)){ //Comprobacion de que el isovalor que buscamos se halla dentro del voxel
				//Cálculo de las densidades en los vértices multiplico la z por -1 ya que en coordenadas de textura las z son positivas
				vec3 vertex[8];
				float density[8];
				vertex[0] = voxelId;
				density[0] = texture(dataTex,vertex[0]*vec3(1,1,-1)/texSize).r;
				vertex[1] = vec3(voxelId.x,voxelId.y,voxelId.z-1);
				density[1] = texture(dataTex,vertex[1]*vec3(1,1,-1)/texSize).r;
				vertex[2] = vec3(voxelId.x,voxelId.y+1,voxelId.z);
				density[2] = texture(dataTex,vertex[2]*vec3(1,1,-1)/texSize).r;
				vertex[3] = vec3(voxelId.x,voxelId.y+1,voxelId.z-1);
				density[3] = texture(dataTex,vertex[3]*vec3(1,1,-1)/texSize).r;
				vertex[4] = vec3(voxelId.x+1,voxelId.y,voxelId.z);
				density[4] = texture(dataTex,vertex[4]*vec3(1,1,-1)/texSize).r;
				vertex[5] = vec3(voxelId.x+1,voxelId.y,voxelId.z-1);
				density[5] = texture(dataTex,vertex[5]*vec3(1,1,-1)/texSize).r; 
				vertex[6] = vec3(voxelId.x+1,voxelId.y+1,voxelId.z);
				density[6] = texture(dataTex,vertex[6]*vec3(1,1,-1)/texSize).r; 
				vertex[7] = vec3(voxelId.x+1,voxelId.y+1,voxelId.z-1);
				density[7] = texture(dataTex,vertex[7]*vec3(1,1,-1)/texSize).r;
				//Calculo de los coeficientes de la interpolacion
				float A   = density[6]+density[0]+density[4]+density[2]+density[7]+density[3]+density[1]+density[5] - 8.0*isoValue;
				float Z   = density[3]+density[5]+density[1]-density[4]-density[2]+density[7]-density[0]-density[6];
				float X   = density[7]+density[6]-density[1]+density[4]-density[0]-density[3]-density[2]+density[5];
				float Y   = density[6]-density[0]-density[1]-density[4]+density[3]-density[5]+density[2]+density[7];
				float XYZ = density[4]+density[2]-density[6]+density[1]-density[0]-density[5]-density[3]+density[7];
				float XY  = density[6]-density[3]+density[1]-density[2]-density[5]-density[4]+density[0]+density[7];
				float YZ  = density[7]+density[0]-density[5]-density[2]+density[3]-density[6]-density[1]+density[4];
				float ZX  = density[0]-density[1]+density[7]+density[2]-density[6]+density[5]-density[4]-density[3];
				
				//Repetir multiplicaciones
				vec3 ptotrans = pto - voxelId; //translacion del pto de la recta para encuadrarlo todo en el origen
				vec4 a;
				a.x = XYZ*v.x*v.y*v.z;
				a.y = XYZ*(ptotrans.z*v.y*v.x+v.z*ptotrans.y*v.x+v.z*v.y*ptotrans.x)+XY*v.y*v.x+ZX*v.z*v.x+YZ*v.z*v.y;
				a.z = XY*(ptotrans.y*v.x+v.y*ptotrans.x)+YZ*(ptotrans.z*v.y+v.z*ptotrans.y)+XYZ*(v.z*ptotrans.y*ptotrans.x+ptotrans.z*ptotrans.y*v.x+ptotrans.z*v.y*ptotrans.x)+ZX*(ptotrans.z*v.x+v.z*ptotrans.x)+Y*v.y+Z*v.z+X*v.x;
				a.w = A+XY*ptotrans.y*ptotrans.x+ZX*ptotrans.z*ptotrans.x+YZ*ptotrans.z*ptotrans.y+XYZ*ptotrans.z*ptotrans.y*ptotrans.x+X*ptotrans.x+Z*ptotrans.z+Y*ptotrans.y;
				vec3 sol;
				int numsoluciones = solveCubic(a,sol);
				//Comprobamos el numero de soluciones válidas de la ecuación y obtenemos la menor
				
				//comprobamos las soluciones obtenidas
				bvec3 GT=greaterThan(sol,vec3(t));
				sol = sol*vec3(GT) + tnext*vec3(not(GT));
				float s=min(min(sol.x,sol.y),sol.z);
				if (s>tnext) discard;

				color = vec4(s,0.0f,0.0f,1.0f);
				
				
			}
			
			voxelId = newvoxelId;
			t = tnext;
			
		}
			
	}
	outColor = color;
}