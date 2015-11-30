#version 430
#define MAX_ITERATIONS 4000

in vec3 entryPoint;

uniform sampler3D dataTex;
uniform sampler2D exitPointTex;
uniform ivec3 texSize;

uniform ivec2      screenSize;
uniform float     stepSize;
float isoValue  = 0.1;

//propiedades de la luz
/*
vec3 Ia = vec3(0.3);
vec3 Id = vec3(1.0);
vec3 Is = vec3(1.0);
vec3 lpos = vec3(0.0);
float alpha = 5000.0;

//propiedades del material (HARDCODEADOS)
vec3 Ka = vec3(1, 0, 0);
vec3 Kd = vec3(0, 1, 0);
vec3 Ks = vec3(1.0);
*/
//matriz para el punto
uniform mat4 model;
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
	res = vec3(0.0);
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

vec4 shade(in vec3 norm, in vec3 colorpos){
	/*vec3 c = vec3(0.0);
	c = Ia * Ka;
	
	vec3 L = normalize(lpos - colorpos);
	vec3 diffuse = Id * Kd * dot(L, norm);
	c += clamp(diffuse, 0.0, 1.0);

	vec3 V = normalize(-colorpos);
	vec3 R = normalize(reflect(-L, norm));
	float factor = max(dot(R, V), 0.001);
	vec3 specular = Is * Ks * pow(factor, alpha);
	c += clamp(specular, 0.0, 1.0);
	
	return vec4(c, 0.0f);*/
	return vec4(norm, 1);
}

void computeDensity(in vec3 voxelId, out float d0, out float d1, out float d2,
	out float d3, out float d4, out float d5, out float d6, out float d7)
{
	ivec3 v0;
	v0 = ivec3(voxelId);
	d0 = texelFetch(dataTex, v0, 0).r;
	d1 = texelFetch(dataTex, v0 + ivec3(0, 0, 1), 0).r;
	d2 = texelFetch(dataTex, v0 + ivec3(1, 0, 1), 0).r;
	d3 = texelFetch(dataTex, v0 + ivec3(1, 0, 0), 0).r;
	d4 = texelFetch(dataTex, v0 + ivec3(0, 1, 0), 0).r;
	d5 = texelFetch(dataTex, v0 + ivec3(0, 1, 1), 0).r;
	d6 = texelFetch(dataTex, v0 + ivec3(1, 1, 1), 0).r;
	d7 = texelFetch(dataTex, v0 + ivec3(1, 1, 0), 0).r;
}

bool checkIsoValue(in float d0, in float d1, in float d2, in float d3, in float d4,
				   in float d5, in float d6, in float d7, in float isovalue )
{
	
	float v1sig = d0-isoValue;
	if ((v1sig*(d1 - isovalue)) < 0.0) return true;
	if ((v1sig*(d2 - isovalue)) < 0.0) return true;
	if ((v1sig*(d3 - isovalue)) < 0.0) return true;
	if ((v1sig*(d4 - isovalue)) < 0.0) return true;
	if ((v1sig*(d5 - isovalue)) < 0.0) return true;
	if ((v1sig*(d6 - isovalue)) < 0.0) return true;
	if ((v1sig*(d7 - isovalue)) < 0.0) return true;
	return false;
	
}


void main()
{
	mat4 vmMatrix = view*model;
	mat4 vmMatrixInverse = inverse(vmMatrix);
	mat4 normalMat = transpose(inverse(vmMatrix));

	vec3 texSizeF = texSize;
	vec3 texSizeInv = 1.0 / texSizeF;
	
	vec4 color = vec4(0,0,0,0);
	vec2 screenCoords = gl_FragCoord.st/screenSize;
    vec3 exitPoint = texture(exitPointTex, screenCoords).xyz;
    
	if (entryPoint.z < exitPoint.z)
	{
		discard;
	}else{
		
		vec3 pto = (vmMatrixInverse * vec4(entryPoint, 1)).xyz;
		vec3 ptofin = (vmMatrixInverse * vec4(exitPoint, 1)).xyz;

		//TODO: Pasar texSize a float
		pto *= texSizeF;
		ptofin *= texSizeF;
		vec3 v = ptofin - pto;
		//TODO: VER LO QUE PASAS CUANDO COMPONENTES DE V SON PROXIMAS A 0
		vec3 vInv = 1.0 / v;

		vec3 incDir = sign(v);
		vec3 singDir = (sing(v) + vec3(1))*0.5;
		vec3 voxelId = floor(pto);
		vec3 voxelCenter;
		vec3 texCoords;
		

        //float tfin = (ptofin.z-pto.z)/(v.z +0.000001f); //esto es para evitar que se cuelgue
		float tfin = 1;
        float t = 0; // t del entrypoint es 0 por la forma en la que se crea la recta---

		int i = 0;
		//aqui va el while
		while (t<tfin&&i<MAX_ITERATIONS){
			i++;

			//calculo del siguiente, esto debe hacerse antes ya que se necesita el valor de tnext
			float tnext;
			vec3 newvoxelId;
			vec3 sig = voxelId + incDir;
			vec3 tsig = (sig - pto)*vInv;
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
			//voxelCenter = voxelId + vec3(0.5);


			float d0, d1, d2, d3, d4, d5, d6, d7;
			computeDensity(voxelId, d0, d1, d2, d3, d4, d5, d6, d7);
			if (checkIsoValue(d0, d1, d2, d3, d4, d5, d6, d7, isoValue)){ //Comprobacion de que el isovalor que buscamos se halla dentro del voxel
				color = vec4(1);
				break;

				//Cálculo de las densidades en los vértices multiplico la z por -1 ya que en coordenadas de textura las z son positivas
				/* con arrays
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
				*/
				//con variables sueltas
				vec3 v0, v1, v2, v3, v4, v5, v6, v7;
				float d0, d1, d2, d3, d4, d5, d6, d7;
				v0 = vec3(voxelId.x - 1, voxelId.y - 1, voxelId.z - 1);
				d0 = texture(dataTex, v0 * vec3(1, 1, -1) / vec3(texSize)).r;
				v1 = vec3(voxelId.x - 1, voxelId.y - 1, voxelId.z);
				d1 = texture(dataTex, v1 * vec3(1, 1, -1) / vec3(texSize)).r;
				v2 = vec3(voxelId.x, voxelId.y - 1, voxelId.z);
				d2 = texture(dataTex, v2 * vec3(1, 1, -1) / vec3(texSize)).r;
				v3 = vec3(voxelId.x, voxelId.y - 1, voxelId.z - 1);
				d3 = texture(dataTex, v3 * vec3(1, 1, -1) / vec3(texSize)).r;
				v4 = vec3(voxelId.x - 1, voxelId.y, voxelId.z - 1);
				d4 = texture(dataTex, v4 * vec3(1, 1, -1) / vec3(texSize)).r;
				v5 = vec3(voxelId.x - 1, voxelId.y, voxelId.z);
				d5 = texture(dataTex, v5 * vec3(1, 1, -1) / vec3(texSize)).r;
				v6 = vec3(voxelId.x, voxelId.y, voxelId.z);
				d6 = texture(dataTex, v6 * vec3(1, 1, -1) / vec3(texSize)).r;
				v7 = vec3(voxelId.x, voxelId.y, voxelId.z - 1);
				d7 = texture(dataTex, v7 * vec3(1, 1, -1) / vec3(texSize)).r;
				//Calculo de los coeficientes de la interpolacion (COEFICIENTES DE MARCOS)
				/*
				float A = d6 + d0 + d4 + d2 + d7 + d3 + d1 + d5 - 8.0*isoValue;
				float Z = d3 + d5 + d1 - d4 - d2 + d7 - d0 - d6;
				float X = d7 + d6 - d1 + d4 - d0 - d3 - d2 + d5;
				float Y = d6 - d0 - d1 - d4 + d3 - d5 + d2 + d7;
				float XYZ = d4 + d2 - d6 + d1 - d0 - d5 - d3 + d7;
				float XY = d6 - d3 + d1 - d2 - d5 - d4 + d0 + d7;
				float YZ = d7 + d0 - d5 - d2 + d3 - d6 - d1 + d4;
				float ZX = d0 - d1 + d7 + d2 - d6 + d5 - d4 - d3;
				*/
				//Calculo de los coeficientes (COEFICIENTES DE LOIC)
				//float A = d0 - 8.0*isoValue;
				float A = d0 - isoValue;
				float X = d3 - d0;
				float Y = d4 - d0;
				float Z = d1 - d0;
				float XY = d7 + d0 - d3 - d4;
				float ZX = d2 + d0 - d4 - d1;
				float YZ = d5 + d0 - d4 - d3;
				float XYZ = d6 + d3 + d4 + d1 - d0 - d5 - d2 - d7;
				
				//Repetir multiplicaciones
				/*vec3 ptotrans = pto - voxelId; //translacion del pto de la recta para encuadrarlo todo en el origen
				vec4 a;
				a.x = XYZ*v.x*v.y*v.z;
				a.y = XYZ*(ptotrans.z*v.y*v.x+v.z*ptotrans.y*v.x+v.z*v.y*ptotrans.x)+XY*v.y*v.x+ZX*v.z*v.x+YZ*v.z*v.y;
				a.z = XY*(ptotrans.y*v.x+v.y*ptotrans.x)+YZ*(ptotrans.z*v.y+v.z*ptotrans.y)+XYZ*(v.z*ptotrans.y*ptotrans.x+ptotrans.z*ptotrans.y*v.x+ptotrans.z*v.y*ptotrans.x)+ZX*(ptotrans.z*v.x+v.z*ptotrans.x)+Y*v.y+Z*v.z+X*v.x;
				a.w = A+XY*ptotrans.y*ptotrans.x+ZX*ptotrans.z*ptotrans.x+YZ*ptotrans.z*ptotrans.y+XYZ*ptotrans.z*ptotrans.y*ptotrans.x+X*ptotrans.x+Z*ptotrans.z+Y*ptotrans.y;
				vec3 sol;
				int numsoluciones = solveCubic(a,sol);
				//Comprobamos el numero de soluciones válidas de la ecuación y obtenemos la menor
				
				//comprobamos las soluciones obtenidas
				bool validsolution = true;
				float s;
				if (numsoluciones == 1){
					s = sol.x;
					if (s >= tnext) validsolution = false;
				}
				else if (numsoluciones == 3){
					bvec3 GT = greaterThan(sol, vec3(t));
					sol = sol*vec3(GT) + tnext*vec3(not(GT));
					s = min(min(sol.x, sol.y), sol.z);
					if (s >= tnext) validsolution = false; //FALLA AQUI EN ESTA COMPROBACIÓN!!!!!!! (Puede ser por la escala numérica?)(Normalizar?)
				}
				if (validsolution){
					//Calculo del punto
					vec3 colorpos = ptotrans + v*s;
					//Devolvemos el punto a coordenadas de la camara (esta en coordenadas del mundo escaladas)
					colorpos = colorpos + voxelId;
					colorpos = colorpos / texSize; //translacion y escalado
					colorpos = (vmMatrix * vec4(colorpos,1)).xyz; //coordenadas de la camara

					//Calculo de normales mediante las derivadas parciales
					vec3 norm = vec3(0.0f);
					norm.x = X + XY * ptotrans.y + ZX * ptotrans.z + XYZ * ptotrans.y * ptotrans.z;
					norm.y = Y + XY * ptotrans.x + YZ * ptotrans.z + XYZ * ptotrans.x * ptotrans.z;
					norm.z = Z + ZX * ptotrans.x + YZ * ptotrans.y + XYZ * ptotrans.x * ptotrans.y;
					norm = norm * -1;
					norm = norm / texSize;// solo escalado porque es un vector y no es homogeneo
					norm = (normalMat * vec4(norm,0)).xyz; //coordenadas de la camara
					norm = normalize(norm);
					color = shade(norm,colorpos);
					break;
				//}*/
			}
			
			voxelId = newvoxelId;
			t = tnext;
			
		}
			
	}
	outColor = color;
}