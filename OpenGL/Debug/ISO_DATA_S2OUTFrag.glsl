#version 420
#define MAX_ITERATIONS 3000
#define EPSILON 0.00000001f


in vec3 entryPoint;

out vec4 outColor;

uniform sampler3D dataTex;
uniform sampler2D exitPointTex;
uniform ivec2    screenSize;

//matriz para el punto
uniform mat4 model;
uniform mat4 view;
uniform mat4 proy;


uniform float isoValue;

/***********PROPIEDADES PARA ILUMINACION (HARDCODEADOS)**************/
//propiedades de la luz
vec3 Ia = vec3(0.6);
vec3 Id = vec3(0.5);
vec3 Is = vec3(0.0);
vec3 lpos = vec3(0,0,0);
//propiedades del material
vec3 Ka = vec3(1, 0.5, 0.0);
vec3 Kd = vec3(1, 1, 1);
vec3 Ks = vec3(1.0);
float alpha = 500.0;



float stepSize = 0.2;
vec4 color = vec4(1,1,1,1);
vec4 backgroundColor = vec4(0, 0, 1, 0);
uniform ivec3 texSize;
vec3 volumeDim = texSize;

float alphaStepFactor = 1.0 / stepSize;

vec3 getVoxelID(in vec3 v, in vec3 p){
	vec3 s = sign(v);
	vec3 pointAux = p + s * 0.01;
	vec3 c = floor(pointAux) + vec3(0.5);
	return c;
}



/* -------------------FUNCIONES AUXILIARES---------------------------*/

bool containsInVec3(in vec3 vector, in float value)
{
	for (int i = 0; i < 3; ++i)
	{
		if (value == vector[i])
			return true;
	}
	return false;
}

int findRepeated(in vec3 vector)
{
	if (vector.y == vector.x) return 1;
	if (vector.z == vector.x) return 2;
	if (vector.z == vector.y) return 2;

	return -1;
}


//-------------SOLVER ECUACIONES CUBICAS------------

float cbrtt(in float x) {
	if (abs(x) < EPSILON) return 0.0;
	if (x > 0.0) return pow(x, 0.33333333333333333333f);
	return -pow(-x, 0.33333333333333333333f);
}

float computeX(in float As, in float Ch, in float Dh, in float delta)
{
	float T0 = -float(sign(Dh)) * abs(As) * sqrt(-delta);
	float T1 = -Dh + T0;
	float p = cbrtt(T1*0.5);
	//TODO: Poner epsilon??
	float q = (T1 == T0) ? -p : -Ch / (p);

	return (Ch <= 0) ? p + q : -Dh / (p*p + q*q + Ch);
}




int solveCubic(in float inA, in float inB, in float inC, in float inD, out vec3 solutions)
{
	float A = inA, B = inB, C = inC, D = inD;

	/* *********************
	/if (A == 0.0)
	{
		// Ecuación cuadrada: BX^2 + CX + D
		float b2 = C*C;
		float menos = 4.0 * B*D;
		float raiz = b2 - menos;
		//if (raiz < EPSILON) raiz = 0.0;
		raiz = sqrt(raiz);
		float s1 = (-C + raiz) / (2.0 * B);
		float s2 = (-C - raiz) / (2.0 * B);
		
		solutions.x = s1;
		solutions.y = s2;
		return 2;
	}
	// *********************** */



	float Baux = B;  // ******
	float Caux = C;	// ******
	B /= 3.0;
	C /= 3.0;

	float sigma1 = A*C - B*B;
	float sigma2 = A*D - B*C;
	float sigma3 = B*D - C*C;
	float delta = 4 * sigma1*sigma3 - sigma2*sigma2;


	if (delta <= 0.0)
	{
		//if (B * B * B * D >= A * C * C * C) 
		if (B * B * B * D > A * C * C * C) // Funciona mejor para rectas?
		{
			//Algortimo A
			
			float x = computeX(A, sigma1, -2.0 * B * sigma1 + A * sigma2, delta);
			//TODO: Controlar el caso de A = 0 (exacto)
			float s1 = (x - B) / A;

			// ************** Caso de 1 simple (arriba) y 1 doble
			float r = s1;
			float e = A;
			float f = Baux + r*e;
			float g = Caux + r*f;
			// Ecuación cuadrada: eX^2 + fX + g
			float f2 = f*f;
			float menos = 4.0 * e*g;
			float raiz = f2 - menos;
			if (raiz < EPSILON) raiz = 0.0;
			raiz = sqrt(raiz);
			float s2 = (-f + raiz) / (2.0 * e);
			

			solutions.x = s1;
			solutions.y = s2;
			return 1; // retornamos 1 en lugar de 2 porque es un caso despreciable de tangencia con la superficie
		}
		else{
			//Algoritmo D
			float x = computeX(D, sigma3, -D * sigma2 + 2.0 * C * sigma3, delta);
			float s1 = -D / (x + C);

			solutions.x = s1;
			return 1;
		}

	}
	else{
		float Ch = sigma1;
		float Dh = -2.0 * B * sigma1 + A * sigma2;
		float theta = 0.33333333333333333 * abs(atan(A*sqrt(delta) , - Dh)); // TODO: en c++ era atan2 ???
		float sqCh2 = 2.0 * sqrt(-Ch);
		float cosTheta = cos(theta);
		float sinTheta = sin(theta);
		float x1 = sqCh2  * cosTheta;
		float x3 = sqCh2 * (-0.5 * cosTheta - 0.86602540378 * sinTheta);
		float xL = (x1 + x3 > 2.0 * B) ? x1 : x3;
		//TODO: Controlar el caso de A = 0 (exacto) Poner infinito
		float s11 = (xL - B) / A;
		float s12 = (x1 - B) / A;
		float s13 = (x3 - B) / A;

		//Sol2 xs
		Ch = sigma3;
		Dh = -D * sigma2 + 2.0 * C * sigma3;
		theta = 0.33333333333333333 * abs(atan(D*sqrt(delta), -Dh)); // TODO: en c++ era atan2 ???
		sqCh2 = 2.0 * sqrt(-Ch);
		cosTheta = cos(theta);
		sinTheta = sin(theta);
		x1 = sqCh2  * cosTheta;
		x3 = sqCh2 * (-0.5 * cosTheta - 0.86602540378 * sinTheta);
		float xS = (x1 + x3 < 2.0 * C) ? x1 : x3;
		float s21 = -D / (xS + C);
		float s22 = -D / (x1 + C);
		float s23 = -D / (x3 + C);

		//The only time this will go wrong is when F^2=EG, which is when the quadratic has a double root. In
		//this case — when the largest and smallest roots of the cubic are equal — it means that we must have a
		//triple-root cubic. But in fact, this will only really cause problems if the triple root is at zero or infinity. We
		//will deal with this another time.
		//https://courses.cs.washington.edu/courses/cse590b/13au/lecture_notes/solvecubic_p5.pdf
		float E = (xS + C)*A;
		float F = -(xL - B)*(xS + C) - A * (-D);
		float G = -D * (xL - B);
		float s3 = (C*F - B*G) / (-B*F + C*E);

		solutions.x = s11;
		solutions.y = s12;
		solutions.z = s13;

		int theresRepeated = findRepeated(solutions);

		while (theresRepeated > -1)
		{
			if (!containsInVec3(solutions, s21))
				solutions[theresRepeated] = s21;
			else if (!containsInVec3(solutions, s22))
				solutions[theresRepeated] = s22;
			else if (!containsInVec3(solutions, s23))
				solutions[theresRepeated] = s23;
			else if (!containsInVec3(solutions, s3))
				solutions[theresRepeated] = s3;

			theresRepeated = findRepeated(solutions);
		}

		return 3;
	}
	return 0;
}

//-----------------------------------------------------------




//---FUNCION DE CALCULO DE COLOR
vec4 shade(in vec3 norm, in vec3 colorpos){
	
	vec3 c = vec3(0.0);
	c = Ia * Ka;

	vec3 L = normalize(lpos - colorpos);
	vec3 diffuse = Id * Kd * dot(L, norm);
	c += clamp(diffuse, 0.0, 1.0);

	vec3 V = normalize(-colorpos);
	vec3 R = normalize(reflect(-L, norm));
	float factor = max(dot(R, V), 0.001);
	vec3 specular = Is * Ks * pow(factor, alpha);
	c += clamp(specular, 0.0, 1.0);

	return vec4(c, 0.0);
	return vec4(norm, 1);
}

//-- OBTIENE EL VECTOR DE DENSIDADES DEL CUBO
void computeDensity(in vec3 voxelId, out float d0, out float d1, out float d2,
	out float d3, out float d4, out float d5, out float d6, out float d7)
{
	
	d0 = texelFetch(dataTex, ivec3(voxelId + vec3(-0.5)), 0).r;
	d1 = texelFetch(dataTex, ivec3(voxelId + vec3(-0.5, -0.5, 0.5)), 0).r;
	d2 = texelFetch(dataTex, ivec3(voxelId + vec3(0.5, -0.5, 0.5)), 0).r;
	d3 = texelFetch(dataTex, ivec3(voxelId + vec3(0.5, -0.5, -0.5)), 0).r;
	d4 = texelFetch(dataTex, ivec3(voxelId + vec3(-0.5, 0.5, -0.5)), 0).r;
	d5 = texelFetch(dataTex, ivec3(voxelId + vec3(-0.5, 0.5, 0.5)), 0).r;
	d6 = texelFetch(dataTex, ivec3(voxelId + vec3(0.5)), 0).r;
	d7 = texelFetch(dataTex, ivec3(voxelId + vec3(0.5, 0.5, -0.5)), 0).r;
}

//-- COMPRUEBA QUE EL ISOVALOR PERTENEZCA AL VOXEL
bool checkIsoValue(in float d0, in float d1, in float d2, in float d3, in float d4,
	in float d5, in float d6, in float d7, in float isovalue)
{

	float v1sig = d0 - isoValue;
	if ((v1sig*(d1 - isovalue)) <= 0.0) return true;
	if ((v1sig*(d2 - isovalue)) <= 0.0) return true;
	if ((v1sig*(d3 - isovalue)) <= 0.0) return true;
	if ((v1sig*(d4 - isovalue)) <= 0.0) return true;
	if ((v1sig*(d5 - isovalue)) <= 0.0) return true;
	if ((v1sig*(d6 - isovalue)) <= 0.0) return true;
	if ((v1sig*(d7 - isovalue)) <= 0.0) return true;
	return false;

}


bool checkSolution(in vec3 solutions, in float t,in float tnext, in int nsol, out float s)
{
	bvec3 validsolution = bvec3(false);
	validsolution.x = solutions.x > t && solutions.x <= tnext;
	if (nsol > 1) validsolution.y = solutions.y > t && solutions.y <= tnext;
	if (nsol > 2) validsolution.z = solutions.z > t && solutions.z <= tnext;
	vec3 valids = solutions * vec3(validsolution);
	vec3 invalid = tnext * vec3(not(validsolution));
	vec3 sum = valids + invalid;
	
	float temp = min(min(sum.x, sum.y), sum.z);
	s = temp;
	return !((s >= tnext - EPSILON) && (s <= tnext + EPSILON)) && 
		   !((s >= t - EPSILON) && (s <= t + EPSILON));
}

bool checkSolution2(in vec3 solutions, in float t, in float tnext, in int nsol, out float s)
{
	bvec3 validsolution = bvec3(false);
	validsolution.x = solutions.x > 0 && solutions.x <= 1;
	if (nsol > 1) validsolution.y = solutions.y > 0 && solutions.y <= 1;
	if (nsol > 2) validsolution.z = solutions.z > 0 && solutions.z <= 1;
	vec3 valids = solutions * vec3(validsolution);
	vec3 invalid = 1 * vec3(not(validsolution));
	vec3 sum = valids + invalid;

	float temp = min(min(sum.x, sum.y), sum.z);
	s = temp;
	return !((s >= 1 - EPSILON) && (s <= 1 + EPSILON));
}


// OBTIENE LA ECUACIÓN DEL ISOVALOR EN FUNCIÓN DE T
vec4 GetTEquation(in float d000, in float d001, in float d101, in float d100, in float d010,
	in float d011, in float d111, in float d110, in float isoValue,
	in vec3 pEntry, in vec3 r)
{
	// Ecuación tras desarrollar la interpolación trilineal de d y sustituyendo x, y, z
	// d000 = a		d100 = b	d010 = c	d110 = e
	// d001 = f		d101 = g	d011 = h	d111 = i
	// x = x0 + rx * t	|		  x0 = j	  y0 = k
	// y = y0 + ry * t	|		  z0 = l	  rx = m
	// z = z0 + rz * t	|		  ry = n	  rz = o

	// Cambio de nomenclatura
	float a = d000; float b = d100; float c = d010; float e = d110;
	float f = d001; float g = d101; float h = d011; float i = d111;
	float j = pEntry.x;
	float k = pEntry.y;
	float l = pEntry.z;
	float m = r.x;
	float n = r.y;
	float o = r.z;

	// Variables auxiliares
	float mno = m*n*o;
	float jno = j*n*o;
	float kmo = k*m*o;
	float klm = k*l*m;
	float lmn = l*m*n;
	float jko = j*k*o;
	float jln = j*l*n;
	float jkl = j*k*l;
	float mn = m*n;
	float mo = m*o;
	float no = n*o;
	float jn = j*n;
	float jo = j*o;
	float jk = j*k;
	float jl = j*l;
	float kl = k*l;
	float km = k*m;
	float ko = k*o;
	float lm = l*m;
	float ln = l*n;

	// Construcción de la ecuación
	vec4 equation = vec4(0.0);	// .x = t^3		.y = t^2	.z = t		.w = 1

	// t^3
	equation.x = (-a*mno) + (b*mno) + (c*mno) + (-e*mno) + (f*mno) + (-g*mno) + (-h*mno) + (i*mno);

	// t^2
	equation.y += a * (-jno - kmo - lmn + mn + mo + no);
	equation.y += b * (jno + kmo + lmn - mn - mo);
	equation.y += c * (jno + kmo + lmn - mn - no);
	equation.y += e * (-jno - kmo - lmn + mn);
	equation.y += f * (jno + kmo + lmn - mo - no);
	equation.y += g * (-jno - kmo - lmn + mo);
	equation.y += h * (-jno - kmo - lmn + no);
	equation.y += i * (jno + kmo + lmn);

	// t
	equation.z += a * (-jko - jln + jn + jo - klm + km + ko + lm + ln - m - n - o);
	equation.z += b * (jko + jln - jn - jo + klm - km - lm + m);
	equation.z += c * (jko + jln - jn + klm - km - ko - ln + n);
	equation.z += e * (-jko - jln + jn - klm + km);
	equation.z += f * (jko + jln - jo + klm - ko - lm - ln + o);
	equation.z += g * (-jko - jln + jo - klm + lm);
	equation.z += h * (-jko - jln - klm + ko + ln);
	equation.z += i * (jko + jln + klm);

	// 1
	equation.w += a * (-jkl + jk + jl - j + kl - k - l + 1);
	equation.w += b * (jkl - jk - jl + j);
	equation.w += c * (jkl - jk - kl + k);
	equation.w += e * (-jkl + jk);
	equation.w += f * (jkl - jl - kl + l);
	equation.w += g * (-jkl + jl);
	equation.w += h * (-jkl + kl);
	equation.w += i * (jkl);

	// Buscamos X*t^3 + Y*t^2 + Z*t + W = isoValue		--->	X*t^3 + Y*t^2 + Z*t + W - isoValue = 0
	equation.w -= isoValue;

	return equation;
}

void main()
{
	mat4 vmMatrix = view*model;
	mat4 vmMatrixInverse = inverse(vmMatrix);
	mat4 vmNMatrix = transpose(inverse(vmMatrix));

	vec2 screenCoords = gl_FragCoord.st / screenSize;
	vec3 exitPoint = texture(exitPointTex, screenCoords).xyz;

	if (entryPoint.z <= exitPoint.z)
	{
		outColor = backgroundColor;
	}
	else{
		//Direccion del rayo. Calculamos el rayo un cubo en el que el 
		volumeDim -= vec3(1.0);
		vec3 volumeDimInv = 1.0 / volumeDim;
		vec3 dir = exitPoint - entryPoint;
		//vector escalado en coordenadas de la textura
		dir = volumeDim*(vmMatrixInverse * vec4(dir, 0)).xyz;
		vec3 pto = volumeDim*(vmMatrixInverse * vec4(entryPoint, 1)).xyz;

		//Direcciones de la recta

		vec3 incDir = sign(dir);

		//variables sobre las que iteramos (t de la recta)
		float t = 0; //la t del entrypoint es 0 (la recta se calcula usando este punto)
		float tfin = 1; //la ultima t será 1 ya que el vector está calculado como exitpoint-entrypoint
		int i = 0;
		vec3 voxelID;
		vec3 currentPoint = pto; //primer punto de corte
		vec4 colorAcum = vec4(0);
		bool end = false;

		while (t <= tfin && i < MAX_ITERATIONS && !end){
			i++;
			//obtenemos el voxelID del punto
			voxelID = getVoxelID(dir, currentPoint);
			
			//calculamos siguiente entrypoint y su t correspondiente
			//primero calculamos que planos se cortan
			vec3 faces = incDir * 0.5 + voxelID;
			vec3 tcortes = (faces - pto) / dir;

			//obtenemos la menor t válida
			bool found = false;
			float tnext;
			if (incDir.x != 0.0)
			{
				tnext = tcortes.x;
				found = true;
			}
			if (incDir.y != 0.0)
			{
				if (found)
					tnext = min(tnext, tcortes.y);
				else
				{
					tnext = tcortes.y;
					found = true;
				}
			}
			if (incDir.z != 0.0)
			{
				if (found)
					tnext = min(tnext, tcortes.z);
				else
				{
					tnext = tcortes.z;
					found = true;
				}
			}
			//float tnext = min(tcortes.x, min(tcortes.y, tcortes.z));
			//calculamos el entrypoint correspondiente
			vec3 pointnext = pto + dir * tnext;

			//AQUI SE OBTIENEN LAS DENSIDADES DE LOS VERTICES
			float d0, d1, d2, d3, d4, d5, d6, d7;
			computeDensity(voxelID, d0, d1, d2, d3, d4, d5, d6, d7);
			if (checkIsoValue(d0, d1, d2, d3, d4, d5, d6, d7, isoValue)){

				float A = d0;
				float X = d3 - d0;
				float Y = d4 - d0;
				float Z = d1 - d0;
				float XY = d7 + d0 - d3 - d4;
				float ZX = d2 + d0 - d3 - d1;
				float YZ = d5 + d0 - d4 - d1;
				float XYZ = d6 + d3 + d4 + d1 - d0 - d5 - d2 - d7;

				vec3 ptotrans = currentPoint - (voxelID - vec3(0.5)); //translacion del voxel al origen
				vec3 pointNextTrans = pointnext - (voxelID - vec3(0.5)); // ***
				vec3 currentDir = pointNextTrans - ptotrans;
				vec4 ec;
				/*
				ec.x = XYZ*dir.x*dir.y*dir.z;
				ec.y = XYZ*(ptotrans.z*dir.y*dir.x + dir.z*ptotrans.y*dir.x + dir.z*dir.y*ptotrans.x) + XY*dir.y*dir.x + ZX*dir.z*dir.x + YZ*dir.z*dir.y;
				ec.z = XY*(ptotrans.y*dir.x + dir.y*ptotrans.x) + YZ*(ptotrans.z*dir.y + dir.z*ptotrans.y) + XYZ*(dir.z*ptotrans.y*ptotrans.x + ptotrans.z*ptotrans.y*dir.x + ptotrans.z*dir.y*ptotrans.x) + ZX*(ptotrans.z*dir.x + dir.z*ptotrans.x) + Y*dir.y + Z*dir.z + X*dir.x;
				ec.w = A - isoValue + XY*ptotrans.y*ptotrans.x + ZX*ptotrans.z*ptotrans.x + YZ*ptotrans.z*ptotrans.y + XYZ*ptotrans.z*ptotrans.y*ptotrans.x + X*ptotrans.x + Z*ptotrans.z + Y*ptotrans.y;
				*/
				//******Irene -> GetTEquation OK, funciona igual con esta línea que comentándola
				ec = GetTEquation(d0, d1, d2, d3, d4, d5, d6, d7, isoValue, ptotrans, currentDir);

				//colorAcum += 0.1*vec4(ec.w);
				//Ecuación cúbica en funcion de T At^3 + Bt^2 + Ct + D -isoValue = 0 
				//Para encontrar los puntos de interseccion con el isovalor
				vec3 sol;
				int nsol;
				//Hay que controlar el caso en el que A=0 (Ecuacion de segundo grado) hecho en el solver
				nsol = solveCubic(ec.x, ec.y, ec.z, ec.w, sol);
				

				
				
				
				//colorAcum = vec4(voxelID/volumeDim, 1);
				//break;

				//Esto saca los t de corte con la isosuperficie
				//Hay que comprobar cuales caen dentro del voxel y dentro de estos, cual es el menor
				// tmin = t, tmax = tnext

				//comprobamos las soluciones obtenidas
				float s;
				bool validsolution;
				validsolution = checkSolution2(sol, t, tnext, nsol, s);



				

				if (validsolution)
				{
					/*
					if (nsol == 1)
						colorAcum = vec4(1, 0, 0, 1);
					else if (nsol == 2)
						colorAcum = vec4(0, 1, 0, 1);
					else if (nsol == 3)
						colorAcum = vec4(0, 0, 1, 1);
					else
						colorAcum = vec4(1, 1, 1, 1);
					break;
					*/
					
					//colorAcum = 2 * vec4(isoValue, 0, 0, 0);
					//calculo del punto
					vec3 ptosol = ptotrans + currentDir*s; //ptotrans + dir*s; //***
					

					//Calculo de normales

					vec3 norm = vec3(0.0);
					norm.x = X + XY * ptosol.y + ZX * ptosol.z + XYZ * ptosol.y * ptosol.z;
					norm.y = Y + XY * ptosol.x + YZ * ptosol.z + XYZ * ptosol.x * ptosol.z;
					norm.z = Z + ZX * ptosol.x + YZ * ptosol.y + XYZ * ptosol.x * ptosol.y;
					norm = norm * -1;


					//deshacemos la translacion de voxelID
					ptosol = ptosol + voxelID - vec3(0.5);
					ptosol = ptosol / volumeDim;
					vec3 texSample = ptosol;    //***
					ptosol = (vmMatrix * vec4(ptosol, 1)).xyz;

					

					norm = norm * volumeDim; // inversa del escalado que pasa a coord cámara M^{-t}*S^{-t}
					norm = (vmNMatrix * vec4(norm, 0)).xyz; //coordenadas de la camara
					norm = normalize(norm);


					colorAcum = shade(norm, ptosol);
					end = true;
				}
			}
			//actualizamos las variables del while
			t = tnext;
			currentPoint = pointnext;
		}
		outColor = colorAcum;
	}
}
