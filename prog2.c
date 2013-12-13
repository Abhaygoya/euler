#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double length = 10.;
double time = .4;
double alpha = 1.;
int N;
double dx, dt, currentTime;


struct vector {
	double mass;
	double mom;
	double ene;
};

double getx(int i){
	double x = length*(i+.5)/N;
	return x-.5*length;
}

//For prob 2, initial conditions are more complicated. First attempt using smooth function f(x)=e^-x^2
void initialize(struct vector * U, double * L){
	int i;
	for(i=0; i<N; i++){
		//Get x value, determine U and F based on x. All calculated from initial vars rho, p, v, and gamma
		double pos = getx(i);
		U[i].mass=1+alpha*exp(-32.*pow(pos,2));
		double v = 5.916*(pow(U[i].mass,.2)-1);
		U[i].mom = U[i].mass*v;
		double P = pow(U[i].mass,1.4);
		L[i] = P;
		U[i].ene = .5*U[i].mass*pow(v,2)+P/.4;


	}
}

void updateF(struct vector * U, struct vector * F){
	int i;
	for(i=0;i<N;i++){
		//double v = 5.916*(pow(U[i].mass,.2)-1);
		//double P = pow(U[i].mass,1.4);
		double v = U[i].mom/U[i].mass;
		double P = (.4)*(U[i].ene-.5*U[i].mass*v*v);
		F[i].mass = U[i].mass*v;
		F[i].mom = U[i].mom*v + P;
		F[i].ene = (U[i].ene+P)*v;
	}
}

void fixed_bcs(struct vector * U, struct vector * F){

	double pos = getx(0);
	U[0].mass=1+alpha*exp(-pow(pos,2));
	double v = 5.916*(pow(U[0].mass,.2)-1);
	U[0].mom = U[0].mass*v;
	double P = pow(U[0].mass,1.4);
	U[0].ene = .5*U[0].mass*pow(v,2)+P/.4;

	pos = getx(N-1);
	U[N-1].mass=1+alpha*exp(-pow(pos,2));
	v = 5.916*(pow(U[N-1].mass,.2)-1);
	U[N-1].mom = U[N-1].mass*v;
	P = pow(U[N-1].mass,1.4);
	U[N-1].ene = .5*U[N-1].mass*pow(v,2)+P/.4;

}

double max(double a, double b, double c){
	if(a>=b&&a>=c) return a;
	if(b>=c&&b>=a) return b;
	else return c;
}

//Calculate alphas at i, plus parameter determines alpha+ for 1 and alpha- for 0
double get_alpha(int i, struct vector * U, struct vector * F, int plus){
	double vL = U[i].mom/U[i].mass;
	double csL = sqrt(1.4*(F[i].mom-U[i].mom*vL)/U[i].mass);
	double vR = U[i+1].mom/U[i+1].mass;
	double csR = sqrt(1.4*(F[i+1].mom-U[i+1].mom*vL)/U[i+1].mass);
	double alpha;
	if(plus==1) alpha = max(0,(vL+csL),(vR+csR));
	else alpha = max(0,(-vL+csL),(-vR+csR));

	return alpha;
}

struct vector getFlux(int i, struct vector * U, struct vector * F){
	double alpha_plus = get_alpha(i,U,F,1);
	double alpha_minus = get_alpha(i,U,F,0);

	struct vector fhll;
	fhll.mass = (alpha_plus*F[i].mass+alpha_minus*F[i+1].mass-alpha_plus*alpha_minus*(U[i+1].mass-U[i].mass))/(alpha_plus+alpha_minus);
	fhll.mom = (alpha_plus*F[i].mom+alpha_minus*F[i+1].mom-alpha_plus*alpha_minus*(U[i+1].mom-U[i].mom))/(alpha_plus+alpha_minus);
	fhll.ene = (alpha_plus*F[i].ene+alpha_minus*F[i+1].ene-alpha_plus*alpha_minus*(U[i+1].ene-U[i].ene))/(alpha_plus+alpha_minus);

	return fhll;

}

void advanceTime(struct vector * U, struct vector * F){
	struct vector Fiph[N-1];
	int i;
	for(i=0;i<N-1;i++) Fiph[i] = getFlux(i,U,F);
	for(i=0;i<N-1;i++){
		U[i].mass -= dt*Fiph[i].mass/dx;
		U[i].mom -= dt*Fiph[i].mom/dx;
		U[i].ene -= dt*Fiph[i].ene/dx;
		U[i+1].mass += dt*Fiph[i].mass/dx;
		U[i+1].mom += dt*Fiph[i].mom/dx;
		U[i+1].ene += dt*Fiph[i].ene/dx;
	}

	fixed_bcs(U,F);
	updateF(U,F);
}


double calculateL1(struct vector * U, double * L){
	double l1error=0.;
	int i;
	for(i=0;i<N;i++) l1error += (pow(U[i].mass,1.4)-L[i])*dx;
	return fabs(l1error);
}


int main(){
	N = 4000;
	while(N<6000){
		dt = time/(N);
		dx = length/N;
		struct vector U[N];
		struct vector F[N];

		double L[N];
		initialize(U,L);
		updateF(U,F);



		currentTime=0.;
		while(currentTime<time){
			advanceTime(U,F);
			currentTime += dt;
			if(isnan(U[N-1].mass)>0){
				printf("%f\n",currentTime);
				break;
			}
		}

		//Get L1 error
		double error = calculateL1(U,L);


		FILE* f;
		char string[15];
		sprintf(string,"%f",time);
		strcat(string,"_seconds.txt");
		f = fopen(string,"w");

		int i;
		for(i=0;i<N;i++) fprintf(f,"%15.15f     %15.15f       %15.15f    %15.15f    %15.15f\n", getx(i),U[i].mass,U[i].mom,U[i].ene,error);
		fclose(f);
		N = 2*N;

	}
	return 0;
}
