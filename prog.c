#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double length = 1.;
double time = .25;
int N;
double dx, dt, currentTime;


struct vector {
	double mass;
	double mom;
	double ene;
};

double getx(int i){
	double x = length*(i+.5)/N;
	return x-.5;
}

void initialize(struct vector * U, struct vector * F){
	int i;
	for(i=0; i<N; i++){
		//Get x value, determine U and F based on x. All calculated from initial vars rho, p, v, and gamma
		double pos = getx(i);
		if(pos<0.){
			U[i].mass = 1.0;
			U[i].mom = 0.;
			U[i].ene = 1.0/.4;
			F[i].mass = 0.;
			F[i].mom = 1.0;
			F[i].ene = 0.;
		}
		else{
			U[i].mass = .125;
			U[i].mom = 0.;
			U[i].ene = .1/.4;
			F[i].mass = 0.;
			F[i].mom = .1;
			F[i].ene = 0.;
		}
		
	}
}

void updateF(struct vector * U, struct vector * F){
	int i;
	for(i=0;i<N;i++){
		double gam = 1.4;
		double rho = U[i].mass;
		double v = U[i].mom/rho;
		double eps = U[i].ene - .5*rho*v*v;
		double P = (gam-1.)*eps;
		F[i].mass = rho*v;
		//We know E(gamma-1)=P from E.O.S, and gamma=1.4 for adiabatic
		F[i].mom = rho*v*v + P;
		F[i].ene = (.5*rho*v*v + eps + P)*v;
	}
}

void fixed_bcs(struct vector * U, struct vector * F){
	U[0].mass = 1.0;
	U[0].mom = 0.;
	U[0].ene = 1.0/.4;
	F[0].mass = 0.;
	F[0].mom = 1.0;
	F[0].ene = 0.;
	
	U[N-1].mass = .125;
	U[N-1].mom = 0.;
	U[N-1].ene = .1/.4;
	F[N-1].mass = 0.;
	F[N-1].mom = .1;
	F[N-1].ene = 0.;
}

double max(double a, double b, double c){
	if(a>=b&&a>=c) return a;
	if(b>=c&&b>=a) return b;
	else return c;
}

//Calculate alphas at i, plus parameter determines alpha+ for 1 and alpha- for 0
double get_alpha(int i, struct vector * U, struct vector * F, int plus){
	double vL = U[i].mom/U[i].mass;
	double csL = sqrt(1.4*(.4*U[i].ene)/U[i].mass);
	double vR = U[i+1].mom/U[i+1].mass;
	double csR = sqrt(1.4*(.4*U[i+1].ene)/U[i+1].mass);
	double alpha;
	if(plus==1) alpha = max(0,(vL+csL),(vR+csR));
	else alpha = max(0,(-vL+csL),(-vR+csR));
	
	return alpha;	
}

struct vector getFlux(int i, struct vector * U, struct vector * F){
	double alpha_plus = get_alpha(i,U,F,1);
	double alpha_minus = get_alpha(i,U,F,0);
	//if(alpha_plus==1./0) printf("Oh no! %f   %d",currentTime,i);
	
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
	updateF(U,F);
	fixed_bcs(U,F);
}
//the array L holds the analytic solution data
double calculateL1(struct vector * U, struct vector * L){
	double l1error=0.;
	int i;
	for(i=0;i<N;i++) l1error += fabs(U[i].mass-L[10000*i/N].mass)*dx;
	return l1error;
}	
		


	
int main(){
	//Read analytic solution from file
	FILE* sol;
	sol = fopen("solution.txt","r");
	struct vector L[10000];
	int i;
	double useless1,useless2;
	for(i=0;i<10000;i++){
		fscanf(sol,"%lf %lf %lf %lf %lf",&useless1, &L[i].mass, &L[i].mom, &L[i].ene, &useless2);
	}
	//Loop through different N's, advancing to time t=.25, and output
	N = 5;
	while(N<15000){
		dt = time/(2*N);
		dx = length/N;
	
		struct vector U[N];
		struct vector F[N];
	
		initialize(U,F);
	
	
	
		currentTime=0.;
		while(currentTime<time){
			advanceTime(U,F);
			currentTime += dt;
		}
		
		
		
		double error = calculateL1(U,L);
		
		FILE* f;
		char string[15];
		sprintf(string,"%d",N);
		strcat(string,".txt");
		f = fopen(string,"w");
		
		for(i=0;i<N;i++){
			double v = U[i].mom/U[i].mass;
			double eps = U[i].ene - .5*U[i].mass*v*v;
			double P = eps*.4;
			fprintf(f,"%15.15f     %15.15f       %15.15f    %15.15f    %15.15f\n", getx(i),U[i].mass,v,P,error);
		}
		
		fclose(f);
		N = 2*N;
		
	}
	return 0;
}
