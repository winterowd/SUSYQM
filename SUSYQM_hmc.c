/* 

Metropolis Monte Carlo for SUSY QM with quartic superpotential

This program generates the configurations of the auxiliary field, from which the 1 and 2 point functions are computed 

v1.0: we compute, ``on the fly'', the 1--point function, <F>
v1.1: we compute, ``on the fly'', the 2--point function, <F[l]F[l+n]>
      we compute, ``on the fly'',  the Binder cumulant, <F^4>-3<F^2>^2
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#define sqr(x) ((x)*(x))
#define sign(x) ((x)>=0.0 ? 1 : -1)
#define abs(x) (sqrt(x*x))

double scale, m2latt, lambdaR; // scale factor,  squared lattice mass, ren.coup.

double *Phi, *Phi_old, *F, *H, *Wpp;

int N; // lattice size
int warms;
int trajecs; // number of trajectories
int meas; // when to measure
double step_size; //step size for integrating EOM
int num_steps; //number of steps in a trajectory
int sign; // sign of the mass term: +1 or -1
int accepts=0; //number of accepts
gsl_rng *r; //random number generator
char start[4];
double sigma_p, sigma_p2;
double lambda;

//debugging variables
double S1, S2, S3, S_mom;
double S1_old, S2_old, S3_old, S_mom_old;
double dPhi;

double Wprime(double phi){// dW/dPhi_n
  double dWdPhin;

  dWdPhin = phi*(1.0L + (sqr(phi)/6.0L));

  return(dWdPhin);
}

double Wprimeprime(double phi) {// dW/dPhi_n
  double d2Wd2Phin;

  d2Wd2Phin = (1.0L + 0.5*sqr(phi));

  return(d2Wd2Phin);
}

double detWprimeprime(double *Phi) { // |detW''|
  int n;
  double det=1.;
  for(n=0; n<N; n++) { 
    //printf("%d %e %e\n", n, Phi[n], Wprimeprime(Phi[n]));
    det = det*Wprimeprime(Phi[n]);
  }
  //printf("det: %e\n", det);
  return(det);
}


double DeltaS(double *Phi, double deltaPhi, int n){
// computes the change in the action, when we propose Phi[n]->Phi[n] + deltaPhi
  double deltaS, deltaSnonloc, deltaSloc, deltaSpoly, Phinew=Phi[n]+deltaPhi;
  double Phiold = Phi[n];

  int np1, nm1; 

  // impose periodic boundary conditions on the scalar

  np1 = (n+1)%N;
  nm1 = (n+N-1)%N;

  deltaSnonloc = -deltaPhi*(Phi[np1]+Phi[nm1]) + deltaPhi*(2.0*Phi[n]+deltaPhi);
  deltaSloc = 0.5*m2latt*( (Wprime(Phinew)-Wprime(Phiold) )*( Wprime(Phinew)+Wprime(Phiold) ));
  deltaSpoly = log((2.0+sqr(Phiold))/(2.0+sqr(Phinew)) );
  deltaS = (1.0/(m2latt*lambdaR))*(deltaSnonloc)+(deltaSloc/lambdaR) + deltaSpoly;
  return(deltaS);
}

double action(double *Phi) {

  int n, np1, nm1;
  double Smom, Snonloc, Sloc, Spoly;
 
  Smom = Snonloc = Sloc = Spoly = 0.0;
  for(n=0; n<N; n++) {
    np1 = (n+1)%N;
    Snonloc += (1.0/(m2latt*lambdaR))*(Phi[n]*Phi[n] - Phi[n]*Phi[np1]); 
    Sloc += 0.5*(m2latt/lambdaR)*( Wprime(Phi[n])*Wprime(Phi[n]) );
    //Spoly += -scale*lambdaR*m2latt*log((1.0/(scale*lambdaR))*(1.0 + 0.5*Phi[n]*Phi[n]));
#ifndef FACTOR
    Spoly += -log(Wprimeprime(Phi[n]));
#endif
    Smom += 0.5*H[n]*H[n]/sigma_p2;
  }
  S1 = Snonloc; S2 = Sloc; S3 = Spoly; S_mom = Smom;
  printf("S1: %e S2: %e S3: %e Smom: %e\n", S1, S2, S3, S_mom);
  return(Smom + Snonloc + Sloc + Spoly);
}//action()

void update_Phi(double eps) {

  int n;
  for(n=0; n<N; n++) {
    Phi[n] += eps*H[n]/sigma_p2;
    //printf("update Phi[%d] = %e\n", n, Phi[n]);
  }

}//update_Phi()

double force(int n) {

  int np1, nm1;
  double force_nonloc, force_loc, force_poly;

  np1 = (n+1)%N;
  nm1 = (n+N-1)%N;
  force_nonloc = force_loc = force_poly = 0.0;
  force_nonloc = (1.0/(m2latt*lambdaR))*(Phi[np1] + Phi[nm1] - 2.0*Phi[n]);
  force_loc = -(m2latt/lambdaR)*Wprime(Phi[n])*Wprimeprime(Phi[n]);
  //force_poly = scale*lambdaR*m2latt*Phi[n]/(1.0 + 0.5*Phi[n]*Phi[n]);
#ifndef FACTOR
  force_poly = Phi[n]/(Wprimeprime(Phi[n]));
#endif
  //printf("F1: %e F2: %e F3: %e\n", c1*force_nonloc, c2*force_loc, c3*force_poly);
  return(force_nonloc + force_loc + force_poly);
}//force

void update_momenta(double eps) {

  int n;
  for(n=0; n<N; n++) {
    H[n] += force(n)*eps;
    //printf("update H[%d] = %e\n", n, H[n]);
  }

}//update_momenta()

void input(double *Phi){// initializes the configuration at random in the interval [-1,1]
  int n;
  double rand;
  if( strncmp("COLD", start, 4) == 0 ) { //cold start
    printf("COLD START!\n");
    for(n=0;n<N;n++){
      Phi_old[n] = Phi[n]=0.0; // cold start
    }
  }
  else { //hot start
    printf("HOT START!\n");
    for(n=0;n<N;n++){
      Phi_old[n] = Phi[n] = 0.00001*(1.0-2.0*gsl_rng_uniform(r)); /* hot start */
      /*rand =drand48();
      printf("rand: %e\n", rand);
      Phi[n] = dPhi*(1.0-2.0*rand);
      printf("Phi[%d]: %e\n", n, Phi[n]);
      Phi_old[n] = Phi[n]; */
      //printf("Phi[%d]: %e\n", n, Phi[n]);
    }
  }
  
}

double vevF1(double *F){//compute <<F[n]>> over the lattice thus obtaining a less noisy estimator for <F>
  
  int n; double vevf1=0.0;
  for(n=0;n<N;n++){
    vevf1 = vevf1 + F[n];
  }
  return(vevf1/(double)N);
}

double potential(double *Phi, int n) {
  return(Wprime(Phi[n])*Wprime(Phi[n]));
}

double new_vevF2(double *F, int n){// smeared estimator for <<F[l]F[l+n]>>
  int l; double vf2 = 0.0;
  for(l=0;l<N;l++){
    vf2 = vf2 + F[l]*F[(l+n)%N];
  }
  vf2 = vf2/(double)N;
  return(vf2);
}



double vevF4(double *F){// smeared estimator for << F[n]^4>>
  int n;
  double vf4 = 0.0;
  for(n=0;n<N;n++){
    vf4 = vf4 + sqr(sqr(F[n]));
  }
  return(vf4/(double)N);
}


double vevF6(double *F){// smeared estimator for << F[n]^6>>
  int n;
  double vf6 = 0.0;
  for(n=0;n<N;n++){
    vf6 = vf6 + sqr(sqr(sqr(F[n])));
  }
  return(vf6/(double)N);
}


double vevF8(double *F){// smeared estimator for << F[n]^8>>
  int n;
  double vf8 = 0.0;
  for(n=0;n<N;n++){
    vf8 = vf8 + sqr(sqr(sqr(sqr(F[n]))));
  }
  return(vf8/(double)N);
}


double SD3a(double *F) { //estimator for <<\Phi[n]^3\Phi[n+1]>>

  int l; double sd3 = 0.0;
  for(l=0;l<N;l++){
    sd3 = sd3 + sqr(F[l])*F[l]*F[(l+1)%N];
  }
  sd3 = sd3/(double)N;
  return(sd3);
}


double SD3b(double *F) { //estimator for <<\Phi[n]^3\Phi[n-1]>>

  int l; double sd3 = 0.0;
  for(l=0;l<N;l++){
    sd3 = sd3 + sqr(F[l])*F[l]*F[(l-1)%N];
  }
  sd3 = sd3/(double)N;
  return(sd3);
}

double SD3c(double *f1, double *f2) { //estimator for <<\Phi[n]^4 Wpp[n]>>

  int l; double sd3 = 0.0;
  for(l=0;l<N;l++){
    sd3 = sd3 + sqr(sqr(f1[l]))*f2[l];
  }
  sd3 = sd3/(double)N;
  return(sd3);

}

double SD1a(double *f1, double *f2) { //estimator for <<\Phi[n]^4 Wpp[n]>>

  int l; double sd1 = 0.0;
  for(l=0;l<N;l++){
    sd1 = sd1 + sqr(f1[l])*f2[l];
  }
  sd1 = sd1/(double)N;
  return(sd1);

}

/* double BiC(double *X){// Binder cumulant: <F^4>-3<F^2>^2
  int n; double vevf2=0.0, vevf4=0.0;

  for(n=0;n<N;n++){
    vevf2 = vevf2 + sqr(X[n]);
    vevf4 = vevf4 + sqr(sqr(X[n]));
  }
  vevf2 = vevf2/(double)N; vevf4 = vevf4/(double)N;
  //double bic = vevf4 - 3.0*sqr(vevf2);
  double bic = vevf4/(3.0*sqr(vevf2));
  return(bic);

}

*/

void ranmom() {

  int n;
  for(n=0;n<N; n++) {
    H[n] = gsl_ran_gaussian(r, sigma_p);
    //H[n] = -H[n]; //test reversibility of HMC
    //printf("HMOM[%d]: %e\n", n, H[n]);
  }

}//init_mom()

void update() {

  int i, n;
  //double S1, S2;
  double deltaPhi;

  /*for(i=0; i<N; i++)
    H[n] = 0.;
  
  for(n=0; n<N; n++) {
    S1 = action(Phi);
    printf("S1: %e\n", S1);
    deltaPhi = dPhi*(1.0L-2.0L*drand48());
    printf("%d deltaPhi: %e\n", n, deltaPhi);
    Phi[n] = Phi[n] + deltaPhi;
    S2 = action(Phi);
    printf("S2: %e\n", S2);
    printf("deltaS: %e\n", S2-S1);
    if(exp(-(S2-S1)) >drand48()) {
      printf("ACCEPT!\n");
    }
    else {
      printf("REJECT!");
      Phi[n] = Phi_old[n];
    }
  }

  return;*/

  ranmom();
  double action_old = action(Phi);
  S1_old = S1; S2_old = S2; S3_old = S3; S_mom_old = S_mom;
  //printf("S1_old: %e S2_old: %e S3_old: %e S_mom_old: %e\n", S1_old, S2_old, S3_old, S_mom_old);
  //printf("S1: %e S2: %e S3: %e S_mom: %e\n", S1, S2, S3, S_mom);
  printf("initial action = %e\n", action_old);
  
  for(i=0; i<num_steps; i++) {

    /* replace 2nd order leap from with Omelyan
    update_Phi(step_size/2.);
    update_momenta(step_size);
    update_Phi(step_size/2.);*/

    update_Phi(lambda*step_size);
    update_momenta(step_size/2.);
    update_Phi((1.-2.*lambda)*step_size);
    update_momenta(step_size/2.);
    update_Phi(lambda*step_size);

  }
  //metropolis step
  printf("final action = %e\n", action(Phi));
  //additional sign flip 
  if(gsl_rng_uniform(r) > .5) {
    for(n=0; n<N; n++)
      Phi[n]=-Phi[n];
  }
  double deltaS = action(Phi) - action_old;
  double rand = gsl_rng_uniform(r);
  if( exp((double)-deltaS) > rand) {
    accepts++;
    printf("ACCEPT: deltaS: %e DS1: %e DS2: %e DS3: %e DSMOM: %e\n", deltaS, (S1-S1_old), 
	   (S2-S2_old), (S3-S3_old), (S_mom-S_mom_old));
    for(n=0; n<N; n++)
      Phi_old[n] = Phi[n];
  }
  else {
    printf("REJECT: deltaS: %e DS1: %e DS2: %e DS3: %e DSMOM: %e\n", deltaS, (S1-S1_old), 
	   (S2-S2_old), (S3-S3_old), (S_mom-S_mom_old));
    for(n=0; n<N; n++)
      Phi[n] = Phi_old[n];
  }

}//update()


int main(int argc, char *argv[]){
  N = atoi(argv[1]);
  lambdaR = atof(argv[2]);
  //  m2latt = 1.0/((double)N*lambdaR*scale);
  m2latt = atof(argv[3]);
  warms = atoi(argv[4]);
  trajecs = atoi(argv[5]);
  meas = atoi(argv[6]);
  step_size = atof(argv[7]);
  num_steps = atoi(argv[8]);
  int seed = atoi(argv[9]);
  strcpy(start, argv[10]);
  sigma_p = atof(argv[11]);
  sigma_p2 = sigma_p*sigma_p;
  lambda = atof(argv[12]);
  printf("N: %d lambdaR: %e m2latt: %e warms: %d trajecs: %d meas: %d step_size: %e num_steps %d seed: %d sigma_p: %e lambda: %e\n", N, lambdaR, m2latt, warms, trajecs, meas, step_size, num_steps, seed, sigma_p, lambda);
  scale = 1.0;
  sign = 1;

#ifdef FACTOR
  printf("FACTOR defined!");
  FILE *vevDetWpp, *vevphiDetWpp, *vevphi2DetWpp, *vevfDetWpp, *vevf2DetWpp;
  vevDetWpp = fopen("vevDetWpp.out", "w");
  vevphiDetWpp = fopen("vevphiDetWpp.out", "w");
  vevphi2DetWpp = fopen("vevphi2DetWpp.out", "w");
  vevfDetWpp = fopen("vevfDetWpp.out", "w");
  vevf2DetWpp = fopen("vevf2DetWpp.out", "w");
#endif
  
  FILE *vevf1, *vevf2, *vevphi, *vevphi2, *vevf4, *vevphi4, *vevphi6, *vevphi8;
  FILE *vevWpp1, *vevWpp2, *vevWpp4, *phiSD3a, *phiSD3b, *phiSD3c, *phiSD1a;
  vevf1   = fopen("vevf1.out", "w");
  vevf2   = fopen("vevf2.out", "w");
  vevphi = fopen("vevphi.out", "w");
  vevphi2 = fopen("vevphi2.out", "w");
  vevf4 = fopen("vevf4.out", "w");
  vevphi4  = fopen("vevphi4.out", "w");
  vevphi6  = fopen("vevphi6.out", "w");
  vevphi8  = fopen("vevphi8.out", "w");
  phiSD3a = fopen("phiSD3a.out", "w");
  phiSD3b = fopen("phiSD3b.out", "w");
  phiSD3c = fopen("phiSD3c.out", "w");
  phiSD1a = fopen("phiSD1a.out", "w");
  vevWpp1 = fopen("vevWpp1.out", "w");
  vevWpp2 = fopen("vevWpp2.out", "w");
  vevWpp4 = fopen("vevWpp4.out", "w");
 
  Phi = malloc((N+1)*sizeof(double));
  Phi_old = malloc((N+1)*sizeof(double));
  F = malloc((N+1)*sizeof(double));
  H = malloc((N+1)*sizeof(double));
  Wpp = malloc((N+1)*sizeof(double));
  
  double bicF, bicPhi, vbicF, vbicPhi;
  int samples;

  int m, n, np1, nm1;

  r = gsl_rng_alloc(gsl_rng_ranlux);
  gsl_rng_set(r,seed);
  //(void)srand48(0);
  //printf("first random %e\n", drand48());
  input(Phi); //initialize the scalar


  for(m=1;m<=warms;m++){
    
    update();
    
  }
  //return;
  printf("WARMUPS FINISHED!\n");
  accepts=0;
  samples=0; bicF=0.0; bicPhi=0.0; vbicF=0.0; vbicPhi=0.0;
  
  for(m=1;m<=trajecs;m++){
    
    update();
    
    
    // record the MC time series 
    if(m%meas==0){
      
      //compute the auxiliary field and exp(-log(Wprimeprime))
      for(n=0;n<N;n++){
	np1 = (n+1)%N; nm1 = (N+n-1)%N;
	F[n] = 0.5*(Phi[np1]-Phi[nm1])+ m2latt*Wprime(Phi[n]);
	Wpp[n] = exp(-log(Wprimeprime(Phi[n])));
      }
      
      
      fprintf(vevf1,"%d\t%g\n",m,vevF1(F));
      for(n=0; n<N/2; n++) {
	fprintf(vevf2,"%d\t%d\t%g\n",m,n,new_vevF2(F,n));
	fprintf(vevphi2,"%d\t%d\t%g\n",m,n,new_vevF2(Phi,n));
	fprintf(vevWpp2,"%d\t%d\t%g\n",m,n,new_vevF2(Wpp,n));
#ifdef FACTOR
	fprintf(vevphi2DetWpp,"%d\t%d\t%g\n",m,n,detWprimeprime(Phi)*new_vevF2(Phi,n));
	fprintf(vevf2DetWpp,"%d\t%d\t%g\n",m,n,detWprimeprime(Phi)*new_vevF2(F,n));
#endif
      }
      for(n=0; n<N; n++) {
	//printf("DEBUG_potential: %d %d %e\n", m, n, potential(Phi,n));
	//printf("DEBUG_phi: %d %d %e\n", m, n, Phi[n]);
      }
      fprintf(vevf4, "%d\t%g\n", m, vevF4(F));
      fprintf(vevphi4, "%d\t%g\n", m, vevF4(Phi));
      fprintf(vevphi, "%d\t%g\n", m, vevF1(Phi));
      fprintf(vevphi6, "%d\t%g\n", m, vevF6(Phi));
      fprintf(vevphi8, "%d\t%g\n", m, vevF8(Phi));
      fprintf(phiSD3a, "%d\t%g\n", m, SD3a(Phi));
      fprintf(phiSD3b, "%d\t%g\n", m, SD3b(Phi));
      fprintf(phiSD3c, "%d\t%g\n", m, SD3c(Phi,Wpp));
      fprintf(phiSD1a, "%d\t%g\n", m, SD1a(Phi,Wpp));
      fprintf(vevWpp1, "%d\t%g\n", m, vevF1(Wpp));
      fprintf(vevWpp4, "%d\t%g\n", m, vevF4(Wpp));
#ifdef FACTOR
      fprintf(vevDetWpp, "%d\t%g\n", m, detWprimeprime(Phi));
      fprintf(vevphiDetWpp, "%d\t%g\n", m, detWprimeprime(Phi)*vevF1(Phi));
      fprintf(vevfDetWpp, "%d\t%g\n", m, detWprimeprime(Phi)*vevF1(F));
#endif
    } 
    
  } 
  printf("TRAJECTORIES FINISHED!\n");
  printf("ACCEPTS: %d\n", accepts);

  fclose(vevf1);
  fclose(vevf2);
  fclose(vevf4);
  fclose(vevphi4);
  fclose(vevphi);
  fclose(vevphi2);
  fclose(vevphi6);
  fclose(vevphi8);
  fclose(phiSD3a);
  fclose(phiSD3b);
  fclose(phiSD3c);
  fclose(phiSD1a);
  fclose(vevWpp1);
  fclose(vevWpp2);
  fclose(vevWpp4);
#ifdef FACTOR
  fclose(vevDetWpp);
  fclose(vevphiDetWpp);
  fclose(vevphi2DetWpp);
  fclose(vevfDetWpp);
  fclose(vevf2DetWpp);
#endif

  free(Phi);
  free(Phi_old);
  free(F);
  free(H);
  free(Wpp);
  gsl_rng_free(r);
  // printf("%g\t%g\n", bicF, vbicF-sqr(bicF));
  //printf("%g\t%g\n", bicPhi, vbicPhi-sqr(bicPhi));
}
