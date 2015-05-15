/* 

HMC for Double Well Potential

*TEST for HMC in QM

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

double f2, lambdaR, mass, a; // well minimum,  ren.coup., mass, lattice spacing

double *Phi, *Phi_old, *H;

int N; // lattice size
int warms;
int hits; //number of hits per site
int trajecs; // number of trajectories
int meas; // when to measure
int sign; // sign of the mass term: +1 or -1
int accepts=0; //number of accepts
gsl_rng *r; //random number generator
char start[4];

//debugging variables
double S1, S2;
double S1_old, S2_old;
double dPhi;

double action(double *Phi) {

  int n, np1;
  double Smom, Snonloc, Sloc;
 
  Smom = Snonloc = Sloc = 0.0;
  for(n=0; n<N; n++) {
    np1 = (n+1)%N;
    Snonloc += mass*(Phi[n]*Phi[n] - Phi[n]*Phi[np1]); 
    Sloc += lambdaR*( Phi[n]*Phi[n] - f2)*( Phi[n]*Phi[n] - f2);
  }
  S1 = Snonloc; S2 = Sloc;
  printf("S1: %e S2: %e\n", S1, S2);
  return(Snonloc/a + a*Sloc);
}//action()

void update_Phi() {

  int n, i;
  double delta;
  for(n=0; n<N; n++) {
    for(i=0; i<hits; i++) {
      delta = dPhi*(1.0-2.0*gsl_rng_uniform(r));
      Phi[n] += delta;
    //printf("update Phi[%d] = %e\n", n, Phi[n]);
    }
  }

}//update_Phi()

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

  double action_old = action(Phi);
  S1_old = S1; S2_old = S2;
  //printf("S1_old: %e S2_old: %e S3_old: %e S_mom_old: %e\n", S1_old, S2_old, S3_old, S_mom_old);
  //printf("S1: %e S2: %e S3: %e S_mom: %e\n", S1, S2, S3, S_mom);
  printf("initial action = %e\n", action_old);
  
  
  update_Phi();

  //metropolis step
  printf("final action = %e\n", action(Phi));
  
  double deltaS = action(Phi) - action_old;
  double rand = gsl_rng_uniform(r);
  if( exp((double)-deltaS) > rand) {
    accepts++;
    printf("ACCEPT: deltaS: %e DS1: %e DS2: %e DSMOM: %e\n", deltaS, (S1-S1_old), 
	   (S2-S2_old), (S_mom-S_mom_old));
    for(n=0; n<N; n++)
      Phi_old[n] = Phi[n];
  }
  else {
    printf("REJECT: deltaS: %e DS1: %e DS2: %e DSMOM: %e\n", deltaS, (S1-S1_old), 
	   (S2-S2_old), (S_mom-S_mom_old));
    for(n=0; n<N; n++)
      Phi[n] = Phi_old[n];
  }

}//update()


int main(int argc, char *argv[]){
  N = atoi(argv[1]);
  lambdaR = atof(argv[2]);
  f2 = atof(argv[3]);
  mass = atof(argv[4]);
  a = atof(argv[5]);
  warms = atoi(argv[6]);
  trajecs = atoi(argv[7]);
  meas = atoi(argv[8]);
  step_size = atof(argv[9]);
  num_steps = atoi(argv[10]);
  int seed = atoi(argv[11]);
  strcpy(start, argv[12]);
  sigma_p = atof(argv[13]);
  sigma_p2 = sigma_p*sigma_p;
  lambda = atof(argv[14]);
  printf("N: %d lambdaR: %e f2: %e mass: %e a: %e warms: %d trajecs: %d meas: %d step_size: %e num_steps %d seed: %d sigma_p: %e lambda: %e\n", N, lambdaR, f2, mass, a, warms, trajecs, meas, step_size, num_steps, seed, sigma_p, lambda);
  
  FILE *vevphi, *vevphi2, *vevphi4, *vevphi6, *vevphi8;
  vevphi = fopen("vevphi.out", "w");
  vevphi2 = fopen("vevphi2.out", "w");
  vevphi4  = fopen("vevphi4.out", "w");
  vevphi6  = fopen("vevphi6.out", "w");
  vevphi8  = fopen("vevphi8.out", "w");
   
  Phi = malloc((N+1)*sizeof(double));
  Phi_old = malloc((N+1)*sizeof(double));
  H = malloc((N+1)*sizeof(double));
  
  int m, n;

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
    
  for(m=1;m<=trajecs;m++){
    
    update();
    
    
    // record the MC time series 
    if(m%meas==0){
      
      for(n=0; n<N/2; n++) {
	fprintf(vevphi2,"%d\t%d\t%g\n",m,n,new_vevF2(Phi,n));
      }
      for(n=0; n<N; n++) {
	printf("DEBUG_phi: %d %d %e\n", m, n, Phi[n]);
      }
      fprintf(vevphi4, "%d\t%g\n", m, vevF4(Phi));
      fprintf(vevphi, "%d\t%g\n", m, vevF1(Phi));
      fprintf(vevphi6, "%d\t%g\n", m, vevF6(Phi));
      fprintf(vevphi8, "%d\t%g\n", m, vevF8(Phi));
      
    } 
    
  } 
  printf("TRAJECTORIES FINISHED!\n");
  printf("ACCEPTS: %d\n", accepts);

  fclose(vevphi4);
  fclose(vevphi);
  fclose(vevphi2);
  fclose(vevphi6);
  fclose(vevphi8);
  
  free(Phi);
  free(Phi_old);
  free(H);
  gsl_rng_free(r);
  // printf("%g\t%g\n", bicF, vbicF-sqr(bicF));
  //printf("%g\t%g\n", bicPhi, vbicPhi-sqr(bicPhi));
}
