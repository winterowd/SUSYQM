/* 

Metropolis Monte Carlo for SUSY QM with quartic superpotential

This program generates the configurations of the auxiliary field, from which the 1 and 2 point functions are computed 

v1.0: we compute, ``on the fly'', the 1--point function, <F>
v1.1: we compute, ``on the fly'', the 2--point function, <F[l]F[l+n]>
      we compute, ``on the fly'',  the Binder cumulant, <F^4>-3<F^2>^2
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define sqr(x) ((x)*(x))
#define sign(x) ((x)>=0.0 ? 1 : -1)


double scale, m2latt, lambdaR; // scale factor,  squared lattice mass, ren.coup.

int N; // lattice size
int MINIT;
int M; // configurations
int P; // packets
int sign; // sign of the mass term: +1 or -1



double Wprime(double phi){// dW/dPhi_n
  double dWdPhin;

  dWdPhin = phi*(1.0L + (sqr(phi)/6.0L));

  return(dWdPhin);
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




void input(double *Phi,double dPhi){// initializes the configuration at random in the interval [-1,1]
  int n;
  double rand;
  for(n=0;n<N;n++){
    rand = drand48();
    //printf("rand: %e\n", rand);
    //Phi[n] = dPhi*(1.0-2.0*drand48()); /* hot start */
    Phi[n] = dPhi*(1.0-2.0*rand);
    //printf("Phi[%d]: %e\n", n, Phi[n]);
    //Phi[n]=0.0; // cold start
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


main(int argc, char *argv[]){
  N = atoi(argv[1]);
  lambdaR = atof(argv[2]);
  double dPhi = atof(argv[3]);
  //  m2latt = 1.0/((double)N*lambdaR*scale);
  m2latt = atof(argv[4]);
  MINIT = atoi(argv[5]);
  M = atoi(argv[6]);
  P = atoi(argv[7]);
  printf("N: %d lambdaR: %e dist: %d dPhi: %e m2latt: %e warms: %d trajecs: %d meas: %d\n", N, lambdaR, dist, dPhi, m2latt, MINIT, M, P);
  scale = 1.0;
  sign = 1;

 

  FILE *acc_rat, *vevf1, *vevf2, *vevphi, *vevphi2, *vevf4, *vevphi4;
  acc_rat = fopen("acc_ratio.out","w");
  vevf1   = fopen("vevf1.out", "w");
  vevf2   = fopen("vevf2.out", "w");
  vevphi = fopen("vevphi.out", "w");
  vevphi2 = fopen("vevphi2.out", "w");
  vevf4 = fopen("vevf4.out", "w");
  vevphi4  = fopen("vevphi4.out", "w");

  double Phi[N+1], F[N+1];
  double deltaPhi;

  double bicF, bicPhi, vbicF, vbicPhi;
  int samples;

  int seed=0;
  int accept; double accept_ratio;
  int m, n, np1, nm1;

  (void)srand48(seed);
  //printf("first random %e\n", drand48());

  input(Phi,dPhi); //initialize the scalar

  

    for(m=1;m<=MINIT;m++){

      accept = 0;

    for(n=0;n<N;n++){
      deltaPhi = dPhi*(1.0L-2.0L*drand48());//choose deltaPhi uniformly
      //printf("%d deltaPhi: %e\n", n, deltaPhi);
      //printf("deltaS: %e\n",DeltaS(Phi, deltaPhi, n));
      if(exp(-DeltaS(Phi,deltaPhi,n))>drand48()){
	//printf("ACCEPT!\n");
	Phi[n]=Phi[n]+deltaPhi; accept++;
      }
      //else 
      //printf("REJECT!\n");
    }


    }

    
samples=0; bicF=0.0; bicPhi=0.0; vbicF=0.0; vbicPhi=0.0;

    for(m=1;m<=M;m++){

      accept = 0;

    for(n=0;n<N;n++){
            deltaPhi = dPhi*(1.0L-2.0L*drand48());//choose deltaPhi uniformly about zero
            if(exp(-DeltaS(Phi,deltaPhi,n))>drand48()){
	Phi[n]=Phi[n]+deltaPhi; accept++;
      }
    }

    //compute the auxiliary field

    for(n=0;n<N;n++){
      np1 = (n+1)%N; nm1 = (N+n-1)%N;
      F[n] = 0.5*(Phi[np1]-Phi[nm1])+ m2latt*Wprime(Phi[n]);
    }

    accept_ratio = (double)accept/(double)N;

    

    // record the MC time series 
    if(m%P==0){
    fprintf(acc_rat,"%d\t%g\n",m,accept_ratio);
    fprintf(vevf1,"%d\t%g\n",m,vevF1(F));
    fprintf(vevf4, "%d\t%g\n", m, vevF4(F));
    fprintf(vevphi4, "%d\t%g\n", m, vevF4(Phi));
    fprintf(vevphi, "%d\t%g\n", m, vevF1(Phi));
    for(n=0; n<N/2; n++) {
      fprintf(vevphi2, "%d\t%d\t%g\n", m, n, new_vevF2(Phi, n));
      fprintf(vevf2,"%d\t%d\t%g\n",m, n, new_vevF2(F, n));
    }
    }
   

   
  }

    

    fclose(acc_rat);
    fclose(vevf1);
    fclose(vevf2);
    fclose(vevf4);
    fclose(vevphi4);
    fclose(vevphi);
    fclose(vevphi2);
    // printf("%g\t%g\n", bicF, vbicF-sqr(bicF));
    //printf("%g\t%g\n", bicPhi, vbicPhi-sqr(bicPhi));
}
