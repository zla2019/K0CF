#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define pi 3.14159265358979323844
#define EULER 0.5772156649015328606
double_complex ci(0.0,1.0);

/* Prototypes */
void bininit(void);
void binner(int inumden,double mom,double corr,double *p1,double *p2);
void printstatements(int onlyzerob);
void coulset(void);
void coulpar(double eta, double_complex *couly,double_complex *cfact1,
	     double_complex *cfact2,double_complex *chype,
	     double_complex *chype1, double_complex *chype2);
double_complex coul(double mom,double r,double zk,double eta,int kk);
double_complex chyper(double_complex a,double_complex b,double_complex cz,
		      int kk);
void partwaveinit(void);
void diffy(int ipart);
double_complex cw2_big_r(int l,double r,double eta);
double_complex cw2_small_r(int l,double r,double eta);
double dgamma(int mm);
void phasespaceinit(int *nb,int *np1,int *np2,double *bweight_num,
		    double *bweight1_den,double *bweight2_den,int *onlyzerob);
void get_b_num(int *ib1,int *ib2,double *bweight_num);
void get_b_den(int *ib1,int *ib2,double *bweight1_den,double *bweight2_den);
void get_pr(int ib1,int j1,double phi1,int ib2,int j2,double phi2,
	    double *p1,double *r1,int *iread1,
	    double *p2,double *r2,int *iread2);
void onepart_filter(int i,double *p,int *ifilter);
void twopart_filter(double *p1,double *p2,int *ifilter);
double get_mom(double *p1,double *p2);
double  corrcalc(double mom,double kdotr,double r);
double bwcalc(double qred,double r);
double triangle(double m0,double m1,double m2);
void get_com_quantities(double *p1,double *r1,double *p2,double *r2,
			double *qred,double *r,double *qdotr,double *mom,
			int *test);
void lorentz(double *u,double *q,double *qprime);
double_complex cgamma(double_complex c);
double get_gamow(double mom);
void momentum_smear(double *p1,double *p2,double *mom,double *qred);
void ranrot(double *p,double *r,double phi);
void store_mixed_pair(int n_mixed_pairs,double *p1,double *pp1,double *p2,
		      double *pp2,int iread1,int iread2,int *iiread1,
		      int *iiread2,double corr,double *ww);
void decompose(double *p1,double *p2,double vs,
	       int long_comoving_boost,
	       int x_along_p,int outwards_boost,
	       double *qmag,double *qx,double *qy,double *qz);
double ran2(void);
void gauss(double *g1,double *g2);
/* end of prototypes */

                    /* CLASS DEFINITIONS */

#define RESULTS_DIR

/* ************* These are used for binning **************** */

/* ********* For binning in one dimension ******************** */
class onedbin_class{ 
public:
  double *numerator,*denominator;
  int *numcounter,*dencounter;
  double *numrms,*denrms;
  int nkmax;
  double maxmom;
  char filename[80];
  onedbin_class(int nk,double biggest_mom,char *fn);
  void bincorr(int inumden,double mom,double corr);
  void printcorr(int onlyzerob);
};
onedbin_class::onedbin_class(int nk,double biggest_mom,char *fn){
  int i;
  nkmax=nk;
  maxmom=biggest_mom;
#ifdef RESULTS_DIR
  strcpy(filename,"results/");
  strcat(filename,fn);
#else
  strcpy(filename,fn);
#endif
  numerator=new double[nkmax];
  denominator=new double[nkmax];
  numcounter=new int[nkmax];
  dencounter=new int[nkmax];
  numrms=new double[nkmax];
  denrms=new double[nkmax];
  
  for(i=0;i<nkmax;i++){
    numerator[i]=0.0;
    denominator[i]=0.0;
    numcounter[i]=0;
    dencounter[i]=0;
    numrms[i]=0.0;
    denrms[i]=0.0;
  }
}

void onedbin_class::bincorr(int inumden,double mom,double corr){
  int kk;
  kk=(int)floor((double)nkmax*mom/maxmom);
  if(kk<nkmax){
    if(inumden==0){
      numerator[kk]=numerator[kk]+corr;
      numcounter[kk]=numcounter[kk]+1;
      numrms[kk]=numrms[kk]+corr*corr;
    }
    else{
      denominator[kk]=denominator[kk]+corr;
      dencounter[kk]=dencounter[kk]+1;
      denrms[kk]=denrms[kk]+corr*corr;
    }
  }
}

void onedbin_class::printcorr(int onlyzerob){
  int kk,counts,no_denominator;
  double mom,error,corr,average_numer,average_denom;
  double normalization,normcount;
  FILE *fptr;
  normalization=1.0;
  no_denominator=onlyzerob;
#ifdef MIXED_PAIRS_FOR_DENOM
  no_denominator=0;
  normalization=0.0;
  normcount=0;
  for(kk=0;kk<nkmax;kk++){
    normcount=normcount+dencounter[kk];
    normalization=normalization+denominator[kk];
  }
  if(normcount>0) normalization=normalization/(double)normcount;
  printf("For %s, normalization=%g \n",filename,normalization);
#endif
  fptr=fopen(filename,"w");
#ifdef REDUCED_MOM
  fprintf(fptr,"!As a function of reduced momentum, k=(p2-p1)/2:\n");
#else
  fprintf(fptr,"!As a function of relative momentum, k=(p2-p1):\n");
#endif
  fprintf(fptr,"!k(MeV/c)  C(k)     +/- numcounts&dencounts\n");
  for(kk=0;kk<nkmax;kk++){
    mom=(kk+0.5)*maxmom/(double)nkmax;

    corr=-1.0;
    if(no_denominator==1){
      if(numcounter[kk]>0) corr=numerator[kk]/numcounter[kk];
    }
    else{
      if(dencounter[kk]>0) corr=numerator[kk]/denominator[kk];
    }

    error=-1.0;
    if(numcounter[kk]!=0 && (no_denominator==1 || dencounter[kk]!=0)){
      average_numer=numerator[kk]/numcounter[kk];
      error=((numrms[kk]/numcounter[kk])
	     -(average_numer*average_numer))/numcounter[kk];
      if(no_denominator!=1) error=error+corr*corr/(double)numcounter[kk];
      if(no_denominator!=1){
	average_denom=denominator[kk]/dencounter[kk];
	error=error+((denrms[kk]/denominator[kk])
		     -(average_denom*average_denom))/dencounter[kk];
	if(no_denominator!=1) error=error+corr*corr/(double)dencounter[kk];
      }
      error=sqrt(fabs(error));
    }
    if(corr>=0.0){
      corr=corr*normalization;
      error=error*normalization;
    }

    fprintf(fptr,"%5.2f    %6.3f   %9g      %d  %d\n",
	    mom,corr,error,numcounter[kk],dencounter[kk]);
  }
  fclose(fptr);
}

/* ********* For binning in three dimensions ******************** */

class threedbin_class{
public:
  double *numerator,*denominator;
  int *numcounter,*dencounter;
  double *numrms,*denrms;
  int nkxmax,nkymax,nkzmax;
  double kxmax,kymax,kzmax;
  char filename[80];
  threedbin_class(int nkx,int nky,int nkz,
		  double kxmax,double kymax,double kzmax,char *fn);
  void bincorr(int inumden,double kx,double ky,double kz,double corr);
  void printcorr(int onlyzerob);
};
threedbin_class::threedbin_class(int nkx,int nky,int nkz,
				 double biggest_kx,double biggest_ky,
				 double biggest_kz,char *fn){
  int ix,iy,iz,i,nmax;
  nkxmax=nkx;
  nkymax=nky;
  nkzmax=nkz;
  kxmax=biggest_kx;
  kymax=biggest_ky;
  kzmax=biggest_kz;
  nmax=nkxmax*nkymax*nkzmax;
#ifdef RESULTS_DIR
  strcpy(filename,"results/");
  strcat(filename,fn);
#else
  strcpy(filename,fn);
#endif
  numerator=new double[nmax];
  denominator=new double[nmax];
  numcounter=new int[nmax];
  dencounter=new int[nmax];
  numrms=new double[nmax];
  denrms=new double[nmax];

  for(i=0;i<nmax;i++){
    numerator[i]=0.0;
    denominator[i]=0.0;
    numcounter[i]=0;
    dencounter[i]=0;
    numrms[i]=0.0;
    denrms[i]=0.0;
  }
}

void threedbin_class::bincorr(int inumden,double kx,double ky,double kz,
			      double corr){
  int kk,ix,iy,iz;
  ix=(int)floor((double)nkxmax*fabs(kx)/kxmax);
  iy=(int)floor((double)nkymax*fabs(ky)/kymax);
  iz=(int)floor((double)nkzmax*fabs(kz)/kzmax);
  if(ix<nkxmax && iy<nkymax && iz<nkzmax){
    kk=ix*(nkymax*nkzmax)+iy*nkzmax+iz;
    if(inumden==0){
      numerator[kk]=numerator[kk]+corr;
      numcounter[kk]=numcounter[kk]+1;
      numrms[kk]=numrms[kk]+corr*corr;
    }
    else{
      denominator[kk]=denominator[kk]+corr;
      dencounter[kk]=dencounter[kk]+1;
      denrms[kk]=denrms[kk]+corr*corr;
    }
  }
}

void threedbin_class::printcorr(int onlyzerob){
  int kk,ix,iy,iz,nmax,counts,no_denominator;
  double kx,ky,kz,error,corr,average_numer,average_denom;
  double normalization,normcount;
  FILE *fptr;
  normalization=1.0;
  no_denominator=onlyzerob;
#ifdef MIXED_PAIRS_FOR_DENOM
  no_denominator=0;
  normalization=0.0;
  normcount=0;
  nmax=nkxmax*nkymax*nkzmax;
  for(kk=0;kk<nmax;kk++){
    normcount=normcount+dencounter[kk];
    normalization=normalization+denominator[kk];
  }
  if(normcount>0) normalization=normalization/(double)normcount;
  printf("For %s, normalization=%g \n",filename,normalization);
#endif
  fptr=fopen(filename,"w");
#ifdef REDUCED_MOM
  fprintf(fptr,"!As a function of reduced momentum, k=(p2-p1)/2:\n");
#else
  fprintf(fptr,"!As a function of relative momentum, k=(p2-p1):\n");
#endif
  fprintf(fptr,"! kx    ky    kz       C(k)    +/-  numcounts&dencounts\n");
  kk=0;
  for(ix=0;ix<nkxmax;ix++){
    kx=(ix+0.5)*kxmax/(double)nkxmax;
    for(iy=0;iy<nkymax;iy++){
      ky=(iy+0.5)*kymax/(double)nkymax;
      for(iz=0;iz<nkzmax;iz++){
	kz=(iz+0.5)*kzmax/(double)nkzmax;

	corr=-1.0;
	if(no_denominator==1){
	  if(numcounter[kk]>0) corr=numerator[kk]/numcounter[kk];
	}
	else{
	  if(dencounter[kk]>0) corr=numerator[kk]/denominator[kk];
	}

	error=-1.0;
	if(numcounter[kk]!=0 && (no_denominator==1 || dencounter[kk]!=0)){
	  average_numer=numerator[kk]/numcounter[kk];
	  error=((numrms[kk]/numcounter[kk])
		 -(average_numer*average_numer))/numcounter[kk];
	  if(no_denominator!=1) error=error+corr*corr/(double)numcounter[kk];
	  if(no_denominator!=1){
	    average_denom=denominator[kk]/dencounter[kk];
	    error=error+((denrms[kk]/denominator[kk])
			 -(average_denom*average_denom))/dencounter[kk];
	    if(no_denominator!=1) error=error+corr*corr/(double)dencounter[kk];
	  }
	  error=sqrt(fabs(error));
	}
	if(corr>=0.0){
	  corr=corr*normalization;
	  error=error*normalization;
	}

	fprintf(fptr,"%5.2f %5.2f %5.2f    %6.3f   %9f      %d  %d\n",
		kx,ky,kz,corr,error,numcounter[kk],dencounter[kk]);
	kk=kk+1;
      }
    }
  }
  fclose(fptr);
}

/* ************* Phasespace Info Class **************** */

class phasespace_class{
public:
  double p[4],r[4];
  int iread;
  phasespace_class(double *pp,double *rr,int iiread);
  /* ~phasespace_class(); */
};
phasespace_class::phasespace_class(double *pp,double *rr,int iiread){
  int alpha;
  for(alpha=0;alpha<4;alpha++){
    p[alpha]=pp[alpha];
    r[alpha]=rr[alpha];
  }
  iread=iiread;
}

/* ************************************************************ */
/* The Following are used for the calculation of the strong-interaction
   corrections to wave functions*/

class partwave_class{
public:
  double_complex *delpsi;
  double_complex delpsifar;
  partwave_class(int nrmax);
};
partwave_class::partwave_class(int nrmax){
  delpsi=new double_complex[nrmax];
}

/* ************************************************************ */
/* The Following are used for the calculation of Coulomb wave functions*/

class coulcalc_class{
public:
  double_complex couly;
  double_complex cfact1;
  double_complex cfact2;
  double_complex chype[41];
  double_complex chype1[6];
  double_complex chype2[6];
};

/* ************************************************************* */
