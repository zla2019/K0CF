/* Set ABSOLUTE_MAXMOM so that program won't bother calculating correlations
   when rel. mom. is greater than ABSOLUTE_MAXMOM. This is ignored if
   "MIXED_PAIRS_FOR_DENOM" is defined in crab.cpp below */
#define ABSOLUTE_MAXMOM 50

/* Four Examples are presented here (A-D). For each binning, modifications
   must be made in four locations in this file as illustrated. One can
   add as many binnings as desired. Binnings are of two types: one-d, where
   the answer will be given as C(x), and three-d, where the answer will
   be given as C(x,y,z). */

/* One of the included routines is 'decompose(p1,p2,vs,
   long_comoving_boost,x_along_p,outwards_boost,
   &kmag,&kx,&ky,&kz)'. This routine returns the 3-dimensions of the
   relative momentum. ( If one only uses Qinv, one can skip calling this
   function since 'mom' equals 'Qinv' and is already fed into 'binner'. )
   When one calls decompose, three additional arguments,
   long_comoving_boost,x_along_p,outwards_boost, are used to determine the
   reference frame, (both relativistically and rotationally).
   1.) If one chooses 'long_comoving_boost' equal to one,
   one boosts along the beam axis to the frame where pz1+pz2=0.
   If 'long_comoving_boost' is not equal to one, the longitudinal frame
   will be set by the value of 'vs'.
   2.) If one chooses 'x_along_p' equal to one,
   one performs a rotation, in the reaction plane defined by the pairs
   momentum, such that the x-axis is parallel to the total momentum.
   3.) If one chooses outwards_boost equal to one,
   one also boost along the direction of p1+p2 to the pair's c.o.m. frame.

/* Example A: Bin the correlation function according to Qinv */
/* Example B: Bin the correlation function according to Qinv
   in 10 different pt bins */
/* Example C: Bin the correlation function according to |Q|, where Q is
   determined in the longitudinally comoving frame. Make three binnings for
   three different directional cuts. */
/* Example D: Save the correlation function in a 3-d binning, according
   to Qx, Qy, and Qz.

/* ********************************************************* */

/* First define the binnings you will do. */

/* A. */
onedbin_class *qinv;

/* C */
onedbin_class *q_plus;
onedbin_class *q_minus;


/* ********************************************************* */

/* Initialize binning objects with desired mesh sizes and momentum ranges */

void bininit(void){
  int ipt,nkmax,nkxmax,nkymax,nkzmax;
  int nptmax;
  double maxmom,kxmax,kymax,kzmax;
  double pt;
  char title[80];


  /* A */
  nkmax=50;
  maxmom=50.0;
  qinv=new onedbin_class(nkmax,maxmom,"qinv.dat");


  /* B */
  nkmax=50;
  maxmom=50.0;
  q_plus=new onedbin_class(nkmax,maxmom,"q_plus.dat");
  q_minus=new onedbin_class(nkmax,maxmom,"q_minus.dat");
}

/* ********************************************************* */

/* Bin |phi|^2 according to p1 and p2 */
/* Use the routine "decompose" to get projections of rel. mom. if desired */

void binner(int inumden,double mom,double corr,double *p1,double *p2){
  double pt,kx,ky,kz,kmag;
  int ipt;
  /* These 3 variables define conv.s used by decompose, 0=no, 1=yes */
  int long_comoving_boost,x_along_p,outwards_boost;
  double vs=0.0;

  /* A */
  qinv->bincorr(inumden,mom,corr);



  /* B */
  kx=p2[1]-p1[1];
  ky=p2[2]-p1[2];
  kz=p2[3]-p1[3];
  if(fabs(ky)<3.0 && fabs(kx)<3.0 && kz>0.0)
    q_plus->bincorr(inumden,mom,corr);
  if(fabs(ky)<3.0 && fabs(kx)<3.0 && kz<0.0)
    q_minus->bincorr(inumden,mom,corr);

}

void printstatements(int onlyzerob){
  int ipt;

  /* A */
  qinv->printcorr(onlyzerob);


  /* B */
  q_plus->printcorr(onlyzerob);
  q_minus->printcorr(onlyzerob);

}






