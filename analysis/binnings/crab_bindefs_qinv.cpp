/* Set ABSOLUTE_MAXMOM so that program won't bother calculating correlations
   when rel. mom. is greater than ABSOLUTE_MAXMOM. This is ignored if
   "MIXED_PAIRS_FOR_DENOM" is defined in crab.cpp below */
#define ABSOLUTE_MAXMOM 100

/* Four Examples are presented here (A-D). For each binning, modifications
   must be made in four locations in this file as illustrated. One can
   add as many binnings as desired. Binnings are of two types:
   1.) one-dimensional, where the answer will be given as C(x).
   2.) three-dimensional, where the answer will be given as C(x,y,z). */

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
   one also boost along the direction of p1+p2 to the pair's c.o.m. frame. */

/* Example A: Bin the correlation function according to Qinv */
/* Example B: Bin the correlation function according to Qinv
   in 10 different pt bins */
/* Example C: Bin the correlation function according to |Q|, where Q is
   determined in the longitudinally comoving frame. Make three binnings for
   three different directional cuts. */
/* Example D: Save the correlation function in a 3-d binning, according
   to Qx, Qy, and Qz. */

/* ********************************************************* */

/* First define the binnings you will do. */

/* A. */
onedbin_class *qinv;

/* B. */
//#define NPTCUTS 10
//#define DELPT 50
//onedbin_class *qinv_ptcut[NPTCUTS];

/* C */
//onedbin_class *q_outwards;
//onedbin_class *q_sidewards;
//onedbin_class *q_beam;

/* D. */
//threedbin_class *qinv3d;

/* ********************************************************* */

/* Initialize binning objects with desired mesh sizes and momentum ranges */

void bininit(void){
  int ipt,nkmax,nkxmax,nkymax,nkzmax;
  int nptmax;
  double maxmom,kxmax,kymax,kzmax;
  double pt;
  char title[80];
//printf("bininit\n");

  /* A */
  nkmax=40;
  maxmom=100;
  qinv=new onedbin_class(nkmax,maxmom,"qinv.dat");

  /* B */
  /*  nkmax=50;
      maxmom=100.0;
      for(ipt=0;ipt<NPTCUTS;ipt++){
      pt=(0.5+ipt)*maxmom/NPTCUTS;
      pt=floor(pt+0.5);
      sprintf(title,"qinv_pt%03g.dat",pt);
      qinv_ptcut[ipt]=new onedbin_class(nkmax,maxmom,title);
      } */

      /* C */
  /* nkmax=50;
     maxmom=100.0;
     q_outwards=new onedbin_class(nkmax,maxmom,"q_outwards.dat");
     q_sidewards=new onedbin_class(nkmax,maxmom,"q_sidewards.dat");
     q_beam=new onedbin_class(nkmax,maxmom,"q_beam.dat"); */

  /* D */
  /* nkxmax=nkymax=nkzmax=10;
     kxmax=kymax=kzmax=50.0;
     qinv3d=new threedbin_class(nkxmax,nkymax,nkzmax,
     kxmax,kymax,kzmax,"qinv3d.dat"); */

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
  /* pt=sqrt(pow(p1[1]+p2[1],2)+pow(p1[2]+p2[2],2));
     ipt=(int)floor(pt/DELPT);
     if(ipt<NPTCUTS){
     qinv_ptcut[ipt]->bincorr(inumden,mom,corr);
     } */

  /* C */
  /* long_comoving_boost=1;
     x_along_p=0;
     outwards_boost=0;
     decompose(p1,p2,vs,
     long_comoving_boost,x_along_p,outwards_boost,
     &kmag,&kx,&ky,&kz);
     if(fabs(ky)<5.0 && fabs(kz)<5.0) q_outwards->bincorr(inumden,kmag,corr);
     if(fabs(kx)<5.0 && fabs(kz)<5.0) q_sidewards->bincorr(inumden,kmag,corr);
     if(fabs(kx)<5.0 && fabs(ky)<5.0) q_beam->bincorr(inumden,kmag,corr); */


  /* D */
  /* long_comoving_boost=1;
     x_along_p=0;
     outwards_boost=0;
     decompose(p1,p2,vs,
     long_comoving_boost,x_along_p,outwards_boost,
     &kmag,&kx,&ky,&kz);
     qinv3d->bincorr(inumden,kx,ky,kz,corr); */

}

void printstatements(int onlyzerob){
  int ipt;

  /* A */
  qinv->printcorr(onlyzerob);

  /* B */
  /* for(ipt=0;ipt<10;ipt++){
     qinv_ptcut[ipt]->printcorr(onlyzerob);
     } */

  /* C */
  /* q_outwards->printcorr(onlyzerob);
     q_sidewards->printcorr(onlyzerob);
     q_beam->printcorr(onlyzerob); */

  /* D */
  /* qinv3d->printcorr(onlyzerob); */

}

