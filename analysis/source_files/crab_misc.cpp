void ranrot(double *p,double *r,double phi){
  double c,s,px0,rx0;
  c=cos(phi);
  s=sin(phi);
  px0=p[1];
  rx0=r[1];
  p[1]=c*px0+s*p[2];
  r[1]=c*rx0+s*r[2];
  p[2]=c*p[2]-s*px0;
  r[2]=c*r[2]-s*rx0;
}

/* ****************************************************************** */

void get_pr(int ib1,int j1,double phi1,int ib2,int j2,double phi2,
	    double *p1,double *r1,int *iread1,
	    double *p2,double *r2,int *iread2){
  int alpha;
  for(alpha=0;alpha<4;alpha++){
    p1[alpha]=phasespace1_ptr[ib1][j1]->p[alpha];
    r1[alpha]=phasespace1_ptr[ib1][j1]->r[alpha];
  }
  *iread1=phasespace1_ptr[ib1][j1]->iread;
  ranrot(p1,r1,phi1);
#if defined IDENTICAL
  for(alpha=0;alpha<4;alpha++){
    p2[alpha]=phasespace1_ptr[ib2][j2]->p[alpha];
    r2[alpha]=phasespace1_ptr[ib2][j2]->r[alpha];
  }
  *iread2=phasespace1_ptr[ib2][j2]->iread;
#else
  for(alpha=0;alpha<4;alpha++){
    p2[alpha]=phasespace2_ptr[ib2][j2]->p[alpha];
    r2[alpha]=phasespace2_ptr[ib2][j2]->r[alpha];
  }
  *iread2=phasespace2_ptr[ib2][j2]->iread;
#endif
  ranrot(p2,r2,phi2);
}

/* ****************************************************************** */

void get_b_num(int *ib1,int *ib2,double *bweight_num){
  int ib;
  double randy;
  randy=ran2();
  ib=0;
  while(randy>bweight_num[ib]){
    ib=ib+1;
  }
  *ib1=ib;
  *ib2=ib;
}

/* ****************************************************************** */

void get_b_den(int *ib1,int *ib2,double *bweight1_den,double *bweight2_den){
  int ib;
  double randy;
  randy=ran2();
//printf("randy = %f\n",randy);
  ib=0;
for (int i = 0; i < 5; i++)
{
//printf("bweight1_den = %f\n",bweight1_den[i]);
//printf("bweight2_den = %f\n",bweight2_den[i]);
}
  while(randy>bweight1_den[ib]){
    ib=ib+1;
  }
  *ib1=ib;
  ib=0;
  while(randy>bweight2_den[ib]){
    ib=ib+1;
  }
  *ib2=ib;
}

/* ****************************************************************** */

void get_com_quantities(double *p1,double *r1,double *p2,double *r2,
			double *qred,double *r,double *qdotr,double *mom,
			int *test){
  int alpha;
  double kdotr;
#if defined REDUCED_MOM
  const double momtest=4.0*ABSOLUTE_MAXMOM*ABSOLUTE_MAXMOM;
#else
  const double momtest=ABSOLUTE_MAXMOM*ABSOLUTE_MAXMOM;
#endif
  double ptot2,pdotr,pp,rr;
#if defined IDENTICAL
  *test=1;
  *mom=-(p2[0]-p1[0])*(p2[0]-p1[0]);
  for(alpha=1;alpha<4;alpha++){
    *mom=*mom+(p2[alpha]-p1[alpha])*(p2[alpha]-p1[alpha]);
  }
#if ! defined MIXED_PAIRS_FOR_DENOM
  if(*mom>momtest){
    *test=0;
    goto FAILED_MOMTEST;
  }
#endif
  pp=(p1[0]+p2[0]);
  rr=(r2[0]-r1[0]);
  pdotr=pp*rr;
  kdotr=(p2[0]-p1[0])*rr;
  ptot2=pp*pp;
  *r=-rr*rr;
  for(alpha=1;alpha<4;alpha++){
    pp=(p1[alpha]+p2[alpha]);
    rr=(r2[alpha]-r1[alpha]);
    pdotr=pdotr-pp*rr;
    kdotr=kdotr-(p2[alpha]-p1[alpha])*rr;
    ptot2=ptot2-pp*pp;
    *r=*r+rr*rr;
  }
  *mom=sqrt(*mom);
  *qred=0.5**mom;
#ifdef REDUCED_MOM
  *mom=*qred;
#endif
  *qdotr=0.5*kdotr;
  *r=sqrt(*r+pdotr*pdotr/ptot2);
#else
  const double  kdotp=MASS2*MASS2-MASS1*MASS1;
  *test=1;
  *mom=-(p2[0]-p1[0])*(p2[0]-p1[0]);
  ptot2=(p1[0]+p2[0])*(p1[0]+p2[0]);
  for(alpha=1;alpha<4;alpha++){
    *mom=*mom+(p2[alpha]-p1[alpha])*(p2[alpha]-p1[alpha]);
    ptot2=ptot2-(p1[alpha]+p2[alpha])*(p1[alpha]+p2[alpha]);
  }
  *mom=*mom+kdotp*kdotp/ptot2;
#if ! defined MIXED_PAIRS_FOR_DENOM
  if(*mom>momtest){
    *test=0;
    goto FAILED_MOMTEST;
  }
#endif
  pp=(p1[0]+p2[0]);
  rr=(r2[0]-r1[0]);
  pdotr=pp*rr;
  kdotr=(p2[0]-p1[0])*rr;
  *r=-rr*rr;
  for(alpha=1;alpha<4;alpha++){
    pp=(p1[alpha]+p2[alpha]);
    rr=(r2[alpha]-r1[alpha]);
    pdotr=pdotr-pp*rr;
    kdotr=kdotr-(p2[alpha]-p1[alpha])*rr;
    *r=*r+rr*rr;
  }
  kdotr=(-kdotr+kdotp*pdotr/ptot2);
  *mom=sqrt(*mom);
  *qred=0.5**mom;
#ifdef REDUCED_MOM
  *mom=*qred;
#endif
  *qdotr=0.5*kdotr;
  *r=sqrt(*r+pdotr*pdotr/ptot2);
#endif
FAILED_MOMTEST:
  return;
}

/* ****************************************************************** */

double get_mom(double *p1,double *p2){
  const double ptotdotk=MASS2*MASS2-MASS1*MASS1;
  double ptot2,mom;
  int alpha;
  ptot2=(p1[0]+p2[0])*(p1[0]+p2[0]);
  mom=(p2[0]-p1[0])*(p2[0]-p1[0]);
  for(alpha=1;alpha<4;alpha++){
    ptot2=ptot2-(p1[alpha]+p2[alpha])*(p1[alpha]+p2[alpha]);
    mom=mom-(p2[alpha]-p1[alpha])*(p2[alpha]-p1[alpha]);
  }
  mom=sqrt(-mom+ptotdotk*ptotdotk/ptot2);
#ifdef REDUCED_MOM
  mom=mom/2.0;
#endif
  return mom;
}

/* ****************************************************************** */

void decompose(double *p1,double *p2,double vs,
	       int long_comoving_boost,
	       int x_along_p,int outwards_boost,
	       double *qmag,double *qx,double *qy,double *qz){
  int alpha;
  double q[4],ptot[4],qprime[4],u[4],minv;
  double oldz,gamma,v_out,oldout;
  double qperp,qside,qbeam,ptot_mag,ptot_perp;
  double q_outwards,q_sidewards,q_inplane;

  for(alpha=0;alpha<4;alpha++){
    q[alpha]=(p2[alpha]-p1[alpha]);
    ptot[alpha]=p1[alpha]+p2[alpha];
  }

  if(fabs(vs)>0.0001 && long_comoving_boost!=1){
    gamma=1.0/sqrt(1.0-vs*vs);
    oldz=ptot[3];
    ptot[3]=gamma*(ptot[3]-vs*ptot[0]);
    ptot[0]=gamma*(ptot[0]-vs*oldz);
    oldz=q[3];
    q[3]=gamma*(q[3]-vs*q[0]);
    q[0]=gamma*(q[0]-vs*oldz);
  }

  if(long_comoving_boost==1){
    vs=ptot[3]/ptot[0];
    gamma=1.0/sqrt(1.0-vs*vs);
    ptot[0]=gamma*(ptot[0]-vs*ptot[3]);
    ptot[3]=0.0;
    oldz=q[3];
    q[3]=gamma*(q[3]-vs*q[0]);
    q[0]=gamma*(q[0]-vs*oldz);
  }

  ptot_perp=sqrt(ptot[1]*ptot[1]+ptot[2]*ptot[2]);
  q_outwards=(q[1]*ptot[1]+q[2]*ptot[2])/ptot_perp;
  q_sidewards=sqrt(q[1]*q[1]+q[2]*q[2]-q_outwards*q_outwards);
  q_inplane=q[3];

  if(x_along_p==1){
    ptot_mag=sqrt(ptot[3]*ptot[3]+ptot_perp*ptot_perp);
    oldout=q_outwards;
    q_outwards=(q_outwards*ptot_perp+q_inplane*ptot[3])/ptot_mag;
    q_inplane=sqrt(oldout*oldout+q[3]*q[3]-q_outwards*q_outwards);
  }

  if(outwards_boost==1){
    if(x_along_p==1 || long_comoving_boost==1){
      v_out=ptot_mag/ptot[0];
      gamma=1.0/sqrt(1.0-v_out*v_out);
      q_outwards=gamma*(q_outwards-v_out*q[0]);
    }
    else{
      q[1]=q_outwards;
      q[2]=q_sidewards;
      q[3]=q_inplane;
      minv=ptot[0]*ptot[0]-ptot[1]*ptot[1]-ptot[2]*ptot[2]-ptot[3]*ptot[3];
      minv=sqrt(minv);
      u[1]=-ptot_perp/minv;
      u[2]=0.0;
      u[3]=-ptot[3]/minv;
      u[0]=ptot[0]/minv;
      lorentz(u,q,qprime);
      q_outwards=qprime[1];
      q_sidewards=qprime[2];
      q_inplane=qprime[3];
    }
  }
#ifdef REDUCED_MOM
  *qx=0.5*q_outwards;
  *qy=0.5*q_sidewards;
  *qz=0.5*q_inplane;
  *qmag=sqrt(*qx**qx+*qy**qy+*qz**qz);
#else
  *qx=q_outwards;
  *qy=q_sidewards;
  *qz=q_inplane;
  *qmag=sqrt(*qx**qx+*qy**qy+*qz**qz);
#endif
}

/* ****************************************************************** */

#ifdef MIXED_PAIRS_FOR_DENOM
void store_mixed_pair(int n_mixed_pairs,double *p1,double *pp1,double *p2, double *pp2,int iread1,int iread2,int *iiread1,int *iiread2,double corr,double *ww){
  int alpha;
  for(alpha=0;alpha<4;alpha++){
    pp1[alpha]=p1[alpha];
    pp2[alpha]=p2[alpha];
  }
  iiread1[n_mixed_pairs]=iread1;
  iiread2[n_mixed_pairs]=iread2;
  ww[n_mixed_pairs]=corr;
}
#endif

/* ****************************************************************** */

void lorentz(double *u,double *q,double *qprime){
  static int n[4]={1,0,0,0};
  static int g[4]={1,-1,-1,-1};
  int mu;
  double udotn=0.0,qdotn=0.0,qdotu=0.0;
  for(mu=0;mu<4;mu++){
    qdotn=qdotn+q[mu]*n[mu]*g[mu];
    qdotu=qdotu+q[mu]*u[mu]*g[mu];
    udotn=udotn+u[mu]*n[mu]*g[mu];
  }
  for(mu=0;mu<4;mu++){
    qprime[mu]=-((qdotu+qdotn)/(1.0+udotn))*(n[mu]+u[mu])
      +2*qdotn*u[mu]+q[mu];
  }
}

/* ****************************************************************** */

double get_gamow(double mom){
  const double rmass=MASS1*MASS2/(MASS1+MASS2);
  const double alpha=1.0/137.036;
  double qred,kmom,eta,gamow;
  int kk;
  kk=(int)floor(mom/INTERACTION_DELK);
  kmom=((double)kk+0.5)*INTERACTION_DELK;
#ifdef REDUCED_MOM
  qred=kmom;
#else
  qred=0.5*kmom;
#endif
  eta=rmass*alpha*(double)Q1Q2/qred;
  gamow=2.0*pi*eta/(exp(2.0*pi*eta)-1.0);
  return gamow;
}
