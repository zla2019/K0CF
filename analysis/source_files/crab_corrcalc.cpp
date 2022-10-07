#define ROOT2 1.41421356237309504880
double  corrcalc(double trueqred,double trueqdotr,double truer){
  double eta,arg,krmax,corr0;
  double xx,xxprime,xxjj,p1,zk;
  int jj,kk,ipart,ipartcount,ispin;
  double pfactor,wsym_leftover,wanti_leftover,wnosym_leftover;
  double qred,qdotr,r;
  const double rmass=MASS1*MASS2/(MASS1+MASS2);
  double_complex cphi1,cphi2,cphis,cphia;
#ifdef STRONG_INTERACTION
  double_complex cfar,cphi[STRONG_NSPINS];
#endif
#if ! defined STRONG_INTERACTION
#if ! defined COULOMB
  arg=trueqdotr/197.323-2.0*pi*floor(trueqdotr/(197.323*2.0*pi));
  cphi1=exp(ci*arg);
  cphis=ROOT2*real(cphi1);
  cphia=ci*ROOT2*imag(cphi1);
  corr0=real(INTERACTION_WSYM*cphis*conj(cphis)
	     +INTERACTION_WANTI*cphia*conj(cphia)
	     +INTERACTION_WNOSYM*cphi1*conj(cphi1));
  goto OUTSIDE_INTERACTION_RANGE;
#endif
#endif

#ifdef REDUCED_MOM
  kk=(int)floor(trueqred/INTERACTION_DELK);
  qred=(0.5+kk)*INTERACTION_DELK;
#else
  kk=(int)floor(2.0*trueqred/INTERACTION_DELK);
  qred=(0.5+kk)*INTERACTION_DELK/2.0;
#endif
  qdotr=trueqdotr*qred/trueqred;
  if(kk>=INTERACTION_NKMAX){
    corr0=1.0;
    goto OUTSIDE_INTERACTION_RANGE;
  }
#ifdef STRONG_INTERACTION
  jj=(int)floor(STRONG_NRMAX*truer/STRONG_RMAX);
  r=(0.5+(double)jj)*STRONG_RMAX/STRONG_NRMAX;
  qdotr=trueqdotr*r/truer;
#else
  r=truer;
#endif

#ifdef COULOMB
  zk=qdotr/qred;
  eta=(double)Q1Q2*(rmass/137.036)/qred;
  cphi1=coul( qred,r,zk,eta,kk);
  cphi2=coul( qred,r,-zk,eta,kk);
//printf("######## end COULOMB #########\n");
#else
  eta=0.0;
  arg=qdotr/197.323-2.0*pi*floor(qdotr/(197.323*2.0*pi));
  cphi1=exp(ci*arg);
  cphi2=conj(cphi1);
#endif

  cphis=(cphi1+cphi2)/ROOT2;
  cphia=(cphi1-cphi2)/ROOT2;
  corr0=0.0;
  /* If there are corrections for strong interactions, add the
     change for each partial wave.  If npartial = 0 then there
     are no strong int. corrections. */
  wsym_leftover=INTERACTION_WSYM;
  wanti_leftover=INTERACTION_WANTI;
  wnosym_leftover=INTERACTION_WNOSYM;

#ifdef STRONG_INTERACTION
//printf("######## start STRONG_INTERACTION #########\n");
  xx=r*qred/197.323;
  p1=qdotr/(r*qred);
  ipart=0;
  for(ispin=0;ispin<STRONG_NSPINS;ispin++){
    if(STRONG_SYMM[ispin]==-1){
      pfactor=ROOT2;
      cphi[ispin]=cphia;
      wanti_leftover=wanti_leftover-STRONG_WEIGHT[ispin];
    }
    if(STRONG_SYMM[ispin]==0){
      pfactor=1.0;
      cphi[ispin]=cphi1;
      wnosym_leftover=wnosym_leftover-STRONG_WEIGHT[ispin];
    }
    if(STRONG_SYMM[ispin]==1){
      pfactor=ROOT2;
      cphi[ispin]=cphis;
      wsym_leftover=wsym_leftover-STRONG_WEIGHT[ispin];
    }
    for(ipartcount=0;ipartcount<STRONG_NPARTIALS[ispin];ipartcount++){
	std::cout << "ipartcount: " << ipartcount << std::endl;
      if(jj<STRONG_NRMAX){
	std::cout << "ispin: " << ispin << " kk: " << kk << " jj: " << jj << std::endl;
	if(STRONG_L[ipart]==0) cphi[ispin]=cphi[ispin]
				 +pfactor*partwave[ipart][kk]->delpsi[jj]/xx;
	if(STRONG_L[ipart]==1) cphi[ispin]=cphi[ispin]
				 +pfactor*p1*3*ci*partwave[ipart][kk]->delpsi[jj]/xx;
	if(STRONG_L[ipart]==2) cphi[ispin]=cphi[ispin]
				 -pfactor*0.5*(3*(p1*p1)-1.0)*
				 5*partwave[ipart][kk]->delpsi[jj]/xx;
	/*  ___________________________ */
      }
      else{
	xxprime=xx-2*pi*floor(xx/(2.0*pi));
	cfar=exp(ci*xxprime-ci*eta*log(2.0*xx));
	if(STRONG_L[ipart]==0) cphi[ispin]=cphi[ispin]
				 -ci*pfactor*partwave[ipart][kk]->delpsifar*cfar/xx;
	if(STRONG_L[ipart]==1) cphi[ispin]=cphi[ispin]
				 -pfactor*p1*3*partwave[ipart][kk]->delpsifar*cfar/xx;
	if(STRONG_L[ipart]==2) cphi[ispin]=cphi[ispin]
				 +ci*pfactor*5*0.5*(3*(p1*p1)-1.0)*partwave[ipart][kk]->delpsifar*cfar/xx;
      }
      ipart=ipart+1;
	std::cout << "ipart: " << ipart << std::endl;
    }
  }
  for(ispin=0;ispin<STRONG_NSPINS;ispin++){
    corr0=corr0+real(STRONG_WEIGHT[ispin]*cphi[ispin]*conj(cphi[ispin]));
  }

//printf("######## finish STRONG_INTERACTION #########\n");
#endif
  corr0=corr0+real(wsym_leftover*cphis*conj(cphis)
		   +wanti_leftover*cphia*conj(cphia)
		   +wnosym_leftover*cphi1*conj(cphi1));
OUTSIDE_INTERACTION_RANGE:
#ifdef BREIT_WIGNER
  corr0=corr0+bwcalc(trueqred,truer);
#endif

  return corr0;
}

#ifdef BREIT_WIGNER

double bwcalc(double qred,double r){
  int ires;
  double e1,e2,v,e,ddelde,gamma,dgammade,delm2,answer,weight,total;
  double bwmom,factor1;
  static double factor2=pow(197.323/(BWINFO_R*sqrt(2.0*pi)),3.0);
  e1=sqrt(MASS1*MASS1+qred*qred);
  e2=sqrt(MASS2*MASS2+qred*qred);
  e=e1+e2;
  v=(qred/e1)+(qred/e2);
  total=0.0;
  for(ires=0;ires<N_BWRESONANCES;ires++){
    bwmom=sqrt(triangle(BWINFO_MINV[ires],MASS1,MASS2));
    factor1=(2*BWINFO_J[ires]+1)/((2*BWINFO_J1[ires]+1)*(2*BWINFO_J2[ires]+1));
    gamma=pow(qred/bwmom,(double)(2*BWINFO_L[ires]+1))
      *BWINFO_GAMMA[ires]*pow(BWINFO_MINV[ires],2.0)/e;
    delm2=BWINFO_MINV[ires]*BWINFO_MINV[ires]-e*e;
    dgammade=-(gamma/e)+3.0*gamma/(qred*v);
    ddelde=(2.0*e*gamma+delm2*dgammade)/(gamma*gamma+delm2*delm2);
    if((r/BWINFO_R)<10.0){
      weight=factor2*exp(-r*r/(2.0*BWINFO_R*BWINFO_R));
    }
    else{
      weight=0.0;
    }
    answer=(2.0*pi/(qred*qred))*factor1*v*ddelde*weight;
    answer=answer*BWINFO_BRANCHING[ires];
    total=total+answer;
  }
  return total;
}

double triangle(double m0,double m1,double m2){
  double answer,m0sq,m1sq,m2sq;
  m0sq=m0*m0;m1sq=m1*m1;m2sq=m2*m2;
  answer=m0sq*m0sq+m1sq*m1sq+m2sq*m2sq;
  answer=answer-2.0*(m0sq*m1sq+m0sq*m2sq+m1sq*m2sq);
  answer=answer/(4.0*m0sq);
  return answer;
}

#endif
