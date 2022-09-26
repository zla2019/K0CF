void partwaveinit(){
  int ipart;
  for(ipart=0;ipart<STRONG_NETPARTIALS;ipart++){
    diffy(ipart);
  }
}
/* *********************************** */
void diffy(int ipart){
  double delr,r00,rmax,mom,phase,phase0,phase1,phase2;
  double vv,r1,r2,r0,vcs;
  double b,eta,qred;
  const double rmass=MASS1*MASS2/(MASS1+MASS2);
  double_complex  cf[STRONG_NRMAX],cf0[STRONG_NRMAX],cfin[2*STRONG_NRMAX+1];
  double_complex cy0,cy1,cy2,ceip;
  int kk,j,mm,nmax,istrong,l;
  FILE *phaseshift_file;
//printf("##########partwaveinit########\n");
  if(ipart==0){
    phaseshift_file=fopen("results/phaseshifts.dat","w");
  }
  else{
    phaseshift_file=fopen("results/phaseshifts.dat","a");
  }

  nmax=2*STRONG_NRMAX;
  l=STRONG_L[ipart];
  fprintf(phaseshift_file,
	  "---- Calculating corrections to interaction for l=%d ----\n",l);
  fprintf(phaseshift_file,"q_red(MeV/c) delta(degrees)\n");

  for(kk=0;kk<INTERACTION_NKMAX;kk++){
    partwave[ipart][kk]=new partwave_class(STRONG_NRMAX);
    mom=((double)kk+0.5)*INTERACTION_DELK;
#ifdef REDUCED_MOM
    qred=mom;
#else
    qred=0.5*mom;
#endif
    rmax=STRONG_RMAX*qred/197.323;
    delr=-rmax/(double)nmax;
    r00=rmax;
#ifdef COULOMB
    eta=(double)Q1Q2*(rmass/137.036)/qred;
#else
    eta=0.0;
#endif
    /* +++++++++++++++ W & WO POTENTIAL LOOP +++++++++++++++++ */
    /* We perform calc with and with out the Reid pot.
       We are only interested in the correction to the wave function
       due to interactions. The first
       time through we turn off the Reid pot.(istrong=1,off) (istrong=2,on) */
    for(istrong=0;istrong<2;istrong++){
      r0=r00;
      r1=r0+delr;
      /* This initiates an incoming coul. wave function at
	 two points at r=10.  See Messiah's appendix on incoming partial
	 Coul. waves for explanation.  We then calc. inward until r=0. */
#if ! defined COULOMB
      cy0=cw2_big_r(l,r0,eta)*exp(-ci*r0);
      cy1=cw2_big_r(l,r1,eta)*exp(-ci*r1);
#else
      if(r0>10.0){
	cy0=cw2_big_r(l,r0,eta)*exp(-ci*r0);
	cy1=cw2_big_r(l,r1,eta)*exp(-ci*r1);
      }
      else{
	cy0=cw2_small_r(l,r0,eta);
	cy1=cw2_small_r(l,r1,eta);
      }
#endif

      cy0=cy0*pow(ci,(double)l);
      cy1=cy1*pow(ci,(double)l);

      cfin[nmax]=cy0;
      cfin[nmax-1]=cy1;
      /* Given the wave func. at two pts the wave func. can be found at
	 third pt. from diff. eq. */
      for(j=nmax-2;j>0;j=j-1){
	r2=r1+delr;
	vcs=1.0;
	if(istrong==0){
	  vv=0.0;
	}
	else{
	  vv=(2.0*rmass/(qred*qred))*POTENTIAL(197.323*r1/qred,ipart);
	}
	/* For screened coulomb */
	/* rvcs=(r1/qred)*197.323;
	   vcs=(rvcs/4.0,3.0) */
	b=1.0-((double)(l*(l+1))/(r1*r1))-(2.0*eta*vcs/r1)-vv;
	cy2=(-b*cy1*(delr*delr)+2*cy1-cy0);
	cfin[j]=cy2;
	r1=r2;
	cy0=cy1;
	cy1=cy2;
      }
      /*phase0 is the phase of the wavefunction at r=0.
	Subtracting the double_complex
	conjugate of the wavefunction *exp(-2i*phase0)
	from the wavefunction
	yields a sol. which satisfies the bound. conditions.Multiply by i
	to correspond to the correct incoming partia-wave phase(
	See Messiah).
	finally, when coming through after doing it for both with and
	withoutinteractions, take the difference. */
      phase1=real(-ci*log(cfin[1]/abs(cfin[1])));
      phase2=real(-ci*log(cfin[2]/abs(cfin[2])));
      phase0=2.0*phase1-phase2;
      ceip=exp(ci*phase0);
      /* Note that we send back to the main prog. the answer with a less
	 dense mesh */
      j=0;
      for(mm=1;mm<2*STRONG_NRMAX;mm=mm+2){
	if(istrong==0){
	  cf0[j]=ci*(cfin[mm]-conj(cfin[mm])*pow(ceip,2))/2;
	}
	if(istrong==1){
	  cf[j]=ci*(cfin[mm]-conj(cfin[mm])*pow(ceip,2))/2;
	  /* For large eta, strong int. correction should be neglible
	     but causes numerical instabilities, so we ignore strong
	     interaction in this case */
	  if(fabs(eta)<pi){
	    partwave[ipart][kk]->delpsi[j]=cf[j]-cf0[j];
	  }
	  else{
	    partwave[ipart][kk]->delpsi[j]=0.0;
	  }
	}
	j=j+1;
      }
      if(istrong==0) phase=phase0;
      if(istrong==1) {
	phase=phase0-phase;
	fprintf(phaseshift_file,"%5.2f    %g\n",qred,phase*180.0/pi);
      }
    }
#ifdef COULOMB
    partwave[ipart][kk]->delpsifar=
      (exp(2.0*ci*phase)-1.0)*cgamma(1.0+(double)l+ci*eta)/2.0;
#else
    partwave[ipart][kk]->delpsifar=
      (exp(2.0*ci*phase)-1.0)*dgamma(l+1)/2.0;
#endif
  }
  fclose(phaseshift_file);
}
/* ********************************** */
double_complex cw2_big_r(int l,double r,double eta){
  double_complex z,a,b,f1,top1,top2,bot,delf,answer;
  double arg;
  int n,j;
  z=-2*r*ci;
  a=(double)l+1.0+ci*eta;
  b=2.0*((double)l+1.0);
  f1=1.0;
  top1=1.0;
  top2=1.0;
  bot=1.0;
  n=5;
  for(j=1;j<=n;j++){
    top1=top1*((double)j-a);
    top2=top2*((double)j+b-a-1.0);
    bot=bot*(double)j;
    delf=(top1*top2)/(bot*pow(z,j));
    f1=f1+delf;
  }
  arg=eta*log(2.0*r);
  arg=arg-2.0*pi*floor(arg/(2.0*pi));
  answer=exp(ci*arg)*f1;
  return answer;
}
/* ************************************* */
double_complex cw2_small_r(int l,double r,double eta){
  /* The notation is like Gr. + R, page 1063.
     The Coulomb wave function is the same as W(i*eta,l+1/2,2*i*rho) */
  double_complex factor,lterm,fact1,fact2;
  double_complex psi1,psi2,psi3,cx,sum1,sum2,delsum1,delp,cdcon1,cdcon2,answer;
  int k;
#ifdef COULOMB
  cdcon1=cgamma(-(double)(l)-ci*eta);
  cdcon2=cgamma((double)(l+1)-ci*eta);
#else
  cdcon1=dgamma(-l);
  cdcon2=dgamma(l+1);
#endif
  factor=pow(-1,(double)(2*l+1))*pow(2*ci*r,(double)(l+1))*exp(-ci*r)/(cdcon1*cdcon2);
  psi1=-EULER;
  psi2=-EULER;
  for(k=1;k<=2*l+1;k++){
    psi2=psi2+1.0/(double)k;
  }
  cx=(double)(l+1)-ci*eta;
  psi3=-EULER-(1.0/cx)+cx*(pi*pi/6.0);
  for(k=1;k<100000;k++){
    delp=-cx*cx/((double)(k*k)*(cx+(double)k));
    psi3=psi3+delp;
    if(abs(delp)<1.0E-12) goto CONVERGE1;
  }
  printf("never escaped loop1 in cw2_small_r!\n");
CONVERGE1:
  lterm=log(2*ci*r);
  fact1=cdcon2/dgamma(2*l+2);
  sum1=fact1*(psi1+psi2-psi3-lterm);
  for(k=1;k<=10000;k++){
    fact1=fact1*(2*ci*r)*((double)(l+k)-ci*eta)/((double)k*(double)(2*l+1+k));
    psi1=psi1+1.0/(double)k;
    psi2=psi2+1.0/(double)(2*l+1+k);
    psi3=psi3+1.0/((double)(k-1)+cx);
    delsum1=fact1*(psi1+psi2-psi3-lterm);
    sum1=sum1+delsum1;
    if(abs(delsum1)<1.0E-15) goto CONVERGE2;
  }
  printf("never escaped loop2 in cw2_small_r!\n");
CONVERGE2:
  fact2=dgamma(2*l+1)*cdcon1/pow(-2*ci*r,2*l+1);
  sum2=fact2;
  for(k=1;k<=2*l;k++){
    fact2=fact2*((double)(k-l-1)-ci*eta)*
      (-2.0*ci*r)/((double)(k)*(double)(2*l-k+1));
    sum2=sum2+fact2;
  }
  sum1=factor*sum1;
  sum2=factor*sum2;
  answer=(sum1+sum2)*exp(pi*eta/2.0);
  return answer;
}
/* ****************************************** */
double dgamma(int mm){
  /* This calc.s gamma functions which are in the form gamma(n)
     where n is an int > 0. */
  double dg;
  int j;
  dg=1.0;
  if(mm<1) {
    for(j=1;j<=-mm+1;j++){
      dg=dg/(1.0+(double)(-j));
    }
  }
  if(mm>1){
    for(j=1;j<=mm-1;j++){
      dg=dg*((double)j);
    }
  }
  return dg;
}
