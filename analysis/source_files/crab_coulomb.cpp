void coulset(){
  double qred,eta,rmass;
  int i0,kk,nkmax;
  double_complex c1phase,c2phase,cgamma;
  double_complex couly,cfact1,cfact2,chype[41],chype1[6],chype2[6];
  rmass=MASS1*MASS2/(MASS1+MASS2);
  for(kk=0;kk<INTERACTION_NKMAX;kk++){
#ifdef REDUCED_MOM
    qred=(0.5+kk)*INTERACTION_DELK;
#else
    qred=(0.5+kk)*INTERACTION_DELK*0.5;
#endif
    eta=(double)Q1Q2*(rmass/137.036)/qred;
    coulpar(eta,&couly,&cfact1,&cfact2,chype,chype1,chype2);
    coulinfo[kk].couly=couly;
    coulinfo[kk].cfact1=cfact1;
    coulinfo[kk].cfact2=cfact2;
    for(i0=0;i0<40;i0++){
      coulinfo[kk].chype[i0]=chype[i0];
    }
    for(i0=0;i0<6;i0++){
      coulinfo[kk].chype1[i0]=chype1[i0];
      coulinfo[kk].chype2[i0]=chype2[i0];
    }
  }
}

/*  ****************************************** */

void coulpar(double eta, double_complex *couly,double_complex *cfact1,double_complex *cfact2,double_complex *chype, double_complex *chype1, double_complex *chype2){
  double_complex a,b,a0,b0;
  double_complex c1top1,c2top2,bot;
  double_complex c1top2,c2top1;
  int j;
  /*See appendix of Messiah.  "cgamma" is the gamma function. "chyper" is
    appropriate hyperbolic function.  This is explained in Messiah's
    section on coul wave func.s. Notation should be explanatory when
    reading the book. */
  *couly=cgamma(1.0+ci*eta);
  *couly=*couly*exp(-0.5*eta*pi);
  a=-ci*eta;
  b=1.0;
  a0=a;
  b0=b;
  for(j=1;j<=40;j++){
    chype[j]=a/(b*(double)j);
    a=a+1.0;
    b=b+1.0;
  }
  c1top1=1.0;
  c2top1=1.0;
  c1top2=1.0;
  c2top2=1.0;
  bot=1.0;
  for(j=1;j<=5;j++){
    c1top1=c1top1*((double)j+a0-1.0);
    c2top1=c2top1*((double)j-a0);
    c1top2=c1top2*((double)j-b0+a0);
    c2top2=c2top2*((double)j+b0-a0-1.0);
    bot=bot*(double)j;
    chype1[j]=(c1top1*c1top2)/(bot);
    chype2[j]=(c2top1*c2top2)/(bot);
  }
  *cfact1=1.0/(cgamma(b0-a0));
  *cfact2=1.0/(cgamma(a0));
}

double_complex coul(double qred,double r,double zk,double eta,int kk){
  double_complex bb,coulguy;
  double arg;
  /* See appendix of Messiah.  "cgamma" is the gamma function. "chyper" is
     appropriate hyperbolic function.  This is explained in Messiah's
     section on coul wave func.s. Notation should be explanatory
     when reading the book. */
  bb=1.0;
  coulguy=coulinfo[kk].couly*chyper(-ci*eta,bb,ci*qred*(r-zk)/197.323,kk);
  arg=zk*qred/197.323;
  arg=arg-2.0*pi*floor(arg/(2.0*pi));
  coulguy=coulguy*exp(ci*arg);
  return coulguy;
}

/* **************************************************** */

double_complex chyper(double_complex a,double_complex b,double_complex cz,
		      int kk){
  double_complex cw1,cw2,cf1,cf2,czarg,czstarj,answer;
#define rcrit 10.0
  double dmag;
  int j;
  dmag=abs(cz);
  if(dmag<rcrit){
    cf1=1.0;
    czstarj=1.0;
    for(j=1;j<=40;j++){
      czstarj=czstarj*cz*coulinfo[kk].chype[j];
      cf1=cf1+czstarj;
      if(abs(czstarj)<0.0001) goto GOOD_ENOUGH;
    }
    printf("chyper not coverging!.\n");
  GOOD_ENOUGH:
    answer=cf1;
  }
  else{
    cf1=1.0;
    cf2=1.0;
    for(j=1;j<=5;j++){
      cf1=cf1+coulinfo[kk].chype1[j]/(pow(-cz,j));
      cf2=cf2+coulinfo[kk].chype2[j]/(pow(cz,j));
    }
    cw1=cf1*=coulinfo[kk].cfact1*(pow(-cz,-a));
    czarg=cz-ci*2.0*pi*floor(real(-ci*cz/(2.0*pi)));
    cw2=cf2*coulinfo[kk].cfact2*pow(cz,a-b)*exp(czarg);
    answer=cw1+cw2;
  }
  return answer;
}

/* ***************************************** */

double_complex cgamma(double_complex c){
  /* This calc.s gamma functions which are in the form gamma(n+i*y)
     where n is an int and y is real. */
  double_complex cg,cphase;
  int mm,j;
  double x,y,phase,delp,cgmag;
  x=real(c);
  y=imag(c);
  phase=-EULER*y;
  for(j=1;j<=100000;j++){
    delp=(y/(double)j)-atan(y/(double)j);
    phase=phase+delp;
    if(fabs(delp)<1E-10) goto CGAMMA_ESCAPE;
  }
  printf("oops not accurate enough, increase jmax\n");
CGAMMA_ESCAPE:
  phase=phase-2.0*pi*floor(phase/(2.0*pi));
  cphase=exp(ci*phase);
  cgmag=sqrt(pi*y/sinh(pi*y));
  mm=(int)floor(x+0.5);
  cg=cgmag*cphase;
  if(mm<1){
    for(j=1;j<=-mm+1;j++){
      cg=cg/(1.0+(double)(-j)+ci*y);
    }
  }
  if(mm>1) {
    for(j=1;j<=mm-1;j++){
      cg=cg*((double)(j)+ci*y);
    }
  }
  return cg;
}
