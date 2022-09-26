#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2() {
  long J;
  long K;
  static long IDUM2=123456789;
  static long IY=0;
  static long IV[NTAB];
  double randy;
  if(IDUM <=0){
    if(-(IDUM)<1) IDUM=1;
    else IDUM = -(IDUM);
    IDUM2=(IDUM);
    for (J=NTAB+7;J>=0;J--) {
      K=(IDUM)/IQ1;
      IDUM=IA1*(IDUM-K*IQ1)-K*IR1;
      if(IDUM<0) IDUM +=IM1;
      if(J<NTAB) IV[J]=IDUM;
    }
    IY=IV[0];
    //printf("Initialization Completed %d\n",IDUM);
  }
  K=(IDUM)/IQ1;
  IDUM=IA1*(IDUM-K*IQ1)-K*IR1;
  if(IDUM<0) IDUM+=IM1;
  K=IDUM2/IQ2;
  IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2;
  if(IDUM2<0) IDUM2+=IM2;
  J=IY/NDIV;
  IY=IV[J]-IDUM2;
  IV[J]=IDUM;
  if(IY<1) IY+=IMM1;
  randy=AM*IY;
  if(randy>RNMX) return RNMX;
  else return randy;
}
//**********************************************************
void gauss(double *g1,double *g2){
  double x,y,r2,r;
  do{
    x=1.0-2.0*ran2();
    y=1.0-2.0*ran2();
    r2=x*x+y*y;
  }while (r2>1.0);
  r=sqrt(r2);
  *g1=(x/r)*sqrt(-2.0*log(r2));
  *g2=(y/x)**g1;
}
void gauss(double *g){
  double x,y,r2,r;
  do{
    x=1.0-2.0*ran2();
    y=1.0-2.0*ran2();
    r2=x*x+y*y;
  }while (r2>1.0);
  r=sqrt(r2);
  *g=(x/r)*sqrt(-2.0*log(r2));
}
