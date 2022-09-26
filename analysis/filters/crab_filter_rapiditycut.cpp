void onepart_filter(int whichpart,double *p,int *ifilter){
  /* This should NOT depend on phi!
     If identical particles, this is only called with whichpart==1,
     If not identical, this gets called with whichpart==1 for the first
     particle and whichpart==2 for the second. */
  double rapidity;
  rapidity=0.5*log((p[0]+p[3])/(p[0]-p[3]));
  if(fabs(rapidity-RAP0)<0.25){
    *ifilter=1;
  }
  else{
    *ifilter=0;
  }
}

/* ********************************************* */

void twopart_filter(double *p1,double *p2,int *ifilter){
  *ifilter=1;
}

/* ********************************************* */

