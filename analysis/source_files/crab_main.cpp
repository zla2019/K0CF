#ifdef COULOMB
coulcalc_class coulinfo[INTERACTION_NKMAX];
#endif

#ifdef STRONG_INTERACTION
partwave_class *partwave[STRONG_NETPARTIALS][INTERACTION_NKMAX];
#endif

phasespace_class *phasespace1_ptr[NBMAX][NPHASEMAX];
#if ! defined IDENTICAL
phasespace_class *phasespace2_ptr[NBMAX][NPHASEMAX];
#endif

int crabmain(){
  double p1[4],p2[4],r1[4],r2[4];
  double ps1[4],ps2[4];

  double corr;
  double q[4],pbar[4],qdotr,qred,r,mom;

  double phi1,phi2;
  int ib1,ib2,nb,onlyzerob;
  int j1,j2,itry,i_numden,n_numden;
  int ii,ifilter,test,iread1,iread2;
  int nmax;
  int np1[NBMAX];
  double bweight_num[NBMAX];
  double bweight1_den[NBMAX];
#ifdef IDENTICAL
  int *np2=np1;
  double *bweight2_den=bweight1_den;
//  double *bweight2_den
//double *bweight1_den;
#else
  int np2[NBMAX];
  double bweight2_den[NBMAX];
#endif
#ifdef MIXED_PAIRS_FOR_DENOM
  int n_mixed_pairs=0,ii1,ii2;
  double ww[NMAX_FOR_MIXING];
  int iiread1[NMAX_FOR_MIXING];
  int iiread2[NMAX_FOR_MIXING];
  double pp1[NMAX_FOR_MIXING][4];
  double pp2[NMAX_FOR_MIXING][4];
#endif
  printf("How many pairs?\n");
  scanf("%d",&nmax);

  phasespaceinit(&nb,np1,np2,bweight_num,bweight1_den,bweight2_den,&onlyzerob);
  bininit();
#ifdef STRONG_INTERACTION
 partwaveinit();
#endif
#ifdef COULOMB
  coulset();
#endif

#ifdef MIXED_PAIRS_FOR_DENOM
  n_numden=1;
#else
  n_numden=2;
  if(onlyzerob==1) n_numden=1;
#endif
  printf("------ Initialization Completed -------\n");

  for(i_numden=0;i_numden<n_numden;i_numden++){
    itry=0;
    for(ii=0;ii<nmax;ii++){
//printf("############## ii = %i ##################\n", ii);
    RETRY:
      itry=itry+1;
      if(i_numden==0){
	get_b_num(&ib1,&ib2,bweight_num); // one proton from impact parameter ib1 and another one from ib2 for same event (ib1=ib2)
      }
      else{
	get_b_den(&ib1,&ib2,bweight1_den,bweight2_den);
//      printf("ib1 = %i\n",ib1); // one proton from impact parameter ib1 and another one from ib2 for mixed event (ib1=ib2 ?)
//      printf("ib2 = %i\n",ib2);
      }
      j1=(int)floor((double)np1[ib1]*ran2());
      j2=(int)floor((double)np2[ib2]*ran2());
//      printf("np1[%i] = %i; np2[%i] = %i\n", ib1, np1[ib1], ib2, np2[ib2]); //number of proton for impact parameter ib1
//printf("iread1 = %i; iread2 = %i \n", iread1, iread2);
#ifdef IDENTICAL
      if(j1==j2) goto RETRY;
#endif
      phi1=2.0*pi*ran2();
      if( onlyzerob==1 || i_numden==1 ){
	phi2=2.0*pi*ran2();
      }
      else{
	phi2=phi1;
      }

      get_pr(ib1,j1,phi1,ib2,j2,phi2,p1,r1,&iread1,p2,r2,&iread2);
      if(iread1==iread2) goto RETRY;
      get_com_quantities(p1,r1,p2,r2,
			 &qred,&r,&qdotr,&mom,&test);
      if(test==0) goto RETRY;//0.5*|p1-p2| < 200 MeV/c
      twopart_filter(p1,p2,&ifilter);
      if(ifilter==0) goto RETRY;
  //    printf("j1 = %i; j2 = %i\n", j1, j2); // random proton id j1 for impact parameter ib1 and random proton id j2 for ib2
      if(i_numden==0){
	corr=corrcalc(qred,qdotr,r);
      }
      else{
	corr=1;
      }
#ifdef SMEAR_MOMENTA
      momentum_smear(p1,p2,&mom,&qred);
#endif
#ifdef MIXED_PAIRS_FOR_DENOM
      if(n_mixed_pairs<NMAX_FOR_MIXING){
	store_mixed_pair(n_mixed_pairs,p1,pp1[n_mixed_pairs],
			 p2,pp2[n_mixed_pairs],
			 iread1,iread2,iiread1,iiread2,corr,ww);
      }
      n_mixed_pairs=n_mixed_pairs+1;
#endif
      binner(i_numden,mom,corr,p1,p2);

      if((10*(ii+1))%nmax==0) printf("%d percent finished\n",
				     (int)(100.0*(double)(ii+1)/(double)nmax));
    }
    if(i_numden==0) printf("# of tries in numerator=%d\n",itry);
    if(i_numden==1) printf("# of tries in denominator=%d\n",itry);
  }

#ifdef MIXED_PAIRS_FOR_DENOM
  itry=0;
//printf("n_mixed_pairs = %i\n", n_mixed_pairs);
  if(n_mixed_pairs>NMAX_FOR_MIXING) n_mixed_pairs=NMAX_FOR_MIXING;
  i_numden=1;
  for(ii=0;ii<nmax;ii++){
  RETRY_MIXED:
    itry=itry+1;
    ii1=(int)floor((double)n_mixed_pairs*ran2());
    ii2=(int)floor((double)n_mixed_pairs*ran2());
    if(iiread1[ii1]==iiread2[ii2]) goto RETRY_MIXED;
    corr=ww[ii1]*ww[ii2];
    twopart_filter(pp1[ii1],pp2[ii2],&ifilter);
    if(ifilter==0) goto RETRY_MIXED;
    mom=get_mom(pp1[ii1],pp2[ii2]);
    binner(i_numden,mom,corr,pp1[ii1],pp2[ii2]);
  }
  printf("no. of tries in denominator=%d\n",itry);
#endif
  printstatements(onlyzerob);
}
