/* This routine reads in phase space pt.s from files and stores them for
   the purposes of calculating correlations. The file 'phasedescription.dat'
   is first read to determine which files have the info about phase space
   points. This file also includes info about impact parameters. See the
   directions to learn how to edit this file. */

/* This routine also calculates the weights used for choosing impact paramters
   when picking pairs. Weights are calculated for both the numerator and
   denominator */

/* This routine assumes info in the phasespace files is in the OSCAR format.
   As long as the files include single lines with the information,
   px,py,pz,x,y,z,t, where these correspond to the final momenta and point of
   last interaction, it should be easy to modify this routine to accomodate
   that format. */

/* If IDENTICAL is defined in the interaction, it means that the particles
   must have the same mass and be chosen according to the same ID numbers.
   CRAB will then only store the points in one array rather than two. */

#define N_BLANK_HEADER_LINES 3 /* OSCAR format has 3 useless lines */

/* define if files have mom. written in MeV/c
   rather than GeV/c (the default) */
/* #define MOM_IN_MEV */

void phasespaceinit(int *nb,int *np1,int *np2,double *bweight_num,double *bweight1_den,double *bweight2_den,int *onlyzerob){
  int ib,i,ident_check,iread=0,alpha;
  char filename[80],cdummy[401];
  double p[4],r[4];
  double b,bmin,bmax,bcut;
  double edummy,mdummy,bdummy,phidummy;
  double bwtot_num,bwtot1_den;
#if ! defined IDENTICAL
  double bwtot2_den;
#endif
  int ident,ifilter,nevents,ievent,n_testparts,idummy,npart,ipart;
  FILE *fptr,*fptr0;

  fptr0=fopen("phasedescription.dat","r");
  fscanf(fptr0,"%d",nb);
  if(*nb>NBMAX){
    printf("Increase NBMAX in crab.cpp\n)");
    exit(1);
  }

  bwtot_num=bwtot1_den=0.0;
#if ! defined IDENTICAL
  bwtot2_den=0.0;
#endif

  for(ib=0;ib<*nb;ib++){
    np1[ib]=0;
#if ! defined IDENTICAL
    np2[ib]=0;
#endif

    fscanf(fptr0,"%s %lf %lf %lf %d %d %lf",
	   filename,&b,&bmin,&bmax,
	   &nevents,&n_testparts,&bcut);
    printf("filename=%s b=%g %g %g  nevents=%d n_testparts=%d bcut=%g\n",
	   filename,b,bmin,bmax,nevents,n_testparts,bcut);

    fptr=fopen(filename,"r");
    /* Read in empty header lines */
    for(i=1;i<=N_BLANK_HEADER_LINES;i++){
      fgets(cdummy,400,fptr);
      printf("%s, LINE #%d: %s",filename,i,cdummy);
    }

    for(ievent=0;ievent<nevents;ievent++){
      /* Read in event header, only npart(# of part.s in that event) is used */
      if(feof(fptr)==1){
	printf("Trying to read past end of file %s\n",filename);
	exit(1);
      }
      fscanf(fptr,"%d %d %lf %lf",&idummy,&npart,&bdummy,&phidummy);
      printf("EVENT HEADER: i_event=%d n_particles=%d b=%g phi=%g\n",
	     idummy,npart,bdummy,phidummy);
      /* Read in phase space pt.s, counter(idummy) and mass (mdummy)
	 are not used. "ident" identifies type of particle */
      for(ipart=0;ipart<npart;ipart++){
	if(feof(fptr)==1){
	  printf("Trying to read past end of file %s, ipart=%d ievent=%d\n",
		 filename,ipart,ievent);
	  exit(1);
	}
	fscanf(fptr,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &idummy,&ident,&p[1],&p[2],&p[3],&p[0],&mdummy,
	       &r[1],&r[2],&r[3],&r[0]);
	iread=iread+1;
	ident_check=0;
//printf("iread = %i \n", iread);
	for(i=0;i<N1TYPES;i++){
	  if(ident==IDENT1[i]) ident_check=1;
	}
	if(ident_check==1){
#if ! defined MOM_IN_MEV
	  for(alpha=1;alpha<4;alpha++) p[alpha]=p[alpha]*1000.0;
#endif
	  p[0]=sqrt(MASS1*MASS1+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
	  onepart_filter(1,p,&ifilter);
	  if(ifilter==1){
	    phasespace1_ptr[ib][np1[ib]]=new phasespace_class(p,r,iread);
	    np1[ib]=np1[ib]+1;
	    if(np1[ib]>NPHASEMAX){
	      printf("Increase NPHASEMAX in crab.cpp\n");
	      exit(1);
	    }
	  }
	}
#if ! defined IDENTICAL
	ident_check=0;
	for(i=0;i<N2TYPES;i++){
	  if(ident==IDENT2[i]) ident_check=1;
	}
	if(ident_check==1){
#if ! defined MOM_IN_MEV
	  for(alpha=1;alpha<4;alpha++) p[alpha]=p[alpha]*1000.0;
#endif
	  p[0]=sqrt(MASS2*MASS2+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
	  onepart_filter(2,p,&ifilter);
	  if(ifilter==1){
	    phasespace2_ptr[ib][np2[ib]]=new phasespace_class(p,r,iread);
	    np2[ib]=np2[ib]+1;
	    if(np2[ib]>NPHASEMAX){
	      printf("Increase NPHASEMAX in crab.cpp\n");
	      exit(1);
	    }
	  }
	}
#endif
      }
    }

    fclose(fptr);
    printf("1: %d particles passed filter for b=%g\n",np1[ib],b);
#if ! defined IDENTICAL
    printf("2: %d particles passed filter for b=%g\n",np2[ib],b);
#endif

    if(*nb!=1){
      /* These are weights used for choosing impact paramters */
      bweight_num[ib]=bwtot_num
	+(bmax*bmax-bmin*bmin)*(double)(np1[ib]*np2[ib])*bcut/
	(double)(n_testparts*nevents*n_testparts*nevents);
      bwtot_num=bweight_num[ib];
      bweight1_den[ib]=bwtot1_den+(bmax*bmax-bmin*bmin)*(double)np1[ib]*bcut/
	(double)(n_testparts*nevents);
      bwtot1_den=bweight1_den[ib];
//printf("bweight_num[%i] = %f\n", ib, bweight_num[ib]);
//printf("bweight1_den[%i] = %f\n", ib, bweight1_den[ib]);
#if ! defined IDENTICAL
      bweight2_den[ib]=bwtot2_den+(bmax*bmax-bmin*bmin)*(double)np2[ib]*bcut/
	(double)(n_testparts*nevents);
      bwtot2_den=bweight2_den[ib];
#endif
    /* **************************************************** */
    }
  }

  fclose(fptr0);

  /* Now we normalize weights used for choosing impact parameters in
     both the numerator and denominator */

  *onlyzerob=0;

  if(*nb!=1){
    for(ib=0;ib<*nb;ib++){
      bweight_num[ib]=bweight_num[ib]/bwtot_num;
      bweight1_den[ib]=bweight1_den[ib]/bwtot1_den;
#if ! defined IDENTICAL
      bweight2_den[ib]=bweight2_den[ib]/bwtot2_den;
#endif
    }
  }
  else{
    bweight_num[0]=1.0;
    bweight1_den[0]=1.0;
#if ! defined IDENTICAL
    bweight2_den[0]=1.0;
#endif
    if(b<0.1) *onlyzerob=1;
  }
}
