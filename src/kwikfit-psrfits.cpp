#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2014 Michael J. Keith, University of Manchester


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <T2toolkit.h>
#include "pfits/pfits.h"
#include "mjk_cmd.h"
#include "mjklog.h"
#include "kwikfit.h"
#include "tempo2.h"

int main (int argc, char** argv){
   const char *template_fname = getS("--template","-t",argc,argv,"example.kwtemplate");
   const char *pgdev = getS("--device","-D",argc,argv,"/xs");
   const char *outfile = getS("--out","-o",argc,argv,"kwikfit.tim");

   writeResiduals=getB("--writeres","",argc,argv,0);
   debugFlag=getB("--debug","-d",argc,argv,0);
   uint64_t nitr = getI("--nitr","-i",argc,argv,3);
   uint64_t nsubbin = getI("--nsub","-s",argc,argv,4);
   logmsg("reading template file: '%s'",template_fname);
   kwikfit_template_t *tmpl = kwikfit_read_template(template_fname);

   FILE *out = fopen(outfile,"w");
   getArgs(&argc,argv);

   fprintf(out,"FORMAT 1\n");
   for(uint64_t ifile = 1; ifile < argc; ifile++){
	  logmsg("Reading '%s' as psrfits",argv[ifile]);
	  fitsfile *fp = openFitsFile(argv[ifile]);

	  dSet *hdr = initialiseDset();

	  loadPrimaryHeader(fp, hdr);
	  const uint64_t nbins=hdr->phead.nbin;
	  //int extractFoldData(fitsfile *fp,dSet *data,float dm,float *fx,float *fy,float *freq_y,float *time_y,float* bpass,int sub0);


	  int colnum;
	  int status=0;
	  if(fits_movnam_hdu(fp,BINARY_TBL,"SUBINT",1,&status)){
		 logerr("error moving to subint table");
	  }
	  float *prof=(float*)calloc(sizeof(float),nbins);
	  double *profile=(double*)calloc(sizeof(double),nbins);

	  if(fits_get_colnum(fp,CASEINSEN,"DATA",&colnum,&status)){
		 logerr("error getting data column");
	  }

	  int initflag=0;
	  if(fits_read_col_flt(fp,colnum,1,1,hdr->phead.nbin,0,prof,&initflag,&status)){
		 logerr("error reading profile");
	  }
	  if(fits_get_colnum(fp,CASEINSEN,"OFFS_SUB",&colnum,&status)){
		 logerr("error getting OFFS_SUB column");
	  }
	  double sub_offs;
	  if(fits_read_col_dbl(fp,colnum,1,1,1,0,&sub_offs,&initflag,&status)){
		 logerr("error reading OFFS_SUB");
	  }

	  char tmp[88];
	  fits_read_key(fp,TSTRING,"TBIN",tmp,NULL,&status);
	  double tbin_alt=atof(tmp);

	  double tbin;
	  if(fits_movnam_hdu(fp,BINARY_TBL,"HISTORY",1,&status)){
		 logerr("error moving to history table");
	  }
	  fits_read_key(fp,TSTRING,"NAXIS2",tmp,NULL,&status);
	  int nrow=atoi(tmp);

	  double freq=hdr->phead.freq;

	  status=0;
	  if(fits_get_colnum(fp,CASEINSEN,"TBIN",&colnum,&status)){
		 logerr("error getting tbin column");
	  }
	  if(fits_read_col_dbl(fp,colnum,nrow,1,1,tbin_alt,&tbin,&initflag,&status)){
		 logerr("error reading tbin");
	  }
	  status=0;
	  if(fits_get_colnum(fp,CASEINSEN,"CTR_FREQ",&colnum,&status)){
		 logerr("error getting ctr_Freq column");
	  }
	  if(fits_read_col_dbl(fp,colnum,nrow,1,1,freq,&freq,&initflag,&status)){
		 logerr("error reading ctr_Freq");
	  }

	  logdbg("status: %d initflag: %d",status,initflag);
	  double min = TKfindMin_f(prof,nbins);
	  for(uint64_t i=0; i < nbins; i++){
		 profile[i]=prof[i]-min;
	  }
	  free(prof);

	  kwikfit_result_t *result = kwikfit_doFit(nbins,profile,tmpl,nitr,nsubbin);
	  kwikfit_plot_result(result,pgdev);
	  printf("%s: %lf \u00b1 %lf\n",argv[ifile],result->phase,result->error);

	  double period=tbin*(double)nbins;
	  long double tstart_s = (long double)hdr->phead.smjd + (long double)hdr->phead.stt_offs;
	  long double t0 = tstart_s + (long double)sub_offs;
	  long double t_off = result->phase * period;
	  long double ToA = (long double)hdr->phead.imjd + (t0 + t_off)/86400.0;

	  logdbg("stt_offs %f" ,hdr->phead.stt_offs);
	  logdbg("imjd %d" ,hdr->phead.imjd);
	  logdbg("smjd %f" ,hdr->phead.smjd);
	  logdbg("sub_offs %lf" ,sub_offs);
	  logdbg("period %lf" ,period);



	  fprintf(out," %s %.6lf \t %.16Lf \t %f 8 -alg KW",argv[ifile],freq,ToA,result->error*period*1e6);
	  fprintf(stdout," %s %.6lf \t %.16Lf \t %f 8 -alg KW\t",argv[ifile],freq,ToA,result->error*period*1e6);
	  double ref=result->amplitudes[0];
	  for(uint64_t i=0; i<result->tmpl->nprof; i++){
		 fprintf(out," -kw%s %.2f",result->tmpl->profs[i].name,result->amplitudes[i]/ref);
	  }
	  fprintf(out,"\n");
	  fprintf(stdout,"\n");
	  fflush(out);
	  kwikfit_free_result(result);
	  freeDset(hdr);
	  closeFitsFile(fp);
   }
   return 0;
}
