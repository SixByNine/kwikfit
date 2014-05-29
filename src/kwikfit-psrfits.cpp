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

   writeResiduals=getB("--writeres","",argc,argv,0);
   debugFlag=getB("--debug","-d",argc,argv,0);
   getArgs(&argc,argv);
   uint64_t nitr = getI("--nitr","-i",argc,argv,3);
   uint64_t nsubbin = getI("--nsub","-s",argc,argv,4);
   logmsg("reading '%s'",template_fname);
   kwikfit_template_t *tmpl = kwikfit_read_template(template_fname);
   kwikfit_write(tmpl,stdout);

   logmsg("Reading '%s' as psrfits",argv[1]);


   fitsfile *fp = openFitsFile(argv[1]);

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
   double tbin=atof(tmp);

   logmsg("status: %d initflag: %d",status,initflag);
   double min = TKfindMin_f(prof,nbins);
   for(uint64_t i=0; i < nbins; i++){
	  profile[i]=prof[i]-min;
   }
   free(prof);
   kwikfit_result_t *result = kwikfit_doFit(nbins,profile,tmpl,nitr,nsubbin);
   kwikfit_plot_result(result,pgdev);
   printf("\n\n");
   printf("%lf \u00b1 %lf\n",result->phase,result->error);
   if(result->phase < 0){
	  printf("%lf \u00b1 %lf\n",result->phase+1.0,result->error);
   }
   result->phase;

   //sub_offs -= tbin*nbins/2; // this is the start of bin 0 relative to the start
   double period=tbin*(double)nbins;
   long double tstart_s = (long double)hdr->phead.smjd + (long double)hdr->phead.stt_offs;
   long double t0 = tstart_s + (long double)sub_offs;
   long double t_off = -result->phase * period;
   long double ToA = (long double)hdr->phead.imjd + (t0 + t_off)/86400.0;

   logmsg("stt_offs %f" ,hdr->phead.stt_offs);
   logmsg("imjd %d" ,hdr->phead.imjd);
   logmsg("smjd %f" ,hdr->phead.smjd);
   logmsg("sub_offs %lf" ,sub_offs);
   logmsg("period %lf" ,period);



   printf("%s 1400 \t %.16Lf \t %f 8\n",argv[1],ToA,result->error*period*1e6);
   return 0;
}
