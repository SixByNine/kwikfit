#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2014 Michael J. Keith, University of Manchester


#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mjk_cmd.h"
#include "mjklog.h"
#include "kwikfit.h"
#include "tempo2.h"

int main (int argc, char** argv){
   const char *template_fname = getS("--template","-t",argc,argv,"example.kwtemplate");
   const char *pgdev = getS("--device","-D",argc,argv,"/xs");

   writeResiduals=getB("--writeres","",argc,argv,0);
   debugFlag=getB("--debug","-d",argc,argv,0);
   char addSin=getB("--sin","-S",argc,argv,0);
   getArgs(&argc,argv);
   uint64_t nitr = getI("--nitr","-i",argc,argv,3);
   uint64_t nsubbin = getI("--nsub","-s",argc,argv,4);
   logmsg("reading '%s'",template_fname);
   kwikfit_template_t *tmpl = kwikfit_read_template(template_fname);
   kwikfit_write(tmpl,stdout);

   FILE *f = fopen(argv[1],"r");
   double profile[4096];
   uint64_t i=0;
   while(!feof(f)){
	  fscanf(f,"%lf",profile+i);
	  if(addSin)profile[i]+=0.3*sin(M_PI*(double)i/512.0);
//	  logmsg("%d %lf",i,profile[i]);
	  i++;
   }
   uint64_t nbins = i;
   kwikfit_result_t *result = kwikfit_doFit(nbins,profile,tmpl,nitr,nsubbin);
   kwikfit_plot_result(result,pgdev);
   printf("\n\n");
   printf("%lf \u00b1 %lf\n",result->phase,result->error);
   if(result->phase < 0){
	  printf("%lf \u00b1 %lf\n",result->phase+1.0,result->error);
   }
   return 0;
}
