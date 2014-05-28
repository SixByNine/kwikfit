#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2014 Michael J. Keith, University of Manchester


#include <stdio.h>
#include <string.h>
#include "mjk_cmd.h"
#include "mjklog.h"
#include "kwikfit.h"

int main (int argc, char** argv){
   debugFlag=1;
   logmsg("TEST");
   logdbg("DEBUG");
   char *template_fname = getS("--template","-t",argc,argv,"example.kwtemplate");
   getArgs(&argc,argv);
   logmsg("reading '%s'",template_fname);
   kwikfit_template_t *tmpl = kwikfit_read_template(template_fname);
   kwikfit_write(tmpl,stdout);

   FILE *f = fopen(argv[1],"r");
   double profile[4096];
   uint64_t i=0;
   while(!feof(f)){
	  fscanf(f,"%lf",profile+i);
	  logmsg("%d %lf",i,profile[i]);
	  i++;
   }
   uint64_t nbins = i;
   kwikfit_doFit(nbins,profile,tmpl);

return 0;
}
