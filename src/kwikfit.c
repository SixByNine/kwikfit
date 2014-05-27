#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2014 Michael J. Keith, University of Manchester

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <TKfit.h>
#include <TKmatrix.h>
#include "mjklog.h"
#include "kwikfit.h"

int debugFlag=0;

void trim(char * s) {
   char * p = s;
   int l = strlen(p);

   while(isspace(p[l - 1])) p[--l] = 0;
   while(* p && isspace(* p)) ++p, --l;

   memmove(s, p, l + 1);
}
kwikfit_template_t *kwikfit_read_template(char* filename){
   FILE* file = fopen(filename,"r");
   char line[1025];
   char str[1024];
   uint64_t iComp=-1;
   uint64_t iProf=-1;
   kwikfit_template_t *result = malloc(sizeof(kwikfit_template_t));
   while (!feof(file)) {
	  fgets(line, 1024, file);
	  trim(line);
	  logdbg("LINE: '%s'",line);
	  if (line[0]=='#' || line[0]=='\0')continue;
	  sscanf(line,"%s",str);
	  if (STREQ(str,"KWTEMPLATE")){
		 strcpy(result->name,line+11);
		 logdbg("KWTEMPLATE '%s'",result->name);
		 result->nprof=0;
		 result->profs=NULL;
		 continue;
	  }
	  if (STREQ(str,"KWPROFILE")){
		 result->nprof++;
		 iProf++;
		 result->profs = realloc(result->profs,
			   sizeof(kwikfit_profile_t)*result->nprof);
		 strcpy(result->profs[iProf].name,line+10);
		 logdbg("KWPROFILE '%s'",result->profs[iProf].name);
		 result->profs[iProf].comps=NULL;
		 result->profs[iProf].ncomp=0;
		 iComp=-1;
		 continue;
	  }
	  if (iProf<0){
		 result->nprof++;
		 iProf++;
		 result->profs = realloc(result->profs,
			   sizeof(kwikfit_profile_t)*result->nprof); 
		 result->profs[iProf].comps=NULL;
		 result->profs[iProf].ncomp=0;
	  }

	  kwikfit_profile_t *prof = result->profs+iProf;
	  prof->comps = realloc(prof->comps,
			sizeof(kwikfit_component_t)*prof->ncomp);
	  prof->ncomp++;
	  iComp++;
	  logdbg("iProf %lu iComp %lu",iProf,iComp);
	  sscanf(line,"%lf %lf %lf",
			&prof->comps[iComp].phase,
			&prof->comps[iComp].concentration,
			&prof->comps[iComp].height);
   }
   return result;
}


void kwikfit_write(kwikfit_template_t* template, FILE* out){
   uint64_t iComp,iProf;
   kwikfit_profile_t *prof;
   fprintf(out,"KWTEMPLATE\n");
   for (iProf=0; iProf < template->nprof; iProf++){
	  prof=template->profs+iProf;
	  fprintf(out,"KWPROFILE\t%s\n",prof->name);
	  for (iComp = 0; iComp < prof->ncomp; iComp++){
		 fprintf(out,"%lf %lf %lf\n",
			   prof->comps[iComp].phase,
			   prof->comps[iComp].concentration,
			   prof->comps[iComp].height);
	  }
	  fprintf(out,"# END %s\n\n",prof->name);
   }
   fprintf(out,"# END %s\n",template->name);
}

double **kwikfit_designMatrix(uint64_t nbins, kwikfit_template_t *tmpl){
   double** matrix = malloc_blas(nbins,tmpl->nprof);
}

void kwikfit_free_designMatrix(double** matrix){
   free_blas(matrix);
}




