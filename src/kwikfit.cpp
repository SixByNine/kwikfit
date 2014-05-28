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
#include <cholesky.h>
#include "kwikfit.h"
#include "fftw3.h"


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
   kwikfit_template_t *result = (kwikfit_template_t*) malloc(sizeof(kwikfit_template_t));
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
		 result->profs = (kwikfit_profile_t*)realloc(result->profs,
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
		 result->profs = (kwikfit_profile_t*)realloc(result->profs,
			   sizeof(kwikfit_profile_t)*result->nprof); 
		 result->profs[iProf].comps=NULL;
		 result->profs[iProf].ncomp=0;
	  }

	  kwikfit_profile_t *prof = result->profs+iProf;
	  prof->comps = (kwikfit_component_t*)realloc(prof->comps,
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


void kwikfit_write(kwikfit_template_t* tmpl, FILE* out){
   uint64_t iComp,iProf;
   kwikfit_profile_t *prof;
   fprintf(out,"KWTEMPLATE\n");
   for (iProf=0; iProf < tmpl->nprof; iProf++){
	  prof=tmpl->profs+iProf;
	  fprintf(out,"KWPROFILE\t%s\n",prof->name);
	  for (iComp = 0; iComp < prof->ncomp; iComp++){
		 fprintf(out,"%lf %lf %lf\n",
			   prof->comps[iComp].phase,
			   prof->comps[iComp].concentration,
			   prof->comps[iComp].height);
	  }
	  fprintf(out,"# END %s\n\n",prof->name);
   }
   fprintf(out,"# END %s\n",tmpl->name);
}

double **kwikfit_designMatrix(const uint64_t nbins, kwikfit_template_t *tmpl, const double ref_phase){
   double** matrix = malloc_blas(nbins,tmpl->nprof);
   uint64_t ibin,iprof,icomp;
   double phase;
   for (ibin = 0; ibin < nbins; ibin++){
	  phase = ref_phase + (double)ibin / (double)nbins;
	  for (iprof=0; iprof < tmpl->nprof; iprof++){
		 kwikfit_profile_t *prof = tmpl->profs+iprof;
		 for (icomp=0; icomp < prof->ncomp; icomp++){
			matrix[ibin][iprof] += vonMises(phase,prof->comps+icomp);
		 }
	  }
   } 
   return matrix;
}

void kwikfit_free_designMatrix(double** matrix){
   free_blas(matrix);
}

double *kwikfit_get_cov(double *profile, uint64_t nbins){
   uint64_t i,j;
   double *cov = (double*)calloc(sizeof(double),nbins);
   double mean=0;
   for(i = 0; i < nbins; i++){
	  mean+=profile[i];
   }
   mean /= (double)nbins;
   for(j = 0; j < nbins; j++){
	  for(i = 0; i < nbins; i++){
		 cov[j] += (profile[i]-mean)*(profile[(i+j)%nbins]-mean);
	  }
   }
   for(i = 0; i < nbins; i++){
	  cov[i] /= (double)nbins;
   }

   return cov;
}

/*
 *
 * double TKleastSquares(double* b, double* white_b,
 *       double** designMatrix, double** white_designMatrix,
 *             int n,int nf, double tol, char rescale_errors,
 *                   double* outP, double* e, double** cvm){
 */
kwikfit_result_t *kwikfit_doFit(uint64_t nbins, double *profile, kwikfit_template_t *tmpl){
   uint64_t i,j;
   int64_t k;

   // first attempt to get the covariance function of the input data
   double * cov = kwikfit_get_cov(profile,nbins);
   
   double **cvm=malloc_uinv(nbins);
   
   for (i=0; i<nbins; i++){
	  for (j=0; j<nbins; j++){
		 k = (int64_t)i-(int64_t)j;
		 while (k<0)k+=nbins;
		 cvm[i][j] = cov[k];
	  }
   }


   double **uinv = malloc_uinv(nbins);
   cholesky_formUinv(uinv,cvm,nbins);
   double phase = 0.;


   free_uinv(uinv);

}
