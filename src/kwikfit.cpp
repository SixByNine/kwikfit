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
	  prof->ncomp++;
	  prof->comps = (kwikfit_component_t*)realloc(prof->comps,
			sizeof(kwikfit_component_t)*prof->ncomp);
	  iComp++;
	  logdbg("iProf %ld iComp %ld",iProf,iComp);
	  sscanf(line,"%lg %lg %lg",
			&(prof->comps[iComp].phase),
			&(prof->comps[iComp].concentration),
			&(prof->comps[iComp].height));
	  prof->comps[iComp].scale_free=0;
	  logdbg("phase=%lf, conc=%lf, h=%lf",prof->comps[iComp].phase,prof->comps[iComp].concentration,prof->comps[iComp].height);
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
   double** matrix = malloc_blas(nbins,tmpl->nprof+1);
   uint64_t ibin,iprof,icomp;
   double phase;
   for (ibin = 0; ibin < nbins; ibin++){
	  phase = ref_phase + (double)ibin / (double)nbins;
	  matrix[ibin][tmpl->nprof]=1;
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
kwikfit_result_t *kwikfit_doFit_INNER(const uint64_t nbins, double *profile, kwikfit_template_t *tmpl, const uint64_t subbin_res,double* best_profile, uint64_t nitr){
   uint64_t i,j;
   int64_t k;
   const double tol = 1.0e-27;  /* Tolerence for singular value decomposition routine */

   // TODO: consider weight the points under the pulse a little higher as the RMS is probably bigger?

   // remove the best fit template so far to get a better estimate of the covariance.
   double *outProf = (double*)calloc(sizeof(double),nbins);
   for (i=0; i<nbins; i++){
	  outProf[i] = profile[i]-best_profile[i];
   }

   // first attempt to get the covariance function of the input data
   double * cov = kwikfit_get_cov(outProf,nbins);

   double **covMatrix=malloc_uinv(nbins);
   cov[0]+=1e-9;

   for (i=0; i<nbins; i++){
	  //printf("%d %lf %d COVAR\n",i,cov[i],nitr);
	  for (j=0; j<nbins; j++){
		 k = (int64_t)i-(int64_t)j;
		 while (k<0)k+=nbins;
		 covMatrix[i][j] = cov[k];
	  }
   }


   double *fit_yvals = (double*) calloc(sizeof(double),nbins);
   double*  white_yvals = (double*) calloc(sizeof(double),nbins);
   double **uinv = malloc_uinv(nbins);

   double phase = 0.;
   double** designMatrix = kwikfit_designMatrix(nbins,tmpl,phase);
   const uint64_t nfit = get_blas_cols(designMatrix);
   logmsg("NFIT = %d NBINS = %d",nfit,nbins);
   double** white_designMatrix = malloc_blas(nbins,nfit);
   cholesky_formUinv(uinv,covMatrix,nbins);
   double *outP = (double*)calloc(sizeof(double),nfit);
   double *err = (double*)calloc(sizeof(double),nfit);
   double **cvm = malloc_uinv(nfit);

   double true_phase;
   double best_chisq=1e99;
   uint64_t ibin,isub;
   for (isub=0; isub<subbin_res; isub++){
	  phase = (double)isub/(double)subbin_res/(double)nbins;
	  designMatrix = kwikfit_designMatrix(nbins,tmpl,phase);
	  TKmultMatrix_sq(uinv,designMatrix,nbins,nfit,white_designMatrix);
	  for (ibin=0; ibin<nbins; ibin++){
		 true_phase=phase + (double)ibin/(double)nbins;
		 memcpy(fit_yvals,profile+ibin,sizeof(double)*(nbins-ibin));
		 if(ibin) memcpy(fit_yvals+nbins-ibin,profile,sizeof(double)*ibin);
		 TKmultMatrixVec_sq(uinv, fit_yvals,nbins,white_yvals);
		 /*
		  * double TKleastSquares(double* b, double* white_b,
		  *       double** designMatrix, double** white_designMatrix,
		  *             int n,int nf, double tol, char rescale_errors,
		  *                   double* outP, double* e, double** cvm)
		  */

		 double chisq = TKleastSquares(fit_yvals,white_yvals,designMatrix,white_designMatrix,nbins,nfit,tol,1,outP,err,cvm);

		 if(chisq < best_chisq){
			best_chisq=chisq;
			logmsg("Phase: %lf Chisq: %lg red-Chisq: %lg",true_phase,chisq,chisq/(nbins-nfit));
			printf("  FIT:\t");
			for (i=0; i<nfit; i++){
			   printf("%lg (%lg)\t",outP[i],err[i]);
			}
			printf("\n");
			TKmultMatrixVec(designMatrix,outP,nbins,nfit,outProf);
			// rotate the profile back
			i = nbins-ibin;
			memcpy(best_profile,outProf+i,sizeof(double)*(nbins-i));
			if(ibin) memcpy(best_profile+nbins-i,outProf,sizeof(double)*i);

		 }

	  }
   }


   for (i=0; i<nbins; i++) printf("%d %lg %lg %lg %lg BEST%d\n",i,
		 best_profile[i],profile[i],profile[i]-best_profile[i],cov[i],nitr); 
   fflush(stdout);
   free_blas(designMatrix);
   free_blas(white_designMatrix);
   free_uinv(uinv);
   free(fit_yvals);
   free(white_yvals);
   free(outP);
   free(outProf);
   free(err);
   free_uinv(cvm);
   if (nitr > 0){
	  kwikfit_doFit_INNER(nbins,profile,tmpl,2*subbin_res,best_profile,nitr-1);
   }
   
}

kwikfit_result_t *kwikfit_doFit(const uint64_t nbins, double *profile, kwikfit_template_t *tmpl, const uint64_t subbin_res){
   double *best_profile=(double*)calloc(sizeof(double),nbins);
   kwikfit_doFit_INNER(nbins,profile,tmpl,subbin_res,best_profile,3);
   free(best_profile);
}

