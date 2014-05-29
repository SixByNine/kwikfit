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
#include <T2toolkit.h>
#include <cholesky.h>
#include "kwikfit.h"
#include "fftw3.h"


void kwikfit_rotate_array(double* in, double* out, const int64_t nbins, const int64_t r){
   int64_t i;
   int64_t j;
   for (i=0; i< nbins; i++){
	  j = i-r;
	  while(j < 0)j+=nbins;
	  j = j%nbins;
	  out[i] = in[j];
   }
}

void trim(char * s) {
   char * p = s;
   int l = strlen(p);

   while(isspace(p[l - 1])) p[--l] = 0;
   while(* p && isspace(* p)) ++p, --l;

   memmove(s, p, l + 1);
}
kwikfit_template_t *kwikfit_read_template(const char* filename){
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
kwikfit_result_t *kwikfit_doFit_INNER(const uint64_t nbins, double *profile, kwikfit_template_t *tmpl, const uint64_t subbin_res,double* best_profile, uint64_t nitr,char do_cvm){
   uint64_t i,j;
   int64_t k;
   const double tol = 1.0e-27;  /* Tolerence for singular value decomposition routine */

   // TODO: consider weight the points under the pulse a little higher as the RMS is probably bigger?

   if(false){
   double delta = (profile[nbins-1]-profile[0])/(double)nbins;
   // remove any jump at the edge of the profile
   for (i=0; i<nbins; i++){
	  profile[i]-=delta*i;
   }
   }

   // remove the best fit template so far to get a better estimate of the covariance.
   double *outProf = (double*)calloc(sizeof(double),nbins);
   for (i=0; i<nbins; i++){
	  outProf[i] = profile[i]-best_profile[i];
   }

   // first attempt to get the covariance function of the input data
   double * cov = kwikfit_get_cov(outProf,nbins);

   double **covMatrix=malloc_uinv(nbins);
   logmsg("cov[0] %lg",cov[0]);
   cov[0]+=1e-9;

   for (i=0; i<nbins; i++){
	  //printf("%d %lf %d COVAR\n",i,cov[i],nitr);
	  for (j=0; j<nbins; j++){
		 k = (int64_t)i-(int64_t)j;
		 while (k<0)k+=nbins;
		 if(do_cvm || i==j)
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
   double *chisq_plot = (double*)calloc(sizeof(double),subbin_res*nbins);
   double *phase_plot = (double*)calloc(sizeof(double),subbin_res*nbins);

   double true_phase;
   double best_phase=0;
   double best_chisq=1e99;
   uint64_t ibin,isub,index;
   for (isub=0; isub<subbin_res; isub++){
	  phase = (double)isub/(double)subbin_res/(double)nbins;
	  designMatrix = kwikfit_designMatrix(nbins,tmpl,phase);
	  TKmultMatrix_sq(uinv,designMatrix,nbins,nfit,white_designMatrix);
	  for (ibin=0; ibin<nbins; ibin++){
		 true_phase=phase + (double)ibin/(double)nbins;
		 if(true_phase >= 0.5)true_phase-=1.0;
		 index = floor((true_phase+0.5)*subbin_res*nbins+0.5);
		 if(phase_plot[index]!=0){
			logmsg("ERROR %d %lg",index,true_phase);
			exit(1);
		 }

		 kwikfit_rotate_array(profile,fit_yvals,nbins,ibin);
		 TKmultMatrixVec_sq(uinv, fit_yvals,nbins,white_yvals);
		 /*
		  * double TKleastSquares(double* b, double* white_b,
		  *       double** designMatrix, double** white_designMatrix,
		  *             int n,int nf, double tol, char rescale_errors,
		  *                   double* outP, double* e, double** cvm)
		  */

		 double chisq = TKleastSquares(fit_yvals,white_yvals,designMatrix,white_designMatrix,nbins,nfit,tol,1,outP,err,cvm);
		 chisq_plot[index] = chisq;
		 phase_plot[index] = true_phase;
		 if(chisq < best_chisq){
			best_chisq=chisq;
			best_phase=true_phase;
			logmsg("Phase: %lf Chisq: %lg red-Chisq: %lg",true_phase,chisq,chisq/(nbins-nfit-1));
			printf("  FIT:\t");
			for (i=0; i<nfit; i++){
			   printf("%lg (%lg)\t",outP[i],err[i]);
			}
			printf("\n");
			TKmultMatrixVec(designMatrix,outP,nbins,nfit,outProf);
			// rotate the profile back
			kwikfit_rotate_array(outProf,best_profile,nbins,-ibin);
		 }
	  }
   }


   fflush(stdout);
   logdbg("Free");
   free_blas(designMatrix);
   free_blas(white_designMatrix);
   free_uinv(uinv);
   free(fit_yvals);
   free(white_yvals);
   if (nitr > 0){
	  free(outP);
	  free(err);
	  free_uinv(cvm);
	  free_uinv(covMatrix);
	  free(outProf);
	  free(phase_plot);
	  free(chisq_plot);
	  return kwikfit_doFit_INNER(nbins,profile,tmpl,2*subbin_res,best_profile,nitr-1,1);
   } else {
	  logmsg("Create output");
	  kwikfit_result_t *result = (kwikfit_result_t*) calloc(sizeof(kwikfit_result_t),1);
	  for (i=0; i<nbins; i++){
		 outProf[i] = profile[i]-best_profile[i];
	  }
	  double polyfit[3];
	  
	  uint64_t nplot = subbin_res*nbins;
	  double min=TKfindMin_d(chisq_plot,nplot);
	  double max=2.*min;
	  double chisq_rot[nplot];
	  const uint64_t centre = ((int64_t)floor(best_phase * nplot + 0.5));

	  logmsg("rotate >>");
	  kwikfit_rotate_array(chisq_plot,chisq_rot,nplot,-centre);
	  int64_t l=0;
	  int64_t r=0;

	  for(l = nplot/2; l > 0; l--){
		 if (chisq_rot[l] > max)break;
	  }
	  for(r=nplot/2; r <nplot; r++){
		 if (chisq_rot[r] > max)break;
	  }
	  uint64_t npolyfit=r-l;
	  double polyfit_x[npolyfit];
	  double polyfit_y[npolyfit];

	  for(i=l; i < r; i++){
		 polyfit_x[i-l] = phase_plot[i];
		 polyfit_y[i-l] = chisq_rot[i];
	  }

	  logmsg("SPAN %ld ] %ld -> %ld \t %lf -> %lf",npolyfit,l,r,polyfit_x[0],polyfit_x[npolyfit-1]);
	  TKfindPoly_d(polyfit_x,polyfit_y,npolyfit,3,polyfit);
	  logmsg("POLYFIT %lg %lg %lg",polyfit[0],polyfit[1],polyfit[2]);

	  min = polyfit[0]-polyfit[1]*polyfit[1]/4.0/polyfit[2];
	  double bpp = -polyfit[1]/2.0/polyfit[2];
	  logmsg("Chisq min=%lf %lf",min,bpp+best_phase);
	  fflush(stdout);
	  double *chifit = (double*)calloc(sizeof(double),nplot);

	  double icept=2*min;
	  double error = (-polyfit[1] + 
		 sqrt(polyfit[1]*polyfit[1] - 4.0 * polyfit[2]*(-icept+polyfit[0])))/2.0/polyfit[2];

	  logmsg("err %lf",error);
	  error-=bpp;
	  logmsg("err %lf",error);


	  double v[3];
	  for (i=0;i<nplot;i++)
	  {
		 TKfitPoly(phase_plot[i],v,3);
		 chisq_rot[i]=0;
		 for (j=0;j<3;j++)
			chisq_rot[i] += v[j]*polyfit[j];
	  }

	  logmsg("rotate <<");
	  kwikfit_rotate_array(chisq_rot,chifit,nplot,centre);

	  result->chifit = chifit;
	  result->phase = best_phase+bpp;
	  result->error = error;
	  result->chisq = best_chisq;
	  result->nfree = nbins-nfit-1;
	  result->nbins = nbins;
	  result->data = profile;
	  result->fit = best_profile;
	  result->residual = outProf;
	  result->data_cov = cov;
	  result->tmpl = tmpl;
	  result->amplitudes = outP;
	  result->amp_cvm = cvm;
	  result->amp_err = err;
	  result->phase_plot = phase_plot;
	  result->chisq_plot = chisq_plot;
	  result->nplot = nplot;
	  return result;
   }
}

kwikfit_result_t *kwikfit_doFit(const uint64_t nbins, double *profile, kwikfit_template_t *tmpl, const uint64_t nitr, const uint64_t subbin_res){
   double *best_profile=(double*)calloc(sizeof(double),nbins);
   return kwikfit_doFit_INNER(nbins,profile,tmpl,subbin_res,best_profile,nitr,0);
   free(best_profile);
}

