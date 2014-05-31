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
   int64_t n=0;
   for (i=0; i< nbins; i++){
	  j = i-r;
	  n=0;
	  while(j < 0){
		 n+=1;
		 j+=nbins;
		 if(n>10)exit(1);
	  }
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
		 sscanf(line+11,"%s",result->name);
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
		 sscanf(line+10,"%s",result->profs[iProf].name);
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
	  phase = (double)ibin / (double)nbins - ref_phase;
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
kwikfit_result_t *kwikfit_doFit_INNER(const uint64_t nbins, double *profile, kwikfit_template_t *tmpl, const uint64_t subbin_res,double* best_profile, const uint64_t nitr,const char do_cvm, const int64_t icount, double cphase, const double cerr){
   uint64_t i,j;
   int64_t k;
   const double tol = 1.0e-27;  /* Tolerence for singular value decomposition routine */
   int64_t left=0;
   int64_t right=nbins;

   //cphase= -cphase;
   if(icount > 2){
	  double p = cphase-2*cerr;
	  while (p < 0){
		 cphase+=1.0;
		 p = cphase-2*cerr;
	  }
	  left = nbins * (cphase-2*cerr);
	  right = nbins * (cphase+2*cerr);
	  left+=nbins;
	  right+=nbins;
	  left-=1;
	  right+=1;
   }
   if(right-left >= nbins){
	  left=0;
	  right=nbins;
   }
   // TODO: consider weight the points under the pulse a little higher as the RMS is probably bigger?

   // remove the best fit template so far to get a better estimate of the covariance.
   double *outProf = (double*)calloc(sizeof(double),nbins);
   for (i=0; i<nbins; i++){
	  if(best_profile[i]<0)best_profile[i]=0;
	  outProf[i] = profile[i]-best_profile[i];
   }

   // first attempt to get the covariance function of the input data
   double * cov = kwikfit_get_cov(outProf,nbins);

   double **covMatrix=malloc_uinv(nbins);
   logmsg("cov[0] %lg",cov[0]);
   cov[0]+=1e-6;

   for (i=0; i<nbins; i++){
	  //printf("%d %lf %d COVAR\n",i,cov[i],nitr);
	  for (j=0; j<nbins; j++){
		 k = (int64_t)i-(int64_t)j;
		 while (k<0)k+=nbins;
		 if(do_cvm || i==j)
			covMatrix[i][j] = cov[k];
	  }
   }

   //double *fit_yvals = (double*) calloc(sizeof(double),nbins);
   //double*  white_yvals = (double*) calloc(sizeof(double),nbins);
   double** fit_yvals=malloc_uinv(nbins);
   double** white_yvals=malloc_uinv(nbins);

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

   uint64_t ibin,isub,index;
   double best_chisq=1e99;
   for (ibin=0; ibin<nbins*subbin_res; ibin++){
	  chisq_plot[ibin]=best_chisq;
   }
   logmsg("Creating whitened profile");
   for (uint64_t iibin=0; iibin<nbins; iibin++){
	  ibin=(iibin+left)%nbins;
	  kwikfit_rotate_array(profile,fit_yvals[ibin],nbins,-ibin);
	  TKmultMatrixVec_sq(uinv, fit_yvals[ibin],nbins,white_yvals[ibin]);
	  if(iibin > right-left)break;
   }
   logmsg("Done creating whitened profile");

   double true_phase;
   double best_phase=0;
   for (isub=0; isub<subbin_res; isub++){
	  phase = (double)isub/(double)subbin_res/(double)nbins;
	  logmsg("sb %d phase %lf",isub,phase);
	  designMatrix = kwikfit_designMatrix(nbins,tmpl,phase);
	  TKmultMatrix_sq(uinv,designMatrix,nbins,nfit,white_designMatrix);
	  for (uint64_t iibin=0; iibin<nbins; iibin++){
		 ibin=(iibin+left)%nbins;
		 true_phase=phase + (double)ibin/(double)nbins;
		 if(true_phase >= 0.5)true_phase-=1.0;
		 index = floor((1.0-true_phase)*subbin_res*nbins+0.5);
		 index = index%(subbin_res*nbins);
		 if(phase_plot[index]!=0){
			logmsg("sb %d left %d, right %d ibin %d iibin %d nbins %d index %d %d",isub,left,right,ibin,iibin,nbins,index,ibin*subbin_res+isub);
			logmsg("ERROR index=%d is reused for phase %lg",index,true_phase);
			exit(1);
		 }
		 phase_plot[index] = true_phase;
		 if(iibin > right-left)continue;

		 /*
		  * double TKleastSquares(double* b, double* white_b,
		  *       double** designMatrix, double** white_designMatrix,
		  *             int n,int nf, double tol, char rescale_errors,
		  *                   double* outP, double* e, double** cvm)
		  */

		 double chisq = TKleastSquares(fit_yvals[ibin],white_yvals[ibin],designMatrix,white_designMatrix,nbins,nfit,tol,1,outP,err,cvm);
		 chisq_plot[index] = chisq;
		 if(chisq < best_chisq){
			best_chisq=chisq;
			best_phase=true_phase;
			logmsg("Phase: %lf Chisq: %lg red-Chisq: %lg index: %ld/%ld",best_phase,chisq,chisq/(nbins-nfit-1),index,subbin_res*nbins);
			printf("  FIT:\t");
			for (i=0; i<nfit; i++){
			   printf("%lg (%lg)\t",outP[i],err[i]);
			}
			printf("\n");
			TKmultMatrixVec(designMatrix,outP,nbins,nfit,outProf);
			// rotate the profile back
			kwikfit_rotate_array(outProf,best_profile,nbins,ibin);
		 }
	  }
   }

   double polyfit[3];

   uint64_t nplot = subbin_res*nbins;
   double min=TKfindMin_d(chisq_plot,nplot);
   double max=2.5*min;
   double chisq_rot[nplot];
   double phase_rot[nplot];
   int64_t centre = ((int64_t)floor((1.0-best_phase) * (double)nplot));
   centre=centre%nplot;

   logmsg("rotate >> %d %lf %lf phs=%lf",centre,chisq_plot[(centre)],min,phase_plot[centre]);
   kwikfit_rotate_array(chisq_plot,chisq_rot,nplot,(nplot/2)-centre);
   kwikfit_rotate_array(phase_plot,phase_rot,nplot,(nplot/2)-centre);
   int64_t l=0;
   int64_t r=0;

   for(l = nplot/2; l > 0; l--){
	  if (chisq_rot[l] > max)break;
   }
   for(r=nplot/2; r <nplot; r++){
	  if (chisq_rot[r] > max)break;
   }


   l-=1;
   r+=1;
   if(l<0)l=0;
   if(r>=nplot)r=nplot-1;
   uint64_t npolyfit=r-l;
   double polyfit_x[npolyfit];
   double polyfit_y[npolyfit];

   for(i=l; i < r; i++){
	  polyfit_x[i-l] = phase_rot[i]-best_phase;
	  if (polyfit_x[i-l] < -0.5)polyfit_x[i-l]+=1;
	  if (polyfit_x[i-l] >  0.5)polyfit_x[i-l]-=1;
	  polyfit_y[i-l] = chisq_rot[i];
	  logmsg("PF %lf %lf",polyfit_x[i-l],polyfit_y[i-l]);
   }

   //   logmsg("rotate <<");
   //   kwikfit_rotate_array(chisq_plot,chisq_rot,nplot,-((nplot/2)-centre));
   //   kwikfit_rotate_array(phase_plot,phase_rot,nplot,-((nplot/2)-centre));

   logmsg("SPAN %ld ] %ld -> %ld \t %lf -> %lf",npolyfit,l,r,polyfit_x[0],polyfit_x[npolyfit-1]);
   TKfindPoly_d(polyfit_x,polyfit_y,npolyfit,3,polyfit);
   logmsg("POLYFIT %lg %lg %lg",polyfit[0],polyfit[1],polyfit[2]);

   min = polyfit[0]-polyfit[1]*polyfit[1]/4.0/polyfit[2];

   double bpp = -polyfit[1]/2.0/polyfit[2];

   logmsg("Chisq min=%lf %lf %lf",min,bpp,bpp+best_phase);
   fflush(stdout);
   double *chifit = (double*)calloc(sizeof(double),nplot);

   double icept=2*min;
   double error = (-polyfit[1] + 
		 sqrt(polyfit[1]*polyfit[1] - 4.0 * polyfit[2]*(-icept+polyfit[0])))/2.0/polyfit[2];

   if(error > 0.5){
	  error=0.5;
	  bpp=0;
   } else {
	  error-=bpp;
   }
   logmsg("err %lf",error);


   double v[3];
   for (i=0;i<nplot;i++)
   {
	  double phase = i/(double)nplot+best_phase;
	  phase = phase_rot[i]-best_phase;
	  if(phase > 0.5)phase-=1.0;
	  TKfitPoly(phase,v,3);
	  chisq_rot[i]=0;
	  for (j=0;j<3;j++)
		 chisq_rot[i] += v[j]*polyfit[j];

	  logmsg("%d %lf %lf",i,phase,chisq_rot[i]);
	  if(phase < -0.5)phase+=1.0;
   }


   kwikfit_rotate_array(chisq_rot,chifit,nplot,-((nplot/2)-centre));



   fflush(stdout);
   logdbg("Free");
   free_blas(designMatrix);
   free_blas(white_designMatrix);
   free_uinv(uinv);
   free_uinv(fit_yvals);
   free_uinv(white_yvals);
   if (nitr > 0){
	  free(outP);
	  free(err);
	  free_uinv(cvm);
	  free_uinv(covMatrix);
	  free(outProf);
	  free(phase_plot);
	  free(chisq_plot);
	  return kwikfit_doFit_INNER(nbins,profile,tmpl,2*subbin_res,best_profile,nitr-1,1,icount+1,(best_phase+bpp), error);
   } else {
	  logmsg("Create output");
	  kwikfit_result_t *result = (kwikfit_result_t*) calloc(sizeof(kwikfit_result_t),1);


	  for (i=0; i<nbins; i++){
		 outProf[i] = profile[i]-best_profile[i];
	  }
	  result->chifit = chifit;
	  result->phase = (best_phase+bpp);
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
   return kwikfit_doFit_INNER(nbins,profile,tmpl,subbin_res,best_profile,nitr,0,0,0,0);
   free(best_profile);
}

void kwikfit_free_result(kwikfit_result_t* res){
   free(res->phase_plot);
   free(res->chisq_plot);
   free(res->data);
   free(res->fit);
   free(res->residual);
   free(res->data_cov);
   free(res->amplitudes);
   free(res->amp_err);
   free_blas(res->amp_cvm);
   free(res->chifit);
   free(res);
}

