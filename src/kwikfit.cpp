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
kwikfit_result_t *kwikfit_doFit_INNER(const uint64_t nbins, double *profile, kwikfit_template_t *tmpl, const uint64_t subbin_res,double* best_profile, const uint64_t nitr,const char do_cvm, const int64_t icount, double cphase, const double cerr, double* old_phase_plot, double* old_chisq_plot){
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
	  left-=5;
	  right+=5;
   }
   if(right-left >= nbins){
	  left=0;
	  right=nbins;
   }
   logdbg("range: %ld",right-left);
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
   logdbg("cov[0] %lg",cov[0]);
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
   logdbg("NFIT = %d NBINS = %d",nfit,nbins);
   double** white_designMatrix = malloc_blas(nbins,nfit);
   cholesky_formUinv(uinv,covMatrix,nbins);
   double *outP = (double*)calloc(sizeof(double),nfit);
   double *err = (double*)calloc(sizeof(double),nfit);
   double **cvm = malloc_uinv(nfit);

   bool ignore[subbin_res*nbins];
   double *chisq_plot = (double*)malloc(sizeof(double)*subbin_res*nbins);
   double *phase_plot = (double*)malloc(sizeof(double)*subbin_res*nbins);

   uint64_t ibin,isub,index;
   double best_chisq=1e100;
   double true_phase;

   for (uint64_t ibin=0; ibin<nbins*subbin_res; ibin++)ignore[ibin]=false;
   if(old_phase_plot!=NULL){
	  for (uint64_t ibin=0; ibin<nbins*subbin_res; ibin++){
		 chisq_plot[ibin]=old_chisq_plot[ibin/2];
		 phase_plot[ibin]=old_phase_plot[ibin/2];
		 ignore[ibin] = (ibin%2 == 1);
		 if(chisq_plot[ibin]<0)ignore[ibin];
	  }
	  free(old_phase_plot);
	  free(old_chisq_plot);
   }

   for(isub=0; isub<subbin_res; isub++){
	  phase = (double)isub/(double)subbin_res/(double)nbins;
	  for (uint64_t iibin=0; iibin<nbins; iibin++){
		 ibin=(iibin+left)%nbins;
		 if(iibin > right-left)continue;
		 true_phase=phase + (double)ibin/(double)nbins;
		 if(true_phase >= 0.5)true_phase-=1.0;
		 index = floor((1.0-true_phase)*subbin_res*nbins+0.5);
		 index = index%(subbin_res*nbins);
		 chisq_plot[index]=best_chisq;
		 phase_plot[index]=0;
		 ignore[index]=false;
		 //logerr("Clearing index=%ld phase %lg (iibin=%ld ibin=%ld isub=%ld)",index,true_phase,iibin,ibin,isub);
	  }
   }

   for (uint64_t ibin=0; ibin<nbins*subbin_res; ibin++){
	  if(ignore[ibin])chisq_plot[ibin]=-1;
   }
   logdbg("Creating whitened profile");
   for (uint64_t iibin=0; iibin<nbins; iibin++){
	  ibin=(iibin+left)%nbins;
	  kwikfit_rotate_array(profile,fit_yvals[ibin],nbins,-ibin);
	  TKmultMatrixVec_sq(uinv, fit_yvals[ibin],nbins,white_yvals[ibin]);
	  if(iibin > right-left)break;
   }
   logdbg("Done creating whitened profile");

   double best_phase=0;
   for (isub=0; isub<subbin_res; isub++){
	  phase = (double)isub/(double)subbin_res/(double)nbins;
	  logdbg("sb %d phase %lf",isub,phase);
	  designMatrix = kwikfit_designMatrix(nbins,tmpl,phase);
	  TKmultMatrix_sq(uinv,designMatrix,nbins,nfit,white_designMatrix);
	  for (uint64_t iibin=0; iibin<nbins; iibin++){
		 if(iibin > right-left)break;
		 ibin=(iibin+left)%nbins;
		 true_phase=phase + (double)ibin/(double)nbins;
		 if(true_phase >= 0.5)true_phase-=1.0;
		 index = floor((1.0-true_phase)*subbin_res*nbins+0.5);
		 index = index%(subbin_res*nbins);
		 if(phase_plot[index]!=0){
			logmsg("sb %d left %d, right %d ibin %d iibin %d nbins %d index %d %d",isub,left,right,ibin,iibin,nbins,index,ibin*subbin_res+isub);
			logerr("ERROR index=%d is reused for phase %lg (iibin=%ld ibin=%ld isub=%ld)",index,true_phase,iibin,ibin,isub);
			exit(1);
		 }
		 phase_plot[index] = true_phase;

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
			logdbg("Phase: %lf Chisq: %lg red-Chisq: %lg index: %ld/%ld",best_phase,chisq,chisq/(nbins-nfit-1),index,subbin_res*nbins);
			if(debugFlag){
			   printf("  FIT:\t");
			   for (i=0; i<nfit; i++){
				  printf("%lg (%lg)\t",outP[i],err[i]);
			   }
			   printf("\n");
			}
			TKmultMatrixVec(designMatrix,outP,nbins,nfit,outProf);
			// rotate the profile back
			kwikfit_rotate_array(outProf,best_profile,nbins,ibin);
		 }
	  }
   }

   uint64_t nplot = subbin_res*nbins;
   for (i=0;i<nplot;i++) {
	  if((chisq_plot[i]) > 1e99){
		 logerr("ERROR: some chisq space not filled!");
		 logmsg("itr = %d left=%ld right=%ld",icount,left,right);
		 logmsg("%d %lf %lg",i,phase_plot[i],chisq_plot[i]);
		 exit(1);
	  }
   }
   double polyfit[3];

   double min=1e100;
   int64_t centre = 0;
   for(i=0;i<nplot;i++){
	  if(ignore[i] && chisq_plot[i] > 0){
		 logerr("Invalid chisq for ignored point");
	  }
	  if(!ignore[i] && chisq_plot[i] < 0){
		 logerr("Invalid chisq for unignored point");
	  }
   }
   for(i=0;i<nplot;i++){
	  if(chisq_plot[i] > 0 && chisq_plot[i] < min){
		 min=chisq_plot[i];
		 centre=i;
	  }
   }
   double max=2.*min;
   double chisq_rot[nplot];
   double phase_rot[nplot];

   logdbg("rotate >> %d %lf %lf phs=%lf",centre,chisq_plot[(centre)],min,phase_plot[centre]);
   kwikfit_rotate_array(chisq_plot,chisq_rot,nplot,(nplot/2)-centre);
   kwikfit_rotate_array(phase_plot,phase_rot,nplot,(nplot/2)-centre);

   int64_t l,r;
   for(l = nplot/2; l > 0; l--){
	  if(!ignore[l]){
		 logdbg("l:%ld %lf %lf",l,max,chisq_rot[l]);
		 if (chisq_rot[l] > max)break;
	  }
   }
   for(r=nplot/2; r <nplot-1; r++){
	  if(!ignore[r]){
		 logdbg("r:%ld",r);
		 if (chisq_rot[r] > max)break;
	  }
   }

   double *chifit = (double*)calloc(sizeof(double),nplot);
   double bpp=0;
   double error = (double)(r-l)/(double)nplot;
   for(i=0;i<nplot;i++)chisq_rot[i]=1e100;
   for(i=l;i<=r;i++)chisq_rot[i]=min;

   /*
	  int64_t l=nplot/2-1;
	  int64_t r=nplot/2+1;
	  for(l = nplot/2; l > 0; l--){
	  logdbg("l:%ld",l);
	  if (chisq_rot[l] > max)break;
	  }
	  for(r=nplot/2; r <nplot; r++){

	  logdbg("r:%ld",r);
	  if (chisq_rot[r] > max)break;
	  }

	  l-=1;
	  r+=2;
	  for(i=l; i < r; i++){
	  if(ignore[i]){
	  r+=1;
	  logmsg("Extend %ld -> %ld",l,r);
	  continue;
	  }
	  }
	  if(l<0)l=0;
	  if(r>=nplot)r=nplot-1;

	  uint64_t npolyfit=r-l;
	  logdbg("l/r %ld %ld N: %ld",l,r,l-r);
	  double polyfit_x[npolyfit];
	  double polyfit_y[npolyfit];

	  j=0;
	  for(i=l; i < r; i++){
	  if(ignore[i]){
	  npolyfit-=1;
	  logerr("Somehow we've lost a bin");
	  continue;
	  }
	  polyfit_x[j] = phase_rot[i]-best_phase;
	  if (polyfit_x[j] < -0.5)polyfit_x[j]+=1;
	  if (polyfit_x[j] >  0.5)polyfit_x[j]-=1;
	  polyfit_y[j] = chisq_rot[i];
	  logdbg("PF %lf %lf",polyfit_x[j],polyfit_y[j]);
	  j+=1;
	  }

	  logdbg("SPAN %ld ] %ld -> %ld \t %lf -> %lf",npolyfit,l,r,polyfit_x[0],polyfit_x[npolyfit-1]);
	  TKfindPoly_d(polyfit_x,polyfit_y,npolyfit,3,polyfit);
	  logdbg("POLYFIT %lg %lg %lg",polyfit[0],polyfit[1],polyfit[2]);

	  double bpp = -polyfit[1]/2.0/polyfit[2];

	  double v[3];
	  min = polyfit[0]-polyfit[1]*polyfit[1]/4.0/polyfit[2];
	  min=0;
	  TKfitPoly(bpp,v,3);
	  double y=0;
	  for (j=0;j<3;j++)
	  min += v[j]*polyfit[j];




	  logdbg("Chisq min=%lf %lf %lf",min,bpp,bpp+best_phase);
	  fflush(stdout);
	  double *chifit = (double*)calloc(sizeof(double),nplot);

	  double icept=2*min;
	  double error = (-polyfit[1] + 
	  sqrt(polyfit[1]*polyfit[1] - 4.0 * polyfit[2]*(-icept+polyfit[0])))/2.0/polyfit[2];

*/
   logdbg("errpre %lf",error);

   if(error > 0.5){
	  error=0.5;
	  bpp=0;
   } else {
	  error-=bpp;
   }
   logdbg("err %lf",error);
   /*
	  for (i=0;i<npolyfit;i++){
	  TKfitPoly(polyfit_x[i],v,3);
	  double y=0;
	  for (j=0;j<3;j++)
	  y += v[j]*polyfit[j];

	  logdbg("%ld %lf %lg %lg PFP%ld\n",i,polyfit_x[i],polyfit_y[i],y,icount);
	  }

	  for (i=0;i<nplot;i++) {
	  double phase = i/(double)nplot+best_phase;
	  phase = phase_rot[i]-best_phase;
	  if(phase > 0.5)phase-=1.0;
	  TKfitPoly(phase,v,3);
	  chisq_rot[i]=0;
	  for (j=0;j<3;j++)
	  chisq_rot[i] += v[j]*polyfit[j];
	  if(phase < -0.5)phase+=1.0;
	  }
	  */

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
	  uint64_t new_subbin_res=subbin_res;
	  if(icount > 1)new_subbin_res*=2;
	  return kwikfit_doFit_INNER(nbins,profile,tmpl,new_subbin_res,best_profile,nitr-1,1,icount+1,(best_phase+bpp), error,phase_plot,chisq_plot);
   } else {
	  logdbg("Create output");
	  kwikfit_result_t *result = (kwikfit_result_t*) calloc(sizeof(kwikfit_result_t),1);

	  double prev=chisq_plot[nplot-1];
	  for (i=0; i<nplot; i++){
		 if(ignore[i])chisq_plot[i]=prev;
		 else prev=chisq_plot[i]; 
	  }
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
   return kwikfit_doFit_INNER(nbins,profile,tmpl,subbin_res,best_profile,nitr,0,0,0,0,NULL,NULL);
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

