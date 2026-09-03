
#include "glip_features.h"


namespace glip{

    namespace Features{


      sFeatures *Create(int n, int nfeats, int style){
	sFeatures *f = NULL;
	int i;
	f = (sFeatures *) calloc(1,sizeof(sFeatures));
	if(f == NULL)
	  glip::Error((char *)MSG1,(char *)"Features::Create");

	f->n = n;
	f->nfeats = nfeats;

	f->label = NULL;
	if(style & glipFeature_ALLOCLABEL)
	  f->label = glip::AllocIntArray(n);

	f->index = NULL;
	if(style & glipFeature_ALLOCINDEX)
	  f->index = glip::AllocIntArray(n);

	f->fv = NULL;
	if(style & glipFeature_ALLOCFV){
	  f->fv = (float **)calloc(n, sizeof(float *));
	  if(f->fv == NULL)
	    glip::Error((char *)MSG1,(char *)"Features::Create");
	  for(i = 0; i < n; i++)
	    f->fv[i] = glip::AllocFloatArray(nfeats);
	}

	f->dist = NULL;
	if(style & glipFeature_ALLOCDIST){
	  f->dist = (float **)calloc(n, sizeof(float *));
	  if(f->dist == NULL)
	    glip::Error((char *)MSG1,(char *)"Features::Create");
	  for(i = 0; i < n; i++)
	    f->dist[i] = glip::AllocFloatArray(n);
	}
	
	return f;
      }

      
      sFeatures *Clone(sFeatures *f){
	sFeatures *c;
	int i, style = 0;
	if(f->label != NULL) style = style | glipFeature_ALLOCLABEL;
	if(f->index != NULL) style = style | glipFeature_ALLOCINDEX;
	if(f->fv != NULL)    style = style | glipFeature_ALLOCFV;
	if(f->dist != NULL)  style = style | glipFeature_ALLOCDIST;
	c = Create(f->n, f->nfeats, style);

	if(style & glipFeature_ALLOCLABEL)
	  memcpy(c->label, f->label, f->n*sizeof(int));
	if(style & glipFeature_ALLOCINDEX)	
	  memcpy(c->index, f->index, f->n*sizeof(int));

	for(i = 0; i < f->n; i++){
	  if(style & glipFeature_ALLOCFV)
	    memcpy(c->fv[i], f->fv[i], f->nfeats*sizeof(float));
	  if(style & glipFeature_ALLOCDIST)
	    memcpy(c->dist[i], f->dist[i], f->n*sizeof(float));
	}
	return c;
      }


      void     Destroy(sFeatures **f){
	sFeatures *tmp;
	int i;
	tmp = *f;
	if(tmp != NULL){
	  if(tmp->label != NULL)
	    glip::FreeIntArray(&tmp->label);
	  if(tmp->index != NULL)
	    glip::FreeIntArray(&tmp->index);
	  if(tmp->fv != NULL){
	    for(i = 0; i < tmp->n; i++)
	      glip::FreeFloatArray(&tmp->fv[i]);
	    free(tmp->fv);
	  }
	  if(tmp->dist != NULL){
	    for(i = 0; i < tmp->n; i++)
	      glip::FreeFloatArray(&tmp->dist[i]);
	    free(tmp->dist);
	  }
	  free(tmp);
	}
	*f = NULL;
      }



      void Randomize(sFeatures *f){
	int i,j,k,tmp;
	float ftmp;
	float *ptmp;
	
	//srand((int)time(NULL));
	for(i = 0; i < f->n; i++){
	  j = glip::RandomInteger(0, f->n-1);

	  if(f->label != NULL){
	    tmp = f->label[i];
	    f->label[i] = f->label[j];
	    f->label[j] = tmp;
	  }

	  if(f->index != NULL){
	    tmp = f->index[i];
	    f->index[i] = f->index[j];
	    f->index[j] = tmp;
	  }

	  if(f->fv != NULL){
	    ptmp = f->fv[i];
	    f->fv[i] = f->fv[j];
	    f->fv[j] = ptmp;
	  }

	  if(f->dist != NULL){
	    ptmp = f->dist[i];
	    f->dist[i] = f->dist[j];
	    f->dist[j] = ptmp;

	    for(k = 0; k < f->n; k++){
	      ftmp = f->dist[k][i];
	      f->dist[k][i] = f->dist[k][j];
	      f->dist[k][j] = ftmp;
	    }
	  }
	  
	}
      }
      
      
      
      sFeatures *RemoveSamples(sFeatures **f,
			       float rate){
	sFeatures *f0 = NULL,*f1 = NULL,*f2 = NULL;
	int *hist,*hist1,*hist2,*index,*pos;
	int Lmax,i,j,l,nsamples,nfeats,p1,p2,style = 0;

	f0 = *f;
	
	if(f0->label != NULL) style = style | glipFeature_ALLOCLABEL;
	if(f0->index != NULL) style = style | glipFeature_ALLOCINDEX;
	if(f0->fv != NULL   ) style = style | glipFeature_ALLOCFV;
	if(f0->dist != NULL ) style = style | glipFeature_ALLOCDIST;

	Lmax = 0;
	for(i = 0; i < f0->n; i++)
	  if(f0->label[i] > Lmax)
	    Lmax = f0->label[i];
	hist  = glip::AllocIntArray(Lmax+1);
	hist1 = glip::AllocIntArray(Lmax+1);
	hist2 = glip::AllocIntArray(Lmax+1);
	pos   = glip::AllocIntArray(Lmax+1);
	
	for(i = 0; i < f0->n; i++)
	  hist[f0->label[i]]++;

	nsamples = 0;
	for(i = 0; i <= Lmax; i++){
	  hist1[i] = ceil(hist[i]*rate);
	  hist2[i] = hist[i] - hist1[i];
	  nsamples += hist1[i];
	}

	pos[0] = 0;
	for(i = 1; i <= Lmax; i++)
	  pos[i] = pos[i-1]+hist[i-1];

	nfeats = f0->nfeats;
	f1 = Create(nsamples, nfeats, style);
       	f2 = Create(f0->n-nsamples, nfeats, style);

	//Index to positions ordered by label.
	index = glip::AllocIntArray(f0->n);
	for(i = 0; i < f0->n; i++){
	  l = f0->label[i];
	  j = pos[l];
	  index[j] = i;
	  pos[l]++;
	}
	
	p1 = p2 = 0;
	for(l = 0; l <= Lmax; l++){

	  for(i = 0; i < hist1[l]; i++){
	    if(f0->fv != NULL)
	      memcpy(f1->fv[p1],
		     f0->fv[index[p1+p2]],
		     nfeats*sizeof(float));
	    if(f0->label != NULL)
	      f1->label[p1] = f0->label[index[p1+p2]];
	    if(f0->index != NULL)
	      f1->index[p1] = f0->index[index[p1+p2]];
	    p1++;
	  }

	  for(j = 0; j < hist2[l]; j++){
	    if(f0->fv != NULL)
	      memcpy(f2->fv[p2],
		     f0->fv[index[p1+p2]],
		     nfeats*sizeof(float));
	    if(f0->label != NULL)
	      f2->label[p2] = f0->label[index[p1+p2]];
	    if(f0->index != NULL)
	      f2->index[p2] = f0->index[index[p1+p2]];
	    p2++;
	  }
	}
	
	//Changes lf e lf2.
	*f = f2;
	Destroy(&f0);
	
	glip::FreeIntArray(&index);
	glip::FreeIntArray(&hist);
	glip::FreeIntArray(&hist1);
	glip::FreeIntArray(&hist2);
	glip::FreeIntArray(&pos);
  
	return f1;
      }



      void Resize(sFeatures **f, int n){
	sFeatures *f0 = NULL;
	sFeatures *c;
	int i, style = 0;
	f0 = *f;
	if(f0->label != NULL) style = style | glipFeature_ALLOCLABEL;
	if(f0->index != NULL) style = style | glipFeature_ALLOCINDEX;
	if(f0->fv != NULL)    style = style | glipFeature_ALLOCFV;
	if(f0->dist != NULL)  style = style | glipFeature_ALLOCDIST;
	c = Create(n, f0->nfeats, style);
	if(style & glipFeature_ALLOCLABEL)
	  memcpy(c->label, f0->label, MIN(n,f0->n)*sizeof(int));

	if(style & glipFeature_ALLOCINDEX)
	  memcpy(c->index, f0->index, MIN(n,f0->n)*sizeof(int));
	
	for(i = 0; i < MIN(n,f0->n); i++){
	  if(style & glipFeature_ALLOCFV)
	    memcpy(c->fv[i], f0->fv[i], f0->nfeats*sizeof(float));
	  if(style & glipFeature_ALLOCDIST)
	    memcpy(c->dist[i], f0->dist[i], MIN(n,f0->n)*sizeof(float));
	}
	*f = c;
	Destroy(&f0);
      }
      
      
      sFeatures *Merge(sFeatures *f1, 
		       sFeatures *f2){
	sFeatures *f;
	int i,n,nfeats,style = 0;
	
	if(f1==NULL && f2==NULL) return NULL;
	if(f1==NULL)
	  return Clone(f2);
	if(f2==NULL)
	  return Clone(f1);
	
	n = f1->n + f2->n;
	nfeats = f1->nfeats;
	if(nfeats != f2->nfeats) return NULL;

	if(f1->label != NULL || f2->label != NULL) style = style | glipFeature_ALLOCLABEL;
	if(f1->index != NULL || f2->index != NULL) style = style | glipFeature_ALLOCINDEX;
	if(f1->fv != NULL    || f2->fv != NULL)    style = style | glipFeature_ALLOCFV;
	if(f1->dist != NULL  || f2->dist != NULL)  style = style | glipFeature_ALLOCDIST;
	
	f = Create(n, nfeats, style);

	if(f1->label != NULL)
	  memcpy(f->label,           f1->label, f1->n*sizeof(int));
	if(f2->label != NULL)
	  memcpy(&(f->label[f1->n]), f2->label, f2->n*sizeof(int));

	if(f1->index != NULL)
	  memcpy(f->index,           f1->index, f1->n*sizeof(int));
	if(f2->index != NULL)
	  memcpy(&(f->index[f1->n]), f2->index, f2->n*sizeof(int));

	if(f1->fv != NULL)
	  for(i = 0; i < f1->n; i++)
	    memcpy(f->fv[i], f1->fv[i], nfeats*sizeof(float));

	if(f2->fv != NULL)
	  for(i = 0; i < f2->n; i++)
	    memcpy(f->fv[i+f1->n], f2->fv[i], nfeats*sizeof(float));

	if(f1->dist != NULL)
	  for(i = 0; i < f1->n; i++)
	    memcpy(f->dist[i], f1->dist[i], f1->n*sizeof(float));

	if(f2->dist != NULL)
	  for(i = 0; i < f2->n; i++)
	    memcpy(&(f->dist[i+f1->n][f1->n]), f2->dist[i], f2->n*sizeof(float));
	
	return f;
      }
      


      int *KNN(sFeatures *f,
	       float *fv,
	       int k){
	int   *best = (int *)calloc(k+1, sizeof(int));
	float *dist = (float *)calloc(k+1, sizeof(float));
	float d;
	int i,j;
	for(i = 0; i < k; i++)
	  dist[i] = FLT_MAX;

	for(i = 0; i < f->n; i++){
	  d = 0.0;
	  for(j = 0; j < f->nfeats; j++)
	    d += SQUARE(fv[j] - f->fv[i][j]);

	  j = k;
	  while(j > 0 && dist[j-1] > d){
	    dist[j] = dist[j-1];
	    best[j] = best[j-1];
	    j--;
	  }
	  best[j] = i;
	  dist[j] = d;
	}
	free(dist);
	return best;
      }



      float KNNMeanDistance(sFeatures *f,
			    float *fv,
			    int k){
	int   *best = (int *)calloc(k+1, sizeof(int));
	float *dist = (float *)calloc(k+1, sizeof(float));
	float d,sum;
	int i,j;
	for(i = 0; i < k; i++)
	  dist[i] = FLT_MAX;

	for(i = 0; i < f->n; i++){
	  d = 0.0;
	  for(j = 0; j < f->nfeats; j++)
	    d += SQUARE(fv[j] - f->fv[i][j]);

	  j = k;
	  while(j > 0 && dist[j-1] > d){
	    dist[j] = dist[j-1];
	    best[j] = best[j-1];
	    j--;
	  }
	  best[j] = i;
	  dist[j] = d;
	}

	sum = 0.0;
	for(i = 0; i < k; i++){
	  sum += sqrtf(dist[i]);
	}
	free(dist);
	free(best);
	return (sum/k);
      }
      


      sFeatures *GetSamples(sImage32 *img,
			    int *S){
	sFeatures *f;
	int i,p;
	f = Create(S[0], 1, glipFeature_ALLOCINDEX | glipFeature_ALLOCFV);
	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  f->index[i-1] = p;
	  f->fv[i-1][0] = img->data[p];
	}
	return f;
      }

      
      sFeatures *GetSamples(sCImage *cimg,
			    int *S){
	sFeatures *f;
	int i,p;
	f = Create(S[0], 3, glipFeature_ALLOCINDEX | glipFeature_ALLOCFV);
	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  f->index[i-1] = p;
	  f->fv[i-1][0] = cimg->C[0]->data[p];
	  f->fv[i-1][1] = cimg->C[1]->data[p];
	  f->fv[i-1][2] = cimg->C[2]->data[p];
	}
	return f;
      }


      sFeatures *GetSamples(sCImage32f *cimg,
			    int *S){
	sFeatures *f;
	int i,p;
	f = Create(S[0], 3, glipFeature_ALLOCINDEX | glipFeature_ALLOCFV);
	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  f->index[i-1] = p;
	  f->fv[i-1][0] = cimg->C[0]->data[p];
	  f->fv[i-1][1] = cimg->C[1]->data[p];
	  f->fv[i-1][2] = cimg->C[2]->data[p];
	}
	return f;
      }



      sFeatures *GetSamples(sScene32 *scn,
			    int *S){
	sFeatures *f;
	int i,p;
	f = Create(S[0], 1, glipFeature_ALLOCINDEX | glipFeature_ALLOCFV);
	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  f->index[i-1] = p;
	  f->fv[i-1][0] = scn->data[p];
	}
	return f;
      }
      
      

      sImage32 *KNNFuzzyClassification(sCImage *cimg,
				       int *S_obj,
				       int *S_bkg,
				       float ds,
				       int k,
				       int maxsamples,
				       int Pmax){
	glip::sFeatures *f_obj = NULL, *f_bkg = NULL;
	glip::sImage32 *membership;
	int p, n;
	float fv[3];
	float d_obj, d_bkg, dmax;
	membership = glip::Image32::Create(cimg->C[0]);
	dmax = sqrtf(255*255*3);
	f_obj = glip::Features::GetSamples(cimg, S_obj);
	f_bkg = glip::Features::GetSamples(cimg, S_bkg);
	glip::Features::Randomize(f_obj);
	glip::Features::Randomize(f_bkg);
	if(f_obj->n > maxsamples) glip::Features::Resize(&f_obj, maxsamples);
	if(f_bkg->n > maxsamples) glip::Features::Resize(&f_bkg, maxsamples);
	n = cimg->C[0]->n;
	for(p = 0; p < n; p++){
	  fv[0] = cimg->C[0]->data[p];
	  fv[1] = cimg->C[1]->data[p];
	  fv[2] = cimg->C[2]->data[p];
	  d_obj = glip::Features::KNNMeanDistance(f_obj, fv, k);
	  d_bkg = glip::Features::KNNMeanDistance(f_bkg, fv, k);
	  membership->data[p] = ROUND(Pmax*((d_bkg+ds)/(d_obj+d_bkg+2.0*ds)));
	}
	glip::Features::Destroy(&f_obj);
	glip::Features::Destroy(&f_bkg);
	return membership;
      }


      sImage32 *KNNFuzzyClassification(sImage32 *img,
				       int *S_obj,
				       int *S_bkg,
				       float ds,
				       int k,
				       int maxsamples,
				       int Pmax){
	glip::sFeatures *f_obj = NULL, *f_bkg = NULL;
	glip::sImage32 *membership;
	int p, n;
	float fv[3];
	float d_obj, d_bkg, dmax;
	membership = glip::Image32::Create(img);
	dmax = glip::Image32::GetMaxVal(img) - glip::Image32::GetMinVal(img);
	f_obj = glip::Features::GetSamples(img, S_obj);
	f_bkg = glip::Features::GetSamples(img, S_bkg);
	glip::Features::Randomize(f_obj);
	glip::Features::Randomize(f_bkg);
	if(f_obj->n > maxsamples) glip::Features::Resize(&f_obj, maxsamples);
	if(f_bkg->n > maxsamples) glip::Features::Resize(&f_bkg, maxsamples);
	n = img->n;
	for(p = 0; p < n; p++){
	  fv[0] = img->data[p];
	  d_obj = glip::Features::KNNMeanDistance(f_obj, fv, k);
	  d_bkg = glip::Features::KNNMeanDistance(f_bkg, fv, k);
	  membership->data[p] = ROUND(Pmax*((d_bkg+ds)/(d_obj+d_bkg+2.0*ds)));
	}
	glip::Features::Destroy(&f_obj);
	glip::Features::Destroy(&f_bkg);
	return membership;
      }


      sImage32 *KNNFuzzyClassification(sCImage32f *cimg,
				       int *S_obj,
				       int *S_bkg,
				       float ds,
				       int k,
				       int maxsamples,
				       int Pmax){
	glip::sFeatures *f_obj = NULL, *f_bkg = NULL;
	glip::sImage32 *membership;
	int p, n;
	float fv[3];
	float d_obj, d_bkg, dmax;
	membership = glip::Image32::Create(cimg->C[0]->ncols,
					  cimg->C[0]->nrows);
	dmax = sqrtf(255*255*3);
	f_obj = glip::Features::GetSamples(cimg, S_obj);
	f_bkg = glip::Features::GetSamples(cimg, S_bkg);
	glip::Features::Randomize(f_obj);
	glip::Features::Randomize(f_bkg);
	if(f_obj->n > maxsamples) glip::Features::Resize(&f_obj, maxsamples);
	if(f_bkg->n > maxsamples) glip::Features::Resize(&f_bkg, maxsamples);
	n = cimg->C[0]->n;
	for(p = 0; p < n; p++){
	  fv[0] = cimg->C[0]->data[p];
	  fv[1] = cimg->C[1]->data[p];
	  fv[2] = cimg->C[2]->data[p];
	  d_obj = glip::Features::KNNMeanDistance(f_obj, fv, k);
	  d_bkg = glip::Features::KNNMeanDistance(f_bkg, fv, k);
	  membership->data[p] = ROUND(Pmax*((d_bkg+ds)/(d_obj+d_bkg+2.0*ds)));
	}
	glip::Features::Destroy(&f_obj);
	glip::Features::Destroy(&f_bkg);
	return membership;
      }



      sScene32 *KNNFuzzyClassification(sScene32 *scn,
				       int *S_obj,
				       int *S_bkg,
				       float ds,
				       int k,
				       int maxsamples,
				       int Pmax){
	glip::sFeatures *f_obj = NULL, *f_bkg = NULL;
	glip::sScene32 *membership;
	int p, n;
	float fv[3];
	float d_obj, d_bkg, dmax;
	membership = glip::Scene32::Create(scn);
	dmax = glip::Scene32::GetMaxVal(scn) - glip::Scene32::GetMinVal(scn);
	f_obj = glip::Features::GetSamples(scn, S_obj);
	f_bkg = glip::Features::GetSamples(scn, S_bkg);
	glip::Features::Randomize(f_obj);
	glip::Features::Randomize(f_bkg);
	if(f_obj->n > maxsamples) glip::Features::Resize(&f_obj, maxsamples);
	if(f_bkg->n > maxsamples) glip::Features::Resize(&f_bkg, maxsamples);
	n = scn->n;
	for(p = 0; p < n; p++){
	  fv[0] = scn->data[p];
	  d_obj = glip::Features::KNNMeanDistance(f_obj, fv, k);
	  d_bkg = glip::Features::KNNMeanDistance(f_bkg, fv, k);
	  membership->data[p] = ROUND(Pmax*((d_bkg+ds)/(d_obj+d_bkg+2.0*ds)));
	}
	glip::Features::Destroy(&f_obj);
	glip::Features::Destroy(&f_bkg);
	return membership;
      }
      
      

      sImage32f *KNNMeanDistanceMap(sCImage32f *cimg,
				    int *S,
				    int k,
				    int maxsamples){
	glip::sFeatures *f = NULL;
	glip::sImage32f *D;
	int p, n;
	float fv[3];
	float d;
	D = glip::Image32f::Create(cimg->C[0]->ncols,
				  cimg->C[0]->nrows);
	f = glip::Features::GetSamples(cimg, S);
	glip::Features::Randomize(f);
	if(f->n > maxsamples) glip::Features::Resize(&f, maxsamples);
	n = cimg->C[0]->n;
	for(p = 0; p < n; p++){	
	  fv[0] = cimg->C[0]->data[p];
	  fv[1] = cimg->C[1]->data[p];
	  fv[2] = cimg->C[2]->data[p];
	  d = glip::Features::KNNMeanDistance(f, fv, k);
	  D->data[p] = d;
	}
	glip::Features::Destroy(&f);
	return D;
      }
      
      
    } //end Features namespace
} //end glip namespace


