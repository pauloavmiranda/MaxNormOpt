
#ifndef _GLIP_FEATURES_H_
#define _GLIP_FEATURES_H_

#include "glip_common.h"
#include "glip_image32.h"
#include "glip_cimage.h"
#include "glip_cimage32f.h"
#include "glip_scene32.h"

/* It supports the following styles: */

#define glipFeature_ALLOCLABEL  0x000001
#define glipFeature_ALLOCINDEX  0x000002
#define glipFeature_ALLOCFV     0x000004
#define glipFeature_ALLOCDIST   0x000008
#define glipFeature_ALLOCALL    0x00000F


namespace glip{

    namespace Features{


      struct sFeatures {
	int n;
	int nfeats;
	int *label; /* For supervised learning from labeled training data */
	int *index; /* Feature vector position in the image */
	float **fv; /* fv[0..n-1][0..nfeats-1]*/

	float **dist; /* dist[0..n-1][0..n-1] */
      };


      sFeatures *Create(int n, int nfeats, int style);
      sFeatures *Clone(sFeatures *f);
      void       Destroy(sFeatures **f);

      /*It shuffles the samples at random.*/
      void Randomize(sFeatures *f);

      /*If n is less than the current size, then 
        it discards the excess samples.*/
      void Resize(sFeatures **f, int n);
      
      /*Returns a subset with samples of all classes in the 
	proportion given by "rate". This subset is removed from the
	original dataset. You should call "Randomize"
	first to ensure the selection of random samples.*/
      sFeatures *RemoveSamples(sFeatures **f,
			       float rate);
      
      /*Concatenate two datasets.*/
      sFeatures *Merge(sFeatures *f1, 
		       sFeatures *f2);

      int *KNN(sFeatures *f,
	       float *fv,
	       int k);

      float KNNMeanDistance(sFeatures *f,
			    float *fv,
			    int k);
      
      sFeatures *GetSamples(sImage32 *img,
			    int *S);
     
      sFeatures *GetSamples(sCImage *cimg,
			    int *S);

      sFeatures *GetSamples(sCImage32f *cimg,
			    int *S);

      sFeatures *GetSamples(sScene32 *scn,
			    int *S);
      
      sImage32 *KNNFuzzyClassification(sCImage *cimg,
				       int *S_obj,
				       int *S_bkg,
				       float ds,
				       int k,
				       int maxsamples,
				       int Pmax);

      sImage32 *KNNFuzzyClassification(sImage32 *img,
				       int *S_obj,
				       int *S_bkg,
				       float ds,
				       int k,
				       int maxsamples,
				       int Pmax);
     
      sImage32 *KNNFuzzyClassification(sCImage32f *cimg,
				       int *S_obj,
				       int *S_bkg,
				       float ds,
				       int k,
				       int maxsamples,
				       int Pmax);

      sScene32 *KNNFuzzyClassification(sScene32 *scn,
				       int *S_obj,
				       int *S_bkg,
				       float ds,
				       int k,
				       int maxsamples,
				       int Pmax);
      
      sImage32f *KNNMeanDistanceMap(sCImage32f *cimg,
				    int *S,
				    int k,
				    int maxsamples);
      
    } //end Features namespace

    typedef Features::sFeatures sFeatures;
    
} //end glip namespace



#endif

