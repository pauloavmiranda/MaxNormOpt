
#ifndef _GLIP_MORPHOLOGY_H_
#define _GLIP_MORPHOLOGY_H_

#include "glip_image32.h"
#include "glip_adjrel.h"
#include "glip_set.h"
#include "glip_pqueue32.h"


namespace glip{
  namespace Image32{

    sImage32 *Dilate(sImage32 *img, sAdjRel *A);
    sImage32 *Erode(sImage32 *img, sAdjRel *A);

    sImage32 *ErodeBin(sImage32 *bin, sSet **seed, float radius);
    
    sImage32 *MorphGrad(sImage32 *img, sAdjRel *A);

    void SupRec_Watershed(sAdjRel *A, 
			  sImage32 *I, sImage32 *J, 
			  sImage32 *L, sImage32 *V);

    /*Removes all background connected components from the stack
      of binary images of I whose area (number of pixels) is <=
      a threshold and outputs a simplified image.*/
    sImage32 *AreaClosing(sAdjRel *A, sImage32 *I, int T);

    sImage32 *VolumeClosing(sAdjRel *A, sImage32 *I, int T);
    
    sImage32 *CloseHoles(sImage32 *img);

    sImage32 *CloseHoles(sImage32 *img, float r);    
    
    sImage32 *CloseHolesFromSeeds(sImage32 *img,
				  sAdjRel *A,
				  int *S);
    
  } //end Image32 namespace
} //end glip namespace

#endif
