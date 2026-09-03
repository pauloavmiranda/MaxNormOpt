#ifndef _GLIP_ANALYSIS_H_
#define _GLIP_ANALYSIS_H_

#include "glip_common.h"
#include "glip_image32.h"
#include "glip_adjrel.h"
#include "glip_pqueue32.h"
#include "glip_stack.h"
#include "glip_queue.h"

namespace glip{
  namespace Image32{

    /**
     * \brief Calculates the Center of Gravity (CoG) for each label. 
     *
     * Returns an array of integers C, where C[x] is the CoG of the label x.
     * Each value p = C[x] is a pixel index, where p = x + y*ncols.
     */
    int *ComputeCoG(sImage32 *label);

    /**
     * \brief Calculates the area for each label. 
     *
     * Returns an array of integers A, where A[x] is the area of the label x.
     */
    int *ComputeArea(sImage32 *label);
    
    int ComputeBinCoG(sImage32 *bin);
    
    sImage32 *RegMin(sImage32 *img, sAdjRel *A);
    sImage32 *LabelBinComp(sImage32 *bin, sAdjRel *A);

    void SelectLargestComp(sImage32 *bin);
    
    //sImage32 *DistTrans(sImage32 *bin, sAdjRel *A, char side); 
    //sImage32 *SignedDistTrans(sImage32 *bin, sAdjRel *A, char side);

    sImage32 *GetObjBorder(sImage32 *bin);
    sImage32 *GetObjBorders(sImage32 *img, sAdjRel *A);
    sImage32 *GetObjBorders(sImage32 *img, sAdjRel *A, bool frame);
    sImage32 *GetBorders(sImage32 *img, sAdjRel *A);
    
    sImage32 *LabelContour(sImage32 *bin);
    sImage32 *LabelContour(sImage32 *bin, sImage32 *contourid);
    
    sImage32 *Mask2EDT(sImage32 *bin, sAdjRel *A,
		       char side, int limit, char sign);
    
    void Mask2EDT(sImage32 *bin, sAdjRel *A,
		  char side, int limit, char sign,
		  sImage32 *cost, sImage32 *root);

    sImage32 *Multiscaleskeletons(sImage32 *bin);


    float PerimeterLength(sImage32 *bin);

    int GetBinArea(sImage32 *bin);
   
    
  } //end Image32 namespace
} //end glip namespace

#endif

