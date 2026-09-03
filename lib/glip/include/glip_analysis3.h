#ifndef _GLIP_ANALYSIS3_H_
#define _GLIP_ANALYSIS3_H_

#include "glip_common.h"
#include "glip_scene32.h"
#include "glip_scene8.h"
#include "glip_adjrel3.h"
#include "glip_pqueue32.h"
#include "glip_stack.h"
#include "glip_queue.h"

namespace glip{
  namespace Scene8{

    int ComputeBinCoG(sScene8 *bin);
    
    sScene8 *GetObjBorders(sScene8 *scn, sAdjRel3 *A);
    sScene8 *GetObjBorders(sScene8 *scn, sAdjRel3 *A, bool frame);
    
    sScene32 *Mask2EDT(sScene8 *bin, sAdjRel3 *A,
		       char side, int limit, char sign);
    
    sScene32 *LabelBinComp(sScene8 *bin, sAdjRel3 *A);
    void SelectLargestComp(sScene8 *bin, sAdjRel3 *A);
    void SelectLargestComp(sScene8 *bin);
    
    int GetBinArea(sScene8 *bin);
    
  } //end Scene8 namespace


} //end glip namespace

#endif

