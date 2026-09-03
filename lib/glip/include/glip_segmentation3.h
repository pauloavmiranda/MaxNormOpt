
#ifndef _GLIP_SEGMENTATION3_H_
#define _GLIP_SEGMENTATION3_H_

#include "glip_common.h"
#include "glip_scene.h"
#include "glip_adjrel3.h"
#include "glip_radiometric3.h"

namespace glip{

  namespace Scene8{

    sScene8 *Threshold(sScene8 *scn, int lower, int higher);
    //sScene8 *GetBoundaries(sScene8 *scn, sAdjRel3 *A);
    sScene8 *GetTransitions(sScene8 *scn, sAdjRel3 *A);

  } //end Scene8 namespace



  namespace Scene16{

    sScene8 *Threshold(sScene16 *scn, int lower, int higher);

  } //end Scene16 namespace



  namespace Scene32{

    sScene8 *Threshold(sScene32 *scn, int lower, int higher);
    int Otsu(sScene32 *scn);
   
  } //end Scene32 namespace


  namespace Scene{

    sScene8 *Threshold(sScene *scn, int lower, int higher);

  } //end Scene namespace


} //end glip namespace


#endif

