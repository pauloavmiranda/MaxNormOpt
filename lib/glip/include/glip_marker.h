
#ifndef _GLIP_MARKER_H_
#define _GLIP_MARKER_H_

#include "glip_common.h"
#include "glip_image32.h"
#include "glip_set.h"
#include "glip_adjrel.h"
#include "glip_analysis.h"
#include "glip_morphology.h"


namespace glip{
#define EXCEPT_BORDER_MAX_X    0b00000001
#define EXCEPT_BORDER_MIN_X    0b00000010
#define EXCEPT_BORDER_MAX_Y    0b00000100
#define EXCEPT_BORDER_MIN_Y    0b00001000
#define EXCEPT_BORDER_MAX_Z    0b00010000
#define EXCEPT_BORDER_MIN_Z    0b00100000
  
  namespace Image32{
    
    float  MaxRadiusByErosion(sImage32 *bin);
    float  MaxObjRadiusByErosion(sImage32 *bin);
    float  MaxBkgRadiusByErosion(sImage32 *bin);
    
    sImage32 *ObjMaskByErosion(sImage32 *bin, float radius);
    sImage32 *BkgMaskByErosion(sImage32 *bin, float radius);

    int  *GetMarkers(sImage32 *label, sAdjRel *A);

    int *GetImageBorder(sImage32 *img);
    int *GetImageBorder(sImage32 *img,
			uchar exception_flags);
    
   } //end Image32 namespace
} //end glip namespace

 
#endif
