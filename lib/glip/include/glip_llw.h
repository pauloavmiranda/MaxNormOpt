
#ifndef _GLIP_LWARCW_H_
#define _GLIP_LWARCW_H_


#include "glip_common.h"
#include "glip_image32.h"
#include "glip_filtering.h"
#include "glip_set.h"
#include "glip_heap64f.h"
#include "glip_analysis.h"


namespace glip{

  namespace LLW{
  
    glip::sImage32 *WeightImage(glip::sImage32 *img,
			       float r);
    
    double *Curvature(glip::sImage32 *img,
		      glip::sImage32 *W,
		      glip::sImage32 *mask,
		      float r,
		      bool draw,
		      double **V_curv_x,
		      double **V_curv_y);
    //----------------------------------------

    double *Curvature(glip::sImage32 *img,
		      glip::sImage32 *spixels,
		      float r,
		      double **V_curv_x,
		      double **V_curv_y);
    
  } //end LLW namespace
 
} //end glip namespace

#endif

