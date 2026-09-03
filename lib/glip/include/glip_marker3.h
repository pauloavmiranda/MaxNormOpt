
#ifndef _GLIP_MARKER3_H_
#define _GLIP_MARKER3_H_

#include "glip_common.h"
#include "glip_scene32.h"
#include "glip_set.h"
#include "glip_adjrel3.h"
#include "glip_analysis3.h"
#include "glip_morphology3.h"
#include "glip_marker.h"

namespace glip{

  namespace Scene32{
    
    int *GetMarkers(sScene32 *label, sAdjRel3 *A);

    void SaveMarkers(char *filename,
		     int *S,
		     sScene32 *label);
    
    int *GetVolumeBorder(sScene32 *scn);
    int *GetVolumeBorder(sScene32 *scn,
			 uchar exception_flags);
    
  } //end Scene32 namespace

  

  namespace Scene8{

    int *Convert2Array(sScene8 *mask);
    int *GetMarkers(sScene8 *label, sAdjRel3 *A);

    void SaveMarkers(char *filename,
		     int *S,
		     sScene8 *label);
    
    int *GetVolumeBorder(sScene8 *scn);
    int *GetVolumeBorder(sScene8 *scn,
			 uchar exception_flags);
    
  } //end Scene8 namespace
  
  
} //end glip namespace

 
#endif
