
#ifndef _GLIP_RADIOMETRIC_H_
#define _GLIP_RADIOMETRIC_H_

#include "glip_common.h"
#include "glip_image32.h"
#include "glip_curve.h"

namespace glip{

  namespace Image32{

    sCurve *Histogram(sImage32 *img);
    sCurve *NormHistogram(sImage32 *img);
    sCurve *NormalizeHistogram(sCurve *hist);
    sCurve *RemoveEmptyBins(sCurve *hist);
    
    sImage32 *LinearStretch(sImage32 *img,
			    int omin, int omax,
			    int nmin, int nmax);
    void LinearStretchinplace(sImage32 *img, 
			      int omin, int omax, 
			      int nmin, int nmax);
    
  } //end Image32 namespace
    

} //end glip namespace

#endif


