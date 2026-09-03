#ifndef _GLIP_ADJREL_H_
#define _GLIP_ADJREL_H_

#include "glip_common.h"

namespace glip{
  namespace AdjRel{

    struct sAdjRel {
      int *dx;
      int *dy;
      int n;
    };

    struct sAdjPxl {
      int *dp;
      int n;
    };
    
    
    sAdjRel *Create(int n);
    void    Destroy(sAdjRel **A);
    sAdjRel *Clone(sAdjRel *A);

    /**
     * \brief Returns the size of the data structure in bytes.
     */
    size_t GetSizeOf(sAdjRel *A);
    
    sAdjRel *Neighborhood_4(); /* 4-neighborhood */
    sAdjRel *Neighborhood_8(); /* 8-neighborhood */
    sAdjRel *Neighborhood_8_counterclockwise();
    sAdjRel *Neighborhood_8_clockwise();
   
    sAdjRel *Circular(float r);
    sAdjRel *Box(int ncols, int nrows);

    //-----------------------------------
    int GetFrameSize(sAdjRel *A);

    int *InverseIndexes(sAdjRel *A);
    
    void Mult(sAdjRel *A, int val);

    sAdjRel *RightSide(sAdjRel *A);
    sAdjRel *LeftSide(sAdjRel *A);

    sAdjRel *RightSide8(sAdjRel *A);
    sAdjRel *LeftSide8(sAdjRel *A);

    void Reduce2Boundary(sAdjRel **A, float r);
    
  } /*end AdjRel namespace*/

  typedef AdjRel::sAdjPxl sAdjPxl;
  typedef AdjRel::sAdjRel sAdjRel;
  
} /*end glip namespace*/


#include "glip_image32.h"

namespace glip{

  namespace Image32{

    sImage32 *Render(sAdjRel *A);

    void DrawAdjRel(sImage32 *img,
		    sAdjRel *A, 
		    int p, int val);

    /*
    sImage32 *GetBoundaries(sImage32 *img, sAdjRel *A);
    sImage32 *GetBoundaries(sImage32 *img, sAdjRel *A, bool frame);
    */
    
  } /*end Image32 namespace*/
 
} /*end glip namespace*/




#include "glip_cimage.h"

namespace glip{
  namespace CImage{

    void DrawAdjRel(sCImage *cimg,
		    sAdjRel *A, 
		    int p, int color);

    
  } //end CImage namespace
} //end glip namespace



namespace glip{
  namespace AdjRel{

    sAdjPxl  *AdjPixels(sAdjRel *A, int ncols);
    sAdjPxl  *AdjPixels(sAdjRel *A, sImage32 *img);

    sAdjPxl  *AdjPixels(sAdjRel *A, sCImage *cimg);

    sAdjPxl  *CreateAdjPxl(int n);
    
    void DestroyAdjPxl(sAdjPxl **N);
    
  } /*end AdjRel namespace*/
} /*end glip namespace*/


#endif

