#ifndef _GLIP_CIMAGE32F_H_
#define _GLIP_CIMAGE32F_H_

#include "glip_image32f.h"
#include "glip_cimage.h"
#include "glip_color.h"

namespace glip{
  namespace CImage32f{

    struct sCImage32f {
      sImage32f *C[3];
    };

    sCImage32f *Create(int ncols, int nrows);
    sCImage32f *Create(sCImage32f *cimg);
    sCImage32f *Create(sImage32f *img);
    void       Destroy(sCImage32f **cimg);
    sCImage32f *Clone(sCImage32f *cimg);
    sCImage32f *Clone(sCImage *cimg);

    void    Set(sCImage32f *cimg, float r, float g, float b);
    
    sCImage32f *RGB2Lab(sCImage *cimg);
 
    sCImage32f *AddFrame(sCImage32f *cimg, int sz, float r, float g, float b);
    sCImage32f *RemFrame(sCImage32f *cimg, int sz);
    
  } //end CImage32f namespace

  typedef CImage32f::sCImage32f sCImage32f;

} //end glip namespace


namespace glip{
  namespace Image32{

    /**
     * \brief Changes an image from RGB to grayscale.
     */
    sImage32 *Lightness(sCImage32f *cimg);

    /**
     * \brief Changes an image from RGB to grayscale.
     *
     * The luminosity method works best overall and is the 
     * default method used if you ask GIMP to change an image 
     * from RGB to grayscale 
     */
    sImage32 *Luminosity(sCImage32f *cimg);

    
  } //end Image32 namespace
} //end glip namespace


#endif

