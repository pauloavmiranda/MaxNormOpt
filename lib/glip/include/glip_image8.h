#ifndef _GLIP_IMAGE8_H_
#define _GLIP_IMAGE8_H_

#include "glip_common.h"
#include "glip_image32.h"

namespace glip{
  namespace Image8{

    struct sImage8 {
      uchar *data;
      uchar **array;
      int nrows; /* numero de linhas (altura) */
      int ncols; /* numero de colunas (largura) */
      int n;     /* numero de pixels */
    };


    sImage8 *Create(int ncols, int nrows);
    void     Destroy(sImage8 **img);

    /**
     * \brief Writes an image to a file on the disc in PGM format. 
     */
    void     Write(sImage8 *img, char *filename);

    /**
     * \brief Returns the maximum pixel value present in a given image.
     */
    int     GetMaxVal(sImage8 *img);
    
  } /*end Image8 namespace*/

  typedef Image8::sImage8 sImage8;

} /*end glip namespace*/

#endif

