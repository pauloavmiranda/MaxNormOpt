#ifndef _GLIP_IMAGE32_H_
#define _GLIP_IMAGE32_H_

#include "glip_common.h"

namespace glip{

  typedef enum {INTER_NONE, INTER_NN, INTER_LINEAR, INTER_CUBIC} InterpolationType;

  /**
   * \brief Common definitions and functions to manipulate a grayscale image of 32 bits per pixel.
   */
  namespace Image32{

    /**
     * \brief Structure representing a digital image. 
     *
     * It supports both linear and two-dimensional access 
     * (i.e., img->data[p] or img->array[y][x] for a pixel
     * (x,y) at address p=x+y*ncols).
     */
    struct sImage32 {
      int *data;
      int **array;
      int nrows; /* numero de linhas (altura) */
      int ncols; /* numero de colunas (largura) */
      int n;     /* numero de pixels */
      int maxval, minval;
      float dx;
      float dy;
    };

    /**
     * \brief A constructor.
     */
    sImage32 *Create(int ncols,int nrows);
    /**
     * \brief A constructor using img as a template.
     */
    sImage32 *Create(sImage32 *img);
    
    /**
     * \brief A destructor.
     */
    void    Destroy(sImage32 **img);
    
    /**
     * \brief A copy constructor.
     */
    sImage32 *Clone(sImage32 *img);
    /**
     * \brief A copy constructor.
     */
    sImage32 *Clone(sImage32 *img, Pixel l, Pixel h);

    /**
     * \brief Copies the contents of one image to another destination image.
     */    
    void     Copy(sImage32 *img,  sImage32 *sub, Pixel l);
    /**
     * \brief Copies the contents of one image to another destination image.
     */
    void     Copy(sImage32 *img,  sImage32 *sub, Pixel l, int bkg);    
    /**
     * \brief Copies the contents of one image to another destination image.
     */
    void     Copy(sImage32 *dest, sImage32 *src);

    /**
     * \brief Returns the size of the data structure in bytes.
     */
    size_t GetSizeOf(sImage32 *img);
    
    /**
     * \brief Calculates the sum of two images pixel by pixel. 
     */
    sImage32 *Add( sImage32 *img1, sImage32 *img2);
    /**
     * \brief Adds a constant to all pixels in the image.
     */
    sImage32 *Add( sImage32 *img,  int value);

    /**
     * \brief Calculates the product of two images pixel by pixel.
     *
     * Gives an image in which each pixel is the product of the 
     * corresponding pixels in img1 and img2.
     */
    sImage32 *Mult(sImage32 *img1, sImage32 *img2);

    /**
     * \brief Reads an image from a file on the disc in PGM format.
     */
    sImage32 *Read(char *filename);
    /**
     * \brief Writes an image to a file on the disc in PGM format. 
     */
    void      Write(sImage32 *img, char *filename);

    /**
     * \brief Changes the amount of gray levels in the image to a given bit depth.
     */
    sImage32 *ConvertToNbits(sImage32 *img, int N);
    /**
     * \brief Changes the amount of gray levels in the image to a given bit depth.
     */
    sImage32 *ConvertToNbits(sImage32 *img, int N, bool Imin);

    /**
     * \brief Calculates the image complement by inverting its values (white becomes black, black becomes white).
     *
     * In the complement of a grayscale image, each pixel value is 
     * subtracted from the maximum pixel value. The difference is 
     * used as the pixel value in the output image. In the output image, 
     * dark areas become lighter and light areas become darker.
     */
    sImage32 *Complement(sImage32 *img);

    /**
     * \brief Returns the minimum pixel value present in a given image.
     */
    int     GetMinVal(sImage32 *img);
    /**
     * \brief Returns the minimum pixel value present in a given image.
     *
     * @param img The input image.
     * @param p Pointer to a variable used to store a pixel address having the minimum value.
     */
    int     GetMinVal(sImage32 *img, int *p);
    /**
     * \brief Returns the maximum pixel value present in a given image.
     */
    int     GetMaxVal(sImage32 *img);
    /**
     * \brief Returns the maximum pixel value present in a given image.
     *
     * @param img The input image.
     * @param p Pointer to a variable used to store a pixel address having the maximum value.
     */
    int     GetMaxVal(sImage32 *img, int *p);
    /**
     * \brief Returns the maximum pixel value present in a given image.
     */
    int     GetMaxVal(sImage32 *img, int ignoredvalue);
    /**
     * \brief Returns the minimum pixel value present in a given image.
     */
    int     GetMinVal(sImage32 *img, int ignoredvalue);
    
    /**
     * \brief Returns the number of occurrences (frequency) of a given value in an image. 
     */
    int     GetFreqVal(sImage32 *img, int val);

    /**
     * \brief Fills an image with a given value.
     */
    void    Set(sImage32 *img, int value);
    /**
     * \brief Replaces a given image value with a new value.
     */
    void    Set(sImage32 *img, int old_value, int new_value);    

    /**
     * \brief Tests whether a pixel is valid for a given image.
     */
    bool    IsValidPixel(sImage32 *img, int x, int y);
    /**
     * \brief Tests whether a pixel is valid for a given image.
     */
    bool    IsValidPixel(sImage32 *img, Pixel u);
    
    /**
     * \brief Calculates the segmentation/binarization by thresholding the image.
     */
    sImage32 *Threshold(sImage32 *img, int L, int H);
    
    //------------------------------------
    /**
     * \brief Adds a frame around an image, with the provided width and fill value.
     */
    sImage32 *AddFrame(sImage32 *img, int sz, int value);
    /**
     * \brief Removes the frame around an image with the given width.
     */
    sImage32 *RemFrame(sImage32 *fimg, int sz);

    //------------------------------------
    /**
     * \brief Resizes the image by the horizontal and vertical factors and interpolation method provided.
     *
     * @param Sx The horizontal resize factor.
     * @param Sy The vertical resize factor.
     */
    sImage32 *Scale(sImage32 *img, float Sx, float Sy,
		    InterpolationType interpolation);
   
    /**
     * \brief Returns a subimage by clipping the minimum bounding box with non-zero values.
     */
    sImage32 *MBB(sImage32 *img);
    /**
     * \brief Calculates the minimum bounding box with non-zero values.
     */
    void      MBB(sImage32 *img, Pixel *l, Pixel *h);


    //---------inline definitions------------------
    inline int GetVal(sImage32 *img, float x, float y,
		      InterpolationType interpolation){
      int i,j;
      if(interpolation == INTER_NONE){
	i = (int)y;
	j = (int)x;
	if(j >= 0 && j < img->ncols &&
	   i >= 0 && i < img->nrows)
	  return img->array[i][j];
	else
	  return 0;
      }
      else if(interpolation == INTER_NN){
	i = ROUND(y);
	j = ROUND(x);
	if(j >= 0 && j < img->ncols &&
	   i >= 0 && i < img->nrows)
	  return img->array[i][j];
	else
	  return 0;
      }
      else if(interpolation == INTER_LINEAR){
	int xl,xr,yb,yt, v1,v2,v3,v4;
	float dx,dy, v12, v34;
	xl = (int)x;
	xr = (int)x + 1;
	yt = (int)y;
	yb = (int)y + 1;
	if(xl >= 0 && xr < img->ncols &&
	   yt >= 0 && yb < img->nrows){
	  v1 = img->array[yt][xl];  v2 = img->array[yt][xr];
	  v3 = img->array[yb][xl];  v4 = img->array[yb][xr];
	  dx = x - xl; dy = y - yt;
	  v12 = (1.0-dx)*v1 + dx*v2;
	  v34 = (1.0-dx)*v3 + dx*v4;
	  return ROUND((1.0-dy)*v12 + dy*v34);
	}
	else
	  return 0;
      }
      else{
	Warning((char *)"Interpolation type not yet implemented",
		(char *)"GetVal");
	return 0;
      }
    }

    
  } //end Image32 namespace

  typedef Image32::sImage32 sImage32;

} //end glip namespace



#include "glip_matrix.h"

namespace glip{
  namespace Image32{

    sImage32 *Create(sMatrix *M);

    sImage32 *Create(sMatrix *M, int Imin, int Imax);
    
    sImage32 *Threshold(sMatrix *M, float L, float H);
    
    sImage32 *WarpPerspective(sImage32 *img,
			      sMatrix *M,
			      int ncols, int nrows,
			      InterpolationType interpolation);
    
    
  } //end Image32 namespace
} //end glip namespace



#endif

