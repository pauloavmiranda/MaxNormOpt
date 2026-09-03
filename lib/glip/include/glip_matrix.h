
#ifndef _GLIP_MATRIX_H_
#define _GLIP_MATRIX_H_

#include "glip_common.h"

namespace glip{
  namespace Matrix{

    /**
     * It supports both linear and two-dimensional access 
     * (i.e., M->data[p] or M->array[i][j] for an entry
     * (i,j) at address p=j+i*ncols).
     */
    struct sMatrix {
      float *data;
      float **array;
      int ncols,nrows;
    };


    sMatrix *Create(int ncols,int nrows);
    sMatrix *Create(int ncols,int nrows,
		    char str[]);

    void     Destroy(sMatrix **mat);
    sMatrix *Clone(sMatrix *mat);
    void     Copy(sMatrix *dest, 
		  sMatrix *src);

    /**
     * \brief Creates an n by n identity matrix.
     */
    sMatrix *Identity(int n);

    //check_equality
    /**
     * \brief The function checks if the matrices have the same dimensions and if all their corresponding elements are identical (A[i,j] == B[i,j]).
     */
    bool Equal(sMatrix *A, sMatrix *B);

    /**
     * \brief The function checks if the matrices have the same dimensions and if all their corresponding elements are identical, up to a tolerance error, in order to deal with variations due to floating-point precision errors.
It checks if the difference between corresponding elements falls within a tiny tolerance range, rather than requiring exact equality (abs(A[i,j] - B[i,j]) < tolerance).
     */
    bool Equal(sMatrix *A, sMatrix *B, float tolerance);

    /**
     * \brief The function solves a system of linear equations Ax = b by the Gaussian elimination method with partial pivoting.
     */
    float *SolveLinearSystem(sMatrix *A, float *b);

    /**
     * \brief Calculates coefficients of perspective transformation which maps (xi,yi) to (ui,vi), (i=0,1,2,3).
     */
    sMatrix *getPerspectiveTransform(float x[4], float y[4],
				     float u[4], float v[4]);

    /**
     * \brief Matrix inversion by Gaussian elimination with partial pivoting.
     */
    sMatrix *Invert(sMatrix *A);

    sMatrix *Transpose(sMatrix *A);
    
    sMatrix *Mult(sMatrix *A, 
		  sMatrix *B);
    sMatrix *MultByScalar(sMatrix *A, float k);

    /**
     * \brief Computes the matrix-vector product M*x.
     *
     * @param out Vector used to store the computed output.
     */
    void Mult(sMatrix *M, float x[], float out[]);
    
    sMatrix *Sub(sMatrix *A, 
		 sMatrix *B);
    sMatrix *Add(sMatrix *A, 
		 sMatrix *B);

    float   GetTrace(sMatrix *M);
    
    void    Print(sMatrix *M);
    void    Print(sMatrix *M, int min_width, int dec_places);
    void    PrintDimension(sMatrix *M);
    
    float   ComputeDistanceL2(sMatrix *Y, 
			      sMatrix *X);
    
    void    Fill(sMatrix *M, float value);
    void    ChangeValue(sMatrix *M, 
			float old_value,
			float new_value);
    
    bool    IsValidEntry(sMatrix *M,
			 int i, int j);
    
    
    sMatrix *Read(char *filename);
    void     Write(sMatrix *M,
		   char *filename);
    
    float   GetMinimumValue(sMatrix *M);
    float   GetMaximumValue(sMatrix *M);

    /**
     * @param axis an option (0->x / 1->y / 2->z).
     */
    sMatrix* RotationMatrix3(int axis, 
			     float th); 
    
    sMatrix* TranslationMatrix3(float dx, float dy, float dz);
    
    sMatrix* TransformVoxel(sMatrix *m, Voxel v);

  } //end Matrix namespace

  typedef Matrix::sMatrix sMatrix;

} //end glip namespace


#include "glip_image32.h"


namespace glip{
  namespace Matrix{

    /**
     * \brief Extracts the two-dimensional binary code from a binary image into a matrix (e.g., binary code of an ArUco marker or QR code).
     */    
    sMatrix *Decode(sImage32 *bin, int ncols, int nrows);

    
  } //end Matrix namespace
} //end glip namespace


#endif


