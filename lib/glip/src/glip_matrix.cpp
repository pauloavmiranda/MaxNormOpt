
#include "glip_matrix.h"

namespace glip{
  namespace Matrix{

    sMatrix *Create(int ncols,int nrows){
      sMatrix *mat=NULL;
      float *aux;
      int i;
      
      mat = (sMatrix *) calloc(1,sizeof(sMatrix));
      if(mat == NULL)
	glip::Error((char *)MSG1,(char *)"Matrix::Create");
      
      aux = (float *)calloc(nrows*ncols, sizeof(float));
      mat->array = (float **) calloc(nrows, sizeof(float *));
      if(mat->array == NULL || aux == NULL)
	glip::Error((char *)MSG1,(char *)"Matrix::Create");

      mat->data = aux;
      mat->array[0] = aux;
      for(i=1; i<nrows; i++) 
	mat->array[i] = mat->array[i-1] + ncols;
    
      mat->ncols = ncols;
      mat->nrows = nrows;
      return(mat);
    }


    sMatrix *Create(int ncols, int nrows,
		    char str[]){
      sMatrix *mat = NULL;
      float val;
      int i,n,p;
      mat = Create(ncols, nrows);
      p = 0;
      i = 0;
      while(sscanf(&str[i], " %f%n", &val,&n) != EOF){
	mat->data[p] = val;
	i += n;
	p++;
      }
      return mat;
    }

    
    void Destroy(sMatrix **mat){
      sMatrix *aux;
      if(mat != NULL){
	aux = *mat;
	if (aux != NULL){
	  if(aux->data !=  NULL) free(aux->data);
	  if(aux->array != NULL) free(aux->array);
	  free(aux);
	  *mat = NULL;
	}
      }
    }
    

    sMatrix *Clone(sMatrix *mat){
      sMatrix *matc;
      matc = Create(mat->ncols, mat->nrows);
      memcpy(matc->array[0], mat->array[0],
	     mat->ncols*mat->nrows*sizeof(float));
      return(matc);
    }


    void     Copy(sMatrix *dest, 
		  sMatrix *src){
      if(dest->ncols!=src->ncols ||
	 dest->nrows!=src->nrows)
	glip::Error((char *)"Incompatible matrices",
		   (char *)"Matrix::Copy");
      
      memcpy(dest->array[0], src->array[0], 
	     src->ncols*src->nrows*sizeof(float));
    }


    sMatrix *Identity(int n){
      sMatrix *I = NULL;
      int i;
      I = Create(n, n);
      for(i = 0; i < n; i++)
	I->array[i][i] = 1.0;
      return I;
    }
    

    bool Equal(sMatrix *A, sMatrix *B){
      int i,j;
      if(A->ncols != B->ncols ||
	 A->nrows != B->nrows){
	return false;
      }
      for(i=0; i<A->nrows; i++){
	for(j=0; j<A->ncols; j++){
	  if( A->array[i][j] != B->array[i][j] )
	    return false;
	}
      }
      return true;
    }


    bool Equal(sMatrix *A, sMatrix *B, float tolerance){
      int i,j;
      if(A->ncols != B->ncols ||
	 A->nrows != B->nrows){
	return false;
      }
      for(i=0; i<A->nrows; i++){
	for(j=0; j<A->ncols; j++){
	  if( fabsf(A->array[i][j] - B->array[i][j]) >= tolerance )
	    return false;
	}
      }
      return true;
    }
    

    void   Print(sMatrix *M){
      int x, y;
  
      printf("\n");
      for(y=0; y<M->nrows; y++){
	for(x=0; x<M->ncols; x++)
	  printf("%12.3f ", M->array[y][x]);
	printf("\n");
      }
      printf("\n");
    }


    void   Print(sMatrix *M, int min_width, int dec_places){
      char str[36];
      int x, y;
      sprintf(str, "%%%d.%df ", min_width, dec_places);
      //printf("str: %s\n", str);
      printf("\n");
      for(y=0; y<M->nrows; y++){
	for(x=0; x<M->ncols; x++)
	  printf(str, M->array[y][x]);
	  //printf("%12.3f ", M->array[y][x]);
	printf("\n");
      }
      printf("\n");
    }
    

    void   PrintDimension(sMatrix *M){
      printf("\n%d x %d\n", M->nrows, M->ncols);
    }
    

    float *SolveLinearSystem(sMatrix *A, float *b){
      sMatrix *B=NULL;
      float *x = NULL;
      float *ptmp;
      float maxmag, tmp;
      int maxrow;
      int i,j,k;
      float m;

      if(A->ncols!=A->nrows)
	glip::Error((char *)"Matrix dimension error",
		    (char *)"Matrix::SolveLinearSystem");
      
      x = (float *)calloc(A->ncols, sizeof(float));
      memcpy(x, b, A->ncols*sizeof(float));
      B = Clone(A);
      
      for(k=0; k<A->nrows; k++){
	//Find the maximum row:
	maxmag = fabsf(B->array[k][k]);
	maxrow = k;
	for(i=k+1; i<A->nrows; i++){
	  if( fabsf(B->array[i][k]) > maxmag ){
	    maxmag = fabsf(B->array[i][k]);
	    maxrow = i;
	  }
	}
	//Exchange the rows:
	ptmp = B->array[maxrow];
	B->array[maxrow] = B->array[k];
	B->array[k] = ptmp;

	tmp = x[maxrow];
	x[maxrow] = x[k];
	x[k] = tmp;
	
	m = B->array[k][k];

	B->array[k][k] = 1.0;
	for(j=k+1; j<A->ncols; j++){
	  B->array[k][j] = B->array[k][j]/m;
	}
	x[k] = x[k]/m;
	
	for(i=0; i<A->nrows; i++){
	  if(i!=k){
	    m = B->array[i][k];
	    
	    B->array[i][k] = 0.0;
	    for(j=k+1; j<A->ncols; j++){
	      B->array[i][j] = B->array[i][j] - m*B->array[k][j];
	    }
	    x[i] = x[i] - m*x[k];
	  }
	}
      }
      Destroy(&B);
      return x;
    }


    sMatrix *getPerspectiveTransform(float x[4], float y[4],
				     float u[4], float v[4]){
      sMatrix *A=NULL, *C=NULL;
      float b[8];
      float *c;
      int i;
      for(i = 0; i <= 3; i++)
	b[i] = u[i];
      for(i = 0; i <= 3; i++)
	b[4+i] = v[i];      

      A = Create(8, 8);
      A->array[0][0]=x[0]; A->array[0][1]=y[0]; A->array[0][2]=1.; A->array[0][6]=-x[0]*u[0]; A->array[0][7]=-y[0]*u[0];
      A->array[1][0]=x[1]; A->array[1][1]=y[1]; A->array[1][2]=1.; A->array[1][6]=-x[1]*u[1]; A->array[1][7]=-y[1]*u[1];
      A->array[2][0]=x[2]; A->array[2][1]=y[2]; A->array[2][2]=1.; A->array[2][6]=-x[2]*u[2]; A->array[2][7]=-y[2]*u[2];
      A->array[3][0]=x[3]; A->array[3][1]=y[3]; A->array[3][2]=1.; A->array[3][6]=-x[3]*u[3]; A->array[3][7]=-y[3]*u[3];

      A->array[4][3]=x[0]; A->array[4][4]=y[0]; A->array[4][5]=1.; A->array[4][6]=-x[0]*v[0]; A->array[4][7]=-y[0]*v[0];
      A->array[5][3]=x[1]; A->array[5][4]=y[1]; A->array[5][5]=1.; A->array[5][6]=-x[1]*v[1]; A->array[5][7]=-y[1]*v[1];
      A->array[6][3]=x[2]; A->array[6][4]=y[2]; A->array[6][5]=1.; A->array[6][6]=-x[2]*v[2]; A->array[6][7]=-y[2]*v[2];
      A->array[7][3]=x[3]; A->array[7][4]=y[3]; A->array[7][5]=1.; A->array[7][6]=-x[3]*v[3]; A->array[7][7]=-y[3]*v[3];
      
      c = SolveLinearSystem(A, b);
      C = Create(3, 3);

      C->array[0][0]=c[0]; C->array[0][1]=c[1]; C->array[0][2]=c[2];
      C->array[1][0]=c[3]; C->array[1][1]=c[4]; C->array[1][2]=c[5];
      C->array[2][0]=c[6]; C->array[2][1]=c[7]; C->array[2][2]=1.;
      
      free(c);
      Destroy(&A);
      return C;
    }
    

    sMatrix *Invert(sMatrix *A){
      //sMatrix *P=NULL; //permutation matrix
      sMatrix *B=NULL;
      sMatrix *I=NULL, *Inv = NULL;
      float *ptmp;
      float maxmag;
      int i,j,k,maxrow;
      float m;
      
      if(A->ncols!=A->nrows)
	glip::Error((char *)"Matrix dimension error",
		    (char *)"Matrix::Invert");
  
      //P = Identity(A->ncols);
      I = Identity(A->ncols);
      B = Clone(A);
      
      for(k=0; k<A->nrows; k++){
	//Find the maximum row:
	maxmag = fabsf(B->array[k][k]);
	maxrow = k;
	for(i=k+1; i<A->nrows; i++){
	  if( fabsf(B->array[i][k]) > maxmag ){
	    maxmag = fabsf(B->array[i][k]);
	    maxrow = i;
	  }
	}
	//Exchange the rows:
	ptmp = B->array[maxrow];
	B->array[maxrow] = B->array[k];
	B->array[k] = ptmp;

	ptmp = I->array[maxrow];
	I->array[maxrow] = I->array[k];
	I->array[k] = ptmp;

	//ptmp = P->array[maxrow];
	//P->array[maxrow] = P->array[k];
	//P->array[k] = ptmp;
	
	m = B->array[k][k];
	
	B->array[k][k] = 1.0;
	for(j=k+1; j<A->ncols; j++)
	  B->array[k][j] = B->array[k][j]/m;
	for(j=0; j<A->ncols; j++)
	  I->array[k][j] = I->array[k][j]/m;

	for(i=0; i<A->nrows; i++){
	  if(i!=k){
	    m = B->array[i][k];
	    
	    B->array[i][k] = 0.0;
	    for(j=k+1; j<A->ncols; j++)
	      B->array[i][j] = B->array[i][j] - m*B->array[k][j];
	    for(j=0; j<A->ncols; j++)
	      I->array[i][j] = I->array[i][j] - m*I->array[k][j];
	  }
	}
      }
      Inv = Create(A->ncols, A->nrows);
      for(i=0; i<A->nrows; i++){
	memcpy(Inv->array[i], I->array[i],
	       A->ncols*sizeof(float));
      }
      //Inv = Mult(I, P);
      //Destroy(&P);
      
      Destroy(&B);
      Destroy(&I);
      return Inv;
    }


    float        GetTrace(sMatrix *M){
      float sum;
      int i;
      
      if(M->ncols!=M->nrows)
	glip::Error((char *)"Matrix dimension error",
		   (char *)"Matrix::GetTrace");
      sum = 0.0;
      for(i=0; i<M->nrows; i++)
	sum += M->array[i][i];
      return sum;
    }
    
    
    sMatrix *Mult(sMatrix *A, 
		  sMatrix *B){
      sMatrix *M = NULL;
      int i,j,k;
      
      if(A->ncols!=B->nrows)
	glip::Error((char *)"Matrix dimension error",
		   (char *)"Matrix::Mult");
  
      M = Create(B->ncols, A->nrows);
      for(i=0; i<M->nrows; i++){
	for(j=0; j<M->ncols; j++){
	  M->array[i][j] = 0.0;
	  for (k=0; k<A->ncols; k++)
	    M->array[i][j] += A->array[i][k]*B->array[k][j];
	}
      }
      return(M);
    }
    
    
    sMatrix *MultByScalar(sMatrix *A, 
			  float k){
      sMatrix *M = NULL;
      int i,j;
      
      M = Create(A->ncols, A->nrows);
      for(i=0; i<M->nrows; i++){
	for(j=0; j<M->ncols; j++){
	  M->array[i][j] = k*A->array[i][j];
	}
      }
      return(M);
    }



    void Mult(sMatrix *M, float x[], float out[]){
      int i,j;
      for(i = 0; i < M->nrows; i++){
	out[i] = 0.0;
	for(j = 0; j < M->ncols; j++)
	  out[i] += M->array[i][j]*x[j];
      }
    }
    

    sMatrix *Sub(sMatrix *A, 
		 sMatrix *B){
      sMatrix *M = NULL;
      int i,j;
      
      if((A->ncols!=B->ncols)||(A->nrows!=B->nrows))
	glip::Error((char *)"Matrix dimension error",
		   (char *)"Matrix::Sub");
      M = Create(A->ncols, A->nrows);
      for(i=0; i<M->nrows; i++){
	for(j=0; j<M->ncols; j++){
	  M->array[i][j] = A->array[i][j] - B->array[i][j];
	}
      }
      return(M);
    }


    sMatrix *Add(sMatrix *A, 
		 sMatrix *B){
      sMatrix *M = NULL;
      int i,j;
      
      if((A->ncols!=B->ncols)||(A->nrows!=B->nrows))
	glip::Error((char *)"Matrix dimension error",
		   (char *)"Matrix::Add");
      
      M = Create(A->ncols, A->nrows);
      for(i=0; i<M->nrows; i++){
	for(j=0; j<M->ncols; j++){
	  M->array[i][j] = A->array[i][j] + B->array[i][j];
	}
      }
      return(M);
    }


    sMatrix *Transpose(sMatrix *A){
      sMatrix *M = NULL;
      int i,j;
      
      M = Create(A->nrows, A->ncols);
      for(i=0; i<M->nrows; i++){
	for(j=0; j<M->ncols; j++){
	  M->array[i][j] = A->array[j][i];
	}
      }
      return(M);
    }


    float ComputeDistanceL2(sMatrix *Y, 
			    sMatrix *X){
      sMatrix *A,*B,*R;
      float d;
  
      A = Sub(X, Y);
      B = Transpose(A);
      R = Mult(A, B);
      d = GetTrace(R);
      d = sqrtf(d);
      Destroy(&A);
      Destroy(&B);
      Destroy(&R);
      return (d);
    }

    
    sMatrix *Read(char *filename){
      sMatrix *M;
      char msg[512];
      int  ncols,nrows,size,n,p;
      double *daux=NULL;
      FILE *fp;

      fp = fopen(filename,"rb");
      if(fp == NULL){
	sprintf(msg,"Cannot open %s",filename);
	glip::Error((char *)msg,(char *)"Matrix::Read");
      }
      fread(&ncols, sizeof(int), 1, fp);
      fread(&nrows, sizeof(int), 1, fp);
      fread(&size,  sizeof(int), 1, fp);

      M = Create(ncols, nrows);
      n = ncols*nrows;
      if(size==sizeof(float))
	fread(M->array[0], sizeof(float), n, fp);
      else if(size==sizeof(double)){
	daux = AllocDoubleArray(n);
	fread(daux, sizeof(double), n, fp);
	for(p=0; p<n; p++)
	  M->array[0][p] = (float)daux[p];
	free(daux);
      }
      else
	glip::Error((char *)"Bad or corrupted file",
		   (char *)"Matrix::Read");
      fclose(fp);
      return M;
    }


    void    Write(sMatrix *M,
		  char *filename){
      char msg[512];
      int n,ncols,nrows,size;
      FILE *fp;
      
      fp = fopen(filename,"wb");
      if(fp == NULL){
	sprintf(msg,"Cannot open %s",filename);
	glip::Error((char *)msg,
		   (char *)"Matrix::Write");
      }
      size  = sizeof(float);
      ncols = M->ncols;
      nrows = M->nrows;
      n = ncols*nrows;
      fwrite(&ncols, sizeof(int),  1, fp);
      fwrite(&nrows, sizeof(int),  1, fp);
      fwrite(&size,  sizeof(int),  1, fp);
      fwrite(M->array[0], sizeof(float), n ,fp);
      fclose(fp);
    }


    bool   IsValidEntry(sMatrix *M,
			int i, int j){
      if((j >= 0)&&(j < M->ncols)&&
	 (i >= 0)&&(i < M->nrows))
	return(true);
      else
	return(false);
    }


    float   GetMaximumValue(sMatrix *M){
      float max;
      int p,n;
      
      n = M->ncols*M->nrows;
      max = M->array[0][0];
      for(p = 1; p < n; p++)
	if(M->array[0][p] > max)
	  max = M->array[0][p];
      return(max);
    }
    
    
    float    GetMinimumValue(sMatrix *M){
      float min;
      int p,n;
      n = M->ncols*M->nrows;
      min = M->array[0][0];
      for(p=1; p<n; p++)
	if(M->array[0][p] < min)
	  min = M->array[0][p];
      return(min);
    }


    void    Fill(sMatrix *M, float value){
      int p,n;
      n = M->ncols*M->nrows;
      for(p=0; p<n; p++)
	M->array[0][p] = value;
    }
    

    void    ChangeValue(sMatrix *M, 
			float old_value,
			float new_value){
      int p,n;
      
      n = M->ncols*M->nrows;
      for(p = 0; p < n; p++)
	if(M->array[0][p] == old_value)
	  M->array[0][p] = new_value;
    }



    // options: 0 (x) / 1 (y) / 2 (z)
    sMatrix* RotationMatrix3(int axis, 
			     float th){
      sMatrix *m;
      m = Create(4,4);
      if (axis==0) {
	m->array[0][0] = 1.0;  m->array[0][1] = 0.0;     m->array[0][2] = 0.0;      m->array[0][3] = 0.0;
	m->array[1][0] = 0.0;  m->array[1][1] = cos(th); m->array[1][2] = -sin(th); m->array[1][3] = 0.0;
	m->array[2][0] = 0.0;  m->array[2][1] = sin(th); m->array[2][2] = cos(th);  m->array[2][3] = 0.0;
	m->array[3][0] = 0.0;  m->array[3][1] = 0.0;     m->array[3][2] = 0.0;      m->array[3][3] = 1.0;
      }
      if (axis==1) {
	m->array[0][0] = cos(th);  m->array[0][1] = 0.0;  m->array[0][2] = sin(th);  m->array[0][3] = 0.0;
	m->array[1][0] = 0.0;      m->array[1][1] = 1;    m->array[1][2] = 0.0;      m->array[1][3] = 0.0;
	m->array[2][0] = -sin(th); m->array[2][1] = 0;    m->array[2][2] = cos(th);  m->array[2][3] = 0.0;
	m->array[3][0] = 0.0;      m->array[3][1] = 0.0;  m->array[3][2] = 0.0;      m->array[3][3] = 1.0;
	
      }
      if (axis==2) {
	m->array[0][0] = cos(th); m->array[0][1] = -sin(th); m->array[0][2] = 0.0;  m->array[0][3] = 0.0;
	m->array[1][0] = sin(th); m->array[1][1] = cos(th);  m->array[1][2] = 0.0;  m->array[1][3] = 0.0;
	m->array[2][0] = 0.0;     m->array[2][1] = 0.0;      m->array[2][2] = 1.0;  m->array[2][3] = 0.0;
	m->array[3][0] = 0.0;     m->array[3][1] = 0.0;      m->array[3][2] = 0.0;  m->array[3][3] = 1.0;
      }
      return m;
    }



    sMatrix* TranslationMatrix3(float dx, float dy, float dz){
      sMatrix *m;
      m = Create(4,4);
      m->array[0][0] = 1.0;  m->array[0][1] = 0.0;  m->array[0][2] = 0.0;  m->array[0][3] = dx;
      m->array[1][0] = 0.0;  m->array[1][1] = 1.0;  m->array[1][2] = 0.0;  m->array[1][3] = dy;
      m->array[2][0] = 0.0;  m->array[2][1] = 0.0;  m->array[2][2] = 1.0;  m->array[2][3] = dz;
      m->array[3][0] = 0.0;  m->array[3][1] = 0.0;  m->array[3][2] = 0.0;  m->array[3][3] = 1.0;
      return m;
    }


    sMatrix* TransformVoxel(sMatrix *m, Voxel v){
      sMatrix *vm,*res;
      vm = Create(1,4);
      vm->array[0][0]=v.c.x;
      vm->array[1][0]=v.c.y;
      vm->array[2][0]=v.c.z;
      vm->array[3][0]=1.0;
      res=Mult(m,vm);
      Destroy(&vm);
      return res;
    }


    sMatrix *Decode(sImage32 *bin, int ncols, int nrows){
      sMatrix *C;
      int i,j,u,v,sum,u1,u2,v1,v2;
      float dx,dy;
      dx = bin->ncols/ncols;
      dy = bin->nrows/nrows;
      C = Create(ncols, nrows);
      for(i = 0; i < nrows; i++){
	v1 = ROUND(i*dy);
	v2 = MIN(ROUND((i+1)*dy), bin->nrows);
	for(j = 0; j < ncols; j++){
	  u1 = ROUND(j*dx);
	  u2 = MIN(ROUND((j+1)*dx), bin->ncols);

	  sum = 0;
	  for(v = v1; v < v2; v++){
	    for(u = u1; u < u2; u++){
	      if(bin->array[v][u] > 0)
		sum += 1;
	    }
	  }
	  if( sum > ((v2-v1)*(u2-u1))/2.0 )
	    C->array[i][j] = 1.0;
	}
      }
      return C;
    }
    

  } //end Matrix namespace
} //end glip namespace

