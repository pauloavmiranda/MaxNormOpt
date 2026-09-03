
//#define APPDEBUG  1

#include "glip.h"

#define ACC_DICE    0
#define ACC_JACCARD 1  //Intersection over Union (IoU)


#ifdef APPDEBUG
size_t mem_size_OIFT_1 = 0;
size_t mem_size_OIFT_2 = 0;
size_t mem_size_ORFC = 0;
#endif


struct sEdge {
  int s;
  int t;
  int index; //ImageGraph data structure index
};


struct sArrayEdges {
  int nedges;
  struct sEdge *E;
};


struct sNeighboringEdges {
  int nnodes;
  int *nnedges; //number of neighboring edges
  int **N;
};


struct sGraphAuxiliaryData {
  struct sArrayEdges AE;
  struct sNeighboringEdges NE;
};


enum AtomType {
  UNARY,
  BINARY,
  REMOVED
};


struct sAtom {
  enum AtomType type;
  int vertex_edge;
  int ls;
  int lt;
};


struct sArrayAtoms {
  int natoms;
  struct sAtom *A;
  int *Phi;
};



struct sArrayBucketsAtoms {
  int nnodes;
  //matrix A_vertex[s][i] = a, where s is a node, i in [0,1], and a is an atom index.
  int **A_vertex; 
  //array n_vertex[s] = n, where s is a node, and n is the number of atoms at vertex s.
  int *n_vertex; 

  int nedges;
  //matrix A_edge[e][i] = a, where e is an edge, i in [0,3], and a is an atom index.
  int **A_edge;
  //array n_edge[e] = n, where e is an edge, and n is the number of atoms at edge e.
  int *n_edge;
};



glip::sImage32 *ReadAnyImage(char *file){
  glip::sImage32 *img;
  char command[512];
  int s;
  s = strlen(file);
  if(strcasecmp(&file[s-3], "pgm") == 0){
    img = glip::Image32::Read(file);
  }
  else{
    sprintf(command, "convert %s image_tmp.pgm", file);
    system(command);
    img = glip::Image32::Read("image_tmp.pgm");
    system("rm image_tmp.pgm");
  }
  return img;
}


glip::sCImage *ReadAnyCImage(char *file){
  glip::sCImage *cimg;
  char command[512];
  int s;

  s = strlen(file);
  if(strcasecmp(&file[s-3], "ppm") == 0){
    cimg = glip::CImage::Read(file);
  }
  else{
    sprintf(command, "convert %s cimage_tmp.ppm", file);
    system(command);
    cimg = glip::CImage::Read("cimage_tmp.ppm");
    system("rm cimage_tmp.ppm");
  }
  return cimg;
}


struct sGraphAuxiliaryData *CreateGraphAuxiliaryData(glip::sImageGraph *graph){
  struct sGraphAuxiliaryData *GAD;
  glip::sAdjRel *A;
  int u_x,u_y,v_x,v_y, s,t,i;
  int ncols, nrows;
  int nne, nedges, e, nnodes;
  A = graph->A;
  ncols = graph->ncols;
  nrows = graph->nrows;
  nnodes = ncols*nrows;
  GAD = (struct sGraphAuxiliaryData *)calloc(1, sizeof(struct sGraphAuxiliaryData));
  if(GAD == NULL){
    printf("Error: CreateGraphAuxiliaryData\n");
    exit(1);
  }

  //allocating memory for graph edge array.
  nedges = glip::ImageGraph::GetNumberOfEdges(graph);
  GAD->AE.nedges = nedges;
  GAD->AE.E = (struct sEdge *)calloc(nedges, sizeof(struct sEdge));
  if(GAD->AE.E == NULL){
    printf("Error: CreateGraphAuxiliaryData\n");
    exit(1);
  }
  
  //allocating memory to store the indices of all neighboring edges of each pixel.
  GAD->NE.nnodes = nnodes;
  GAD->NE.nnedges = (int *)calloc(nnodes, sizeof(int));
  GAD->NE.N = (int **)calloc(nnodes, sizeof(int *));
  if(GAD->NE.nnedges == NULL || GAD->NE.N == NULL){
    printf("Error: CreateGraphAuxiliaryData\n");
    exit(1);
  }
  for(s = 0; s < nnodes; s++){
    GAD->NE.N[s] = (int *)calloc(A->n-1, sizeof(int));
    if(GAD->NE.N[s] == NULL){
      printf("Error: CreateGraphAuxiliaryData\n");
      exit(1);
    }
  }

  //filling the data structures.
  e = 0;
  for(s = 0; s < nnodes; s++){
    u_x = s%ncols;
    u_y = s/ncols;
    for(i=1; i<A->n; i++){
      v_x = u_x + A->dx[i];
      v_y = u_y + A->dy[i];
      if(v_x >= 0 && v_x < ncols &&
	 v_y >= 0 && v_y < nrows){
	t = v_x + ncols*v_y;
	if(s < t){
	  GAD->AE.E[e].s = s;
	  GAD->AE.E[e].t = t;
	  GAD->AE.E[e].index = i;
	  
	  nne = GAD->NE.nnedges[s];
	  GAD->NE.N[s][nne] = e;
	  GAD->NE.nnedges[s] += 1;

	  nne = GAD->NE.nnedges[t];
	  GAD->NE.N[t][nne] = e;
	  GAD->NE.nnedges[t] += 1;
	  
	  e++;
	}
      }
    }
  }
  return GAD;
}


void DestroyGraphAuxiliaryData(struct sGraphAuxiliaryData **pGAD){
  struct sGraphAuxiliaryData *GAD;
  int s;
  if(pGAD != NULL){
    GAD = *pGAD;
    if(GAD != NULL){
      for(s = 0; s < GAD->NE.nnodes; s++){
	if(GAD->NE.N[s] != NULL)
	  free(GAD->NE.N[s]);
      }
      if(GAD->AE.E != NULL)
	free(GAD->AE.E);
      if(GAD->NE.nnedges != NULL)
	free(GAD->NE.nnedges);
      if(GAD->NE.N != NULL)
	free(GAD->NE.N);
      if(GAD != NULL)
	free(GAD);
      *pGAD = NULL;
    }
  }
}


size_t GetSizeOf(struct sGraphAuxiliaryData *GAD){
  size_t size = 0;
  int s;
  if(GAD == NULL)
    return 0;
  size = sizeof(struct sGraphAuxiliaryData); //GAD
  size += GAD->NE.nnodes*sizeof(int*); //GAD->NE.N
  size += GAD->NE.nnodes*sizeof(int); //GAD->NE.nnedges
  size += GAD->AE.nedges*sizeof(struct sEdge); //GAD->AE.E
  for(s = 0; s < GAD->NE.nnodes; s++){
    size += GAD->NE.nnedges[s]*sizeof(int); //GAD->NE.N[s]
  }
  return size;
}


struct sArrayAtoms *CreateOrderedArrayAtoms(glip::sImageGraph *graph,
					    struct sGraphAuxiliaryData *GAD,
					    int *S_bkg,
					    int *S_obj){
  struct sArrayAtoms *AA = NULL, *AA_sorted = NULL;
  glip::sPQueue32 *Q=NULL;
  glip::sAdjRel *A;
  int a,s,t,e,i,j,k,w;
  int *i_inv;
  A = graph->A;
  AA_sorted = (struct sArrayAtoms *)calloc(1, sizeof(struct sArrayAtoms));
  AA        = (struct sArrayAtoms *)calloc(1, sizeof(struct sArrayAtoms));
  if(AA == NULL || AA_sorted == NULL){
    printf("Error: CreateOrderedArrayAtoms\n");
    exit(1);
  }

  AA->natoms = 2*GAD->NE.nnodes + 4*GAD->AE.nedges;
  AA->A = (struct sAtom *)calloc(AA->natoms, sizeof(struct sAtom));
  AA->Phi = (int *)calloc(AA->natoms, sizeof(int));
  if(AA->A == NULL || AA->Phi == NULL){
    printf("Error: CreateOrderedArrayAtoms\n");
    exit(1);
  }

  AA_sorted->natoms = AA->natoms;
  AA_sorted->A = (struct sAtom *)calloc(AA->natoms, sizeof(struct sAtom));
  AA_sorted->Phi = (int *)calloc(AA->natoms, sizeof(int));
  if(AA_sorted->A == NULL || AA_sorted->Phi == NULL){
    printf("Error: CreateOrderedArrayAtoms\n");
    exit(1);
  }

  //filling the array with atoms derived from nodes.
  a = 0;
  for(s = 0; s < GAD->NE.nnodes; s++){
    AA->A[a].type = UNARY;
    AA->A[a].vertex_edge = s;
    AA->A[a].ls = 0;
    AA->A[a].lt = NIL;
    a++;

    AA->A[a].type = UNARY;
    AA->A[a].vertex_edge = s;
    AA->A[a].ls = 1;
    AA->A[a].lt = NIL;
    a++;
  }
  for(k = 1; k <= S_bkg[0]; k++){
    s = S_bkg[k];
    AA->Phi[2*s]   = 0;
    AA->Phi[2*s+1] = graph->Wmax+1;
  }
  for(k = 1; k <= S_obj[0]; k++){
    s = S_obj[k];
    AA->Phi[2*s]   = graph->Wmax+1;
    AA->Phi[2*s+1] = 0;
  }

  //filling the array with atoms derived from edges.
  i_inv = glip::AdjRel::InverseIndexes(A);
  for(e = 0; e < GAD->AE.nedges; e++){
    AA->Phi[a] = 0;
    AA->A[a].type = BINARY;
    AA->A[a].vertex_edge = e;
    AA->A[a].ls = 0;
    AA->A[a].lt = 0;
    a++;
    
    AA->Phi[a] = 0;
    AA->A[a].type = BINARY;
    AA->A[a].vertex_edge = e;
    AA->A[a].ls = 1;
    AA->A[a].lt = 1;
    a++;

    s = GAD->AE.E[e].s;
    i = GAD->AE.E[e].index;
    w = (graph->n_link[s])[i];
    AA->Phi[a] = w;
    AA->A[a].type = BINARY;
    AA->A[a].vertex_edge = e;
    AA->A[a].ls = 1;
    AA->A[a].lt = 0;
    a++;

    t = GAD->AE.E[e].t;
    j = i_inv[i];
    w = (graph->n_link[t])[j];
    AA->Phi[a] = w;
    AA->A[a].type = BINARY;
    AA->A[a].vertex_edge = e;
    AA->A[a].ls = 0;
    AA->A[a].lt = 1;
    a++;
  }
  free(i_inv);

  //The ordering of atoms in decreasing order of Phi.
  /*
  Q = glip::PQueue32::Create(graph->Wmax+2, AA->natoms, AA->Phi);
  for(a = 0; a < AA->natoms; a++){
    if(AA->A[a].type == BINARY && AA->A[a].ls == AA->A[a].lt)
      glip::PQueue32::FastInsertElem(Q, a); //Insert elem as last.
    else
      glip::PQueue32::FastInsertElemAsFirst(Q, a);
  }

  k = 0;
  while(!glip::PQueue32::IsEmpty(Q)) {
    a = glip::PQueue32::FastRemoveMaxFIFO(Q);
    AA_sorted->Phi[k] = AA->Phi[a];
    AA_sorted->A[k] = AA->A[a];
    k++;
  }
  glip::PQueue32::Destroy(&Q);
  */
  //----------------------------------
  int *count, *start, *end;
  count = (int *)calloc(graph->Wmax+2, sizeof(int));
  start = (int *)calloc(graph->Wmax+2, sizeof(int));
  end   = (int *)calloc(graph->Wmax+2, sizeof(int));
  for(a = 0; a < AA->natoms; a++){
    count[ AA->Phi[a] ] += 1;
  }
  start[0] = 0;
  end[0] = count[0]-1;
  for(i = 1; i <= graph->Wmax+1; i++){
    start[i] = start[i-1] + count[i-1];
    end[i] = start[i] + count[i]-1;
  }
  for(a = AA->natoms-1; a >= 0; a--){
    w = AA->Phi[a];
    if(AA->A[a].type == BINARY && AA->A[a].ls == AA->A[a].lt){ //Insert elem as last.
      k = AA->natoms-1 - start[w];
      AA_sorted->Phi[k] = AA->Phi[a];
      AA_sorted->A[k] = AA->A[a];
      start[w] += 1;
    }
    else{
      k = AA->natoms-1 - end[w];
      AA_sorted->Phi[k] = AA->Phi[a];
      AA_sorted->A[k] = AA->A[a];
      end[w] -= 1;
    }
  }
  free(count);
  free(start);
  free(end);
  //----------------------------------
  free(AA->Phi);
  free(AA->A);
  free(AA);
 
  return AA_sorted;
}



void DestroyArrayAtoms(struct sArrayAtoms **pAA){
  struct sArrayAtoms *AA;
  if(pAA != NULL){
    AA = *pAA;
    if(AA != NULL){
      if(AA->Phi != NULL)
	free(AA->Phi);
      if(AA->A != NULL)
	free(AA->A);
      if(AA != NULL)
	free(AA);
      *pAA = NULL;
    }
  }
}



size_t GetSizeOf(struct sArrayAtoms *AA){
  size_t size = 0;
  if(AA == NULL)
    return 0;
  size = sizeof(struct sArrayAtoms); //AA
  size += AA->natoms*sizeof(struct sAtom); //AA->A
  size += AA->natoms*sizeof(int); //AA->Phi
  return size;
}



int *LoadSeedsTxt(char *filename,
		  int label,
		  int ncols, int nrows){
  int *S = NULL;
  FILE *fp;
  int i,j,nseeds;
  int x,y,id,lb,p;
  fp = fopen(filename, "r");
  if (fp == NULL) {
    printf("Error reading seeds.\n");
    exit(1);
  }
  fscanf(fp, " %d", &nseeds);
  S = (int *)calloc((nseeds + 1), sizeof(int));
  S[0] = nseeds;
  j = 0;
  for (i = 0; i < nseeds; i++) {
    fscanf(fp, " %d %d %d %d", &x, &y, &id, &lb);
    if(x >= 0 && x < ncols &&
       y >= 0 && y < nrows &&
       lb == label){
      p = x + y*ncols;
      j++;
      S[j] = p;
    }
  }
  S[0] = j;
  fclose(fp);
  S = (int *)realloc(S, (j+1)*sizeof(int));
  return S;
}



struct sArrayBucketsAtoms *CreateArrayBucketsAtoms(struct sArrayAtoms *AA,
						   int nnodes,
						   int nedges){
  struct sArrayBucketsAtoms *ABA = NULL;
  int s,e,a,i;
  ABA = (struct sArrayBucketsAtoms *)calloc(1, sizeof(struct sArrayBucketsAtoms));
  if(ABA == NULL){
    printf("Error: CreateArrayBucketsAtoms\n");
    exit(1);
  }
  ABA->A_vertex = (int **)calloc(nnodes, sizeof(int *));
  ABA->A_edge   = (int **)calloc(nedges, sizeof(int *));
  ABA->n_vertex = (int *)calloc(nnodes, sizeof(int));
  ABA->n_edge   = (int *)calloc(nedges, sizeof(int));
  if(ABA->A_vertex == NULL || ABA->A_edge == NULL ||
     ABA->n_vertex == NULL || ABA->n_edge == NULL){
    printf("Error: CreateArrayBucketsAtoms\n");
    exit(1);
  }
  for(s = 0; s < nnodes; s++){
    ABA->A_vertex[s] = (int *)calloc(2, sizeof(int));
    if(ABA->A_vertex[s] == NULL){
      printf("Error: CreateArrayBucketsAtoms\n");
      exit(1);
    }
  }
  for(e = 0; e < nedges; e++){
    ABA->A_edge[e] = (int *)calloc(4, sizeof(int));
    if(ABA->A_edge[e] == NULL){
      printf("Error: CreateArrayBucketsAtoms\n");
      exit(1);
    }
  }

  ABA->nnodes = nnodes;
  ABA->nedges = nedges;
  for(a = 0; a < AA->natoms; a++){
    //creates the initial array A_vertex of buckets of atoms, indexed by nodes.
    if(AA->A[a].type == UNARY){
      s = AA->A[a].vertex_edge;
      i = ABA->n_vertex[s];
      ABA->A_vertex[s][i] = a;
      ABA->n_vertex[s] += 1;
    }
    //creates the initial array A_edge of buckets of atoms, indexed by edges.
    else{ //BINARY
      e = AA->A[a].vertex_edge;
      i = ABA->n_edge[e];
      ABA->A_edge[e][i] = a;
      ABA->n_edge[e] += 1;
    }
  }
  return ABA;
}



void DestroyArrayBucketsAtoms(struct sArrayBucketsAtoms **pABA){
  struct sArrayBucketsAtoms *ABA = NULL;
  int s,e;
  if(pABA != NULL){
    ABA = *pABA;
    if(ABA != NULL){
      for(s = 0; s < ABA->nnodes; s++)
	if(ABA->A_vertex[s] != NULL)
	  free(ABA->A_vertex[s]);
      for(e = 0; e < ABA->nedges; e++)
	if(ABA->A_edge[e] != NULL)
	  free(ABA->A_edge[e]);
      if(ABA->A_vertex != NULL)
	free(ABA->A_vertex);
      if(ABA->A_edge != NULL)
	free(ABA->A_edge);
      if(ABA->n_vertex != NULL)
	free(ABA->n_vertex);
      if(ABA->n_edge != NULL)
	free(ABA->n_edge);
      if(ABA != NULL)
	free(ABA);
      *pABA = NULL;
    }
  }  
}


size_t GetSizeOf(struct sArrayBucketsAtoms *ABA){
  size_t size = 0;
  int e,s;
  if(ABA == NULL)
    return 0;
  size = sizeof(struct sArrayBucketsAtoms); //ABA
  size += ABA->nedges*sizeof(int); //ABA->n_edge
  size += ABA->nnodes*sizeof(int); //ABA->n_vertex
  size += ABA->nedges*sizeof(int*); //ABA->A_edge
  size += ABA->nnodes*sizeof(int*); //ABA->A_vertex
  for(e = 0; e < ABA->nedges; e++)
    size += 4*sizeof(int); //ABA->A_edge[e]
  for(s = 0; s < ABA->nnodes; s++)
    size += 2*sizeof(int); //ABA->A_vertex[s]
  return size;
}


inline bool remove_incompatible_unary_atom(int s,
					   int e, //adjacent edge with atoms removed.
					   struct sGraphAuxiliaryData *GAD,
					   struct sArrayBucketsAtoms *ABA,
					   struct sArrayAtoms *H){
  bool v_ls[2];
  int k,a,lb,n;
  v_ls[0] = false; //Is label 0 for s valid?
  v_ls[1] = false; //Is label 1 for s valid?
  for(k = 0; k < ABA->n_edge[e]; k++){
    a = ABA->A_edge[e][k];
    if(GAD->AE.E[e].s == s)
      lb = H->A[a].ls;
    else //GAD->AE.E[e].t == s
      lb = H->A[a].lt;
    v_ls[lb] = true;
  }
#ifdef APPDEBUG
  if(v_ls[0] == false && v_ls[1] == false)
    printf("Error: there is no valid unary atom.\n");
#endif
  for(k = 0; k < ABA->n_vertex[s]; k++){
    a = ABA->A_vertex[s][k];
    lb = H->A[a].ls;
    if(v_ls[lb] == false){
      n = ABA->n_vertex[s];
      ABA->A_vertex[s][k] = ABA->A_vertex[s][n-1];
      ABA->n_vertex[s] -= 1;
      H->A[a].type = REMOVED;
      return true;
    }
  }
  return false;
}


inline bool remove_incompatible_binary_atom(int e,
					    int s, //adjacent node with atoms removed.
					    struct sGraphAuxiliaryData *GAD,
					    struct sArrayBucketsAtoms *ABA,
					    struct sArrayAtoms *H){
  int a,inv_lb,lb,k,n;
  bool removed = false;
#ifdef APPDEBUG
  if(ABA->n_vertex[s] != 1)
    printf("Error: wrong amount of atoms remaining.\n");
#endif
  a = ABA->A_vertex[s][0];
  inv_lb = 1 - H->A[a].ls;
  for(k = 0; k < ABA->n_edge[e]; k++){
    a = ABA->A_edge[e][k];
    if(GAD->AE.E[e].s == s)
      lb = H->A[a].ls;
    else //GAD->AE.E[e].t == s
      lb = H->A[a].lt;

    if(lb == inv_lb){
      n = ABA->n_edge[e];
      ABA->A_edge[e][k] = ABA->A_edge[e][n-1];
      ABA->n_edge[e] -= 1;
      k--;
      H->A[a].type = REMOVED;
      removed = true;
    }
  }
  return removed;
}




//=================ORFC=================================
glip::sImage32 *IFT_fmin_transpose(glip::sImageGraph *sg,
				   int *S){
  glip::sPQueue32 *Q = NULL; 
  glip::sImage32 *V;
  int i,j,p,q,n, edge,tmp;
  glip::Pixel u,v;
  glip::sAdjRel *A;
  int *i_inv;

  A = sg->A;
  i_inv = glip::AdjRel::InverseIndexes(A);
  n = sg->ncols*sg->nrows;
  V = glip::Image32::Create(sg->ncols, sg->nrows);
  Q = glip::PQueue32::Create(sg->Wmax+2, n, V->data);

#ifdef APPDEBUG
  mem_size_ORFC = 0;
  mem_size_ORFC += sizeof(int)*A->n; //i_inv
  mem_size_ORFC += glip::Image32::GetSizeOf(V);
  mem_size_ORFC += glip::PQueue32::GetSizeOf(Q);
#endif
  
  glip::Image32::Set(V, INT_MIN);
  for(i = 1; i <= S[0]; i++){
    V->data[S[i]] = sg->Wmax+1;
    glip::PQueue32::FastInsertElem(Q, S[i]);
  }
  while(!glip::PQueue32::IsEmpty(Q)){
    p = glip::PQueue32::FastRemoveMaxFIFO(Q);
    u.x = p%V->ncols; 
    u.y = p/V->ncols; 
    for(i=1; i<A->n; i++){
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if(v.x >= 0 && v.x < V->ncols &&
	 v.y >= 0 && v.y < V->nrows){
	q = v.x + V->ncols*v.y;
	if(Q->L.elem[q].color != BLACK){
	  j = i_inv[i];
	  edge = (sg->n_link[q])[j];
	  tmp  = MIN(V->data[p], edge);
	  if(tmp > V->data[q]){
	    if(Q->L.elem[q].color == GRAY)
	      glip::PQueue32::FastRemoveElem(Q, q);
	    V->data[q] = tmp;
	    glip::PQueue32::FastInsertElem(Q, q);
	  }
	}
      }
    }
  }
  glip::PQueue32::Destroy(&Q);
  free(i_inv);
  return V;
}



void ORFC_MinMax(glip::sImageGraph *sg,
		 int *S_bkg,
		 int *S_obj,
		 glip::sImage32 *label){
  glip::sImage32 *value_e;
  glip::sPQueue32 *QS=NULL;
  glip::sQueue *Q=NULL;
  glip::sAdjRel *A = sg->A;
  int p,q,j,i,k,energy,nsi;
  glip::Pixel u,v;
  int *s_energy,*s_pixel;

  Q = glip::Queue::Create(sg->ncols*sg->nrows);
  nsi = S_obj[0];
  s_energy = glip::AllocIntArray(nsi);
  s_pixel  = glip::AllocIntArray(nsi);
  QS = glip::PQueue32::Create(sg->Wmax+2, nsi, s_energy);
  value_e = IFT_fmin_transpose(sg, S_bkg);
  glip::Image32::Set(label, 0);

#ifdef APPDEBUG
  size_t size_data = 0;
  size_data += glip::Queue::GetSizeOf(Q);
  size_data += sizeof(int)*nsi; //s_energy
  size_data += sizeof(int)*nsi; //s_pixel
  size_data += glip::PQueue32::GetSizeOf(QS);  
  size_data += glip::Image32::GetSizeOf(value_e);
  mem_size_ORFC = MAX(size_data, mem_size_ORFC);
#endif
  
  k = 0;
  for(i = 1; i <= S_obj[0]; i++){
    p = S_obj[i];
    s_energy[k] = value_e->data[p];
    s_pixel[k] = p;
    glip::PQueue32::InsertElem(&QS, k);
    k++;
  }
  
  while(!glip::PQueue32::IsEmpty(QS)){
    j = glip::PQueue32::RemoveMinFIFO(QS);
    p = s_pixel[j];
    energy = s_energy[j];
    
    if(label->data[p] != 0) continue;
    
    glip::Queue::Reset(Q);
    glip::Queue::Push(Q, p);
    label->data[p] = 1;
    while(!glip::Queue::IsEmpty(Q)){
      p = glip::Queue::Pop(Q);
      u.x = p%sg->ncols;
      u.y = p/sg->ncols;      
      for(i = 1; i < A->n; i++){
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	if(v.x >= 0 && v.x < sg->ncols &&
	   v.y >= 0 && v.y < sg->nrows){	    
	  q = v.x + v.y*sg->ncols;
	  if(energy < (sg->n_link[p])[i] &&
	     label->data[q] == 0){
	    label->data[q] = 1;
	    glip::Queue::Push(Q, q);
	  }
	}
      }
    }
  }
  glip::Queue::Destroy(&Q);
  glip::Image32::Destroy(&value_e);
  glip::PQueue32::Destroy(&QS);
  glip::FreeIntArray(&s_energy);
  glip::FreeIntArray(&s_pixel);      
}


//=================OIFT=================================

void OIFT_MinMax_1st(glip::sImageGraph *sg,
		     int *S_bkg,
		     int *S_obj,
		     glip::sImage32 *label){
  glip::sImage32 *value;
  glip::sPQueue32 *Q=NULL;
  int i,j,p,q,n,w,tmp;
  int u_x,u_y,v_x,v_y;
  glip::sAdjRel *A;
  int *i_inv;
  
  value = glip::Image32::Create(label);  
  n = label->n;
  Q = glip::PQueue32::Create(sg->Wmax*2+1+2,n,value->data);
  A = sg->A;
  i_inv = glip::AdjRel::InverseIndexes(A);

#ifdef APPDEBUG
  mem_size_OIFT_1 = 0;
  mem_size_OIFT_1 += glip::Image32::GetSizeOf(value);
  mem_size_OIFT_1 += glip::PQueue32::GetSizeOf(Q);
  mem_size_OIFT_1 += sizeof(int)*A->n; //i_inv
#endif
  
  glip::Image32::Set(label, 0);
  for(i=1; i<=S_bkg[0]; i++){
    p = S_bkg[i];
    value->data[p] = sg->Wmax*2+1+1;
    label->data[p] = 0;
    glip::PQueue32::FastInsertElem(Q, p);
  }
  for(i=1; i<=S_obj[0]; i++){
    p = S_obj[i];
    value->data[p] = sg->Wmax*2+1+1;
    label->data[p] = 1;
    glip::PQueue32::FastInsertElem(Q, p);
  }
  
  while(!glip::PQueue32::IsEmpty(Q)) {
    p = glip::PQueue32::FastRemoveMaxFIFO(Q);
    u_x = p%sg->ncols;
    u_y = p/sg->ncols;
    for(i=1; i<A->n; i++){
      v_x = u_x + A->dx[i];
      v_y = u_y + A->dy[i];
      if(v_x >= 0 && v_x < sg->ncols &&
	 v_y >= 0 && v_y < sg->nrows){
	q = v_x + sg->ncols*v_y;
	if(Q->L.elem[q].color != BLACK){
	  if(label->data[p]==0){
	    j = i_inv[i];
	    w = (sg->n_link[q])[j]*2 + 1;
	  }
	  else
	    w = (sg->n_link[p])[i]*2;

	  tmp = MIN(value->data[p], w);
	  if(tmp > value->data[q]){
	    if(Q->L.elem[q].color == GRAY)
	      glip::PQueue32::FastRemoveElem(Q, q);
	    value->data[q] = tmp;
	    label->data[q] = label->data[p];
	    glip::PQueue32::FastInsertElem(Q, q);
	  }
	}
      }
    }
  }
  glip::Image32::Destroy(&value);
  glip::PQueue32::Destroy(&Q);
  free(i_inv);
}



void OIFT_MinMax_2nd(glip::sImageGraph *sg,
		     int *S_bkg,
		     int *S_obj,
		     glip::sImage32 *label){
  glip::sImage32 *value;
  glip::sPQueue32 *Q=NULL;
  int i,j,p,q,n,w;
  int u_x,u_y,v_x,v_y;
  glip::sAdjRel *A;
  int *i_inv;
  
  value = glip::Image32::Create(label);  
  n = label->n;
  Q = glip::PQueue32::Create(sg->Wmax+2,n,value->data);
  A = sg->A;
  i_inv = glip::AdjRel::InverseIndexes(A);

#ifdef APPDEBUG
  mem_size_OIFT_2 = 0;
  mem_size_OIFT_2 += glip::Image32::GetSizeOf(value);
  mem_size_OIFT_2 += glip::PQueue32::GetSizeOf(Q);
  mem_size_OIFT_2 += sizeof(int)*A->n; //i_inv
#endif
  
  glip::Image32::Set(label, 0);
  for(i=1; i<=S_bkg[0]; i++){
    p = S_bkg[i];
    value->data[p] = sg->Wmax+1;
    label->data[p] = 0;
    glip::PQueue32::FastInsertElem(Q, p);
  }
  for(i=1; i<=S_obj[0]; i++){
    p = S_obj[i];
    value->data[p] = sg->Wmax+1;
    label->data[p] = 1;
    glip::PQueue32::FastInsertElem(Q, p);
  }
  
  while(!glip::PQueue32::IsEmpty(Q)) {
    p = glip::PQueue32::FastRemoveMaxFIFO(Q);
    u_x = p%sg->ncols;
    u_y = p/sg->ncols;
    for(i=1; i<A->n; i++){
      v_x = u_x + A->dx[i];
      v_y = u_y + A->dy[i];
      if(v_x >= 0 && v_x < sg->ncols &&
	 v_y >= 0 && v_y < sg->nrows){
	q = v_x + sg->ncols*v_y;
	if(Q->L.elem[q].color != BLACK){
	  if(label->data[p]==0){
	    j = i_inv[i];
	    w = (sg->n_link[q])[j];
	  }
	  else
	    w = (sg->n_link[p])[i];
	  
	  if(w > value->data[q]){
	    if(Q->L.elem[q].color == GRAY)
	      glip::PQueue32::FastRemoveElem(Q, q);
	    value->data[q] = w;
	    label->data[q] = label->data[p];
	    glip::PQueue32::FastInsertElem(Q, q);
	  }
	}
      }
    }
  }
  glip::Image32::Destroy(&value);
  glip::PQueue32::Destroy(&Q);
  free(i_inv);
}


void PrintTimeAcc(struct timeval tic,
		  struct timeval toc,
		  glip::sImage32 *label,
		  glip::sImage32 *gtruth,
		  char *message,
		  FILE *report,
		  int acc_type,
		  bool last){
  double totaltime, acc = 0.0;
  totaltime = ((toc.tv_sec-tic.tv_sec)*1000.0 + 
	       (toc.tv_usec-tic.tv_usec)*0.001);
  printf("%s\n\t time: %f ms\n", message, totaltime);
  fprintf(report, "%lf;", totaltime);
  if(gtruth != NULL){
    if(acc_type == ACC_DICE){
      acc = glip::Image32::DiceSimilarity(gtruth, label);
      printf("\t dice: %f\n", acc);
    }
    else{ //acc_type == ACC_JACCARD
      acc = glip::Image32::JaccardSimilarity(gtruth, label);
      printf("\t jaccard: %f\n", acc);
    }
  }
  if(last)
    fprintf(report, "%lf\n", acc);
  else
    fprintf(report, "%lf;", acc);
}


int main(int argc, char **argv){
  glip::sImage32 *img = NULL, *tmp = NULL;
  glip::sCImage *cimg = NULL, *ctmp = NULL;
  glip::sImage32 *label = NULL, *gtruth = NULL;
  glip::sImageGraph *graph = NULL;
  glip::sQueue *K_vertex = NULL; //*K_edge = NULL;
  char *IsInK = NULL;
  struct sGraphAuxiliaryData *GAD = NULL;
  struct sArrayBucketsAtoms *ABA = NULL;
  struct sArrayAtoms *H = NULL;
  char filename[512];
  int *S_bkg, *S_obj;
  int a,s,t,e,k,nv,ne,energy,i;
  float polarity = 0.5;
  struct timeval tic,toc;
  FILE *report;
  
  // check number of parameters
  if(argc < 3){
    fprintf(stdout,"usage:\n");
    fprintf(stdout,"OptMaxNorm <image> <file_seeds> [polarity] [gtruth]\n");
    fprintf(stdout,"image .... input image (.pgm, .png, .jpg),\n");
    fprintf(stdout,"file_seeds ........ the file with seeds (.txt),\n");
    fprintf(stdout,"polarity .......... boundary polarity in the interval [-1, 1] (default = 0.5),\n");
    fprintf(stdout,"gtruth ............ ground truth image.\n");
    exit(0);
  }

  report = fopen("report.csv", "a");
  if(report == NULL){
    printf("Error creating report file.\n");
    exit(0);
  }
  if( ftell(report) == 0 ){ //start of file
    fprintf(report, "OMN time; OMN Acc; OIFT-1 time; OIFT-1 Acc; OIFT-2 time; OIFT-2 Acc; ORFC time; ORFC Acc; ORFC-ch time; ORFC-ch Acc\n");
  }
  
  strcpy(filename, argv[1]);
  img = ReadAnyImage(filename);
  cimg = ReadAnyCImage(filename);

  int niter = ROUND(img->ncols/600);
  printf("niter: %d\n", niter);
  for(i = 1; i <= niter; i++){
    ctmp = glip::CImage::GaussianBlur(cimg);
    glip::CImage::Destroy(&cimg);
    cimg = ctmp;
    tmp = glip::Image32::GaussianBlur(img);
    glip::Image32::Destroy(&img);
    img = tmp;
  }
  
  S_bkg = LoadSeedsTxt(argv[2], 0, img->ncols, img->nrows);
  S_obj = LoadSeedsTxt(argv[2], 1, img->ncols, img->nrows);

  //graph = glip::ImageGraph::ByHomogeneityAffinity(img, 1.5);
  if(argc >= 4)
    polarity = atof(argv[3]);
  printf("polarity: %f\n", polarity);

  if(argc >= 5)
    gtruth = ReadAnyImage(argv[4]);

  if( glip::CImage::IsGrayscale(cimg) ){
    graph = glip::ImageGraph::ByEuclideanDistance(img, 1.5, niter);
  }
  else{
    glip::sCImage32f *cimg_lab;
    cimg_lab = glip::CImage32f::RGB2Lab(cimg);
    graph = glip::ImageGraph::ByEuclideanDistance(cimg_lab, 1.5, niter);
    glip::CImage32f::Destroy(&cimg_lab);
  }

  glip::ImageGraph::Orient2Digraph(graph, img, 100.0*polarity);
  glip::ImageGraph::ChangeType(graph, CAPACITY);

  //------------------------------------------------
  gettimeofday(&tic,NULL);
  //------------------------------------------------
  GAD = CreateGraphAuxiliaryData(graph);

  H = CreateOrderedArrayAtoms(graph, GAD, S_bkg, S_obj);

  ABA = CreateArrayBucketsAtoms(H, GAD->NE.nnodes, GAD->AE.nedges);

  K_vertex = glip::Queue::Create(ABA->nnodes);
  //K_edge   = glip::Queue::Create(ABA->nedges);

  IsInK = (char *)calloc(ABA->nnodes, sizeof(char));
  if(IsInK == NULL){
    printf("Error: insufficient memory\n");
    exit(1);
  }

#ifdef APPDEBUG
  size_t mem_size;
  //mem_size  = GetSizeOf(GAD);
  mem_size  = GetSizeOf(H);
  mem_size += GetSizeOf(ABA);
  mem_size += glip::Queue::GetSizeOf(K_vertex);
  mem_size += ABA->nnodes*sizeof(char); //IsInK
  printf("Memory consumption of MNORA: %zu bytes.\n", mem_size);
#endif
  
  for(a = 0; a < H->natoms; a++){
    if(H->A[a].type == REMOVED) //skip removed atom.
      continue;
    else if(H->A[a].type == UNARY){
      s = H->A[a].vertex_edge;
      nv = ABA->n_vertex[s];
      if(nv > 1){
	k = 0;
	while(ABA->A_vertex[s][k] != a)
	  k++;

	ABA->A_vertex[s][k] = ABA->A_vertex[s][nv-1];
	ABA->n_vertex[s] -= 1;
	IsInK[s] = 1;
	glip::Queue::Push(K_vertex, s);
      }
    }
    else{ //BINARY
      e = H->A[a].vertex_edge;
      ne = ABA->n_edge[e];
      if(ne > 1){
	k = 0;
	while(ABA->A_edge[e][k] != a)
	  k++;

	ABA->A_edge[e][k] = ABA->A_edge[e][ne-1];
	ABA->n_edge[e] -= 1;

	s = GAD->AE.E[e].s;
	t = GAD->AE.E[e].t;
	if(remove_incompatible_unary_atom(s, e, GAD, ABA, H)){
	  IsInK[s] = 1;
	  glip::Queue::Push(K_vertex, s);
	}
	if(remove_incompatible_unary_atom(t, e, GAD, ABA, H)){
	  IsInK[t] = 1;
	  glip::Queue::Push(K_vertex, t);
	}
      }
    }
    while(!glip::Queue::IsEmpty(K_vertex)){
      s = glip::Queue::Pop(K_vertex);
      IsInK[s] = 0;
      for(k = 0; k < GAD->NE.nnedges[s]; k++){
	e = GAD->NE.N[s][k];
	//Remove inconsistentes de ABA->A_edge[e]
	//Se aconteceu remocao:
	if(remove_incompatible_binary_atom(e, s, GAD, ABA, H)){
	  if(GAD->AE.E[e].s != s)
	    t = GAD->AE.E[e].s;
	  else
	    t = GAD->AE.E[e].t;

	  //somente insere t em K se tiver ocorrido alguma
	  //remocao de inconsistência e se não for repetido.
	  if(remove_incompatible_unary_atom(t, e, GAD, ABA, H)){
	    if(IsInK[t] == 0){
	      IsInK[t] = 1;
	      glip::Queue::Push(K_vertex, t);
	    }
	  }
	}
	
      }
    }
  }

  //Para construir segmentação final com base nos átomos que sobraram.
  label = glip::Image32::Create(img);
  for(s = 0; s < label->n; s++){
#ifdef APPDEBUG
    if(ABA->n_vertex[s] != 1)
      printf("Error: wrong amount of atoms remaining.\n");
#endif
    a = ABA->A_vertex[s][0];
    label->data[s] = H->A[a].ls;
  }
  //------------------------------------------------  
  gettimeofday(&toc,NULL);
  PrintTimeAcc(tic, toc, label, gtruth, (char *)"OMN segmentation:",
	       report, ACC_DICE, false);
  //------------------------------------------------
  glip::Image32::Write(label, (char*)"label_omn.pgm");
#ifdef APPDEBUG
  energy = glip::ift::GetEnergy_Max(graph, label, 1);
  printf("\t energy: %d\n", energy);
#endif
  
  free(IsInK);
  glip::Queue::Destroy(&K_vertex);
  //glip::Queue::Destroy(&K_edge);
  DestroyArrayBucketsAtoms(&ABA);
  DestroyArrayAtoms(&H);
  DestroyGraphAuxiliaryData(&GAD);

  //------------------------------------------------
  gettimeofday(&tic,NULL);
  //------------------------------------------------
  OIFT_MinMax_1st(graph, S_bkg, S_obj, label);
  //------------------------------------------------  
  gettimeofday(&toc,NULL);
  PrintTimeAcc(tic, toc, label, gtruth, (char *)"OIFT-1 segmentation:",
	       report, ACC_DICE, false);
  //------------------------------------------------
  glip::Image32::Write(label, (char*)"label_oift-1.pgm");
#ifdef APPDEBUG
  energy = glip::ift::GetEnergy_Max(graph, label, 1);
  printf("\t energy: %d\n", energy);
#endif

#ifdef APPDEBUG
  printf("Memory consumption of OIFT-1: %zu bytes.\n", mem_size_OIFT_1);
#endif
  
  //------------------------------------------------
  gettimeofday(&tic,NULL);
  //------------------------------------------------
  OIFT_MinMax_2nd(graph, S_bkg, S_obj, label);
  //------------------------------------------------  
  gettimeofday(&toc,NULL);
  PrintTimeAcc(tic, toc, label, gtruth, (char *)"OIFT-2 segmentation:",
	       report, ACC_DICE, false);
  //------------------------------------------------
  glip::Image32::Write(label, (char*)"label_oift-2.pgm");
#ifdef APPDEBUG
  energy = glip::ift::GetEnergy_Max(graph, label, 1);
  printf("\t energy: %d\n", energy);
#endif

#ifdef APPDEBUG
  printf("Memory consumption of OIFT-2: %zu bytes.\n", mem_size_OIFT_2);
#endif
  
  //------------------------------------------------
  gettimeofday(&tic,NULL);
  //------------------------------------------------
  ORFC_MinMax(graph, S_bkg, S_obj, label);
  //------------------------------------------------  
  gettimeofday(&toc,NULL);
  PrintTimeAcc(tic, toc, label, gtruth, (char *)"ORFC segmentation:",
	       report, ACC_DICE, false);
  //------------------------------------------------
  glip::Image32::Write(label, (char*)"label_orfc.pgm");
#ifdef APPDEBUG
  energy = glip::ift::GetEnergy_Max(graph, label, 1);
  printf("\t energy: %d\n", energy);
#endif

#ifdef APPDEBUG
  printf("Memory consumption of ORFC: %zu bytes.\n", mem_size_ORFC);
#endif
  
  //------------------------------------------------
  gettimeofday(&tic,NULL);
  //------------------------------------------------
  ORFC_MinMax(graph, S_bkg, S_obj, label);
  tmp = glip::Image32::CloseHoles(label, 1.5);
  glip::Image32::Destroy(&label);
  label = tmp;
  //------------------------------------------------  
  gettimeofday(&toc,NULL);
  PrintTimeAcc(tic, toc, label, gtruth, (char *)"ORFC-ch segmentation:",
	       report, ACC_DICE, true);
  //------------------------------------------------
  glip::Image32::Write(label, (char*)"label_orfc-ch.pgm");
#ifdef APPDEBUG
  energy = glip::ift::GetEnergy_Max(graph, label, 1);
  printf("\t energy: %d\n", energy);
#endif
  
  free(S_bkg);
  free(S_obj);
  glip::ImageGraph::Destroy(&graph);
  glip::Image32::Destroy(&label);
  glip::Image32::Destroy(&img);
  glip::CImage::Destroy(&cimg);  
  if(gtruth != NULL)
    glip::Image32::Destroy(&gtruth);
  fclose(report);
  
  return 0;
}

