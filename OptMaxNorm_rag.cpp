
//#define APPDEBUG  1

#include "glip.h"

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



struct sGraphAuxiliaryData *CreateGraphAuxiliaryData(glip::sGraph *graph){
  struct sGraphAuxiliaryData *GAD;
  int s,t,i;
  int nne, nedges, e;
  
  GAD = (struct sGraphAuxiliaryData *)calloc(1, sizeof(struct sGraphAuxiliaryData));
  if(GAD == NULL){
    printf("Error: CreateGraphAuxiliaryData\n");
    exit(1);
  }

  //allocating memory for graph edge array.
  nedges = glip::Graph::GetNumberOfArcs(graph)/2;
  GAD->AE.nedges = nedges;
  GAD->AE.E = (struct sEdge *)calloc(nedges, sizeof(struct sEdge));
  if(GAD->AE.E == NULL){
    printf("Error: CreateGraphAuxiliaryData\n");
    exit(1);
  }
  
  //allocating memory to store the indices of all neighboring edges of each pixel.
  GAD->NE.nnodes = graph->nnodes;
  GAD->NE.nnedges = (int *)calloc(graph->nnodes, sizeof(int));
  GAD->NE.N = (int **)calloc(graph->nnodes, sizeof(int *));
  if(GAD->NE.nnedges == NULL || GAD->NE.N == NULL){
    printf("Error: CreateGraphAuxiliaryData\n");
    exit(1);
  }

  for(s = 0; s < graph->nnodes; s++){
    GAD->NE.N[s] = (int *)calloc(graph->nodes[s].outdegree, sizeof(int));
    if(GAD->NE.N[s] == NULL){
      printf("Error: CreateGraphAuxiliaryData\n");
      exit(1);
    }
  }
  
  //filling the data structures.
  e = 0;
  for(s = 0; s < graph->nnodes; s++){
    for(i = 0; i < graph->nodes[s].outdegree; i++){
      t = graph->nodes[s].adjList[i];
      if(s < t){
	GAD->AE.E[e].s = s;
	GAD->AE.E[e].t = t;
	GAD->AE.E[e].index = NIL; //not used

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



struct sArrayAtoms *CreateOrderedArrayAtoms(glip::sGraph *graph,
					    glip::sImage32 *spixels,
					    char *filename,
					    struct sGraphAuxiliaryData *GAD){
  struct sArrayAtoms *AA = NULL, *AA_sorted = NULL;
  glip::sPQueue32 *Q=NULL;
  int a,px,py,p,qx,qy,q,s,t,e,i,j,k,n,nedges,phi_max,tmp,nne;
  int phi_0, phi_1, phi_00, phi_01, phi_10, phi_11, max_phi_01_10;
  FILE *fp;

  fp = fopen(filename, "r");
  if(fp == NULL){
    printf("Error opening file.\n");
    exit(1);
  }
  
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
  phi_max = 0;
  a = 0;
  for(s = 0; s < GAD->NE.nnodes; s++){
    AA->Phi[a] = 0;
    AA->A[a].type = UNARY;
    AA->A[a].vertex_edge = s;
    AA->A[a].ls = 0;
    AA->A[a].lt = NIL;
    a++;

    AA->Phi[a] = 0;
    AA->A[a].type = UNARY;
    AA->A[a].vertex_edge = s;
    AA->A[a].ls = 1;
    AA->A[a].lt = NIL;
    a++;
  }

  fscanf(fp,"%d", &n);
  for(i = 0; i < n; i++){
    fscanf(fp,"%d %d %d %d", &px, &py, &phi_0, &phi_1);
    p = px + py*spixels->ncols;
    s = spixels->data[p];
    if(s < 0 || s >= graph->nnodes){
      printf("Error: invalid node in graph_construction_MaxNorm\n");
      exit(1);
    }
    phi_max = MAX(phi_max, phi_0);
    phi_max = MAX(phi_max, phi_1);
    AA->Phi[s*2]   = phi_0;
    AA->Phi[s*2+1] = phi_1;
  }

  //filling the array with atoms derived from edges.
  e = 0;
  for(s = 0; s < graph->nnodes; s++){
    for(i = 0; i < graph->nodes[s].outdegree; i++){
      t = graph->nodes[s].adjList[i];
      if(s < t){

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

	phi_10 = graph->nodes[s].Warcs[i];
	phi_max = MAX(phi_max, phi_10);
	AA->Phi[a] = phi_10;
	AA->A[a].type = BINARY;
	AA->A[a].vertex_edge = e;
	AA->A[a].ls = 1;
	AA->A[a].lt = 0;
	a++;
	
	phi_01 = glip::Graph::GetArcWeight(graph, t, s);
	phi_max = MAX(phi_max, phi_01);
	AA->Phi[a] = phi_01;
	AA->A[a].type = BINARY;
	AA->A[a].vertex_edge = e;
	AA->A[a].ls = 0;
	AA->A[a].lt = 1;
	a++;
	
	e++;
      }
    }
  }

  fscanf(fp,"%d", &n);
  for(i = 0; i < n; i++){
    fscanf(fp,"%d %d %d %d %d %d %d %d", &px, &py, &qx, &qy, &phi_00, &phi_01, &phi_10, &phi_11);
    p = px + py*spixels->ncols;
    q = qx + qy*spixels->ncols;
    s = spixels->data[p];
    t = spixels->data[q];
    if(s == t)
      continue;
    if(s > t){
      tmp = s;      s = t;           t = tmp;
      tmp = phi_01; phi_01 = phi_10; phi_10 = tmp;
    }

    k = 0;
    e = GAD->NE.N[s][0];
    while(GAD->AE.E[e].t != t){
      k++;
      if(k >= GAD->NE.nnedges[s])
	break;
      e = GAD->NE.N[s][k];
    }
    if(k >= GAD->NE.nnedges[s]){
      printf("Warning: penalties could not be added to the graph.\n");
      continue;
    }
    a = 2*GAD->NE.nnodes + e*4;

    if(phi_10 >= 0){
      phi_max = MAX(phi_max, phi_10);
      AA->Phi[a+2] = phi_10;
    }
    
    if(phi_01 >= 0){
      phi_max = MAX(phi_max, phi_01);
      AA->Phi[a+3] = phi_01;
    }
    
    max_phi_01_10 = MAX(AA->Phi[a+2], AA->Phi[a+3]);

    printf("max_phi_01_10: %d\n", max_phi_01_10);
   
    if(phi_00 >= 0){
      if(phi_00 > max_phi_01_10){
	printf("Warning: violation of condition MAX(A,D) <= MAX(C,B) detected.\nA = %d has been reduced to MAX(C,B) = %d.\n",phi_00, max_phi_01_10);
	phi_00 = max_phi_01_10;
      }
      phi_max = MAX(phi_max, phi_00);
      AA->Phi[a] = phi_00;
    }
    
    if(phi_11 >= 0){
      if(phi_11 > max_phi_01_10){
	printf("Warning: violation of condition MAX(A,D) <= MAX(C,B) detected.\nD = %d has been reduced to MAX(C,B) = %d.\n",phi_11, max_phi_01_10);
	phi_11 = max_phi_01_10;
      }
      phi_max = MAX(phi_max, phi_11);
      AA->Phi[a+1] = phi_11;
    }
    
  }
  
  //The ordering of atoms in decreasing order of Phi.
  Q = glip::PQueue32::Create(phi_max+2, AA->natoms, AA->Phi);
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
  free(AA->Phi);
  free(AA->A);
  free(AA);
  fclose(fp);
  
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


//=============OIFT=========================
glip::sGraph *graph_construction_MaxNorm(glip::sGraph *graph,
					 glip::sImage32 *spixels,
					 char *filename){
  glip::sGraph *G;
  int i,e,px,py,p,qx,qy,q,s,t,w,n,node_s,node_t;
  int phi_0, phi_1, phi_00, phi_01, phi_10, phi_11, max_phi_01_10;
  FILE *fp;

  fp = fopen(filename, "r");
  if(fp == NULL){
    printf("Error opening file.\n");
    exit(1);
  }

  G = glip::Graph::Create(graph->nnodes+2, 10, NULL);
  for(p = 0; p < graph->nnodes; p++){
    for(i = 0; i < graph->nodes[p].outdegree; i++){
      q = graph->nodes[p].adjList[i];
      w = graph->nodes[p].Warcs[i];
      glip::Graph::AddDirectedEdge(G, p, q, w);
    }
  }

  node_s = graph->nnodes;
  node_t = graph->nnodes+1;
  fscanf(fp,"%d", &n); //number of unary relations.
  for(i = 0; i < n; i++){
    fscanf(fp,"%d %d %d %d", &px, &py, &phi_0, &phi_1);
    p = px + py*spixels->ncols;
    s = spixels->data[p];
    if(s < 0 || s >= graph->nnodes){
      printf("Error: invalid node in graph_construction_MaxNorm\n");
      exit(1);
    }
    glip::Graph::AddDirectedEdge(G, node_s, s, phi_0);
    glip::Graph::AddDirectedEdge(G, s, node_s,     0);
    
    glip::Graph::AddDirectedEdge(G, s, node_t, phi_1);
    glip::Graph::AddDirectedEdge(G, node_t, s,     0);
  }

  fscanf(fp,"%d", &n);
  for(e = 0; e < n; e++){
    fscanf(fp,"%d %d %d %d %d %d %d %d", &px, &py, &qx, &qy, &phi_00, &phi_01, &phi_10, &phi_11); //number of binary relations with equal label.
    p = px + py*spixels->ncols;
    q = qx + qy*spixels->ncols;
    s = spixels->data[p];
    t = spixels->data[q];
    if(s == t)
      continue;

    if(phi_10 >= 0)
      glip::Graph::UpdateDirectedEdge(G, s, t, phi_10);
    if(phi_01 >= 0)
      glip::Graph::UpdateDirectedEdge(G, t, s, phi_01);

    max_phi_01_10 = MAX(glip::Graph::GetArcWeight(G, s, t),
			glip::Graph::GetArcWeight(G, t, s));
    
    if(phi_00 > max_phi_01_10){
      printf("Warning: violation of condition MAX(A,D) <= MAX(C,B) detected.\nA = %d has been reduced to MAX(C,B) = %d.\n",phi_00, max_phi_01_10);
      phi_00 = max_phi_01_10;
    }
    if(phi_11 > max_phi_01_10){
      printf("Warning: violation of condition MAX(A,D) <= MAX(C,B) detected.\nD = %d has been reduced to MAX(C,B) = %d.\n",phi_11, max_phi_01_10);
      phi_11 = max_phi_01_10;
    }

    if(phi_01 >= phi_10){
      if(phi_00 >= 0){
	glip::Graph::UpdateDirectedEdgeIfHigher(G, node_s, s, phi_00);
	glip::Graph::UpdateDirectedEdgeIfHigher(G, s, node_s,      0);
      }
      if(phi_11 >= 0){
	glip::Graph::UpdateDirectedEdgeIfHigher(G, t, node_t, phi_11);
	glip::Graph::UpdateDirectedEdgeIfHigher(G, node_t, t,      0);
      }
    }
    else{
      if(phi_00 >= 0){
	glip::Graph::UpdateDirectedEdgeIfHigher(G, node_s, t, phi_00);
	glip::Graph::UpdateDirectedEdgeIfHigher(G, t, node_s,      0);
      }
      if(phi_11 >= 0){
	glip::Graph::UpdateDirectedEdgeIfHigher(G, s, node_t, phi_11);
	glip::Graph::UpdateDirectedEdgeIfHigher(G, node_t, s,      0);
      }
    }
  }
  fclose(fp);
  return G;
}


int *Segmentation_OIFT_MinMax_1st(glip::sGraph *graph,
				  int *S_bkg,
				  int *S_obj){
  glip::sGraph *transpose;
  glip::sPQueue32 *Q=NULL;
  glip::sGraph *g;
  int i,j,p,q,n;
  int w,d;
  int *label;
  int *value;
  int Wmax;
  Wmax = glip::Graph::GetMaximumArc(graph);
  n = graph->nnodes;
  transpose = glip::Graph::Transpose(graph);
  label = (int *)calloc(n, sizeof(int));
  value = (int *)calloc(n, sizeof(int));
  Q = glip::PQueue32::Create(Wmax*2+3, n, value);

  for(i=1; i<=S_bkg[0]; i++){
    p = S_bkg[i];
    value[p] = Wmax*2+1+1;
    label[p] = 0;
    glip::PQueue32::FastInsertElem(Q, p);
  }
  for(i=1; i<=S_obj[0]; i++){
    p = S_obj[i];
    value[p] = Wmax*2+1+1;
    label[p] = 1;
    glip::PQueue32::FastInsertElem(Q, p);
  }

  while(!glip::PQueue32::IsEmpty(Q)) {
    p = glip::PQueue32::FastRemoveMaxFIFO(Q);
    
    if(label[p]==0){ g = transpose; d = 1; }
    else{   	     g = graph;     d = 0; }
    
    for(i = 0; i < g->nodes[p].outdegree; i++){
      q = g->nodes[p].adjList[i];
      
      if(Q->L.elem[q].color != BLACK){
	w = g->nodes[p].Warcs[i] * 2 + d;
	w = MIN(w, value[p]);
	if(w > value[q]){
	  if(Q->L.elem[q].color == GRAY)
	    glip::PQueue32::FastRemoveElem(Q, q);
	  value[q] = w;
	  label[q] = label[p];
	  glip::PQueue32::FastInsertElem(Q, q);
	}
      }
    }
  }

  free(value);
  glip::PQueue32::Destroy(&Q);
  glip::Graph::Destroy(&transpose);
  return label;
}



int *Segmentation_OIFT_MinMax_2nd(glip::sGraph *graph,
				  int *S_bkg,
				  int *S_obj){
  glip::sGraph *transpose;
  glip::sPQueue32 *Q=NULL;
  glip::sGraph *g;
  int i,j,p,q,n;
  int w;
  int *label;
  int *value;
  int Wmax;
  Wmax = glip::Graph::GetMaximumArc(graph);
  n = graph->nnodes;
  transpose = glip::Graph::Transpose(graph);
  label = (int *)calloc(n, sizeof(int));
  value = (int *)calloc(n, sizeof(int));
  Q = glip::PQueue32::Create(Wmax+2, n, value);

  for(i=1; i<=S_bkg[0]; i++){
    p = S_bkg[i];
    value[p] = Wmax+1;
    label[p] = 0;
    glip::PQueue32::FastInsertElem(Q, p);
  }
  for(i=1; i<=S_obj[0]; i++){
    p = S_obj[i];
    value[p] = Wmax+1;
    label[p] = 1;
    glip::PQueue32::FastInsertElem(Q, p);
  }

  while(!glip::PQueue32::IsEmpty(Q)) {
    p = glip::PQueue32::FastRemoveMaxFIFO(Q);
    
    if(label[p]==0){ g = transpose; }
    else{   	     g = graph;     }
    
    for(i = 0; i < g->nodes[p].outdegree; i++){
      q = g->nodes[p].adjList[i];
      
      if(Q->L.elem[q].color != BLACK){
	w = g->nodes[p].Warcs[i];
	if(w > value[q]){
	  if(Q->L.elem[q].color == GRAY)
	    glip::PQueue32::FastRemoveElem(Q, q);
	  value[q] = w;
	  label[q] = label[p];
	  glip::PQueue32::FastInsertElem(Q, q);
	}
      }
    }
  }

  free(value);
  glip::PQueue32::Destroy(&Q);
  glip::Graph::Destroy(&transpose);
  return label;
}



int *Segmentation_OIFT_MinMax_1st(glip::sGraph *graph,
				  glip::sImage32 *spixels,
				  char *filename){
  glip::sGraph *G;
  int S_bkg[2];
  int S_obj[2];
  int *label;
  int n;
  G = graph_construction_MaxNorm(graph, spixels, filename);
  n = G->nnodes;
  S_bkg[0] = 1; S_bkg[1] = n-1;
  S_obj[0] = 1; S_obj[1] = n-2;
  label = Segmentation_OIFT_MinMax_1st(G, S_bkg, S_obj);
#ifdef APPDEBUG
  int energy;
  energy = glip::ift::GetEnergy_Max(G, label, 1);
  printf("energy: %d\n", energy);
#endif
  glip::Graph::Destroy(&G);
  return label;
}


int *Segmentation_OIFT_MinMax_2nd(glip::sGraph *graph,
				  glip::sImage32 *spixels,
				  char *filename){
  glip::sGraph *G;
  int S_bkg[2];
  int S_obj[2];
  int *label;
  int n;
  G = graph_construction_MaxNorm(graph, spixels, filename);
  n = G->nnodes;
  S_bkg[0] = 1; S_bkg[1] = n-1;
  S_obj[0] = 1; S_obj[1] = n-2;
  label = Segmentation_OIFT_MinMax_2nd(G, S_bkg, S_obj);
#ifdef APPDEBUG
  int energy;
  energy = glip::ift::GetEnergy_Max(G, label, 1);
  printf("energy: %d\n", energy);
#endif
  glip::Graph::Destroy(&G);
  return label;
}

//=============ORFC=========================
void IFT_fmin(glip::sGraph *graph,
	      int *S,
	      int *cost){
  glip::sPQueue32 *Q;
  int tmp, w, Wmax;
  int n,p,q,i;
  n = graph->nnodes;
  Wmax = glip::Graph::GetMaximumArc(graph);
  Q = glip::PQueue32::Create(Wmax+2, n, cost);
  
  for(p = 0; p < n; p++)
    cost[p] = INT_MIN;
  
  for(i = 1; i <= S[0]; i++){
    cost[S[i]] = Wmax+1;
    glip::PQueue32::FastInsertElem(Q, S[i]);
  }
  
  while(!glip::PQueue32::IsEmpty(Q)){
    p = glip::PQueue32::FastRemoveMaxFIFO(Q);
    
    for(i = 0; i < graph->nodes[p].outdegree; i++){
      q = graph->nodes[p].adjList[i];
      if(Q->L.elem[q].color != BLACK){
	w = graph->nodes[p].Warcs[i];
	tmp = MIN(cost[p], w);
	
	if(tmp > cost[q]){
	  if(Q->L.elem[q].color == GRAY)
	    glip::PQueue32::FastRemoveElem(Q, q);
	  cost[q] = tmp;
	  glip::PQueue32::FastInsertElem(Q, q);
	}
      }
    }
  }
  glip::PQueue32::Destroy(&Q);
}


int *Segmentation_ORFC_MinMax(glip::sGraph *graph,
			      int *S_bkg,
			      int *S_obj){
  glip::sGraph *transpose;
  int *L, *E;
  glip::sPQueue32 *QS=NULL;
  glip::sQueue *Q = glip::Queue::Create(graph->nnodes);
  int p,q,j,i,k,energy,nsi,Wmax;
  int *s_energy,*s_pixel;

  transpose = glip::Graph::Transpose(graph);
  L = (int *)calloc(graph->nnodes, sizeof(int));
  E = glip::AllocIntArray(graph->nnodes);
  IFT_fmin(transpose, S_bkg, E);

  nsi = S_obj[0];
  s_energy = glip::AllocIntArray(nsi);
  s_pixel  = glip::AllocIntArray(nsi);
  Wmax = glip::Graph::GetMaximumArc(graph);
  QS = glip::PQueue32::Create(Wmax+2, nsi, s_energy);
  
  k = 0;
  for(i = 1; i <= S_obj[0]; i++){
    p = S_obj[i];
    s_energy[k] = E[p];
    s_pixel[k] = p;
    glip::PQueue32::InsertElem(&QS, k);
    k++;
    L[p] = 0;
  }
  
  while(!glip::PQueue32::IsEmpty(QS)){
    j = glip::PQueue32::RemoveMinFIFO(QS);
    p = s_pixel[j];
    energy = s_energy[j];
    
    if(L[p] != 0) continue;
    
    glip::Queue::Reset(Q);
    glip::Queue::Push(Q, p);
    L[p] = 1;
    
    while(!glip::Queue::IsEmpty(Q)){
      p = glip::Queue::Pop(Q);
      for(i = 0; i < graph->nodes[p].outdegree; i++){
	q = graph->nodes[p].adjList[i];
	if(energy < graph->nodes[p].Warcs[i] && L[q] == 0){
	  L[q] = 1;
	  glip::Queue::Push(Q, q);
	}
      }
    }
  }
  glip::Queue::Destroy(&Q);
  glip::PQueue32::Destroy(&QS);
  glip::FreeIntArray(&s_energy);
  glip::FreeIntArray(&s_pixel);
  glip::FreeIntArray(&E);
  glip::Graph::Destroy(&transpose);
  return L;
}


int *Segmentation_ORFC_MinMax(glip::sGraph *graph,
			      glip::sImage32 *spixels,
			      char *filename){
  glip::sGraph *G;
  int S_bkg[2];
  int S_obj[2];
  int *label;
  int n;
  G = graph_construction_MaxNorm(graph, spixels, filename);
  n = G->nnodes;
  S_bkg[0] = 1; S_bkg[1] = n-1;
  S_obj[0] = 1; S_obj[1] = n-2;
  label = Segmentation_ORFC_MinMax(G, S_bkg, S_obj);
#ifdef APPDEBUG
  int energy;
  energy = glip::ift::GetEnergy_Max(G, label, 1);
  printf("energy: %d\n", energy);
#endif
  glip::Graph::Destroy(&G);
  return label;
}



void PrintTime(struct timeval tic,
	       struct timeval toc,
	       char *message){
  double totaltime;
  totaltime = ((toc.tv_sec-tic.tv_sec)*1000.0 + 
	       (toc.tv_usec-tic.tv_usec)*0.001);
  printf("%s\n\t time: %f ms\n", message, totaltime);
}



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



glip::sImage32 *SPixelsLabels2Image(int *L,
				   glip::sImage32 *spixels){
  glip::sImage32 *label;
  int p;
  label = glip::Image32::Create(spixels);
  for(p = 0; p < spixels->n; p++)
    label->data[p] = L[spixels->data[p]];
  return label;
}


int main(int argc, char **argv){
  glip::sQueue *K_vertex = NULL;
  char *IsInK = NULL;
  struct sGraphAuxiliaryData *GAD = NULL;
  struct sArrayBucketsAtoms *ABA = NULL;
  struct sArrayAtoms *H = NULL;
  char filename[512];
  int a,s,t,e,k,nv,ne,energy,lb;
  struct timeval tic,toc;
  glip::sImage32 *img = NULL, *spixels = NULL, *gtruth = NULL;
  int k_sp, SPsize = 100;
  float polarity = 0.5, alpha = 0.06;
  glip::sGraph *graph;
  double dice = 0.0;
  glip::sImage32 *label = NULL;
  int *L = NULL;
  
  // check number of parameters
  if(argc < 3){
    fprintf(stdout,"usage:\n");
    fprintf(stdout,"OptMaxNorm_rag <image> <penalties> [polarity] [SPsize] [gtruth]\n");
    fprintf(stdout,"image ....... image (.pgm),\n");
    fprintf(stdout,"penalties ... text file with phi penalties (.txt),\n");
    fprintf(stdout,"polarity .......... boundary polarity in the interval [-1, 1] (default = 0.5),\n");
    fprintf(stdout,"SPsize ............ average size of a superpixel in number of pixels (default = 100),\n");
    fprintf(stdout,"gtruth ............ ground truth image.\n");
    exit(0);
  }

  img = ReadAnyImage(argv[1]);
  if(argc >= 5)
    SPsize = atoi(argv[4]);
  printf("SPsize: %d\n", SPsize);
  k_sp = ROUND(img->n/((float)SPsize));
  //------------------------------------------------
  gettimeofday(&tic,NULL);
  //------------------------------------------------
  spixels = glip::Superpixels::IFT_SLIC(img, k_sp, alpha,
				       0.5, 0.5, 10);
  //------------------------------------------------  
  gettimeofday(&toc,NULL);
  PrintTime(tic, toc, (char *)"Superpixel time:");
  //------------------------------------------------ 
#ifdef APPDEBUG
  glip::sCImage *ctmp;
  glip::Image32::Write(spixels, (char *)"spixels.pgm");
  ctmp = glip::CImage::RandomColorize(spixels);
  glip::CImage::Write(ctmp, (char *)"spixels.ppm");
  glip::CImage::Destroy(&ctmp);
  printf("spixel: min=%d, max=%d\n",
	 glip::Image32::GetMinVal(spixels),
	 glip::Image32::GetMaxVal(spixels));
#endif
  strcpy(filename, argv[2]);

  if(argc >= 4)
    polarity = atof(argv[3]);
  printf("polarity: %f\n", polarity);
  graph = glip::Graph::RAG(img, spixels, polarity, 1.0);
  glip::Graph::ChangeType(graph, CAPACITY);
  
  if(argc >= 6)
    gtruth = ReadAnyImage(argv[5]);
  
  
  //------------------------------------------------
  printf("----------------------------\n");
  gettimeofday(&tic,NULL);
  //------------------------------------------------
  GAD = CreateGraphAuxiliaryData(graph);

  H = CreateOrderedArrayAtoms(graph, spixels, filename, GAD);

  ABA = CreateArrayBucketsAtoms(H, GAD->NE.nnodes, GAD->AE.nedges);

  K_vertex = glip::Queue::Create(ABA->nnodes);

  IsInK = (char *)calloc(ABA->nnodes, sizeof(char));
  if(IsInK == NULL){
    printf("Error: insufficient memory\n");
    exit(1);
  }
  
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

  //------------------------------------------------  
  gettimeofday(&toc,NULL);
  PrintTime(tic, toc, (char *)"MNORA segmentation:");
  //------------------------------------------------

  L = (int *)calloc(GAD->NE.nnodes, sizeof(int));
  if(L == NULL){
    printf("Error: insufficient memory\n");
    exit(1);
  }
  
  //Para montar a segmentação final com base nos átomos que sobraram.
  for(s = 0; s < GAD->NE.nnodes; s++){
#ifdef APPDEBUG
    if(ABA->n_vertex[s] != 1)
      printf("Error: wrong amount of atoms remaining.\n");
#endif
    a = ABA->A_vertex[s][0];
    lb = H->A[a].ls;
    L[s] = lb;
  }

#ifdef APPDEBUG
  //computing energy:
  energy = INT_MIN;
  for(s = 0; s < GAD->NE.nnodes; s++){
    a = ABA->A_vertex[s][0];
    lb = H->A[a].ls;
    energy = MAX(energy, H->Phi[a]);
  }
  for(e = 0; e < GAD->AE.nedges; e++){
    if(ABA->n_edge[e] != 1)
      printf("Error: wrong amount of atoms remaining.\n");

    a = ABA->A_edge[e][0];
    energy = MAX(energy, H->Phi[a]);
  }
  printf("energy: %d\n", energy);
#endif

  free(IsInK);
  glip::Queue::Destroy(&K_vertex);
  DestroyArrayBucketsAtoms(&ABA);
  DestroyArrayAtoms(&H);
  DestroyGraphAuxiliaryData(&GAD);

  label = SPixelsLabels2Image(L, spixels);
  glip::Image32::Write(label, (char*)"label_mnora.pgm");
  if(gtruth != NULL){
    dice = glip::Image32::DiceSimilarity(gtruth, label);
    printf("Dice: %lf\n",dice);
  }
  glip::Image32::Destroy(&label);
  free(L);
  
  //------------------------------------------------
  printf("----------------------------\n");
  gettimeofday(&tic,NULL);
  //------------------------------------------------

  L = Segmentation_OIFT_MinMax_1st(graph, spixels, filename);
 
  //------------------------------------------------  
  gettimeofday(&toc,NULL);
  PrintTime(tic, toc, (char *)"OIFT 1st segmentation:");
  //------------------------------------------------

  label = SPixelsLabels2Image(L, spixels);
  glip::Image32::Write(label, (char*)"label_oift-1.pgm");
  if(gtruth != NULL){
    dice = glip::Image32::DiceSimilarity(gtruth, label);
    printf("Dice: %lf\n",dice);
  }
  glip::Image32::Destroy(&label);
  free(L);

  
  //------------------------------------------------
  printf("----------------------------\n");
  gettimeofday(&tic,NULL);
  //------------------------------------------------

  L = Segmentation_OIFT_MinMax_2nd(graph, spixels, filename);
 
  //------------------------------------------------  
  gettimeofday(&toc,NULL);
  PrintTime(tic, toc, (char *)"OIFT 2nd segmentation:");
  //------------------------------------------------

  label = SPixelsLabels2Image(L, spixels);
  glip::Image32::Write(label, (char*)"label_oift-2.pgm");
  if(gtruth != NULL){
    dice = glip::Image32::DiceSimilarity(gtruth, label);
    printf("Dice: %lf\n",dice);
  }
  glip::Image32::Destroy(&label);
  free(L);
  
  //------------------------------------------------
  printf("----------------------------\n");
  gettimeofday(&tic,NULL);
  //------------------------------------------------

  L = Segmentation_ORFC_MinMax(graph, spixels, filename);
 
  //------------------------------------------------  
  gettimeofday(&toc,NULL);
  PrintTime(tic, toc, (char *)"ORFC segmentation:");
  //------------------------------------------------

  label = SPixelsLabels2Image(L, spixels);
  glip::Image32::Write(label, (char*)"label_orfc.pgm");
  if(gtruth != NULL){
    dice = glip::Image32::DiceSimilarity(gtruth, label);
    printf("Dice: %lf\n",dice);
  }
  glip::Image32::Destroy(&label);
  free(L);
  
  return 0;
}

