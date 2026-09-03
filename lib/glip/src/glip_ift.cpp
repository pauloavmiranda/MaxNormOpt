
#include "glip_ift.h"
#include <queue>

namespace glip{
  namespace ift{

    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_conquer_path(int p,
			 sImageGraph *sg,
			 sImage32 *P_sum, 
			 sImage32 *V,
			 sPQueue32 *Q,
			 sImage32 *label);
    
    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_prune_tree(int p,
		       sImageGraph *sg,
		       sImage32 *P_sum, 
		       sImage32 *V,
		       sPQueue32 *Q,
		       sQueue *Qfifo,
		       sImage32 *label);
    
    sSet *COIFT_new_seeds(sImageGraph *sg,
			  sSet *Si1,
			  sSet *Si2,
			  sImage32 *E,
			  sImage32 *label_oift);

    //View Boundary Band
    /*
    Image32::Image32 *BB_OIFT_View(ImageGraph::ImageGraph *G,
				   Image32::Image32 *L,
				   Image32::Image32 *C,
				   float delta,
				   int band_id);
    */ 
   
    void BB_OIFT_Propagate(int p,
			   sImageGraph *sg,
			   sPQueue32 *Q,
			   sImage32 *V,
			   sImage32 *L,
			   int *i_inv);
    
    void BB_OIFT_Propagate_bkg(int p,
			       sImageGraph *sg,
			       sPQueue32 *Q,
			       sPQueue32 *Qe,
			       sImage32 *V,
			       sImage32 *L,
			       int *i_inv);
    
    void BB_OIFT_Propagate_obj(int p,
			       sImageGraph *sg,
			       sPQueue32 *Q,
			       sPQueue32 *Qi,
			       sImage32 *V,
			       sImage32 *L,
			       int *i_inv);

    sImage32 *BB_OIFT_GetLeafNodes(sImage32 *pred,
				   sAdjRel *A);
    
    //----------------------------------

    void IFT_fmax(sImageGraph *g,
		  int *S,
		  sImage32 *label,
		  sImage32 *cost,
		  sImage32 *pred){
      sPQueue32 *Q=NULL;
      int i,p,q,n,cst;
      Pixel u,v;
      sAdjRel *A = g->A;

      n = g->ncols*g->nrows;
      Q = PQueue32::Create(2+g->Wmax, n, cost->data);

      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) cost->data[p] = INT_MAX;
	else                    cost->data[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if (glip::Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){	  
	    q = v.x + v.y*label->ncols;

	    cst = MAX( (g->n_link[p])[i], cost->data[p] );
	    
	    if(cst < cost->data[q]){	    
	      if(Q->L.elem[q].color == GRAY)
		glip::PQueue32::FastRemoveElem(Q, q);
	      cost->data[q]  = cst;
	      pred->data[q]  = p;
	      label->data[q] = label->data[p];
	      glip::PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      glip::PQueue32::Destroy(&Q);
    }

    

    void IFT_feuc(sImageGraph *g,
		  int *S,
		  sImage32 *label,
		  sImage32 *cost,
		  sImage32 *pred){
      sImage32 *Dx=NULL,*Dy=NULL;
      sPQueue32 *Q=NULL;
      int i,p,q,n,cst,dx,dy,sq_diag;
      Pixel u,v;
      sAdjRel *A = g->A;

      Dx = Image32::Create(cost->ncols, cost->nrows);
      Dy = Image32::Create(cost->ncols, cost->nrows);

      n = g->ncols*g->nrows;
      sq_diag = g->ncols*g->ncols + g->nrows*g->nrows;
      Q = PQueue32::Create(sq_diag+1, n, cost->data);

      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) cost->data[p] = INT_MAX;
	else                    cost->data[p] = 0;
      }

      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if (glip::Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){	  
	    q = v.x + v.y*label->ncols;

	    dx  = Dx->data[p] + abs(A->dx[i]);
	    dy  = Dy->data[p] + abs(A->dy[i]);
	    cst = dx*dx + dy*dy; //(g->n_link[p])[i];
	    
	    if(cst < cost->data[q]){	    
	      if(Q->L.elem[q].color == GRAY)
		glip::PQueue32::FastRemoveElem(Q, q);
	      Dx->data[q] = dx;
	      Dy->data[q] = dy;
	      cost->data[q]  = cst;
	      pred->data[q]  = p;
	      label->data[q] = label->data[p];
	      glip::PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      
      glip::PQueue32::Destroy(&Q);
      Image32::Destroy(&Dx);
      Image32::Destroy(&Dy);
    }




    void IFT_fsum(sImageGraph *g,
		  int *S,
		  sImage32 *label,
		  sImage32 *cost,
		  sImage32 *pred){
      sPQueue32 *Q=NULL;
      int i,p,q,n,cst;
      Pixel u,v;
      sAdjRel *A = g->A;

      n = g->ncols*g->nrows;
      Q = PQueue32::Create(g->Wmax+2, n, cost->data);

      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) cost->data[p] = INT_MAX;
	else                    cost->data[p] = 0;
      }

      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::InsertElem(&Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::InsertElem(&Q, p);
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::RemoveMinFIFO(Q);
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if (glip::Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){	  
	    q = v.x + v.y*label->ncols;

	    cst = cost->data[p] + (g->n_link[p])[i];
	    
	    if(cst < cost->data[q]){	    
	      if(Q->L.elem[q].color == GRAY)
		glip::PQueue32::RemoveElem(Q, q);
	      cost->data[q]  = cst;
	      pred->data[q]  = p;
	      label->data[q] = label->data[p];
	      glip::PQueue32::InsertElem(&Q, q);
	    }
	  }
	}
      }
      
      glip::PQueue32::Destroy(&Q);
    }

    
    
    void IFT_fw(sImageGraph *g,
		int *S,
		sImage32 *label,
		sImage32 *cost,
		sImage32 *pred){
      sPQueue32 *Q=NULL;
      int i,p,q,n,cst;
      Pixel u,v;
      sAdjRel *A = g->A;
      
      n = g->ncols*g->nrows;
      Q = PQueue32::Create(g->Wmax+2, n, cost->data);

      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) cost->data[p] = INT_MAX;
	else                    cost->data[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if (glip::Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){	  
	    q = v.x + v.y*label->ncols;
	    
	    cst = (g->n_link[p])[i];
	    
	    if(cst < cost->data[q]){	    
	      if(Q->L.elem[q].color == GRAY)
		glip::PQueue32::FastRemoveElem(Q, q);
	      cost->data[q]  = cst;
	      pred->data[q]  = p;
	      label->data[q] = label->data[p];
	      glip::PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      glip::PQueue32::Destroy(&Q);
    }
    
    //----------------------------------
    
    int GetEnergy_Min(sImageGraph *sg,
		      sImage32 *label,
		      int lb){
      sAdjRel *A;
      int u_x,u_y,v_x,v_y,p,q,n,i;
      int energy,w;

      energy = INT_MAX;
      A = sg->A;
      n = sg->ncols*sg->nrows;
      for(p = 0; p < n; p++){
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(label->data[p] == lb && label->data[q] != lb){
	      w = (sg->n_link[p])[i];
	      energy = MIN(energy, w);
	    }
	  }
	}
      }
      return energy;
    }


    int GetEnergy_Max(sImageGraph *sg,
		      sImage32 *label,
		      int lb){
      sAdjRel *A;
      int u_x,u_y,v_x,v_y,p,q,n,i;
      int energy,w;

      energy = INT_MIN;
      A = sg->A;
      n = sg->ncols*sg->nrows;
      for(p = 0; p < n; p++){
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(label->data[p] == lb && label->data[q] != lb){
	      w = (sg->n_link[p])[i];
	      energy = MAX(energy, w);
	    }
	  }
	}
      }
      return energy;
    }
    

    long long GetEnergy_Sum(sImageGraph *sg,
			    sImage32 *label,
			    int lb){
      sAdjRel *A;
      int u_x,u_y,v_x,v_y,p,q,n,i;
      long long energy,w;

      energy = 0;
      A = sg->A;
      n = sg->ncols*sg->nrows;
      for(p = 0; p < n; p++){
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(label->data[p] == lb && label->data[q] != lb){
	      w = (sg->n_link[p])[i];
	      energy += w;
	    }
	  }
	}
      }
      return energy;
    }



    int GetEnergy_Min(sGraph *graph,
		      int *label,
		      int lb){
      int p,q,i;
      int energy,w;

      energy = INT_MAX;
      for(p = 0; p < graph->nnodes; p++){
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  w = graph->nodes[p].Warcs[i];
	  if(label[p] == lb && label[q] != lb){
	    energy = MIN(energy, w);
	  }
	}
      }
      return energy;
    }


    int GetEnergy_Max(sGraph *graph,
		      int *label,
		      int lb){
      int p,q,i;
      int energy,w;

      energy = INT_MIN;
      for(p = 0; p < graph->nnodes; p++){
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  w = graph->nodes[p].Warcs[i];
	  if(label[p] == lb && label[q] != lb){
	    energy = MAX(energy, w);
	  }
	}
      }
      return energy;
    }
    

    float GetEnergy_Mean(sGraph *graph,
			 int *label,
			 int lb){
      int p,q,i;
      int w, n = 0;
      float sum = 0.0;
      for(p = 0; p < graph->nnodes; p++){
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  w = graph->nodes[p].Warcs[i];
	  if(label[p] == lb && label[q] != lb){
	    sum += w;
	    n++;
	  }
	}
      }
      return (sum/n);
    }

    
    
    
    /*Weighted Distance Transform.*/
    sImage32 *SC_Pred_fsum(sImageGraph *sg,
			   int *S,
			   float power){
      sHeap *Q=NULL;
      int i,p,q,n;
      float edge,tmp;
      float *cost=NULL;
      int u_x,u_y,v_x,v_y;
      sImage32 *pred;
      sAdjRel *A;
      float *Dpq;
      
      n    = sg->ncols*sg->nrows;
      pred = Image32::Create(sg->ncols, sg->nrows);
      cost = glip::AllocFloatArray(n);
      Q = Heap::Create(n, cost);
      A = sg->A;

      //--------------------
      Dpq = (float *)malloc(A->n*sizeof(float));
      for(i=1; i<A->n; i++){
	Dpq[i] = sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]);
      }
      //--------------------
      
      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++)
	cost[p] = FLT_MAX;
      
      for(i=1; i<=S[0]; i++){
	cost[S[i]] = 0.0;
	Heap::Insert_MinPolicy(Q, S[i]);
      }
	
      while(!Heap::IsEmpty(Q)){
	Heap::Remove_MinPolicy(Q, &p);
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(Q->color[q] != BLACK){
	      
	      edge = (sg->n_link[p])[i];
	      tmp  = cost[p] + powf(MAX(edge,1.0), power) - 1.0 + Dpq[i];
	      
	      if(tmp < cost[q]){
		Heap::Update_MinPolicy(Q, q, tmp);
		pred->data[q] = p;
	      }
	    }
	  }
	}
      }
      free(Dpq);
      glip::FreeFloatArray(&cost);
      Heap::Destroy(&Q);
      return pred;
    }


    
    sScene32 *SC_Pred_fsum(sGraph *graph,
			   sScene32 *scn,
			   int *S,
			   float power){
      sHeap *Q=NULL;
      int i,p,q,n;
      float edge,tmp;
      float *cost=NULL;
      int u_x,u_y,u_z,v_x,v_y,v_z;
      glip::sScene32 *pred;
      float Dpq,dx2,dy2,dz2;
      dx2 = scn->dx*scn->dx;
      dy2 = scn->dy*scn->dy;
      dz2 = scn->dz*scn->dz;
      n    = graph->nnodes;
      pred = glip::Scene32::Create(scn);
      cost = glip::AllocFloatArray(n);
      Q = Heap::Create(n, cost);
      
      Scene32::Fill(pred, NIL);
      for(p = 0; p < n; p++)
	cost[p] = FLT_MAX;
      
      for(i=1; i<=S[0]; i++){
	cost[S[i]] = 0.0;
	Heap::Insert_MinPolicy(Q, S[i]);
      }
	
      while(!Heap::IsEmpty(Q)){
	Heap::Remove_MinPolicy(Q, &p);
	u_x = glip::Scene32::GetAddressX(scn, p);
	u_y = glip::Scene32::GetAddressY(scn, p);
	u_z = glip::Scene32::GetAddressZ(scn, p);	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->color[q] != BLACK){
	    v_x = glip::Scene32::GetAddressX(scn, q);
	    v_y = glip::Scene32::GetAddressY(scn, q);
	    v_z = glip::Scene32::GetAddressZ(scn, q);
	    Dpq = sqrtf((u_x-v_x)*(u_x-v_x)*dx2 + 
			(u_y-v_y)*(u_y-v_y)*dy2 + 
			(u_z-v_z)*(u_z-v_z)*dz2);
	    edge = graph->nodes[p].Warcs[i];
	    tmp  = cost[p] + powf(MAX(edge,1.0), power) - 1.0 + Dpq;
	    
	    if(tmp < cost[q]){
	      Heap::Update_MinPolicy(Q, q, tmp);
	      pred->data[q] = p;
	    }
	  }
	}
      }
      glip::FreeFloatArray(&cost);
      Heap::Destroy(&Q);
      return pred;
    }
    


    sScene32 *SC_Pred_fsum(sScene32 *W,
			   sAdjRel3 *A,
			   int *S, float power){
      sHeap *Q=NULL;
      int i,p,q,n;
      float edge,tmp;
      float *cost=NULL;
      glip::Voxel u,v;
      sScene32 *pred;
      float *Dpq;
      
      n    = W->n;
      pred = Scene32::Create(W);
      cost = glip::AllocFloatArray(n);
      Q = Heap::Create(n, cost);

      //--------------------
      Dpq = (float *)malloc(A->n*sizeof(float));
      for(i=1; i<A->n; i++){
	Dpq[i] = sqrtf(A->d[i].axis.x*A->d[i].axis.x*W->dx*W->dx + 
		       A->d[i].axis.y*A->d[i].axis.y*W->dy*W->dy + 
		       A->d[i].axis.z*A->d[i].axis.z*W->dz*W->dz);
      }
      //--------------------
      
      Scene32::Fill(pred, NIL);
      for(p = 0; p < n; p++)
	cost[p] = FLT_MAX;
      
      for(i=1; i<=S[0]; i++){
	cost[S[i]] = 0.0;
	Heap::Insert_MinPolicy(Q, S[i]);
      }
	
      while(!Heap::IsEmpty(Q)){
	Heap::Remove_MinPolicy(Q, &p);
	u.c.x = glip::Scene32::GetAddressX(W, p);
	u.c.y = glip::Scene32::GetAddressY(W, p);
	u.c.z = glip::Scene32::GetAddressZ(W, p);
	for(i=1; i<A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(glip::Scene32::IsValidVoxel(W,v)){
	    q = glip::Scene32::GetVoxelAddress(W, v);
	    if(Q->color[q] != BLACK){
	      edge = W->data[p] + W->data[q];
	      tmp  = cost[p] + powf(MAX(edge,1.0), power) - 1.0 + Dpq[i];
	      
	      if(tmp < cost[q]){
		Heap::Update_MinPolicy(Q, q, tmp);
		pred->data[q] = p;
	      }
	    }
	  }
	}
      }
      free(Dpq);
      glip::FreeFloatArray(&cost);
      Heap::Destroy(&Q);
      return pred;
    }
    

    //---------------------------------------

    void ORFC(sImageGraph *sg,
	      int *S,
	      sImage32 *label){
      glip::sImage32 *tmp,*value_e;
      glip::sPQueue32 *QS=NULL;
      glip::sQueue *Q = glip::Queue::Create(sg->ncols*sg->nrows);
      glip::sAdjRel *A = sg->A;
      int p,q,j,i,k,energy,nsi;
      glip::Pixel u,v;
      int *s_energy,*s_pixel;

      nsi = 0;
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(label->data[p] > 0)
	  nsi++;
      }
      s_energy = glip::AllocIntArray(nsi);
      s_pixel  = glip::AllocIntArray(nsi);
      QS = glip::PQueue32::Create(sg->Wmax+2, nsi, s_energy);
  
      glip::ImageGraph::Transpose(sg);
      
      tmp = glip::Image32::Clone(label);
      value_e = glip::ift::Cost_fmin(sg, S, 0, tmp);
      glip::Image32::Destroy(&tmp);

      glip::ImageGraph::Transpose(sg);
      
      k = 0;
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(label->data[p] > 0){
	  s_energy[k] = value_e->data[p];
	  s_pixel[k] = p;
	  glip::PQueue32::InsertElem(&QS, k);
	  k++;
	  label->data[p] = NIL;
	}
      }

      while(!glip::PQueue32::IsEmpty(QS)){
	j = glip::PQueue32::RemoveMinFIFO(QS);
	p = s_pixel[j];
	energy = s_energy[j];
	
	if(label->data[p] != NIL) continue;
	
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
	    //if(glip::Image32::IsValidPixel(label, v.x, v.y)){
	    if(v.x >= 0 && v.x < label->ncols &&
	       v.y >= 0 && v.y < label->nrows){	    
	      q = v.x + v.y*sg->ncols;
	      if(energy < (sg->n_link[p])[i] &&
		 label->data[q] == NIL){
		label->data[q] = 1;
		glip::Queue::Push(Q, q);
	      }
	    }
	  }
	}
      }
      for(p=0; p<sg->ncols*sg->nrows; p++){
	if(label->data[p] == NIL)
	  label->data[p] = 0;
      }

      glip::Queue::Destroy(&Q);
      glip::Image32::Destroy(&value_e);
      glip::PQueue32::Destroy(&QS);
      glip::FreeIntArray(&s_energy);
      glip::FreeIntArray(&s_pixel);      
    }



    //Outer Cut (the polarity must be embedded in the graph):
    void ORFC_MaxMin(sGraph *graph,
		     sGraph *transpose,
		     int *S,
		     int *label){
      int *tmp,*value_e;
      glip::sPQueue32 *QS=NULL;
      glip::sQueue *Q = glip::Queue::Create(graph->nnodes);
      int p,q,j,i,k,energy,nsi,Wmax;
      int *s_energy,*s_pixel;
      
      nsi = 0;
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(label[p] > 0)
	  nsi++;
      }
      s_energy = glip::AllocIntArray(nsi);
      s_pixel  = glip::AllocIntArray(nsi);
      Wmax = Graph::GetMaximumArc(graph);
      QS = glip::PQueue32::Create(Wmax+2, nsi, s_energy);

      tmp = (int *)calloc(graph->nnodes, sizeof(int));
      memcpy(tmp, label, graph->nnodes*sizeof(int));
      value_e = glip::ift::Cost_fmax(transpose, S, 0, tmp);
      free(tmp);

      k = 0;
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(label[p] > 0){
	  s_energy[k] = value_e[p];
	  s_pixel[k] = p;
	  glip::PQueue32::InsertElem(&QS, k);
	  k++;
	  label[p] = NIL;
	}
      }

      while(!glip::PQueue32::IsEmpty(QS)){
	j = glip::PQueue32::RemoveMaxFIFO(QS);
	p = s_pixel[j];
	energy = s_energy[j];

	if(label[p] != NIL) continue;
	
	glip::Queue::Reset(Q);
	glip::Queue::Push(Q, p);
	label[p] = 1;
	
	while(!glip::Queue::IsEmpty(Q)){
	  p = glip::Queue::Pop(Q);
	  for(i = 0; i < graph->nodes[p].outdegree; i++){
	    q = graph->nodes[p].adjList[i];
	    if(energy > graph->nodes[p].Warcs[i] &&
	       label[q] == NIL){
	      label[q] = 1;
	      glip::Queue::Push(Q, q);
	    }
	  }
	}
      }
      for(p = 0; p < graph->nnodes; p++){
	if(label[p] == NIL)
	  label[p] = 0;
      }
      free(value_e);
      glip::Queue::Destroy(&Q);
      glip::PQueue32::Destroy(&QS);
      glip::FreeIntArray(&s_energy);
      glip::FreeIntArray(&s_pixel);
    }
    

    //---------------------------------------

    //Outer Cut:
    void OIFT(sImage32 *W,
	      sAdjRel *A,
	      sImage32 *img,
	      float per,
	      int *S,
	      sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w,Wmax;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      float per_pq;
      
      value = Image32::Create(W->ncols,
			      W->nrows);
      n = label->n;
      Wmax = glip::Image32::GetMaxVal(W)*2;
      Wmax *= (1.0 + fabsf(per)/100.0);     
      Q = PQueue32::Create(Wmax+2,n,value->data);

      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){

	      w = W->data[p] + W->data[q];
	      if(label->data[p] > 0)
		per_pq = per;
	      else
		per_pq = -per;

	      if(img->data[p] > img->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(img->data[p] < img->data[q])
		w *= (1.0 - per_pq/100.0);
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
    }




    //Outer Cut:
    void OIFT(sAdjRel3 *A,
	      sScene32 *scn,
	      float per,
	      int *S,
	      sScene32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w,Wmax;
      sScene32 *value;
      glip::Voxel u,v;
      float per_pq;
      
      value = Scene32::Create(scn);
      n = label->n;
      Wmax = glip::Scene32::GetMaximumValue(scn);
      Wmax *= (1.0 + fabsf(per)/100.0);     
      Q = PQueue32::Create(Wmax+2,n,value->data);

      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.c.x = glip::Scene32::GetAddressX(label, p);
	u.c.y = glip::Scene32::GetAddressY(label, p);
	u.c.z = glip::Scene32::GetAddressZ(label, p);	
	
	for(i=1; i<A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(glip::Scene32::IsValidVoxel(label, v)){
	    q = glip::Scene32::GetVoxelAddress(label,v);
	    if(Q->L.elem[q].color != BLACK){

	      w = abs(scn->data[p] - scn->data[q]);
	      if(label->data[p] > 0)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Scene32::Destroy(&value);
      PQueue32::Destroy(&Q);
    }

    
    
    
    //Outer Cut:
    void OIFT(sScene32 *W,
	      sAdjRel3 *A,
	      sScene32 *scn,
	      float per,
	      int *S,
	      sScene32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w,Wmax;
      sScene32 *value;
      glip::Voxel u,v;
      float per_pq;
      
      value = Scene32::Create(W);
      n = label->n;
      Wmax = glip::Scene32::GetMaximumValue(W)*2;
      Wmax *= (1.0 + fabsf(per)/100.0);     
      Q = PQueue32::Create(Wmax+2,n,value->data);

      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.c.x = glip::Scene32::GetAddressX(label, p);
	u.c.y = glip::Scene32::GetAddressY(label, p);
	u.c.z = glip::Scene32::GetAddressZ(label, p);	
	
	for(i=1; i<A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(glip::Scene32::IsValidVoxel(label, v)){
	    q = glip::Scene32::GetVoxelAddress(label,v);
	    if(Q->L.elem[q].color != BLACK){

	      w = W->data[p] + W->data[q];
	      if(label->data[p] > 0)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Scene32::Destroy(&value);
      PQueue32::Destroy(&Q);
    }



    //Outer Cut:
    void OIFT(sScene32 *Wx,
	      sScene32 *Wy,
	      sScene32 *Wz,
	      sScene32 *scn,
	      float per,
	      int *S,
	      sScene32 *label){
      sPQueue32 *Q=NULL;
      sAdjRel3 *A;
      int i,j,p,q,n,t;
      int w,Wmax;
      sScene32 *value;
      sScene32 *W[3];
      glip::Voxel u,v;
      float per_pq;
      int T[7];
      W[0] = Wx;
      W[1] = Wy;
      W[2] = Wz;
      A = glip::AdjRel3::Spheric(1.0);
      for(i=1; i<A->n; i++){
	if(A->d[i].axis.x != 0) T[i] = 0;
	if(A->d[i].axis.y != 0) T[i] = 1;
	if(A->d[i].axis.z != 0) T[i] = 2;
      }
      value = Scene32::Create(scn);
      n = label->n;
      Wmax = MAX(glip::Scene32::GetMaximumValue(Wx),
		 MAX(glip::Scene32::GetMaximumValue(Wy),
		     glip::Scene32::GetMaximumValue(Wz)))*2;
      Wmax *= (1.0 + fabsf(per)/100.0);     
      Q = PQueue32::Create(Wmax+2,n,value->data);

      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.c.x = glip::Scene32::GetAddressX(label, p);
	u.c.y = glip::Scene32::GetAddressY(label, p);
	u.c.z = glip::Scene32::GetAddressZ(label, p);	
	
	for(i=1; i<A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(glip::Scene32::IsValidVoxel(label, v)){
	    q = glip::Scene32::GetVoxelAddress(label,v);
	    if(Q->L.elem[q].color != BLACK){
	      t = T[i];
	      w = W[t]->data[p] + W[t]->data[q];
	      if(label->data[p] > 0)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Scene32::Destroy(&value);
      PQueue32::Destroy(&Q);
      AdjRel3::Destroy(&A);
    }

    
    //---------------------------------------
    
    void OIFT_in(sImageGraph *sg,
		 int *S,
		 sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->n;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);
      A = sg->A;

      i_inv = glip::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }

      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){	  
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p] != 0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }
    
    
    void OIFT(sImageGraph *sg,
	      int *S,
	      sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);
      A = sg->A;

      i_inv = glip::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }


    void OIFT_MaxMin(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *value,
		     int niter){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,it = 0;
      int w;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      Image32::Set(pred, NIL);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax*2+2+1,n,value->data);
      A = sg->A;

      i_inv = glip::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = sg->Wmax*2+2; //sg->Wmax+1; //INT_MAX;
	else                    value->data[p] = 0; 
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  PQueue32::FastInsertElem(Q, S[i]);
	  //label->data[S[i]] += 2;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p]!=NIL){
	    PQueue32::FastInsertElem(Q, p);
	    //label->data[p] += 2;
	  }
	}
      }

      while(!PQueue32::IsEmpty(Q)) {
	if(it == niter)
	  break;

	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	//label->data[p] -= 2;
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
		w = MAX(value->data[p], w*2);
	      }
	      else{
		w = (sg->n_link[p])[i];
		w = MAX(value->data[p], w*2+1);
	      }
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		pred->data[q] = p;
		PQueue32::FastInsertElem(Q, q);
		//label->data[q] += 2;
	      }
	    }
	  }
	}

	it++;
      }
      /*
      for(u_y = 0; u_y < label->nrows; u_y++){
	for(u_x = 0; u_x < label->ncols; u_x++)
	  printf("%d ", label->array[u_y][u_x]);
	printf("\n");
      }
      */      
      PQueue32::Destroy(&Q);
      free(i_inv);
    }

    

    void OIFT_MinMax(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *value,
		     int niter){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,it = 0;
      int w;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      Image32::Set(pred, NIL);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);
      A = sg->A;

      i_inv = glip::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = 0;
	else                    value->data[p] = sg->Wmax+1; //INT_MAX;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  PQueue32::FastInsertElem(Q, S[i]);
	  label->data[S[i]] += 2;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p]!=NIL){
	    PQueue32::FastInsertElem(Q, p);
	    label->data[p] += 2;
	  }
	}
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMaxFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	label->data[p] -= 2;
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w > value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		pred->data[q] = p;
		PQueue32::FastInsertElem(Q, q);
		label->data[q] += 2;
	      }
	    }
	  }
	}

	it++;
	if(it == niter)
	  break;
      }
      PQueue32::Destroy(&Q);
      free(i_inv);
    }



    void OIFT_TZ2Bkg(sImageGraph *sg,
		     int *S,
		     sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax*2+3,n,value->data);
      A = sg->A;

      i_inv = glip::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j] * 2;
	      }
	      else
		w = (sg->n_link[p])[i] * 2 + 1;

	      w = MAX(w, value->data[p]);
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }



    void OIFT_TZ2Obj(sImageGraph *sg,
		     int *S,
		     sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax*2+3,n,value->data);
      A = sg->A;

      i_inv = glip::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j] * 2 + 1;
	      }
	      else
		w = (sg->n_link[p])[i] * 2;

	      w = MAX(w, value->data[p]);
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }



    

    void OIFT_TZ(sImageGraph *sg,
		 int *S,
		 sImage32 *label){
      sImage32 *lb_tzbkg, *lb_tzobj;
      int p;
      lb_tzbkg = glip::Image32::Clone(label);
      OIFT_TZ2Bkg(sg, S, lb_tzbkg);
      glip::ImageGraph::Transpose(sg);
      lb_tzobj = glip::Image32::Create(label);
      for(p = 0; p < label->n; p++){
	if(label->data[p] == NIL)
	  lb_tzobj->data[p] = NIL;
	else if(label->data[p] == 0)
	  lb_tzobj->data[p] = 1;
	else
	  lb_tzobj->data[p] = 0;
      }
      OIFT_TZ2Bkg(sg, S, lb_tzobj);
      for(p = 0; p < label->n; p++){
	if(lb_tzbkg->data[p] == 1)
	  label->data[p] = 1;
	else if(lb_tzobj->data[p] == 1)
	  label->data[p] = 0;
	else
	  label->data[p] = 2; /*Tie-zone*/
      }
      glip::ImageGraph::Transpose(sg);
      glip::Image32::Destroy(&lb_tzbkg);
      glip::Image32::Destroy(&lb_tzobj);
    }



    bool isOIFT(sImageGraph *sg,
		int *S,
		sImage32 *Slabel,
		sImage32 *label,
		sImage32 *pred,
		sImage32 *ord,
		bool complete_check){
      bool seg, forest;
      seg = isOIFT_Segmentation(sg,
				S,
				Slabel,
				label);
      if(complete_check)
	forest = isOIFT_Forest(sg,
			       S,
			       Slabel,
			       pred,
			       ord);
      else
	forest = isForest(sg,
			  pred);
      return (seg && forest);
    }

    
    
    bool isOIFT_Segmentation(sImageGraph *sg,
			     int *S,
			     sImage32 *Slabel,
			     sImage32 *label){
      //*************
      //static int i = 0;
      char filename[512];
      //*************
      glip::sImage32 *tz;
      bool flag = true;
      bool energy_test;
      int p;
      tz = glip::Image32::Clone(Slabel);
      OIFT_TZ(sg, S, tz);
      //*************
      //sprintf(filename, "tiezone.pgm");
      //Image32::Write(tz, filename);
      //*************
      for(p = 0; p < label->n; p++){
	if(label->data[p] != tz->data[p] &&
	   tz->data[p] != 2){
	  flag = false;
	  break;
	}
      }
      energy_test = (GetEnergy_Min(sg, tz, 1) == GetEnergy_Min(sg, label, 1));
      printf("within tie-zone: %d, energy test: %d\n", flag, energy_test);
      glip::Image32::Destroy(&tz);
      return (flag && energy_test);
    }


    

    bool isOIFT_Forest_tmp(sImageGraph *sg,
			   sImage32  *label,
			   sImage32  *Tpred,
			   sImage32  *value,
			   sPQueue32 *Q,
			   int *i_inv,
			   sImage32  *pred,
			   sImage32 *ord){
      struct node_oift_info { int label; int pred; int value; int color; };
      struct node_oift_info *backup;
      int u_x,u_y,v_x,v_y;
      int i,j,k,p,q,n,np,x;
      int w,bucket;
      bool flag = true;
      sAdjRel *A;
      int *F = NULL;
      A = sg->A;
      n = label->n;
      if(PQueue32::IsEmpty(Q))
	return true;
      
      bucket = Q->C.minvalue;
      while(Q->C.first[bucket] == NIL)
	bucket++;
      
      np = 0;
      p = Q->C.first[bucket];
      do{
	if(Tpred->data[p] == pred->data[p])
	  np++;
	p = Q->L.elem[p].next;
      }while(p != NIL);

      //printf("np: %d\n", np);
      
      if(np == 0)
	return false;
      
      F = (int *)calloc(np, sizeof(int));	
      np = 0;
      p = Q->C.first[bucket];
      w = value->data[p];
      do{
	if(Tpred->data[p] == pred->data[p]){
	  F[np] = p;
	  np++;
	}
	p = Q->L.elem[p].next;
      }while(p != NIL);	

      //ordenacao por insercao:
      for(i = 0; i < np-1; i++){
	// Insere F[i+1] em F[0],...,F[i].
	x = F[i+1];
	j = i;
	while( j >= 0 && ord->data[F[j]] > ord->data[x]){
	  F[j+1] = F[j];
	  j -= 1;
	}
	F[j+1] = x;
      }
      
      backup = (struct node_oift_info *)calloc(A->n,
					       sizeof(struct node_oift_info));
      
      for(k = 0; k < np; k++){
	p = F[k];
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	//Backup:
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    backup[i].value = value->data[q];
	    backup[i].label = label->data[q];
	    backup[i].pred  = Tpred->data[q];
	    backup[i].color = Q->L.elem[q].color;
	  }
	}

	PQueue32::FastRemoveElem(Q, p);
	flag = true;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		Tpred->data[q] = p;
		PQueue32::FastInsertElem(Q, q);
	      }
	      
	      if(Tpred->data[q] != p && p == pred->data[q]){
		flag = false;
		break;
	      }
	      
	    }
	  }
	}

	if(flag)
	  flag = isOIFT_Forest_tmp(sg, label, Tpred, value,
				   Q, i_inv, pred, ord);
	
	if(flag)
	  break;

	//Restore:
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    label->data[q] = backup[i].label;
	    Tpred->data[q] = backup[i].pred;
	    if(backup[i].color == WHITE){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      value->data[q] = backup[i].value;
	    }
	    else if(backup[i].color == GRAY){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      value->data[q] = backup[i].value;
	      PQueue32::FastInsertElem(Q, q);
	    }
	    Q->L.elem[q].color = backup[i].color;
	  }
	}
	PQueue32::FastInsertElem(Q, p);
      }
      free(F);
      free(backup);
      
      return flag;
    }

    

    
    bool isOIFT_Forest(sImageGraph *sg,
		       int *S,
		       sImage32 *Slabel,
		       sImage32 *pred,
		       sImage32 *ord){
      sPQueue32 *Q=NULL;
      int i,p,n;
      sImage32 *value, *label, *Tpred;
      int *i_inv;
      bool flag = true;

      if(!isForest(sg, pred))
	return false;
      
      label = Image32::Clone(Slabel);
      Tpred = Image32::Create(sg->ncols,
			      sg->nrows);
      Image32::Set(Tpred, NIL);
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->n;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);

      i_inv = glip::AdjRel::InverseIndexes(sg->A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      flag = isOIFT_Forest_tmp(sg,
			       label,
			       Tpred,
			       value,
			       Q,
			       i_inv,
			       pred,
			       ord);

      Image32::Destroy(&value);
      Image32::Destroy(&label);
      Image32::Destroy(&Tpred);
      PQueue32::Destroy(&Q);
      free(i_inv);      
      return flag;
    }



    bool isForest(sImageGraph *sg,
		  sImage32 *pred){
      glip::sBMap *path,*acyclic;
      bool forest = true;
      int p,q;
      acyclic = glip::BMap::Create(pred->n);
      glip::BMap::Fill(acyclic, 0);
      path = glip::BMap::Create(pred->n);
      glip::BMap::Fill(path, 0);      
      for(p = 0; p < pred->n; p++){
	q = p;
	do{
	  if(glip::BMap::Get(path, q) == 1){
	    forest = false;
	    break;
	  }
	  else if(glip::BMap::Get(acyclic, q) == 1)
	    break;

	  glip::BMap::Set1(path, q);
	  
	  q = pred->data[q];
	}while(q != NIL);

	if(!forest) break;

	q = p;
	do{
	  if(glip::BMap::Get(acyclic, q) == 1)
	    break;
	  glip::BMap::Set0(path, q);
	  glip::BMap::Set1(acyclic, q);	  
	  q = pred->data[q];
	}while(q != NIL);
      }      
      glip::BMap::Destroy(&acyclic);
      glip::BMap::Destroy(&path);

      if(!forest)
	printf("Cycle detected.\n");
      
      return forest;
    }




    void OIFT_guided(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *ord){
      sHeap32fi_lex *Q=NULL;
      int i,j,p,q,n;
      int w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = (float *)calloc(label->n, sizeof(float));
      n = label->ncols*label->nrows;
      Q = glip::Heap32fi_lex::Create(label->n, value, ord->data);
      A = sg->A;

      i_inv = glip::AdjRel::InverseIndexes(A);

      Image32::Set(pred, NIL);
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value[p] = FLT_MAX;
	else                    value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  glip::Heap32fi_lex::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    glip::Heap32fi_lex::Insert_MinPolicy(Q, p);
      }

      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->color[q] != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i];
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value[q]){
		//value[q] = w;
		label->data[q] = label->data[p];
		pred->data[q] = p;
		glip::Heap32fi_lex::Update_MinPolicy(Q, q, w, ord->data[q]);
	      }
	    }
	  }
	}
      }
      free(value);
      Heap32fi_lex::Destroy(&Q);
      free(i_inv);
    }


    
    
    
    void OIFT(sGraph *graph,
	      sGraph *transpose,
	      int *S,
	      int *label){
      sPQueue32 *Q=NULL;
      sGraph *g;
      int i,j,p,q,n;
      int w;
      int *value;
      int Wmax;
      Wmax = Graph::GetMaximumArc(graph);
      n = graph->nnodes;
      value = (int *)malloc(n*sizeof(int));
      Q = PQueue32::Create(Wmax+2, n, value);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) value[p] = INT_MAX;
	else              value[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);

	if(label[p]==0) g = transpose;
	else   	        g = graph;

	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  if(Q->L.elem[q].color != BLACK){

	    /*
	    if(label[p]==0)
	      w = Graph::GetArcWeight(graph, q, p);
	    else
            */
	    w = g->nodes[p].Warcs[i];
	    
	    if(w < value[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      value[q] = w;
	      label[q] = label[p];
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      free(value);
      PQueue32::Destroy(&Q);
    }
    



    void OIFT_TZ2Bkg(sGraph *graph,
		     sGraph *transpose,
		     int *S,
		     int *label){
      sPQueue32 *Q=NULL;
      sGraph *g;
      int i,j,p,q,n;
      int w,d;
      int *value;
      int Wmax;
      Wmax = Graph::GetMaximumArc(graph);
      n = graph->nnodes;
      value = (int *)malloc(n*sizeof(int));
      Q = PQueue32::Create(Wmax*2+3, n, value);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) value[p] = INT_MAX;
	else              value[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);

	if(label[p]==0){ g = transpose; d = 0; }
	else{   	 g = graph;     d = 1; }

	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  if(Q->L.elem[q].color != BLACK){

	    /*
	    if(label[p]==0)
	      w = Graph::GetArcWeight(graph, q, p);
	    else
            */
	    w = g->nodes[p].Warcs[i] * 2 + d;

	    w = MAX(w, value[p]);
	    
	    if(w < value[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      value[q] = w;
	      label[q] = label[p];
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      free(value);
      PQueue32::Destroy(&Q);
    }




    void OIFT_Heap(sImageGraph *sg,
		   int *S,
		   sImage32 *label){
      sHeap *Q=NULL;
      int i,j,p,q,n;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      n = label->ncols*label->nrows;
      value = glip::AllocFloatArray(n);
      Q = Heap::Create(n, value);
      A = sg->A;

      i_inv = glip::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value[p] = FLT_MAX;
	else                    value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    Heap::Insert_MinPolicy(Q, p);
      }

      while(!Heap::IsEmpty(Q)) {
	Heap::Remove_MinPolicy(Q, &p);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->color[q] != BLACK){
      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		wi = (sg->n_link[q])[j];
	      }
	      else
		wi = (sg->n_link[p])[i];

	      if(wi == INT_MAX)
		w = FLT_MAX;
	      else
		w = (float)wi;
	      
	      if(w < value[q]){
		label->data[q] = label->data[p];
		Heap::Update_MinPolicy(Q, q, w);
	      }
	    }
	  }
	}
      }
      glip::FreeFloatArray(&value);
      Heap::Destroy(&Q);
      free(i_inv);
    }



    void OIFT_Heap(sGraph *graph,
		   sGraph *transpose,
		   int *S,
		   int *label){
      sHeap *Q=NULL;
      sGraph *g;
      int i,j,p,q,n;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;

      n = graph->nnodes;
      value = glip::AllocFloatArray(n);
      Q = Heap::Create(n, value);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) value[p] = FLT_MAX;
	else              value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    Heap::Insert_MinPolicy(Q, p);
      }

      while(!Heap::IsEmpty(Q)) {
	Heap::Remove_MinPolicy(Q, &p);

	if(label[p]==0) g = transpose;
	else   	        g = graph;
	
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  if(Q->color[q] != BLACK){

	    /*
	    if(label[p]==0)
	      wi = Graph::GetArcWeight(graph, q, p);
	    else
	    */
	    wi = g->nodes[p].Warcs[i];
	    
	    if(wi == INT_MAX)
	      w = FLT_MAX;
	    else
	      w = (float)wi;
	    
	    if(w < value[q]){
	      label[q] = label[p];
	      Heap::Update_MinPolicy(Q, q, w);
	    }
	  }
	}
      }
      glip::FreeFloatArray(&value);
      Heap::Destroy(&Q);
    }
    
    

    /*
    void EOIFT(sImageGraph *sg,
	       int *S,
	       sImage32 *label){
      sPQueue32 *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n;
      int l_ant,e_obj,e_bkg,e_max;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Qobj = PQueue32::Create(sg->Wmax+2,n,value->data);
      Qbkg = PQueue32::Create(sg->Wmax+2,n,value->data);
      Q	= glip::AllocIntArray(n);
      A = sg->A;

      i_inv = glip::AdjRel::InverseIndexes(A);

      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value->data[S[i]] = 0;
	  if(label->data[S[i]] == 0)
	    PQueue32::FastInsertElem(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    PQueue32::FastInsertElem(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value->data[p] = 0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qobj, p);
	  }
	  else{
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(!PQueue32::IsEmpty(Qobj) && !PQueue32::IsEmpty(Qbkg)) {
	p_obj = PQueue32::FastGetMinFIFO(Qobj);
	p_bkg = PQueue32::FastGetMinFIFO(Qbkg);

	e_obj = value->data[p_obj];
	e_bkg = value->data[p_bkg];
	if(e_obj < e_bkg){
	  e_max = e_bkg;
	  p = p_obj;
	  PQueue32::FastRemoveElem(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = e_obj;
	  p = p_bkg;
	  PQueue32::FastRemoveElem(Qbkg, p);
	}
	else{
	  e_max = e_obj;
          if(l_ant == 0){
	    p = p_obj;
	    PQueue32::FastRemoveElem(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    PQueue32::FastRemoveElem(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}

	Qobj->L.elem[p].color = BLACK;
	Qbkg->L.elem[p].color = BLACK;

	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){
		
		if(label->data[p]==0){
		  j = i_inv[i];
		  w = (sg->n_link[q])[j];
		}
		else
		  w = (sg->n_link[p])[i];

                if(w < e_max){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  else if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = BLACK;
		  Qbkg->L.elem[q].color = BLACK;
		  label->data[q] = label->data[p];
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value->data[q]){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = WHITE;
		  Qbkg->L.elem[q].color = WHITE;
		  value->data[q] = w;
		  label->data[q] = label->data[p];
		  if(label->data[q] > 0)
		    PQueue32::FastInsertElem(Qobj, q);
		  else
		    PQueue32::FastInsertElem(Qbkg, q);
		}
	      }
	    }
	  }
	}
      }

      while(!PQueue32::IsEmpty(Qobj)){
	p = PQueue32::FastRemoveMinFIFO(Qobj);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){

		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		Qobj->L.elem[q].color = BLACK;
		label->data[q] = label->data[p];
		Qtop++;
		Q[Qtop] = q;
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Qobj);
      PQueue32::Destroy(&Qbkg);
      glip::FreeIntArray(&Q);
      free(i_inv);
    }
    */

    void EOIFT(sImageGraph *sg,
	       int *S,
	       sImage32 *label,
	       int e_max){
      sPQueue32 *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant,e_obj,e_bkg;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Qobj = PQueue32::Create(sg->Wmax+2,n,value->data);
      Qbkg = PQueue32::Create(sg->Wmax+2,n,value->data);
      Q	= glip::AllocIntArray(n);
      A = sg->A;

      i_inv = glip::AdjRel::InverseIndexes(A);

      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value->data[S[i]] = 0;
	  if(label->data[S[i]] == 0)
	    PQueue32::FastInsertElem(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    PQueue32::FastInsertElem(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value->data[p] = 0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qobj, p);
	  }
	  else{
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      //while(!PQueue32::IsEmpty(Qobj) && !PQueue32::IsEmpty(Qbkg)){
      while(Qobj->nadded != 0 && Qbkg->nadded != 0){
	p_obj = PQueue32::FastGetMinFIFO(Qobj);
	p_bkg = PQueue32::FastGetMinFIFO(Qbkg);

	e_obj = value->data[p_obj];
	e_bkg = value->data[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  p = p_obj;
	  PQueue32::FastRemoveElem(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  p = p_bkg;
	  PQueue32::FastRemoveElem(Qbkg, p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0){
	    p = p_obj;
	    PQueue32::FastRemoveElem(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    PQueue32::FastRemoveElem(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}

	Qobj->L.elem[p].color = BLACK;
	Qbkg->L.elem[p].color = BLACK;

	lp = label->data[p];

	goto label04;
	
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;

	label04:

	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  w = (sg->n_link[q])[j];
		}
		else
		  w = (sg->n_link[p])[i];

                if(w < e_max){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  else if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = BLACK;
		  Qbkg->L.elem[q].color = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value->data[q]){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = WHITE;
		  Qbkg->L.elem[q].color = WHITE;
		  value->data[q] = w;
		  label->data[q] = lp;
		  if(lp > 0)
		    PQueue32::FastInsertElem(Qobj, q);
		  else
		    PQueue32::FastInsertElem(Qbkg, q);
		}
	      }
	    }
	  }
	}
      }

      //while(!PQueue32::IsEmpty(Qobj)){
      while(Qobj->nadded != 0){
	p = PQueue32::FastRemoveMinFIFO(Qobj);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){
		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		Qobj->L.elem[q].color = BLACK;
		label->data[q] = label->data[p];
		Qtop++;
		Q[Qtop] = q;
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Qobj);
      PQueue32::Destroy(&Qbkg);
      glip::FreeIntArray(&Q);
      free(i_inv);
    }
    


    //------------

    void EOIFT(sGraph *graph,
	       sGraph *transpose,
	       int *S,
	       int *label,
	       int e_max){
      sPQueue32 *Qobj=NULL, *Qbkg=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant,e_obj,e_bkg;
      int w;
      int *value;
      int *Q=NULL;
      int Qtop = -1;
      int Wmax;
      Wmax = Graph::GetMaximumArc(graph);
      n = graph->nnodes;
      value = (int *)malloc(n*sizeof(int));

      Qobj = PQueue32::Create(Wmax+2, n, value);
      Qbkg = PQueue32::Create(Wmax+2, n, value);
      Q	= glip::AllocIntArray(n);

      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0;
	  if(label[S[i]] == 0)
	    PQueue32::FastInsertElem(Qbkg, S[i]);
	  else if(label[S[i]] != NIL)
	    PQueue32::FastInsertElem(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    value[p] = INT_MAX;
	    label[p] = 0;
	  }
	  else
	    value[p] = 0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    value[p] = 0;
	    PQueue32::FastInsertElem(Qbkg, p);
	  }
	  else if(label[p] != NIL){
	    value[p] = 0;
	    PQueue32::FastInsertElem(Qobj, p);
	  }
	  else{
	    value[p] = INT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      //while(!PQueue32::IsEmpty(Qobj) && !PQueue32::IsEmpty(Qbkg)){
      while(Qobj->nadded != 0 && Qbkg->nadded != 0){
	p_obj = PQueue32::FastGetMinFIFO(Qobj);
	p_bkg = PQueue32::FastGetMinFIFO(Qbkg);

	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  p = p_obj;
	  PQueue32::FastRemoveElem(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  p = p_bkg;
	  PQueue32::FastRemoveElem(Qbkg, p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0){
	    p = p_obj;
	    PQueue32::FastRemoveElem(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    PQueue32::FastRemoveElem(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}
	
	Qobj->L.elem[p].color = BLACK;
	Qbkg->L.elem[p].color = BLACK;

	lp = label[p];

	if(lp==0) g = transpose;
	else      g = graph;
	
	goto label05;

	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;

	label05:
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];
	    if(Qobj->L.elem[q].color != BLACK){
	      w = g->nodes[p].Warcs[i];
	      
	      if(w < e_max){
		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		Qobj->L.elem[q].color = BLACK;
		Qbkg->L.elem[q].color = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < value[q]){
		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		if(Qbkg->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qbkg, q);
		Qobj->L.elem[q].color = WHITE;
		Qbkg->L.elem[q].color = WHITE;
		value[q] = w;
		label[q] = lp;
		if(lp > 0)
		  PQueue32::FastInsertElem(Qobj, q);
		else
		  PQueue32::FastInsertElem(Qbkg, q);		  
	      }
	      
	    }
	  }
	}
      }
      
      //while(!PQueue32::IsEmpty(Qobj)){
      while(Qobj->nadded != 0){
	p = PQueue32::FastRemoveMinFIFO(Qobj);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  for(i = 0; i < graph->nodes[p].outdegree; i++){
	    q = graph->nodes[p].adjList[i];
	    
	    if(Qobj->L.elem[q].color != BLACK){
	      
	      if(Qobj->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Qobj, q);
	      Qobj->L.elem[q].color = BLACK;
	      label[q] = label[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      free(value);
      PQueue32::Destroy(&Qobj);
      PQueue32::Destroy(&Qbkg);
      glip::FreeIntArray(&Q);
    }
    

    //----------------------------------------
    /*
    void EOIFT_Heap_2(sImageGraph *sg,
		      int *S,
		      sImage32 *label){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n;
      int l_ant;
      float e_obj,e_bkg,e_max;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      //----------
      n = label->n;
      value = glip::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= glip::AllocIntArray(n);
      A = sg->A;
      
      i_inv = glip::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qobj, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(!Heap::IsEmpty(Qobj) && !Heap::IsEmpty(Qbkg)) {
	Heap::Get_MinPolicy(Qobj, &p_obj);
	Heap::Get_MinPolicy(Qbkg, &p_bkg);

	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = e_bkg;
	  p = p_obj;
	  Heap::Delete_MinPolicy(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = e_obj;
	  p = p_bkg;
	  Heap::Delete_MinPolicy(Qbkg, p);
	}
	else{
	  e_max = e_obj;
          if(l_ant == 0){
	    p = p_obj;
	    Heap::Delete_MinPolicy(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    Heap::Delete_MinPolicy(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}

	//-----------------
	//printf("e_obj: %5d, e_bkg: %5d, ", ROUND(e_obj), ROUND(e_bkg));
	//printf("x: %4d, y: %4d\n", p%label->ncols, p/label->ncols);
	//-----------------
	
	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	//Stack::Push(Q, p);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  //p = Stack::Pop(Q);
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols; //PixelX(label, p);
	  u_y = p/label->ncols; //PixelY(label, p);
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->color[q] != BLACK){
		
		if(label->data[p]==0){
		  j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX)
		  w = FLT_MAX;
		else
		  w = (float)wi;
		
                if(w < e_max){
		  if(Qobj->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qobj, q);
		  else if(Qbkg->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qbkg, q);
		  Qobj->color[q] = BLACK;
		  Qbkg->color[q] = BLACK;
		  label->data[q] = label->data[p];
		  //Stack::Push(Q, q);
		  Qtop++;
		  Q[Qtop] = q;
		  //---------
		  //NbyQ++;
		}
		else if(w < value[q]){
		  label->data[q] = label->data[p];
		  if(label->data[q] > 0){
		    if(Qbkg->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qbkg, q);
		    Heap::Update_MinPolicy(Qobj, q, w);
		  }
		  else{
		    if(Qobj->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qobj, q);
		    Heap::Update_MinPolicy(Qbkg, q, w);
		  }
		}
	      }
	    }
	  }
	}
      }

      while(!Heap::IsEmpty(Qobj)){
	Heap::Remove_MinPolicy(Qobj, &p);
	//Stack::Push(Q, p);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  //p = Stack::Pop(Q);
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols; //PixelX(label, p);
	  u_y = p/label->ncols; //PixelY(label, p);
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->color[q] != BLACK){

		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		Qobj->color[q] = BLACK;
		label->data[q] = label->data[p];
		//Stack::Push(Q, q);
		Qtop++;
		Q[Qtop] = q;
		//--------
		//NbyQ++;
	      }
	    }
	  }
	}
      }

      //printf("NbyQ: %d -> %f\n",NbyQ, (float)NbyQ/(float)n);
      
      glip::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      glip::FreeIntArray(&Q);
      free(i_inv);
    }
    */
    //------------------------------------

    /*
    void EOIFT_Heap(sImageGraph *sg,
		    int *S,
		    sImage32 *label,
		    float e_max){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      //----------
      n = label->n;
      value = glip::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= glip::AllocIntArray(n);

      A = sg->A;
      
      i_inv = glip::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qobj, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(Qobj->last > 0 && Qbkg->last > 0){
      //---------------------
	p_obj = Qobj->pixel[1];
	p_bkg = Qbkg->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  Heap::Remove_MinPolicy(Qobj, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  Heap::Remove_MinPolicy(Qbkg, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0){
	    Heap::Remove_MinPolicy(Qobj, &p);
	  }
          else{
	    Heap::Remove_MinPolicy(Qbkg, &p);
	  }
          l_ant = 1 - l_ant;
	}

	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	lp = label->data[p];
	
	goto label01;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label01:
	  
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->color[q] != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX) continue;

		w = (float)wi;
		
                if(w < e_max){
		  if(Qobj->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qobj, q);
		  else if(Qbkg->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qbkg, q);
		  Qobj->color[q] = BLACK;
		  Qbkg->color[q] = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value[q]){
		  label->data[q] = lp;
		  if(lp > 0){
		    if(Qbkg->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qbkg, q);
		    Heap::Update_MinPolicy(Qobj, q, w);
		  }
		  else{
		    if(Qobj->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qobj, q);
		    Heap::Update_MinPolicy(Qbkg, q, w);
		  }
		}
	      }
	    }
	  }
	}
      }

      for(i = 1; i <= Qobj->last; i++){
	p = Qobj->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	Qobj->color[p] = BLACK;
      }
      
      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;
	u_x = p%label->ncols;
	u_y = p/label->ncols;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Qobj->color[q] != BLACK){
	      Qobj->color[q] = BLACK;
	      label->data[q] = label->data[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      
      glip::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      glip::FreeIntArray(&Q);
      free(i_inv);
    }
    */
    //-----------------------------------------------------


    void EOIFT_Heap(sImageGraph *sg,
		    int *S,
		    sImage32 *label,
		    float e_max){
      sHeapPair *QP=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      n = label->n;
      value = glip::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= glip::AllocIntArray(n);
      A = sg->A;
      
      i_inv = glip::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label->data[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label->data[p] != NIL){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0)
	    HeapPair::Remove_MinPolicy_1(QP, &p);
          else
	    HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	}

	QP->color[p] = BLACK;

	lp = label->data[p];
	
	goto label03;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label03:
	  
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(QP->color[q] != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX) continue;

		w = (float)wi;
		
                if(w < e_max){
		  if(QP->color[q] == GRAY)
		    HeapPair::Delete_MinPolicy(QP, q);
		  QP->color[q] = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value[q]){
		  label->data[q] = lp;
		  if(lp > 0)
		    HeapPair::Update_MinPolicy_1(QP, q, w);
		  else
		    HeapPair::Update_MinPolicy_0(QP, q, w);
		}
	      }
	    }
	  }
	}
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }
      
      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;
	u_x = p%label->ncols;
	u_y = p/label->ncols;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(QP->color[q] != BLACK){
	      QP->color[q] = BLACK;
	      label->data[q] = label->data[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      
      glip::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      glip::FreeIntArray(&Q);
      free(i_inv);
    }
    

    //-----------------------------------------------------
    
    void EOIFT_Heap_2(sImageGraph *sg,
		      int *S,
		      sImage32 *label,
		      float e_max){
      sHeapPair *QP=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      int *T=NULL;
      float *Tv = NULL;
      int Ttop = -1;
      n = label->n;
      value = glip::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= glip::AllocIntArray(n);
      T	= glip::AllocIntArray(n);
      Tv = glip::AllocFloatArray(n);

      A = sg->A;
      i_inv = glip::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  Tv[S[i]] = value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label->data[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    Tv[p] = value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    Tv[p] = value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label->data[p] != NIL){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    Tv[p] = value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	  //------------------
	  //p_obj = QP->pixel[QP->n];
	  //while(QP->last_1 <= QP->n && value[p_obj] < e_max){
	  //  HeapPair::Remove_MinPolicy_1(QP, &p_obj);
	  //  Qtop++;
	  //  Q[Qtop] = p_obj;
	  //  p_obj = QP->pixel[QP->n];
	  //}
	  //------------------
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	  //------------------
	  //p_bkg = QP->pixel[1];
	  //while(QP->last_0 > 0 && value[p_bkg] < e_max){
	  //  HeapPair::Remove_MinPolicy_0(QP, &p_bkg);
	  //  Qtop++;
	  //  Q[Qtop] = p_bkg;
	  //  p_bkg = QP->pixel[1];
	  //}
	  //------------------
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0) HeapPair::Remove_MinPolicy_1(QP, &p);
          else           HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	  //------------------
	  //if(e_obj < e_max){
	  //  p_obj = QP->pixel[QP->n];
	  //  while(QP->last_1 <= QP->n && value[p_obj] < e_max){
	  //    HeapPair::Remove_MinPolicy_1(QP, &p_obj);
	  //    Qtop++;
	  //    Q[Qtop] = p_obj;
	  //    p_obj = QP->pixel[QP->n];
	  //  }
	  //  p_bkg = QP->pixel[1];
	  //  while(QP->last_0 > 0 && value[p_bkg] < e_max){
	  //    HeapPair::Remove_MinPolicy_0(QP, &p_bkg);
	  //    Qtop++;
	  //    Q[Qtop] = p_bkg;
	  //    p_bkg = QP->pixel[1];
	  //  }
	  //}
	  //else e_max = e_obj;
	  //------------------
	}

	//QP->color[p] = BLACK;

	lp = label->data[p];
	
	goto label01;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label01:
	  
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(QP->color[q] != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX) continue;

		w = (float)wi;
		
                if(w < e_max){
		  if(QP->color[q] == GRAY)
		    HeapPair::Delete_MinPolicy(QP, q);
		  QP->color[q] = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < Tv[q]){ //if(w < value[q]){
		  if(Tv[q] >= value[q]){
		    label->data[q] = lp;
		    Ttop++;
		    T[Ttop] = q;
		  }
		  Tv[q] = w;
		}
	      }
	    }
	  }
	}
	while(Ttop > -1){
	  p = T[Ttop];
	  w = Tv[p];
	  Ttop--;
	  if(QP->color[p] != BLACK){
	    if(label->data[p] > 0)
	      HeapPair::Update_MinPolicy_1(QP, p, w);
	    else
	      HeapPair::Update_MinPolicy_0(QP, p, w);
	  }
	}
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }
      
      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;
	u_x = p%label->ncols;
	u_y = p/label->ncols;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(QP->color[q] != BLACK){
	      QP->color[q] = BLACK;
	      label->data[q] = label->data[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      
      glip::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      glip::FreeIntArray(&Q);
      glip::FreeIntArray(&T);
      glip::FreeFloatArray(&Tv);
      free(i_inv);
    }


    //------------------------------------
    /*
    void EOIFT_Heap(sGraph *graph,
		    sGraph *transpose,
		    int *S,
		    int *label){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n;
      int l_ant;
      float e_obj,e_bkg,e_max;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sStack *Q=NULL;
      
      n = graph->nnodes;
      value = glip::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= Stack::Create(n);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) value[p] = FLT_MAX;
	else              value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  if(label[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p] == 0)
	    Heap::Insert_MinPolicy(Qbkg, p);
	  else if(label[p] != NIL)
	    Heap::Insert_MinPolicy(Qobj, p);
      }

      l_ant = 0;
      while(!Heap::IsEmpty(Qobj) && !Heap::IsEmpty(Qbkg)) {
	Heap::Get_MinPolicy(Qobj, &p_obj);
	Heap::Get_MinPolicy(Qbkg, &p_bkg);

	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = e_bkg;
	  p = p_obj;
	}
	else if(e_obj > e_bkg){
	  e_max = e_obj;
	  p = p_bkg;
	}
	else{
	  e_max = e_obj;
          if(l_ant == 0)
	    p = p_obj;
          else
	    p = p_bkg;
          l_ant = 1 - l_ant;
	}

	if(Qobj->color[p] == GRAY)
	  Heap::Delete_MinPolicy(Qobj, p);
	if(Qbkg->color[p] == GRAY)
	  Heap::Delete_MinPolicy(Qbkg, p);
	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	Stack::Push(Q, p);
	while(!Stack::IsEmpty(Q)){
	  p = Stack::Pop(Q);
	  
	  if(label[p]==0) g = transpose;
	  else   	  g = graph;
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];

	    if(Qobj->color[q] != BLACK){

	      
	      //if(label[p]==0)
		//wi = Graph::GetArcWeight(graph, q, p);
	      //else
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX)
		w = FLT_MAX;
	      else
		w = (float)wi;

	      if(w < e_max){
		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		if(Qbkg->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qbkg, q);
		Qobj->color[q] = BLACK;
		Qbkg->color[q] = BLACK;
		
		label[q] = label[p];
		Stack::Push(Q, q);
	      }
	      else if(w < value[q]){
		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		if(Qbkg->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qbkg, q);
		value[q] = w;
		label[q] = label[p];
		if(label[q] > 0)
		  Heap::Insert_MinPolicy(Qobj, q);
		else
		  Heap::Insert_MinPolicy(Qbkg, q);		  
	      }
	      
	    }
	  }

	}
	
      }

      while(!Heap::IsEmpty(Qobj)){
	Heap::Remove_MinPolicy(Qobj, &p);
	Stack::Push(Q, p);
	while(!Stack::IsEmpty(Q)){
	  p = Stack::Pop(Q);

	  for(i = 0; i < graph->nodes[p].outdegree; i++){
	    q = graph->nodes[p].adjList[i];
	    
	    if(Qobj->color[q] != BLACK){

	      if(Qobj->color[q] == GRAY)
		Heap::Delete_MinPolicy(Qobj, q);
	      Qobj->color[q] = BLACK;
	      label[q] = label[p];
	      Stack::Push(Q, q);
	    }
	  }
	}
      }

      while(!Heap::IsEmpty(Qbkg)){
	Heap::Remove_MinPolicy(Qbkg, &p);
	Stack::Push(Q, p);
	Qobj->color[p] = BLACK;
	while(!Stack::IsEmpty(Q)){
	  p = Stack::Pop(Q);

	  for(i = 0; i < transpose->nodes[p].outdegree; i++){
	    q = transpose->nodes[p].adjList[i];
	    
	    if(Qobj->color[q] != BLACK){
	      
	      if(Qbkg->color[q] == GRAY)
		Heap::Delete_MinPolicy(Qbkg, q);
	      Qobj->color[q] = BLACK;
	      label[q] = label[p];
	      Stack::Push(Q, q);
	    }
	  }
	}
      }

      glip::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      Stack::Destroy(&Q);
    }
    */



    /*    
    void EOIFT_Heap(sGraph *graph,
		    sGraph *transpose,
		    int *S,
		    int *label,
		    float e_max){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      int *Q=NULL;
      int Qtop = -1;
      
      n = graph->nnodes;
      value = glip::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= glip::AllocIntArray(n);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qbkg, p);
	  }
	  else if(label[p] != NIL){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qobj, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(Qobj->last > 0 && Qbkg->last > 0){
	p_obj = Qobj->pixel[1];
	p_bkg = Qbkg->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  Heap::Remove_MinPolicy(Qobj, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  Heap::Remove_MinPolicy(Qbkg, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0)
	    Heap::Remove_MinPolicy(Qobj, &p);
          else
	    Heap::Remove_MinPolicy(Qbkg, &p);
          l_ant = 1 - l_ant;
	}

	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	lp = label[p];
	
	if(lp==0) g = transpose;
	else   	  g = graph;
	
	goto label02;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];

	    if(Qobj->color[q] != BLACK){
	      
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX) continue;

	      w = (float)wi;

	      if(w < e_max){
		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		else if(Qbkg->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qbkg, q);
		Qobj->color[q] = BLACK;
		Qbkg->color[q] = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < value[q]){
		label[q] = lp;
		if(lp > 0){
		  if(Qbkg->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qbkg, q);
		  Heap::Update_MinPolicy(Qobj, q, w);
		}
		else{
		  if(Qobj->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qobj, q);
		  Heap::Update_MinPolicy(Qbkg, q, w);
		}
	      }
	      
	    }
	  }

	}
	
      }

      for(i = 1; i <= Qobj->last; i++){
	p = Qobj->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	Qobj->color[p] = BLACK;
      }

      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;      
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Qobj->color[q] != BLACK){
	    Qobj->color[q] = BLACK;
	    label[q] = label[p];
	    Qtop++;
	    Q[Qtop] = q;
	  }
	}
      }
	
      glip::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      glip::FreeIntArray(&Q);
    }
    */



    void EOIFT_Heap(sGraph *graph,
		    sGraph *transpose,
		    int *S,
		    int *label,
		    float e_max){
      sHeapPair *QP=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int *Q=NULL;
      int Qtop = -1;
      n = graph->nnodes;
      value = glip::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= glip::AllocIntArray(n);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label[p] != NIL){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0) HeapPair::Remove_MinPolicy_1(QP, &p);
          else   	 HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	}

	QP->color[p] = BLACK;

	lp = label[p];
	
	if(lp==0) g = transpose;
	else   	  g = graph;
	
	goto label06;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label06:
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];
	    if(QP->color[q] != BLACK){
	      
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX) continue;

	      w = (float)wi;

	      if(w < e_max){
		if(QP->color[q] == GRAY)
		  HeapPair::Delete_MinPolicy(QP, q);
		QP->color[q] = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < value[q]){
		label[q] = lp;
		if(lp > 0)
		  HeapPair::Update_MinPolicy_1(QP, q, w);
		else
		  HeapPair::Update_MinPolicy_0(QP, q, w);
	      }
	    }
	  }
	}
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }

      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;      
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(QP->color[q] != BLACK){
	    QP->color[q] = BLACK;
	    label[q] = label[p];
	    Qtop++;
	    Q[Qtop] = q;
	  }
	}
      }
	
      glip::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      glip::FreeIntArray(&Q);
    }
    

    
    void EOIFT_Heap_2(sGraph *graph,
		      sGraph *transpose,
		      int *S,
		      int *label,
		      float e_max){
      sHeapPair *QP=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int *Q=NULL;
      int Qtop = -1;
      int *T=NULL;
      float *Tv = NULL;
      int Ttop = -1;
      n = graph->nnodes;
      value = glip::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= glip::AllocIntArray(n);
      T	= glip::AllocIntArray(n);
      Tv = glip::AllocFloatArray(n);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  Tv[S[i]] = value[S[i]] = 0.0;
	  if(label[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    Tv[p] = value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	  else
	    Tv[p] = value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label[p] != NIL){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    Tv[p] = value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0) HeapPair::Remove_MinPolicy_1(QP, &p);
          else   	 HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	}

	//QP->color[p] = BLACK;

	lp = label[p];
	
	if(lp==0) g = transpose;
	else   	  g = graph;
	
	goto label02;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label02:
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];
	    if(QP->color[q] != BLACK){
	      
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX) continue;

	      w = (float)wi;

	      if(w < e_max){
		if(QP->color[q] == GRAY)
		  HeapPair::Delete_MinPolicy(QP, q);
		QP->color[q] = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < Tv[q]){
		if(Tv[q] >= value[q]){
		  label[q] = lp;
		  Ttop++;
		  T[Ttop] = q;
		}
		Tv[q] = w;
	      }
	    }
	  }
	}
	while(Ttop > -1){
	  p = T[Ttop];
	  w = Tv[p];
	  Ttop--;
	  if(QP->color[p] != BLACK){
	    if(label[p] > 0)
	      HeapPair::Update_MinPolicy_1(QP, p, w);
	    else
	      HeapPair::Update_MinPolicy_0(QP, p, w);
	  }
	}	
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }

      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;      
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(QP->color[q] != BLACK){
	    QP->color[q] = BLACK;
	    label[q] = label[p];
	    Qtop++;
	    Q[Qtop] = q;
	  }
	}
      }
	
      glip::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      glip::FreeIntArray(&Q);
      glip::FreeIntArray(&T);
      glip::FreeFloatArray(&Tv);
    }
    
    
    
    
    void IFT_fmax_Heap(sGraph *graph,
		       int *S,
		       int *label,
		       float *cost){
      sHeap *Q;
      float tmp, w;
      int n,p,q,i;
      n = graph->nnodes;
      Q = glip::Heap::Create(n, cost);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) cost[p] = FLT_MAX;
	else              cost[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  glip::Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    glip::Heap::Insert_MinPolicy(Q, p);
      }
      
      while(!glip::Heap::IsEmpty(Q)){
	glip::Heap::Remove_MinPolicy(Q, &p);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->color[q] != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = MAX(cost[p], w);

	    if(tmp < cost[q]){
	      glip::Heap::Update_MinPolicy(Q, q, tmp);
	      label[q] = label[p];
	    }
	  }
	}
      }
      glip::Heap::Destroy(&Q);
    }




    void IFT_fmax(sGraph *graph,
		  int *S,
		  int *label,
		  int *cost){
      sPQueue32 *Q;
      int tmp, w, Wmax;
      int n,p,q,i;
      n = graph->nnodes;
      Wmax = Graph::GetMaximumArc(graph);
      Q = PQueue32::Create(Wmax+2, n, cost);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) cost[p] = INT_MAX;
	else              cost[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->L.elem[q].color != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = MAX(cost[p], w);

	    if(tmp < cost[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      cost[q] = tmp;
	      label[q] = label[p];
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
    }



    void IFT_fw(sGraph *graph,
		int *S,
		int *label,
		int *cost){
      sPQueue32 *Q;
      int tmp, w, Wmax;
      int n,p,q,i;
      n = graph->nnodes;
      Wmax = Graph::GetMaximumArc(graph);
      Q = PQueue32::Create(Wmax+2, n, cost);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) cost[p] = INT_MAX;
	else              cost[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->L.elem[q].color != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = w;

	    if(tmp < cost[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      cost[q] = tmp;
	      label[q] = label[p];
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
    }
    
    

    
    void IFT_fw_Heap(sGraph *graph,
		     int *S,
		     int *label,
		     float *cost,
		     int *pred){
      sHeap *Q;
      float tmp, w;
      int n,p,q,i;
      n = graph->nnodes;
      Q = glip::Heap::Create(n, cost);

      for(p = 0; p < n; p++){
	pred[p] = NIL;
	if(label[p]==NIL) cost[p] = FLT_MAX;
	else              cost[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  glip::Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    glip::Heap::Insert_MinPolicy(Q, p);
      }
      
      while(!glip::Heap::IsEmpty(Q)){
	glip::Heap::Remove_MinPolicy(Q, &p);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->color[q] != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = w;

	    if(tmp < cost[q]){
	      glip::Heap::Update_MinPolicy(Q, q, tmp);
	      label[q] = label[p];
	      pred[q] = p;
	    }
	  }
	}
      }
      glip::Heap::Destroy(&Q);
    }
   


    void IFT_fw(sGraph *graph,
		int *S,
		int *label,
		int *cost,
		int *pred){
      sPQueue32 *Q;
      int tmp, w, Wmax;
      int n,p,q,i;
      n = graph->nnodes;
      Wmax = Graph::GetMaximumArc(graph);
      Q = PQueue32::Create(Wmax+2, n, cost);
      
      for(p = 0; p < n; p++){
	pred[p] = NIL;
	if(label[p]==NIL) cost[p] = INT_MAX;
	else              cost[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->L.elem[q].color != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = w;

	    if(tmp < cost[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      cost[q] = tmp;
	      label[q] = label[p];
	      pred[q] = p;
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
    }
   
    
    //---------------------------------------
    // Convex IFT:

    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_conquer_path(int p,
			 sImageGraph *sg,
			 sImage32 *P_sum, 
			 sImage32 *V,
			 sPQueue32 *Q,
			 sImage32 *label){
      int i,q,edge;
      Pixel u,v;
      sAdjRel *A;
      
      A = sg->A;
      do{
	if(Q->L.elem[p].color == GRAY)
	  PQueue32::FastRemoveElem(Q, p);
	Q->L.elem[p].color = BLACK;
	
	label->data[p] = 1;
	
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols*v.y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      edge = (sg->n_link[p])[i];
	      
	      if(edge < V->data[q] && q != P_sum->data[p]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		V->data[q] = edge;
		label->data[q] = 1; //label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
	p = P_sum->data[p];
	if(p == NIL) break;
      }while(Q->L.elem[p].color != BLACK);
    }

    
    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_prune_tree(int p,
		       sImageGraph *sg,
		       sImage32 *P_sum, 
		       sImage32 *V,
		       sPQueue32 *Q,
		       sQueue *Qfifo,
		       sImage32 *label){
      Pixel u,v;
      int i,q,edge;
      sAdjRel *A = sg->A;
      
      if(Q->L.elem[p].color == GRAY)
	PQueue32::FastRemoveElem(Q, p);
      Q->L.elem[p].color = BLACK;
      
      label->data[p] = 0;
      
      Queue::Push(Qfifo, p);
      
      //printf("Prune tree\n");
      
      while(!Queue::IsEmpty(Qfifo)){
	p = Queue::Pop(Qfifo);
	u.x = p%label->ncols; 
	u.y = p/label->ncols; 
	
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v.x,v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols*v.y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(P_sum->data[q] == p){
		Queue::Push(Qfifo, q);
		
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		Q->L.elem[q].color = BLACK;
		
		label->data[q] = 0;
	      }
	      else{
		edge = (sg->n_link[p])[i];
		if(edge < V->data[q]){
		  if(Q->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Q, q);
		  V->data[q] = edge;
		  label->data[q] = 0;
		  PQueue32::FastInsertElem(Q, q);
		}
	      }
	      
	    }
	  }
	}
	
      }
    }

    
    void SC_IFT(sImageGraph *sg,
		int *S,
		sImage32 *label,
		sImage32 *P_sum){
      sPQueue32 *Q=NULL;
      sQueue *Qfifo=NULL;
      sImage32 *V;
      int p,n,i;
      
      n = sg->ncols*sg->nrows;
      V = Image32::Create(sg->ncols, sg->nrows);
      Q = PQueue32::Create(sg->Wmax+2, n, V->data);
      Qfifo = Queue::Create(n);
      
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) V->data[p] = INT_MAX;
	else                    V->data[p] = 0;
      }

      for(i = 1; i <= S[0]; i++)
	PQueue32::FastInsertElem(Q, S[i]);
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	if(label->data[p] > 0)
	  SC_conquer_path(p, sg, P_sum, V, Q, label);
	else if(label->data[p] == 0)
	  SC_prune_tree(p, sg, P_sum, V, Q, Qfifo, label);
      }
      
      Image32::Destroy(&V);
      Queue::Destroy(&Qfifo);
      PQueue32::Destroy(&Q);
    }


    //----------------------------------------

    sImage32 *Cost_fmin(sImageGraph *sg,
			int *S, int lb,
			sImage32 *label){
      sPQueue32 *Q = NULL; 
      sImage32 *V;
      int i,p,q,n, edge,tmp;
      Pixel u,v;
      sAdjRel *A;
      
      n = sg->ncols*sg->nrows;
      V = Image32::Create(sg->ncols, sg->nrows);
      Q = PQueue32::Create(sg->Wmax+2, n, V->data);
      A = sg->A;
      
      for(p=0; p<n; p++){
	if(label->data[p] == lb) V->data[p] = sg->Wmax+1;
	else                     V->data[p] = INT_MIN;
      }

      if(S != NULL){
	for(i = 1; i <= S[0]; i++){
	  if(label->data[S[i]] == lb)
	    PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p] == lb)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMaxFIFO(Q);
	u.x = p%label->ncols; 
	u.y = p/label->ncols; 
	
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols*v.y;
	    if(Q->L.elem[q].color != BLACK){
	      edge = (sg->n_link[p])[i];          
	      tmp  = MIN(V->data[p], edge);
	      if(tmp > V->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		V->data[q] = tmp; //mapa de custos
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
      return V;
    }



    int *Cost_fmax(sGraph *graph,
		   int *S, int lb,
		   int *label){
      sPQueue32 *Q = NULL; 
      int *V;
      int i,p,q,n, edge,tmp,Wmax;
      
      n = graph->nnodes;
      V = (int *)calloc(graph->nnodes, sizeof(int));
      Wmax = Graph::GetMaximumArc(graph);
      Q = PQueue32::Create(Wmax+2, n, V);
      
      for(p=0; p<n; p++){
	if(label[p] == lb) V[p] = 0;
	else               V[p] = INT_MAX;
      }

      if(S != NULL){
	for(i = 1; i <= S[0]; i++){
	  if(label[S[i]] == lb)
	    PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(label[p] == lb)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);

	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  
	  if(Q->L.elem[q].color != BLACK){
	    edge = graph->nodes[p].Warcs[i];
	    tmp  = MAX(V[p], edge);
	    if(tmp < V[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      V[q] = tmp; //mapa de custos
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
      return V;
    }
    
    //----------------------------------------
    
    sImage32 *Cost_fmax(sImageGraph *sg,
			int *S, int lb,
			sImage32 *label){
      sPQueue32 *Q = NULL; 
      sImage32 *V;
      int i,p,q,n, edge,tmp;
      Pixel u,v;
      sAdjRel *A;
      
      n = sg->ncols*sg->nrows;
      V = Image32::Create(sg->ncols, sg->nrows);
      Q = PQueue32::Create(sg->Wmax+2, n, V->data);
      A = sg->A;
      
      for(p=0; p<n; p++){
	if(label->data[p] == lb) V->data[p] = 0;
	else                     V->data[p] = INT_MAX;
      }

      if(S != NULL){
	for(i = 1; i <= S[0]; i++){
	  if(label->data[S[i]] == lb)
	    PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p] == lb)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.x = p%label->ncols; 
	u.y = p/label->ncols; 
	
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols*v.y;
	    if(Q->L.elem[q].color != BLACK){
	      edge = (sg->n_link[p])[i];          
	      tmp  = MAX(V->data[p], edge);
	      if(tmp < V->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		V->data[q] = tmp; //mapa de custos
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
      return V;
    }
    
    
    sSet *COIFT_new_seeds(sImageGraph *sg,
			  int *S,
			  sImage32 *E,
			  sImage32 *label_oift){
      sPQueue32 *Q = NULL;
      sImage32 *V, *pred, *label;
      Pixel u,v;
      sAdjRel *A;
      sSet *newSi = NULL;
      int i,j,p,q,n, tmp, Emax, Emin, e_q, e_p, w_pq; //Wmin;
      
      //Wmin = MinimumWeight(sg);
      Emax = Image32::GetMaxVal(E);
      Emin = E->data[S[1]];
      
      V = Image32::Create(sg->ncols, sg->nrows);
      label = Image32::Create(sg->ncols, sg->nrows);
      pred = Image32::Create(sg->ncols, sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create((Emax+1)+2, n, V->data);
      A = sg->A;
      
      Image32::Set(pred,  NIL);
      Image32::Set(V, INT_MAX);
      Image32::Set(label, NIL);

      for(i = 1; i <= S[0]; i++){
	if(label_oift->data[S[i]] == 1)
	  label->data[S[i]] = 1;
      }

      V->data[S[1]] = 0;
      PQueue32::FastInsertElem(Q, S[1]);
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	
	if(E->data[p] < Emin)
	  Emin = E->data[p];
	
	for(i = 1; i < A->n; i++){
          v.x = u.x + A->dx[i];
          v.y = u.y + A->dy[i];
	  
          //if(Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols * v.y;
            
	    if(Q->L.elem[q].color != BLACK){
	      //------------------------------------
	      e_p = E->data[p];
	      e_q = E->data[q];
	      w_pq = (sg->n_link[p])[i];
	      if(label->data[q] == 1){ //Seed pixel in Si
		tmp = 0;
	      }
	      else if(e_q >= Emin && label_oift->data[q] > 0){
		tmp = 1;
	      }
	      /*
		else if(e_q == e_p){
		tmp = 1;
		}
	      */
	      /*
		else if(e_q > w_pq){
		tmp = 0;
		}
		else if(e_q == w_pq && w_pq > Wmin){
		tmp = 1;
		}
	      */
	      else{
		tmp = Emax - e_q + 2;
	      }
	      //------------------------------------
	      if(tmp < V->data[q]){
		if(Q->L.elem[q].color == GRAY) 
		  PQueue32::FastRemoveElem(Q, q);                     
		V->data[q] = tmp;
		pred->data[q] = p;
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
          }
	}
      }
      
      // calculando os predecessores das sementes Si2, para obter as novas sementes internas
      for(i = 1; i <= S[0]; i++){
	if(label->data[S[i]] == 1){
	  p = S[i];
	  q = pred->data[p];
	  while(q != NIL && (label->data[q] != 1)){ //enquanto o predecessor nao for raiz(i.e.,Si1) e seu rótulo nao for marcado antes(i.e., NIL)
	    Set::Insert(&newSi, q);
	    //printf("newSi->elem: %i\n",newSi->elem);
	    //printf("p: %i\n",q);
	    label->data[q] = 1;
	    q = pred->data[q];
	  }
	}
      }
      
      Image32::Destroy(&V);
      Image32::Destroy(&label);
      Image32::Destroy(&pred);
      PQueue32::Destroy(&Q);
      return newSi;
    }
   

    void COIFT(sImageGraph *sg,
	       int *S,
	       sImage32 *label){
      sImage32 *V_bkg, *Vc_bkg, *ero_V_bkg, *label_oift;
      sSet *Snew=NULL;
      int i_worst,i,n,p,energy,e;
      sAdjRel *A = AdjRel::Circular(2.5);
      
      if(S == NULL) return;
      
      label_oift = Image32::Clone(label);
      OIFT(sg, S, label_oift);
      
      /*1. IFT_max com semente do fundo*/
      ImageGraph::Transpose(sg);

      // mapa de custos via IFT_MAX usando só sementes do fundo
      V_bkg = Cost_fmax(sg, S, 0, label);
      ImageGraph::Transpose(sg);
      
      ero_V_bkg = Image32::Erode(V_bkg, A);
      //Vc_bkg = Image32::Complement(ero_V_bkg);
      
      //Image32::Write(V_bkg, "./out/V_bkg.pgm");
      //Image32::Write(ero_V_bkg, "./out/ero_V_bkg.pgm");
      //Image32::Write(Vc_bkg, "./out/Vc_bkg.pgm");
      
      /*2. IFT com f(pi_s.<s,t>)= Vc_bkg(t) */
      /*
	Si1 = Si;
	Si2 = Si->next;
	Si->next = NULL;
      */
      
      energy = INT_MAX;
      i_worst = 1;
      for(i = 1; i <= S[0]; i++){
	e = ero_V_bkg->data[S[i]];
	if(label->data[S[i]] == 1 && e < energy){
	  energy = e;
	  i_worst = i;
	}
      }

      p = S[1];
      S[1] = S[i_worst];
      S[i_worst] = p;
      
      Snew = COIFT_new_seeds(sg, S, ero_V_bkg, label_oift);
      
      /*3. IFT com novas sementes*/
      Image32::DrawSet(label, Snew,   1);

      OIFT(sg, NULL, label);
      
      AdjRel::Destroy(&A);
      Image32::Destroy(&ero_V_bkg);
      Image32::Destroy(&V_bkg);
      Image32::Destroy(&label_oift);
      Set::Destroy(&Snew);
    }
    

    //----------------------------------------

    sImage32 *BB_Geodesic_Cost(sImage32 *pred,
			       sAdjRel *A){
      sQueue *Qfifo=NULL;
      sImage32 *cost;
      Pixel u,v;
      int n,p,q,i;
      int *Dpq;
      
      Dpq = glip::AllocIntArray(A->n);
      for(i = 1; i < A->n; i++){
	Dpq[i] = ROUND(10*sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]));
      }
      
      n = pred->ncols*pred->nrows;
      cost  = Image32::Create(pred->ncols, pred->nrows);
      Image32::Set(cost, INT_MAX);
      Qfifo = Queue::Create(n);
      
      for(p = 0; p < n; p++){
	if(pred->data[p] == NIL){
	  cost->data[p] = 0;
	  Queue::Push(Qfifo, p);
	}
      }
      
      while(!Queue::IsEmpty(Qfifo)){
	p = Queue::Pop(Qfifo);
	u.x = p%pred->ncols; 
	u.y = p/pred->ncols; 
	
	for(i = 1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(pred, v.x, v.y)){
	  if(v.x >= 0 && v.x < pred->ncols &&
	     v.y >= 0 && v.y < pred->nrows){
	    q = v.x + pred->ncols*v.y;
	    if(p == pred->data[q]){
	      cost->data[q] = cost->data[p] + Dpq[i];
	      Queue::Push(Qfifo, q);
	    }
	  }
	}
      }
      glip::FreeIntArray(&Dpq);
      Queue::Destroy(&Qfifo);
      return cost;
    }



    sScene32 *BB_Geodesic_Cost(sScene32 *pred,
			       sAdjRel3 *A){
      sQueue *Qfifo=NULL;
      sScene32 *cost;
      Voxel u,v;
      int n,p,q,i;
      int *Dpq;
      
      Dpq = glip::AllocIntArray(A->n);
      for(i = 1; i < A->n; i++){
	Dpq[i] = ROUND(10*sqrtf(A->d[i].axis.x*A->d[i].axis.x +
				A->d[i].axis.y*A->d[i].axis.y +
				A->d[i].axis.z*A->d[i].axis.z));
      }
      
      n = pred->n;
      cost  = Scene32::Create(pred);
      Scene32::Fill(cost, INT_MAX);
      Qfifo = Queue::Create(n);
      
      for(p = 0; p < n; p++){
	if(pred->data[p] == NIL){
	  cost->data[p] = 0;
	  Queue::Push(Qfifo, p);
	}
      }
      
      while(!Queue::IsEmpty(Qfifo)){
	p = Queue::Pop(Qfifo);
	u.c.x = glip::Scene32::GetAddressX(pred, p);
	u.c.y = glip::Scene32::GetAddressY(pred, p);
	u.c.z = glip::Scene32::GetAddressZ(pred, p);	
	
	for(i = 1; i < A->n; i++){
	  v.c.x = u.c.x + A->d[i].axis.x;
	  v.c.y = u.c.y + A->d[i].axis.y;
	  v.c.z = u.c.z + A->d[i].axis.z;
	  if(Scene32::IsValidVoxel(pred, v)){
	    q = Scene32::GetVoxelAddress(pred, v);
	    if(p == pred->data[q]){
	      cost->data[q] = cost->data[p] + Dpq[i];
	      Queue::Push(Qfifo, q);
	    }
	  }
	}
      }
      glip::FreeIntArray(&Dpq);
      Queue::Destroy(&Qfifo);
      return cost;
    }


    
    //Always assuming that raw_map is at least double
    //the size on each dimension
    sImage32 *BB_CropTemplate(sImage32 *cost_template,
			      int *S,
			      sImage32 *label,
			      int lb){
      sImage32 *cropped;
      int x_center = 0, y_center = 0, tx_center, ty_center;
      int start_x, start_y, x, y, Imin;
      int i, p,q, ns = 0;
      //Calcua centro de massa das sementes internas:
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(label->data[p] == lb){
	  x_center += p%label->ncols;
	  y_center += p/label->ncols;
	  ns++;
	}
      }
      if(ns != 0){
	x_center /= ns;
	y_center /= ns;
      }
      //------------------------
      cropped = Image32::Create(label->ncols, label->nrows);
      tx_center = cost_template->ncols/2;
      ty_center = cost_template->nrows/2;
      if(cost_template->array[ty_center][tx_center] != 0){
	Imin = glip::Image32::GetMinVal(cost_template, &p);
	tx_center = p%cost_template->ncols;
	ty_center = p/cost_template->ncols;
      }
      start_x = tx_center - x_center;
      start_y = ty_center - y_center;
      p = 0;
      for (y = start_y; y < start_y + label->nrows; y++) {  
	for (x = start_x; x < start_x + label->ncols; x++) {      
	  if(x >= 0 && x < cost_template->ncols &&
	     y >= 0 && y < cost_template->nrows){
	    q = y * cost_template->ncols + x;
	    cropped->data[p] = cost_template->data[q];
	    p++;
	  }
	  else
	    cropped->data[p] = INT_MAX;
	}
      }
      return cropped;
    }

    /*
    sImage32 *BB_OIFT_View(sImageGraph *G,
                           sImage32 *L,
			   sImage32 *C,
			   float delta,
			   int band_id) {
      sAdjRel *A = G->A;
      int p, q, i, n, Cmin;
      Pixel u, v;
      sImage32 *B;

      n = G->ncols * G->nrows;
      B = Image32::Create(G->ncols, G->nrows);
      
      Cmin = INT_MAX;
      for (p = 0; p < n; p++) {
	if (L->data[p] == 0) continue;
	u.x = p % G->ncols;  
	u.y = p / G->ncols;
	for (i = 1; i < A->n; i++) {
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  if (Image32::IsValidPixel(L, v.x, v.y)) {
	    q = v.x + G->ncols * v.y;
	    if (L->data[q] == band_id) 
	      if (C->data[p] < Cmin) Cmin = C->data[p];
	  }
	}
      }
      
      for (p = 0; p < n; p++) {
	if (C->data[p] <= Cmin + delta && C->data[p] >= Cmin)
	  B->data[p] = 1;
      }
      return B;
    }
    */
    
    void BB_OIFT_Propagate(int s,
			   sImageGraph *G,
			   sPQueue32 *Q,
			   sImage32 *V,
			   sImage32 *L,
			   int *i_inv) {
      Pixel u, v;
      sAdjRel *A = G->A;
      int j, i, t, tmp;
      
      u.x = s % L->ncols;  
      u.y = s / L->ncols;  
      
      for (i = 1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	//if (Image32::IsValidPixel(L, v.x, v.y)) {
	if(v.x >= 0 && v.x < L->ncols &&
	   v.y >= 0 && v.y < L->nrows){
	  t = v.x + L->ncols * v.y;
	  if (Q->L.elem[t].color != BLACK) {

	    if(L->data[s]==0){
	      j = i_inv[i];
	      tmp = (G->n_link[t])[j];
	    }
	    else
	      tmp = (G->n_link[s])[i];
	    
	    if (tmp < V->data[t]) {
	      if (Q->L.elem[t].color == GRAY)
		PQueue32::FastRemoveElem(Q, t);
	      V->data[t] = tmp;
	      L->data[t] = L->data[s];
	      PQueue32::FastInsertElem(Q, t);
	    }
	  }
	}
      }
    }

    
    void BB_OIFT_Propagate_bkg(int s,
			       sImageGraph *G,
			       sPQueue32 *Q,
			       sPQueue32 *Qe,
			       sImage32 *V,
			       sImage32 *L,
			       int *i_inv) {
      sAdjRel *A = G->A;
      Pixel u, v;
      int i, t;
      
      BB_OIFT_Propagate(s, G, Q, V, L, i_inv);

      u.x = s % L->ncols;  
      u.y = s / L->ncols;  
      for (i = 1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	//if (Image32::IsValidPixel(L, v.x, v.y)) {
	if(v.x >= 0 && v.x < L->ncols &&
	   v.y >= 0 && v.y < L->nrows){
	  t = v.x + L->ncols * v.y;
	  if (Qe->L.elem[t].color == WHITE) {
	    if ((L->data[t] != 0) || (Q->L.elem[t].color != BLACK))
	      PQueue32::FastInsertElem(Qe, t);
	  }
	}
      }
    }

    
    void BB_OIFT_Propagate_obj(int s,
			       sImageGraph *G,
			       sPQueue32 *Q,
			       sPQueue32 *Qi,
			       sImage32 *V,
			       sImage32 *L,
			       int *i_inv) {
      sAdjRel *A = G->A;
      Pixel u, v;
      int i, t;
      int white, black;
      white = black = 0;
      
      u.x = s % L->ncols;
      u.y = s / L->ncols;
      for (i = 1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	//if (Image32::IsValidPixel(L, v.x, v.y)) {
	if(v.x >= 0 && v.x < L->ncols &&
	   v.y >= 0 && v.y < L->nrows){
	  t = v.x + L->ncols * v.y;
	  if (Qi->L.elem[t].color == WHITE) {
	    if (Q->L.elem[t].color == GRAY) { // mais um daqueles erros que não acontece nunca...
	      L->data[t] = 1;
	      PQueue32::FastRemoveElem(Q, t);
	      PQueue32::FastInsertElem(Qi, t);
	      BB_OIFT_Propagate(t, G, Q, V, L, i_inv);
	    }
	    else {
	      //printf("Erro 2 - %d not GRAY on queue Q - Color = %d\n", t, Q->L.elem[t].color);
	      if (Q->L.elem[t].color == BLACK) black++;
	      else if (Q->L.elem[t].color == WHITE) white++;
	      else printf("Error!\n");
	    }
	  }
	}
      }
    }

    
    sImage32 *BB_OIFT_GetLeafNodes(sImage32 *pred,
				   sAdjRel *A) {
      int n, p, q, i, isleaf;
      Pixel u, v;
      sImage32 *bin;
      
      bin = Image32::Create(pred->ncols, pred->nrows);
      n = pred->ncols * pred->nrows;
      for (p = 0; p < n; p++) {
	isleaf = 1;
	u.x = p % pred->ncols;  
	u.y = p / pred->ncols;  
	for (i = 1; i < A->n; i++) {
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if (Image32::IsValidPixel(pred, v.x, v.y)) {
	  if(v.x >= 0 && v.x < pred->ncols &&
	     v.y >= 0 && v.y < pred->nrows){
	    q = v.x + pred->ncols*v.y;
	    if(pred->data[q] == p) {
	      isleaf = 0;
	      break;
	    }
	  }
	}
	if (isleaf) bin->data[p] = 1;
	else        bin->data[p] = 0;
      }
      return bin;
    }


    
    void BB_OIFT(sImageGraph *G,
		 int *S,
		 sImage32 *L,
		 sImage32 *C,
		 sImage32 *P,
		 float delta) {
      sPQueue32 *Q  = NULL; // Fila "principal"
      sPQueue32 *Qi = NULL; // Candidatos a borda do objeto (interna - OBJ)
      sPQueue32 *Qe = NULL; // Candidatos a borda do fundo  (externa - BKG)
      sPQueue32 *Qt = NULL; // Leaves
      sImage32 *V, *Leaf, *B;
      int i, j, s, t, n, Cmin, Cmax, p, q, change_Qt;
      sSet *seed = NULL;
      sAdjRel *A = G->A;
      float tmp;
      Pixel u, v;
      int *i_inv = NULL;

      i_inv = glip::AdjRel::InverseIndexes(A);
      
      n = G->ncols * G->nrows;
      V = Image32::Create(G->ncols, G->nrows);
      
      Cmax = Image32::GetMaxVal(C);
     
      Qt = PQueue32::Create(Cmax + 2, n, C->data);
      if(P != NULL){
	Leaf = BB_OIFT_GetLeafNodes(P, G->A);
	for (t = 0; t < n; t++)
	  if (Leaf->data[t])
	    PQueue32::FastInsertElem(Qt, t);
	Image32::Destroy(&Leaf);
      }
      
      Q  = PQueue32::Create(G->Wmax + 2, n, V->data);
      Qi = PQueue32::Create(Cmax + 2, n, C->data);
      Qe = PQueue32::Create(Cmax + 2, n, C->data);
      
      for (t = 0; t < n; t++)
	V->data[t] = INT_MAX;

      Cmin = INT_MAX;
      for(i = 1; i <= S[0]; i++){
	t = S[i];
	if(L->data[t] == 1){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
	else if(L->data[t] == 0){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	  PQueue32::FastInsertElem(Qe, t);
	  if (C->data[t] < Cmin)
	    Cmin = C->data[t];	  
	}
      }
      
      for (t = 0; t < n; t++) {
	if ((C->data[t] - Cmin >= delta) && (V->data[t] != 0)) {
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	  PQueue32::FastInsertElem(Qe, t);
	}
      }
      
      while (!PQueue32::IsEmpty(Q)) {
	s = PQueue32::FastRemoveMinFIFO(Q);
	if (L->data[s] == 0) { 
	  if (Qe->L.elem[s].color == GRAY)
	    PQueue32::FastRemoveElem(Qe, s);
	  BB_OIFT_Propagate_bkg(s, G, Q, Qe, V, L, i_inv);
	  
	  change_Qt = true;
	  while(change_Qt == true) {
	    change_Qt = false;
	    while (!PQueue32::IsEmpty(Qe) &&
		   (PQueue32::FastGetMaxVal(Qe) - PQueue32::FastGetMinVal(Qe) >=
		    delta)) {
	      t = PQueue32::FastRemoveMaxFIFO(Qe);
	      if (Q->L.elem[t].color == GRAY) { 
		L->data[t] = 0;
		PQueue32::FastRemoveElem(Q, t);
		BB_OIFT_Propagate_bkg(t, G, Q, Qe, V, L, i_inv);
	      }
	      else {
		printf("Erro 1: C(t): %d\n", C->data[t]);
		/*
		printf("Bmin: %d, Bmax: %d\n", Qe->C.minvalue, Qe->C.maxvalue);
		printf("Omin: %d, Omax: %d\n", Qi->C.minvalue, Qi->C.maxvalue);
		printf("Ocolor: %d\n", Qi->L.elem[t].color);
		printf("q: %d, (x,y)=(%d,%d), Qcolor: %d, ", t, t % L->ncols, t/L->ncols, Q->L.elem[t].color);
		printf("L: %d\n", L->data[t]);
		*/
	      }
	    }

	    while (!PQueue32::IsEmpty(Qt) &&
		   (PQueue32::FastGetMaxVal(Qt) - PQueue32::FastGetMinVal(Qe) >=
		    delta)) {
	      t = PQueue32::FastRemoveMaxFIFO(Qt);
	      change_Qt = true;
	      if (Q->L.elem[t].color != BLACK) {
		L->data[t] = 0;
		if(Q->L.elem[t].color == GRAY)  PQueue32::FastRemoveElem(Q, t);
		if(Qe->L.elem[t].color == GRAY) PQueue32::FastRemoveElem(Qe, t);
		BB_OIFT_Propagate_bkg(t, G, Q, Qe, V, L, i_inv);
	      }
	      else if (L->data[t] == 1) printf("Erro 3:\n");
	    }
	  }
	}
	else if (L->data[s] == 1) {
	  PQueue32::FastInsertElem(Qi, s);
	  BB_OIFT_Propagate(s, G, Q, V, L, i_inv);
	  while (PQueue32::FastGetMaxVal(Qi) - PQueue32::FastGetMinVal(Qi) >=
		 delta) {
	    t = PQueue32::FastRemoveMinFIFO(Qi);
	    BB_OIFT_Propagate_obj(t, G, Q, Qi, V, L, i_inv);
	  }
	}
      }
      Image32::Destroy(&V);
      PQueue32::Destroy(&Qt);
      PQueue32::Destroy(&Qi);
      PQueue32::Destroy(&Qe);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }

    
    
    void RBB_OIFT(sImageGraph *G,
		  int *S,
		  sImage32 *L,
		  sImage32 *C,
		  sImage32 *P,
		  float delta) {
      sPQueue32 *Q  = NULL; // Fila "principal"
      sPQueue32 *Qi = NULL; // Candidatos a borda do objeto (interna - OBJ)
      sPQueue32 *Qe = NULL; // Candidatos a borda do fundo  (externa - BKG)
      sPQueue32 *Qt = NULL; // Leaves
      sImage32 *V, *Leaf, *B;
      int i, j, s, t, n, Cmin, Cmax, p, q, change_Qt;
      sSet *seed = NULL;
      sAdjRel *A = G->A;
      float tmp;
      Pixel u, v;
      int *i_inv = NULL;

      i_inv = glip::AdjRel::InverseIndexes(A);
      
      n = G->ncols * G->nrows;
      V = Image32::Create(G->ncols, G->nrows);
      
      Cmax = Image32::GetMaxVal(C);
     
      Qt = PQueue32::Create(Cmax + 2, n, C->data);
      if(P != NULL){
	Leaf = BB_OIFT_GetLeafNodes(P, G->A);
	for (t = 0; t < n; t++)
	  if (Leaf->data[t])
	    PQueue32::FastInsertElem(Qt, t);
	Image32::Destroy(&Leaf);
      }
      
      Q  = PQueue32::Create(G->Wmax + 2, n, V->data);
      Qi = PQueue32::Create(Cmax + 2, n, C->data);
      Qe = PQueue32::Create(Cmax + 2, n, C->data);
      
      for (t = 0; t < n; t++)
	V->data[t] = INT_MAX;

      Cmin = INT_MAX;
      for(i = 1; i <= S[0]; i++){
	t = S[i];
	if(L->data[t] == 1){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
	else if(L->data[t] == 0){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	  PQueue32::FastInsertElem(Qe, t);
	  if (C->data[t] < Cmin)
	    Cmin = C->data[t];	  
	}
      }
      
      for (t = 0; t < n; t++) {
	if ((C->data[t] - Cmin >= delta * Cmin) && (V->data[t] != 0)) {
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	  PQueue32::FastInsertElem(Qe, t);
	}
      }
      
      while (!PQueue32::IsEmpty(Q)) {
	s = PQueue32::FastRemoveMinFIFO(Q);
	if (L->data[s] == 0) { 
	  if (Qe->L.elem[s].color == GRAY)
	    PQueue32::FastRemoveElem(Qe, s);
	  BB_OIFT_Propagate_bkg(s, G, Q, Qe, V, L, i_inv);
	  
	  change_Qt = true;
	  while(change_Qt == true) {
	    change_Qt = false;
	    while (!PQueue32::IsEmpty(Qe) &&
		   (PQueue32::FastGetMaxVal(Qe) - PQueue32::FastGetMinVal(Qe) >=
		    delta * PQueue32::FastGetMinVal(Qe))) {
	      t = PQueue32::FastRemoveMaxFIFO(Qe);
	      if (Q->L.elem[t].color == GRAY) { 
		L->data[t] = 0;
		PQueue32::FastRemoveElem(Q, t);
		BB_OIFT_Propagate_bkg(t, G, Q, Qe, V, L, i_inv);
	      }
	      else {
		printf("Erro 1: C(t): %d\n", C->data[t]);
		/*
		printf("Bmin: %d, Bmax: %d\n", Qe->C.minvalue, Qe->C.maxvalue);
		printf("Omin: %d, Omax: %d\n", Qi->C.minvalue, Qi->C.maxvalue);
		printf("Ocolor: %d\n", Qi->L.elem[t].color);
		printf("q: %d, (x,y)=(%d,%d), Qcolor: %d, ", t, t % L->ncols, t/L->ncols, Q->L.elem[t].color);
		printf("L: %d\n", L->data[t]);
		*/
	      }
	    }
	    
	    while (!PQueue32::IsEmpty(Qt) &&
		   (PQueue32::FastGetMaxVal(Qt) - PQueue32::FastGetMinVal(Qe) >=
		    delta * PQueue32::FastGetMinVal(Qe))) {
	      t = PQueue32::FastRemoveMaxFIFO(Qt);
	      change_Qt = true;
	      if (Q->L.elem[t].color != BLACK) {
		L->data[t] = 0;
		if(Q->L.elem[t].color == GRAY)  PQueue32::FastRemoveElem(Q, t);
		if(Qe->L.elem[t].color == GRAY) PQueue32::FastRemoveElem(Qe, t);
		BB_OIFT_Propagate_bkg(t, G, Q, Qe, V, L, i_inv);
	      }
	      else if (L->data[t] == 1) printf("Erro 3:\n");
	    }
	  }
	}
	else if (L->data[s] == 1) {
	  PQueue32::FastInsertElem(Qi, s);
	  BB_OIFT_Propagate(s, G, Q, V, L, i_inv);
	  while (PQueue32::FastGetMaxVal(Qi) - PQueue32::FastGetMinVal(Qi) >=
		 delta * PQueue32::FastGetMinVal(Qi)) {
	    t = PQueue32::FastRemoveMinFIFO(Qi);
	    BB_OIFT_Propagate_obj(t, G, Q, Qi, V, L, i_inv);
	  }
	}
      }
      Image32::Destroy(&V);
      PQueue32::Destroy(&Qt);
      PQueue32::Destroy(&Qi);
      PQueue32::Destroy(&Qe);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }

    //----------------------------------------

    void B_OIFT(sImageGraph *G,
		int *S,
		sImage32 *L,
		sImage32 *C,
		float delta) {
      sPQueue32 *Q  = NULL; // Main queue.
      sPQueue32 *Qo = NULL; // Object pixels.
      sPQueue32 *Qb = NULL; // Background pixels.
      sPQueue32 *Qx = NULL; // Undefined pixels.
      sImage32 *V;
      int i, j, s, t, n, Cmin, Cmax, p, q;
      sSet *seed = NULL;
      sAdjRel *A = G->A;
      float tmp;
      Pixel u, v;
      int *i_inv = NULL;

      i_inv = glip::AdjRel::InverseIndexes(A);
      
      n = G->ncols * G->nrows;
      V = Image32::Create(G->ncols, G->nrows);
      
      Cmax = Image32::GetMaxVal(C);
     
      Q  = PQueue32::Create(G->Wmax + 2, n, V->data);
      Qo = PQueue32::Create(Cmax + 2, n, C->data);
      Qb = PQueue32::Create(Cmax + 2, n, C->data);
      Qx = PQueue32::Create(Cmax + 2, n, C->data);
      
      for (t = 0; t < n; t++){
	V->data[t] = INT_MAX;
	PQueue32::FastInsertElem(Qx, t);
      }

      for(i = 1; i <= S[0]; i++){
	t = S[i];
	if(L->data[t] == 1){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
	else if(L->data[t] == 0){
	  V->data[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
      }
      
      while (!PQueue32::IsEmpty(Q)) {
	s = PQueue32::FastRemoveMinFIFO(Q);
	BB_OIFT_Propagate(s, G, Q, V, L, i_inv);
	PQueue32::FastRemoveElem(Qx, s);
	
	if(L->data[s] == 0){ 
	  PQueue32::FastInsertElem(Qb, s);

	  while (!PQueue32::IsEmpty(Qx) &&
		 (PQueue32::FastGetMaxVal(Qx) - PQueue32::FastGetMinVal(Qb) >= delta)){
	    t = PQueue32::FastRemoveMaxFIFO(Qx);
	    PQueue32::FastInsertElem(Qb, t);
	    L->data[t] = 0;
	    if (Q->L.elem[t].color == GRAY)
	      PQueue32::FastRemoveElem(Q, t);
	    else
	      Q->L.elem[t].color = BLACK;
	    BB_OIFT_Propagate(t, G, Q, V, L, i_inv);
	  }
	}
	else if(L->data[s] == 1){
	  PQueue32::FastInsertElem(Qo, s);
	  
	  while (!PQueue32::IsEmpty(Qx) &&
		 (PQueue32::FastGetMaxVal(Qo) - PQueue32::FastGetMinVal(Qx) >= delta)){
	    t = PQueue32::FastRemoveMinFIFO(Qx);
	    PQueue32::FastInsertElem(Qo, t);
	    L->data[t] = 1;
	    if (Q->L.elem[t].color == GRAY)
	      PQueue32::FastRemoveElem(Q, t);
	    else
	      Q->L.elem[t].color = BLACK;
	    BB_OIFT_Propagate(t, G, Q, V, L, i_inv);
	  }
	}
      }
      Image32::Destroy(&V);
      PQueue32::Destroy(&Qo);
      PQueue32::Destroy(&Qb);
      PQueue32::Destroy(&Qx);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }



    void B_OIFT_Propagate(int s,
			  sGraph *graph,
			  sGraph *transpose,
			  sPQueue32 *Q,
			  int *V,
			  int *L){
      int i, t, tmp;
      sGraph *g;

      if(L[s]==0) g = transpose;
      else   	  g = graph;
      
      for(i = 0; i < g->nodes[s].outdegree; i++){
	t = g->nodes[s].adjList[i];

	if (Q->L.elem[t].color != BLACK) {
	  tmp = g->nodes[s].Warcs[i];
  
	  if (tmp < V[t]) {
	    if (Q->L.elem[t].color == GRAY)
	      PQueue32::FastRemoveElem(Q, t);
	    V[t] = tmp;
	    L[t] = L[s];
	    PQueue32::FastInsertElem(Q, t);
	  }
	}
      }
    }
    


    void B_OIFT(sGraph *graph,
		sGraph *transpose,
		int *S,
		int *L,
		int *C,
		float delta) {
      sPQueue32 *Q  = NULL; // Main queue.
      sPQueue32 *Qo = NULL; // Object pixels.
      sPQueue32 *Qb = NULL; // Background pixels.
      sPQueue32 *Qx = NULL; // Undefined pixels.
      int *V;
      int i, j, s, t, n, Cmin, Cmax, p, q, Wmax;
      sSet *seed = NULL;
      float tmp;
      
      n = graph->nnodes;
      V = glip::AllocIntArray(n);

      Cmax = INT_MIN;
      for(p=0; p<n; p++)
	if(C[p] > Cmax)
	  Cmax = C[p];

      Wmax = glip::Graph::GetMaximumArc(graph);
      
      Q  = PQueue32::Create(Wmax + 2, n, V);
      Qo = PQueue32::Create(Cmax + 2, n, C);
      Qb = PQueue32::Create(Cmax + 2, n, C);
      Qx = PQueue32::Create(Cmax + 2, n, C);
      
      for (t = 0; t < n; t++){
	V[t] = INT_MAX;
	PQueue32::FastInsertElem(Qx, t);
      }

      for(i = 1; i <= S[0]; i++){
	t = S[i];
	if(L[t] == 1){
	  V[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
	else if(L[t] == 0){
	  V[t] = 0;
	  PQueue32::FastInsertElem(Q, t);
	}
      }
      
      while (!PQueue32::IsEmpty(Q)) {
	s = PQueue32::FastRemoveMinFIFO(Q);
	B_OIFT_Propagate(s, graph, transpose, Q, V, L);
	PQueue32::FastRemoveElem(Qx, s);
	
	if(L[s] == 0){ 
	  PQueue32::FastInsertElem(Qb, s);

	  while (!PQueue32::IsEmpty(Qx) &&
		 (PQueue32::FastGetMaxVal(Qx) - PQueue32::FastGetMinVal(Qb) >= delta)){
	    t = PQueue32::FastRemoveMaxFIFO(Qx);
	    PQueue32::FastInsertElem(Qb, t);
	    L[t] = 0;
	    if (Q->L.elem[t].color == GRAY)
	      PQueue32::FastRemoveElem(Q, t);
	    else
	      Q->L.elem[t].color = BLACK;
	    B_OIFT_Propagate(t, graph, transpose, Q, V, L);
	  }
	}
	else if(L[s] == 1){
	  PQueue32::FastInsertElem(Qo, s);
	  
	  while (!PQueue32::IsEmpty(Qx) &&
		 (PQueue32::FastGetMaxVal(Qo) - PQueue32::FastGetMinVal(Qx) >= delta)){
	    t = PQueue32::FastRemoveMinFIFO(Qx);
	    PQueue32::FastInsertElem(Qo, t);
	    L[t] = 1;
	    if (Q->L.elem[t].color == GRAY)
	      PQueue32::FastRemoveElem(Q, t);
	    else
	      Q->L.elem[t].color = BLACK;
	    B_OIFT_Propagate(t, graph, transpose, Q, V, L);
	  }
	}
      }
      glip::FreeIntArray(&V);
      PQueue32::Destroy(&Qo);
      PQueue32::Destroy(&Qb);
      PQueue32::Destroy(&Qx);
      PQueue32::Destroy(&Q);
    }

    
    
    //----------------------------------------

    int *GetSeedsByLabel(int *S,
			 sImage32 *label,
			 int lb){
      int *Slb = NULL;
      int l, s, i, j, size = 0;
      for(i = 1; i <= S[0]; i++){
	s = S[i];
	l = label->data[s];
	if(l == lb) size++;
      }
      Slb = glip::AllocIntArray(size+1);
      Slb[0] = size;
      i = 1;
      for(j = 1; j <= S[0]; j++){
	s = S[j];
	l = label->data[s];
	if(l == lb){
	  Slb[i] = s;
	  i++;
	}
      }
      return Slb;
    }
    
    
    int *GetAllInternalSeedsByLabel(int *S,
				    sImage32 *label,
				    int lb,
				    int *hr,
				    int nlayers){
      int *Sall, *Sl, *Stmp;
      int l,size,i,j;
      Sall = GetSeedsByLabel(S, label, lb);
      for(l=0; l < nlayers; l++){
	if(hr[l] == lb-1){
	  //Sl = GetSeedsByLabel(S, label, l+1);
	  Sl = GetAllInternalSeedsByLabel(S, label, l+1,
					  hr, nlayers);
	  size = Sl[0] + Sall[0];
	  Stmp = glip::AllocIntArray(size+1);
	  Stmp[0] = size;
	  j =  1;
	  for(i = 1; i <= Sall[0]; i++){
	    Stmp[j] = Sall[i];
	    j++;
	  }
	  for(i = 1; i <= Sl[0]; i++){
	    Stmp[j] = Sl[i];
	    j++;
	  }
	  glip::FreeIntArray(&Sall);
	  glip::FreeIntArray(&Sl);
	  Sall = Stmp;
	}
      }
      return Sall;
    }
    

    sImageGraph *GetPolarityGraph(sImageGraph *graph,
				  sCImage *cimg,
				  sImage32 *img,
				  int pol){
      sImageGraph *sg_oriented;
      sg_oriented = glip::ImageGraph::Clone(graph);
      if(cimg != NULL){
	glip::sImage32 *lumi;
	lumi = glip::Image32::Luminosity(cimg);
	glip::ImageGraph::Orient2Digraph(sg_oriented, lumi, pol);
	glip::Image32::Destroy(&lumi);
      }
      else
	glip::ImageGraph::Orient2Digraph(sg_oriented, img, pol);
      return sg_oriented;
    }
    

    
    void HL_OIFT(sImageGraph *graph,
		 sCImage *cimg,
		 sImage32 *img,
		 float radius,
		 char *hierarchy,
		 int *S,
		 sImage32 *label){
      int p;
      //-----------------------
      glip::sLayeredGraph *lg;
      glip::sImageGraph *sg, *sg_oriented;
      //------------------------
      FILE *fpHierarchy;
      int nlayers, typeL, nRelations, relationT, layer1, layer2;
      int i, j, n, nSeeds, k, s, Wmax, pol;
      int x, y, lb, ncols, nrows;
      int *lg_label;
      
      ncols = label->ncols;
      nrows = label->nrows;
      glip::ImageGraph::ChangeType(graph, DISSIMILARITY);
      
      fpHierarchy = fopen((char *)hierarchy, "r");
      if (!fpHierarchy){
	glip::Error((char *)MSG1,(char *)"Error opening file of Hierarchy Relations!");
      }
      
      /*------------- Set number of layers -------------*/
      (void)fscanf(fpHierarchy,"%d\n",&nlayers);
      //printf("NLayers =  %d\n",nlayers);
      
      /*-------------GIVEN BY USER: TYPE OF LAYER -------------*/
      /*typeLayer:
	0 = normal,
	1 = GSC,
      */
      int *typeLayer = glip::AllocIntArray(nlayers);
      int *polLayer = glip::AllocIntArray(nlayers);
      
      for(i = 0; i < nlayers; i++) {
	(void)fscanf(fpHierarchy,"%d %d\n",&typeL, &pol);
	typeLayer[i] = typeL;
	polLayer[i] = pol;
      }
      
      /*------------- Creating HIERARCHY -------------*/
      /*Relations is given by inclusion(=1) or exclusion(=2) binary relation between: (son,father) and (brothers) that defines a hierarchy, in a .txt*/
      /*Example: 1|0|1 = inclusion of layer 0 in layer 1 
	or 2|0|1 = exclusion between layer 0 and layer 1*/
      
      /*Initialization of Hierarchy*/
      int *hr = glip::AllocIntArray(nlayers);
      for(i = 0; i < nlayers; i++){
	hr[i]= -1; 
      }
      
      /*Reading a .txt fill, and set in "hr" with only inclusion cases*/
      // hr[0] = 1;  // i(0,1) = 1 0 1  // means: "1 is father of 0" or "0 is included in 1"
      (void)fscanf(fpHierarchy,"%d\n",&nRelations);
      for(i = 0; i < nRelations; i++) {
	(void)fscanf(fpHierarchy,"%d %d %d\n",&relationT,&layer1,&layer2);
	if(relationT == 1){ // it's inclusion relation
	  hr[layer1] = layer2;
	}
      }
      
      fclose(fpHierarchy);
      
      /* ------------- Create label image -------------*/
      n = ncols*nrows*nlayers;
      lg_label = glip::AllocIntArray(n);
      /*Set NIL (-1) in all matrix positions*/
      for(i = 0; i < n; i++){
	lg_label[i] = NIL;
      }
      
      nSeeds = S[0];
      int* Seeds = glip::AllocIntArray(nSeeds*nlayers + 1);
      
      k = 1;
      for(i = 1; i <= S[0]; i++){
	s = S[i];
	x = s%ncols;
	y = s/ncols;
	lb = label->data[s];
	
	if(lb != 0){
	  p = x + y*ncols + nrows*ncols*(lb-1); 
	  lg_label[p] = 1;
	  Seeds[k] = p;
	  k++;
	  /*Essa parte faltou no algoritmo do artigo:*/
	  for(j = 0; j < nlayers ; j++){
	    /*
	    if(j != lb-1 && hr[lb-1] == j){
	      p = x + y*ncols + nrows*ncols*(j);
	      lg_label[p] = 1;
	      Seeds[k] = p;
	      k++;
	    }
	    else if(j != lb-1){
	      p = x + y*ncols + nrows*ncols*(j);
	      lg_label[p] = 0;
	      Seeds[k] = p;
	      k++;
	    }
	    */
	    if(j != lb-1 && hr[j] == lb-1){
	      p = x + y*ncols + nrows*ncols*(j);
	      lg_label[p] = 0;
	      Seeds[k] = p;
	      k++;
	    }
	  }
	}
	else{ // it is background seed, then transfer for all layers
	  for(j = 0; j < nlayers ; j++){
	    p = x + y*ncols + nrows*ncols*(j);
	    lg_label[p] = 0;
	    Seeds[k] = p;
	    k++;
	  }
	}
      }
      Seeds[0] = k-1;

      /* ------------- Create LAYERED GRAPH -------------*/
      lg = glip::LayeredGraph::Create(nlayers, ncols*nrows);
      
      /*Create a graph image with 8-neighborhood*/
      sg = graph;
      //glip::SparseGraph::SuppressZeroWeightedArcs(sg);
      Wmax = sg->Wmax;

      /*Set cost of arcs for each layer*/
      for(i = 0; i < nlayers; i++){
	if(typeLayer[i] == 0){ /* 0 = normal*/
	  if(polLayer[i] == 0)
	    glip::LayeredGraph::SetArcs(lg, sg, i);
	  else{
	    sg_oriented = GetPolarityGraph(sg, cimg, img, polLayer[i]);
	    if(sg_oriented->Wmax > Wmax) Wmax = sg_oriented->Wmax;
	    glip::LayeredGraph::SetArcs(lg, sg_oriented, i);
	    glip::ImageGraph::Destroy(&sg_oriented);
	  }
	}
	else if(typeLayer[i] == 1){ /* 1 = GSC */
	  glip::sImageGraph *sg_GSC;
	  glip::sImage32 *P_sum;
	  int *Slb = NULL;
	  if(polLayer[i] == 0)
	    sg_GSC = glip::ImageGraph::Clone(sg);
	  else{
	    sg_GSC = GetPolarityGraph(sg, cimg, img, polLayer[i]);
	    if(sg_GSC->Wmax > Wmax) Wmax = sg_GSC->Wmax;
	  }
	  Slb = GetAllInternalSeedsByLabel(S, label, i+1, hr, nlayers);
	  //glip::SparseGraph::SuppressZeroWeightedArcs(sg_GSC);
	  P_sum = glip::ift::SC_Pred_fsum(sg_GSC, Slb, 0.1);
	  glip::ImageGraph::Orient2DigraphOuter(sg_GSC, P_sum);
	  glip::LayeredGraph::SetArcs(lg, sg_GSC, i);  
	  glip::Image32::Destroy(&P_sum);
	  glip::ImageGraph::Destroy(&sg_GSC);
	  glip::FreeIntArray(&Slb);
	}
      }

      
      for(i = 0; i < nlayers ; i++){
	if(hr[i] != -1){ //son and father
	  glip::LayeredGraph::SetArcs(lg, i, hr[i], ncols, 0.0, radius);
	  glip::LayeredGraph::SetArcs(lg, hr[i], i, ncols, (float)(Wmax+1), radius);
	}
      }
      

      //Second exclusion
      for(i = 0; i < nlayers-1; i++){
	for(j = i+1; j < nlayers; j++){
	  if( hr[i] == hr[j] ){ //same father -> brothers -> exclusion
	    glip::LayeredGraph::SetArcs(lg, i, j, ncols, (float)(Wmax+1), radius);
	    glip::LayeredGraph::SetArcs(lg, j, i, ncols, 0.0, radius);
	  }
	}
      }
      

      /* ------------- Executa a Hierarchical Layered OIFT -------------*/ 
      glip::ift::HL_OIFT(lg, Wmax, Seeds, lg_label, hr);

      
      /* ------------- OUTPUT RESULT-------------*/
      
      glip::sQueue *FIFO;
      FIFO = glip::Queue::Create(nlayers);
      int *depth = glip::AllocIntArray(nlayers);
      int dmax = 0;
      for(i = 0; i < nlayers ; i++){
	if(hr[i] == -1){
	  depth[i] = 0;
	  glip::Queue::Push(FIFO, i);
	}
	else
	  depth[i] = -1;
      }
      while(!glip::Queue::IsEmpty(FIFO)){
	j = glip::Queue::Pop(FIFO);
	for(i = 0; i < nlayers ; i++){
	  if(hr[i] == j){
	    depth[i] = depth[j] + 1;
	    if(depth[i] > dmax) dmax = depth[i];
	    glip::Queue::Push(FIFO, i);
	  }
	}
      }
      
      int sizeImg = ncols*nrows;
      
      glip::Image32::Set(label, 0);
      int d;
      for(d = 0; d <= dmax; d++){
	for(i = 0; i < nlayers; i++){
	  if(depth[i] == d){
	    for(j = 0; j < sizeImg; j++){
	      if(lg_label[j+ sizeImg*i] != 0)
		label->data[j] = i+1;
	    }
	  }
	}
      }

      glip::LayeredGraph::Destroy(&lg);
      
      glip::FreeIntArray(&typeLayer);
      glip::FreeIntArray(&polLayer);
      glip::FreeIntArray(&hr);
      glip::FreeIntArray(&lg_label);
      glip::FreeIntArray(&Seeds);
      glip::FreeIntArray(&depth);
      glip::Queue::Destroy(&FIFO);
    }


    //----------------------------------------

    
    void HL_OIFT(sLayeredGraph *lg,
		 int Wmax, int *S, int *L, int *hr){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,layer_p,layer_q,size_layer,tmp2=0;
      int w_pq=0, w_qp=0, tmp=0;
      int* value;
      sGraphNode *A;
      int exclusionCase = 0;
      int *i_inv = NULL;
      int i_inv_size;
      bool flag;

      A = lg->graph->nodes;
      p = 0;
      i_inv_size = A[p].outdegree;
      i_inv = glip::AllocIntArray(i_inv_size);
      flag = true;
      while(flag){
	flag = false;
	for(i = 0; i < A[p].outdegree && i < i_inv_size; i++){
	  q = A[p].adjList[i];
	  if(p == q){ flag = true; break; }
	  for(j = 0; j < A[q].outdegree; j++)
	    if(A[q].adjList[j] == p)
	      i_inv[i] = j;
	}
	p++;
      }
      
      size_layer = lg->nnodesperlayer;
      n = size_layer*lg->nlayers;
      
      //printf("Size Layer: %d\n",size_layer);
      //printf("N pixels in the Graph: %d\n",n);
      
      /*Initialization*/
      value = glip::AllocIntArray(n);
      //Q = glip::PQueue32::Create(Wmax*lg->nlayers+2,n,value);
      Q = glip::PQueue32::Create(Wmax+3,n,value);
      
      /*Insert in value*/
      for(p=0; p<n; p++){
	if(L[p]==NIL) value[p] = INT_MAX; //(int)floor(FLT_MAX+0.5);
	else          value[p] = 0;
      }
      
      /*Insert in PQueue*/
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  glip::PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(L[p]!=NIL)
	    glip::PQueue32::FastInsertElem(Q, p);	    
      }
      
      /*
	clock_t t;
	t = clock();
      */
      
      /*Starting OIFT*/
      while(!glip::PQueue32::IsEmpty(Q)){
	p = glip::PQueue32::FastRemoveMinFIFO(Q);
	
	for(i=0; i<A[p].outdegree; i++){
	  q = A[p].adjList[i];
	  //if(q != NIL){
	  if(Q->L.elem[q].color != BLACK){
	    //w_pq = ROUND(glip::Graph::GetArcWeight(lg->graph,p,q));
	    w_pq = ROUND(A[p].Warcs[i]);
	    //if(w_pq < 0) continue;
	    
	    layer_p = p/size_layer;
	    layer_q = q/size_layer;
	    exclusionCase = 0;
	    /* Analize each relation type*/
	    if (layer_p == layer_q){  /*SAME LAYER*/
	      // Get w_qp
	      //w_qp = ROUND(glip::Graph::GetArcWeight(lg->graph,q,p));
	      w_qp = ROUND(A[q].Warcs[i_inv[i]]);
	      exclusionCase = 0;
	      if(L[p] != 0)
		tmp = w_pq;
	      else 
		tmp = w_qp; 
	    }
	    /*INTER LAYERS: relations where defined*/
	    else if(layer_p != layer_q){  
	      // Get w_qp
	      if(w_pq == 0)
		w_qp = Wmax+1; //INT_MAX;
	      else
		w_qp = 0;
	      // inclusion
	      /*layer_p is "INCLUDED" in layer_q /OR/ layer_q is "INCLUDED" in layer_p */
	      if((hr[layer_p] == layer_q) || (hr[layer_q] == layer_p)){  
		exclusionCase = 0;
		if(L[p] != 0)
		  tmp = w_pq;
		else 
		  tmp = w_qp; 
	      }
	      // exclusion
	      /*layer_q is EXCLUDED from layer_p */
	      else if(hr[layer_p] == hr[layer_q]){  
		exclusionCase = 1;
		/* It was defined to change "inputs"/"outputs", in/from the higher layer for exclusion case*/
		if (layer_p < layer_q){ 
		  if(L[p] != 0)
		    tmp = w_qp;
		  else 
		    tmp = w_pq;
		}
		else{
		  if(L[p] != 0)
		    tmp = w_pq;
		  else 
		    tmp = w_qp;
		}
	      }
	    }
	    
	    if(tmp < value[q]){
	      if(Q->L.elem[q].color == GRAY)
		glip::PQueue32::FastRemoveElem(Q, q);
	      
	      value[q] = tmp;
	      
	      if(exclusionCase == 1){ /*Treat exclusion case*/
		if(L[p] != 0)
		  L[q] = 0;
	      }
	      else{ /*exclusion == 0*/
		L[q] = L[p];
	      }
	      
	      glip::PQueue32::FastInsertElem(Q, q);
	    } 
	  }
	  //}
	}
      }

      /*
	t = clock() - t;
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
	printf("HLOIFT took %d seconds(t) to execute \n", ((int)t));
	printf("HLOIFT took %f seconds(t/C_P_S) to execute \n", time_taken);      
      */
      
      glip::PQueue32::Destroy(&Q);
      glip::FreeIntArray(&value);
      glip::FreeIntArray(&i_inv);
    }


    //--------------------------------------------------------
 
    void HL_OIFT_2(sLayeredGraph *lg,
		   int Wmax, int *S, int *L, int *hr){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,layer_p,layer_q,size_layer,tmp2=0;
      int w_pq=0, w_qp=0, tmp=0;
      int* value;
      sGraphNode *A;
      int exclusionCase = 0;

      A = lg->graph->nodes;
      size_layer = lg->nnodesperlayer;
      n = size_layer*lg->nlayers;

      /*Initialization*/
      value = glip::AllocIntArray(n);
      Q = glip::PQueue32::Create(Wmax+3,n,value);
      
      /*Insert in value*/
      for(p=0; p<n; p++){
	if(L[p]==NIL) value[p] = INT_MAX; //(int)floor(FLT_MAX+0.5);
	else          value[p] = 0;
      }
      
      /*Insert in PQueue*/
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  glip::PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(L[p]!=NIL)
	    glip::PQueue32::FastInsertElem(Q, p);	    
      }
      
      /*Starting OIFT*/
      while(!glip::PQueue32::IsEmpty(Q)){
	p = glip::PQueue32::FastRemoveMinFIFO(Q);
	
	for(i=0; i<A[p].outdegree; i++){
	  q = A[p].adjList[i];
	  //if(q != NIL){
	  if(Q->L.elem[q].color != BLACK){
	    //w_pq = ROUND(glip::Graph::GetArcWeight(lg->graph,p,q));
	    w_pq = ROUND(A[p].Warcs[i]);
	    //if(w_pq < 0) continue;
	    
	    layer_p = p/size_layer;
	    layer_q = q/size_layer;
	    exclusionCase = 0;
	    /* Analize each relation type*/
	    if (layer_p == layer_q){  /*SAME LAYER*/
	      // Get w_qp
	      w_qp = ROUND(glip::Graph::GetArcWeight(lg->graph,q,p));
	      exclusionCase = 0;
	      if(L[p] != 0)
		tmp = w_pq;
	      else 
		tmp = w_qp; 
	    }
	    /*INTER LAYERS: relations where defined*/
	    else if(layer_p != layer_q){  
	      // Get w_qp
	      if(w_pq == 0)
		w_qp = Wmax+1; //INT_MAX;
	      else
		w_qp = 0;
	      // inclusion
	      /*layer_p is "INCLUDED" in layer_q /OR/ layer_q is "INCLUDED" in layer_p */
	      if((hr[layer_p] == layer_q) || (hr[layer_q] == layer_p)){  
		exclusionCase = 0;
		if(L[p] != 0)
		  tmp = w_pq;
		else 
		  tmp = w_qp; 
	      }
	      // exclusion
	      /*layer_q is EXCLUDED from layer_p */
	      else if(hr[layer_p] == hr[layer_q]){  
		exclusionCase = 1;
		/* It was defined to change "inputs"/"outputs", in/from the higher layer for exclusion case*/
		if (layer_p < layer_q){ 
		  if(L[p] != 0)
		    tmp = w_qp;
		  else 
		    tmp = w_pq;
		}
		else{
		  if(L[p] != 0)
		    tmp = w_pq;
		  else 
		    tmp = w_qp;
		}
	      }
	    }
	    
	    if(tmp < value[q]){
	      if(Q->L.elem[q].color == GRAY)
		glip::PQueue32::FastRemoveElem(Q, q);
	      
	      value[q] = tmp;
	      
	      if(exclusionCase == 1){ /*Treat exclusion case*/
		if(L[p] != 0)
		  L[q] = 0;
	      }
	      else{ /*exclusion == 0*/
		L[q] = L[p];
	      }
	      
	      glip::PQueue32::FastInsertElem(Q, q);
	    } 
	  }
	}
      }

      glip::PQueue32::Destroy(&Q);
      glip::FreeIntArray(&value);
    }
   


    //Compute watershed by fpeak from markers
    void IFT_fpeak(sImage32 *grad,
		   sAdjRel *A,
		   sImage32 *label){
      sPQueue32 *Q;
      sImage32 *cost; //*Rmin;
      int tmp, w, Wmax;
      int n,p,q,i,px,py,qx,qy;
      n = grad->n;
      Wmax = glip::Image32::GetMaxVal(grad);
      //Rmin = glip::Image32::RegMin(grad, A);
      cost = glip::Image32::Create(grad);
      Q = glip::PQueue32::Create(Wmax+2, n, cost->data);

      for(p = 0; p < n; p++){
	//if(Rmin->data[p] == 0){
	if(label->data[p] == NIL){
	  cost->data[p] = INT_MAX;
	  //label->data[p] = NIL;
	}
	else{
	  cost->data[p] = grad->data[p];
	  //label->data[p] = Rmin->data[p]-1;
	  PQueue32::FastInsertElem(Q, p);
	}
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	px = p%grad->ncols;
	py = p/grad->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(glip::Image32::IsValidPixel(grad, qx,qy)){
	  if(qx >= 0 && qx < grad->ncols &&
	     qy >= 0 && qy < grad->nrows){
	    q = qx + qy*grad->ncols;
	      
	    if(Q->L.elem[q].color != BLACK){
	      w = grad->data[q];
	      tmp = MAX(cost->data[p], w);
	      
	      if(tmp < cost->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		cost->data[q] = tmp;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      glip::PQueue32::Destroy(&Q);
      glip::Image32::Destroy(&cost);
      //glip::Image32::Destroy(&Rmin);
    }

    
    //Compute watershed by fwv from markers
    void IFT_fwv(sImage32 *grad,
		 sAdjRel *A,
		 sImage32 *label){
      sPQueue32 *Q;
      sImage32 *cost; //*Rmin;
      int tmp, Wmax;
      int n,p,q,i,px,py,qx,qy;
      n = grad->n;
      Wmax = glip::Image32::GetMaxVal(grad);
      //Rmin = glip::Image32::RegMin(grad, A);
      cost = glip::Image32::Create(grad);
      Q = glip::PQueue32::Create(Wmax+2, n, cost->data);

      for(p = 0; p < n; p++){
	//if(Rmin->data[p] == 0){
	if(label->data[p] == NIL){
	  cost->data[p] = INT_MAX;
	  //label->data[p] = NIL;
	}
	else{
	  cost->data[p] = grad->data[p];
	  //label->data[p] = Rmin->data[p]-1;
	  PQueue32::FastInsertElem(Q, p);
	}
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	px = p%grad->ncols;
	py = p/grad->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(glip::Image32::IsValidPixel(grad, qx,qy)){
	  if(qx >= 0 && qx < grad->ncols &&
	     qy >= 0 && qy < grad->nrows){
	      q = qx + qy*grad->ncols;
	      
	    if(Q->L.elem[q].color != BLACK){
	      tmp = grad->data[q];
	      
	      if(tmp < cost->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		cost->data[q] = tmp;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      glip::PQueue32::Destroy(&Q);
      glip::Image32::Destroy(&cost);
      //glip::Image32::Destroy(&Rmin);
    }




    void RelaxMobj(sImageGraph *sg,
		   int *S,
		   sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1[100],*flabel_2[100],*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic,Lmax,l,lb,lm;
      float sw,w;

      Lmax = glip::Image32::GetMaxVal(label);
      glip::ImageGraph::ChangeType(sg, DISSIMILARITY);
      
      A = sg->A;
      ninic = 1;
      n = label->n;
      for(l = 0; l <= Lmax; l++){
	flabel_1[l] = (float *)calloc(n, sizeof(float));
	flabel_2[l] = (float *)calloc(n, sizeof(float));
      }
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = glip::Image32::Create(label);
      for(p = 0; p < n; p++){
	lb = label->data[p];
	flabel_1[lb][p] = 1.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(glip::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	for(l = 0; l <= Lmax; l++)
	  memcpy(flabel_2[l], flabel_1[l], n*sizeof(float));

	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  for(l = 0; l <= Lmax; l++)
	    flabel_2[l][p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(glip::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      w = (sg->n_link[p])[i];
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      for(l = 0; l <= Lmax; l++)
		flabel_2[l][p] += w*flabel_1[l][q];
	    }
	  }
	  for(l = 0; l <= Lmax; l++)
	    flabel_2[l][p] /= sw;
	}
	for(l = 0; l <= Lmax; l++){
	  tmp = flabel_1[l];
	  flabel_1[l] = flabel_2[l];
	  flabel_2[l] = tmp;
	}
	
	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  lb = label->data[p];
	  for(l = 0; l <= Lmax; l++){
	    if(lb == l)
	      flabel_1[l][p] = 1.0;
	    else
	      flabel_1[l][p] = 0.0;
	  }
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(glip::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	lm = 0;
	for(l = 1; l <= Lmax; l++){
	  if(flabel_1[l][p] > flabel_1[lm][p])
	    lm = l;
	}
	label->data[p] = lm;
      }

      for(l = 0; l <= Lmax; l++){
	free(flabel_1[l]);
	free(flabel_2[l]);
      }
      free(mask_nodes); 
      glip::Image32::Destroy(&mask);
    }
    


    float *Relax_dual(sImageGraph *sg,
		      int *S,
		      sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;

      A = sg->A;
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = glip::Image32::Create(label);
      for(p = 0; p < n; p++){
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(glip::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(glip::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      w = (sg->n_link[p])[i];

	      //w = w*w;
	      //w = w*w;
	      //w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(glip::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      //free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      glip::Image32::Destroy(&mask);

      return flabel_1;
    }
    

    void Relax(sImageGraph *sg,
	       int *S,
	       sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;

      glip::ImageGraph::ChangeType(sg, DISSIMILARITY);
      
      A = sg->A;
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = glip::Image32::Create(label);
      for(p = 0; p < n; p++){
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(glip::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(glip::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      w = (sg->n_link[p])[i];
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(glip::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      glip::Image32::Destroy(&mask);
    }


    void ORelax_1(sImageGraph *sg,
		  int *S,
		  sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;
      int *i_inv = NULL;
      int Imax;
      Imax = glip::Image32::GetMaxVal(label);
      glip::ImageGraph::ChangeType(sg, DISSIMILARITY);
      //printf("Wmax: %d\n", sg->Wmax);
      
      A = sg->A;

      i_inv = glip::AdjRel::InverseIndexes(A);
      
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = glip::Image32::Create(label);
      for(p = 0; p < n; p++){
	/*
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	*/
	flabel_1[p] = (float)label->data[p]/(float)Imax;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(glip::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(glip::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      
	      if(flabel_1[q] < 0.5){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  flabel_1[p] = (float)label->data[p]/(float)Imax;
	  /*
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	  */
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(glip::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      glip::Image32::Destroy(&mask);
      free(i_inv);
    }



    float *ORelax_dual(sImageGraph *sg,
		       int *S,
		       sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;
      int *i_inv = NULL;

      A = sg->A;
      i_inv = glip::AdjRel::InverseIndexes(A);
      
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = glip::Image32::Create(label);
      for(p = 0; p < n; p++){
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(glip::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(glip::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;

	      /*
	      if(flabel_1[p] >= flabel_1[q]){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      */
	      if(flabel_1[q] < 0.5){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      
	      //w = w*w;
	      //w = w*w;
	      //w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(glip::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      //free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      glip::Image32::Destroy(&mask);
      free(i_inv);

      return flabel_1;
    }


    

    void ORelax(sImageGraph *sg,
		int *S,
		sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;
      int *i_inv = NULL;
      int Imax;
      Imax = glip::Image32::GetMaxVal(label);
      glip::ImageGraph::ChangeType(sg, DISSIMILARITY);
      
      A = sg->A;
      i_inv = glip::AdjRel::InverseIndexes(A);
      
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = glip::Image32::Create(label);
      for(p = 0; p < n; p++){
	/*
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	*/
	flabel_1[p] = (float)label->data[p]/(float)Imax;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(glip::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(glip::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      
	      if(flabel_1[p] >= flabel_1[q]){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  flabel_1[p] = (float)label->data[p]/(float)Imax;
	  /*
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	  */
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(glip::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      glip::Image32::Destroy(&mask);
      free(i_inv);
    }




    void ORelax_i(sImageGraph *sg,
		  int *S,
		  sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;
      int *i_inv = NULL;

      glip::ImageGraph::ChangeType(sg, DISSIMILARITY);
      
      A = sg->A;
      i_inv = glip::AdjRel::InverseIndexes(A);
      
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = glip::Image32::Create(label);
      for(p = 0; p < n; p++){
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(glip::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(glip::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;
	      
	      if(flabel_1[p] > flabel_1[q]){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(glip::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      glip::Image32::Destroy(&mask);
      free(i_inv);
    }



    void ORelax_s(sImageGraph *sg,
		  int *S,
		  sImage32 *label, int ntimes){
      sImage32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      sAdjRel *A;
      int *mask_nodes;
      int n,p,q,i,j,k,px,py,qx,qy,nlast,ninic;
      float sw,w;
      int *i_inv = NULL;
      //printf("Entrou\n");
      
      glip::ImageGraph::ChangeType(sg, DISSIMILARITY);
      
      A = sg->A;
      i_inv = glip::AdjRel::InverseIndexes(A);
      
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = glip::Image32::Create(label);
      for(p = 0; p < n; p++){
	if(label->data[p] > 0)
	  flabel_1[p] = 1.0;
	else
	  flabel_1[p] = 0.0;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(glip::Image32::IsValidPixel(label, qx,qy)){
	  if(qx >= 0 && qx < label->ncols &&
	     qy >= 0 && qy < label->nrows){
	    q = qx + qy*label->ncols;
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  px = p%label->ncols;
	  py = p/label->ncols;
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    //if(glip::Image32::IsValidPixel(label, qx,qy)){
	    if(qx >= 0 && qx < label->ncols &&
	       qy >= 0 && qy < label->nrows){
	      q = qx + qy*label->ncols;

	      if(flabel_1[p] == flabel_1[q]){
		k = i_inv[i];
		w = ((sg->n_link[p])[i] + (sg->n_link[q])[k])/2.;
	      }
	      else if(flabel_1[p] > flabel_1[q]){
		w = (sg->n_link[p])[i];
	      }
	      else{
		k = i_inv[i];
		w = (sg->n_link[q])[k];
	      }
	      
	      w = sg->Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  if(label->data[p] > 0)
	    flabel_1[p] = 1.0;
	  else
	    flabel_1[p] = 0.0;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    px = p%label->ncols;
	    py = p/label->ncols;
	    for(i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      //if(glip::Image32::IsValidPixel(label, qx,qy)){
	      if(qx >= 0 && qx < label->ncols &&
		 qy >= 0 && qy < label->nrows){
		q = qx + qy*label->ncols;
		
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }


      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      glip::Image32::Destroy(&mask);
      free(i_inv);
    }




    void ORelax_1(sAdjRel3 *A,
		  sScene32 *scn,
		  float per,
		  int *S,
		  sScene32 *label,
		  int ntimes){
      sScene32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      int *mask_nodes;
      int n,p,q,i,j,k,nlast,ninic;
      Voxel u,v;
      float sw,w,per_pq,dmin;
      int Imax,Wmax;
      float *Dpq;

      //--------------------
      dmin = MIN(scn->dx, MIN(scn->dy, scn->dz));
      Dpq = (float *)malloc(A->n*sizeof(float));
      for(i=1; i<A->n; i++){
	Dpq[i] = sqrtf(A->d[i].axis.x*A->d[i].axis.x*scn->dx*scn->dx + 
		       A->d[i].axis.y*A->d[i].axis.y*scn->dy*scn->dy + 
		       A->d[i].axis.z*A->d[i].axis.z*scn->dz*scn->dz)/dmin;
      }
      //--------------------
      Imax = glip::Scene32::GetMaximumValue(label);
      Wmax = glip::Scene32::GetMaximumValue(scn);
      Wmax *= (1.0 + fabsf(per)/100.0);
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = glip::Scene32::Create(label);
      for(p = 0; p < n; p++){
	flabel_1[p] = (float)label->data[p]/(float)Imax;
	u.c.x = glip::Scene32::GetAddressX(label,p);
	u.c.y = glip::Scene32::GetAddressY(label,p);
	u.c.z = glip::Scene32::GetAddressZ(label,p);
	for(i=1; i < A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(glip::Scene32::IsValidVoxel(label, v)){
	    q = glip::Scene32::GetVoxelAddress(label,v);
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  u.c.x = glip::Scene32::GetAddressX(label,p);
	  u.c.y = glip::Scene32::GetAddressY(label,p);
	  u.c.z = glip::Scene32::GetAddressZ(label,p);
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    v.v = u.v + A->d[i].v;
	    if(glip::Scene32::IsValidVoxel(label, v)){
	      q = glip::Scene32::GetVoxelAddress(label,v);
	      w = abs(scn->data[p] - scn->data[q]);
	      if(flabel_1[q] < 0.5)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
		
	      w = Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;
	      w /= Dpq[i];
	      
	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  flabel_1[p] = (float)label->data[p]/(float)Imax;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    u.c.x = glip::Scene32::GetAddressX(label,p);
	    u.c.y = glip::Scene32::GetAddressY(label,p);
	    u.c.z = glip::Scene32::GetAddressZ(label,p);
	    for(i=1; i < A->n; i++){
	      v.v = u.v + A->d[i].v;
	      if(glip::Scene32::IsValidVoxel(label, v)){
		q = glip::Scene32::GetVoxelAddress(label,v);
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }

      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      glip::Scene32::Destroy(&mask);
    }

    


    void ORelax_1(sScene32 *W,
		  sAdjRel3 *A,
		  sScene32 *scn,
		  float per,
		  int *S,
		  sScene32 *label,
		  int ntimes){
      sScene32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      int *mask_nodes;
      int n,p,q,i,j,k,nlast,ninic;
      Voxel u,v;
      float sw,w,per_pq;
      int Imax,Wmax;
      Imax = glip::Scene32::GetMaximumValue(label);
      Wmax = glip::Scene32::GetMaximumValue(W)*2;
      Wmax *= (1.0 + fabsf(per)/100.0);
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = glip::Scene32::Create(label);
      for(p = 0; p < n; p++){
	flabel_1[p] = (float)label->data[p]/(float)Imax;
	u.c.x = glip::Scene32::GetAddressX(label,p);
	u.c.y = glip::Scene32::GetAddressY(label,p);
	u.c.z = glip::Scene32::GetAddressZ(label,p);
	for(i=1; i < A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(glip::Scene32::IsValidVoxel(label, v)){
	    q = glip::Scene32::GetVoxelAddress(label,v);
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  u.c.x = glip::Scene32::GetAddressX(label,p);
	  u.c.y = glip::Scene32::GetAddressY(label,p);
	  u.c.z = glip::Scene32::GetAddressZ(label,p);
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    v.v = u.v + A->d[i].v;
	    if(glip::Scene32::IsValidVoxel(label, v)){
	      q = glip::Scene32::GetVoxelAddress(label,v);
	      w = W->data[p] + W->data[q];
	      if(flabel_1[q] < 0.5)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
		
	      w = Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  flabel_1[p] = (float)label->data[p]/(float)Imax;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    u.c.x = glip::Scene32::GetAddressX(label,p);
	    u.c.y = glip::Scene32::GetAddressY(label,p);
	    u.c.z = glip::Scene32::GetAddressZ(label,p);
	    for(i=1; i < A->n; i++){
	      v.v = u.v + A->d[i].v;
	      if(glip::Scene32::IsValidVoxel(label, v)){
		q = glip::Scene32::GetVoxelAddress(label,v);
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }

      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      glip::Scene32::Destroy(&mask);
    }



    

    void ORelax_1(sScene32 *Wx,
		  sScene32 *Wy,
		  sScene32 *Wz,
		  sScene32 *scn,
		  float per,
		  int *S,
		  sScene32 *label,
		  int ntimes){
      sScene32 *mask;
      float *flabel_1,*flabel_2,*tmp;
      int *mask_nodes;
      int n,p,q,i,j,k,nlast,ninic,t;
      Voxel u,v;
      float sw,w,per_pq;
      int Imax,Wmax;
      sAdjRel3 *A;
      sScene32 *W[3];
      int T[7];
      W[0] = Wx;
      W[1] = Wy;
      W[2] = Wz;
      A = glip::AdjRel3::Spheric(1.0);
      for(i=1; i<A->n; i++){
	if(A->d[i].axis.x != 0) T[i] = 0;
	if(A->d[i].axis.y != 0) T[i] = 1;
	if(A->d[i].axis.z != 0) T[i] = 2;
      }
      Wmax = MAX(glip::Scene32::GetMaximumValue(Wx),
		 MAX(glip::Scene32::GetMaximumValue(Wy),
		     glip::Scene32::GetMaximumValue(Wz)))*2;
      Imax = glip::Scene32::GetMaximumValue(label);

      Wmax *= (1.0 + fabsf(per)/100.0);
      ninic = 1;
      n = label->n;
      flabel_1 = (float *)malloc(sizeof(float)*n);
      flabel_2 = (float *)malloc(sizeof(float)*n);
      mask_nodes = (int *)malloc(sizeof(int)*(n+1));
      mask_nodes[0] = 0;
      mask = glip::Scene32::Create(label);
      for(p = 0; p < n; p++){
	flabel_1[p] = (float)label->data[p]/(float)Imax;
	u.c.x = glip::Scene32::GetAddressX(label,p);
	u.c.y = glip::Scene32::GetAddressY(label,p);
	u.c.z = glip::Scene32::GetAddressZ(label,p);
	for(i=1; i < A->n; i++){
	  v.v = u.v + A->d[i].v;
	  if(glip::Scene32::IsValidVoxel(label, v)){
	    q = glip::Scene32::GetVoxelAddress(label,v);
	    if(label->data[p] != label->data[q]){
	      mask->data[p] = 1;
	      mask_nodes[0]++;
	      mask_nodes[mask_nodes[0]] = p;
	      break;
	    }
	  }
	}
      }

      while(ntimes > 0){
	//Update flabel:
	memcpy(flabel_2, flabel_1, n*sizeof(float));
	for(j = 1; j <= mask_nodes[0]; j++){
	  p = mask_nodes[j];
	  u.c.x = glip::Scene32::GetAddressX(label,p);
	  u.c.y = glip::Scene32::GetAddressY(label,p);
	  u.c.z = glip::Scene32::GetAddressZ(label,p);
	  flabel_2[p] = 0.0;
	  sw = 0.0;
	  for(i=1; i < A->n; i++){
	    v.v = u.v + A->d[i].v;
	    if(glip::Scene32::IsValidVoxel(label, v)){
	      q = glip::Scene32::GetVoxelAddress(label,v);
	      t = T[i];
	      w = W[t]->data[p] + W[t]->data[q];
	      if(flabel_1[q] < 0.5)
		per_pq = per;
	      else
		per_pq = -per;

	      if(scn->data[p] > scn->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(scn->data[p] < scn->data[q])
		w *= (1.0 - per_pq/100.0);
		
	      w = Wmax - w;
	      w = w*w;
	      w = w*w;
	      w = w*w;

	      sw += w;
	      flabel_2[p] += w*flabel_1[q];
	    }
	  }
	  flabel_2[p] /= sw;
	}
	tmp = flabel_1;
	flabel_1 = flabel_2;
	flabel_2 = tmp;

	for(i = 1; i <= S[0]; i++){
	  p = S[i];
	  flabel_1[p] = (float)label->data[p]/(float)Imax;
	}
	
	ntimes--;
	
	if(ntimes > 0){
	  //Dilate mask:
	  nlast = mask_nodes[0];
	  for(j = ninic; j <= mask_nodes[0]; j++){
	    p = mask_nodes[j];
	    u.c.x = glip::Scene32::GetAddressX(label,p);
	    u.c.y = glip::Scene32::GetAddressY(label,p);
	    u.c.z = glip::Scene32::GetAddressZ(label,p);
	    for(i=1; i < A->n; i++){
	      v.v = u.v + A->d[i].v;
	      if(glip::Scene32::IsValidVoxel(label, v)){
		q = glip::Scene32::GetVoxelAddress(label,v);
		if(mask->data[q] == 0){
		  mask->data[q] = 1;
		  nlast++;
		  mask_nodes[nlast] = q;
		}
	      }
	    }
	  }
	  ninic = mask_nodes[0] + 1;
	  mask_nodes[0] = nlast;
	}
      }

      for(p = 0; p < n; p++){
	if(flabel_1[p] < 0.5)
	  label->data[p] = 0;
	else
	  label->data[p] = 1;
      }
      
      free(flabel_1);
      free(flabel_2);
      free(mask_nodes); 
      glip::Scene32::Destroy(&mask);
      glip::AdjRel3::Destroy(&A);
    }


    
    
    void GGC_maxmin(sImageGraph *graph,
		    sCImage *cimg,
		    sImage32 *img,
		    int method,
		    float power_geodesic,
		    float delta,
		    float theta_hedgehog,
		    float R,
		    int postproc,
		    int niterations,
		    int conn,
		    int pol,
		    int shapepriors,
		    int costtemplate,
		    int *S,
		    sImage32 *label){
      int *S1=NULL;
      sGQueue *Q=NULL;
      sImageGraph *sg = NULL;
      sGraph *g = NULL, *transpose = NULL;
      sImage32 *P_sum = NULL, *C_sum = NULL, *raw_map = NULL;
      sImage32 *pred, *cost, *Slabel;
      int i,p,n,ns1;
      //-------------------------
      glip::timer tic,toc;
      float totaltime;
      FILE *fp;
      int reduction = 0;
      //-------------------------      
      float A=750.0, B=75.0, n1=2.0, n2=1.0, n3=1.0;
      int m=2;
      
      Slabel = glip::Image32::Clone(label);
      pred = Image32::Create(label);
      cost = Image32::Create(label);
      Image32::Set(pred, NIL);
      Image32::Set(cost, INT_MAX);
      
      S1 = (int *)malloc((S[0]+1)*sizeof(int));
      ns1 = 0;
      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(label->data[p] == 1){
	  ns1++;
	  S1[ns1] = p;
	}
      }
      S1[0] = ns1;
      S1 = (int *)realloc(S1, (S1[0]+1)*sizeof(int));

      ImageGraph::ChangeType(graph, DISSIMILARITY);
      
      //------------------------
      if(pol == 0 && shapepriors <= 1 && !conn){
	switch(shapepriors){
	case 0:
	  IFT_fw(graph, S, label, cost, pred);
	  break;
	case 1:
	  P_sum = SC_Pred_fsum(graph, S1, power_geodesic);
	  SC_IFT(graph, S, label, P_sum);
	  Image32::Destroy(&P_sum);
	  break;
	}
      }
      else{
	sg = ImageGraph::Clone(graph);
	if(cimg != NULL){ //COLOR_IMAGE
	  sImage32 *lumi;
	  lumi = Image32::Luminosity(cimg);
	  ImageGraph::Orient2Digraph(sg, lumi, pol);
	  Image32::Destroy(&lumi);
	}
	else
	  ImageGraph::Orient2Digraph(sg, img, pol);
	
	if(shapepriors != 0){
	  if(costtemplate == 0){
	    P_sum = SC_Pred_fsum(sg, S1, power_geodesic);
	    C_sum = BB_Geodesic_Cost(P_sum, graph->A);
	  }
	  else{
	    switch(costtemplate){
	    case 1: //Circle
	      raw_map = Image32::Read("./templates/circle.pgm"); break;
	    case 2: //Square
	      raw_map = Image32::Read("./templates/square.pgm"); break;
	    case 3:
	      raw_map = Image32::Read("./templates/ellipse2:1.pgm"); break;
	    case 4:
	      raw_map = Image32::Read("./templates/ellipse3:1.pgm"); break;
	    case 5:
	      raw_map = Image32::Read("./templates/rectangle2:1.pgm"); break;
	    case 6:
	      raw_map = Image32::Read("./templates/rectangle3:1.pgm"); break;
	    case 7:
	      fp = fopen("gielis.txt", "r");
	      if(fp != NULL){
		fscanf(fp," %f %f %f %f %f %d", &A, &B, &n1, &n2, &n3, &m);
		fclose(fp);
	      }
	      raw_map = glip::Band::GielisEquationTemplate(2001, A, B,
							  n1, n2, n3, m);
	      break;	      
	    }
	    
	    C_sum = BB_CropTemplate(raw_map, S, label, 1);
	    P_sum = NULL;
	    Image32::Destroy(&raw_map);
	  }
	}
	
	if(shapepriors == 1 && costtemplate == 0){ //Star Convexity
	  ImageGraph::Orient2DigraphOuter(sg, P_sum);
	}
	else if(shapepriors == 2){ //Hedgehog
	  if(R <= 1.6 && sg->A->n == 9)
	    ImageGraph::Convert2HedgehogDigraph(sg, C_sum, theta_hedgehog);
	  else{
	    g = Graph::Clone(sg);
	    Graph::HedgehogDigraph(g, C_sum, theta_hedgehog, R);
	    transpose = Graph::Transpose(g);
	  }
	}
	else if(shapepriors == 3){ //Local Band
	  g = Graph::Clone(sg);
	  printf("NumberOfArcs: %d (image graph)\n", Graph::GetNumberOfArcs(g));
	  Graph::LocalBandConstraint(g, C_sum, ROUND(delta), R);
	  printf("NumberOfArcs: %d (local band)\n", Graph::GetNumberOfArcs(g));
	  fp = fopen("tmp.txt", "r");
	  if(fp !=NULL){
	    fscanf(fp, "%d", &reduction);
	    fclose(fp);
	  }
	  if(reduction){
	    gettimeofday(&tic,NULL);
	    Graph::LocalBandReduction(&g, C_sum, ROUND(delta));
	    gettimeofday(&toc,NULL);
	    totaltime = (toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001;
	    printf("\nReduction Time: %f ms\n",totaltime);	    
	    printf("NumberOfArcs: %d (reduction)\n", Graph::GetNumberOfArcs(g));
	  }
	  transpose = Graph::Transpose(g);
	}
	
	if(conn){
	  COIFT(sg, S, label);
	}
	else if(shapepriors < 4 && g == NULL){
	  gettimeofday(&tic,NULL);
	  switch(method){
	  case 0:
	    OIFT(sg, NULL, label); break;
	  case 1:
	    EOIFT(sg, NULL, label, 0); break;
	  case 2:
	    OIFT_Heap(sg, NULL, label);	break;
	  case 3:
	    EOIFT_Heap(sg, NULL, label, 0.0); break;
	  case 4:
	    OIFT_TZ2Bkg(sg, NULL, label); break;
	  case 5:
	    OIFT_TZ2Obj(sg, NULL, label); break;	    
	  case 6:
	    EOIFT_Heap_2(sg, NULL, label, 0.0); break;
	  }
	  gettimeofday(&toc,NULL);
	  totaltime = (toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001;
	  printf("\nInner Time: %f ms\n",totaltime);

	  /****** REMOVER *****/

	  printf("energy: %d\n", GetEnergy_Min(sg, label, 1));
	  
	  if( isOIFT_Segmentation(sg, S, Slabel, label) )
	    printf("isOIFT: OK\n");
	  else
	    printf("isOIFT: Error\n");
	  
	  /********************/
	}
	else if(shapepriors < 4 && g != NULL){
	  gettimeofday(&tic,NULL);
	  switch(method){
	  case 0:
	    OIFT(g, transpose, NULL, label->data); break;
	  case 1:
	    EOIFT(g, transpose, NULL, label->data, 0); break;
	  case 2:
	    OIFT_Heap(g, transpose, NULL, label->data); break;
	  case 3:
	    EOIFT_Heap(g, transpose, NULL, label->data, 0.0); break;
	  case 6:
	    EOIFT_Heap_2(g, transpose, NULL, label->data, 0.0); break;
	  }
	  gettimeofday(&toc,NULL);
	  totaltime = (toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001;
	  printf("\nInner Time: %f ms\n",totaltime);
	  printf("Energy: %d\n", GetEnergy_Min(g, label->data, 1));
	}
	else{
	  if(shapepriors == 4)
	    BB_OIFT(sg, S, label,
		    C_sum, P_sum, delta);
	  else if(shapepriors == 5){
	    //---code included to fix errors with zero costs---
	    //---because delta*C(t) cannot be zero-------------
	    for(p = 0; p < C_sum->n; p++) 
	      C_sum->data[p] += 1;
	    //-------------------------------------------------
	    RBB_OIFT(sg, S, label,
		     C_sum, P_sum, delta);
	  }
	  else if(shapepriors == 6){
	    B_OIFT(sg, S, label, C_sum, delta);
	  }
	}
	//printf("energy: %d\n", glip::ift::GetEnergy_Min(sg, APP->Data.label, 1));
      }
      
      //--------REMOVER----------------
      /*      
      if(sg != NULL)
	ImageGraph::Destroy(&sg);
      sg = ImageGraph::Clone(graph);
      ImageGraph::Orient2Digraph(sg, img, pol);
      */
      //-------------------------------
      
      switch(postproc){
      case 1:
	if(sg != NULL)
	  Relax(sg, S, label, niterations);
	else
	  Relax(graph, S, label, niterations);
	break;
      case 2:
	if(sg != NULL)
	  ORelax(sg, S, label, niterations);
	else
	  ORelax(graph, S, label, niterations);
	break;
      case 3:
	if(sg != NULL)
	  ORelax_1(sg, S, label, niterations);
	else
	  ORelax_1(graph, S, label, niterations);
	break;
      case 4:
	Image32::ModeFilterLabel(label, (float)niterations);
	break;
      case 5:
	if(sg != NULL)
	  ORelax_i(sg, S, label, niterations);
	else
	  ORelax_i(graph, S, label, niterations);
	break;
      case 6:
	if(sg != NULL)
	  ORelax_s(sg, S, label, niterations);
	else
	  ORelax_s(graph, S, label, niterations);
	break;
      case 7:
	Relax(graph, S, label, niterations);
	break;
      }

      /****** REMOVER *****/
      /*
      if(sg != NULL){
	if( isOIFT(sg, S, Slabel, label) )
	  printf("pos isOIFT: OK\n");
	else
	  printf("pos isOIFT: Error\n");
      }
      */
      /********************/
      
      //------------------------

      if(Slabel != NULL) Image32::Destroy(&Slabel);
      
      if(sg != NULL) ImageGraph::Destroy(&sg);
      if(g  != NULL) Graph::Destroy(&g);
      if(transpose != NULL) Graph::Destroy(&transpose);
      if(C_sum != NULL) Image32::Destroy(&C_sum);
      if(P_sum != NULL) Image32::Destroy(&P_sum);
      if(pred != NULL) Image32::Destroy(&pred);
      if(cost != NULL) Image32::Destroy(&cost);
      free(S1);
    }



    //--------------------------------------

    
    void DOIFT_removeSubTree(int q_in,
			     sHeap *Q,
			     sImage32 *pred,
			     sImage32 *root,
			     sAdjRel *A){
      std::queue<int> path, frontier_path;
      int i, p, q;
      int px,py,qx,qy;
      
      path.push(q_in);
      
      while (!path.empty()) {
	p = path.front();
	path.pop();
	px = p%root->ncols;
	py = p/root->ncols;
	
	pred->data[p] = NIL;
	root->data[p] = p;
	
	if (Q->color[p] == GRAY)
	  glip::Heap::Delete_MinPolicy(Q, p);
	Q->color[p] = WHITE;
	Q->cost[p] = FLT_MAX;
	
	for (i = 1; i < A->n; i++) {
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  if(glip::Image32::IsValidPixel(root, qx, qy)){
	    q = qx + qy*root->ncols;
	    if (p == pred->data[q])
	      path.push(q);
	    else if(Q->cost[q] != FLT_MAX)
	      frontier_path.push(q);
	  }
	}
      }
      
      while (!frontier_path.empty()) {
	p = frontier_path.front();
	frontier_path.pop();
	if (Q->cost[p] != FLT_MAX) {
	  if (Q->color[p] != BLACK) {
	    glip::Heap::Update_MinPolicy(Q, p, Q->cost[p]);
	  } else {
	    Q->color[p] = WHITE;
	    glip::Heap::Insert_MinPolicy(Q, p);
	  }
	}
      }
    }




    void DOIFT_removeSubTree(int q_in,
			     sHeap *Q,
			     int *pred,
			     int *root,
			     sGraph *graph){
      std::queue<int> path, frontier_path;
      int i, p, q;
      
      path.push(q_in);
      
      while (!path.empty()) {
	p = path.front();
	path.pop();
	
	pred[p] = NIL;
	root[p] = p;
	
	if (Q->color[p] == GRAY)
	  glip::Heap::Delete_MinPolicy(Q, p);
	Q->color[p] = WHITE;
	Q->cost[p] = FLT_MAX;

	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if (p == pred[q])
	    path.push(q);
	  else if(Q->cost[q] != FLT_MAX)
	    frontier_path.push(q);
	}
      }
      
      while (!frontier_path.empty()) {
	p = frontier_path.front();
	frontier_path.pop();
	if (Q->cost[p] != FLT_MAX) {
	  if (Q->color[p] != BLACK) {
	    glip::Heap::Update_MinPolicy(Q, p, Q->cost[p]);
	  } else {
	    Q->color[p] = WHITE;
	    glip::Heap::Insert_MinPolicy(Q, p);
	  }
	}
      }
    }

    

    //1st from DOIFT paper.
    void DOIFT(sImageGraph *sg,
	       int *S, //new seeds only.
	       int *R, //new removal only.
	       sImage32 *label,
	       sHeap *Q,
	       sImage32 *pred,
	       sImage32 *root){
      //glip::sQueue *frontier = NULL;
      int *i_inv;
      int u_x,u_y,v_x,v_y;
      int p,q,i,j,rq,rp;
      float w,tmp;
      
      i_inv = glip::AdjRel::InverseIndexes(sg->A);

      if(R != NULL && R[0] > 0){
	DOIFT_TreeRemoval(R, sg, Q,
			  pred, root);

	/*
	while(!glip::Queue::IsEmpty(frontier)){
	  p = glip::Queue::Pop(frontier);
	  if (Q->color[p] != BLACK) {
	    glip::Heap::Update_MinPolicy(Q, p, Q->cost[p]);
	  } else {
	    Q->color[p] = WHITE;
	    glip::Heap::Insert_MinPolicy(Q, p);
	  }
	  Q->color[p] = GRAY;
	}
	*/
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(Q->cost[p] > 0.0){
	  pred->data[p] = NIL;
	  root->data[p] = p;
	  
	  if (Q->color[p] != BLACK) {
	    glip::Heap::Update_MinPolicy(Q, p, 0.0);
	  } else {
	    Q->color[p] = WHITE;
	    Q->cost[p] = 0.0;
	    glip::Heap::Insert_MinPolicy(Q, p);
	  }
	  Q->color[p] = GRAY;
	}
      }

      //----------
      while(!Heap::IsEmpty(Q)) {
	Heap::Remove_MinPolicy(Q, &p);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	for(i=1; i<sg->A->n; i++){
	  v_x = u_x + sg->A->dx[i];
	  v_y = u_y + sg->A->dy[i];
	  if(Image32::IsValidPixel(label,v_x,v_y)){
	    q = v_x + label->ncols*v_y;
	    
	    //if (Q->color[q] != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i];
		w = 2.0*(sg->n_link[q])[j];
	      }
	      else
		w = 2.0*(sg->n_link[p])[i] + 1.0;
	      
	      tmp = MAX(w, Q->cost[p]);
	      
	      rq = root->data[q];
	      rp = root->data[p];
	      
	      if(tmp < Q->cost[q]){
		label->data[q] = label->data[p];
		root->data[q] = root->data[p];
		pred->data[q] = p;
		
		if (Q->color[q] != BLACK) {
		  glip::Heap::Update_MinPolicy(Q, q, tmp);
		} else {
		  Q->color[q] = WHITE;
		  Q->cost[q] = tmp;
		  glip::Heap::Insert_MinPolicy(Q, q);
		}
	      }
	      else if(p == pred->data[q]){
		if(rp != rq || tmp > Q->cost[q]){
		  DOIFT_removeSubTree(q,
				      Q,
				      pred,
				      root,
				      sg->A);
		  break;
		}
	      }
	   //}
	  }
	}
      }
      //----------
      //glip::Queue::Destroy(&frontier);
      free(i_inv);
    }




    //1st from DOIFT paper.
    void DOIFT(sGraph *graph,
	       sGraph *transpose,
	       int *S, //new seeds only.
	       int *R, //new removal only.
	       int *label,
	       sHeap *Q,
	       int *pred,
	       int *root){
      int p,q,i,j,rq,rp;
      float w,tmp;
      sGraph *g;
      
      if(R != NULL && R[0] > 0){
	DOIFT_TreeRemoval(R, graph, transpose, Q,
			  pred, root);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	if(Q->cost[p] > 0.0){
	  pred[p] = NIL;
	  root[p] = p;
	  
	  if (Q->color[p] != BLACK) {
	    glip::Heap::Update_MinPolicy(Q, p, 0.0);
	  } else {
	    Q->color[p] = WHITE;
	    Q->cost[p] = 0.0;
	    glip::Heap::Insert_MinPolicy(Q, p);
	  }
	  Q->color[p] = GRAY;
	}
      }

      //----------
      while(!Heap::IsEmpty(Q)) {
	Heap::Remove_MinPolicy(Q, &p);

	if(label[p]==0) g = transpose;
	else   	        g = graph;
	
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  w = g->nodes[p].Warcs[i];
	  if(label[p] == 0) w = 2.0*w;
	  else              w = 2.0*w + 1.0;
	    
	  tmp = MAX(w, Q->cost[p]);
	  
	  rq = root[q];
	  rp = root[p];
	    
	  if(tmp < Q->cost[q]){
	    label[q] = label[p];
	    root[q] = root[p];
	    pred[q] = p;
	    
	    if (Q->color[q] != BLACK) {
	      glip::Heap::Update_MinPolicy(Q, q, tmp);
	    } else {
	      Q->color[q] = WHITE;
	      Q->cost[q] = tmp;
	      glip::Heap::Insert_MinPolicy(Q, q);
	    }
	  }
	  else if(p == pred[q]){
	    if(rp != rq || tmp > Q->cost[q]){
	      DOIFT_removeSubTree(q,
				  Q,
				  pred,
				  root,
				  graph);
	      break;
	    }
	  }
	}
      }
    }
    


    void DOIFT_TreeRemoval(int *R,
			   sImageGraph *sg,
			   sHeap *Q,
			   sImage32 *pred,
			   sImage32 *root){
      glip::sQueue *T;
      int i,k,r,p,q;
      int px,py,qx,qy;
      
      T = glip::Queue::Create(root->n);

      for(i = 1; i <= R[0]; i++){
	r = root->data[R[i]];
	if(Q->cost[r] != FLT_MAX){
	  if (Q->color[r] == GRAY)
	    glip::Heap::Delete_MinPolicy(Q, r);
	  Q->cost[r] = FLT_MAX;
	  Q->color[r] = WHITE;
	  pred->data[r] = NIL;
	  root->data[r] = r;
	  glip::Queue::Push(T, r);
	}
      }
      while (!glip::Queue::IsEmpty(T)) {
	p = glip::Queue::Pop(T);
	px = p%root->ncols;
	py = p/root->ncols;

	for (i = 1; i < (sg->A)->n; i++) {
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  if(qx >= 0 && qx < root->ncols &&
	     qy >= 0 && qy < root->nrows){
	    q = qx + qy*root->ncols;
	    if(p == pred->data[q]){
	      if (Q->color[q] == GRAY)
		glip::Heap::Delete_MinPolicy(Q, q);
	      Q->cost[q] = FLT_MAX;
	      Q->color[q] = WHITE;
	      pred->data[q] = NIL;
	      root->data[q] = q;
	      glip::Queue::Push(T, q);
	    }
	    else if(Q->cost[root->data[q]] != FLT_MAX){
	      if(Q->color[q] != GRAY){
		Q->color[q] = WHITE;
		glip::Heap::Insert_MinPolicy(Q, q);
	      }
	    }
	  }
	}
      }
      glip::Queue::Destroy(&T);
    }




    void DOIFT_TreeRemoval(int *R,
			   sGraph *graph,
			   sGraph *transpose,
			   sHeap *Q,
			   int *pred,
			   int *root){
      glip::sQueue *T;
      int i,k,r,p,q;
      
      T = glip::Queue::Create(graph->nnodes);

      for(i = 1; i <= R[0]; i++){
	r = root[R[i]];
	if(Q->cost[r] != FLT_MAX){
	  if (Q->color[r] == GRAY)
	    glip::Heap::Delete_MinPolicy(Q, r);
	  Q->cost[r] = FLT_MAX;
	  Q->color[r] = WHITE;
	  pred[r] = NIL;
	  root[r] = r;
	  glip::Queue::Push(T, r);
	}
      }
      while (!glip::Queue::IsEmpty(T)) {
	p = glip::Queue::Pop(T);

	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];

	  if(p == pred[q]){
	    if (Q->color[q] == GRAY)
	      glip::Heap::Delete_MinPolicy(Q, q);
	    Q->cost[q] = FLT_MAX;
	    Q->color[q] = WHITE;
	    pred[q] = NIL;
	    root[q] = q;
	    glip::Queue::Push(T, q);
	  }
	  else if(Q->cost[root[q]] != FLT_MAX){
	    if(Q->color[q] != GRAY){
	      Q->color[q] = WHITE;
	      glip::Heap::Insert_MinPolicy(Q, q);
	    }
	  }
	}
      }
      glip::Queue::Destroy(&T);
    }

    
    //----------------------------------------------------
    
    void DOIFT_TreeRemoval_v2_v3(int *R,
				 sImageGraph *sg,
				 sHeap32fi_lex *Q,
				 sImage32 *pred,
				 sImage32 *maxorder,
				 int *T){
      glip::sBMap *Fcolor;
      int i,j,k,r,p,q;
      int px,py,qx,qy;
      int top = pred->n+1;
      
      T[0] = 0;
      for(j = 1; j <= R[0]; j++){
	r = R[j];
	T[0]++;
	T[T[0]] = r;

	if (Q->color[r] == GRAY)
	  glip::Heap32fi_lex::Delete_MinPolicy(Q, r);
	Q->cost1[r] = FLT_MAX;
	maxorder->data[r] = INT_MAX;
	Q->color[r] = WHITE;
	pred->data[r] = NIL;
      }

      //---------------------------------------
      while(T[0] > 0){
	p = T[T[0]];
	T[0]--;
	
	px = p%pred->ncols;
	py = p/pred->ncols;
	for(i = 1; i < (sg->A)->n; i++){
	  qx = px + (sg->A)->dx[i];
	  qy = py + (sg->A)->dy[i];
	  //if(glip::Image32::IsValidPixel(pred, qx, qy)){
	  if(qx >= 0 && qx < pred->ncols &&
	     qy >= 0 && qy < pred->nrows){
	    q = qx + qy*pred->ncols;
	    if(p == pred->data[q]){
	      T[0]++;
	      T[T[0]] = q;

	      if (Q->color[q] == GRAY)
		glip::Heap32fi_lex::Delete_MinPolicy(Q, q);
	      Q->cost1[q] = FLT_MAX;
	      maxorder->data[q] = INT_MAX;
	      Q->color[q] = WHITE;
	      pred->data[q] = NIL;
	    }
	    else if(Q->cost1[q] != FLT_MAX &&
		    Q->color[q] != GRAY){
	      top--;
	      T[top] = q;
	    }
	  }
	}
      }
      //---------------------------------------	
      //printf("per: %f\n",(float)(pred->n+1-top)/((float)(pred->n)));
      while(top < pred->n+1){
	q = T[top];
	top++;
	if(Q->cost1[q] != FLT_MAX && Q->color[q] != GRAY){
	  Q->color[q] = WHITE;
	  glip::Heap32fi_lex::Insert_MinPolicy(Q, q);	  
	}
      }
    }





    void DOIFT_TreeRemoval_v2_v3(int *R,
				 sGraph *graph,
				 sHeap32fi_lex *Q,
				 int *pred,
				 int *maxorder,
				 int *T){
      glip::sBMap *Fcolor;
      int i,j,k,r,p,q;
      int top = graph->nnodes+1;
      
      T[0] = 0;
      for(j = 1; j <= R[0]; j++){
	r = R[j];
	T[0]++;
	T[T[0]] = r;

	if (Q->color[r] == GRAY)
	  glip::Heap32fi_lex::Delete_MinPolicy(Q, r);
	Q->cost1[r] = FLT_MAX;
	maxorder[r] = INT_MAX;
	Q->color[r] = WHITE;
	pred[r] = NIL;
      }

      //---------------------------------------
      while(T[0] > 0){
	p = T[T[0]];
	T[0]--;
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];

	  if(p == pred[q]){
	    T[0]++;
	    T[T[0]] = q;

	    if (Q->color[q] == GRAY)
	      glip::Heap32fi_lex::Delete_MinPolicy(Q, q);
	    Q->cost1[q] = FLT_MAX;
	    maxorder[q] = INT_MAX;
	    Q->color[q] = WHITE;
	    pred[q] = NIL;
	  }
	  else if(Q->cost1[q] != FLT_MAX &&
		  Q->color[q] != GRAY){
	    top--;
	    T[top] = q;
	  }
	}
      }
      //---------------------------------------	
      while(top < graph->nnodes+1){
	q = T[top];
	top++;
	if(Q->cost1[q] != FLT_MAX && Q->color[q] != GRAY){
	  Q->color[q] = WHITE;
	  glip::Heap32fi_lex::Insert_MinPolicy(Q, q);	  
	}
      }
    }

    ////DOIFT__2 -> DOIFT_v2    2nd from DOIFT paper
    void DOIFT_v2(sImageGraph *sg,
		  int *S, //new seeds only.
		  int *R, //new removal only.
		  sImage32 *label,
		  sHeap32fi_lex *Q,
		  sImage32 *pred,
		  sImage32 *maxorder,
		  int *iter){
      int *i_inv;
      int u_x,u_y,v_x,v_y;
      int p,q,i,j,tq,cm;
      float w,Wm;
      int Removal[2];
      int maxord;
      int *T = NULL;
      
      i_inv = glip::AdjRel::InverseIndexes(sg->A);

      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((pred->n+1)*sizeof(int));
	DOIFT_TreeRemoval_v2_v3(R, sg, Q, pred,
				maxorder, T);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred->data[p] = NIL;
	if(label->data[p]==0)
	  cm = 1; //1;
	else
	  cm = 2; //2;
	  
	if (Q->color[p] != BLACK) {
	  glip::Heap32fi_lex::Update_MinPolicy(Q, p, 0.0, cm); //INT_MAX);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder->data[p] = cm; //INT_MAX;
	  glip::Heap32fi_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);

	/*
	if(maxorder->data[p] == INT_MAX){
	  (*iter)++;
	  maxorder->data[p] = *iter;
	}
	*/
	//-------------------------------------------------------
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<sg->A->n; i++){
	  v_x = u_x + sg->A->dx[i];
	  v_y = u_y + sg->A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){	  
	    q = v_x + label->ncols*v_y;
	    
	    if(label->data[p]==0){
	      j = i_inv[i];
	      w = (sg->n_link[q])[j];
	      cm = 1; //1;
	    }
	    else{
	      w = (sg->n_link[p])[i];
	      cm = 2; //2;
	    }
	      
	    Wm = MAX(w, Q->cost1[p]);
	    
	    if(Wm != w)
	      maxord = maxorder->data[p];
	    else if( Q->cost1[p] == Wm )
	      maxord = maxorder->data[p] + 2;
	    else
	      maxord = cm;
	      
	    if((Wm < Q->cost1[q] || (Wm == Q->cost1[q] && maxord < maxorder->data[q]))){
	      label->data[q] = label->data[p];
	      pred->data[q] = p;
	      
	      if (Q->color[q] != BLACK) {
		glip::Heap32fi_lex::Update_MinPolicy(Q, q, Wm, maxord);
	      } else {
		Q->color[q] = WHITE;
		Q->cost1[q] = Wm;
		maxorder->data[q] = maxord;
		glip::Heap32fi_lex::Insert_MinPolicy(Q, q);
	      }
	    }
	    else if(p == pred->data[q]){
	      if((label->data[q] != label->data[p])||
		 (Wm != Q->cost1[q]) ||
		 (maxord > maxorder->data[q])){
		if(T == NULL)
		  T = (int *)malloc((pred->n+1)*sizeof(int));

		Removal[0] = 1; Removal[1] = q;
		DOIFT_TreeRemoval_v2_v3(Removal, sg, Q, pred,
					maxorder, T);
		break;
	      }
	    }
	  }
	}
	//-------------------------------------------------------
      }
      //----------
      free(i_inv);
      if(T != NULL)
	free(T);
    }


    //DOIFT__2 -> DOIFT_v2    2nd from DOIFT paper
    void DOIFT_v2(sGraph *graph,
		  sGraph *transpose,
		  int *S,
		  int *R,
		  int *label,
		  sHeap32fi_lex *Q,
		  int *pred,
		  int *maxorder,
		  int *iter){
      int p,q,i,j,tq,cm;
      float w,Wm;
      int Removal[2];
      int maxord;
      int *T = NULL;
      sGraph *g;

      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	DOIFT_TreeRemoval_v2_v3(R, graph, Q, pred,
				maxorder, T);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred[p] = NIL;
	if(label[p]==0)
	  cm = 1; //1;
	else
	  cm = 2; //2;
	  
	if (Q->color[p] != BLACK) {
	  glip::Heap32fi_lex::Update_MinPolicy(Q, p, 0.0, cm);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder[p] = cm;
	  glip::Heap32fi_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);
	
	if(label[p]==0){ g = transpose; cm = 1; } //1;
	else{   	 g = graph;     cm = 2; } //2;
	//-------------------------------------------------------
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  w = g->nodes[p].Warcs[i];
	  
	  Wm = MAX(w, Q->cost1[p]);
	    
	  if(Wm != w)
	    maxord = maxorder[p];
	  else if( Q->cost1[p] == Wm )
	    maxord = maxorder[p] + 2;
	  else
	    maxord = cm;
	      
	  if((Wm < Q->cost1[q] || (Wm == Q->cost1[q] && maxord < maxorder[q]))){
	    label[q] = label[p];
	    pred[q] = p;
	      
	    if (Q->color[q] != BLACK) {
	      glip::Heap32fi_lex::Update_MinPolicy(Q, q, Wm, maxord);
	    } else {
	      Q->color[q] = WHITE;
	      Q->cost1[q] = Wm;
	      maxorder[q] = maxord;
	      glip::Heap32fi_lex::Insert_MinPolicy(Q, q);
	    }
	  }
	  else if(p == pred[q]){
	    if((label[q] != label[p])||
	       (Wm != Q->cost1[q]) ||
	       (maxord > maxorder[q])){
	      if(T == NULL)
		T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	      
	      Removal[0] = 1; Removal[1] = q;
	      DOIFT_TreeRemoval_v2_v3(Removal, graph, Q, pred,
				      maxorder, T);
	      break;
	    }
	  }
	}
      
	//-------------------------------------------------------
      }
      //----------
      if(T != NULL)
	free(T);
    }



    //DOIFT__2p -> DOIFT_v3   3th from DOIFT paper
    void DOIFT_v3(sImageGraph *sg,
		  int *S, //new seeds only.
		  int *R, //new removal only.
		  sImage32 *label,
		  sHeap32fi_lex *Q,
		  sImage32 *pred,
		  sImage32 *maxorder,
		  int *iter){
      int *i_inv;
      int u_x,u_y,v_x,v_y;
      int p,q,i,j,tq,cm;
      float w,Wm;
      int Removal[2];
      int maxord;
      int *T = NULL;
      bool seenfirst_pq, seenfirst_qp, aux;
      
      i_inv = glip::AdjRel::InverseIndexes(sg->A);

      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((pred->n+1)*sizeof(int));
	DOIFT_TreeRemoval_v2_v3(R, sg, Q, pred,
				maxorder, T);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred->data[p] = NIL;
	if(label->data[p]==0)
	  cm = 1; //1;
	else
	  cm = 2; //2;
	  
	if (Q->color[p] != BLACK) {
	  glip::Heap32fi_lex::Update_MinPolicy(Q, p, 0.0, cm);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder->data[p] = cm;
	  glip::Heap32fi_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);
	//-------------------------------------------------------
	u_x = p%label->ncols; 
	u_y = p/label->ncols; 
	
	for(i=1; i<sg->A->n; i++){
	  v_x = u_x + sg->A->dx[i];
	  v_y = u_y + sg->A->dy[i];
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){	  
	    q = v_x + label->ncols*v_y;
	    
	    if(label->data[p]==0){
	      j = i_inv[i];
	      w = (sg->n_link[q])[j];
	      cm = 1; //1;
	    }
	    else{
	      w = (sg->n_link[p])[i];
	      cm = 2; //2;
	    }
	      
	    Wm = MAX(w, Q->cost1[p]);
	    
	    if(Wm != w)
	      maxord = maxorder->data[p];
	    else if( Q->cost1[p] == Wm )
	      maxord = maxorder->data[p] + 2;
	    else
	      maxord = cm;

	    seenfirst_pq = seenfirst_qp = false;
	    if(pred->data[q] != NIL){
	      aux = (w > Q->cost1[p] && Q->cost1[q] > Q->cost1[pred->data[q]] && w == Q->cost1[q]);
	      seenfirst_pq = aux && (Q->cost1[p] < Q->cost1[pred->data[q]] ||
				     (Q->cost1[p] == Q->cost1[pred->data[q]] &&
				      maxorder->data[p] < maxorder->data[pred->data[q]]));
	      seenfirst_qp = aux && (Q->cost1[p] > Q->cost1[pred->data[q]] ||
				     (Q->cost1[p] == Q->cost1[pred->data[q]] &&
				      maxorder->data[p] > maxorder->data[pred->data[q]]));
	    }
	    
	    if(Wm < Q->cost1[q] ||
	       (Wm == Q->cost1[q] && ((maxord < maxorder->data[q] && !seenfirst_qp) || seenfirst_pq))){
	      label->data[q] = label->data[p];
	      pred->data[q] = p;
	      if(Q->color[q] == GRAY)
		glip::Heap32fi_lex::Delete_MinPolicy(Q, q);
	      Q->color[q] = WHITE;
	      Q->cost1[q] = Wm;
	      maxorder->data[q] = maxord;
	      glip::Heap32fi_lex::Insert_MinPolicy(Q, q);
	    }
	    else if(p == pred->data[q]){
	      if((label->data[q] != label->data[p])||
		 (Wm != Q->cost1[q]) ||
		 (maxord > maxorder->data[q])){
		if(T == NULL)
		  T = (int *)malloc((pred->n+1)*sizeof(int));

		Removal[0] = 1; Removal[1] = q;
		DOIFT_TreeRemoval_v2_v3(Removal, sg, Q, pred,
					maxorder, T);
		break;
	      }
	    }
	  }
	}
	//-------------------------------------------------------
      }
      //----------
      free(i_inv);
      if(T != NULL)
	free(T);
    }



    //DOIFT__2p -> DOIFT_v3   3th from DOIFT paper
    void DOIFT_v3(sGraph *graph,
		  sGraph *transpose,
		  int *S,
		  int *R,
		  int *label,
		  sHeap32fi_lex *Q,
		  int *pred,
		  int *maxorder,
		  int *iter){
      int p,q,i,j,tq,cm;
      float w,Wm;
      int Removal[2];
      int maxord;
      int *T = NULL;
      sGraph *g;
      bool seenfirst_pq, seenfirst_qp, aux;
      
      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	DOIFT_TreeRemoval_v2_v3(R, graph, Q, pred,
				maxorder, T);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred[p] = NIL;
	if(label[p]==0)
	  cm = 1; //1;
	else
	  cm = 2; //2;
	  
	if (Q->color[p] != BLACK) {
	  glip::Heap32fi_lex::Update_MinPolicy(Q, p, 0.0, cm);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder[p] = cm;
	  glip::Heap32fi_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);
	
	if(label[p]==0){ g = transpose; cm = 1; } //1;
	else{   	 g = graph;     cm = 2; } //2;
	//-------------------------------------------------------
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  w = g->nodes[p].Warcs[i];
	  
	  Wm = MAX(w, Q->cost1[p]);
	    
	  if(Wm != w)
	    maxord = maxorder[p];
	  else if( Q->cost1[p] == Wm )
	    maxord = maxorder[p] + 2;
	  else
	    maxord = cm;

	  seenfirst_pq = seenfirst_qp = false;
	  if(pred[q] != NIL){
	    aux = (w > Q->cost1[p] && Q->cost1[q] > Q->cost1[pred[q]] && w == Q->cost1[q]);
	    seenfirst_pq = aux && (Q->cost1[p] < Q->cost1[pred[q]] ||
				   (Q->cost1[p] == Q->cost1[pred[q]] &&
				    maxorder[p] < maxorder[pred[q]]));
	    seenfirst_qp = aux && (Q->cost1[p] > Q->cost1[pred[q]] ||
				   (Q->cost1[p] == Q->cost1[pred[q]] &&
				    maxorder[p] > maxorder[pred[q]]));
	  }
	  
	  if(Wm < Q->cost1[q] ||
	     (Wm == Q->cost1[q] && ((maxord < maxorder[q] && !seenfirst_qp) || seenfirst_pq))){
	    label[q] = label[p];
	    pred[q] = p;
	    if(Q->color[q] == GRAY)
	      glip::Heap32fi_lex::Delete_MinPolicy(Q, q);
	    Q->color[q] = WHITE;
	    Q->cost1[q] = Wm;
	    maxorder[q] = maxord;
	    glip::Heap32fi_lex::Insert_MinPolicy(Q, q);
	  }
	  else if(p == pred[q]){
	    if((label[q] != label[p])||
	       (Wm != Q->cost1[q]) ||
	       (maxord > maxorder[q])){
	      if(T == NULL)
		T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	      
	      Removal[0] = 1; Removal[1] = q;
	      DOIFT_TreeRemoval_v2_v3(Removal, graph, Q, pred,
				      maxorder, T);
	      break;
	    }
	  }
	}
      
	//-------------------------------------------------------
      }
      //----------
      if(T != NULL)
	free(T);
    }



    //DOIFT__2p_hierarchy -> DOIFT_v3_hierarchy
    void DOIFT_v3_hierarchy(sGraph *graph,
			    sGraph *transpose,
			    int *S,
			    int *R,
			    int *label,
			    sHeap32fi_lex *Q,
			    int *pred,
			    int *maxorder,
			    sHeap32 *E,
			    int *H,
			    int *obj_area){
      int p,q,b,i,j,tq,cm;
      float w,Wm;
      int Removal[2];
      int maxord, T_area;
      int *T = NULL;
      sGraph *g;
      bool seenfirst_pq, seenfirst_qp, aux;

      T_area = *obj_area;
      
      if(R != NULL && R[0] > 0){
	if(T == NULL)
	  T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	DOIFT_TreeRemoval_v2_v3(R, graph, Q, pred,
				maxorder, T);
      }

      for(i = 1; i <= S[0]; i++){
	p = S[i];
	pred[p] = NIL;
	if(label[p]==0)
	  cm = 1; //1;
	else
	  cm = 2; //2;
	  
	if (Q->color[p] != BLACK) {
	  glip::Heap32fi_lex::Update_MinPolicy(Q, p, 0.0, cm);
	} else {
	  Q->color[p] = WHITE;
	  Q->cost1[p] = 0.0;
	  maxorder[p] = cm;
	  glip::Heap32fi_lex::Insert_MinPolicy(Q, p);
	}
	Q->color[p] = GRAY;
      }

      //----------
      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);
	
	if(label[p]==0){
	  g = transpose; cm = 1;
	  if(H[p] == 0){
	    (*obj_area)--;
	    H[p] = T_area;
	  }
	  if(E->color[p] == GRAY)
	    glip::Heap32::Delete_MaxPolicy(E, p);
	} //1;
	else{   	 g = graph;     cm = 2; } //2;
	//-------------------------------------------------------
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  w = g->nodes[p].Warcs[i];
	  
	  Wm = MAX(w, Q->cost1[p]);
	    
	  if(Wm != w)
	    maxord = maxorder[p];
	  else if( Q->cost1[p] == Wm )
	    maxord = maxorder[p] + 2;
	  else
	    maxord = cm;

	  seenfirst_pq = seenfirst_qp = false;
	  if(pred[q] != NIL){
	    aux = (w > Q->cost1[p] && Q->cost1[q] > Q->cost1[pred[q]] && w == Q->cost1[q]);
	    seenfirst_pq = aux && (Q->cost1[p] < Q->cost1[pred[q]] ||
				   (Q->cost1[p] == Q->cost1[pred[q]] &&
				    maxorder[p] < maxorder[pred[q]]));
	    seenfirst_qp = aux && (Q->cost1[p] > Q->cost1[pred[q]] ||
				   (Q->cost1[p] == Q->cost1[pred[q]] &&
				    maxorder[p] > maxorder[pred[q]]));
	  }
	  
	  if(Wm < Q->cost1[q] ||
	     (Wm == Q->cost1[q] && ((maxord < maxorder[q] && !seenfirst_qp) || seenfirst_pq))){
	    label[q] = label[p];
	    pred[q] = p;
	    if(Q->color[q] == GRAY)
	      glip::Heap32fi_lex::Delete_MinPolicy(Q, q);
	    Q->color[q] = WHITE;
	    Q->cost1[q] = Wm;
	    maxorder[q] = maxord;
	    glip::Heap32fi_lex::Insert_MinPolicy(Q, q);
	  }
	  else if((p == pred[q]) &&
		  ((label[q] != label[p])||
		   (Wm != Q->cost1[q]) ||
		   (maxord > maxorder[q]))){

	    if(label[q] == label[p])
	      printf("Ops...\n");
	    
	    if(T == NULL)
	      T = (int *)malloc((graph->nnodes+1)*sizeof(int));
	      
	    Removal[0] = 1; Removal[1] = q;
	    DOIFT_TreeRemoval_v2_v3(Removal, graph, Q, pred,
				    maxorder, T);
	    break;
	  }
	  else if(label[p] == 0 && label[q] == 1){
	    if(E->color[q] != GRAY)
	      glip::Heap32::Insert_MaxPolicy(E, q);
	  }
	  else if(label[p] == 1 && label[q] == 0 && H[q] > 0){
	    if(E->color[p] != GRAY)
	      glip::Heap32::Insert_MaxPolicy(E, p);
	  }
	}
      
	//-------------------------------------------------------
      }
      //----------
      if(T != NULL)
	free(T);
    }
    



    int *OIFT_area_hierarchy(sGraph *graph,
			     sGraph *transpose,
			     int *S0,
			     int *S1){
			     //glip::sImage32 *SP){
      /*
      char filename[512];
      glip::sImage32 *temp = NULL;
      int sp,it = 0;
      temp = glip::Image32::Create(SP);
      */
      //-------------------------------------
      int *H = NULL, *Rset = NULL, *S = NULL;
      glip::sHeap32fi_lex *Q = NULL;
      glip::sHeap32 *E = NULL;
      float *cost1=NULL;
      int *energy=NULL, *maxorder=NULL, *label=NULL, *pred=NULL;
      int p,p_max,i, obj_area;
      H = glip::AllocIntArray(graph->nnodes);
      S = glip::AllocIntArray(graph->nnodes+1);
      cost1  = glip::AllocFloatArray(graph->nnodes);
      energy   = glip::AllocIntArray(graph->nnodes);
      maxorder = glip::AllocIntArray(graph->nnodes);
      label    = glip::AllocIntArray(graph->nnodes);
      pred     = glip::AllocIntArray(graph->nnodes);
      for(p=0; p < graph->nnodes; p++){
	H[p] = 0;
	cost1[p] = FLT_MAX;
	maxorder[p] = INT_MAX;
	pred[p] = NIL;
	label[p] = NIL;
      }
      Q = glip::Heap32fi_lex::Create(graph->nnodes,
				    cost1, maxorder);
      
      for(i=1; i<=S1[0]; i++){
	p = S1[i];
	label[p] = 1;
      }
      glip::ift::IFT_fmax(graph, S1, label, energy);
      //------------------------------------------
      /*
      for(p = 0; p < temp->n; p++){
	sp = SP->data[p];
	temp->data[p] = energy[sp];
      }
      glip::Image32::Write(temp, (char*)"energy.pgm");
      */
      //------------------------------------------
      for(i=1; i<=S1[0]; i++){
	p = S1[i];
	energy[p] = -1;
      }
      E = glip::Heap32::Create(graph->nnodes, energy);
      
      for(p=0; p < graph->nnodes; p++){
	label[p] = NIL;
      }
      for(i=1; i<=S1[0]; i++){
	p = S1[i];
	label[p] = 1;
      }
      for(i=1; i<=S0[0]; i++){
	p = S0[i];
	label[p] = 0;
      }

      S[0] = 0;
      for(i=1; i<=S1[0]; i++){
	S[0]++;
	S[S[0]] = S1[i];
      }
      for(i=1; i<=S0[0]; i++){
	S[0]++;
	S[S[0]] = S0[i];
      }
      //------------------------------
      //printf("#Internal seeds: %d\n", S1[0]);
      //------------------------------
      obj_area = graph->nnodes;
      while(obj_area > S1[0]){
      
	glip::ift::DOIFT_v3_hierarchy(graph,
				      transpose,
				      S,
				      Rset,
				      label,
				      Q,
				      pred,
				      maxorder,
				      E,
				      H,
				      &obj_area);
	//------------------------------------------
	/*
	it++;
	printf("it: %d, obj_area: %d\n", it, obj_area);
	for(p = 0; p < temp->n; p++){
	  sp = SP->data[p];
	  temp->data[p] = (E->color[sp] == GRAY)*255;
	}
	sprintf(filename, "objboundary_%02d.pgm", it);
	glip::Image32::Write(temp, filename);
	for(p = 0; p < temp->n; p++){
	  sp = SP->data[p];
	  temp->data[p] = label[sp]*255;
	}
	sprintf(filename, "label_%02d.pgm", it);
	glip::Image32::Write(temp, filename);
	*/
	//------------------------------------------
	
	Heap32::Remove_MaxPolicy(E, &p_max);
	
	S[0] = 1;
	S[1] = p_max;
	label[p_max] = 0;

      }

      for(i=1; i<=S1[0]; i++){
	p = S1[i];
	H[p] = S1[0];
      }
      
      glip::Heap32::Destroy(&E);      
      glip::Heap32fi_lex::Destroy(&Q);
      if(S  != NULL) glip::FreeIntArray(&S);
      if(energy   != NULL) glip::FreeIntArray(&energy);
      if(label    != NULL) glip::FreeIntArray(&label);
      if(pred     != NULL) glip::FreeIntArray(&pred);
      if(maxorder != NULL) glip::FreeIntArray(&maxorder);
      if(cost1 != NULL) glip::FreeFloatArray(&cost1);
      //--------------------------------------------
      //glip::Image32::Destroy(&temp);
      //--------------------------------------------
      return H;
    }


    

  } /*end ift namespace*/
} /*end glip namespace*/


