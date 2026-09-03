
#include "glip_stack.h"

namespace glip{
  namespace Stack{

    sStack *Create(int n) {
      sStack *S;
      S = (sStack *) malloc(sizeof(sStack));
      if(S==NULL) glip::Error((char *)MSG1,
			     (char *)"Stack::Create");
      S->n   = n;
      S->top = -1;
      S->data = glip::AllocIntArray(n);
      return S;
    }

    void    Destroy(sStack **S) {
      sStack *aux = *S;
      if(aux) {
	if(aux->data) glip::FreeIntArray(&aux->data);
	free(aux);
	*S = NULL;
      }
    }


    size_t GetSizeOf(sStack *S){
      size_t size = 0;
      if(S == NULL)
	return 0;
      size = sizeof(sStack);
      size += S->n*sizeof(int);
      //printf("%zu\n", size);
      return size;
    }

    
  } //end Stack namespace
} //end glip namespace

