
#ifndef _GLIP_STACK_H
#define _GLIP_STACK_H 1

#include "glip_common.h"

namespace glip{
  namespace Stack{

    struct sStack {
      int *data;
      int top;
      int n;
    };
    
    sStack *Create(int n);
    void    Destroy(sStack **S);

    /**
     * \brief Returns the size of the data structure in bytes.
     */
    size_t GetSizeOf(sStack *S);
    
    //---------inline definitions------------------
    inline void    Push(sStack *S, int p){
      S->data[++(S->top)] = p;
    }
    
    inline int    Pop(sStack *S){
      //if(S->top == -1) return -1;
      return(S->data[S->top--]);
    }
    
    inline void Clear(sStack *S){
      S->top = -1;
    }
    
    inline int IsEmpty(sStack *S){
      return(S->top == -1);
    }

  } //end Stack namespace

  typedef Stack::sStack sStack;

} //end glip namespace

#endif
