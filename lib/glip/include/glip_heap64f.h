#ifndef _GLIP_HEAP64F_H_
#define _GLIP_HEAP64F_H_

#include "glip_common.h"
#include "glip_gpqueue_by_Falcao.h"

#include "glip_heap.h"

namespace glip{
  namespace Heap64f{

    struct sHeap64f {
      double *cost;
      char *color;
      int *pixel;
      int *pos;
      int last;
      int n;
    };

    /* Auxiliary Functions */

    char IsFull(sHeap64f *H);
    char IsEmpty(sHeap64f *H);
    sHeap64f *Create(int n, double *cost);
    void Destroy(sHeap64f **H);

    void Insert_MaxPolicy(sHeap64f *H, int pixel);
    void Insert_MinPolicy(sHeap64f *H, int pixel);
    
    void Remove_MaxPolicy(sHeap64f *H, int *pixel);
    void Remove_MinPolicy(sHeap64f *H, int *pixel);

    void Get_MaxPolicy(sHeap64f *H, int *pixel);
    void Get_MinPolicy(sHeap64f *H, int *pixel);
    
    void Update_MaxPolicy(sHeap64f *H, int p, double value);
    void Update_MinPolicy(sHeap64f *H, int p, double value);

    void GoUp_MaxPolicy(sHeap64f *H, int i);
    void GoUp_MinPolicy(sHeap64f *H, int i);

    void GoDown_MaxPolicy(sHeap64f *H, int i);
    void GoDown_MinPolicy(sHeap64f *H, int i);

    void Reset(sHeap64f *H);

    void Delete_MaxPolicy(sHeap64f *H, int pixel);
    void Delete_MinPolicy(sHeap64f *H, int pixel);
    

  } //end Heap64f namespace

  typedef Heap64f::sHeap64f sHeap64f;

} //end glip namespace

#endif

