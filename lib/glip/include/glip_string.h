
#ifndef _GLIP_STRING_H_
#define _GLIP_STRING_H_

#include "glip_common.h"

namespace glip{
  namespace String{

    //Removes the leading and trailing white space.
    void      Trim(char *str);
    void      SubString(char *str,
			int beginIndex,
			int endIndex);
    void      ReplaceCharacter(char *str,
			       char old_c,
			       char new_c);

  } //end String namespace
} //end glip namespace


#endif



