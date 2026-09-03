
#ifndef _GLIP_CONFIGFILE_H_
#define _GLIP_CONFIGFILE_H_

#include "glip_attributelist.h"

namespace glip{
  namespace ConfigFile{

    struct sConfigFile{
      sAttributeList **lines;
      int nlines;
    };
    
    sConfigFile *Create();
    void Destroy(sConfigFile **cf);
    
    sConfigFile *Read(char filename[]);

  } //end ConfigFile namespace

  typedef ConfigFile::sConfigFile sConfigFile;

} //end glip namespace

 
#endif

