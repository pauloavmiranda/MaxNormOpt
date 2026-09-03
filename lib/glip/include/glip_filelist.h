// version 00.00.04

#ifndef _GLIP_FILELIST_H_
#define _GLIP_FILELIST_H_

#include "glip_common.h"
#include "glip_arraylist.h"
#include "glip_string.h"

extern "C" {
#include <glob.h>
}


namespace glip{
  namespace FileList{

    struct sFileList {
      sArrayList *A;
      int n;   //Number of files added.
    };

    sFileList *Create(int cap);
    void       Destroy(sFileList **L);
    
    sFileList *Read(char *filename);
    void       Write(sFileList *L,
		     char *filename);

    void      AddFile(sFileList *L, char *file);
    char     *GetFile(sFileList *L, int index);
    bool      HasFile(sFileList *L, char *file);
    
    
    void AddFilesInDir(sFileList *L,
		       char *dir);
    void AddFilesInDirRec(sFileList *L,
			  char *dir);

    //It shuffles the files at random.
    //Should use "void srand(unsigned int seed);" before calling.
    void  Randomize(sFileList *L);
    
    void  Resize(sFileList *L, int n);
    
    //Trims the capacity of this FileList instance 
    //to be the list's current size.
    void  Trim2Size(sFileList *L);

    void  DeleteFilesInFileList(sFileList *L);
    
    bool      FileExists(char *file);
    void      RemoveFileDirectory(char *file);
    void      RemoveFileExtension(char *file);
    void      MergeRelativePath(char *dir, char *rel);

  } /*end FileList namespace*/

  typedef FileList::sFileList sFileList;

} /*end glip namespace*/

#endif

