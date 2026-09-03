
/*-------------------------------------------------------------*/
/*--- Public header file for the library.                   ---*/
/*---                                               bzlib.h ---*/
/*-------------------------------------------------------------*/

/* ------------------------------------------------------------------
   This file is part of bzip2/libbzip2, a program and library for
   lossless, block-sorting data compression.

   bzip2/libbzip2 version 1.0.4 of 20 December 2006
   Copyright (C) 1996-2006 Julian Seward <jseward@bzip.org>

   Please read the WARNING, DISCLAIMER and PATENTS sections in the 
   README file.

   This program is released under the terms of the license contained
   in the file LICENSE.
   ------------------------------------------------------------------ */


#ifndef _glip_BZLIB_H
#define _glip_BZLIB_H

#ifdef __cplusplus
extern "C" {
#endif

#define glip_BZ_RUN               0
#define glip_BZ_FLUSH             1
#define glip_BZ_FINISH            2

#define glip_BZ_OK                0
#define glip_BZ_RUN_OK            1
#define glip_BZ_FLUSH_OK          2
#define glip_BZ_FINISH_OK         3
#define glip_BZ_STREAM_END        4
#define glip_BZ_SEQUENCE_ERROR    (-1)
#define glip_BZ_PARAM_ERROR       (-2)
#define glip_BZ_MEM_ERROR         (-3)
#define glip_BZ_DATA_ERROR        (-4)
#define glip_BZ_DATA_ERROR_MAGIC  (-5)
#define glip_BZ_IO_ERROR          (-6)
#define glip_BZ_UNEXPECTED_EOF    (-7)
#define glip_BZ_OUTBUFF_FULL      (-8)
#define glip_BZ_CONFIG_ERROR      (-9)

typedef 
   struct {
      char *next_in;
      unsigned int avail_in;
      unsigned int total_in_lo32;
      unsigned int total_in_hi32;

      char *next_out;
      unsigned int avail_out;
      unsigned int total_out_lo32;
      unsigned int total_out_hi32;

      void *state;

      void *(*bzalloc)(void *,int,int);
      void (*bzfree)(void *,void *);
      void *opaque;
   } 
   bz_stream;


#ifndef glip_BZ_IMPORT
#define glip_BZ_EXPORT
#endif

#ifndef glip_BZ_NO_STDIO
/* Need a definitition for FILE */
#include <stdio.h>
#endif

//#ifdef _WIN32
//#   include <windows.h>
//#   ifdef small
//      /* windows.h define small to char */
//#      undef small
//#   endif
//#   ifdef glip_BZ_EXPORT
//#   define glip_BZ_API(func) WINAPI func
//#   define glip_BZ_EXTERN extern
//#   else
//   /* import windows dll dynamically */
//#   define glip_BZ_API(func) (WINAPI * func)
//#   define glip_BZ_EXTERN
//#   endif
//#else
//#   define glip_BZ_API(func) func
//#   define glip_BZ_EXTERN extern
//#endif

#   define glip_BZ_API(func) func
#   define glip_BZ_EXTERN extern


/*-- Core (low-level) library functions --*/

glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzCompressInit) (
      bz_stream* strm, 
      int        blockSize100k, 
      int        verbosity, 
      int        workFactor 
   );

glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzCompress) ( 
      bz_stream* strm, 
      int action 
   );

glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzCompressEnd) ( 
      bz_stream* strm 
   );

glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzDecompressInit) ( 
      bz_stream *strm, 
      int       verbosity, 
	  int small
   );

glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzDecompress) ( 
      bz_stream* strm 
   );

glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzDecompressEnd) ( 
      bz_stream *strm 
   );



/*-- High(er) level library functions --*/

#ifndef glip_BZ_NO_STDIO
#define glip_BZ_MAX_UNUSED 5000

typedef void glip_BZFILE;

glip_BZ_EXTERN glip_BZFILE* glip_BZ_API(glip_BZ2_bzReadOpen) ( 
      int*  bzerror,   
      FILE* f, 
      int   verbosity, 
      int   small,
      void* unused,    
	  int nUnused
   );

glip_BZ_EXTERN void glip_BZ_API(glip_BZ2_bzReadClose) ( 
      int*    bzerror, 
      glip_BZFILE* b 
   );

glip_BZ_EXTERN void glip_BZ_API(glip_BZ2_bzReadGetUnused) ( 
      int*    bzerror, 
      glip_BZFILE* b, 
      void**  unused,  
      int*    nUnused 
   );

glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzRead) ( 
      int*    bzerror, 
      glip_BZFILE* b, 
      void*   buf, 
      int     len 
   );

glip_BZ_EXTERN glip_BZFILE* glip_BZ_API(glip_BZ2_bzWriteOpen) ( 
      int*  bzerror,      
      FILE* f, 
      int   blockSize100k, 
      int   verbosity, 
      int   workFactor 
   );

glip_BZ_EXTERN void glip_BZ_API(glip_BZ2_bzWrite) ( 
      int*    bzerror, 
      glip_BZFILE* b, 
      void*   buf, 
      int     len 
   );

glip_BZ_EXTERN void glip_BZ_API(glip_BZ2_bzWriteClose) ( 
      int*          bzerror, 
      glip_BZFILE*       b, 
      int           abandon, 
      unsigned int* nbytes_in, 
      unsigned int* nbytes_out 
   );

glip_BZ_EXTERN void glip_BZ_API(glip_BZ2_bzWriteClose64) ( 
      int*          bzerror, 
      glip_BZFILE*       b, 
      int           abandon, 
      unsigned int* nbytes_in_lo32, 
      unsigned int* nbytes_in_hi32, 
      unsigned int* nbytes_out_lo32, 
      unsigned int* nbytes_out_hi32
   );
#endif


/*-- Utility functions --*/

glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzBuffToBuffCompress) ( 
      char*         dest, 
      unsigned int* destLen,
      char*         source, 
      unsigned int  sourceLen,
      int           blockSize100k, 
      int           verbosity, 
      int           workFactor 
   );

glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzBuffToBuffDecompress) ( 
      char*         dest, 
      unsigned int* destLen,
      char*         source, 
      unsigned int  sourceLen,
      int           small, 
	  int verbosity
   );

/*--
   Code contributed by Yoshioka Tsuneo (tsuneo@rr.iij4u.or.jp)
   to support better zlib compatibility.
   This code is not _officially_ part of libbzip2 (yet);
   I haven't tested it, documented it, or considered the
   threading-safeness of it.
   If this code breaks, please contact both Yoshioka and me.
--*/

glip_BZ_EXTERN const char * glip_BZ_API(glip_BZ2_bzlibVersion) (
      void
   );

#ifndef glip_BZ_NO_STDIO
glip_BZ_EXTERN glip_BZFILE * glip_BZ_API(glip_BZ2_bzopen) (
      const char *path,
      const char *mode
   );

glip_BZ_EXTERN glip_BZFILE * glip_BZ_API(glip_BZ2_bzdopen) (
      int        fd,
      const char *mode
   );
         
glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzread) (
      glip_BZFILE* b, 
      void* buf, 
      int len 
   );

glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzwrite) (
      glip_BZFILE* b, 
      void*   buf, 
      int     len 
   );

glip_BZ_EXTERN int glip_BZ_API(glip_BZ2_bzflush) (
      glip_BZFILE* b
   );

glip_BZ_EXTERN void glip_BZ_API(glip_BZ2_bzclose) (
      glip_BZFILE* b
   );

glip_BZ_EXTERN const char * glip_BZ_API(glip_BZ2_bzerror) (
      glip_BZFILE *b, 
      int    *errnum
   );
#endif

#ifdef __cplusplus
}
#endif

#endif

/*-------------------------------------------------------------*/
/*--- end                                           bzlib.h ---*/
/*-------------------------------------------------------------*/
