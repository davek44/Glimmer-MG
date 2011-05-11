//  A. L. Delcher
//
//  File:  extract.hh
//
//  Last Modified:  13 May 2005
//
//  Declarations for  extract.cc



#ifndef  __EXTRACT_HH_INCLUDED
#define  __EXTRACT_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "gene.hh"


static void  Output_Subsequence
     (const string & seq, long int i, long int len, int incr,
      const char * tag, long int start, long int end);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Usage
    (void);

#endif
