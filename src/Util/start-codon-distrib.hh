//  A. L. Delcher
//
//  File:  start-codon-distrib.hh
//
//  Last Modified:  25 July 2005
//
//  Declarations for  start-codon-distrib.cc



#ifndef  __START_CODON_DISTRIB_HH_INCLUDED
#define  __START_CODON_DISTRIB_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "gene.hh"


struct  Count_Entry_t
  {
   char  codon [4];
   int  count;
  };


static bool  Count_Entry_Cmp
    (Count_Entry_t const & a, Count_Entry_t const & b);



static void  Incr
    (vector <Count_Entry_t> & entry, const char s [4]);
static void  Parse_Command_Line
    (int argc, char * argv []);
static long int  Seq_Sub
    (long int sub, long int seq_len);
static void  Usage
    (void);

#endif
