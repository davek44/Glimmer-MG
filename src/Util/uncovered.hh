//  A. L. Delcher
//
//  File:  uncovered.hh
//
//  Last Modified:  10 June 2005
//
//  Declarations for  uncovered.cc



#ifndef  __UNCOVERED_HH_INCLUDED
#define  __UNCOVERED_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "gene.hh"


struct  Region_t
  {
   long int  lo, hi;  // Range of region in "between" coordinates
  };


static bool  Region_Cmp
    (Region_t const & a, Region_t const & b)
  { return  (a . lo < b . lo); }



static void  Coalesce_Regions
    (vector <Region_t> & list);
static void  Output_Subsequence
     (const string & seq, long int i, long int len,
      const char * tag, long int start, long int end);
static void  Output_Uncovered
    (const string & seq, long int seq_len, const vector <Region_t> & list);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Usage
    (void);

#endif
