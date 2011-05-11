//  A. L. Delcher
//
//  File:  multi-extract.hh
//
//  Last Modified:  Tue Aug  9 09:30:18 EDT 2005
//
//  Declarations for  multi-extract.cc



#ifndef  __MULTI_EXTRACT_HH_INCLUDED
#define  __MULTI_EXTRACT_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "gene.hh"


struct  Coord_Info_t
  {
   char  * id, * tag;
   long int  start, end;
   int  dir;
  };


static bool  By_Tag
    (Coord_Info_t const & a, Coord_Info_t const & b)
  {
   return  (strcmp (a . tag, b . tag) < 0);
  }



static void  Find_Matches
    (char * p, const vector <Coord_Info_t> & list, int & sub, int & num);
static void  Output_Subsequence
    (const string & seq, long int i, long int len, int incr,
     const char * id, const char * tag, long int start, long int end);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Usage
    (void);

#endif
