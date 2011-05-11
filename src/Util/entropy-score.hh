//  A. L. Delcher
//
//  File:  entropy-score.hh
//
//  Last Modified:  29 July 2005
//
//  Declarations for  entropy-score.cc



#ifndef  __ENTROPY_SCORE_HH_INCLUDED
#define  __ENTROPY_SCORE_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "gene.hh"


static double  Entropy_Distance_Ratio
    (const string & buff);
static int  On_Seq
    (int i, int seq_len);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Read_Entropy_Profiles
    (const char * fn, bool & errflg);
static void  Usage
    (void);

#endif
