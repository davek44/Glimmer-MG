//  Programmers:   Arthur L. Delcher
//                 Doug Harmon
//  File:          build-icm.hh
//  Last Updated:  Mon May  3 08:07:22 EDT 2004
//
//  Declarations for  build-icm.cc


#ifndef _BUILD_ICM_HH
#define _BUILD_ICM_HH


#include  "icm.hh"


// The 16 possible base pair sequences
const int  AA_PAIR = 0;
const int  AC_PAIR = 1;
const int  AG_PAIR = 2;
const int  AT_PAIR = 3;
const int  CA_PAIR = 4;
const int  CC_PAIR = 5;
const int  CG_PAIR = 6;
const int  CT_PAIR = 7;
const int  GA_PAIR = 8;
const int  GC_PAIR = 9;
const int  GG_PAIR = 10;
const int  GT_PAIR = 11;
const int  TA_PAIR = 12;
const int  TC_PAIR = 13;
const int  TG_PAIR = 14;
const int  TT_PAIR = 15;


// The 4 possible bases
const int  A_BASE = 0;
const int  C_BASE = 1;
const int  G_BASE = 2;
const int  T_BASE = 3;


// Function prototypes

static void  Parse_Command_Line
    (int argc, char * argv []);
static int  Read_String
    (FILE * fp, char * & s, long int & s_size, char * & tag, long int & tag_size);
static void  Read_Training_Data
    (FILE * fp);
static void  Set_Stop_Codons
    (void);
static void  Usage
    (char * command);

#endif
