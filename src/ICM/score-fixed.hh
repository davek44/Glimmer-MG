//    Programmer:  Arthur L. Delcher
//          File:  score-fixed.hh
//  Last Updated:  Mon Jun  7 10:36:48 EDT 2004
//
//  Declarations for  score-fixed.cc


#ifndef _SCORE_FIXED_HH
#define _SCORE_FIXED_HH


#include  "icm.hh"


static void  Parse_Command_Line
    (int argc, char * argv []);
static int  Read_String
    (FILE * fp, char * & s, long int & s_size, char * & tag,
     long int & tag_size);
static void  Usage
    (char * command);


#endif
