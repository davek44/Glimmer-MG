//    Programmer:  Arthur L. Delcher
//          File:  build-fixed.hh
//  Last Updated:  Fri Jun  4 16:33:12 EDT 2004
//
//  Declarations for  build-fixed.cc


#ifndef _BUILD_FIXED_HH
#define _BUILD_FIXED_HH


#include  "icm.hh"


static void  Parse_Command_Line
    (int argc, char * argv []);
static int  Read_String
    (FILE * fp, char * & s, long int & s_size, char * & tag,
     long int & tag_size);
static void  Read_Training_Data
    (FILE * fp);
static void  Usage
    (char * command);



#endif
