//  A. L. Delcher
//
//  File:  window-acgt.hh
//
//  Last Modified:  21 June 2005
//
//  Declarations for  window-acgt.cc



#ifndef  __WINDOW_ACGT_HH_INCLUDED
#define  __WINDOW_ACGT_HH_INCLUDED


#include  "delcher.hh"
#include  "gene.hh"


static void  Finish
    (int win_pos, int win_size, int win_next, int count [5],
     const vector <char> & window, int win_sub);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Print_Line
    (int win_pos, int win_size, int count [5]);
int static  Subscript
    (char ch);
static void  Usage
    (void);

#endif
