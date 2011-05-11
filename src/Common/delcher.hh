//  A. L. Delcher
//
//  File:  delcher.hh
//
//  Last Modified:  23 October 2003
//
//  Common generic routines declarations


#ifndef  __DELCHER_HH_INCLUDED
#define  __DELCHER_HH_INCLUDED

#define  ALLOW_LONG_OPTIONS  1

#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <ctype.h>
#include  <float.h>
#include  <time.h>
#include  <assert.h>
#include  <errno.h>
#include  <getopt.h>
#include  <limits.h>
#include  <algorithm>
#include  <string>
#include  <new>
#include  <cstdlib>
#include  <iostream>
#include  <iomanip>
#include  <fstream>
#include  <vector>
#include  <string>
#include  <cstring>

#include  "exceptions.hh"


using namespace  std;


#define  TRUE  1
#define  FALSE  0
#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif


#define  SAFE_CALLOC(x,y)  Safe_calloc (x, y, __FILE__, __LINE__)
#define  SAFE_MALLOC(x)  Safe_malloc (x, __FILE__, __LINE__)
#define  SAFE_REALLOC(x,y)  Safe_realloc (x, y, __FILE__, __LINE__)


const int  MAX_ERROR_MSG_LEN = 1000;
  // Length of longest possible error message


extern char  Clean_Exit_Msg_Line [MAX_ERROR_MSG_LEN];
  // String to write error messages to before exiting
extern int  Verbose;
  // Flag to determine level of debugging output


const char *  Commatize
    (long int  n);
void  Clean_Exit
    (const char * msg, const char * src_fname = NULL, size_t line_num = 0);
FILE *  File_Open
    (const string & fname, const string & mode, const char * src_fname = NULL,
     size_t line_num = 0);
char  First_Non_Blank
    (const char * s);
int  Int_Power
    (int a, int b);
void  Make_Lower_Case
    (char * s);
void  Make_Upper_Case
    (char * s);
const char *  Num_Or_Max
    (int n, int mx = INT_MAX);
double  Percent
    (double a, double b);
const char *  Printable
    (bool b);
const char *  Printable
    (char * s);
double  Pseudo_Normal
    (void);
double  Ratio
    (double a, double b);
void  Reverse_String
    (char * s);
void  Reverse_String
    (string & s);
void *  Safe_calloc
    (size_t n, size_t len, const char * src_fname = NULL,
     size_t line_num = 0);
void *  Safe_malloc
    (size_t len, const char * src_fname = NULL, size_t line_num = 0);
void *  Safe_realloc
    (void * q, size_t len, const char * src_fname = NULL,
     size_t line_num = 0);
char *  Strip_Trailing
    (char * s, char ch);


template <class DT>  void  Incr_Limited
    (DT & A, DT limit);
template <class DT>  DT  Max
    (DT, DT);
template <class DT>  DT  Min
    (DT, DT);
template <class DT>  void  Swap
    (DT &, DT &);



template <class DT>  void  Incr_Limited
    (DT & A, DT limit)

// Increment  A  by 1, but only if it's less than  limit .

  {
   if  (A < limit)
       A ++;

   return;
  }



template <class DT>  DT  Max
    (DT A, DT B)

// Return the larger of  A  and  B .

  {
   if  (A > B)
       return  A;
     else
       return  B;
  }



template <class DT>  DT  Min
    (DT A, DT B)

// Return the smaller of  A  and  B .

  {
   if  (A < B)
       return  A;
     else
       return  B;
  }



template <class DT>  void  Reverse
    (vector <DT> & v)

// Reverse the order of entries in  v .

  {
   DT  s;
   int  i, j, n;

   n = v . size ();
   for  (i = 0, j = n - 1;  i < j;  i ++, j --)
     {
      s = v [i];
      v [i] = v [j];
      v [j] = s;
     }

   return;
  }



template <class DT>  void  Swap
    (DT & A, DT & B)

// Swap the values in  A  and  B .

  {
   DT  Save;

   Save = A;
   A = B;
   B = Save;

   return;
  }



#endif
