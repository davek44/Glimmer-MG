//  A. L. Delcher
//
//  File:  delcher.cc
//
//  Last Modified:  23 October 2003
//
//  Common generic routines.


#include  "delcher.hh"

const int  COMMATIZE_BUFF_LEN = 50;
  // Length of buffer for creating string with commas


char  Clean_Exit_Msg_Line [MAX_ERROR_MSG_LEN];
  // String to write error messages to before exiting
int  Global_Debug_Flag = 0;
  // Used to set different debugging levels
int  Verbose = 0;
  // Flag to determine level of debugging output


const char *  Commatize
    (long int  n)

//  Return a string representing the value of  n  with commas
//  every three digits.

  {
   static char  buff [COMMATIZE_BUFF_LEN];
   bool  is_negative = false;
   int  i, ct;

   buff [COMMATIZE_BUFF_LEN - 1] = '\0';

   if  (n == 0)
       {
        buff [COMMATIZE_BUFF_LEN - 2] = '0';
        return  buff + 48;
       }

   i = COMMATIZE_BUFF_LEN - 2;
   if  (n < 0)
       {
        is_negative = true;
        n *= -1;
       }

   for  (ct = 0;  n > 0;  ct ++)
     {
      if  (ct == 3)
          {
           buff [i --] = ',';
           ct = 0;
          }
      buff [i --] = char ('0' + n % 10);
      n /= 10;
     }

   if  (is_negative)
       buff [i --] = '-';

   return  buff + i + 1;
  }



void  Clean_Exit
    (const char * msg, const char * src_fname, size_t line_num)

//  Write string  msg  to  stderr  and also a line indicating
//  the error happen in source file  src_fname  at line  line_num
//  if they are not  NULL  and  0  respectively.
//  Then exit with an error condition.

  {
   fprintf (stderr, "%s\n", msg);
   if  (src_fname != NULL)
       fprintf (stderr, "  in file  %s", src_fname);
   if  (line_num != 0)
       fprintf (stderr, "  at line  %lu", (long unsigned) (line_num));
   fprintf (stderr, "  errno = %d\n", errno);

   exit (EXIT_FAILURE);
  }



FILE *  File_Open
    (const string & fname, const string & mode, const char * src_fname,
     size_t line_num)

//  Open  fname  in  mode  and return a pointer to its control
//  block.  If fail, print a message and exit, assuming the call came from
//  source file  src_fname  at line  line_num .

  {
   FILE  *  fp;

   fp = fopen (fname . c_str (), mode . c_str ());
   if  (fp == NULL)
       {
        sprintf (Clean_Exit_Msg_Line,
                 "ERROR:  Could not open file  %s", fname . c_str ());
        Clean_Exit (Clean_Exit_Msg_Line, src_fname, line_num);
       }

   return  fp;
  }



char  First_Non_Blank
    (const char * s)

//  Return the first non-white-space character in  s .
//  If none, return  ' ' .

  {
   const char  * p;

   for  (p = s;  * p != '\0';  p ++)
     if  (! isspace (* p))
         return  * p;

   return  ' ';
  }



void  Make_Lower_Case
    (char * s)

//  Convert all letters in  s  to lower-case.

  {
   char  * p;

   for  (p = s;  * p != '\0';  p ++)
     * p = tolower (* p);
  }



void  Make_Upper_Case
    (char * s)

//  Convert all letters in  s  to upper-case.

  {
   char  * p;

   for  (p = s;  * p != '\0';  p ++)
     * p = toupper (* p);
  }



int  Int_Power
    (int a, int b)

//  Return  a  raised to the power  b .  Both are assumed
//  to be non-negative;

  {
   int  result = 1;
   int  p = a;

   while  (b > 0)
     {
      if  (b & 1)
          result *= p;
      p = p * p;
      b >>= 1;
     }

   return  result;
  }



const char *  Num_Or_Max
    (int n, int mx)

//  Return a string representation of  n  or else "max"
//  if  n >= mx .

  {
   static char  buff [20];

   if  (n >= mx)
       return  "max";

   sprintf (buff, "%d", n);
   return  buff;
  }



double  Percent
    (double a, double b)

//  Return  a / b  as a percentage.  Return  0.0  if  b = 0.0 .

  {
   if  (b == 0.0)
       return  0.0;

   return  100.0 * (a / b);
  }



const char *  Printable
    (bool b)

//  Return a string representing the value of  b .

  {
   if  (b)
       return  "true";
     else
       return  "false";
  }



const char *  Printable
    (char * s)

//  Return "none" if  s  is  NULL ; otherwise, return  s  itself.

  {
   if  (s == NULL)
       return  "none";
     else
       return  s;
  }



double  Pseudo_Normal
    (void)

//  Return a simple approximation to a standard normally distributed value,
//  i.e., mean = 0.0 and s.d. = 1.0

  {
   double  sum = 0.0;
   int  i;

   for  (i = 0;  i < 12;  i ++)
     sum += drand48 ();

   return  sum - 6.0;
  }



double  Ratio
    (double a, double b)

//  Return  a / b , or  0.0  if  b = 0.0 .

  {
   if  (b == 0.0)
       return  0.0;

   return  a / b;
  }



void  Reverse_String
    (char * s)

//  Reverse the order of characters in string  s .

  {
   int  i, j, n;

   n = strlen (s);
   for  (i = 0, j = n - 1;  i < j;  i ++, j --)
     {
      char  ch;

      ch = s [i];
      s [i] = s [j];
      s [j] = ch;
     }

   return;
  }



void  Reverse_String
    (string & s)

//  Reverse the order of characters in string  s .

  {
   int  i, j, n;

   n = s . length ();
   for  (i = 0, j = n - 1;  i < j;  i ++, j --)
     {
      char  ch;

      ch = s [i];
      s [i] = s [j];
      s [j] = ch;
     }

   return;
  }



void *  Safe_calloc
    (size_t n, size_t len, const char * src_fname, size_t line_num)

//  Allocate and return a pointer to enough memory to hold an
//  array with  n  entries of  len  bytes each.  All memory is
//  cleared to 0.  If fail, print a message and exit, assuming the
//  call came from  source file  src_fname  at line  line_num .

  {
   void  * p;

   p = calloc (n, len);
   if  (p == NULL)
       {
        sprintf (Clean_Exit_Msg_Line,
                 "ERROR:  calloc failed  %lu x %lu",
                 (long unsigned) (n), (long unsigned) (len));
        Clean_Exit (Clean_Exit_Msg_Line, src_fname, line_num);
       }

   return  p;
  }



void *  Safe_malloc
    (size_t len, const char * src_fname, size_t line_num)

//  Allocate and return a pointer to  len  bytes of memory.
//  If fail, print a message and exit, assuming the call came from
//  source file  src_fname  at line  line_num .

  {
   void  * p;

   p = malloc (len);
   if  (p == NULL)
       {
        sprintf (Clean_Exit_Msg_Line,
                 "ERROR:  malloc failed  %lu  bytes",
                 (long unsigned) (len));
        Clean_Exit (Clean_Exit_Msg_Line, src_fname, line_num);
       }

   return  p;
  }



void *  Safe_realloc
    (void * q, size_t len, const char * src_fname, size_t line_num)

//  Reallocate memory for  q  to  len  bytes and return a
//  pointer to the new memory.  If fail, print a message and exit,
//  assuming the call came from source file  src_fname  at line  line_num .

  {
   void  * p;

   p = realloc (q, len);
   if  (p == NULL)
       {
        sprintf (Clean_Exit_Msg_Line,
                 "ERROR:  realloc failed  %lu  bytes",
                 (long unsigned) (len));
        Clean_Exit (Clean_Exit_Msg_Line, src_fname, line_num);
       }

   return  p;
  }



char *  Strip_Trailing
    (char * s, char ch)

//  Remove all occurrences of  ch  at the end of  s  by writing
//  '\0's over them.  Return  s  so this can be used as a function.

  {
   int  i, len;

   len = strlen (s);

   for  (i = len - 1;  i >= 0 && s [i] == ch;  i --)
     s [i] = '\0';

   return  s;
  }



