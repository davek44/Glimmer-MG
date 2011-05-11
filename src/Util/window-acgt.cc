//  A. L. Delcher
//
//  File:  window-acgt.cc
//
//  Last Modified:  Tue May 23 10:30:35 EDT 2006
//
//  Read a multifasta file from stdin and report the acgt content
//  of windows in it.  Command line arguments specify the
//  length of windows and their separation.
//  Output goes to stdout in the format
//  window-start  window-len  A's  C's  G's  T's  #other  %GC
//  Note that the last window can be shorter than the specified
//  length.



#include  "window-acgt.hh"


// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

static bool  Output_Percents = false;
  // If set true (by the -p option) then output percentages instead
  // of counts
static int  Window_Len;
  // Width of window to process; specified on command line
static int  Window_Skip;
  // Number of characters to slide window before reporting the next
  // result



int  main
    (int argc, char * argv [])

  {
   vector <char>  window;  // the actual window characters (a ring buffer)
   char  line [MAX_LINE];
   int  win_pos;    // position of the first character in the current window
   int  win_next;   // next window position to be printed
   int  win_size;   // number of characters in the current window
   int  win_sub;    // subscript of next position in  window
   int  last_pos;   // last window position printed
   int  count [5] = {0};
   int  i;
   

   Verbose = 0;

   Parse_Command_Line (argc, argv);

   window . resize (Window_Len);

   win_pos = win_next = 1;
   win_sub = win_size = last_pos = 0;

   while  (fgets (line, MAX_LINE, stdin) != NULL)
     {
      if  (First_Non_Blank (line) == '>')
          {
           if  (win_pos != last_pos)
               Finish (win_pos, win_size, win_next, count, window, win_sub);
           fputs (line, stdout);
           printf ("%8s %7s %6s %6s %6s %6s %6s %6s\n", "Position", "Length",
                "As", "Cs", "Gs", "Ts", "Other", "%GC");
           win_pos = win_next = 1;
           win_sub = win_size = last_pos = 0;
           for  (i = 0;  i < 5;  i ++)
             count [i] = 0;
          }
        else
          {
           char  * p;

           for  (p = line;  * p != '\0';  p ++)
             if  (! isspace (* p))
                 {
                  if  (win_size == Window_Len)
                      {
                       count [Subscript (window [win_sub])] --;
                         // Substract for character sliding out of the window
                       win_pos ++;
                      }
                    else
                      win_size ++;
                  count [Subscript (* p)] ++;
                  window [win_sub] = * p;
                  win_sub = (win_sub + 1) % Window_Len;

                  if  (win_size == Window_Len && win_pos == win_next)
                      {
                       Print_Line (win_pos, win_size, count);
                       last_pos = win_pos;
                       win_next += Window_Skip;
                      }
                 }
          }
     }

   if  (win_pos != last_pos)
       Finish (win_pos, win_size, win_next, count, window, win_sub);

   return  0;
  }



static void  Finish
    (int win_pos, int win_size, int win_next, int count [5],
     const vector <char> & window, int win_sub)

//  Print the final line for the information in the window
//  beginning at position  win_pos  and  containing  win_size
//  characters.  The counts to be printed are in  count .
//  The ring buffer of window characters is  window  and the current
//  position in it is  win_sub .

  {
   while  (win_pos < win_next && win_size > 0)
     {
      count [Subscript (window [win_sub])] --;
      win_pos ++;
      win_size --;
      win_sub = (win_sub + 1) % Window_Len;
     }

   if  (win_size > 0)
       Print_Line (win_pos, win_size, count);

   return;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   bool  errflg = false;
   int  ch;

   optarg = NULL;

#if  ALLOW_LONG_OPTIONS
   int  option_index = 0;
   static struct option  long_options [] = {
        {"help", 0, 0, 'h'},
        {"percen", 0, 0, 'p'},
        {0, 0, 0, 0}
      };

   while  (! errflg
        && ((ch = getopt_long (argc, argv, "hp",
                     long_options, & option_index)) != EOF))
#else
   while  (! errflg
        && ((ch = getopt (argc, argv, "hp")) != EOF))
#endif

     switch  (ch)
       {
        case  'h' :
          Usage ();
          exit (EXIT_SUCCESS);

        case  'p' :
          Output_Percents = true;
          break;

        default :
          errflg = true;
       }

   if  (errflg || optind > argc - 2)
       {
        Usage ();
        exit (EXIT_FAILURE);
       }

   Window_Len = int (strtol (argv [optind ++], NULL, 10));
   Window_Skip = int (strtol (argv [optind ++], NULL, 10));

   if  (Window_Len < 1)
       {
        sprintf (Clean_Exit_Msg_Line, "ERROR:  Bad window length = %d", Window_Len);
        SIMPLE_THROW (Clean_Exit_Msg_Line);
       }
   if  (Window_Skip < 1)
       {
        sprintf (Clean_Exit_Msg_Line, "ERROR:  Bad window skip = %d", Window_Skip);
        SIMPLE_THROW (Clean_Exit_Msg_Line);
       }

   return;
  }



static void  Print_Line
    (int win_pos, int win_size, int count [5])

//  Print the output line for window beginning at position
//   win_pos  containing  win_size  characters and counts in
//   count .

  {
   int  i;

   printf ("%8d %7d", win_pos, win_size);
   if  (Output_Percents)
       for  (i = 0;  i < 5;  i ++)
         printf (" %6.1f", Percent (count [i], win_size));
     else
       for  (i = 0;  i < 5;  i ++)
         printf (" %6d", count [i]);
   printf (" %6.1f", Percent (count [1] + count [2], win_size));
   putchar ('\n');

   return;
  }



int static  Subscript
    (char ch)

//  Return the subscript (in  count ) corresponding to  ch .

  {
   switch (tolower (ch))
     {
      case  'a' :
        return  0;
      case  'c' :
        return  1;
      case  'g' :
        return  2;
      case  't' :
        return  3;
      default :
        return  4;
     }
  }



static void  Usage
    (void)

//  Print to stderr description of options and command line for
//  this program.

  {
   fprintf (stderr,
       "USAGE:  window-acgt [options] window-len window-skip < input-file\n"
       "\n"
       "Read multi-fasta-format file from standard input.\n"
       "Print the acgt-content of windows in each sequence.\n"
       "The width of windows is <window-len>.  The number of\n"
       "positions between windows to report is <window-skip>\n"
       "So the overlap between consecutive windows is\n"
       "<window-len> - <window-skip> positions\n"
       "Output goes to standard output in the format:\n"
       "  window-start  window-len  A's C's G's T's #other %%GC\n"
       "Note the last window in the sequence can be shorter than\n"
       "<window-len> if the sequence ends prematurely\n"
       "\n"
       "Options:\n"
       " -h  or  --help\n"
       "    Print this message\n"
       " -p  or  --percent\n"
       "    Output percentages instead of counts\n"
       "\n");

   return;
  }



