//  A. L. Delcher
//
//  File:  extract.cc
//
//  Last Modified:  5 Aug 2005
//
//  This program reads the sequence in the file named as the
//  first command-line argument and then extracts from it
//  sequences specified by coordinates.  The resulting sequences
//  are output (in multifasta or two-string format) to stdout.



#include  "extract.hh"


// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

static char  * Coord_Info = NULL;
  // Name of file with list of coordinates (or a single coordinate
  // specification).
static bool  Fasta_Output_Format = true;
  // Determines format of output.
static bool  Is_Circular = true;
  // Determines whether the input sequence is regarded as a circular
  // genome.
static long int  Min_Len = 0;
  // Output sequences shorter than this are not printed.
static char  * Sequence_File_Name = NULL;
  // Name of file with input sequence
static bool  Skip_Start = false;
  // If set true (by -s option) then omit the first 3 letters
  // of each output sequence.
static bool  Skip_Stop = false;
  // If set true (by -t option) then omit the last 3 letters
  // of each output sequence.
static bool  Use_Direction = false;
  // If set true (by -d option) then use 4th field of coords
  // to determine direction in a circular genome.



int  main
    (int argc, char * argv [])

  {
   FILE  * fp;
   string  sequence, hdr;
   char  line [MAX_LINE], tag [MAX_LINE];
   long int  seq_len;

   Verbose = 0;

   Parse_Command_Line (argc, argv);

   fp = File_Open (Sequence_File_Name, "r");

   if  (! Fasta_Read (fp, sequence, hdr))
       {
        sprintf (Clean_Exit_Msg_Line, "ERROR:  Failed to read file %s",
             Sequence_File_Name);
        Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
       }
   fclose (fp);

   seq_len = sequence . length ();

   if  (strcmp (Coord_Info, "-") == 0)
       fp = stdin;
     else
       fp = File_Open (Coord_Info, "r");

   while  (fgets (line, MAX_LINE, fp) != NULL)
     {
      long int  i, start, end, extract_len;
      int  dir;

      if  (Use_Direction)
          {
           if  (sscanf (line, "%s %ld %ld %d", tag, & start, & end, & dir) != 4)
               {
                fprintf (stderr, "ERROR:  Skipped following coord line\n");
                fputs (line, stderr);
                continue;
               }
          }
        else
          {
           if  (sscanf (line, "%s %ld %ld", tag, & start, & end) != 3)
               {
                fprintf (stderr, "ERROR:  Skipped following coord line\n");
                fputs (line, stderr);
                continue;
               }
           if  ((start < end && (! Is_Circular || end - start <= seq_len / 2))
                    || (Is_Circular && start - end > seq_len / 2))
               dir = +1;
             else
               dir = -1;
          }

      if  (dir > 0)
          {
           extract_len = 1 + end - start;
           if  (extract_len < 0)
               extract_len += seq_len;
           if  (extract_len < Min_Len)
               continue;

           i = start - 1;
           if  (Skip_Start)
               {
                i += 3;
                extract_len -= 3;
                start += 3;
               }
           if  (Skip_Stop)
               extract_len -= 3;

           if  (extract_len >= Min_Len)
               Output_Subsequence (sequence, i, extract_len, 1, tag, start, end);
          }
        else
          {
           extract_len = 1 + start - end;
           if  (extract_len < 0)
               extract_len += seq_len;
           if  (extract_len < Min_Len)
               continue;

           i = start - 1;
           if  (Skip_Start)
               {
                i -= 3;
                extract_len -= 3;
                start -= 3;
               }
           if  (Skip_Stop)
               extract_len -= 3;

           if  (extract_len >= Min_Len)
               Output_Subsequence (sequence, i, extract_len, -1, tag, start, end);
          }
        
     }

   return  0;
  }



static void  Output_Subsequence
     (const string & seq, long int i, long int len, int incr,
      const char * tag, long int start, long int end)

//  Print to stdout the subsequence of  seq  beginning at subscript
//   i  with length  len .  If  incr  is positive, output the forward
//  strand; otherwise, output the reverse complement strand.  Use  tag ,
//   start  and  end  to label the output.

  {
   const int  fasta_width = 60;
   long int  seq_len;
   long int  ct = 0, line_ct = 0;

   if  (Fasta_Output_Format)
       printf (">%s  %ld %ld  len=%ld\n", tag, start, end, len);
     else
       printf ("%-10s ", tag);

   if  (incr > 0)
       incr = 1;
     else
       incr = -1;
   seq_len = seq . length ();
        
   while  (ct < len)
     {
      if  (i < 0)
          i += seq_len;
      else if  (i >= seq_len)
          i -= seq_len;

      if  (Fasta_Output_Format && line_ct == fasta_width)
          {
           putchar ('\n');
           line_ct = 0;
          }
      if  (incr > 0)
          putchar (seq [i]);
        else
          putchar (Complement (seq [i]));

      i += incr;
      ct ++;
      line_ct ++;
     }
   putchar ('\n');

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
        {"2_fields", 0, 0, 'd'},
        {"dir", 0, 0, 'd'},
        {"help", 0, 0, 'h'},
        {"minlen", 1, 0, 'l'},
        {"nostart", 0, 0, 's'},
        {"nostop", 0, 0, 't'},
        {"nowrap", 0, 0, 'w'},
        {0, 0, 0, 0}
      };

   while  (! errflg
        && ((ch = getopt_long (argc, argv, "2dhl:stw",
                     long_options, & option_index)) != EOF))
#else
   while  (! errflg
        && ((ch = getopt (argc, argv, "2dhl:stw")) != EOF))
#endif

     switch  (ch)
       {
        case  '2' :
          Fasta_Output_Format = false;
          break;

        case  'd' :
          Use_Direction = true;
          break;

        case  'h' :
          Usage ();
          exit (EXIT_SUCCESS);

        case  'l' :
          Min_Len = strtol (optarg, NULL, 10);
          break;

        case  's' :
          Skip_Start = true;
          break;

        case  't' :
          Skip_Stop = true;
          break;

        case  'w' :
          Is_Circular = false;
          break;

        default :
          errflg = true;
       }

   if  (errflg || optind > argc - 2)
       {
        Usage ();
        exit (EXIT_FAILURE);
       }

   Sequence_File_Name = argv [optind ++];
   Coord_Info = argv [optind ++];

   return;
  }



static void  Usage
    (void)

//  Print to stderr description of options and command line for
//  this program.

  {
   fprintf (stderr,
       "USAGE:  extract [options] <sequence-file> <coords>\n"
       "\n"
       "Read fasta-format <sequence-file> and extract from it the\n"
       "subsequences specified by <coords>.  By default, <coords>\n"
       "is the name of a file containing lines of the form\n"
       "  <tag>  <start>  <stop>  [<frame>] ...\n"
       "Coordinates are inclusive counting from 1, e.g., \"1 3\"\n"
       "represents the 1st 3 characters of the sequence.\n"
       "For each line the corresponding region of <sequence-file>\n"
       "is extracted and output (after reverse-complementing if necessary)\n"
       "Output goes to stdout in multi-fasta format, unless the -2 option\n"
       "is specified\n"
       "\n"
       "Options:\n"
       " -2\n"
       " --2_fields\n"
       "    Output each sequence as 2 fields (tag and sequence) on a single line\n"
       " -d\n"
       " --dir\n"
       "    Use the 4th column of each input line to specify the direction\n"
       "    of the sequence.  Positive is forward, negative is reverse\n"
       "    The input sequence is assumed to be circular\n"
       " -h\n"
       " --help\n"
       "    Print this message\n"
       " -l <n>\n"
       " --minlen <n>\n"
       "    Don't output any sequence shorter than <n> characters\n"
       " -s\n"
       " --nostart\n"
       "    Omit the first 3 characters of each output string\n"
       " -t\n"
       " --nostop\n"
       "    Omit the last 3 characters of each output string\n"
       " -w\n"
       " --nowrap\n"
       "    Use the actual input coordinates without any wraparound\n"
       "    that would be needed by a circular genome.  Without this\n"
       "    option, the output sequence is the shorter of the two ways\n"
       "    around the circle.  Use the -d option to specify direction\n"
       "    explicitly.\n"
       "\n");

   return;
  }



