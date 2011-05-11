//  A. L. Delcher
//
//  File:  uncovered.cc
//
//  Last Modified:  10 June 2005
//
//  This program reads the sequence in the file named as the
//  first command-line argument and then extracts from it
//  sequences that are not covered by the list of regions
//  specified in the file named as the second command-line argument.
//  Output is a multifasta file sent to stdout.



#include  "uncovered.hh"


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
   vector <Region_t>  region_list;
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
      Region_t  region;
      long int  i, j, start, end, extract_len;
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

           i = start - 1;
           if  (Skip_Start)
               {
                i += 3;
                extract_len -= 3;
               }
           if  (Skip_Stop)
               extract_len -= 3;

           j = i + extract_len;
           if  (j <= seq_len)
               {
                region . lo = i;
                region . hi = j;
                region_list . push_back (region);
               }
             else  // wraparound end
               {
                region . lo = i;
                region . hi = seq_len;
                region_list . push_back (region);
                region . lo = 0;
                region . hi = j - seq_len;
                region_list . push_back (region);
               }
          }
        else
          {
           extract_len = 1 + start - end;
           if  (extract_len < 0)
               extract_len += seq_len;

           i = start;
           if  (Skip_Start)
               {
                i -= 3;
                extract_len -= 3;
               }
           if  (Skip_Stop)
               extract_len -= 3;

           j = i - extract_len;
           if  (j >= 0)
               {
                region . lo = j;
                region . hi = i;
                region_list . push_back (region);
               }
             else  // wraparound beginning
               {
                region . lo = 0;
                region . hi = i;
                region_list . push_back (region);
                region . lo = seq_len + j;
                region . hi = seq_len;
                region_list . push_back (region);
               }
          }
     }

   Coalesce_Regions (region_list);

   Output_Uncovered (sequence, seq_len, region_list);

   return  0;
  }



static void  Coalesce_Regions
    (vector <Region_t> & list)

//  Merge overlapping regions in  list  changing list to
//  a sorted list of the merged regions.

  {
   int  i, j, n;

   n = list . size ();
   if  (n == 0)
       return;

   sort (list . begin (), list . end (), Region_Cmp);

   j = 0;
   for  (i = 1;  i < n;  i ++)
     if  (list [i] . lo <= list [j] . hi)
         {
          if  (list [j] . hi < list [i] . hi)
              list [j] . hi = list [i] . hi;
         }
       else
         list [++ j] = list [i];

   list . resize (j + 1);

   return;
  }



static void  Output_Subsequence
     (const string & seq, long int i, long int len,
      const char * tag, long int start, long int end)

//  Print to stdout the subsequence of  seq  beginning at subscript
//   i  with length  len .  Use  tag ,  start  and  end  to label the output.

  {
   const int  fasta_width = 60;
   long int  ct = 0, line_ct = 0;

   if  (Fasta_Output_Format)
       printf (">%s  %ld %ld  len=%ld\n", tag, start, end, len);
     else
       printf ("%-10s ", tag);

   while  (ct < len)
     {
      if  (Fasta_Output_Format && line_ct == fasta_width)
          {
           putchar ('\n');
           line_ct = 0;
          }

      putchar (seq [i]);
      i ++;
      ct ++;
      line_ct ++;
     }
   putchar ('\n');

   return;
  }



static void  Output_Uncovered
    (const string & seq, long int seq_len, const vector <Region_t> & list)

//  Output the portions of  seq  that are not in  list  if they
//  are long enough.

  {
   long int  a, b, len;
   char  tag [100];
   int  i, n, ct = 0;

   n = list . size ();
   a = 0;

   for  (i = 0;  i < n;  i ++)
     {
      b = list [i] . lo;
      len = b - a;
      if  (len > 0 && len >= Min_Len)
          {
           sprintf (tag, "seq%05d", ++ ct);
           Output_Subsequence (seq, a, len, tag, a + 1, b);
          }
      a = list [i] . hi;
     }

   len = seq_len - a;
   if  (len > 0 && len >= Min_Len)
       {
        sprintf (tag, "seq%05d", ++ ct);
        Output_Subsequence (seq, a, len, tag, a + 1, seq_len);
       }

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
        {"dir", 0, 0, 'd'},
        {"minlen", 1, 0, 'l'},
        {"nostart", 0, 0, 's'},
        {"nostop", 0, 0, 't'},
        {"nowrap", 0, 0, 'w'},
        {0, 0, 0, 0}
      };

   while  (! errflg
        && ((ch = getopt_long (argc, argv, "2dhl:sw",
                     long_options, & option_index)) != EOF))
#else
   while  (! errflg
        && ((ch = getopt (argc, argv, "2dhl:sw")) != EOF))
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
       "USAGE:  uncovered [options] <sequence-file> <coords>\n"
       "\n"
       "Read fasta-format <sequence-file> and extract from it the\n"
       "subsequences not covered by regions specified in <coords>.\n"
       "By default, <coords>\n"
       "is the name of a file containing lines of the form\n"
       "  <tag>  <start>  <stop>  [<frame>] ...\n"
       "Coordinates are inclusive counting from 1, e.g., \"1 3\"\n"
       "represents the 1st 3 characters of the sequence.\n"
       "Output goes to stdout in multi-fasta format, unless the -2 option\n"
       "is specified\n"
       "\n"
       "Options:\n"
       " -2"
       "    Output each sequence as 2 fields (tag and sequence) on a single line\n"
       " -d\n"
       " --dir\n"
       "    Use the 4th column of each input line to specify the direction\n"
       "    of the sequence.  Positive is forward, negative is reverse\n"
       "    The input sequence is assumed to be circular\n"
       " -h\n"
       "    Print this message\n"
       " -l <n>\n"
       " --minlen <n>\n"
       "    Don't output any sequence shorter than <n> characters\n"
       " -s\n"
       " --nostart\n"
       "    Omit the first 3 characters of each <coords> region\n"
       " -t\n"
       " --nostop\n"
       "    Omit the last 3 characters of each <coords> region\n"
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



