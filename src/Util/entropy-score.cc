//  A. L. Delcher
//
//  File:  entropy-score.cc
//
//  Last Modified:  29 July 2005
//
//  This program reads the sequence in the file named as the
//  first command-line argument and then scores specified
//  regions in it by entropy distance.  Results are output
//  to  stdout .



#include  "entropy-score.hh"


// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

static char  * Coord_Info = NULL;
  // Name of file with list of coordinates (or a single coordinate
  // specification).
static bool  Is_Circular = true;
  // Determines whether the input sequence is regarded as a circular
  // genome.
static long int  Min_Len = 0;
  // Output sequences shorter than this are not printed.
static double  Neg_Entropy_Profile [20] = DEFAULT_NEG_ENTROPY_PROF;
  // Entropy distribution of amino-acids in non-genes
static double  Pos_Entropy_Profile [20] = DEFAULT_POS_ENTROPY_PROF;
  // Entropy distribution of amino-acids in genes
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
      string  buff;
      long int  i, start, end, extract_len;
      int  dir, n;

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
                start += 3;
               }
           if  (Skip_Stop)
               extract_len -= 3;

           if  (extract_len < Min_Len)
               continue;

           Forward_Strand_Transfer (buff, sequence, On_Seq (i, seq_len),
                extract_len);
          }
        else
          {
           extract_len = 1 + start - end;
           if  (extract_len < 0)
               extract_len += seq_len;

           i = start - 1;
           if  (Skip_Start)
               {
                i -= 3;
                extract_len -= 3;
                start -= 3;
               }
           if  (Skip_Stop)
               extract_len -= 3;

           if  (extract_len < Min_Len)
               continue;

           Reverse_Strand_Transfer (buff, sequence, On_Seq (i, seq_len),
                extract_len);
          }

      n = strlen (line);
      if  (line [n - 1] == '\n');
          line [n - 1] = '\0';
      printf ("%s \t%5.3f\n", line, Entropy_Distance_Ratio (buff));
     }

   return  0;
  }



static double  Entropy_Distance_Ratio
    (const string & buff)

//  Return the distance ratio for the entropy profile for the
//  gene in  buff .  The ratio is the distance to global
//   Pos_Entropy_Profile  over the distance to global  Neg_Entropy_Profile .

  {
   int  count [26] = {0};
   double  ep [20];
   double  pos_dist, neg_dist, ratio;
   char  aa;
   int  i, len;

   len = buff . length ();
   for  (i = 0; i < len;  i += 3)
     {
      aa = Codon_Translation (buff . c_str () + i);
      if  (aa != '*')
          count [aa - 'A'] ++;
     }
   Counts_To_Entropy_Profile (count, ep);

   pos_dist = neg_dist = 0.0;
   for  (i = 0;  i < 20;  i ++)
     {
      pos_dist += pow (ep [i] - Pos_Entropy_Profile [i], 2);
      neg_dist += pow (ep [i] - Neg_Entropy_Profile [i], 2);
     }

   pos_dist = sqrt (pos_dist);
   neg_dist = sqrt (neg_dist);
   if  (neg_dist == 0.0)
       {
        if  (pos_dist == 0.0)
            ratio = 1.0;
          else
            ratio = 1e3;
       }
     else
       ratio = pos_dist / neg_dist;

   return  ratio;
  }



static int  On_Seq
    (int i, int seq_len)

//  Return the subscript equivalent to  i  on a sequence of
//  length  seq_len  (with subscripts starting at 0)
//  assuming circular wraparounds.

  {
   while  (i < 0)
     i += seq_len;
   while  (seq_len <= i)
     i -= seq_len;

   return  i;
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
        {"entropy", 1, 0, 'E'},
        {"help", 0, 0, 'h'},
        {"minlen", 1, 0, 'l'},
        {"nostart", 0, 0, 's'},
        {"nostop", 0, 0, 't'},
        {"nowrap", 0, 0, 'w'},
        {0, 0, 0, 0}
      };

   while  (! errflg
        && ((ch = getopt_long (argc, argv, "2dE:hl:sw",
                     long_options, & option_index)) != EOF))
#else
   while  (! errflg
        && ((ch = getopt (argc, argv, "2dE:hl:sw")) != EOF))
#endif

     switch  (ch)
       {
        case  'd' :
          Use_Direction = true;
          break;

        case  'E' :
          Read_Entropy_Profiles (optarg, errflg);
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



static void  Read_Entropy_Profiles
    (const char * fn, bool & errflg)

//  Read positive and negative entropy profiles from the
//  file name  fn .  If not successful, set  errflg  to  true .
//  Save the entropy profiles in globals  Pos_Entropy_Profile
//  and  Neg_Entropy_Profile .

  {
   FILE  * fp;
   char  line [MAX_LINE];
   int  i;

   fp = File_Open (fn, "r");
   fgets (line, MAX_LINE, fp);  // skip header line
   for  (i = 0;  i < 20;  i ++)
     if  (fscanf (fp, "%s %lf %lf\n", line, Pos_Entropy_Profile + i,
             Neg_Entropy_Profile + i) != 3)
         {
          errflg = true;
          return;
         }

   fclose (fp);

   return;
  }



static void  Usage
    (void)

//  Print to stderr description of options and command line for
//  this program.

  {
   fprintf (stderr,
       "USAGE:  entropy-score [options] <sequence-file> <coords>\n"
       "\n"
       "Read fasta-format <sequence-file> and then score regions in\n"
       "it specified by <coords>.  By default, <coords>\n"
       "is the name of a file containing lines of the form\n"
       "  <tag>  <start>  <stop>  [<frame>] ...\n"
       "Coordinates are inclusive counting from 1, e.g., \"1 3\"\n"
       "represents the 1st 3 characters of the sequence.\n"
       "Output is the same format as <coords> put with the entropy\n"
       "distance score appended to each line\n"
       "Output goes to  stdout .\n"
       "\n"
       "Options:\n"
       " -d\n"
       " --dir\n"
       "    Use the 4th column of each input line to specify the direction\n"
       "    of the sequence.  Positive is forward, negative is reverse\n"
       "    The input sequence is assumed to be circular\n"
       " -E <filename>\n"
       " --entropy <filename>\n"
       "    Read entropy profiles from <filename>.  Format is one header\n"
       "    line, then 20 lines of 3 columns each.  Columns are amino acid,\n"
       "    positive entropy, negative entropy.  Rows must be in order\n"
       "    by amino acid code letter\n"
       " -h\n"
       " --help\n"
       "    Print this message\n"
       " -l <n>\n"
       " --minlen <n>\n"
       "    Don't output any sequence shorter than <n> characters\n"
       " -s\n"
       " --nostart\n"
       "    Omit the first 3 characters of each specified string\n"
       " -t\n"
       " --nostop\n"
       "    Omit the last 3 characters of each specified string\n"
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



