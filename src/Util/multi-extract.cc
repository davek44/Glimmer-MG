//  A. L. Delcher
//
//  File:  multi-extract.cc
//
//  Last Modified:  Tue Aug  9 09:30:18 EDT 2005
//
//  This program reads the sequences in the file named as the
//  first command-line argument and then extracts from it
//  subsequences specified by the second command-line file.  The
//  resulting sequences are output (in multifasta or two-string format)
//  to stdout.



#include  "multi-extract.hh"


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
   vector <Coord_Info_t>  coord_list;
   Coord_Info_t  coord;
   string  sequence, hdr;
   char  line [MAX_LINE], id [MAX_LINE], tag [MAX_LINE];
   long int  seq_len;

   Verbose = 0;

   Parse_Command_Line (argc, argv);

   // Read and store the coordinate information
   if  (strcmp (Coord_Info, "-") == 0)
       fp = stdin;
     else
       fp = File_Open (Coord_Info, "r");

   while  (fgets (line, MAX_LINE, fp) != NULL)
     {
      long int  start, end;
      int  dir;

      if  (Use_Direction)
          {
           if  (sscanf (line, "%s %s %ld %ld %d", id, tag, & start, & end, & dir) != 5)
               {
                fprintf (stderr, "ERROR:  Skipped following coord line\n");
                fputs (line, stderr);
                continue;
               }
          }
        else
          {
           if  (sscanf (line, "%s %s %ld %ld", id, tag, & start, & end) != 4)
               {
                fprintf (stderr, "ERROR:  Skipped following coord line\n");
                fputs (line, stderr);
                continue;
               }
          }
      coord . id = strdup (id);
      coord . tag = strdup (tag);
      coord . start = start;
      coord . end = end;
      coord . dir = (Use_Direction ? dir : 0);
      coord_list . push_back (coord);
     }

   fclose (fp);

   // Sort the coordinate information by tag

   sort (coord_list . begin (), coord_list . end (), By_Tag);


   fp = File_Open (Sequence_File_Name, "r");

   while  (Fasta_Read (fp, sequence, hdr))
     {
      char  * p;
      long int  extract_len, loc;
      int  i, dir, sub, num;

      strcpy (tag, hdr . c_str ());
      p = strtok (tag, " \t\n");

      Find_Matches (p, coord_list, sub, num);

      if  (num == 0)
          continue;

      seq_len = sequence . length ();

      for  (i = sub;  i - sub < num;  i ++)
        {
         Coord_Info_t  * cp = & coord_list [i];

         if  (Use_Direction)
             dir = cp -> dir;
         else if  ((cp -> start < cp -> end
                      && (! Is_Circular || cp -> end - cp -> start <= seq_len / 2))
                       || (Is_Circular && cp -> start - cp -> end > seq_len / 2))
             dir = +1;
           else
             dir = -1;

         if  (dir > 0)
             {
              extract_len = 1 + cp -> end - cp -> start;
              if  (extract_len < 0)
                  extract_len += seq_len;

              loc = cp -> start - 1;
              if  (Skip_Start)
                  {
                   loc += 3;
                   extract_len -= 3;
                  }
              if  (Skip_Stop)
                  extract_len -= 3;
             }
           else
             {
              extract_len = 1 + cp -> start - cp -> end;
              if  (extract_len < 0)
                  extract_len += seq_len;

              loc = cp -> start - 1;
              if  (Skip_Start)
                  {
                   loc -= 3;
                   extract_len -= 3;
                  }
              if  (Skip_Stop)
                  extract_len -= 3;
             }

         if  (extract_len >= Min_Len)
             Output_Subsequence (sequence, loc, extract_len, dir, cp -> id,
                  cp -> tag, cp -> start, cp -> end);
        }
     }

   fclose (fp);

   return  0;
  }



static void  Find_Matches
    (char * p, const vector <Coord_Info_t> & list, int & sub, int & num)

//  Find in  list  the entries whose  tag  matches  p .  Set  sub  to
//  the subscript of the first match and  num  to the number of matches.
//  If no matches are found, set  num  to zero.   list  must be
//  sorted into ascending order by  tag  so that the matches will
//  be consecutive in  list .

  {
   int  i, n;

   n = list . size ();
   for  (i = 0;  i < n && strcmp (p, list [i] . tag) > 0;  i ++)
     ;

   sub = i;
   num = 0;

   if  (i == n || strcmp (p, list [i] . tag) != 0)
       return;


   for  ( ;  i < n && strcmp (p, list [i] . tag) == 0;  i ++)
     num ++;

   return;
  }



static void  Output_Subsequence
     (const string & seq, long int i, long int len, int incr,
      const char * id, const char * tag, long int start, long int end)

//  Print to stdout the subsequence of  seq  beginning at subscript
//   i  with length  len .  If  incr  is positive, output the forward
//  strand; otherwise, output the reverse complement strand.  Use  id ,
//   tag ,  start  and  end  to label the output.

  {
   const int  fasta_width = 60;
   long int  seq_len;
   long int  ct = 0, line_ct = 0;

   if  (Fasta_Output_Format)
       printf (">%s  %s  %ld %ld  len=%ld\n", id, tag, start, end, len);
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
       "USAGE:  multi-extract [options] <sequence-file> <coords>\n"
       "\n"
       "Read multi-fasta-format <sequence-file> and extract from it the\n"
       "subsequences specified by <coords>.  By default, <coords>\n"
       "is the name of a file containing lines of the form\n"
       "  <id>  <tag>  <start>  <stop>  [<frame>] ...\n"
       "<id> is the identifier for the subsequence\n"
       "<tag> is the tag of the sequence in <sequence-file> from which\n"
       "to extract the entry\n"
       "Coordinates are inclusive counting from 1, e.g., \"1 3\"\n"
       "represents the 1st 3 characters of the sequence.\n"
       "Specify \"-\" for <coords> to read it from standard input\n"
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



