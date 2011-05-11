//  A. L. Delcher
//
//  File:  entropy-profile.cc
//
//  Last Modified:  Tue May 23 10:30:35 EDT 2006
//
//  This program reads a multifasta file of gene sequences.
//  It translates each sequence to it protein product and
//  computes and prints the natural and entropy distributions
//  of the amino acids.  It also translates the reverse-complement
//  of each sequence, and prints the same distributions for
//  those sequences.



#include  "entropy-profile.hh"


// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

static bool  Brief_Output = false;
  // Determines output format
static long int  Min_Len = 0;
  // Sequences shorter than this are ignored



int  main
    (int argc, char * argv [])

  {
   string  sequence, hdr;
   string  rev_sequence;
   double  entropy [26], prob [26];
   double  ep [20];
   double  rev_entropy [26], rev_prob [26];
   double  rev_ep [20];
   int  count [26] = {0}, rev_count [26] = {0};
   int  total = 0, rev_total = 0;;
   int  i, j;

   Verbose = 0;

   Parse_Command_Line (argc, argv);

   while  (Fasta_Read (stdin, sequence, hdr))
     {
      const char  * seq, * rev_seq;
      char  aa;
      int  i, n;

      n = sequence . length ();
      if  (n < Min_Len || n % 3 != 0)
          continue;

      rev_sequence = seq;
      Reverse_Complement (rev_sequence);

      seq = sequence . c_str ();
      rev_seq = rev_sequence . c_str ();
      for  (i = 0;  i < n;  i += 3)
        {
         aa = Codon_Translation (seq + i);
         if  (aa != '*')
             count [aa - 'A'] ++;
         aa = Codon_Translation (rev_seq + i);
         if  (aa != '*')
             rev_count [aa - 'A'] ++;
        }

     }

   for  (i = 0;  i < 26;  i ++)
     if  (IS_AMINO [i])
         {
          total += count [i];
          rev_total += rev_count [i];
         }

   Counts_To_Entropy_Profile (count, ep);
   Counts_To_Entropy_Profile (rev_count, rev_ep);

   if  (Brief_Output)
       {
        printf ("AA  %8s  %8s\n", "Positive", "Negative");
        for  (i = j = 0;  i < 26;  i ++)
          if  (IS_AMINO [i])
              {
               printf (" %c  %8.5f  %8.5f\n", 'A' + i, ep [j], rev_ep [j]);
               j ++;
              }
       }
     else
       {
        printf ("%2s %29s   %29s\n", "", "--- Forward Translation ----",
             "--- Reverse Translation ----");
        printf ("%2s %6s %6s  %6s  %6s   %6s %6s  %6s  %6s\n",
             "AA", "Count", "Percen", "Entrpy", "EFrac",
             "Count", "Percen", "Entrpy", "EFrac");
        for  (i = 0;  i < 26;  i ++)
          if  (IS_AMINO [i])
              {
               prob [i] = double (count [i]) / total;
               entropy [i] = -1.0 * prob [i] * log (prob [i]);
               rev_prob [i] = double (rev_count [i]) / rev_total;
               rev_entropy [i] = -1.0 * rev_prob [i] * log (rev_prob [i]);
              }

        for  (i = j = 0;  i < 26;  i ++)
          if  (IS_AMINO [i])
              {
               printf ("%c: %6d %5.1f%%  %6.3f  %6.3f   %6d %5.1f%%  %6.3f  %6.3f\n",
                    'A' + i, count [i], Percent (count [i], total),
                    entropy [i], ep [j],
                    rev_count [i], Percent (rev_count [i], rev_total),
                    rev_entropy [i], rev_ep [j]);
               j ++;
              }
       }

   return  0;
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
        {"brief", 0, 0, 'b'},
        {"help", 0, 0, 'h'},
        {"minlen", 1, 0, 'l'},
        {0, 0, 0, 0}
      };

   while  (! errflg
        && ((ch = getopt_long (argc, argv, "bhl:",
                     long_options, & option_index)) != EOF))
#else
   while  (! errflg
        && ((ch = getopt (argc, argv, "bhl:")) != EOF))
#endif

     switch  (ch)
       {
        case  'b' :
          Brief_Output = true;
          break;

        case  'h' :
          Usage ();
          exit (EXIT_SUCCESS);

        case  'l' :
          Min_Len = strtol (optarg, NULL, 10);
          break;

        default :
          errflg = true;
       }

   if  (errflg || optind > argc - 0)
       {
        Usage ();
        exit (EXIT_FAILURE);
       }

   return;
  }



static void  Usage
    (void)

//  Print to stderr description of options and command line for
//  this program.

  {
   fprintf (stderr,
       "USAGE:  entropy-profile [options] < input-file\n"
       "\n"
       "Read multi-fasta-format gene sequences from stdin.\n"
       "Translate each to its protein product and then print\n"
       "the natural and entropy distributions of the amino acids\n"
       "Output goes to stdout\n"
       "\n"
       "Options:\n"
       " -b\n"
       " --brief\n"
       "    Brief output:  3 columns with header line\n"
       " -h\n"
       " --help\n"
       "    Print this message\n"
       " -l <n>\n"
       " --minlen <n>\n"
       "    Don't output any sequence shorter than <n> characters\n"
       "\n");

   return;
  }



