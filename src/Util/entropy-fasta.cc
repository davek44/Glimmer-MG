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



#include  "entropy-fasta.hh"


// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

static bool  Brief_Output = false;
  // Determines output format
static long int  Min_Len = 0;
  // Sequences shorter than this are ignored
static double  Pos_Entropy_Profile [20] = DEFAULT_POS_ENTROPY_PROF;
  // Entropy distribution of amino-acids in genes
static double  Neg_Entropy_Profile [20] = DEFAULT_NEG_ENTROPY_PROF;
  // Entropy distribution of amino-acids in non-genes


int  main
    (int argc, char * argv [])

  {
   string  sequence, hdr;
   double  entropy [26], prob [26];
   double  ep [20];
   int  i, j;

   Verbose = 0;

   while  (Fasta_Read (stdin, sequence, hdr))
     {
      if  (sequence . length () % 3 != 0)
      {
	   cerr << hdr << " not divisible by 3\n";
	   exit(EXIT_FAILURE);
      }

      cout << ">" << hdr << "\t" << Entropy_Distance_Ratio(sequence) << "\n" << sequence << "\n";
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



