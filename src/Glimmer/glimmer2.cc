//
//  Programmer:  A. Delcher
//
//        File:  ~delcher/Glimmer2/glimmer2.cc
//
//     Version:  2.01  31 Jul 98
//                 Change probability model
//                 Simplify wraparounds
//                 Move start codons to eliminate overlaps
//                 Discount independent model scores when
//                    there are no overlaps
//                 Uses Harmon's model
//
//     Version:  2.03  9 Dec 2002
//               Include raw scores in output
//               Add strict option to use independent intergenic
//                 model that discounts stop codons
//               Add option to score each entry from a list of coordinates
//                 separately, without overlapping/voting rules
//
//     Version:  2.10  5 Feb 2003
//               Strict option to use independent intergenic
//                 model that discounts stop codons is only behaviour
//
//    Copyright (c) 1999 by Arthur Delcher, Steven Salzberg, Simon
//    Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  This program finds open reading frames in the file named
//  on the command line and scores them using the probability
//  model in the file indicated by the second command-line
//  parameter.
//


#include  "delcher.h"
#include  "gene.h"


const int  DEFAULT_MIN_GENE_LEN = 90;
const double  DEFAULT_MIN_OLAP_PERCENT = 0.10;
const int  DEFAULT_MIN_WEAK_LEN = INT_MAX;
const int  DEFAULT_SHOW_FRAME_AND_LENGTH = TRUE;
const int  DEFAULT_THRESHOLD_SCORE = 90;
const int  DEFAULT_MIN_OLAP = 30;
const int  DEFAULT_CHOOSE_FIRST_START_CODON = TRUE;
const int  DEFAULT_USE_INDEPENDENT = TRUE;
const int  DEFAULT_USE_STRICT_INDEPENDENT = TRUE;
const int  DEFAULT_VOTE_THRESHOLD = 150;
const int  DEFAULT_IGNORE_OPTION = FALSE;
const int  MAX_ITERATIONS = 10;
const int  MAX_RIBOSOME_PATTERN_LEN = 15;
const int  MAX_START_SHIFT = 150;
const int  MIN_LONG_GENE_LEN = 90;                   // Was 240
const double  MIN_PERCENT_LEN_DIFF = 0.10;
const int  MIN_SCORABLE_LEN = 60;
const int  OLAP_THRESHOLD_SCORE = 50;
const double  ATG_THRESHOLD = 7.0;
const double  GTG_THRESHOLD = 8.0;
const double  TTG_THRESHOLD = 9.0;
const double  BIG_NEGATIVE = -1000.0;
const int  MAX_FREE_LEN = 3;
const int  NUM_FRAME_MODELS = 3;
const int  ORF_SIZE_INCR = 1000;
const double  SMALL_OLAP_PERCENT = 0.10;
const int  UPSTREAM_LEN = 15;
const int  UPSTREAM_OFFSET = 2;

                         //  Tail of 16s ribosomal RNA in *reverse* order
// #define  DEFAULT_RIBOSOME_PATTERN  "tctttcctccac"   // For SangerTB
// #define  DEFAULT_RIBOSOME_PATTERN  "tcctccactagg"   // For H.pylori
#define  DEFAULT_RIBOSOME_PATTERN  "tcctcca"   // For H.pylori

const int  MODEL_LEN = 12;
const int  ALPHABET_SIZE = 4;
const int  MAX_NAME_LEN = 256;
const double  MAX_LOG_DIFF = -46.0;    // approx 1e-20 difference in probs


#include  "icm.h"

const unsigned int  OK = 0x0;
const unsigned int  REJECTED = 0x1;
const unsigned int  SHADOWED_BY = 0x2;
const unsigned int  SHADOWS_ANOTHER = 0x4;
const unsigned int  RBS_START_SHIFT = 0x8;
const unsigned int  MIGHT_CHANGE = 0x10;
const unsigned int  JUNK_ORF = 0x20;
const unsigned int  GIVEN_GENE = 0x40;
const unsigned int  OFF_LIMITS = 0x80;

const char  CONTAINS_CHAR = '>';
const char  ELIM_OLAP_CHAR = 'E';
const char  NEAR_REJECT_CHAR = 'N';
const char  NOPROB_CHAR = ' ';
const char  REJECT_CHAR = 'R';
const char  SCORES_WORSE_CHAR = 'B';
const char  SHADOWED_CHAR = '<';
const char  SHORTER_CHAR = 'S';

enum  Gene_Category_Type
         {NONE, REGULAR, VOTED, WEAK};
enum  Olap_Fix_Type
         {NEITHER_CAN_MOVE, ONLY_A_CAN_MOVE, ONLY_B_CAN_MOVE, BOTH_CAN_MOVE};

struct  Overlap_Node
  {
   char  Problem_Code;
   long int  From, Olap, Delay, Lo, Hi;
   int  Other_Frame, Score;
   Overlap_Node  * Next;
  };

struct  Gene_Ref
  {
   long int  Lo, Hi, Len, Max_Hi, Min_Lo, Top;
   long int  Delay_Len, Delay_Cause;
   int  Frame;
   unsigned int  Status;
   Gene_Category_Type  Category;
   int  Score;
   double  Raw_Score;
   Overlap_Node  * Overlap_List;
   void  Clear_Status  (unsigned int S)      // Make Status not include S
     {
      Status &= ~ S;
     }
   int  Has_Status  (unsigned int S)         // Return whether Status matches S
     {
      return  ((Status & S) == S);
     }
   void  Set_Status  (unsigned int S)        // Set Status to include S
     {
      Status |= S;
     }
  };

struct  ED_Struct
  {
   float  Free_i, Free_j, Both_Free, Match;
   int  Free_i_Len, Free_j_Len, Both_i_Len, Both_j_Len;
  }  ED_Score [1 + MAX_RIBOSOME_PATTERN_LEN] [1 + UPSTREAM_LEN];

struct Ignore_Node
  {
    int Low_Address;
    int High_Address;
    int Frame;
  };
typedef  Ignore_Node  * Ignore_Ptr;


void  Add_Overlap
    (Gene_Ref &, long int, long int, long int, char, int);
void  Append_Gene_Ref
    (Gene_Ref * &, long int, long int, long int, int, int, Gene_Category_Type,
     double raw_score);
double  Bulge_Cost
    (int);
int  Can_Delay_Start
    (Gene_Ref &, Overlap_Node *, char, long int);
int  Choose_Score
    (int [7], int);
long int  Choose_Start
    (long int, long int);
void  Determine_Changes
    (Gene_Ref * [], Gene_Ref [], long int);
double  Doublet_Score
    (char, char, char, char);
double  Edit_Distance
    (const char *, const char *);
void  Evaluate_Overlap
    (Overlap_Node *, long int, Gene_Ref []);
long int  Extend_Data
    (char * & Data, long int Len, long int Max_Extend);
void  Find_Overlaps
    (Gene_Ref *, long int);
int  Gene_Ref_Cmp
    (const void *, const void *);
Olap_Fix_Type  Get_Olap_Fix
    (long int, long int, int, long int, long int, int);
double  Loop_Cost
    (int, int);
void  Make_Codon_Log_Prob
    (double codon_log_prob [64], double base_prob [4]);
int  Make_Final_Changes
    (Gene_Ref * [], Gene_Ref [], long int);
int  Make_Sure_Changes
    (Gene_Ref [], long int);
int  Match
    (char, char);
void  Permute
    (int [], int);
static void  Print_Category
    (Gene_Category_Type Category);
static void  Print_Separate_Score_Headings
    (void);
void  Process_Ignore ( Ignore_Ptr & );
void  Process_Options
    (int, char * []);
void  Process_Orflist
    (void);
void  Read_Probability_Model
    (char *);
void  RNAbin
    (char * Data, long int genomeLen, long int coords[][2],
     long int newStartsArr [], long int N);
  // Olga's new function
void  Score_Multifasta_Orfs
    (FILE * fp);
void  Score_Olap_Region
    (long int, long int, int, int, int &, int &);
void  Score_String
    (long int, long int, double [ALPHABET_SIZE], int [7],
     int, int &, double &);
static void  Set_Ignore_Indep_Len
    (void);
  // Set  Ignore_Indep_Len  based on  Ch_Ct  if it's not
  // already set by command-line option
void  Set_Indep_Probs_From_Data
    (double Ch_Ct [], FILE * fp);
void  Show_Gene_Info
    (Gene_Ref [], long int);
void  Simple_Score
    (char [], long int, int, double [ALPHABET_SIZE], int [7],
     int, int &, double &);
void  Slide_Both_Starts
    (Gene_Ref &, Overlap_Node *, char, Gene_Ref &, Overlap_Node *,
     long int);
void  Slide_One_Start
    (Gene_Ref &, Overlap_Node *, long int &);
void  Transfer
    (char *, long int, int);
static void  Usage
    (char * command);



static int  Allow_Partial_Orfs = FALSE;
  // If set true by -X option, then score orfs that
  // extend to the end of the sequence
double  Ch_Ct [ALPHABET_SIZE] = {0.0};
int  Choose_First_Start_Codon = DEFAULT_CHOOSE_FIRST_START_CODON;
static double  Codon_Log_Prob [64] = {0.0};
  // Log of probability of non-stop codons using independent probabilities
  // of each base
char  * Data;
long int  Data_Len;
int  Ignore_Option = DEFAULT_IGNORE_OPTION;
char  * Ignore_File_Name;
long int  Extended_Len;
static bool  GC_From_Parameter = false;
  // If true, then GC content for independent model comes from
  // the value specified by the -C option
int  Genome_Is_Circular = TRUE;
long int  Ignore_Indep_Len = LONG_MAX;
static long int  Input_Size = INIT_SIZE;
  // Size of input string buffer
long int  Min_Gene_Len = DEFAULT_MIN_GENE_LEN;
long int  Min_Olap = DEFAULT_MIN_OLAP;
double  Min_Olap_Percent = DEFAULT_MIN_OLAP_PERCENT;
long int  Min_Weak_Len = DEFAULT_MIN_WEAK_LEN;
static char  Name [MAX_LINE];
  // First ID string on fasta header line
char  * Orf_Buffer;
long int  Orf_Buffer_Len;
int  Orflist_Option = FALSE;
char  * Orflist_File_Name;
long int  (* Original_Coord) [2];
long int  * Revised_Start;
char  Ribosome_Pattern [1 + MAX_RIBOSOME_PATTERN_LEN] = DEFAULT_RIBOSOME_PATTERN;
static bool  Separate_Multifasta_Orfs = false;
  // If set true by -M option then input is multifasta file
  // of orfs to be scored separately (like Orflist_Option)
int  Show_Frame_And_Length = DEFAULT_SHOW_FRAME_AND_LENGTH;
int  Threshold_Score = DEFAULT_THRESHOLD_SCORE;
int  Use_Independent = DEFAULT_USE_INDEPENDENT;
int  Use_Strict_Independent = DEFAULT_USE_STRICT_INDEPENDENT;
int  Vote_Threshold = DEFAULT_VOTE_THRESHOLD;

int  Scoring_Overlap = FALSE;                       // Temporary
long int  Big_Diff_Ct = 0, Small_Diff_Ct = 0;       // Temporary
double  Diff_Sum = 0.0;                             // Temporary
int  Global_Show_Details = FALSE;                   // Temporary
int  Global_Check = FALSE;                          // Temporary



int  main
    (int argc, char * argv [])

  {
   FILE  * fp;
   Gene_Ref  * Gene, * * Ptr;
   Ignore_Ptr  Ignore;
   long int  Ignore_Seg, Initial_Base, Max_Extend;
   Overlap_Node  * P;
   int  Frame, In_Frame_Score, Is_Tentative_Gene, Score [7], Weak_Score;
   double  Raw_Score;
   long int  Votes [6] = {0};
   unsigned  Codon;
   Gene_Category_Type  Category;
   long int  For_Prev [3] = {LONG_MAX, LONG_MAX, LONG_MAX};
   long int  Rev_Prev [3] = {LONG_MAX, LONG_MAX, LONG_MAX};
   long int  For_Start [3] = {0};
   long int  Rev_Start [3] = {0};
   long int  ID_Num = 0;
   long int  Changes_Made;
   long int  Lo, Hi;
   long int  i, j, Len, Gene_Len, Start;

   if  (argc < 3)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   Process_Options (argc, argv);

   Process_Ignore (Ignore);

   Data = (char *) Safe_malloc (Input_Size);

   fp = File_Open (argv [1], "r");

   Read_Probability_Model (argv [2]);

   if  (Separate_Multifasta_Orfs)
       {
        if  (! GC_From_Parameter)
            Set_Indep_Probs_From_Data (Ch_Ct, fp);

        Set_Ignore_Indep_Len ();

        Make_Codon_Log_Prob (Codon_Log_Prob, Ch_Ct);

        for  (i = 0;  i < ALPHABET_SIZE;  i ++)
          Ch_Ct [i] = log (Ch_Ct [i]);

        Print_Separate_Score_Headings ();

        Score_Multifasta_Orfs (fp);

        return  0;
       }

   Read_String (fp, Data, Input_Size, Name, FALSE);
   fclose (fp);

   Data_Len = strlen (Data + 1);
   for  (i = 1;  i <= Data_Len;  i ++)
     {
      Data [i] = Filter (tolower (Data [i]));
      if  (! GC_From_Parameter)
          {
           switch  (Data [i])
             {
              case  'a' :
              case  't' :
                Ch_Ct [0] += 1.0;
                break;
              case  'c' :
              case  'g' :
                Ch_Ct [1] += 1.0;
                break;
             }
          }
     }

   if  (! GC_From_Parameter)
       {
        Ch_Ct [2] = Ch_Ct [1];
        Ch_Ct [3] = Ch_Ct [0];
        for  (i = 0;  i < ALPHABET_SIZE;  i ++)
          Ch_Ct [i] = Ch_Ct [i] / (2.0 * Data_Len);
       }

   Set_Ignore_Indep_Len ();

   Make_Codon_Log_Prob (Codon_Log_Prob, Ch_Ct);

   for  (i = 0;  i < ALPHABET_SIZE;  i ++)
     Ch_Ct [i] = log (Ch_Ct [i]);
   
   Orf_Buffer_Len = ORF_SIZE_INCR;
   Orf_Buffer = (char *) Safe_malloc (Orf_Buffer_Len);
   Orf_Buffer [0] = ' ';

   Gene = (Gene_Ref *) Safe_malloc (sizeof (Gene_Ref));

   if  (Orflist_Option)
       Print_Separate_Score_Headings ();
     else
       {
        printf ("Minimum gene length = %ld\n", Min_Gene_Len);
        printf ("Minimum overlap length = %ld\n", Min_Olap);
        printf ("Minimum overlap percent = %.1f%%\n", 100.0 * Min_Olap_Percent);
        printf ("Threshold score = %d\n", Threshold_Score);
        printf ("Use independent scores = %s\n", Use_Independent ? "True" : "False");
        if  (Use_Independent)
            printf ("Ignore independent score on orfs longer than %ld\n",
                       Ignore_Indep_Len - 1);
        printf ("Use strict independent model = %s\n",
                Use_Strict_Independent ? "True" : "False");
        printf ("Use first start codon = %s\n",
                       Choose_First_Start_Codon ? "True" : "False");
        if  (! Choose_First_Start_Codon)
            printf ("   Ribosome pattern = %s\n", Ribosome_Pattern);

        if( Ignore_Option )
          {
            printf("The ignore file was: \"%s\" and ", Ignore_File_Name);  
            printf("the following regions were skipped:\n");
            for(i = 0; Ignore[i].Low_Address > -1; i++)
              printf("     start: %10d   |   end: %10d \n",
                     Ignore[i].Low_Address, Ignore[i].High_Address );
          }

        putchar ('\n');
        printf ("              Orf     Gene                 Lengths"
                "     Gene    -- Frame Scores -");
        if  (Use_Independent)
            printf ("  Indep");
        putchar ('\n');
        printf ("  ID#  Fr    Start    Start      End      Orf  Gene"
                "    Score   F1 F2 F3 R1 R2 R3");
        if  (Use_Independent)
            printf ("  Score");
        putchar ('\n');
       }

   if  (Genome_Is_Circular)
       {
	 if ( Ignore_Option )
	   Max_Extend = Ignore[0].Low_Address;
	 else
	   Max_Extend = Data_Len;

        Extended_Len = Extend_Data  (Data, Data_Len, Max_Extend);
        Codon = Ch_Mask (Data [Data_Len - 1]) << 4 | Ch_Mask (Data [Data_Len]);
       }
     else
       {
        Extended_Len = Data_Len;
        Codon = Ch_Mask ('g') << 4 | Ch_Mask ('g');
        for  (i = 0;  i < 3;  i ++)
          For_Prev [i] = i;
       }

   Frame = 0;
   Initial_Base = 1;  
   Ignore_Seg = 0;

   if ( Ignore_Option )
     {
       if (Ignore[0].Low_Address == 0 && Ignore[0].High_Address != 0)
	 Initial_Base = Ignore[0].High_Address;
     }

   if  (Orflist_Option)
       {
        Process_Orflist ();

        return  0;
       }

   if  (Allow_Partial_Orfs)
       {
        if  (Genome_Is_Circular)
            {
             fprintf (stderr, "ERROR:  Must use -l option with -X option\n");
             exit (EXIT_FAILURE);
            }
        
        for  (i = 0;  i < 3;  i ++)
          {
           For_Start [i] = i + 1;
           Rev_Prev [i] = i;
          }
        Extended_Len = Data_Len + 3;
        Data = (char *) Safe_realloc (Data, Extended_Len + 2);
        strcat (Data, "tag");
       }


   for  (i = Initial_Base;  i <= Extended_Len;  i ++)
     {
      Frame = (Frame + 1) % 3;

      if  (i == Ignore[Ignore_Seg].Low_Address)
        {
          i = Ignore[Ignore_Seg].High_Address;
          for ( j = 0; j < 3; j++ )
            {
              For_Prev  [j] = LONG_MAX;
	      Rev_Prev  [j] = LONG_MAX;
              For_Start [j] = 0;
              Rev_Start [j] = 0;
            }
	  for ( j = 0; j < 6; j++ )
	    Votes[j] = 0;
          Ignore_Seg++;
	  Codon = Ch_Mask (Data [i - 1]) << 4 | Ch_Mask (Data [i]);
          i++;
        }
       
      Codon = (Codon & SHIFT_MASK) << 4;
      Codon |= Ch_Mask (Data [i]);

      if  (Is_Forward_Stop (Codon)
             || (Allow_Partial_Orfs && i > Data_Len))
          {
           Len = i - For_Prev [Frame] - 3;
           if  (Len >= Min_Gene_Len  && For_Prev [Frame] <= Data_Len)
                if  (For_Start [Frame] != 0)
                    {
                     Gene_Len = 1 + i - 3 - For_Start [Frame];
                     Start = Choose_Start (For_Start [Frame], Gene_Len);
                     Gene_Len = 1 + i - 3 - Start;
                     if  (Gene_Len >= Min_Gene_Len)
                         {
                          Score_String (Start, i - 3, Ch_Ct, Score,
                                        Use_Independent && Len < Ignore_Indep_Len,
                                        Weak_Score, Raw_Score);
                          In_Frame_Score = Score [0];
                          Is_Tentative_Gene = (In_Frame_Score >= Threshold_Score
                                 || Votes [Frame] + In_Frame_Score
                                            >= Vote_Threshold
                                 || (Weak_Score >= Threshold_Score
                                       && Gene_Len >= Min_Weak_Len));
                          if  (In_Frame_Score >= Threshold_Score)
                              Category = REGULAR;
                          else if  (Votes [Frame] + In_Frame_Score
                                            >= Vote_Threshold)
                              Category = VOTED;
                          else if  (Weak_Score >= Threshold_Score
                                       && Gene_Len >= Min_Weak_Len)
                              Category = WEAK;
                            else
                              Category = NONE;
                          if  (Is_Tentative_Gene)
                              printf ("%5ld ", ++ ID_Num);
                            else
                              printf ("%5s ", "");
                          printf (" F%1d %8ld %8ld %8ld %8ld %5ld   ",
                                 Frame + 1, For_Prev [Frame] + 1,
                                 Start, i - 3,
                                 Len, Gene_Len);
                          printf (" %4d   ", In_Frame_Score);
                          Permute (Score, Frame + 1);
                          for  (j = 0;  j < 6;  j ++)
                            if  (Score [j] < 0)
                                printf ("  _");
                              else
                                {
                                 printf (" %2d", Score [j]);
                                 Votes [j] += Score [j];
                                }
                          if  (Use_Independent)
                              printf ("   %2d", Score [6]);
                          printf ("  %4ld", Votes [Frame]);
                          printf ("  %6.3f", Raw_Score);
                          Print_Category (Category);
                          putchar ('\n');
                          if  (Is_Tentative_Gene)
                              Append_Gene_Ref (Gene, ID_Num, Start,
                                   i - 3, Frame + 1, In_Frame_Score, Category,
                                   Raw_Score);
                         }
                    }
           For_Prev [Frame] = i;
           For_Start [Frame] = 0;
           Votes [Frame] = 0;
          }
      if  (Is_Forward_Start (Codon) && For_Start [Frame] == 0)
          For_Start [Frame] = i - 2;

      if  (Is_Reverse_Stop (Codon)
             || (Allow_Partial_Orfs && i > Data_Len))
          {
           if  (Allow_Partial_Orfs && i > Data_Len)
               Rev_Start [Frame] = i - 3;
           Len = i - Rev_Prev [Frame] - 3;
           if  (Len >= Min_Gene_Len && Rev_Prev [Frame] <= Data_Len)
                if  (Rev_Start [Frame] != 0)
                    {
                     Gene_Len = Rev_Start [Frame] - Rev_Prev [Frame];
                     Start = Choose_Start (Rev_Start [Frame], - Gene_Len);
                     Gene_Len = Start - Rev_Prev [Frame];
                     if  (Gene_Len >= Min_Gene_Len)
                         {
                          Score_String (Start,
                                    Rev_Prev [Frame] + 1, Ch_Ct, Score,
                                    Use_Independent && Len < Ignore_Indep_Len,
                                    Weak_Score, Raw_Score);
                          In_Frame_Score = Score [0];
                          Is_Tentative_Gene = (In_Frame_Score >= Threshold_Score
                                 || Votes [3 + (3 - Frame) % 3] + In_Frame_Score
                                        >= Vote_Threshold
                                 || (Weak_Score >= Threshold_Score
                                       && Gene_Len >= Min_Weak_Len));
                          if  (In_Frame_Score >= Threshold_Score)
                              Category = REGULAR;
                          else if  (Votes [3 + (3 - Frame) % 3] + In_Frame_Score
                                        >= Vote_Threshold)
                              Category = VOTED;
                          else if  (Weak_Score >= Threshold_Score
                                       && Gene_Len >= Min_Weak_Len)
                              Category = WEAK;
                            else
                              Category = NONE;
                          if  (Is_Tentative_Gene)
                              printf ("%5ld ", ++ ID_Num);
                            else
                              printf ("%5s ", "");
                          printf (" R%1d %8ld %8ld %8ld %8ld %5ld   ",
                                 1 + ((2 - Frame) + 1) % 3,
                                 i - 3, Start,
                                 Rev_Prev [Frame] + 1, Len,
                                 Gene_Len);
                          printf (" %4d   ", In_Frame_Score);
                          Permute (Score, - Frame - 1);
                          for  (j = 0;  j < 6;  j ++)
                            if  (Score [j] < 0)
                                printf ("  _");
                              else
                                {
                                 printf (" %2d", Score [j]);
                                 if  ((j < 3 && For_Prev [j] < Rev_Start [Frame])
                                       || (j >= 3
                                          && Rev_Prev [(6 - j) % 3] < Rev_Start [Frame]))
                                     Votes [j] += Score [j];
                                }
                          if  (Use_Independent)
                              printf ("   %2d", Score [6]);
                          printf ("  %4ld", Votes [3 + (3 - Frame) % 3]);
                          printf ("  %6.3f", Raw_Score);
                          Print_Category (Category);
                          putchar ('\n');
                          if  (Is_Tentative_Gene)
                              Append_Gene_Ref (Gene, ID_Num, Rev_Prev [Frame] + 1,
                                   Start , - Frame - 1, In_Frame_Score, Category,
                                   Raw_Score);
                         }
                    }
           Rev_Prev [Frame] = i;
           Rev_Start [Frame] = 0;
           Votes [3 + (3 - Frame) % 3] = 0;
          }
      if  (Is_Reverse_Start (Codon))
          Rev_Start [Frame] = i;
     }

   Initial_Base = Data_Len;
   Ignore_Seg = 0;

   if ( Ignore_Option )
     {
       for(Ignore_Seg = 0; Ignore[Ignore_Seg].Low_Address > -1; Ignore_Seg++)
	 {}
       Ignore_Seg--;
       if (Ignore[Ignore_Seg].High_Address == Data_Len)
	 Initial_Base = Ignore[i].Low_Address;
     }

   if  (! Genome_Is_Circular)
       for  (i = Initial_Base;  Initial_Base - i < 3;  i --)
         {
	   Frame = i % 3;

	   if  (i == Ignore[Ignore_Seg].High_Address)
	     {
	       i = Ignore[Ignore_Seg].Low_Address;
	       for ( j = 0; j < 3; j++ )
		 {
		   For_Prev  [j] = j;
		   Rev_Prev  [j] = j;
		   For_Start [j] = 0;
		   Rev_Start [j] = 0;
		 }
	       for ( j = 0; j < 6; j++ )
		 Votes[j] = 0;
	       Ignore_Seg--;
	       Codon = Ch_Mask (Data [i - 1]) << 4 | Ch_Mask (Data [i]);
	       i--;
	     }
	   
          if  (Rev_Start [Frame] != 0
                 && (Gene_Len = Rev_Start [Frame] - Rev_Prev [Frame])
                        >= Min_Gene_Len)
              {
               Len = i - Rev_Prev [Frame];
               Start = Choose_Start (Rev_Start [Frame], - Gene_Len);
               Gene_Len = Start - Rev_Prev [Frame];
               if  (Gene_Len >= Min_Gene_Len)
                   {
                    Score_String (Start,
                              Rev_Prev [Frame] + 1, Ch_Ct, Score,
                              Use_Independent && Len < Ignore_Indep_Len,
                              Weak_Score, Raw_Score);
                    In_Frame_Score = Score [0];
                    Is_Tentative_Gene = (In_Frame_Score >= Threshold_Score
                           || Votes [3 + (3 - Frame) % 3] + In_Frame_Score
                                  >= Vote_Threshold
                           || (Weak_Score >= Threshold_Score
                                 && Gene_Len >= Min_Weak_Len));
                    if  (In_Frame_Score >= Threshold_Score)
                        Category = REGULAR;
                    else if  (Votes [3 + (3 - Frame) % 3] + In_Frame_Score
                                  >= Vote_Threshold)
                        Category = VOTED;
                    else if  (Weak_Score >= Threshold_Score
                                 && Gene_Len >= Min_Weak_Len)
                        Category = WEAK;
                      else
                        Category = NONE;
                    if  (Is_Tentative_Gene)
                        printf ("%5ld ", ++ ID_Num);
                      else
                        printf ("%5s ", "");
                    printf (" R%1d %8ld %8ld %8ld %8ld %5ld   ",
                           1 + ((2 - Frame) + 1) % 3,
                           i, Start,
                           Rev_Prev [Frame] + 1, Len,
                           Gene_Len);
                    printf (" %4d   ", In_Frame_Score);
                    Permute (Score, - Frame - 1);
                    for  (j = 0;  j < 6;  j ++)
                      if  (Score [j] < 0)
                          printf ("  _");
                        else
                          {
                           printf (" %2d", Score [j]);
                           if  ((j < 3 && For_Prev [j] < Rev_Start [Frame])
                                 || (j >= 3
                                    && Rev_Prev [(6 - j) % 3] < Rev_Start [Frame]))
                               Votes [j] += Score [j];
                          }
                    if  (Use_Independent)
                        printf ("   %2d", Score [6]);
                    printf ("  %4ld", Votes [3 + (3 - Frame) % 3]);
                    printf ("  %6.3f", Raw_Score);
                    Print_Category (Category);
                    putchar ('\n');
                    if  (Is_Tentative_Gene)
                        Append_Gene_Ref (Gene, ID_Num, Rev_Prev [Frame] + 1,
                             Start , - Frame - 1, In_Frame_Score, Category,
                             Raw_Score);
                   }
              }
         }

   printf ("End = %ld\n", Data_Len);

   if  (! Choose_First_Start_Codon)
       {
//  Find likely ribosome binding sites and shift some start sites based
//  on it.  This is Olga's code.

        Original_Coord = (long int (*) [2])
                            Safe_malloc ((1 + ID_Num) * sizeof (long int [2]));
        Revised_Start = (long int *) Safe_malloc ((1 + ID_Num) * sizeof (long int));

        for  (i = 1;  i <= ID_Num;  i ++)
         if  (Gene [i] . Frame > 0)
             {
              Original_Coord [i] [0] = Gene [i] . Lo;
              Original_Coord [i] [1] = Gene [i] . Hi;
             }
           else
             {
              Original_Coord [i] [0] = Gene [i] . Hi;
              Original_Coord [i] [1] = Gene [i] . Lo;
             }

        RNAbin (Data, Extended_Len, Original_Coord + 1, Revised_Start + 1, ID_Num);

#if  0
{
 double  Sum = 0.0;
 long int  Ct = 0, Shift;

   for  (i = 1;  i <= ID_Num;  i ++)
     if  (Original_Coord [i] [0] != Revised_Start [i])
         {
          Shift = labs (Revised_Start [i] - Original_Coord [i] [0]); 
          printf ("%5ld:  %7ld  %7ld  %7ld  %4ld  %4ld\n", i, Original_Coord [i] [0],
                 Original_Coord [i] [1], Revised_Start [i], Shift,
                 1 + labs (Original_Coord [i] [1] - Original_Coord [i] [0]));
          Sum += Shift;
          Ct ++;
         }
   printf ("Moved %ld starts by average of %.1f bases\n",
               Ct, Sum / Ct);
}
#endif

        for  (i = 1;  i <= ID_Num;  i ++)
          if  (Original_Coord [i] [0] != Revised_Start [i])
              {
               if  (Gene [i] . Frame > 0)
                   Gene [i] . Lo = Revised_Start [i];
                 else
                   Gene [i] . Hi = Revised_Start [i];
               Gene [i] . Len = 1 + Gene [i] . Hi - Gene [i] . Lo;
               Gene [i] . Set_Status (RBS_START_SHIFT);
              }
       }

   Ptr = (Gene_Ref * *) Safe_malloc ((1 + ID_Num) * sizeof (Gene_Ref *));
   for  (i = 1;  i <= ID_Num;  i ++)
     Ptr [i] = Gene + i;


   Changes_Made = 1;
   for  (i = 1;  i <= MAX_ITERATIONS && Changes_Made > 0;  i ++)
     {
      Find_Overlaps (Gene, ID_Num);

      qsort (Ptr + 1, ID_Num, sizeof (Gene_Ref *), Gene_Ref_Cmp);
                            //  Sort potential genes by descending length

      Determine_Changes (Ptr, Gene, ID_Num);

// Global_Show_Details = TRUE;
Show_Gene_Info (Gene, ID_Num);

      Changes_Made = Make_Sure_Changes (Gene, ID_Num);
fprintf (stderr, "Changes_Made = %ld\n", Changes_Made);
fprintf (stderr, "Done iteration %2ld\n", i);

     }


   Changes_Made = 1;
   for  (i = 1;  i <= MAX_ITERATIONS && Changes_Made > 0;  i ++)
     {
      Find_Overlaps (Gene, ID_Num);

      qsort (Ptr + 1, ID_Num, sizeof (Gene_Ref *), Gene_Ref_Cmp);
                            //  Sort potential genes by descending length

      Determine_Changes (Ptr, Gene, ID_Num);

// Global_Show_Details = TRUE;
Show_Gene_Info (Gene, ID_Num);

      Changes_Made = Make_Final_Changes (Ptr, Gene, ID_Num);
fprintf (stderr, "Changes_Made = %ld\n", Changes_Made);
fprintf (stderr, "Done iteration %2ld\n", i);
     }

               //  One last time

   Find_Overlaps (Gene, ID_Num);

   qsort (Ptr + 1, ID_Num, sizeof (Gene_Ref *), Gene_Ref_Cmp);
                         //  Sort potential genes by descending length

   Determine_Changes (Ptr, Gene, ID_Num);

   printf ("\n\nPutative Genes:\n");
   for  (i = 1;  i <= ID_Num;  i ++)
     if  (! Gene [i] . Has_Status (REJECTED))
         {
          Lo = Gene [i] . Lo;
          if  (Lo > Data_Len)
              Lo -= Data_Len;
          Hi = Gene [i] . Hi;
          if  (Hi > Data_Len)
              Hi -= Data_Len;
          if  (Gene [i] . Frame > 0)
              printf ("%5ld %8ld %8ld", i, Lo, Hi);
            else
              printf ("%5ld %8ld %8ld", i, Hi, Lo);
          if  (Show_Frame_And_Length)
              printf ("  [%+2d L=%4ld r=%5.3f]", Gene [i] . Frame,
                             Gene [i] . Len, Gene [i] . Raw_Score);
          switch  (Gene [i] . Category)
            {
             case  VOTED :
               printf ("  [Vote]");
               break;
             case  WEAK :
               printf ("  [Weak]");
               break;
             default :
               ; // Nothing
            }
          for  (P = Gene [i] . Overlap_List;  P != NULL;  P = P -> Next)
            {
             if  (Gene [P -> From] . Has_Status (REJECTED))
                 {
                  fprintf (stderr, "ERROR:  Unexpected reject\n");
                  assert (FALSE);
                  exit (-1);
                 }
             switch  (P -> Problem_Code)
               {
                case  REJECT_CHAR :
                  fprintf (stderr, "ERROR:  Unexpected reject\n");
                  assert (FALSE);
                  exit (-1);
                  break;
                case  NEAR_REJECT_CHAR :
                  printf ("  [NearRejectBy #%ld L=%ld S=%d]", P -> From, P -> Olap,
                                 P -> Score);
                  break;
                case  SCORES_WORSE_CHAR :
                  printf ("  [LowScoreBy #%ld L=%ld S=%d]", P -> From, P -> Olap,
                                 P -> Score);
                  break;
                case  ELIM_OLAP_CHAR :
                  printf ("  [OlapWith #%ld L=%ld S=%d]", P -> From, P -> Olap,
                                 P -> Score);
                  break;
                case  SHORTER_CHAR :
                  printf ("  [ShorterThan #%ld L=%ld S=%d]", P -> From, P -> Olap,
                                 P -> Score);
                  break;
                case  CONTAINS_CHAR :
                  printf ("  [Contains #%ld]", P -> From);
                  break;
                case  SHADOWED_CHAR :
                  printf ("  [ShadowedBy #%ld]", P -> From);
                  break;
                case  NOPROB_CHAR :
                  break;
                default :
                  assert (P -> Problem_Code == REJECT_CHAR);
               }
            }
          if  (Gene [i] . Has_Status (RBS_START_SHIFT))
              printf ("  [RBS Start Move]");
          if  (Gene [i] . Delay_Len > 0)
              printf ("  [DelayedBy #%ld L=%ld]", Gene [i] . Delay_Cause,
                              Gene [i] . Delay_Len);
          putchar ('\n');
         }

   return  0;
  }



void  Add_Overlap
    (Gene_Ref & R, long int Other, long int Lo, long int Hi,
     char Problem, int Other_Frame)

//  Add a node to  R 's overlap list for the overlap with gene number
//  Other  at positions  Lo .. Hi .  Store  Problem  and  Other_Frame
//  in this node.  Use a simple forward list and maintain the nodes
//  in order of decreasing overlap length.

  {
   Overlap_Node  * New_Node, * P, * * Attach;

   New_Node = (Overlap_Node *) Safe_malloc (sizeof (Overlap_Node));

   New_Node -> From = Other;
   New_Node -> Lo = Lo;
   New_Node -> Hi = Hi;
   New_Node -> Olap = 1 + Hi - Lo;
   New_Node -> Score = 0;
   New_Node -> Delay = 0;
   New_Node -> Problem_Code = Problem;
   New_Node -> Other_Frame = Other_Frame;

   Attach = & (R . Overlap_List);
   for  (P = R . Overlap_List;  P != NULL && New_Node -> Olap <= P -> Olap;
                   P = P -> Next)
     Attach = & (P -> Next);

   New_Node -> Next = P;
   * Attach = New_Node;

   return;
  }



void  Append_Gene_Ref
    (Gene_Ref * & Gene, long int ID, long int Lo, long int Hi,
     int Frame, int Score, Gene_Category_Type Category, double raw_score)

//  Add reference for new gene  ID  at positions  Lo .. Hi  in
//  frame  Frame  with  Score  to  Gene .

  {
   Gene = (Gene_Ref *) Safe_realloc (Gene,
                                    (1 + ID) * sizeof (Gene_Ref));

   Gene [ID] . Lo = Lo;
   Gene [ID] . Hi = Hi;
   if  (ID > 1 && Gene [ID - 1] . Max_Hi > Hi)
       Gene [ID] . Max_Hi = Gene [ID - 1] . Max_Hi;
     else
       Gene [ID] . Max_Hi = Hi;      // Hi's may be out of order because
                                         // of start codon position.
   Gene [ID] . Len = 1 + Hi - Lo;
   Gene [ID] . Delay_Len = 0;
   Gene [ID] . Delay_Cause = 0;
   Gene [ID] . Frame = Frame;
   Gene [ID] . Status = OK;
   Gene [ID] . Overlap_List = NULL;
   Gene [ID] . Score = Score;
   Gene [ID] . Category = Category;
   Gene [ID] . Raw_Score = raw_score;
  }



double  Bulge_Cost  (int N)

/* Return the energy cost of a bulge on one side of  N  bases. */

  {
   if  (N <= 0)
       return  BIG_NEGATIVE;

   if  (N < 4)
       return  -3.3;

   return  -3.3 - (N - 3) * (15.8 - 3.3) / 27.0;
  }



int  Can_Delay_Start
    (Gene_Ref & Gene, Overlap_Node * P, char Code, long int  Min_Delay)

//  Determine whether the start of gene  Gene  can be moved downstream
//  by at least  Min_Delay  bases to eliminate the overlap in  P .
//  If so, indicate by how much, mark  Gene  as potentially changeable,
//  and return  TRUE .  Otherwise, return  FALSE .  In any case,
//  save  Code  in  * P .

  {
   char  Codon [4];
   int  Score [7], Weak_Score;
   double  Raw_Score;
   long int  i, j, Len;

   P -> Problem_Code = Code;
   if  (Gene . Frame > 0)
       {
        if  (P -> Lo != Gene . Lo)
            return  FALSE;
        i = Gene . Lo + 3 * ((Min_Delay - 1) / 3);
        Len = 1 + Gene . Hi - i;
        do
          {
           i += 3;
           Len -= 3;
           Transfer (Codon, i, 3);
          }  while  (! Is_Start (Codon) && ! Is_Stop (Codon)
                            && Len > Min_Gene_Len + 2);
        if  (! Is_Start (Codon) || Len < Min_Gene_Len)
            return  FALSE;

        j = Choose_Start (i, Len);
        Len -= j - i;
        if  (Len < Min_Gene_Len)
            return  FALSE;

        Score_String (j, Gene . Hi, Ch_Ct, Score,
                         Use_Independent && Len < Ignore_Indep_Len,
                         Weak_Score, Raw_Score);
        if  (Score [0] < Threshold_Score)
            return  FALSE;

        P -> Delay = j - Gene . Lo;
        Gene . Set_Status (MIGHT_CHANGE);
       }
     else
       {
        if  (P -> Hi != Gene . Hi)
            return  FALSE;
        i = Gene . Hi - 3 * ((Min_Delay - 1) / 3);
        Len = 1 + i - Gene . Lo;
        do
          {
           i -= 3;
           Len -= 3;
           Transfer (Codon, i, -3);
          }  while  (! Is_Start (Codon) && ! Is_Stop (Codon)
                            && Len > Min_Gene_Len + 2);
        if  (! Is_Start (Codon) || Len < Min_Gene_Len)
            return  FALSE;

        j = Choose_Start (i, - Len);
        Len -= i - j;
        if  (Len < Min_Gene_Len)
            return  FALSE;

        Score_String (i, Gene . Lo, Ch_Ct, Score,
                         Use_Independent && Len < Ignore_Indep_Len,
                         Weak_Score, Raw_Score);
        if  (Score [0] < Threshold_Score)
            return  FALSE;

        P -> Delay = Gene . Hi - j;
        Gene . Set_Status (MIGHT_CHANGE);
       }
   return  TRUE;
  }



int  Choose_Score  (int S [7], int F)

//  Return the score in  S  corresponding to frame  F .

  {
   switch  (F)
     {
      case  1 :
        return  S [0];
      case  2 :
        return  S [1];
      case  3 :
        return  S [2];
      case  -1 :
        return  S [3];
      case  -2 :
        return  S [5];
      case  -3 :
        return  S [4];
      default :
        fprintf (stderr, "ERROR:  Bad frame value  %d  in  Choose_Score ()\n", F);
     }

   return  S [0];
  }



long int  Choose_Start  (long int Begin, long int Len)

/* Return the position of the first base of the most likely start
*  codon for the gene whose first start codon is at base  Begin
*  and has length  Len , which is positive in the forward direction
*  and negative in the reverse complement direction. */

  {
   double  Adjustment, Score;
   long int  Min_Len, Original_Start;
   char  Buffer [1 + UPSTREAM_LEN], Codon [4];

   if  (Choose_First_Start_Codon)
       return  Begin;

   Original_Start = Begin;
   Codon [3] = '\0';
   Min_Len = Max (Min_Gene_Len, Len - MAX_START_SHIFT);

   if  (Len > 0)
       {
        Transfer (Codon, Begin, 3);
        switch  (Codon [0])
          {
           case  'a' :
             Adjustment = 1.0;
             break;
           case  'g' :
             Adjustment = 0.0;
             break;
           case  't' :
             Adjustment = -1.0;
             break;
           default :
             Adjustment = -1.0;
          }
        do
          {
           Transfer (Buffer, Begin - UPSTREAM_LEN - UPSTREAM_OFFSET, UPSTREAM_LEN);
           Score = Edit_Distance (Ribosome_Pattern - 1, Buffer - 1);
           switch  (Codon [0])
             {
              case  'a' :
                if  (Score >= ATG_THRESHOLD + Adjustment)
                    return  Begin;
                break;
              case  'g' :
                if  (Score >= GTG_THRESHOLD + Adjustment)
                    return  Begin;
                break;
              case  't' :
                if  (Score >= TTG_THRESHOLD + Adjustment)
                    return  Begin;
                break;
             }
           do
             {
              Begin += 3;
              Len -= 3;
              if  (Begin > Data_Len)
                  Begin -= Data_Len;
              Transfer (Codon, Begin, 3);
             }  while  (! Is_Start (Codon) && ! Is_Stop (Codon)
                               && Len > Min_Len);
          }  while  (! Is_Stop (Codon) && Len > Min_Len);
       }
     else
       {
        Transfer (Codon, Begin, -3);
        switch  (Codon [0])
          {
           case  'a' :
             Adjustment = 1.0;
             break;
           case  'g' :
             Adjustment = 0.0;
             break;
           case  't' :
             Adjustment = -1.0;
             break;
           default :
             Adjustment = -1.0;
          }
        do
          {
           Transfer (Buffer, Begin + UPSTREAM_LEN + UPSTREAM_OFFSET, - UPSTREAM_LEN);
           Score = Edit_Distance (Ribosome_Pattern - 1, Buffer - 1);
           switch  (Codon [0])
             {
              case  'a' :
                if  (Score >= ATG_THRESHOLD + Adjustment)
                    return  Begin;
                break;
              case  'g' :
                if  (Score >= GTG_THRESHOLD + Adjustment)
                    return  Begin;
                break;
              case  't' :
                if  (Score >= TTG_THRESHOLD + Adjustment)
                    return  Begin;
                break;
             }
           do
             {
              Begin -= 3;
              Len -= 3;
              if  (Begin < 1)
                  Begin += Data_Len;
              Transfer (Codon, Begin, -3);
             }  while  (! Is_Start (Codon) && ! Is_Stop (Codon)
                               && Len > Min_Len);
          }  while  (! Is_Stop (Codon) && Len > Min_Len);
       }

   return  Original_Start;
  }


int  Cmp_Ignore (const void * A, const void * B)

// Comparison function for qsort to sort Ignore_Nodes by
// assending base pair number in the lower address

  {
   
    if  ( ((Ignore_Ptr) A) -> Low_Address < ((Ignore_Ptr) B) -> Low_Address)
      return  -1;
    else if  ( ((Ignore_Ptr) A) -> Low_Address > ((Ignore_Ptr) B) -> Low_Address)
      return  1;
    else
      return  0;
  }


void  Determine_Changes
    (Gene_Ref * Ptr [], Gene_Ref Gene [], long int N)

//  Consider all overlaps in  Gene [1 .. N]  in order by
//  pointers in  Ptr [1 .. N] .  Based on them reject or
//  annotate entries in  Gene .  If a gene's start site is
//  delayed, then possibly shift its entry in  Ptr .

  {
   Overlap_Node  * P;
   long int  i, Sub;

   for  (i = 1;  i <= N;  i ++)
     {
      if  (Ptr [i] -> Has_Status (REJECTED))
          continue;
      Sub = Ptr [i] - Gene;
      for  (P = Ptr [i] -> Overlap_List;  P != NULL;  P = P -> Next)
        {
         if  (Gene [P -> From] . Has_Status (REJECTED))
             continue;
         if  (P -> Problem_Code == SHADOWED_CHAR
               || P -> Problem_Code == CONTAINS_CHAR)
             continue;   // Analyze subwindows here???
         if  (Gene [Sub] . Len < Gene [P -> From] . Len
                || (Gene [Sub] . Len == Gene [P -> From] . Len
                       && Sub < P -> From))
             continue;
         Evaluate_Overlap (P, Sub, Gene);
        }
     }

   return;
  }



double  Doublet_Score  (char A, char B, char P, char Q)

/* Return the energy released by consecutive bases  AB  binding
*  to pair  PQ  where  PQ  is in  5'-3' order.   Values from
*  Lewin's Genes V, p. 115. */

  {
   B = tolower (B);
   P = tolower (P);
   Q = tolower (Q);

   switch  (tolower (A))
     {
      case  'a' :
        switch  (B)
          {
           case  'a' :
             if  (P == 't' && Q == 't')
                 return  0.9;
             goto  Error;
           case  'c' :
             if  (P == 't' && Q == 'g')
                 return  1.8;
             goto  Error;
           case  'g' :
             if  (P == 't' && Q == 'c')
                 return  2.3;
             else if  (P == 't' && Q == 't')
                 return  1.15;                    /* Guess */
             goto  Error;
           case  't' :
             if  (P == 't' && Q == 'a')
                 return  1.1;
             else if  (P == 't' && Q == 'g')
                 return  0.65;                    /* Guess */
             goto  Error;
          }

      case  'c' :
        switch  (B)
          {
           case  'a' :
             if  (P == 'g' && Q == 't')
                 return  2.1;
             goto  Error;
           case  'c' :
             if  (P == 'g' && Q == 'g')
                 return  2.9;
             goto  Error;
           case  'g' :
             if  (P == 'g' && Q == 'c')
                 return  3.4;
             else if  (P == 'g' && Q == 't')
                 return  1.50;                    /* Guess */
             goto  Error;
           case  't' :
             if  (P == 'g' && Q == 'a')
                 return  2.3;
             else if  (P == 'g' && Q == 'g')
                 return  1.15;                    /* Guess */
             goto  Error;
          }

      case  'g' :
        switch  (B)
          {
           case  'a' :
             if  (P == 'c' && Q == 't')
                 return  1.7;
             else if  (P == 't' && Q == 't')
                 return  0.85;                    /* Guess */
             goto  Error;
           case  'c' :
             if  (P == 'c' && Q == 'g')
                 return  2.0;
             else if  (P == 't' && Q == 'g')
                 return  1.00;                    /* Guess */
             goto  Error;
           case  'g' :
             if  (P == 'c' && Q == 'c')
                 return  2.9;
             else if  (P == 't' && Q == 'c')
                 return  1.45;                    /* Guess */
             else if  (P == 'c' && Q == 't')
                 return  1.45;                    /* Guess */
             else if  (P == 't' && Q == 't')
                 return  0.5;                     /* Guess */
             goto  Error;
           case  't' :
             if  (P == 'c' && Q == 'a')
                 return  1.8;
             else if  (P == 't' && Q == 'a')
                 return  0.9;                     /* Guess */
             else if  (P == 'c' && Q == 'g')
                 return  0.9;                     /* Guess */
             else if  (P == 't' && Q == 'g')
                 return  0.5;                     /* Guess */
             goto  Error;
          }

      case  't' :
        switch  (B)
          {
           case  'a' :
             if  (P == 'a' && Q == 't')
                 return  0.9;
             else if  (P == 'g' && Q == 't')
                 return  0.5;                     /* Guess */
             goto  Error;
           case  'c' :
             if  (P == 'a' && Q == 'g')
                 return  1.7;
             else if  (P == 'g' && Q == 'g')
                 return  0.85;                    /* Guess */
             goto  Error;
           case  'g' :
             if  (P == 'a' && Q == 'c')
                 return  2.1;
             else if  (P == 'g' && Q == 'c')
                 return  1.05;                    /* Guess */
             else if  (P == 'a' && Q == 't')
                 return  1.05;                    /* Guess */
             else if  (P == 'g' && Q == 't')
                 return  0.5;                     /* Guess */
             goto  Error;
           case  't' :
             if  (P == 'a' && Q == 'a')
                 return  0.9;
             else if  (P == 'g' && Q == 'a')
                 return  0.5;                     /* Guess */
             else if  (P == 'a' && Q == 'g')
                 return  0.5;                     /* Guess */
             else if  (P == 'g' && Q == 'g')
                 return  0.5;                     /* Guess */
             goto  Error;
          }

     }
   
  Error:
   fprintf (stderr, "ERROR:  Bad doublet pair  %c%c  and  %c%c\n",
              A, B, P, Q);

   return  0;
  }



double  Edit_Distance  (const char * P, const char * T)

/* Find and return the highest enery match of string  P [1 ...]  within
*  string  T [1 ...] . */

  {
   double  Max, Best, Final_Score, X;
   int  Best_i = 0, Best_j = 0;
   int  a, b, i, j, M, N, Len, Len1, Max_i = 0;

   M = strlen (P + 1);
   N = strlen (T + 1);

   for  (i = 1;  i <= M;  i ++)
     {
      ED_Score [i] [0] . Free_i = 0.0;
      ED_Score [i] [0] . Free_j = BIG_NEGATIVE;
      ED_Score [i] [0] . Both_Free = BIG_NEGATIVE;
      ED_Score [i] [0] . Match = BIG_NEGATIVE;
      ED_Score [i] [0] . Free_i_Len = i;
      ED_Score [i] [0] . Free_j_Len = 0;
      ED_Score [i] [0] . Both_i_Len = 0;
      ED_Score [i] [0] . Both_j_Len = 0;
     }
   for  (j = 1;  j <= N;  j ++)
     {
      ED_Score [0] [j] . Free_i = BIG_NEGATIVE;
      ED_Score [0] [j] . Free_j = 0.0;
      ED_Score [0] [j] . Both_Free = BIG_NEGATIVE;
      ED_Score [0] [j] . Match = BIG_NEGATIVE;
      ED_Score [0] [j] . Free_i_Len = 0;
      ED_Score [0] [j] . Free_j_Len = j;
      ED_Score [0] [j] . Both_i_Len = 0;
      ED_Score [0] [j] . Both_j_Len = 0;
     }
   ED_Score [0] [0] . Free_i = 0.0;
   ED_Score [0] [0] . Free_j = 0.0;
   ED_Score [0] [0] . Both_Free = 0.0;
   ED_Score [0] [0] . Match = BIG_NEGATIVE;
   ED_Score [0] [0] . Free_i_Len = 0;
   ED_Score [0] [0] . Free_j_Len = 0;
   ED_Score [0] [0] . Both_i_Len = 0;
   ED_Score [0] [0] . Both_j_Len = 0;

   Final_Score = BIG_NEGATIVE;
   for  (j = 1;  j <= N;  j ++)
     {
      Max = BIG_NEGATIVE;
      for  (i = 1;  i <= M;  i ++)
        {
         Best = BIG_NEGATIVE;
         Len = 0;
         X = ED_Score [i - 1] [j] . Match;
         if  (X > Best)
             {
              Best = X;
              Len = 1;
             }
         X = ED_Score [i - 1] [j] . Free_i;
         a = ED_Score [i - 1] [j] . Free_i_Len;
         if  (X > Best && a < MAX_FREE_LEN)
             {
              Best = X;
              Len = 1 + a;
             }
         ED_Score [i] [j] . Free_i = Best;
         ED_Score [i] [j] . Free_i_Len = Len;
         if  (Best > Max)
             {
              Max = Best;
              Max_i = i;
             }
         
         Best = BIG_NEGATIVE;
         Len = 0;
         X = ED_Score [i] [j - 1] . Match;
         if  (X > Best)
             {
              Best = X;
              Len = 1;
             }
         X = ED_Score [i] [j - 1] . Free_j;
         a = ED_Score [i] [j - 1] . Free_j_Len;
         if  (X > Best && a < MAX_FREE_LEN)
             {
              Best = X;
              Len = 1 + a;
             }
         ED_Score [i] [j] . Free_j = Best;
         if  (Best > BIG_NEGATIVE)
             ED_Score [i] [j] . Free_j_Len = Len;
           else
             ED_Score [i] [j] . Free_j_Len = 0;

         Best = BIG_NEGATIVE;
         Len = Len1 = 0;
         X = ED_Score [i - 1] [j - 1] . Match;
         if  (X > Best)
             {
              Best = X;
              Len = Len1 = 1;
             }
         X = ED_Score [i - 1] [j] . Free_j;
         if  (X > Best)
             {
              Best = X;
              Len = 1;
              Len1 = ED_Score [i - 1] [j] . Free_j_Len;
             }
         X = ED_Score [i] [j - 1] . Free_i;
         if  (X > Best)
             {
              Best = X;
              Len = ED_Score [i] [j - 1] . Free_i_Len;
              Len1 = 1;
             }
         X = ED_Score [i - 1] [j] . Both_Free;
         a = ED_Score [i - 1] [j] . Both_i_Len;
         if  (X > Best && a < MAX_FREE_LEN)
             {
              Best = X;
              Len = a + 1;
              Len1 = ED_Score [i - 1] [j] . Both_j_Len;
             }
         X = ED_Score [i] [j - 1] . Both_Free;
         a = ED_Score [i] [j - 1] . Both_j_Len;
         if  (X > Best && a < MAX_FREE_LEN)
             {
              Best = X;
              Len = ED_Score [i] [j - 1] . Both_i_Len;
              Len1 = a + 1;
             }
         X = ED_Score [i - 1] [j - 1] . Both_Free;
         a = ED_Score [i - 1] [j - 1] . Both_i_Len;
         b = ED_Score [i - 1] [j - 1] . Both_j_Len;
         if  (X > Best && a < MAX_FREE_LEN && b < MAX_FREE_LEN)
             {
              Best = X;
              Len = a + 1;
              Len1 = b + 1;
             }
         ED_Score [i] [j] . Both_Free = Best;
         ED_Score [i] [j] . Both_i_Len = Len;
         ED_Score [i] [j] . Both_j_Len = Len1;

         if  (! Match (P [i], T [j]))
             ED_Score [i] [j] . Match = BIG_NEGATIVE;
           else
             {
              Best = 0.0;
              X = ED_Score [i - 1] [j - 1] . Free_i
                      + Bulge_Cost (ED_Score [i - 1] [j - 1] . Free_i_Len);
              if  (X > Best)
                  Best = X;
              X = ED_Score [i - 1] [j - 1] . Free_j
                      + Bulge_Cost (ED_Score [i - 1] [j - 1] . Free_j_Len);
              if  (X > Best)
                  Best = X;
              X = ED_Score [i - 1] [j - 1] . Both_Free
                      + Loop_Cost (ED_Score [i - 1] [j - 1] . Both_i_Len,
                                      ED_Score [i - 1] [j - 1] . Both_j_Len);
              if  (X > Best)
                  Best = X;
              X = ED_Score [i - 1] [j - 1] . Match;
              if  (X != BIG_NEGATIVE)
                  {
                   X += Doublet_Score (P [i - 1], P [i], T [j - 1], T [j]);
                   if  (X > Best)
                       Best = X;
                  }
              ED_Score [i] [j] . Match = Best;
              if  (Best > Max)
                  {
                   Max = Best;
                   Max_i = i;
                  }
             }
        }
      if  (Max > Final_Score)
          {
           Final_Score = Max;
           Best_i = Max_i;
           Best_j = j;
          }
     }

   return  Final_Score;
  }



void  Evaluate_Overlap  (Overlap_Node * A_Olap_Node, long int A, Gene_Ref Gene [])

//  Analyze the overlap pointed to by  A_Olap_Node  with  Gene [A]  and
//  save the results in  (* A_Olap_Node)  and  Gene [A] .

  {
   Overlap_Node  * B_Olap_Node;
   Olap_Fix_Type  Olap_Fix;
   long int  A_Lo, A_Hi, B_Lo, B_Hi;
   long int  B, Min_Shift, Percent_Min;
   int  A_Score, B_Score;

   B = A_Olap_Node -> From;
   Score_Olap_Region (A_Olap_Node -> Lo, A_Olap_Node -> Hi,
                         Gene [A] . Frame, A_Olap_Node -> Other_Frame,
                         A_Score, B_Score);

   for  (B_Olap_Node = Gene [B] . Overlap_List;
            B_Olap_Node != NULL && B_Olap_Node -> From != A;
               B_Olap_Node = B_Olap_Node -> Next)
     ;
   assert (B_Olap_Node != NULL);

   A_Olap_Node -> Score = A_Score;
   B_Olap_Node -> Score = B_Score;

   Percent_Min = (long int) ceil (Min_Olap_Percent *
        Min (Gene [A] . Len, Gene [B] . Len));
   Min_Shift = 1 + A_Olap_Node -> Olap - Max (Min_Olap, Percent_Min);
   
   if  (Min_Shift <= 0 || Min_Shift > A_Olap_Node -> Olap)
       {
        fprintf (stderr, "ERROR:  Shouldn't be here.  Overlap is ignorable\n");
        fprintf (stderr, "Min_Shift = %ld for genes %ld and %ld   Olap = %ld\n",
                            Min_Shift, A, B, A_Olap_Node -> Olap);
        fprintf (stderr, "Percent_Min = %ld  A_Len = %ld   B_Len = %ld\n",
                            Percent_Min, Gene [A] . Len, Gene [B] . Len);
        exit (-3);
       }

   assert (Gene [A] . Len >= Gene [B] . Len);

   if  (Gene [A] . Category == REGULAR
          && (Gene [B] . Category == WEAK
                || Gene [B] . Category == VOTED))
       {
        Gene [B] . Set_Status (MIGHT_CHANGE);
        B_Olap_Node -> Problem_Code = REJECT_CHAR;
        return;
       }
   if  (Gene [B] . Category == REGULAR
          && (Gene [A] . Category == WEAK
                || Gene [A] . Category == VOTED))
       {
        Gene [A] . Set_Status (MIGHT_CHANGE);
        A_Olap_Node -> Problem_Code = REJECT_CHAR;
        return;
       }

   A_Lo = Gene [A] . Lo;
   A_Hi = Gene [A] . Hi;
   B_Lo = Gene [B] . Lo;
   B_Hi = Gene [B] . Hi;
   if  (A_Hi > Data_Len && B_Hi < A_Lo)
       {
        B_Lo += Data_Len;
        B_Hi += Data_Len;
       }
   else if  (B_Hi > Data_Len && A_Hi < B_Lo)
       {
        A_Lo += Data_Len;
        A_Hi += Data_Len;
       }

   Olap_Fix = Get_Olap_Fix (A_Lo, A_Hi, Gene [A] . Frame,
                            B_Lo, B_Hi, Gene [B] . Frame);

// Approximately equal lengths

   if  (Gene [A] . Len - Gene [B] . Len < MIN_PERCENT_LEN_DIFF * Gene [B] . Len)
       {
        if  (A_Score > B_Score && A_Score >= OLAP_THRESHOLD_SCORE)
            {
             switch  (Olap_Fix)
               {
                case  NEITHER_CAN_MOVE :
                  A_Olap_Node -> Problem_Code = ELIM_OLAP_CHAR;
                  B_Olap_Node -> Problem_Code = SCORES_WORSE_CHAR;
                  return;
                case  ONLY_A_CAN_MOVE :
                  if  (Min_Shift <= SMALL_OLAP_PERCENT * Gene [A] . Len
                         && Can_Delay_Start (Gene [A], A_Olap_Node,
                                             ELIM_OLAP_CHAR, Min_Shift))
                      return;
                  A_Olap_Node -> Problem_Code = ELIM_OLAP_CHAR;
                  B_Olap_Node -> Problem_Code = SCORES_WORSE_CHAR;
                  return;
                case  ONLY_B_CAN_MOVE :
                  if  (Can_Delay_Start (Gene [B], B_Olap_Node,
                                             SCORES_WORSE_CHAR, Min_Shift))
                      return;
                  A_Olap_Node -> Problem_Code = ELIM_OLAP_CHAR;
                  B_Olap_Node -> Problem_Code = NEAR_REJECT_CHAR;
                  return;
                case  BOTH_CAN_MOVE :
                  Slide_Both_Starts (Gene [B], B_Olap_Node, NEAR_REJECT_CHAR,
                                     Gene [A], A_Olap_Node, Min_Shift);
                  return;
                default :
                  assert (FALSE);
               }
            }
        else if  (B_Score > A_Score && B_Score >= OLAP_THRESHOLD_SCORE)
            {
             switch  (Olap_Fix)
               {
                case  NEITHER_CAN_MOVE :
                  A_Olap_Node -> Problem_Code = SCORES_WORSE_CHAR;
                  B_Olap_Node -> Problem_Code = ELIM_OLAP_CHAR;
                  return;
                case  ONLY_A_CAN_MOVE :
                  if  (Can_Delay_Start (Gene [A], A_Olap_Node,
                                             SCORES_WORSE_CHAR, Min_Shift))
                      return;
                  A_Olap_Node -> Problem_Code = SCORES_WORSE_CHAR;
                  B_Olap_Node -> Problem_Code = ELIM_OLAP_CHAR;
                  return;
                case  ONLY_B_CAN_MOVE :
                  if  (Min_Shift <= SMALL_OLAP_PERCENT * Gene [B] . Len
                         && Can_Delay_Start (Gene [B], B_Olap_Node,
                                             ELIM_OLAP_CHAR, Min_Shift))
                      return;
                  A_Olap_Node -> Problem_Code = NEAR_REJECT_CHAR;
                  B_Olap_Node -> Problem_Code = ELIM_OLAP_CHAR;
                  return;
                case  BOTH_CAN_MOVE :
                  Slide_Both_Starts (Gene [A], A_Olap_Node, NEAR_REJECT_CHAR,
                                     Gene [B], B_Olap_Node, Min_Shift);
                  return;
                default :
                  assert (FALSE);
               }
            }
          else    // Scores are essentially tied
            {
             A_Olap_Node -> Problem_Code = ELIM_OLAP_CHAR;
             B_Olap_Node -> Problem_Code = ELIM_OLAP_CHAR;
             switch  (Olap_Fix)
               {
                case  NEITHER_CAN_MOVE :
                  return;
                case  ONLY_A_CAN_MOVE :
                  if  (Can_Delay_Start (Gene [A], A_Olap_Node,
                                             ELIM_OLAP_CHAR, Min_Shift))
                      return;
                  return;
                case  ONLY_B_CAN_MOVE :
                  Can_Delay_Start (Gene [B], B_Olap_Node,
                                   ELIM_OLAP_CHAR, Min_Shift);
                  return;
                case  BOTH_CAN_MOVE :
                  Slide_Both_Starts (Gene [B], B_Olap_Node, ELIM_OLAP_CHAR,
                                     Gene [A], A_Olap_Node, Min_Shift);
                  return;
                default :
                  assert (FALSE);
               }
            }
        return;
       }

//  A is significantly longer

   if  (A_Score > B_Score && A_Score >= OLAP_THRESHOLD_SCORE)
       {
        Gene [B] . Set_Status (MIGHT_CHANGE);
        switch  (Olap_Fix)
          {
           case  NEITHER_CAN_MOVE :
             B_Olap_Node -> Problem_Code = REJECT_CHAR;
             return;
           case  ONLY_A_CAN_MOVE :
             if  (Min_Shift <= SMALL_OLAP_PERCENT * Gene [A] . Len
                    && Gene [B] . Len >= MIN_LONG_GENE_LEN
                    && Can_Delay_Start (Gene [A], A_Olap_Node,
                                        ELIM_OLAP_CHAR, Min_Shift))
                 return;
             B_Olap_Node -> Problem_Code = REJECT_CHAR;
             return;
           case  ONLY_B_CAN_MOVE :
             if  (Can_Delay_Start (Gene [B], B_Olap_Node,
                                   ELIM_OLAP_CHAR, Min_Shift))
                 return;
             B_Olap_Node -> Problem_Code = REJECT_CHAR;
             return;
           case  BOTH_CAN_MOVE :
             Slide_Both_Starts (Gene [B], B_Olap_Node, REJECT_CHAR,
                                Gene [A], A_Olap_Node, Min_Shift);
             return;
           default :
             assert (FALSE);
          }
       }
   else if  (B_Score > A_Score && B_Score >= OLAP_THRESHOLD_SCORE)
       {
        A_Olap_Node -> Problem_Code = SCORES_WORSE_CHAR;
        B_Olap_Node -> Problem_Code = SHORTER_CHAR;
        switch  (Olap_Fix)
          {
           case  NEITHER_CAN_MOVE :
             return;
           case  ONLY_A_CAN_MOVE :
             Can_Delay_Start (Gene [A], A_Olap_Node,
                              SCORES_WORSE_CHAR, Min_Shift);
             return;
           case  ONLY_B_CAN_MOVE :
             if  (Min_Shift <= SMALL_OLAP_PERCENT * Gene [B] . Len)
                 Can_Delay_Start (Gene [B], B_Olap_Node,
                                  SHORTER_CHAR, Min_Shift);
             return;
           case  BOTH_CAN_MOVE :
             Slide_Both_Starts (Gene [A], A_Olap_Node, SCORES_WORSE_CHAR,
                                Gene [B], B_Olap_Node, Min_Shift);
             return;
           default :
             assert (FALSE);
          }
       }
     else
       {
        A_Olap_Node -> Problem_Code = ELIM_OLAP_CHAR;
        B_Olap_Node -> Problem_Code = SHORTER_CHAR;
        switch  (Olap_Fix)
          {
           case  NEITHER_CAN_MOVE :
             return;
           case  ONLY_A_CAN_MOVE :
             Can_Delay_Start (Gene [A], A_Olap_Node,
                              ELIM_OLAP_CHAR, Min_Shift);
             return;
           case  ONLY_B_CAN_MOVE :
             if  (Min_Shift <= SMALL_OLAP_PERCENT * Gene [B] . Len)
                 Can_Delay_Start (Gene [B], B_Olap_Node,
                                  SHORTER_CHAR, Min_Shift);
             return;
           case  BOTH_CAN_MOVE :
             Slide_Both_Starts (Gene [A], A_Olap_Node, ELIM_OLAP_CHAR,
                                Gene [B], B_Olap_Node, Min_Shift);
             return;
           default :
             assert (FALSE);
          }
       }

   return;
  }



long int  Extend_Data  (char * & Data, long int Data_Len, long int Max_Extend)

//  Allocate additional memory at the end of  Data [1 .. Len]
//  and duplicate a long enough prefix of it to ensure a
//  stop codon in all 6 reading frames.  Return the length of
//  the extended version.  This will allow processing without
//  calculating wraparounds.

  {
   unsigned  Codon;
   int  Frame, Found_Stop [6] = {0}, Num_Stops_Found = 0;
   long int  i, j, New_Len;

   Codon = Ch_Mask (Data [Data_Len - 1]) << 4 | Ch_Mask (Data [Data_Len]);
   Frame = 0;

   for  (i = 1;  Num_Stops_Found < 6 && i <= Data_Len;  i ++)
     {
      Codon = (Codon & SHIFT_MASK) << 4;
      Codon |= Ch_Mask (Data [i]);
      Frame = (Frame + 1) % 3;
      if  (Is_Forward_Stop (Codon))
          {
           if  (! Found_Stop [Frame])
               {
                Found_Stop [Frame] = TRUE;
                Num_Stops_Found ++;
               }
          }
      if  (Is_Reverse_Stop (Codon))
          {
           if  (! Found_Stop [3 + Frame])
               {
                Found_Stop [3 + Frame] = TRUE;
                Num_Stops_Found ++;
               }
          }
     }

   if (i - 1 > Max_Extend)
     New_Len = Data_Len + Max_Extend;
   else
     New_Len = Data_Len + i - 1;

   Data = (char *) Safe_realloc (Data, 2 + New_Len);
   for  (j = 1;  j < i;  j ++)
     Data [Data_Len + j] = Data [j];
   Data [Data_Len + i] = '\0';

   if  (i == Data_Len)
       {
        //  No stops found in string at all.  Shouldn't really happen.

        New_Len += 11;
        Data = (char *) Safe_realloc (Data, 2 + New_Len);
        strcat (Data, "tagctagctag");  // Ensure a stop codon in each frame.
       }

   return  New_Len;
  }



void  Find_Overlaps  (Gene_Ref * Gene, long int N)

//  Find all overlaps between genes in  Gene [1 .. N]  and
//  put entries for them on the  Overlap_List 's for each gene.

  {
   Overlap_Node  * Save;
   char  Prob_i, Prob_j;
   int  Other_Frame;
   long int  i, j, Len_i, Len_j, Olap_Hi, Olap_Lo, Olap_Len;

   Gene [N] . Min_Lo = Gene [N] . Lo;
   for  (i = N - 1;  i > 0;  i --)
     if  (Gene [i] . Lo < Gene [i + 1] . Min_Lo)
         Gene [i] . Min_Lo = Gene [i] . Lo;
       else
         Gene [i] . Min_Lo = Gene [i + 1] . Min_Lo;

   for  (i = 1;  i <= N;  i ++)
     while  (Gene [i] . Overlap_List != NULL)
       {
        Save = Gene [i] . Overlap_List;
        Gene [i] . Overlap_List = Save -> Next;
        free (Save);
       }

   for  (i = 2;  i <= N;  i ++)
     {
      if  (Gene [i] . Has_Status (REJECTED))
          continue;
      Gene [i] . Status = OK;
      Len_i = 1 + Gene [i] . Hi - Gene [i] . Lo;
      for  (j = i - 1;  j > 0 && Gene [i] . Lo <= Gene [j] . Max_Hi;  j --)
        if  (! Gene [j] . Has_Status (REJECTED)
                 && Gene [j] . Lo <= Gene [i] . Hi
                 && Gene [i] . Lo <= Gene [j] . Hi)
            {
             Len_j = 1 + Gene [j] . Hi - Gene [j] . Lo;
             Olap_Lo = Max (Gene [i] . Lo, Gene [j] . Lo);
             Olap_Hi = Min (Gene [i] . Hi, Gene [j] . Hi);
             Olap_Len = 1 + Olap_Hi - Olap_Lo;
             if  (Olap_Lo == Gene [i] . Lo && Olap_Hi == Gene [i] . Hi)
                 {
                  Prob_i = SHADOWED_CHAR;
                  Prob_j = CONTAINS_CHAR;
                 }
             else if  (Olap_Lo == Gene [j] . Lo && Olap_Hi == Gene [j] . Hi)
                 {
                  Prob_i = CONTAINS_CHAR;
                  Prob_j = SHADOWED_CHAR;
                 }
               else
                 Prob_i = Prob_j = NOPROB_CHAR;
             if  (Olap_Len >= Min_Olap
                   && (Olap_Len >= Min_Olap_Percent * Len_i
                         || Olap_Len >= Min_Olap_Percent * Len_j))
                 {
                  Add_Overlap (Gene [j], i, Olap_Lo, Olap_Hi, Prob_j,
                               Gene [i] . Frame);
                  Add_Overlap (Gene [i], j, Olap_Lo, Olap_Hi, Prob_i,
                               Gene [j] . Frame);
                 }
            }
      if  (Gene [i] . Hi > Data_Len)    // Check wraparounds
          for  (j = 1;  j < i
                    && Data_Len + Gene [j] . Min_Lo <= Gene [i] . Hi;  j ++)
            if  (! Gene [j] . Has_Status (REJECTED)
                     && Data_Len + Gene [j] . Lo <= Gene [i] . Hi
                     && Gene [i] . Lo <= Data_Len + Gene [j] . Hi)
                {
                 Len_j = 1 + Gene [j] . Hi - Gene [j] . Lo;
                 Olap_Lo = Max (Gene [i] . Lo, Data_Len + Gene [j] . Lo);
                 Olap_Hi = Min (Gene [i] . Hi, Data_Len + Gene [j] . Hi);
                 Olap_Len = 1 + Olap_Hi - Olap_Lo;
                 if  (Olap_Lo == Gene [i] . Lo && Olap_Hi == Gene [i] . Hi)
                     {
                      Prob_i = SHADOWED_CHAR;
                      Prob_j = CONTAINS_CHAR;
                     }
                 else if  (Olap_Lo == Data_Len + Gene [j] . Lo
                             && Olap_Hi == Data_Len + Gene [j] . Hi)
                     {
                      Prob_i = CONTAINS_CHAR;
                      Prob_j = SHADOWED_CHAR;
                     }
                   else
                     Prob_i = Prob_j = NOPROB_CHAR;
                 if  (Olap_Len >= Min_Olap
                       && (Olap_Len >= Min_Olap_Percent * Len_i
                             || Olap_Len >= Min_Olap_Percent * Len_j))
                     {
                      Other_Frame = (Gene [i] . Hi - Data_Len) % 3;
                      if  (Gene [i] . Frame > 0)
                          Other_Frame = 1 + Other_Frame;
                        else
                          Other_Frame = - Other_Frame - 1;
                      Add_Overlap (Gene [j], i, Olap_Lo - Data_Len,
                                   Olap_Hi - Data_Len, Prob_j,
                                   Other_Frame);
                      Other_Frame = (Gene [j] . Hi + Data_Len) % 3;
                      if  (Gene [j] . Frame > 0)
                          Other_Frame = 1 + Other_Frame;
                        else
                          Other_Frame = - Other_Frame - 1;
                      Add_Overlap (Gene [i], j, Olap_Lo, Olap_Hi,
                                   Prob_i, Other_Frame);
                     }
                }
     }

   return;
  }



int  Gene_Ref_Cmp  (const void * A, const void * B)

/* Comparison function for  qsort  to sort  Gene_Refs  by
*  descending length. */

  {
   long int  A_Len, B_Len;

   A_Len = 1 + (* ((Gene_Ref * *) A)) -> Hi - (* ((Gene_Ref * *) A)) -> Lo;
   B_Len = 1 + (* ((Gene_Ref * *) B)) -> Hi - (* ((Gene_Ref * *) B)) -> Lo;

   if  (A_Len > B_Len)
       return  -1;
   else if  (A_Len < B_Len)
       return  1;
     else
       return  0;
  }



void  Indep_Eval  (char X [], int T, double P [], double & Prob_X)

//  Set  Prob_X  to the log of the probability of generating DNA string
//  X [1 .. T]  using the independent logs of probabilities of single
//  characters in  P [] .

  {
   const double  MIN_LOG_PROB_FACTOR = -6.0;
   double  fwd_score, rev_score, new_score;
   int  i, j, frame;

   if  (Use_Strict_Independent)
       {  // New, stronger independent model, fewer predictions
        new_score = MIN_LOG_PROB_FACTOR * T;
        for  (frame = 0;  frame < 3;  frame ++)
          {
           fwd_score = rev_score = 0.0;
           for  (j = 1;  j <= frame;  j ++)
             {
              fwd_score += P [Nucleotide_To_Subscript (X [j])];
              rev_score += P [Rev_Nucleotide_To_Subscript (X [j])];
             }
           for  (i = 1 + frame;  i <= T - 2;  i += 3)
             {
              fwd_score += Codon_Log_Prob [Codon_To_Subscript (X + i)];
              rev_score += Codon_Log_Prob [Rev_Codon_To_Subscript (X + i)];
             }
           for  (j = i;  j <= T;  j ++)
             {
              fwd_score += P [Nucleotide_To_Subscript (X [j])];
              rev_score += P [Rev_Nucleotide_To_Subscript (X [j])];
             }
           if  (fwd_score > new_score)
               new_score = fwd_score;
           if  (rev_score > new_score)
               new_score = rev_score;
          }
        Prob_X = new_score;
       }
     else
       {  // Old, weaker independent model, more predictions
        Prob_X = 0.0;
        for  (i = 1;  i <= T;  i ++)
          switch  (X [i])
            {
             case  'a' :
               Prob_X += P [0];
               break;
             case  'c' :
               Prob_X += P [1];
               break;
             case  'g' :
               Prob_X += P [2];
               break;
             case  't' :
               Prob_X += P [3];
               break;
            }
        Prob_X = Max (Prob_X, MIN_LOG_PROB_FACTOR * T);
       }

   return;
  }



Olap_Fix_Type  Get_Olap_Fix
    (long int A_Lo, long int A_Hi, int A_Frame,
     long int B_Lo, long int B_Hi, int B_Frame)

//  Return the allowable start moves for the overlap between
//  gene A at positions  A_Lo .. A_Hi  in frame  A_Frame
//  and gene B at positions  B_Lo .. B_Hi  in frame  B_Frame .
//  Only the sign of the frame matters.  The actual frame
//  maybe wrong for "wraparound" overlaps.


  {
   if  (A_Lo < B_Lo)                      // A is on the left
       {
//        assert (
if  (! (B_Lo <= A_Hi && A_Hi < B_Hi))
    {
     printf ("ERROR:  Unexpected overlap\n");
     printf ("A_Lo = %7ld  A_Hi = %7ld  A_Frame = %2d\n",
                A_Lo, A_Hi, A_Frame);
     printf ("B_Lo = %7ld  B_Hi = %7ld  B_Frame = %2d\n",
                B_Lo, B_Hi, B_Frame);
     exit (-2);
    }
        if  (A_Frame > 0)
            {
             if  (B_Frame > 0)
                 return  ONLY_B_CAN_MOVE;
               else
                 return  NEITHER_CAN_MOVE;
            }
          else
            {
             if  (B_Frame > 0)
                 return  BOTH_CAN_MOVE;
               else
                 return  ONLY_A_CAN_MOVE;
            }
       }
     else                                // A is on the right
       {
//        assert (A_Lo <= B_Hi && B_Hi < A_Hi);
if  (! (A_Lo <= B_Hi && B_Hi < A_Hi))
    {
     printf ("ERROR:  Unexpected overlap\n");
     printf ("A_Lo = %7ld  A_Hi = %7ld  A_Frame = %2d\n",
                A_Lo, A_Hi, A_Frame);
     printf ("B_Lo = %7ld  B_Hi = %7ld  B_Frame = %2d\n",
                B_Lo, B_Hi, B_Frame);
     exit (-2);
    }
        if  (A_Frame > 0)
            {
             if  (B_Frame > 0)
                 return  ONLY_A_CAN_MOVE;
               else
                 return  BOTH_CAN_MOVE;
            }
          else
            {
             if  (B_Frame > 0)
                 return  NEITHER_CAN_MOVE;
               else
                 return  ONLY_B_CAN_MOVE;
            }
       }
  }



double  Loop_Cost  (int M, int N)

/* Return the energy cost of a loop with  M  bases on one side and
*  N  bases on the other. */

  {
   double  Cost;
   int  Min;

   if  (M <= 0 || N <= 0)
       return  BIG_NEGATIVE;

   if  (M < N)
       Min = M;
     else
       Min = N;

   if  (Min < 4)
       Cost = -0.8;
     else
       Cost = -0.8 - (Min - 3) * (8.4 - 0.8) / 27.0;

   if  (M == N)
       return  Cost;

   return  Cost - abs (M - N) * 1.1;
  }



void  Make_Codon_Log_Prob
    (double clp [64], double bp [4])

//  Set entries in  clp  to the log of the probability of
//  each codon using the indepedent base probabilities in  bp
//  but normalized so that stop codons have zero probability.

  {
   char  alphabet [] = "acgt";
   char  codon [4] = "aaa";
   double  sum = 0.0;
   int  i, j, k, sub;

   for  (i = 0;  i < 4;  i ++)
     {
      codon [0] = alphabet [i];
      for  (j = 0;  j < 4;  j ++)
        {
         codon [1] = alphabet [j];
         for (k = 0;  k < 4;  k ++)
           {
            codon [2] = alphabet [k];

            sub = 16 * i + 4 * j + k;
            if  (Is_Stop (codon))
                clp [sub] = 0.0;
              else
                {
                 clp [sub] = bp [i] * bp [j] * bp [k];
                 sum += clp [sub];
                }
           }
        }
     }

   if  (sum <= 0.0 || sum > 1.0)
       {
        fprintf (stderr, "ERROR:  Bad sum = %f in  Make_Codon_Log_Prob\n",
                 sum);
        exit (EXIT_FAILURE);
       }

   for  (i = 0;  i < 64;  i ++)
     if  (clp [i] == 0.0)
         clp [i] = -1000.0;
       else
         clp [i] = log (clp [i] / sum);

   return;
  }



int  Make_Final_Changes
    (Gene_Ref * Ptr [], Gene_Ref Gene [], long int N)

//  Make all remaining changes to entries in  Gene [1 .. N]  that
//  are caused by other entries.  Make changes in order of pointers
//  in  Ptr [1 .. N] .  Return the number of changes made.

  {
   Overlap_Node  * P;
   long int  i, Change_Ct, Delay_Cause, Delay_Len;

               //  First make all rejects

   Change_Ct = 0;
   for  (i = 1;  i <= N;  i ++)
     if  (Ptr [i] -> Has_Status (MIGHT_CHANGE))
         {
          for  (P = Ptr [i] -> Overlap_List;  P != NULL;  P = P -> Next)
            if  (P -> Problem_Code == REJECT_CHAR
                   && ! Gene [P -> From] . Has_Status (REJECTED))
                {
                 Ptr [i] -> Clear_Status (MIGHT_CHANGE);
                 Ptr [i] -> Set_Status (REJECTED);
                 Change_Ct ++;
                 break;
                }
         }

               //  Make longest delay for each remaining gene

   for  (i = 1;  i <= N;  i ++)
     if  (Ptr [i] -> Has_Status (MIGHT_CHANGE))
         {
          Delay_Len = 0;
          Delay_Cause = 0;
          for  (P = Ptr [i] -> Overlap_List;  P != NULL;  P = P -> Next)
            if  (! Gene [P -> From] . Has_Status (REJECTED))
                switch  (P -> Problem_Code)
                  {
                   case  REJECT_CHAR :
                     fprintf (stderr, "ERROR:  Unexpected reject\n");
                     assert (FALSE);
                     exit (-1);
                   case  SCORES_WORSE_CHAR :
                   case  ELIM_OLAP_CHAR :
                   case  SHORTER_CHAR :
                     if  (P -> Delay == 0)
                         break;
                     if  (P -> Delay > Delay_Len)
                         {
                          Delay_Len = P -> Delay;
                          Delay_Cause = P -> From;
                         }
                     break;
                   case  CONTAINS_CHAR :
                   case  SHADOWED_CHAR :
                   case  NOPROB_CHAR :
                     break;
                  }
          if  (Delay_Len > 0)
              {
               Change_Ct ++;
               Ptr [i] -> Delay_Len += Delay_Len;
               Ptr [i] -> Delay_Cause = Delay_Cause;
               if  (Ptr [i] -> Frame > 0)
                   Ptr [i] -> Lo += Delay_Len;
                 else
                   Ptr [i] -> Hi -= Delay_Len;
               Ptr [i] -> Len -= Delay_Len;
              }
         }

   return  Change_Ct;
  }



int  Make_Sure_Changes
    (Gene_Ref Gene [], long int N)

//  Make any changes to entries in  Gene [1 .. N]  that are caused
//  by other entries that are *not* going to change.  Return the
//  number of changes made.

  {
   Overlap_Node  * P;
   int  Still_Might_Change;
   long int  i, Change_Ct, Delay_Cause, Delay_Len;

   Change_Ct = 0;
   for  (i = 1;  i <= N;  i ++)
     if  (Gene [i] . Has_Status (MIGHT_CHANGE))
         {
          Still_Might_Change = FALSE;
          Delay_Len = 0;
          Delay_Cause = 0;
          for  (P = Gene [i] . Overlap_List;  P != NULL;  P = P -> Next)
            switch  (P -> Problem_Code)
              {
               case  REJECT_CHAR :
                 if  (Gene [P -> From] . Has_Status (MIGHT_CHANGE))
                     Still_Might_Change = TRUE;
                 else if  (! Gene [P -> From] . Has_Status (REJECTED))
                     Gene [i] . Set_Status (REJECTED);
                 break;
               case  SCORES_WORSE_CHAR :
               case  ELIM_OLAP_CHAR :
               case  SHORTER_CHAR :
                 if  (P -> Delay == 0)
                     break;
                 if  (Gene [P -> From] . Has_Status (MIGHT_CHANGE))
                     Still_Might_Change = TRUE;
                 else if  (! Gene [P -> From] . Has_Status (REJECTED))
                     {
                      if  (P -> Delay > Delay_Len)
                          {
                           Delay_Len = P -> Delay;
                           Delay_Cause = P -> From;
                          }
                     }
                 break;
               case  CONTAINS_CHAR :
               case  SHADOWED_CHAR :
               case  NOPROB_CHAR :
                 break;
              }
          if  (Gene [i] . Has_Status (REJECTED))
              {
               Gene [i] . Clear_Status (MIGHT_CHANGE);
               Change_Ct ++;
              }
          else if  (Delay_Len > 0)
              {
               Change_Ct ++;
               Gene [i] . Delay_Len += Delay_Len;
               Gene [i] . Delay_Cause = Delay_Cause;
               if  (Gene [i] . Frame > 0)
                   Gene [i] . Lo += Delay_Len;
                 else
                   Gene [i] . Hi -= Delay_Len;
               Gene [i] . Len -= Delay_Len;
               Still_Might_Change = FALSE;
               for  (P = Gene [i] . Overlap_List;  P != NULL;  P = P -> Next)
                 if  (P -> Delay > Delay_Len
                        && ! Gene [P -> From] . Has_Status (REJECTED))
                     {
                      Still_Might_Change = TRUE;
                      break;
                     }
               if  (! Still_Might_Change)
                   Gene [i] . Clear_Status (MIGHT_CHANGE);
              }
         }

   return  Change_Ct;
  }



int  Match  (char P, char Q)

/* Return true iff bases  P  and  Q  bind together. */

  {
   Q = tolower (Q);
   switch  (tolower (P))
     {
      case  'a' :
        return  (Q == 't');
      case  'c' :
        return  (Q == 'g');
      case  'g' :
        return  (Q == 'c' || Q == 't');
      case  't' :
        return  (Q == 'a' || Q == 'g');
     }

   return  0;
  }



void  Permute  (int S [], int F)

/* Permute the values in  S  which are in frame  F  to make them
*  consistent with scores from frame  +1 . */

  {
   int  Save;

   switch  (F)
     {
      case  1 :
        Save = S [4];
        S [4] = S [5];
        S [5] = Save;
        return;
      case  2 :
        Save = S [0];
        S [0] = S [2];
        S [2] = S [1];
        S [1] = Save;
        Save = S [3];
        S [3] = S [5];
        S [5] = Save;
        return;
      case  3 :
        Save = S [0];
        S [0] = S [1];
        S [1] = S [2];
        S [2] = Save;
        Save = S [3];
        S [3] = S [4];
        S [4] = Save;
        return;
      case  -1 :
        Save = S [0];
        S [0] = S [3];
        S [3] = Save;
        Save = S [1];
        S [1] = S [5];
        S [5] = S [2];
        S [2] = S [4];
        S [4] = Save;
        return;
      case  -2 :
        Save = S [0];
        S [0] = S [4];
        S [4] = S [2];
        S [2] = S [5];
        S [5] = Save;
        Save = S [1];
        S [1] = S [3];
        S [3] = Save;
        return;
      case  -3 :
        Save = S [0];
        S [0] = S [5];
        S [5] = S [1];
        S [1] = S [4];
        S [4] = Save;
        Save = S [2];
        S [2] = S [3];
        S [3] = Save;
        return;
     }
   return;
  }


static void  Print_Category
    (Gene_Category_Type Category)

//  Print to  stdout  a representation of  Category .

  {
   switch  (Category)
     {
      case  NONE :
      case  REGULAR :
        printf ("     ");
        break;
      case  VOTED :
        printf (" vote");
        break;
      case  WEAK :
        printf (" weak");
        break;
      default :
        fprintf (stderr, "ERROR:  Unexpected category = %d\n",
                 (int) Category);
        exit (-1);
     }

   return;
  }


static void  Print_Separate_Score_Headings
    (void)

// Print headings for scoring of separate orfs without overlap/voting
// rules.

  {
   printf ("Threshold score = %d\n", Threshold_Score);
   printf ("Ignore independent score on orfs longer than %ld\n",
           Ignore_Indep_Len - 1);
   putchar ('\n');
   printf ("                             ---- Frame Scores ----"
           "  Indep    Raw\n");
   printf ("   Tag      Start      End   +1  +2  +3  -1  -2  -3"
           "  Score   Score\n");

   return;
  }



void  Process_Ignore ( Ignore_Ptr &Ignore_Array )

// Function to take a file of ignore regions to an array of pointers 
// of type Ignore_Node. The input filename obtained by Process_Options
// is used, and the address pairs are sorted in assending order
// afterwhich all overlapping entries are combined.

   {
     FILE * Ignore_File;
     char Ignore_Line [MAX_INPUT], *Token, Delim[] = "\t\n\r\f\x20";
     int  i, j, Temp[2], Current, Line_No, Mem_Hunk_Incr;
     int  Ignore_Regions = 0;
     int  Shrink_Ignore;     

     if ( !Ignore_Option )
      {
	Ignore_Array = (Ignore_Ptr) Safe_malloc (sizeof(Ignore_Node));
	Ignore_Array[0].Low_Address  = -1;
	Ignore_Array[0].High_Address = -1;
	return;
      }

     Mem_Hunk_Incr = 10;
     assert (Ignore_File_Name != NULL);
     Ignore_File = File_Open (Ignore_File_Name, "r");
     Ignore_Array = (Ignore_Ptr) Safe_malloc (Mem_Hunk_Incr*sizeof(Ignore_Node));
     
     i = 0;
     Line_No = 0;
     while ( fgets (Ignore_Line, MAX_INPUT, Ignore_File) )
       {
	 Line_No++;
	 if ( !(i % Mem_Hunk_Incr ) && i > 0)
	   Ignore_Array = (Ignore_Ptr) Safe_realloc  
	     (Ignore_Array, (i + Mem_Hunk_Incr + 1)*sizeof(Ignore_Node));
	 
	 Token  = strtok (Ignore_Line, Delim);	
	 j = -1;
	 while ( Token )
	   {
	     if( isdigit(Token[0]) )
	       Temp[++j] = atoi( Token );
	     
	     if( j == 1 )
	       break;
	     
	     Token = strtok (NULL, Delim);
	   }
	 
	 switch ( j )
	   {
	   case -1 :   // no integers on the line get another line
	     break;
	     
	   case  0 :   // one iteger on the line this is a probably an error
	     fprintf(stderr, "ERROR: Error on line %d of ignore file ", Line_No );
	     fprintf(stderr, "it must contain either two integers or none.\n");
	     fprintf(stderr, "       Line %d skipped.\n", Line_No );
	     break;
	     
	   case 1  :
	     if ( Temp[0] < Temp[1] )
	       {
		 Shrink_Ignore = (Temp[1] - Temp[0] + 1) % 3;
		 if ( Temp[1] - Temp[0]  >  3 )
		    Temp[1] -= Shrink_Ignore;
		 Ignore_Array[i].Low_Address  = Temp[0];
		 Ignore_Array[i].High_Address = Temp[1];
	       }
	     else
	       {
		 Shrink_Ignore = (Temp[0] - Temp[1] + 1) % 3;
		 if ( Temp[0] - Temp[1]  >  3 )
		    Temp[0] -= Shrink_Ignore;
		 Ignore_Array[i].Low_Address  = Temp[1];
		 Ignore_Array[i].High_Address = Temp[0];
	       }

	     Ignore_Regions++;
	     i++;
	     break;
	     
	   default :
	     break;
	   }
       }
     fclose (Ignore_File);
     
     qsort(Ignore_Array, Ignore_Regions, sizeof(Ignore_Node), Cmp_Ignore);
     
     // reduce the overlapping ignore regions if any.
     Current = 0;
     for ( i = 1; i < Ignore_Regions; i++ )
       {
	 if( Ignore_Array[Current].High_Address > Ignore_Array[i].High_Address )
	   {
	     // do nothing, i is completely shadowed
	   }
	 else if ( Ignore_Array[Current].High_Address > Ignore_Array[i].Low_Address )
	   {      //  arrays overlap combine both
	     Ignore_Array[Current].High_Address = Ignore_Array[i].High_Address;
	   }
	 else if ( ++Current < i )
	   {     
	     Ignore_Array[Current].Low_Address  = Ignore_Array[i].Low_Address;
	     Ignore_Array[Current].High_Address = Ignore_Array[i].High_Address;
	   }
       }
     Ignore_Regions = Current + 1;
     Ignore_Array[Ignore_Regions].Low_Address  = -1;   // Signal end of data
     Ignore_Array[Ignore_Regions].High_Address = -1;   //  

   }


void  Process_Options  (int argc, char * argv [])

//  Process command-line options and set corresponding global switches
//  and parameters.
//
//    -C n   Use n as GC percentage for independent model
//    -f     Use ribosome-binding energy to choose start codon
//    +f     Use first codon in orf as start codon
//    -g n   Set minimum gene length to n
//    -i <filename> Use <filename> to select regions of bases that are off 
//           limits, so that no bases within that area will be examined
//    -l     Assume linear rather than circular genome, i.e., no wraparound
//    -L <filename> Use <filename> to specify a list of orfs that should
//           be scored separately, with no overlap rules
//    -M     Input is a multifasta file of separate genes to be scored
//           separately, with no overlap rules
//    -o n   Set minimum overlap length to n.  Overlaps shorter than this
//           are ignored.
//    -p n   Set minimum overlap percentage to n%.  Overlaps shorter than
//           this percentage of *both* strings are ignored.
//    -q n   Set the maximum length orf that can be rejected because of
//           the independent probability score column to n
//    -r     Don't use independent probability score column
//    +r     Use independent probability score column
//    -s s   Use string s as the ribosome binding pattern to find start codons.
//    -S     Don't use stricter independent intergenic model
//    +S     Do use stricter independent intergenic model that doesn't
//           give probabilities to in-frame stop codons.
//    -t n   Set threshold score for calling as gene to n.  If the in-frame
//           score >= n, then the region is given a number and considered
//           a potential gene.
//    -w n   Use "weak" scores on tentative genes n or longer.  Weak
//           scores ignore the independent probability score.
//    -X     Allow orfs that extend off the ends of the sequence

  {
   char  * P;
   long int  W;
   double  X;
   bool  errflg = false;
   int  i, L;

   for  (i = 3;  i < argc;  i ++)
     {
      switch  (argv [i] [0])
        {
         case  '-' :
           switch  (argv [i] [1])
             {
              case  'C' :       // use value as GC percentage of independent model
                GC_From_Parameter = true;
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                X = strtod (P, NULL);
                if  (errno == ERANGE || X < 0.0 || X > 100.0)
                    fprintf (stderr, "ERROR:  Bad GC percentage %s\n", P);
                  else
                    {
                     Ch_Ct [1] = Ch_Ct [2] = X / 200.0;
                     Ch_Ct [0] = Ch_Ct [3] = 0.50 - Ch_Ct [1];
                    }
                break;
              case  'f' :       // use function to choose start codon in gene
                Choose_First_Start_Codon = FALSE;
                break;
              case  'g' :       // minimum gene length
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                W = strtol (P, NULL, 10);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad minimum gene length %s\n", P);
                  else
                    Min_Gene_Len = W;
                assert (Min_Gene_Len > 0);
                break;
	      case  'i' :       // ignore regions of the genome
		if  (argv [i] [2] != '\0')
		  P = argv [i] + 2;
		else
		  P = argv [++ i];
                P = strtok (P, " \t\n\"'");
                L = strlen (P);
		Ignore_File_Name = (char *) Safe_malloc(MAX_LINE);
		strcpy( Ignore_File_Name, P );
		Ignore_Option = TRUE;
                break;
              case  'l' :       // linear, not circular genome
                Genome_Is_Circular = FALSE;
                break;
	      case  'L' :       // list of separate orfs to score
		if  (argv [i] [2] != '\0')
		  P = argv [i] + 2;
		else
		  P = argv [++ i];
                P = strtok (P, " \t\n\"'");
                L = strlen (P);
		Orflist_File_Name = (char *) Safe_malloc(MAX_LINE);
		strcpy( Orflist_File_Name, P );
		Orflist_Option = TRUE;
                break;
              case  'M' :       // input is multifasta list of genes to score separately
                Separate_Multifasta_Orfs = true;
                break;
              case  'o' :       // minimum overlap length
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                W = strtol (P, NULL, 10);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad minimum overlap length %s\n", P);
                  else
                    Min_Olap = W;
                assert (Min_Olap >= 0);
                if  (Min_Olap == 0)
                    Min_Olap = 1;
                break;
              case  'p' :       // minimum overlap percent
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                X = strtod (P, NULL);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad minimum overlap percent %s\n", P);
                  else
                    Min_Olap_Percent = X / 100.0;
                assert (Min_Olap_Percent >= 0.0 && Min_Olap_Percent <= 100.0);
                break;
              case  'q' :       // max length of orf rejected by independent score
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                W = strtol (P, NULL, 10);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad reject orf length %s\n", P);
                  else
                    Ignore_Indep_Len = W;
                assert (Ignore_Indep_Len >= 0 && Ignore_Indep_Len < LONG_MAX);
                break;
              case  'r' :       // don't use random/independent score column
                Use_Independent = FALSE;
                break;
              case  's' :       // string for ribosome binding pattern
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                P = strtok (P, " \t\n\"'");
                L = strlen (P);
                assert (L <= MAX_RIBOSOME_PATTERN_LEN);
                for  (i = 0;  i < L;  i ++)
                  Ribosome_Pattern [i] = Filter (P [i]);
                Ribosome_Pattern [L] = '\0';
                break;
#if  0
              case  'S' :       // don't use strict independent model
                Use_Strict_Independent = FALSE;
                break;
#endif
              case  't' :       // threshold score for calling as gene
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                W = strtol (P, NULL, 10);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad threshold score %s\n", P);
                  else
                    Threshold_Score = W;
                assert (Threshold_Score > 0 && Threshold_Score < 100);
                break;
              case  'w' :       // length of genes to use weak score on
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                W = strtol (P, NULL, 10);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad weak score length %s\n", P);
                  else
                    Min_Weak_Len = W;
                assert (Min_Weak_Len >= 0 && Min_Weak_Len < INT_MAX);
                break;
              case  'X' :       // process orfs that extend off sequence ends
                Allow_Partial_Orfs = TRUE;
                break;
              default :
                fprintf (stderr, "Unrecognized option %s\n", argv [i]);
                errflg = true;
             }
           break;
         case  '+' :
           switch  (argv [i] [1])
             {
              case  'f' :       // automatically use first start codon in gene
                Choose_First_Start_Codon = TRUE;
                break;
              case  'r' :       // use random/independent score column
                Use_Independent = TRUE;
                break;
              case  'S' :       // use strict version of independent model
                Use_Strict_Independent = TRUE;
                break;
              default :
                fprintf (stderr, "Unrecognized option %s\n", argv [i]);
                errflg = true;
             }
           break;
         default :
           fprintf (stderr, "Unrecognized option %s\n", argv [i]);
        }
     }

   if  (errflg)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   return;
  }



void  Process_Orflist
    (void)

//  Read the list of orf coordinates from file name  Orflist_File_Name
//  and score each separately, printing the results to stdout.
//  No overlap rules are used and the orfs are not added to the
//  list for any further processing.

  {
   FILE  * fp;
   char  line [MAX_LINE];
   char  tag [MAX_LINE];
   int  score [7];
   int  weak_score;
   double  raw_score;
   int  total, bad_ct;
   int  i, start, stop;

   fp = File_Open (Orflist_File_Name, "r");

   total = bad_ct = 0;

   while  (fgets (line, MAX_LINE, fp) != NULL)
     {
      if  (sscanf (line, "%s %d %d", tag, & start, & stop) != 3)
          continue;

      total ++;

      Score_String (start, stop, Ch_Ct, score, TRUE,
                    weak_score, raw_score);
      printf ("%8s %8d %8d ", tag, start, stop);
      for  (i = 0;  i < 6;  i ++)
        printf (" %3d", score [i]);
      printf (" %6d", score [6]);
      printf (" %7.3f", raw_score);

      if  (score [0] < Threshold_Score)
          {
           if  (1 + abs (stop - start) >= Ignore_Indep_Len)
               printf ("  ignore--too long");
             else
               {
                bad_ct ++;
                printf ("  reject");
               }
          }

      putchar ('\n');
     }

   fclose (fp);

   printf ("\nReject = %d (%.1f%%)  Keep = %d (%.1f%%)  Total = %d\n",
           bad_ct, Percent (bad_ct, total), total - bad_ct,
           Percent (total - bad_ct, total), total);

   return;
  }



void  Read_Probability_Model  (char * Param)

//  Read in the probability model indicated by  Param .

  {
   FILE  * fp;

   fp = File_Open (Param, "r");   // maybe rb ?

   Read_Scoring_Model (fp);

   fclose (fp);

   return;
  }



void  Score_Multifasta_Orfs
    (FILE * fp)

//  Read each fasta record in  fp  and score it as a separate
//  gene.  Print its scores to  stdout  and the string
//  "reject" at the end if its main score is below global
//   Threshold_Score .

  {
   int  total, bad_ct;

   total = bad_ct = 0;

   while  (Read_String (fp, Data, Input_Size, Name, FALSE))
     {
      int  score [7], has_stop [7];
      int  weak_score;
      double  raw_score;
      bool  all_stops;
      int  i, len;

      len = strlen (Data + 1);
      total ++;

      printf ("%8s %8d %8d ", Name, 1, len);
      if  (len % 3 != 0)
          {
           printf (" len = %d mod 3  frameshift?\n", len % 3);
           continue;
          }

      Find_Stop_Codons (Data, len, has_stop);
      all_stops = true;
      for  (i = 0;  i < 6 && all_stops;  i ++)
        all_stops = has_stop [i];
      if  (all_stops)
          {
           printf (" stop in all 6 reading frames   reject\n");
           continue;
          }

      Score_String (1, len, Ch_Ct, score, TRUE,
                    weak_score, raw_score);
      for  (i = 0;  i < 6;  i ++)
        printf (" %3d", score [i]);
      printf (" %6d", score [6]);
      printf (" %7.3f", raw_score);

      if  (score [0] < Threshold_Score)
          {
           if  (len >= Ignore_Indep_Len)
               printf ("  ignore--too long");
             else
               {
                bad_ct ++;
                printf ("  reject");
               }
          }

      putchar ('\n');
     }

   fclose (fp);

   printf ("\nReject = %d (%.1f%%)  Keep = %d (%.1f%%)  Total = %d\n",
           bad_ct, Percent (bad_ct, total), total - bad_ct,
           Percent (total - bad_ct, total), total);

   return;
  }



void  Score_Olap_Region
    (long int Lo, long int Hi, int A_Frame, int B_Frame,
     int & A_Score, int & B_Score)

//  Score the region in  Data [Lo .. Hi]  and set  A_Score
//  and  B_Score  to the scores in frames  A_Frame  and  B_Frame
//  respectively.

  {
   int  Frame, Score [7], Weak_Score;
   double  Raw_Score;

   Score_String (Lo, Hi, Ch_Ct, Score, FALSE, Weak_Score, Raw_Score);

   Frame = Lo % 3;
   if  (Frame == 0)
       Frame = 3;
   Permute (Score, Frame);

   A_Score = Choose_Score (Score, A_Frame);
   B_Score = Choose_Score (Score, B_Frame);

   return;
  }



void  Score_String
    (long int Start, long int Stop, double Ch_Ct [ALPHABET_SIZE],
     int Score [7], int Use_Independent, int & Weak_Score,
     double & Raw_Score)

//  Set  Score []  to integer scores from 0 to 99 for the string in
//  global  Data [Start .. Stop]  using the simple Markov model in
//  global  Context_Delta  and the simple independent log probabilities
//  in  Ch_Ct .  Data_Len  is the last position in  Data  to compute
//  wraparounds.
//  Scores in  Score [0 .. 5]  represent scores in each of 6 reading
//  frames.  If  Use_Independent  is true, also use the independent model
//  and put its score in  Score [6] ; otherwise, set  Score [6]  to an
//  artificially low value.  Set  Weak_Score  to what  Score [0]
//  would be if the independent model were not used.  Set  Raw_Score
//  to the log probability per base of  Score [0]  without normalizing
//  to 0..99 with the other scores.

  {
   long int  Len;

   if  (Start < Stop)
       {
        Len = 1 + Stop - Start;
        if  (2 + Len > Orf_Buffer_Len)
            {
             Orf_Buffer_Len = Max (2 + Len, Orf_Buffer_Len + ORF_SIZE_INCR);
             Orf_Buffer = (char *) Safe_realloc (Orf_Buffer, Orf_Buffer_Len);
            }
        Transfer (Orf_Buffer + 1, Start, Len);
       }
     else
       {
        Len = 1 + Start - Stop;
        if  (2 + Len > Orf_Buffer_Len)
            {
             Orf_Buffer_Len = Max (2 + Len, Orf_Buffer_Len + ORF_SIZE_INCR);
             Orf_Buffer = (char *) Safe_realloc (Orf_Buffer, Orf_Buffer_Len);
            }
        Transfer (Orf_Buffer + 1, Start, - Len);
       }

   Simple_Score (Orf_Buffer, Len, MODEL_LEN, Ch_Ct,
                   Score, Use_Independent, Weak_Score, Raw_Score);

   return;
  }



static void  Set_Ignore_Indep_Len
    (void)

//  Print GC content from global  Ch_Ct  and then set global
//   Ignore_Indep_Len  to length of orf would expect to occur
//  once at random in a million bases.

  {
   double  poisson_lambda;

   printf ("GC Proportion = %.1f%%\n", 100.0 * (Ch_Ct [1] + Ch_Ct [2]));

   if  (Ignore_Indep_Len == LONG_MAX)
       {
        poisson_lambda = Ch_Ct [3] * Ch_Ct [0]
                           * (2.0 * Ch_Ct [2] + Ch_Ct [0]);
        assert (poisson_lambda != 0.0);
        Ignore_Indep_Len
             = (long int) floor (3.0 * log (2.0 * 1000000 * poisson_lambda)
                 / poisson_lambda);
       }

   return;
  }



void  Set_Indep_Probs_From_Data
    (double Ch_Ct [], FILE * fp)

//  Set entries in  Ch_Ct  to probabilities of each letter in
//  multifasta file  fp  (which must already be opened).
//  Rewind  fp  when finished.

  {
   long int  total = 0L;
   int  i;

   for  (i = 0;  i < ALPHABET_SIZE;  i ++)
     Ch_Ct [i] = 0.0;

   while  (Read_String (fp, Data, Input_Size, Name, FALSE))
     {
      int  len;

      len = strlen (Data + 1);

      for  (i = 1;  i <= len;  i ++)
        {
         switch (tolower (Data [i]))
           {
            case  'a' :
            case  't' :
              Ch_Ct [0] += 1.0;
              total ++;
              break;
            case  'c' :
            case  'g' :
              Ch_Ct [1] += 1.0;
              total ++;
              break;
           }
        }
     }

   Ch_Ct [2] = Ch_Ct [1];
   Ch_Ct [3] = Ch_Ct [0];
   for  (i = 0;  i < ALPHABET_SIZE;  i ++)
     Ch_Ct [i] = Ch_Ct [i] / (2.0 * total);

   rewind (fp);

   return;
  }



void  Show_Gene_Info
    (Gene_Ref Gene [], long int N)

//  Print out summary information on what's in  Gene [1 .. N] .

  {
   Overlap_Node  * P;
   int  Is_Maybe_Reject, Is_Sure_Reject;
   long int  Change_Ct = 0L, Reject_Ct = 0L, Sure_Reject_Ct = 0L;
   long int  Olap_Ct = 0L, Potential_Ct = 0L;
   long int  i;

   for  (i = 1;  i <= N;  i ++)
     {
      if  (Global_Show_Details)
          printf ("Gene %5ld  %+2d  %7ld  %7ld  %4ld  %03o  Sco = %2d\n",
                  i, Gene [i] . Frame, Gene [i] . Lo, Gene [i] . Hi,
                  Gene [i] . Len, Gene [i] . Status, Gene [i] . Score);
      if  (Gene [i] . Has_Status (REJECTED))
          continue;
      Potential_Ct ++;
      if  (Gene [i] . Has_Status (MIGHT_CHANGE))
          Change_Ct ++;
      Is_Maybe_Reject = Is_Sure_Reject = FALSE;
      for  (P = Gene [i] . Overlap_List;  P != NULL;  P = P -> Next)
        {
         if  (Gene [P -> From] . Has_Status (REJECTED))
             continue;
         Olap_Ct ++;
         if  (P -> Problem_Code == REJECT_CHAR)
             {
              Is_Maybe_Reject = TRUE;
              if  (! Gene [P -> From] . Has_Status (MIGHT_CHANGE))
                  {
                   Is_Sure_Reject = TRUE;
                   if  (Global_Show_Details)
                       printf ("    Gene %5ld (%03o) is sure reject from %5ld (%03o)\n",
                               i, Gene [i] . Status, P -> From,
                               Gene [P -> From] . Status);
                  }
             }
         if  (Global_Show_Details)
             printf ("    From %5ld  %7ld  %7ld  %4ld  Del = %3ld  Sco = %2d  Cod = `%c'  OF = %+2d\n",
                     P -> From, P -> Lo, P -> Hi, P -> Olap, P -> Delay,
                     P -> Score, P -> Problem_Code, P -> Other_Frame);
        }
      if  (Is_Maybe_Reject)
          Reject_Ct ++;
      if  (Is_Sure_Reject)
          Sure_Reject_Ct ++;
     }

   putchar ('\n');
   printf ("   Original Genes = %5ld\n", N);
   printf ("  Potential Genes = %5ld\n", Potential_Ct);
   printf ("        Avg Olaps = %5.1f\n",
           Potential_Ct == 0 ? 0.0 : double (Olap_Ct) / Potential_Ct);
   printf ("Potential Changes = %5ld\n", Change_Ct);
   printf ("Potential Rejects = %5ld\n", Reject_Ct);
   printf ("     Sure Rejects = %5ld\n", Sure_Reject_Ct);

   return;
  }



void  Simple_Score
    (char X [], long int T, int Model_Len, double Ch_Ct [ALPHABET_SIZE],
     int Score [], int Use_Independent, int & Weak_Score,
     double & Raw_Score)

//  Set  Score  to the probabilites of string  X [1 .. T]  being
//  generated in each of the 3 forward and 3 reverse reading frames
//  using simple nonhomogeneous Markov models in global  Context_Delta []
//  with model length equal to  Model_Len .   If  Use_Independent  is true,
//  also use the independent model and put its score in  Score [6] ;
//  otherwise, set  Score [6]  to an artificially low value.
//  Set  Weak_Score  to what  Score [0]  would be in the independent
//  model is not used.  Set  Raw_Score  to the log probability per base
//  of  Score [0]  without normalizing to 0..99 with the other scores.

  {
   double  Max, Min, Sum, S [7], W [7];
   double  Weak_Max, Weak_Min, Weak_Sum;
   int  i, Has_Stop [7], Offset;

   Find_Stop_Codons (X, T, Has_Stop);

   Max = - DBL_MAX;
   Min = DBL_MAX;

   if  (! Has_Stop [0])
       {
        Fast_Evaluate (X + 1, T, Model_Len, MODEL [0],
                         MODEL [1], MODEL [2], S [0]);
        if  (S [0] > Max)
            Max = S [0];
        if  (S [0] < Min)
            Min = S [0];
       }
   if  (! Has_Stop [1])
       {
        Fast_Evaluate (X + 1, T, Model_Len, MODEL [2],
                         MODEL [0], MODEL [1], S [1]);
        if  (S [1] > Max)
            Max = S [1];
        if  (S [1] < Min)
            Min = S [1];
       }
   if  (! Has_Stop [2])
       {
        Fast_Evaluate (X + 1, T, Model_Len, MODEL [1],
                         MODEL [2], MODEL [0], S [2]);
        if  (S [2] > Max)
            Max = S [2];
        if  (S [2] < Min)
            Min = S [2];
       }

   Offset = T % 3;
   Reverse_Complement (X, T);
   if  (! Has_Stop [3 + Offset])
       {
        Fast_Evaluate (X + 1, T, Model_Len, MODEL [0],
                         MODEL [1], MODEL [2], S [3 + Offset]);
        if  (S [3 + Offset] > Max)
            Max = S [3 + Offset];
        if  (S [3 + Offset] < Min)
            Min = S [3 + Offset];
       }
   Offset = (Offset + 1) % 3;
   if  (! Has_Stop [3 + Offset])
       {
        Fast_Evaluate (X + 1, T, Model_Len, MODEL [1],
                         MODEL [2], MODEL [0], S [3 + Offset]);
        if  (S [3 + Offset] > Max)
            Max = S [3 + Offset];
        if  (S [3 + Offset] < Min)
            Min = S [3 + Offset];
       }
   Offset = (Offset + 1) % 3;
   if  (! Has_Stop [3 + Offset])
       {
        Fast_Evaluate (X + 1, T, Model_Len, MODEL [2],
                         MODEL [0], MODEL [1], S [3 + Offset]);
        if  (S [3 + Offset] > Max)
            Max = S [3 + Offset];
        if  (S [3 + Offset] < Min)
            Min = S [3 + Offset];
       }

   Weak_Max = Max;
   Weak_Min = Min;

   Has_Stop [6] = ! Use_Independent;
   if  (Use_Independent)
       {
        Indep_Eval (X, T, Ch_Ct, S [6]);
        if  (S [6] > Max)
            Max = S [6];
        if  (S [6] < Min)
            Min = S [6];
       }

   assert (Max != - DBL_MAX && Min != DBL_MAX);

   if  (Min < Max + MAX_LOG_DIFF)
       Min = Max + MAX_LOG_DIFF;
   if  (Weak_Min < Weak_Max + MAX_LOG_DIFF)
       Weak_Min = Weak_Max + MAX_LOG_DIFF;

   Sum = 0.0;
   for  (i = 0;  i < 7;  i ++)
     if  (Has_Stop [i])
         W [i] = -1.0;
     else if  (S [i] >= Min)
         {
          W [i] = exp (S [i] - Min);
          Sum += W [i];
         }
       else
         W [i] = 0.0;
   assert (Sum > 0.0);

   for  (i = 0;  i < 7;  i ++)
     if  (Has_Stop [i])
         Score [i] = -1;
       else
         {
          Score [i] = int (100.0 * (W [i] / Sum));
          if  (Score [i] > 99)
              Score [i] = 99;
         }

   if  (! Use_Independent)
       Weak_Score = Score [0];
     else
       {
        Weak_Sum = 0.0;
        for  (i = 0;  i < 6;  i ++)
          if  (! Has_Stop [i] && S [i] >= Weak_Min)
              {
               W [i] = exp (S [i] - Weak_Min);
               Weak_Sum += W [i];
              }
        assert (Weak_Sum > 0.0);
        Weak_Score = int (100.0 * (W [0] / Weak_Sum));
        if  (Weak_Score > 99)
            Weak_Score = 99;
       }

   Raw_Score = S [0] / T;

   if  (Score [6] < 0)
       Score [6] = 0;

   return;
  }



void  Slide_One_Start
    (Gene_Ref & Gene, Overlap_Node * P, long int & Distance_Moved)

//  Try to move the start of gene  Gene  until it scores better
//  than the other frame in the remaining portion of the
//  overlap pointed to by  P .  It also must meet the minimum
//  length and overall score criteria for being a gene.
//  If successful, set  Distance_Moved  to the number of bases the
//  start shifted; otherwise, set  Distance_Moved  to  0 .

  {
   char  Codon [4];
   int  Olap_Is_OK, Score [7], This_Score, Other_Score, Weak_Score;
   double  Raw_Score;
   long int  i, Len;

   if  (Gene . Frame > 0)
       {
        i = Gene . Lo + P -> Delay;
        Len = 1 + Gene . Hi - i;
        do
          {
           do
             {
              i += 3;
              Len -= 3;
              Transfer (Codon, i, 3);
             }  while  (! Is_Start (Codon) && ! Is_Stop (Codon)
                               && Len > Min_Gene_Len + 2);
           if  (! Is_Start (Codon) || Len < Min_Gene_Len)
               {
                Distance_Moved = 0;
                return;
               }
           if  (1 + P -> Hi - i < MIN_SCORABLE_LEN)
               Olap_Is_OK = TRUE;
             else
               {
                Score_Olap_Region (i, P -> Hi,
                                   Gene . Frame, P -> Other_Frame,
                                   This_Score, Other_Score);
                Olap_Is_OK = (Other_Score <= This_Score);
               }
           Score_String (i, Gene . Hi, Ch_Ct, Score,
                            Use_Independent && Len < Ignore_Indep_Len,
                            Weak_Score, Raw_Score);
          }  while  (Score [0] < Threshold_Score || ! Olap_Is_OK);
        Distance_Moved = i - Gene . Lo - P -> Delay;
       }
     else
       {
        i = Gene . Hi - P -> Delay;
        Len = 1 + i - Gene . Lo;
        do
          {
           do
             {
              i -= 3;
              Len -= 3;
              Transfer (Codon, i, -3);
             }  while  (! Is_Start (Codon) && ! Is_Stop (Codon)
                               && Len > Min_Gene_Len + 2);
           if  (! Is_Start (Codon) || Len < Min_Gene_Len)
               {
                Distance_Moved = 0;
                return;
               }
           if  (1 + i - P -> Lo < MIN_SCORABLE_LEN)
               Olap_Is_OK = TRUE;
             else
               {
                Score_Olap_Region (P -> Lo, i,
                                   Gene . Frame, P -> Other_Frame,
                                   This_Score, Other_Score);
                Olap_Is_OK = (Other_Score <= This_Score);
               }
           Score_String (i, Gene . Lo, Ch_Ct, Score,
                            Use_Independent && Len < Ignore_Indep_Len,
                            Weak_Score, Raw_Score);
          }  while  (Score [0] < Threshold_Score || ! Olap_Is_OK);
        Distance_Moved = Gene . Hi - P -> Delay - i;
       }

   Gene . Set_Status (MIGHT_CHANGE);

   return;
  }



void  Slide_Both_Starts
    (Gene_Ref & Gene_A, Overlap_Node * A_Node, char Failure_Code,
     Gene_Ref & Gene_B, Overlap_Node * B_Node, long int Min_Shift)

//  Try to move the starts of genes  Gene_A  and  Gene_B  to
//  resolve the overlap indicated in their respective overlap
//  nodes  A_Node  and  B_Node .  Try to move  A  first.
//  If unsuccessful, set  A_Node -> Problem_Code  to  Failure_Code .
//  Otherwise, successively alternate between the genes.
//  If ultimately succeed set problem codes to  ELIM_OLAP_CHAR .
//  Otherwise, set the problem code of the gene that fails to
//  move to  SCORES_WORSE_CHAR .  The combination of all moves
//  must be at least  Min_Shift .

  {
   const int  MAX_SLIDE_STEPS = 50;
   long int  i, A_Move, B_Move, Total_Move;

   assert (A_Node -> Delay == 0 && B_Node -> Delay == 0);

   Total_Move = 0;
   for  (i = 0;  i < MAX_SLIDE_STEPS;  i ++)
     {
      Slide_One_Start (Gene_A, A_Node, A_Move);
      if  (A_Move == 0)
          {
           if  (Total_Move == 0)
               A_Node -> Problem_Code = Failure_Code;
             else
               A_Node -> Problem_Code = SCORES_WORSE_CHAR;
           B_Node -> Problem_Code = ELIM_OLAP_CHAR;
           return;
          }
      A_Node -> Delay += A_Move;
      Total_Move += A_Move;

      if  (Total_Move >= Min_Shift)
          {
           A_Node -> Problem_Code = ELIM_OLAP_CHAR;
           B_Node -> Problem_Code = ELIM_OLAP_CHAR;
           return;
          }
      
      Slide_One_Start (Gene_B, B_Node, B_Move);
      if  (B_Move == 0)
          {
           B_Node -> Problem_Code = SCORES_WORSE_CHAR;
           A_Node -> Problem_Code = ELIM_OLAP_CHAR;
           return;
          }
      B_Node -> Delay += B_Move;
      Total_Move += B_Move;

      if  (Total_Move >= Min_Shift)
          {
           A_Node -> Problem_Code = ELIM_OLAP_CHAR;
           B_Node -> Problem_Code = ELIM_OLAP_CHAR;
           return;
          }
     }

   assert (FALSE);

   return;
  }



void  Transfer  (char * S, long int Start, int Len)

//  Transfer  |Len|  characters from  Data [Start ...]  to
//   S  and add null terminator.  If  Len > 0 go in forward direction;
//  otherwise, go in reverse direction and use complements.

  {
   long int  i, j;

   if  (Len > 0)
       {
        for  (i = 0;  i < Len;  i ++)
          {
           j = Start + i;
           if  (j > Data_Len)
               j -= Data_Len;
           else if  (j < 1)
               j += Data_Len;
           S [i] = Filter (tolower (Data [j]));
          }
        S [i] = '\0';
       }
     else
       {
        for  (i = 0;  i < - Len;  i ++)
          {
           j = Start - i;
           if  (j > Data_Len)
               j -= Data_Len;
           else if  (j < 1)
               j += Data_Len;
           S [i] = Filter (tolower (Complement (Data [j])));
          }
        S [i] = '\0';
       }

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
       "USAGE:  %s <genome-file> <icm-file> [options] \n"
       "\n"
       "Find/Score potential genes in <genome-file> using\n"
       "the probability model in <icm-file>\n"
       "\n"
       "Options:\n"
       " -C n   Use n as GC percentage of independent model\n"
       "        Note:  n should be a percentage, e.g., -C 45.2\n"
       " -f     Use ribosome-binding energy to choose start codon\n"
       " +f     Use first codon in orf as start codon\n"
       " -g n   Set minimum gene length to n\n"
       " -i <filename> Use <filename> to select regions of bases that are off \n"
       "        limits, so that no bases within that area will be examined\n"
       " -l     Assume linear rather than circular genome, i.e., no wraparound\n"
       " -L <filename> Use <filename> to specify a list of orfs that should\n"
       "        be scored separately, with no overlap rules\n"
       " -M     Input is a multifasta file of separate genes to be scored\n"
       "        separately, with no overlap rules\n"
       " -o n   Set minimum overlap length to n.  Overlaps shorter than this\n"
       "        are ignored.\n"
       " -p n   Set minimum overlap percentage to n%%.  Overlaps shorter than\n"
       "        this percentage of *both* strings are ignored.\n"
       " -q n   Set the maximum length orf that can be rejected because of\n"
       "        the independent probability score column to (n - 1)\n"
       " -r     Don't use independent probability score column\n"
       " +r     Use independent probability score column\n"
       " -s s   Use string s as the ribosome binding pattern to find start codons.\n"
       " +S     Do use stricter independent intergenic model that doesn't\n"
       "        give probabilities to in-frame stop codons.  (Option is obsolete\n"
       "        since this is now the only behaviour\n"
       " -t n   Set threshold score for calling as gene to n.  If the in-frame\n"
       "        score >= n, then the region is given a number and considered\n"
       "        a potential gene.\n"
       " -w n   Use \"weak\" scores on tentative genes n or longer.  Weak\n"
       "        scores ignore the independent probability score.\n"
       " -X     Allow orfs extending off ends of sequence to be scored\n"
       "\n",
       command);

   return;
  }



