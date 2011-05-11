//  A. L. Delcher
//
//  File:  long-orfs.hh
//
//  Last Modified:  2 Aug 2005
//
//  Declarations for  long-orfs.cc



#ifndef  __LONG_ORFS_HH_INCLUDED
#define  __LONG_ORFS_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "gene.hh"
#include  <queue>


// Default values of global variables

static const bool  DEFAULT_GENOME_IS_CIRCULAR = true;
static const int  DEFAULT_MIN_GENE_LEN = 90;
static const int  DEFAULT_MAX_OLAP_BASES = 30;
static const double  LONG_ORF_SCORE_PER_BASE = 0.03;
  // artificially good score value for sufficiently long orfs


struct  Orf_Interval_t
  {
   int  lo, hi;
   int  frame;
   bool  deleted;
  };


static bool  Orf_Interval_Cmp
    (const Orf_Interval_t & a, const Orf_Interval_t & b)
  { return  (a . lo < b . lo || (a . lo == b . lo && a . hi < b . hi)); }


struct  Range_t
  {
   int  lo, hi;
  };


static bool  Range_Cmp
    (const Range_t & a, const Range_t & b)
  { return  (a . lo < b . lo || (a . lo == b . lo && a . hi < b . hi)); }


struct  Position_t
  {
   int  lo, hi, max_prev;
  };


struct  Start_t
  {
   int  j, pos, which;
   double  score, rate;
  };


// Function prototypes

static void  Echo_General_Settings
    (FILE * fp);
static void  Echo_Specific_Settings
    (FILE * fp, int len);
static void  Eliminate_Overlapping
    (vector <Orf_Interval_t> & orf_list, int max_olap);
static double  Entropy_Distance_Ratio
    (int start, int len, int fr);
static void  Entropy_Filter
    (vector <Orf_t> & orf_list, double cutoff);
static int  Find_Optimal_Len
    (const vector <Orf_Interval_t> & interval);
static void  Find_Orfs
    (vector <Orf_t> & orf_list);
static void  Finish_Orfs
    (bool use_wraparound, const int prev_rev_stop [3],
     const int last_rev_start [3], int last_position,
     vector <Orf_t> & orf_list);
static void  Get_Ignore_Regions
    (void);
static void  Get_Intervals
    (vector <Orf_Interval_t> & interval_list, const vector <Orf_t> & orf_list);
static void  Handle_First_Forward_Stop
     (int fr, int pos, int start_pos, int first_base, int & gene_len,
      int & orf_len, bool use_wraparound);
static void  Handle_Last_Reverse_Stop
     (int fr, const int prev_rev_stop [3], const int last_rev_start [3],
      int & gene_len, int & orf_len, bool use_wraparound, int last_position);
static int  Intersect_Size
    (int a, int b, int c, int d);
static int  On_Seq_0
    (int i);
static int  On_Seq_1
    (int i);
static void  Output_Orfs
    (FILE * fp, const vector <Orf_Interval_t> & interval, int & total_len);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Print_Comma_Separated_Strings
    (const vector <const char *> & v, FILE * fp);
static void  Read_Entropy_Profiles
    (const char * fn, bool & errflg);
static void  Remove_Shorter
    (vector <Orf_Interval_t> & interval, int len);
static void  Set_Start_And_Stop_Codons
    (void);
static void  Usage
    (void);
static void  Wrap_Around_Back
    (int wfr, int pos, int & gene_len, int & orf_len);
static void  Wrap_Through_Front
    (int fr, int pos, int & gene_len, int & orf_len);

#endif
