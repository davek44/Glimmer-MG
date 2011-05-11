//  A. L. Delcher
//
//  File:  icm.hh
//
//  Last Modified:  3 May 2004
//
//  Declarations for the Interpolated Context Model (ICM) used by Glimmer



#ifndef _ICM_H_INCLUDED
#define _ICM_H_INCLUDED


#include  "delcher.hh"
#include  "gene.hh"

using namespace std;


#define  STORE_MUT_INFO  1
  // If true, save best mutual information in training node and include it
  // in ascii dump of model


const int  ALPHABET_SIZE = 4;
  // Number of characters is strings being modelled
const int  ALPHA_SQUARED = (ALPHABET_SIZE * ALPHABET_SIZE);
  // Square of the preceding
const char  ALPHA_STRING [] = "acgt";
  // String of allowed characters in alphabet, used for converting
  // characters to subscripts

const int  NUM_CHI2_ENTRIES = 7;
  // Number of entries in following tables
const float  CHI2_VAL [NUM_CHI2_ENTRIES]
    = {2.37, 4.11, 6.25, 7.81, 9.35, 11.3, 12.8};
  // Table of chi-squared values
const float  CHI2_SIGNIFICANCE [NUM_CHI2_ENTRIES]
    = {0.50, 0.75, 0.90, 0.95, 0.975, 0.99, 0.995};
  // Corresponding significance values

const int  DEFAULT_MODEL_LEN = 12;
  // Number of bases in window used for scoring.  The last base
  // in this window is Markov dependent on the preceding bases
const int  ID_STRING_LEN = 150;
  // Length of string written at start of icm file
const int  DEFAULT_MODEL_DEPTH = 7;
  // The max number of bases to restrict (use) in the window
const unsigned  NUM_FIXED_LENGTH_PARAMS = 6;
  // The number of binary integer parameters at the start of the
  // fixed-length icm file
const double  MUT_INFO_BIAS = 0.03;
  // A mutual information position on the left (i.e., earlier in the
  // context window) must be at least this much better than one
  // to its right in order to be used.
const int  DEFAULT_PERIODICITY = 3;
  // Number of models to cycle through in each string.  Must be at least 1

const char  MAX_MI_CHAR = '*';
  // Indicates position in context of current node's max mutual info
  // Used for labels of ascii version of output
const char  SEPARATOR_CHAR = '|';
  // Used to separate groups of periodicity characters in ascii version
  // of output

const double  MAX_LOG_DIFF = -46.0;
  // Treat differences larger than this as essentially infinite
  // to avoid over/underflow problems.  This is approximately
  // 1e-20 difference in probs
const double  MUT_INFO_EPSILON = 1e-4;
  // Smaller than this is regarded as zero mutual information
const double  PSEUDO_COUNT = 0.001;
  // Part of a count inherited in proportion to the parent's probabilities
  // to prevent 0 counts.
const int  SAMPLE_SIZE_BOUND = 400;
  // If the counts for the bases sum to less than this value,
  // use chi^2 formula to weight the probability

const int  ICM_VERSION_ID = 200;
  // Integer version number stored in model prefix for compatibility check


#define  PARENT(x) ((int) ((x) - 1) / ALPHABET_SIZE) 
//  Macro that returns the parent of x in our tree



// Temporary so it will compile

enum  ICM_Model_t
    {UNKNOWN_TYPE};


struct  ICM_Training_Node_t
  {
   short int  mut_info_seq;
     // The base that is to be restricted in mut_info_pos
   int  (* count) [ALPHA_SQUARED];
     // count [x] [y] (where y is a pair of letters (e.g., aa, ac, ag, .. tt))
     // is the number of times the first base of y occurs at position x
     // and the last base of y occurs at position (model_len - 1)
  };


struct  ICM_Score_Node_t
  {
   short int  mut_info_pos;
#if  STORE_MUT_INFO
   float  mut_info;
#endif
   float  prob [ALPHABET_SIZE];
  };


class  ICM_t
  {
  protected:
   bool  empty;
   int  model_len;
     // Number of bases in window
   int  model_depth;
     // Most levels in model tree, i.e., most number of positions
     // to be dependent on
   int  periodicity;
     // Number of different models, alternating cyclicly
   int  num_nodes;
     // Number of nodes in tree
   ICM_Score_Node_t  * * score;

  public:
   ICM_t
       (int m = DEFAULT_MODEL_LEN,
        int d = DEFAULT_MODEL_DEPTH,
        int p = DEFAULT_PERIODICITY);
   ~ ICM_t
       ();

   int  Get_Model_Len
       (void)
     { return  model_len; }
   int  Get_Periodicity
       (void)
     { return  periodicity; }

   void  Build_Indep_WO_Stops
       (double gc_frac, const vector <const char *> & stop_codon);
   void  Build_Reverse_Codon_WO_Stops
    (double codon_prob [64], const vector <const char *> & stop_codon);
   void  Cumulative_Score
       (const string & s, vector <double> & score, int frame)  const;
   void  Cumulative_Score_String
       (char * string, int len, int frame, double * score);
   void  Display
       (FILE * fp);
   void  Frame_Score
       (const string & s, vector <double> & score, int frame)  const;
   void  Full_Window_Distrib
       (char * string, int frame, float * dist);
   double  Full_Window_Prob
       (const char * string, int frame)  const;
   void  Input
       (FILE * fp);
   void  Output
       (FILE * fp, bool binary_form);
   void  Output_Node
       (FILE * fp, ICM_Score_Node_t * node, int id, int frame, bool binary_form);
   double  Partial_Window_Prob
       (int predict_pos, const char * string, int frame)  const;
   void  Read
       (char * path);
   double  Score_String
       (const char * string, int len, int frame)  const;
   void  Set_Label_String
       (char * label, int id, int frame);
   void  Write_Header
       (FILE * fp, bool binary_form);
   void Copy
       (ICM_t & icm);
  };


class  ICM_Training_t  :  public ICM_t
  {
  private:
   ICM_Training_Node_t  * * train;
   int  (* count_memory) [ALPHA_SQUARED];
     // Holds dynamically allocated block for all counts to avoid
     // loss at tail from separate callocs.  Saved here to allow
     // freeing it in the destructor.

   void  Complete_Tree
       (const vector <char *> & data);
   void  Count_Char_Pairs_Restricted
       (const char * string, int level);
   ICM_Training_Node_t *  Get_Training_Node
       (const char * w, int frame, int level);
   void  Interpolate_Probs
       (int frame, int sub, int ct []);
   void  Take_Logs
       (void);

  public:
   ICM_Training_t
       (int m = DEFAULT_MODEL_LEN,
        int d = DEFAULT_MODEL_DEPTH,
        int p = DEFAULT_PERIODICITY);
   ~ ICM_Training_t
       ();

   void  Train_Model
       (const vector <char *> & data);
  };


class  Fixed_Length_ICM_t
  {
  private:
   int  length;
   int  max_depth;
   int  special_position;
   ICM_Model_t  model_type;
   int  * permutation;
   vector <ICM_t *> sub_model;

  public:
   Fixed_Length_ICM_t
       (int len = 1, int sp = 0, int * perm = NULL, ICM_Model_t = UNKNOWN_TYPE);
   ~ Fixed_Length_ICM_t
       ();

   int  Get_Length
       (void)
     {
      return  length;
     }
   double  Score_Window
       (char * w);

   // Match MDD_t interface
   int  getModelLength
       (void)
     {
      return  length;
     }
   int  getModelType
       (void)
     {
      return  model_type;
     }
   int  getSpecialPosition
       (void)
     {
      return  special_position;
     }
   void  read
       (const char * path);
   double  subrange_score
       (char * w, int lo, int hi);
   double  score
       (char * w)
     {
      return  Score_Window (w);
     }
  };


class  Fixed_Length_ICM_Training_t
  {
  private:
   int  length;
   int  max_depth;
   int  special_position;
   ICM_Model_t  model_type;
   int  * permutation;
   vector <ICM_Training_t *> sub_model;

  public:
   Fixed_Length_ICM_Training_t
       (int len, int md, int sp, int * perm, ICM_Model_t mt = UNKNOWN_TYPE);
   ~ Fixed_Length_ICM_Training_t
       ();
   void  Output
       (FILE * fp, bool binary_form);
   void  Train_Model
       (vector <char *> & data);
   void  Write_Header
       (FILE * fp, bool binary_form);
  };



void  Count_Char_Pairs
    (int ct [] [ALPHA_SQUARED], char * string, int w, int period);
void  Count_Single_Chars
    (int ct [ALPHABET_SIZE], char * string, int w, int period);
double  Get_Mutual_Info
    (int ct [], int n, int sum);
void  Permute_Data
    (vector <char *> & data, int * perm);
void  Permute_String
    (char * s, int * perm, int n);
int  Subscript
    (char ch);

#endif

