//  A. L. Delcher
//  D. R. Kelley
//
//  File:  glimmer_base.hh
//
//  Declarations for  Glimmer Base


#ifndef  __GLIMMER_BASE_HH_INCLUDED
#define  __GLIMMER_BASE_HH_INCLUDED


#include  "delcher.hh"
#include  "kelley.hh"
#include  "fasta.hh"
#include  "gene.hh"
#include  "icm.hh"
#include <sys/stat.h>
#include <algorithm>
#include <map>

// Default values of global variables

static const int  DEFAULT_MIN_GENE_LEN = 75;
static const int  DEFAULT_MAX_OLAP_BASES = 30;
static const int  DEFAULT_RIBOSOME_WINDOW_SIZE = 20;
static const double  DEFAULT_START_PROB []
     = {0.60, 0.30, 0.10};
static const int  DEFAULT_THRESHOLD_SCORE = 30;
static const float  DEFAULT_PRIOR = -1;
static const int  DEFAULT_USE_FIRST_START_CODON = false;
static const int  DEFAULT_USE_INDEPENDENT_SCORE = true;
static const int  HI_SCORE = 100;
  // the highest possible ICM score for an orf
static const double  LONG_ORF_SCORE_PER_BASE = 0.01; // 0.03;
  // artificially good score value for sufficiently long orfs
  //**ALD Should maybe change to a lower value like 0.01 ??

struct  Event_Node_t
  {
   int  id : 24;
   int  frame : 3;
   unsigned  is_first_start : 1;
   unsigned  disqualified : 1;
   unsigned  truncated : 1;
   Event_t  e_type;
   int  pos, pwm_sep;
     // pos is the last base of the codon, numbered starting at 1
   double  score, pwm_score;
   vector<Error_t> errors;
   Event_Node_t  * frame_pred;
   Event_Node_t  * best_pred;

   Event_Node_t  ()   // default constructor
     { is_first_start = disqualified = truncated = 0; }

   void  Set_Frame_From_Pos
       (void);
  };

static bool  Event_Pos_Cmp
    (Event_Node_t * const & a, Event_Node_t * const & b)
  { return  (a -> pos < b -> pos); }

struct  Range_t
  {
   int  lo, hi;
  };

struct  Orf_Pos_t
  {
   int  start, stop, dir;
   char  * tag;
  };

static bool  Range_Cmp
    (const Range_t & a, const Range_t & b)
  { return  (a . lo < b . lo); }

struct  Start_t
{
     int  j, pos;
     double  score, rate;
     int  which : 8;
     unsigned  truncated : 1;  
     bool  first;
     vector<Error_t> errors;
};

static bool Start_Cmp
	(const Start_t & a, const Start_t & b)
{ return a.pos < b.pos; };


struct vec_error_cmp
{
     bool operator()(const vector<Error_t> & v1, const vector<Error_t> & v2) const 
	  {
	       unsigned int v1_size = v1.size();
	       unsigned int v2_size = v2.size();

	       if(v1_size != v2_size)
		    return v1_size < v2_size;
	       else {	       
		    for(unsigned int i = 0; i < v1_size; i++) {
			 if(v1[i].pos != v2[i].pos)
			      return v1[i].pos < v2[i].pos;
			 else
			      if(v1[i].type != v2[i].type)
				   return v1[i].type < v2[i].type;
		    }

		    return 0;
	       }
	  }
};


extern int  Verbose;
extern bool Dave_Log;
extern bool Length_Log;
extern bool Detail_Log;
extern bool Sequence_Log;

extern bool Allow_Truncated_Orfs;
extern Event_Node_t * Best_Event[6];
extern string  Command_Line;
extern vector <double>  Cumulative_Score [6];
extern const char  * Fasta_Header;
extern Event_Node_t  First_Event, Final_Event;
extern vector <Codon_t>  Fwd_Start_Pattern;
extern vector <Codon_t>  Fwd_Stop_Pattern;
extern bool  GC_Frac_Set;
extern int  Genbank_Xlate_Code;
extern ICM_t  Gene_ICM;
extern int  Gene_ID_Ct;
extern bool  Genome_Is_Circular;
extern char  * ICM_File_Name;
extern char  * Ignore_File_Name;
extern int  Ignore_Score_Len;
extern vector <Range_t>  Ignore_Region;
extern double  Indep_GC_Frac;
extern ICM_t  Indep_Model;
extern Event_Node_t  * Last_Event [6];
extern PWM_t  LogOdds_PWM;
extern int  Min_Gene_Len;
extern int  Max_Olap_Bases;
extern double  Neg_Entropy_Profile [20];
extern int  Num_Start_Codons;
extern int  Num_Stop_Codons;
extern char  * Orflist_File_Name;
extern char  * Output_Tag;
extern double  Pos_Entropy_Profile [20];
extern vector <Codon_t>  Rev_Start_Pattern;
extern vector <Codon_t>  Rev_Stop_Pattern;
extern PWM_t  Ribosome_PWM;
extern int  Ribosome_Window_Size;
extern bool  Separate_Orf_Input;
extern vector <Orf_Pos_t>  Orf_Pos_List;
extern string  Sequence;
extern int  Sequence_Ct;
extern char  * Sequence_File_Name;
extern int  Sequence_Len;
extern vector <const char *>  Start_Codon;
extern vector <double>  Start_Prob;
extern vector <const char *>  Stop_Codon;
extern string  Tag;
extern int  Threshold_Score;
extern bool  Use_Entropy_Profiles;
extern bool  Use_First_Start_Codon;
extern bool  Use_Independent_Score;

extern bool Allow_Indels;
extern int Dist_Max_Overlap;
extern double Event_Threshold;
extern char * Feature_File;
extern AdjDist_Dist_t LogOdds_AdjDist;
extern AdjOr_Dist_t LogOdds_AdjOr;
extern float LogOdds_Fudge;
extern Length_Dist_t LogOdds_Length;
extern float LogOdds_Prior;
extern Start_Dist_t LogOdds_Start;
extern vector< pair<double,int> > Meta_PWM_Save;
extern vector<PWM_t> Meta_Ribosome_PWMs;
extern int Min_Indel_ORF_Len;
extern double Start_Threshold;
extern bool User_Adj;
extern bool User_Length;
extern bool User_RBS;
extern bool User_Start;



void  Add_Events_Fwd
    (const Orf_t & orf, vector <Start_t> & start_list, int & id);
void  Add_Events_Rev
    (const Orf_t & orf, vector <Start_t> & start_list, int & id);
void  Add_PWM_Score
    (Event_Node_t * p);
void Blend_Length
    (vector<double> & length_dist, vector<double> & par_length_dist, vector<double> & nonpar_length_dist, double par_cumprob);
void  Clear_Events
    (void);
void  Complement_Transfer
    (string & buff, const string & s, int lo, int hi);
void  Disqualify
    (Event_Node_t * p, int cutoff);
void  Do_Fwd_Stop_Codon
    (int i, int frame, int prev_fwd_stop [3], int first_fwd_start [3],
     int first_fwd_stop [3], int first_base, bool hit_ignore,
     vector <Orf_t> & orf_list);
void  Do_Rev_Stop_Codon
    (int i, int frame, int prev_rev_stop [3], int last_rev_start [3],
     bool hit_ignore, vector <Orf_t> & orf_list);
void  Echo_Specific_Settings
    (FILE * fp, int len);
int  Find_Uncovered_Position
    (vector <Event_Node_t *> ep);
void  Find_Orfs
    (vector <Orf_t> & orf_list);
void  Finish_Orfs
    (bool use_wraparound, const int prev_rev_stop [3],
     const int last_rev_start [3], int last_position,
     vector <Orf_t> & orf_list);
int  Frame_To_Sub
    (int f);
void  Get_Ignore_Regions
    (void);
vector<int> Get_Sequence_Lengths
    ();
void  Handle_First_Forward_Stop
    (int fr, int pos, int start_pos, int first_base, int & gene_len,
     int & orf_len, bool use_wraparound);
void  Handle_First_Reverse_Stop
    (int pos, int last_start, int & gene_len, int & orf_stop, bool hit_ignore);
void  Handle_Last_Reverse_Stop
    (int fr, const int prev_rev_stop [3], const int last_rev_start [3],
     int & gene_len, int & orf_len, int & orf_stop, bool use_wraparound, int last_position);
void  Initialize_Terminal_Events
    (Event_Node_t & first_event, Event_Node_t & final_event,
     Event_Node_t * best_event [6], Event_Node_t * last_event [6]);
double  Olap_Score_Adjustment
    (int lo, int hi, int f1, int f2);
int  On_Seq_1
    (int i);
void Parse_Features
    (const char* feature_file);
void Parse_Features_Err
    (string line);
int  Position_To_Frame
    (int p);
void  Print_Comma_Separated_Strings
    (const vector <const char *> & v, FILE * fp);
void  Print_Headings
    (FILE * fp);
void Print_Start
    (Start_t & start, int stop_pos, Event_Node_t * ne);
const char  * Print_String
    (Event_t e);
void  Prob_To_Logs
    (vector <double> & v);
void  Process_Events
    (void);
void  Process_Fwd_Start_Rev_Stop_Event
    (Event_Node_t * ep);
void  Process_Fwd_Stop_Rev_Start_Event
    (Event_Node_t * ep);
void  Process_Initial_Event
    (Event_Node_t * ep);
void  PWM_Meta_Score_Fwd_Start
    (int pos, double & score, int & separation);
void  PWM_Meta_Score_Rev_Start
    (int pos, double & score, int & separation);
void  PWM_Score_Fwd_Start
    (int pos, const PWM_t & pwm, int window, double & score, int & separation);
void  PWM_Score_Rev_Start
    (int pos, const PWM_t & pwm, int window, double & score, int & separation);
void Read_Dist_Dist
    (ifstream & ff, vector<float> & dist_dist);
float Read_Length_Dist
    (ifstream & ff, vector<double> & length_dist);
void Read_Orient_Dist
    (ifstream & ff, vector<float> & orients);
void Read_Start_Dist
    (ifstream & ff, vector<float> & starts);
void  Requalify
    (Event_Node_t * p, int cutoff);
void  Reverse_Complement_Transfer
    (string & buff, const string & s, int lo, int hi);
void  Reverse_Transfer
    (string & buff, const string & s, int start, int len);
void  Set_Final_Event
    (Event_Node_t & fe, Event_Node_t * best_event [6],
     int seq_len);
void  Set_GC_Fraction
    ();
void  Set_Ignore_Score_Len
    (void);
void  Set_Start_And_Stop_Codons
    (void);
void  Shift_Events
    (vector <Event_Node_t *> & ep, int reference_pos);
void  Show_Events
    (FILE * fp);
void  Wrap_Around_Back
    (int wfr, int pos, int & gene_len, int & orf_len);
void  Wrap_Through_Front
    (int fr, int pos, int & gene_len, int & orf_len);

#endif
