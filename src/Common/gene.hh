//  A. L. Delcher
//
//  File:  gene.hh
//
//  Last Modified:  23 October 2003
//
//  DNA- and gene-related routines delcarations



#ifndef  __GENE_HH_INCLUDED
#define  __GENE_HH_INCLUDED

#include  "xlate_tables.hh"


const unsigned  ATG_MASK = 0x184;
const unsigned  CAA_MASK = 0x211;
const unsigned  CAC_MASK = 0x212;
const unsigned  CAG_MASK = 0x214;
const unsigned  CAT_MASK = 0x218;
const unsigned  CAY_MASK = 0x21a;
const unsigned  CTA_MASK = 0x281;
const unsigned  CTG_MASK = 0x284;
const unsigned  GTG_MASK = 0x484;
const unsigned  RTG_MASK = 0x584;
const unsigned  TAA_MASK = 0x811;
const unsigned  TAG_MASK = 0x814;
const unsigned  TAR_MASK = 0x815;
const unsigned  TCA_MASK = 0x821;
const unsigned  TGA_MASK = 0x841;
const unsigned  TRA_MASK = 0x851;
const unsigned  TTA_MASK = 0x881;
const unsigned  TTG_MASK = 0x884;
const unsigned  TYA_MASK = 0x8a1;
const unsigned  YTA_MASK = 0xa81;
const unsigned  SHIFT_MASK = 0xFF;

const unsigned  DELETE_FLAG = 0x01;
const unsigned  TRUNCATED_START_FLAG = 0x02;
const unsigned  TRUNCATED_STOP_FLAG = 0x04;

const long int  INCR_SIZE = 10000;
const long int  INIT_SIZE = 10000;
const int  MAX_LINE = 300;

#define  DEFAULT_POS_ENTROPY_PROF  {0.08468,0.01606,0.05739,0.05752,0.04328,\
  0.07042,0.02942,0.05624,0.04442,0.05620,0.03029,0.03975,0.05116,0.04098,\
  0.05989,0.08224,0.05660,0.06991,0.02044,0.03310}
#define  DEFAULT_NEG_ENTROPY_PROF  {0.07434,0.03035,0.05936,0.04729,0.05662,\
  0.07704,0.05777,0.05328,0.03360,0.05581,0.01457,0.03718,0.04594,0.05977,\
  0.08489,0.05990,0.04978,0.07227,0.01050,0.01974}
const char  * const DEFAULT_START_CODON []
     = {"atg", "gtg", "ttg"};
const char  * const DEFAULT_STOP_CODON []
     = {"taa", "tag", "tga"};


enum  Event_t
  {INITIAL, FWD_START, FWD_STOP, REV_START, REV_STOP, TERMINAL};

class  Codon_t
  {
  private:
   static const unsigned  shift_mask = 0xff;
   static const unsigned  reverse_shift_mask = 0xff0;

   unsigned int  data;
     // Represent the codon as a 12-bit string.  Each character
     // is 4 bits, representing whether it can be a, c, g or t.
     // a is 1, c is 2, g is 4 and t is 8.
     // E.g., 'a' is 0001; IUPAC character 's' (which is 'c' or 'g')
     // is 0110.
   

  public:
   Codon_t ()
     { data = 0x0; }

   void  Clear
       ()
     { data = 0x0; }
   bool  Can_Be
       (const vector <Codon_t> & a, int & which);
   bool  Must_Be
       (const vector <Codon_t> & a, int & which);
   void  Print
       (FILE * fp)
     { fprintf (fp, "%03x", data); }
   void  Reverse_Complement
       (void);
   void  Reverse_Shift_In
       (char ch);
   void  Set_From
       (const char * s);
   void  Shift_In
       (char ch);
  };


class  Orf_t
  {
  protected:
   int  stop_position;
     // first base (i.e., lowest subscript) counting positions
     // starting at 1
   int  frame;
     // is determined by the leftmost position of the stop codon,
     // positions starting at 1, positive for forward, negative for
     // reverse
   int  orf_len;
   int  gene_len;

  public:
   Orf_t  ()
     { stop_position = 0;  frame = 0; }

   int  Get_Frame  (void)  const
     { return  frame; }
   int  Get_Gene_Len  (void)  const
     { return  gene_len; }
   int  Get_Orf_Len  (void)  const
     { return  orf_len; }
   int  Get_Stop_Position  (void)  const
     { return  stop_position; }

   void  Set_Frame  (int i)
     { frame = i; }
   void  Set_Gene_Len  (int i)
     { gene_len = i; }
   void  Set_Orf_Len  (int i)
     { orf_len = i; }
   void  Set_Stop_Position  (int i)
     { stop_position = i; }
     
  };

class Error_t
{
public:
     Error_t(int p, int t)
	  { pos = p; type = t; }

     int pos;
     int type; // 0 insertion, 1 deletion, 2 substitution
};

struct  DNA_vect_t
  {
   double  p [4];
  };


class  PWM_t
  {
  private:
   vector <DNA_vect_t>  col;

  public:
   PWM_t  ()
     {}

   void  Check  (void)
     { cerr << "PWM_t Check:  size = " << col . size () << endl; }
   void  Counts_To_Prob
       (void);
   double  Column_Score
       (char ch, int col)  const;
   bool  Is_Empty  (void)  const
     { return  col . empty (); }
   void  Make_Log_Odds_WRT_GC
       (double gc_frac);
   void  Print
       (FILE * fp);
   void  Probs_To_Logs
    (void);
   bool  Read
       (FILE * fp);
   int  Width  (void)  const
     { return   int (col . size ()); }

   PWM_t &  operator =
       (const PWM_t & src);
  };


class  Length_Dist_t
{
private:
     vector< vector <double> > Full_Log_Odds;
     vector< vector <double> > Trunc_Log_Odds;
     vector< vector <double> > Trunc2_Log_Odds;
     int min_aa_len;
     unsigned int full_trunc_merge[3];
     vector<double> fragment_lengths;

     void Choose_Frags(vector<int> & Frag_Lengths);
     int Choose_Frag_Dist(int frag_length) const;
     double Map_Length(int length) const;

public:
     Length_Dist_t  ();
     void  Check  (void)
	  { cerr << "Length_Dist_t Check:  size = " << Full_Log_Odds . size () << endl; }     
     double  Score (unsigned int length, bool truncated_5p, bool truncated_3p, unsigned int fragment_length)  const;
     double Huge_Score(unsigned int length, const vector<double> & My_Log_Odds) const;
     bool  Is_Empty  (void)  const
	  { return  Full_Log_Odds . empty (); }
     void  Make_Log_Odds (vector<double> & Gene_Lengths, vector<double> & Non_Lengths, vector<int> & Frag_Lengths, unsigned int min_gene_len);
     void  Print (const char * file_name) const;
     int  Width  (void)  const
	  { return   int (Full_Log_Odds . size ()); }
};


class  Start_Dist_t
{
private:
     vector <float> Log_Odds;
     const double* default_start_prob;

public:
     Start_Dist_t  (const double*);
     float  Score (int start)  const
	  { return Log_Odds[start]; }
     bool  Is_Empty  (void)  const
	  { return  Log_Odds . empty (); }
     void  Make_Log_Odds (vector<float> & Gene_Starts, vector<float> & Non_Starts);
     void  Print (const char * file_name) const;
};


class  AdjOr_Dist_t
{
private:
     float Log_Odds_Fwd_Fwd;
     float Log_Odds_Fwd_Rev;
     float Log_Odds_Rev_Fwd;
     float Log_Odds_Rev_Rev;

public:
     AdjOr_Dist_t  ();
     float  Score (int or1, int or2)  const;
     float  Score (Event_t & e1, Event_t & e2)  const;
     void  Make_Log_Odds (vector<float> & Gene_AdjOr, vector<float> & Non_AdjOr);
     void  Print (const char * file_name) const;
};


class  AdjDist_Dist_t
{
private:
     int max_overlap;
     vector<float> Log_Odds_Fwd_Fwd;
     vector<float> Log_Odds_Fwd_Rev;
     vector<float> Log_Odds_Rev_Fwd;
     void Make_Log_Odds(vector<float> & Log_Odds, vector<float> & Gene_AdjDist, vector<float> & Non_AdjDist);
     void Overlap_Smooth(vector<float> & dists);

public:
     AdjDist_Dist_t ();
     void Set_MaxOverlap(int mo);
     float  Score (Event_t & e1, Event_t & e2, int length)  const;
     void  Make_Log_Odds_Fwd_Fwd (vector<float> & Gene_AdjDist, vector<float> & Non_AdjDist);
     void  Make_Log_Odds_Fwd_Rev (vector<float> & Gene_AdjDist, vector<float> & Non_AdjDist);
     void  Make_Log_Odds_Rev_Fwd (vector<float> & Gene_AdjDist, vector<float> & Non_AdjDist);
     void  Print (const char * file_name) const;
};


class  Gene_t  :  public Orf_t
  {
  private:
   unsigned int  status;
   int  id;
   double  score;
   vector<Error_t> errors;

  public:
   Gene_t  ()
     { status = 0; }
   Gene_t  (const Orf_t & orf) : Orf_t (orf)
     { status = 0; }

   int  Get_ID  (void)  const
     { return  id; }
   double  Get_Score  (void)  const
     { return  score; }
   unsigned int  Get_Status  (void)  const
     { return  status; }
   unsigned int  Get_Status_Bit
       (unsigned int u)  const;
   vector<Error_t> Get_Errors (void)
     { return errors; }

   void  Set_ID  (int i)
     { id = i; }
   void  Set_Score  (double d)
     { score = d; }
   void  Set_Status  (unsigned int u)
     { status = u; }
   void  Set_Status_Bit  (unsigned int u)
     { status |= u; }
   void Set_Errors (vector<Error_t> ers)
     { errors = ers; }

   void  Clear_Status  (void)
     { status = 0; }
  };



bool  By_ID
    (const Gene_t & a, const Gene_t & b);
unsigned  Ch_Mask
    (char);
int  Char_Sub
    (char ch);
char  Codon_Translation
    (const char * c, int transl_tabl = 1);
char  Complement
    (char ch);
void  Counts_To_Entropy_Profile
    (int count [26], double ep [20]);
int  Filter
    (char ch);
void  Find_Stop_Codons
    (const char * s, int t, int stop []);
int  First_In_Frame_Stop
    (char * s, int frame);
void  Forward_Strand_Transfer
    (string & t, const string & s, int start, int len);
int  Is_Forward_Start
    (unsigned codon);
int  Is_Forward_Stop
    (unsigned codon);
int  Is_Reverse_Start
    (unsigned codon);
int  Is_Reverse_Stop
    (unsigned codon);
int  Is_Start
    (const char * s);
int  Is_Stop
    (const char * s);
int  Nucleotide_To_Subscript
    (char ch);
int  Read_String
    (FILE * fp, char * & t, long int & size, char name [], int partial);
void  Reverse_Complement
    (char * s);
void  Reverse_Complement
    (string & s);
void  Reverse_Strand_Transfer
    (string & t, const string & s, int start, int len);
void  Set_Stop_Codons_By_Code
    (vector <const char *> & stop_codon, int code, bool & errflg);


#endif
