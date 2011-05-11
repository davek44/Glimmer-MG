//  A. L. Delcher
//
//  File:  glimmer3.hh
//
//  Last Modified:  Tue May  9 10:25:40 EDT 2006
//
//  Declarations for  Glimmer3


#ifndef  __GLIMMER3_HH_INCLUDED
#define  __GLIMMER3_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "gene.hh"
#include  "icm.hh"
#include <sys/stat.h>
#include <algorithm>

static void  All_Frame_Score
    (const string & s, int offset, int frame, vector <double> & af);
static void  Echo_General_Settings
    (FILE * fp);
static double  Entropy_Distance_Ratio
    (int start, int len, int fr);
static void  Find_Stops_Reverse
    (const string & s, int len, vector <bool> & has_stop);
static void  Fix_Wrap
    (int & p, const int n);
static void  Get_Orf_Pos_List
    (void);
static void  Integerize_Scores
    (const vector <double> ds, int hi_score, const vector <bool> set_zero,
    vector <int> & is);
static int  On_Seq_0
    (int i);
static void  Output_Extra_Start_Info
    (FILE * fp, int i, int lo, int hi, int frame,
     vector <Start_t> & start_list);
static void  Parse_Command_Line
    (int argc, char * argv []);
template  <class DT>
static void  Permute_By_Frame
    (vector <DT> & v, int frame);
static void  Print_Orflist_Headings
    (FILE * fp);
static void  Read_Entropy_Profiles
    (const char * fn, bool & errflg);
static void  Read_Sequences
    (FILE * fp, vector <string> & seq_list, vector <string> & hdr_list,
     int & seq_ct);
static void  Score_Orflist
    (FILE * detail_fp, FILE * summary_fp);
static void  Score_Orfs
    (vector <Orf_t> & orf_list, vector <Gene_t> & gene_list, FILE * fp);
static void  Score_Separate_Input
    (const string & seq, const string & hdr, int seq_num, FILE * detail_fp,
     FILE * predict_fp);
static void  Trace_Back
    (FILE * fp, const Event_Node_t & final_event);
static void  Usage
    (void);

#endif
