//  D. R. Kelley
//
//  File:  glimmer-mg.hh
//
//  Last Modified:  Tue May  9 10:25:40 EDT 2006
//
//  Declarations for  Glimmer-MG


#ifndef  __GLIMMERMG_HH_INCLUDED
#define  __GLIMMERMG_HH_INCLUDED


#include  "glimmer_base.hh"
#include  "delcher.hh"
#include  "fasta.hh"
#include  "gene.hh"
#include  "icm.hh"
#include <sys/stat.h>
#include <algorithm>
#include <map>

//-- Include hash_map
#ifdef __GNUC__
#if __GNUC__ < 3
  #include <hash_map.h>
  #include <hash_set.h>
  namespace Sgi { using ::hash_map; }; // inherit globals
  #define HASHMAP std
#elif __GNUC__ == 3
  #include <ext/hash_map>
  #include <ext/hash_set>
  #if __GNUC_MINOR__ == 0
    namespace Sgi = std;               // GCC 3.0
    #define HASHMAP std
  #else
    namespace Sgi = ::__gnu_cxx;       // GCC 3.1 and later
    #define HASHMAP __gnu_cxx
  #endif
#elif __GNUC__ > 3
  #include <ext/hash_map>
  #include <ext/hash_set>
  namespace Sgi = ::__gnu_cxx;         // GCC 4.0 and later
  #define HASHMAP __gnu_cxx
#endif
#else      // ...  there are other compilers, right?
  namespace Sgi = std;
  #define HASHMAP std
#endif

static string Classes_ICM_File
    (vector<string> & seq_classes);
static void Clean_Quality_454
    ();
static void Complement_Transfer_Qual
   (vector<int> & buff, int start, int len);
static void Cumulative_Frame_Score
    (int  frame, int lo, int hi, vector<double> & score, vector<double> & indep_score);
static void  Echo_General_Settings
    (FILE * fp);
static int Fwd_Prev_Stop
    (int end_point);
static void Parse_Classes
    (const char* class_file);
static void  Parse_Command_Line
    (int argc, char * argv []);
static double Pass_Stop_Penalty
    (int frame, int lo, int hi);
static void Read_Meta_ICMs
    ();
static void Read_Meta_RBS
    ();
static void Read_Meta_Lengths
    ();
static void Read_Meta_Starts
    ();
static void Read_Meta_Stops
    ();
static void Read_Meta_AdjOr
    ();
static void Read_Meta_AdjDist
    ();
static void Read_Meta_GC
    ();
static void  Reverse_Transfer_Qual
    (vector<int> & buff, int start, int len);
static int   Rev_Next_Stop
    (int end_point);
static void Save_Prev_Stops
    ();
static void  Score_All_Frames
    ();
static void Score_Indels
(Orf_t & orf, vector<Start_t> & start_list, vector<Error_t> & errors, double suffix_score, int suffix_j, vector<double> & score, vector<double> & indep_score, int q, int k, int j);
static void  Score_Orfs_Errors
    (vector <Orf_t> & orf_list, FILE * fp);
static void  Score_Orf_Starts
    (Orf_t & orf, vector<Start_t> & start_list, int end_point, double suffix_score, int suffix_j, vector<Error_t> & indels);
static void Set_Quality_454
    ();
static void  Trace_Back
    (FILE * fp, const Event_Node_t & final_event);
static void Update_Meta_Adj
    ();
static void Update_Meta_Null_ICM
    ();
static void Update_Meta_Length
    ();
static void Update_Meta_RBS
    ();
static void Update_Meta_Start
    ();
static void Update_Meta_Stop
    ();
static void  Usage
    (void);

#endif
