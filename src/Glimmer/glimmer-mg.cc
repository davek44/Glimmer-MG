//  D. R. Kelley
//
//  File:  glimmer_mg.cc
//
//  Copyright (c) 2011 University of Maryland Center for Bioinformatics
//  & Computational Biology

#include "glimmer_base.hh"
#include "glimmer-mg.hh"


bool  Allow_Truncated_Orfs = true;
  // If set true by -X option, then score orfs that
  // extend to the end of the sequence
Event_Node_t  * Best_Event [6];
  // Best parse event up to the current point in each reading frame
string  Command_Line;
  // Command, options and parameters that invoked the program
const char  * Fasta_Header;
  // Header on first line of fasta input file
Event_Node_t  First_Event, Final_Event;
  // First and last nodes in DAG of possible parse events
vector <Codon_t>  Fwd_Start_Pattern;
  // Bit patterns representing possible forward start codons
vector <Codon_t>  Fwd_Stop_Pattern;
  // Bit patterns representing possible forward stop codons
bool  GC_Frac_Set = false;
  // If true, then  Indep_GC_Frac  is set by -C option; otherwise,
  // it is determined from the input sequence data.
int  Genbank_Xlate_Code = 0;
  // Holds the Genbank translation table number that determines
  // stop codons and codon translation.
ICM_t  Gene_ICM;
  // The interpolated context model (ICM) of the coding
  // part of genes.
int  Gene_ID_Ct = 0;
  // Counter used to assign ID numbers to tentative genes
bool  Genome_Is_Circular = false;
  // If true, input sequences are assumed to be circularly connected
  // so genes will be allowed to wrap around the end
char  * ICM_File_Name = NULL;
  // Name of the file containing the probability model
int  Ignore_Score_Len = INT_MAX;
  // Genes at least this long do not count the independent model
  // in their score
double  Indep_GC_Frac = -1.0;
  // GC proportion used in simple independent model.
  // Set from counts of input sequences or by -C option
ICM_t  Indep_Model (3, 2, 3);
  // The ICM for an independent model of bases, based on GC-percentage
  // but without in-frame stop codons
Event_Node_t  * Last_Event [6];
  // Last parse event up to the current point in each reading frame
PWM_t  LogOdds_PWM;
  // Log odds wrt background gc-fraction of  Ribosome_PWM.
int  Min_Gene_Len = DEFAULT_MIN_GENE_LEN;
  // Shortest (in nucleotides) gene that will be considered for scoring
int  Max_Olap_Bases = DEFAULT_MAX_OLAP_BASES;
  // Overlaps of this many or fewer bases are allowed between adjacent
  // genes
double  Neg_Entropy_Profile [20] = DEFAULT_NEG_ENTROPY_PROF;
  // Entropy distribution of amino-acids in non-genes
int  Num_Start_Codons;
  // Number of different start codon patterns
int  Num_Stop_Codons;
  // Number of different stop codon patterns
char  * Output_Tag = NULL;
  // Prefix used for output files
double  Pos_Entropy_Profile [20] = DEFAULT_POS_ENTROPY_PROF;
  // Entropy distribution of amino-acids in genes
vector <Codon_t>  Rev_Start_Pattern;
  // Bit patterns representing possible reverse start codons
vector <Codon_t>  Rev_Stop_Pattern;
  // Bit patterns representing possible reverse stop codons
PWM_t  Ribosome_PWM;
  // Position weight matrix for the ribosome binding pattern
int  Ribosome_Window_Size = DEFAULT_RIBOSOME_WINDOW_SIZE;
  // Width of window before starts in which to look for matches to
  //  Ribosome_PWM .
vector <Orf_Pos_t>  Orf_Pos_List;
  // List of orfs specified by the -L option to be scored separatedly
string  Sequence;
  // The input sequence to be scored.
int  Sequence_Ct;
  // The number of sequences in the input fasta file
char  * Sequence_File_Name = NULL;
  // Name of the input sequence file
int  Sequence_Len;
  // Length of genomic sequence string being processed.
vector <const char *>  Start_Codon;
  // Sequences assumed to be start codons
vector <const char *>  Stop_Codon;
  // Sequences assumed to be stop codons
string  Tag;
  // The fasta-header lines of the sequence in  Sequence

////////////////////////////////////////////
// Dave's variables
////////////////////////////////////////////
bool Allow_Indels = false;
  // Allow Glimmer to shift the frame of a gene at low quality bases
int Dist_Max_Overlap = -1;
  // Overlaps of this many or fewer bases are allowed between adjacent
  // genes (as defined by the adjacent distance distributions
double Event_Threshold = -3;
  // Minimum score for an ORF to be added as an event
char * Feature_File = NULL;
  // Name of file used to read in other models
AdjDist_Dist_t LogOdds_AdjDist;
  // Log odds ratio between gene adjacent distances and noncoding orf adjacent distances
AdjOr_Dist_t LogOdds_AdjOr;
  // Log odds ratio between gene adjacent orientation and noncoding orf adjacent orientation
float LogOdds_Fudge = 1.0;
  // Fudge factor to be added to the LogOdds_Prior
Length_Dist_t LogOdds_Length;
  // Log odds ratio between gene length and noncoding orf length
float LogOdds_Prior = DEFAULT_PRIOR;
  // Log odds ratio between genes and noncoding orf counts  
Start_Dist_t LogOdds_Start(DEFAULT_START_PROB);
  // Log odds ratio between gene start codon and noncoding orf start codon
double Start_Threshold = -6;
  // Minimum score for an ORF to be considered by Add_Events
bool User_Adj = false;
bool User_Length = false;
bool User_RBS = false;
bool User_Start = false;

static int Chunk_Sequences = 500000;
  // Number of sequences to allow in a single chunk
static char * Quality_File_Name = NULL;
  // Quality value file name
static double Indel_Suffix_Score_Threshold = -12;
  // Minimum score for an ORF to consider allowing indels
static bool Allow_Subs = false;
  // Allow Glimmer to pass through (error-predicted) stop codons
static int Indel_Quality_Threshold = 18;
  // Threshold for a base's quality value to allow an indel
static int Indel_Max = 2;
  // Number of indels per ORF to allow
vector< vector<double> > Frame_Scores(6);
  // ICM scores for all read bases in all read frames
static vector<int> Quality_Values;
static bool User_ICM = false;
static bool User_Stop = false;
  // Indicator variables to decide what the user has specified versus
  // what should be taken from classifications
static string ICM_dir = "/fs/szasmg3/dakelley/phymm3/.genomeData";
  // ICM directories
static vector<int> Fwd_Prev_Stops;
  // Saved positions of the previous forward stop codon for each sequence position
static vector<int> Rev_Next_Stops;
  // Saved positions of the next reverse stop codon for each sequence position


namespace HASHMAP
{
  template<> struct hash< std::string >
  {
    size_t operator()( const std::string& x ) const
    {
      return HASHMAP::hash< const char* >()( x.c_str() );
    }
  };
}

typedef HASHMAP::hash_map< string, vector<string> > class_hash;
static class_hash classifications;
  // Phymm classification for each sequence that maps to Phymm genomes used to build models
  // Path to directory with gene and noncoding ICM's
typedef HASHMAP::hash_map< string, ICM_t* > icm_hash;
static icm_hash Sequence_ICMs;
  // ICM's corresponding to classifications of sequences
typedef HASHMAP::hash_map< string, PWM_t* > rbs_hash;
static rbs_hash Sequence_RBSs;
  // RBS PWMs corresponding to classifications of sequences
typedef HASHMAP::hash_map< string, vector<double> > length_hash;
static length_hash Sequence_Lengths_Gene;
static length_hash Sequence_Lengths_Non;
  // Lengths corresponding to classifications of sequences
typedef HASHMAP::hash_map< string, float > prior_hash;
static prior_hash Sequence_Prior;
  // Prior probability log likelihood ratio to start every gene with
typedef HASHMAP::hash_map< string, vector<float> > start_hash;
static start_hash Sequence_Starts_Gene;
static start_hash Sequence_Starts_Non;
  // Start codons corresponding to classifications of sequences
typedef HASHMAP::hash_map< string, vector<float> > adjor_hash;
static adjor_hash Sequence_AdjOr_Gene;
static adjor_hash Sequence_AdjOr_Non;
  // Adjacent orientations corresponding to classifications of sequences
typedef HASHMAP::hash_map< string, vector<float> > adjdist_hash;
static adjdist_hash Sequence_AdjDist_ff_Gene;
static adjdist_hash Sequence_AdjDist_ff_Non;
static adjdist_hash Sequence_AdjDist_fr_Gene;
static adjdist_hash Sequence_AdjDist_fr_Non;
static adjdist_hash Sequence_AdjDist_rf_Gene;
static adjdist_hash Sequence_AdjDist_rf_Non;
  // Adjacent distances corresponding to classifications of sequences
typedef HASHMAP::hash_map< string, float > gc_hash;
static gc_hash Sequence_GC;
  // GC% corresponding to classifications of sequences
typedef HASHMAP::hash_map< string, int > transl_hash;
static transl_hash Sequence_Transl;
  // GenBank trans_table code corresponding to classifications of sequences

typedef HASHMAP::hash_set< string > seq_set;
typedef HASHMAP::hash_map< string, seq_set > icm_set_hash;
static icm_set_hash ICM_Sequences;
  // Maps ICM names to sets of reads that will use them


int  main
(int argc, char * argv [])

{
     FILE  * sequence_fp, * quality_fp, * detail_fp, * predict_fp;
     vector <string>  seq_list, hdr_list;
     vector < vector<int> > qual_list;
     vector <Orf_t>  orf_list;
     vector <Gene_t>  gene_list;
     string  hdr, filename;
     time_t  now;
     int  i, chunk_i;
     icm_set_hash::const_iterator icm_map_it;
     seq_set::const_iterator seq_set_it;
     string header, header_prefix;
     int seqs_analyzed = 0;     

     try
     {
	  now = time (NULL);
	  cerr << "Starting at " << ctime (& now) << endl;

	  Verbose = 0;

	  Parse_Command_Line (argc, argv);

	  Set_Start_And_Stop_Codons ();

	  // open output files
	  if(Detail_Log) {
	       filename = Output_Tag;
	       filename . append (".detail");
	       detail_fp = File_Open (filename, "w", __FILE__, __LINE__);

	       Echo_General_Settings (stderr);
	       fprintf (detail_fp, "Command:  %s\n\n", Command_Line . c_str ());
	       Echo_General_Settings (detail_fp);
	  } else
	       detail_fp = NULL;

	  filename = Output_Tag;
	  filename . append (".predict");
	  predict_fp = File_Open (filename, "w", __FILE__, __LINE__);

	  // prepare other models  
	  if(Feature_File != NULL)
	       Parse_Features(Feature_File);

	  if(!classifications.empty()) {
	       if(!User_Length)
		    Read_Meta_Lengths();
	       if(!User_Start)
		    Read_Meta_Starts();
	       if(!User_Adj) {
		    Read_Meta_AdjOr();
		    Read_Meta_AdjDist();
	       }
	       if(!User_Stop)
		    Read_Meta_Stops();
	  }
	  
	  // prepare ICM(s)
	  if(User_ICM) {
	       // get GC
	       if  (! GC_Frac_Set)
		    Set_GC_Fraction ();
	       // make null ICM
	       Indep_Model . Build_Indep_WO_Stops (Indep_GC_Frac, Stop_Codon);
	       Set_Ignore_Score_Len ();

	       // get gene ICM
	       //Gene_ICM . Read (ICM_File_Name);

	       // set up empty hash map to read it
	       seq_set dummy;
	       ICM_Sequences[string(ICM_File_Name)] = dummy;

	  } else if(!classifications.empty()) {
	       // use classification ICMs
	       Read_Meta_ICMs();
	       Read_Meta_GC();

	       // Gene_ICM and Indep_Model will be updated on the fly
	  } else {
	       sprintf (Clean_Exit_Msg_Line, "ERROR:  Must specify ICM with -m or sequence classifications with -c\n");
	       Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
	  }

	  // prepare RBS PWM(s)
	  if(User_RBS) {
	       // get GC
	       if  (! GC_Frac_Set)
		    Set_GC_Fraction ();

	       LogOdds_PWM = Ribosome_PWM;
	       LogOdds_PWM.Make_Log_Odds_WRT_GC (Indep_GC_Frac);
	  } else if(!classifications.empty()){
	       // use classification RBS PWMs
	       Read_Meta_RBS();
	  }

	  // open files
	  sequence_fp = File_Open (Sequence_File_Name, "r", __FILE__, __LINE__);
	  if(Quality_File_Name != NULL)
	       quality_fp = File_Open (Quality_File_Name, "r", __FILE__, __LINE__);

	  bool end_of_reads = false;
	  while (!end_of_reads)
	  {
	       // read in chunk of sequences
	       for(chunk_i = 0; !end_of_reads && chunk_i < Chunk_Sequences; chunk_i++)
	       {
		    // make lists big enough
		    if(chunk_i >= (int)seq_list.size())
		    {
			 seq_list.push_back("");
			 hdr_list.push_back("");
			 vector<int> q;
			 qual_list.push_back(q);
		    }

		    if(!Fasta_Read(sequence_fp, seq_list[chunk_i], hdr_list[chunk_i]))
			 end_of_reads = true;
		    if(Quality_File_Name != NULL)
			 Fasta_Qual_Vec_Read(quality_fp, qual_list[chunk_i], header);
	       }
	       if(end_of_reads)
		    chunk_i--;

	       // iterate over ICMs to be loaded
	       for(icm_map_it = ICM_Sequences.begin(); icm_map_it != ICM_Sequences.end(); icm_map_it++)
	       {
		    // read it here meta or not
		    Gene_ICM.Read((char*)icm_map_it->first.c_str());

		    for(i = 0; i < chunk_i; i++)
		    {
			 Fasta_Header = hdr_list[i].c_str();
			 header_prefix = split(string(Fasta_Header))[0];

			 // does this sequence match this ICM
			 if(!User_ICM)
			      seq_set_it = icm_map_it->second.find(header_prefix);

			 if(User_ICM || seq_set_it != icm_map_it->second.end())
			 {
			      // prepare sequence
			      Sequence = seq_list[i];
			      Sequence_Len = Sequence . length ();
			      for(int seq_i = 0; seq_i < Sequence_Len; seq_i++)
				   Sequence[seq_i] = tolower(Filter(Sequence[seq_i]));

			      // prepare quality values
			      if(Allow_Indels)
			      {
				   Quality_Values = qual_list[i];
				   if(Quality_File_Name == NULL)
					Set_Quality_454();
				   else
					Clean_Quality_454();
			      }

			      if(Detail_Log){
				   fprintf (detail_fp, "\n\n>%s\n", Fasta_Header);
				   Echo_Specific_Settings (detail_fp, Sequence_Len);
			      }
			      fprintf (predict_fp, ">%s\n", Fasta_Header);
			 
			      // update classification-based models	      
			      if(!classifications.empty())
			      {
				   if(!User_RBS)
					Update_Meta_RBS();
				   if(!User_Length)
					Update_Meta_Length();
				   if(!User_Start)
					Update_Meta_Start();
				   if(!User_Adj)
					Update_Meta_Adj();
				   if(!User_Stop)
					Update_Meta_Stop();
				   if(!User_ICM)
					Update_Meta_Null_ICM();
			      }

			      Initialize_Terminal_Events (First_Event, Final_Event, Best_Event, Last_Event);

			      if(Detail_Log)
				   Print_Headings (detail_fp);

			      if(Sequence_Log) {
				   cerr << "Analyzing Sequence #" << ++seqs_analyzed << " " << Fasta_Header << endl;
				   cerr << "Start Find_Orfs" << endl;
			      }
			      Find_Orfs (orf_list);

			      if(Sequence_Log)
				   cerr << "Start Score_Orfs" << endl;
			      //Score_Orfs (orf_list, gene_list, detail_fp);
			      Score_Orfs_Errors (orf_list, detail_fp);

			      if  (Verbose > 1)
				   Show_Events (stdout);

			      if(Sequence_Log)
				   cerr << "Start Process_Events" << endl;
			      Process_Events ();
			      Set_Final_Event (Final_Event, Best_Event, Sequence_Len);

			      if(Sequence_Log)
				   cerr << "Start Trace_Back" << endl;
			      Trace_Back (predict_fp, Final_Event);

			      gene_list . clear ();
			      orf_list . clear ();

			      Clear_Events ();
			 }			      
		    }
	       }
	  }

	  fclose (sequence_fp);
	  if(Quality_File_Name != NULL)
	       fclose(quality_fp);

	  if(Detail_Log)
	       fclose (detail_fp);
	  fclose (predict_fp);	  
     }
     catch (std :: exception & e)
     {
	  cerr << "** Standard Exception **" << endl;
	  cerr << e << endl;
	  exit (EXIT_FAILURE);
     }
     
     return  0;
}


static string Classes_ICM_File(vector<string> & seq_classes)

// Given the sequence's classifications, find the file name
// for the best ICM.

{
     struct stat st_file_info;
     vector<string> strain_nc1, strain_nc2;
     string icm_file;

     if(seq_classes.size() >= 2) {

	  // find best existing double
	  bool icm2_found = false;
	  for(unsigned int i = 1; !icm2_found && i < seq_classes.size(); i++) {	
	       if(seq_classes[0].compare(seq_classes[i]) < 0) {
		    strain_nc1 = split(seq_classes[0], '|');
		    strain_nc2 = split(seq_classes[i], '|');
	       } else {
		    strain_nc2 = split(seq_classes[0], '|');
		    strain_nc1 = split(seq_classes[i], '|');
	       }
	       
	       icm_file = ICM_dir + "/" + strain_nc1[0] + "/" + strain_nc1[1] + "_2/" + strain_nc2[0] + "/" + strain_nc2[1] + ".gicm";

	       if(stat(icm_file.c_str(), &st_file_info) == 0)
		    icm2_found = true;
	  }

	  if(!icm2_found) {
	       // or just use single
	       strain_nc1 = split(seq_classes[0], '|');
	       icm_file = ICM_dir + "/" + strain_nc1[0] + "/" + strain_nc1[1] + ".gicm"; 
	  }

     } else {

	  // construct file name
	  strain_nc1 = split(seq_classes[0], '|');
	  icm_file = ICM_dir + "/" + strain_nc1[0] + "/" + strain_nc1[1] + ".gicm"; 
     }

     return icm_file;
}


static void Clean_Quality_454()

// "Clean" the quality values for a 454 sequence
// where all quality values in a homopolymer run 
// are the same so that only the last nt has the
// quality value of interest.

{
     unsigned int i;

     // Can't be zero
     for(i = 0; i < Quality_Values.size(); i++)
	  if(Quality_Values[i] <= 0)
	       Quality_Values[i] = 1;

     if(Sequence_Len != (int)Quality_Values.size()) {
	  sprintf (Clean_Exit_Msg_Line, "ERROR:  %s sequence length does not match quality values length\n",Fasta_Header);
	  Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
     }

     // Use only the last qv of a homopolymer run
     for(i = 1; i < (unsigned int)Sequence_Len; i++) {
	  if(Sequence[i] == Sequence[i-1]) {
	       // set middling nt's to high quality
	       Quality_Values[i-1] = Max(Quality_Values[i-1], Indel_Quality_Threshold+1);
	  }
     }
}


static void Complement_Transfer_Qual (vector<int> & buff, int start, int len)

//  Copy to string  buff  the substring of Quality_Values  starting at subscript
//   start  and going to the right for a length of  len .

{
   buff . resize (len);
   for (int j = 0;  j < len;  j ++, start ++)
	buff [j] = Quality_Values [start];
}


static void Cumulative_Frame_Score
    (int  frame, int lo, int hi, vector<double> & score, vector<double> & indep_score)

// Fill in score and indep_score using the precomputed Frame_Scores from lo
// to hi in the frame given

{
     double cum_score = 0;
     int f = 1;
     int len = hi - lo;
     int si;

     if(frame > 0) {
	  score.resize(len);
	  indep_score.resize(len);
	  si = hi-1;
	  for(int i = 0; i < len; i++) {
	       score[i] = cum_score + Frame_Scores[f][si];
	       cum_score = score[i];
	       indep_score[i] = 0;

	       si--;
	       if(f == 2)
		    f = 0;
	       else
		    f++;
	  }
     } else {
	  score.resize(len);
	  indep_score.resize(len);
	  si = lo-1;
	  for(int i = 0; i < len; i++) {
	       score[i] = cum_score + Frame_Scores[3+f][si];
	       cum_score = score[i];
	       indep_score[i] = 0;

	       si++;
	       if(f == 2)
		    f = 0;
	       else
		    f++;
	  }
     }
}


static void  Echo_General_Settings
    (FILE * fp)

//  Output values of global variables and parameter settings
//  to  fp .

  {
   int  i, n;

   fprintf (fp, "Sequence file = %s\n", Sequence_File_Name);
   fprintf (fp, "Number of sequences = %d\n", Sequence_Ct);
   fprintf (fp, "ICM model file = %s\n", ICM_File_Name);

   fprintf (fp, "Circular genome = %s\n", Printable (Genome_Is_Circular));

   fprintf (fp, "Truncated orfs = %s\n", Printable (Allow_Truncated_Orfs));
   fprintf (fp, "Minimum gene length = %d bp\n", Min_Gene_Len);
   fprintf (fp, "Maximum overlap bases = %d\n", Max_Olap_Bases);
   if  (Genbank_Xlate_Code != 0)
	fprintf (fp, "Translation table = %d\n", Genbank_Xlate_Code);
   fprintf (fp, "Start codons = ");
   Print_Comma_Separated_Strings (Start_Codon, fp);
   fputc ('\n', fp);
   fprintf (fp, "Stop codons = ");
   Print_Comma_Separated_Strings (Stop_Codon, fp);
   fputc ('\n', fp);

   fprintf (fp, "GC percentage = %.1f%%\n", 100.0 * Indep_GC_Frac);
   fprintf (fp, "Ignore score on orfs longer than %s\n",
            Num_Or_Max (Ignore_Score_Len));

   return;
  }


static int Fwd_Prev_Stop(int end_point)

// Walk back in the Sequence from the given point
// until we hit a forward stop codon or the end
// of the sequence

{
     if(end_point >= 0 && end_point < Sequence_Len)
	  return Fwd_Prev_Stops[end_point];
     else
	  return end_point;

/*	  
	  Codon_t codon;
	  int i = end_point;
	  int c, which;

	  while(i >= 2) {
	       // get next codon
	       for(c = 0; c < 3; c++) {
		    codon.Reverse_Shift_In(Sequence[i]);
		    i--;
	       }

	       if(codon. Must_Be (Fwd_Stop_Pattern, which))
		    return i+3;
	  }

	  return i;
*/
}


static void Save_Prev_Stops()

// Fill in vectors mapping sequence positions to their
// previous stop codons

{
     // Forward

     Fwd_Prev_Stops.resize(Sequence_Len);
     
     int last_stops[3] = {0,1,-1};
     int frame = 0;
     Codon_t codon;
     int which;

     for(int i = 0; i < Sequence_Len; i++) {
	  // add next nt
	  codon.Shift_In(Sequence[i]);

          // if stop codon
	  if(i >= 2 && codon.Must_Be(Fwd_Stop_Pattern, which))
	       last_stops[frame] = i;

	  // set stop
	  Fwd_Prev_Stops[i] = last_stops[frame];

	  // update frame
	  frame = (frame+1)%3;
     }

     // Reverse

     Rev_Next_Stops.resize(Sequence_Len);

     last_stops[0] = Sequence_Len-1;
     last_stops[1] = Sequence_Len-2;
     last_stops[2] = Sequence_Len;
     
     frame = 0;

     for(int i = Sequence_Len-1; i >= 0; i--) {
	  // add next nt
	  codon.Shift_In(Complement(Sequence[i]));

	  // if stop codon
	  if(i <= Sequence_Len-3 && codon.Must_Be(Fwd_Stop_Pattern, which))
	       last_stops[frame] = i;

	  // set stop
	  Rev_Next_Stops[i] = last_stops[frame];

	  // update frame
	  frame = (frame+1)%3;
     }
}


static void Parse_Classes(const char* class_file)

// Parse metagenomic sequence classifications

{
     string line;
     vector<string> a;

     ifstream class_in(class_file);
     if(!class_in.good()) {
	  sprintf (Clean_Exit_Msg_Line, "ERROR:  Cannot open classification file %s\n",class_file);
	  Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
     }

     while(getline(class_in, line)) {
	  a = split(line);

	  vector<string> v;
	  for(unsigned int i = 1; i < a.size(); i++)
	       v.push_back(a[i]);

	  classifications[a[0]] = v;
     }
     class_in.close();
}



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   FILE  * fp;
   char  * p, * q;
   bool  errflg = false;
   int  i, ch;

   optarg = NULL;
   Command_Line = argv [0];

#if  ALLOW_LONG_OPTIONS
   int  option_index = 0;
   static struct option  long_options [] = {
        {"rbs_pwm", 1, 0, 'b'},
	{"class", 1, 0, 'c'},
	{"features", 1, 0, 'F'},
        {"gene_len", 1, 0, 'g'},
        {"help", 0, 0, 'h'},
	{"indel", 0, 0, 'i'},
	{"icm", 1, 0, 'm'},
        {"max_olap", 1, 0, 'o'},
	{"quality", 1, 0, 'q'},
	{"circular", 0, 0, 'r'},
	{"sub", 0, 0, 's'},
	{"fudge", 1, 0, 'u'},
        {"trans_table", 1, 0, 'z'},
        {"stop_codons", 1, 0, 'Z'},
        {0, 0, 0, 0}
      };

   while  (! errflg && ((ch = getopt_long (argc, argv,
        "b:c:f:g:him:o:P:q:rsu:z:Z:",
        long_options, & option_index)) != EOF))
#else
   while  (! errflg && ((ch = getopt (argc, argv,
        "b:c:f:g:him:o:P:q:rsu:z:Z:")) != EOF))
#endif

     switch  (ch)
       {
        case  'b' :
          Command_Line . append (" -b ");
          Command_Line . append (optarg);
          fp = File_Open (optarg, "r", __FILE__, __LINE__);
          Ribosome_PWM . Read (fp);
	  fclose(fp);
          Ribosome_PWM . Counts_To_Prob ();
          Ribosome_PWM . Probs_To_Logs ();
          if  (Verbose > 1)
              Ribosome_PWM . Print (stderr);
          User_RBS = true;
          break;

       case 'c' :
	    Command_Line.append(" -c");
	    Command_Line.append(optarg);	    
	    Parse_Classes(optarg);
	    break;

       case 'f' :
	    Command_Line.append(" -f");
	    Command_Line.append(optarg);
	    Feature_File = optarg;
	    break;

       case  'g' :
          Command_Line . append (" -g ");
          Command_Line . append (optarg);
          Min_Gene_Len = strtol (optarg, & p, 10);
          if  (p == optarg || Min_Gene_Len <= 0)
              {
               fprintf (stderr, "ERROR:  Bad minimum gene length (-g option)\n"
                    "  value = \"%s\"", optarg);
               errflg = true;
              }
          break;

       case  'h' :
          Command_Line . append (" -h");
          errflg = true;
          break;

       case  'i' :
	  Command_Line . append (" -i");
	  Allow_Indels = true;
	  break;
	 
       case  'm' :
	    Command_Line . append (" -m ");
	    Command_Line . append (optarg);
	    ICM_File_Name = optarg;
	    User_ICM = true;
	    break;

       case  'o' :
          Command_Line . append (" -o ");
          Command_Line . append (optarg);
          Max_Olap_Bases = strtol (optarg, & p, 10);
          if  (p == optarg || Max_Olap_Bases < 0)
              {
               fprintf (stderr, "ERROR:  Bad max overlap bases (-o option)\n"
                    "  value = \"%s\"", optarg);
               errflg = true;
              }
          break;

       case 'q' :
	    Command_Line . append (" -Q ");
	    Command_Line . append (optarg);
	    Allow_Indels = true;
	    Quality_File_Name = optarg;
	    break;

       case 'r' :
	    Command_Line . append (" -r ");
	    Genome_Is_Circular = true;
	    Allow_Truncated_Orfs = false;
	    break;

       case 's' :
	    Command_Line . append (" -S");
	    Allow_Subs = true;
	    break;

       case  'u' :
	    Command_Line . append (" -u ");
	    Command_Line . append (optarg);
	    LogOdds_Fudge = strtod (optarg, & p);
	    if  (p == optarg) {
		 fprintf (stderr, "ERROR:  Bad value for fudge factor (-u option)\n"
			  "  value = \"%s\"", optarg);
		 errflg = true;
	    }
	    LogOdds_Prior += LogOdds_Fudge;
	    break;

        case  'z' :
          Command_Line . append (" -z ");
          Command_Line . append (optarg);
	  User_Stop = true;
          Genbank_Xlate_Code = strtol (optarg, & p, 10);
          Set_Stop_Codons_By_Code (Stop_Codon, Genbank_Xlate_Code, errflg);
          break;

        case  'Z' :
          Command_Line . append (" -Z ");
          Command_Line . append (optarg);
	  User_Stop = true;
          Stop_Codon . clear ();
          for  (p = strtok (optarg, ",");  p != NULL;  p = strtok (NULL, ","))
            {
             q = strdup (p);
             Make_Lower_Case (q);
             Stop_Codon . push_back (q);
            }
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = true;
       }

   if  (errflg)
       {
        Usage ();
        exit (EXIT_FAILURE);
       }

   for  (i = optind;  i < argc;  i ++)
   {
	Command_Line . append (" ");
	Command_Line . append (argv [i]);
   }


   if  (optind > argc - 2)
   {
	Usage ();
	exit (EXIT_FAILURE);
   }	
   Sequence_File_Name = argv [optind ++];
   Output_Tag = argv [optind ++];

   ////////////////////////////////////////////
   // check for bad input
   ////////////////////////////////////////////
   if(Allow_Indels && Allow_Subs) {
	sprintf (Clean_Exit_Msg_Line, "ERROR:  Cannot use --indel and --sub simultaneously\n");
	Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
   }

   return;
  }


static double Pass_Stop_Penalty(int frame, int lo, int hi)

// Compute the log likelihood ratio penalty for passing 
// through this ORF's previous stop codon.

{
     double default_p = 0.999;
     double codon_p[3] = {default_p, default_p, default_p};
     int stop_i[3] = {lo-3, lo-2, lo-1};
     if(frame < 0) {
	  stop_i[0] = hi+1;
	  stop_i[1] = hi;
	  stop_i[2] = hi-1;
     }

     // use quality values if available
     if(Quality_File_Name != NULL) {
	  codon_p[0] = 1.0 - pow(10.0, -(double)Quality_Values[stop_i[0]]/10.0);
	  codon_p[1] = 1.0 - pow(10.0, -(double)Quality_Values[stop_i[1]]/10.0);
	  codon_p[2] = 1.0 - pow(10.0, -(double)Quality_Values[stop_i[2]]/10.0);
     }
     
     // compute probability stop codon is correct trying to account for different stops
     double p_stop = codon_p[0];
     if ((frame > 0 && Sequence[stop_i[1]] == 'a') || (frame < 0 && Sequence[stop_i[1]] == 't'))
	  p_stop *= 2.0/3.0*codon_p[1] + 1.0/3.0;
     else
	  p_stop *= codon_p[1];
     if ((frame > 0 && Sequence[stop_i[2]] == 'a') || (frame < 0 && Sequence[stop_i[2]] == 't'))
	  p_stop *= 2.0/3.0*codon_p[2] + 1.0/3.0;
     else
	  p_stop *= codon_p[2];

     return log(1.0-p_stop) - log(p_stop);
}


static void Read_Meta_ICMs()

// Using classifications, determine all ICMs that we will want to score with
// and make a mapping between ICMs and the reads they will apply to

{
     string seq_header, icm_file;
     vector<string> seq_classes;
     icm_set_hash::iterator isi;

     class_hash::const_iterator ci;
     for(ci = classifications.begin(); ci != classifications.end(); ci++) {
	  seq_header = ci->first;
	  seq_classes = ci->second;
	  
	  // determine ICM file name
	  icm_file = Classes_ICM_File(seq_classes);

	  // check hash table
	  isi = ICM_Sequences.find(icm_file);
	  if(isi == ICM_Sequences.end()) {
	       // if unfound, make a new set and add
	       seq_set icm_seqs;
	       icm_seqs.insert(seq_header);
	       ICM_Sequences[icm_file] = icm_seqs;
	  } else {
	       // if found, just add
	       isi->second.insert(seq_header);
	  }
     }
}


static void Read_Meta_RBS()

// Using classifications, determine all RBS PWM's that we will want to score
// with (which is currently only the top classification for each sequence) and
// save in a hash_map.  Incorporate null model via classification GC%.

{
     string train_str, rbs_file, gc_file, line;
     vector<string> seq_classes, strain_nc;
     //float class_gc;
     rbs_hash::const_iterator pi;
     FILE  * fp;

     class_hash::const_iterator ci;
     for(ci = classifications.begin(); ci != classifications.end(); ci++)
     {
	  seq_classes = ci->second;
	  for(unsigned int i = 0; i < seq_classes.size(); i++) {
	       train_str = seq_classes[i];

	       pi = Sequence_RBSs.find(train_str);
	       if(pi == Sequence_RBSs.end()) {
		    strain_nc = split(train_str, '|');
		    rbs_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".motif";

		    // get RBS PWM
		    PWM_t * class_rbs = new PWM_t();
		    fp = File_Open (rbs_file.c_str(), "r", __FILE__, __LINE__);
		    class_rbs->Read(fp);
		    fclose(fp);
		    class_rbs->Counts_To_Prob();
		    /*
		      class_rbs->Probs_To_Logs();
	       
		      // get GC%
		      gc_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".gc.txt";
		      ifstream gc_open(gc_file.c_str());
		      if(gc_open.good()) {
		      getline(gc_open, line);
		      class_gc = strtod(line.c_str(), NULL);
		      } else {
		      cerr << "WARNING: GC classification file unavailable " << gc_file << endl;
		      class_gc = 0.5;
		      }
		      gc_open.close();
	       
		      // make log odds
		      class_rbs->Make_Log_Odds_WRT_GC(class_gc);
		    */

		    // save
		    Sequence_RBSs[train_str] = class_rbs;
	       }
	  }
     }
}


static void Read_Meta_Lengths()

// Using classifications, determine all length distributions that we will
// want to score with and save the counts in a hash_map.

{
     string train_str, length_file;
     vector<string> seq_classes, strain_nc;
     length_hash::const_iterator li;

     class_hash::const_iterator ci;
     for(ci = classifications.begin(); ci != classifications.end(); ci++) {
	  seq_classes = ci->second;
	  for(unsigned int i = 0; i < seq_classes.size(); i++) {
	       train_str = seq_classes[i];

	       li = Sequence_Lengths_Gene.find(train_str);
	       if(li == Sequence_Lengths_Gene.end()) {
		    strain_nc = split(train_str, '|');

		    // gene length dist
		    length_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".lengths.genes.txt";
		    //length_file = "/fs/szasmg/dakelley/research/gene_prediction/results/2-23/all.lengths.genes.txt";
		    ifstream length_gene_in(length_file.c_str());
		    vector<double> lengths_gene;
		    float gene_count = 0;
		    if(length_gene_in.good()) {
			 //cout << length_file << endl;
			 gene_count = Read_Length_Dist(length_gene_in, lengths_gene);
			 length_gene_in.close();
		    } else
			 cerr << "ERROR:  Cannot open gene length file " << length_file << endl;
		    Sequence_Lengths_Gene[train_str] = lengths_gene;

		    // non length dist
		    length_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".lengths.non.txt";
		    ifstream length_non_in(length_file.c_str());
		    vector<double> lengths_non;
		    float nonorf_count = 0;
		    if(length_non_in.good()) {
			 nonorf_count = Read_Length_Dist(length_non_in, lengths_non);
			 length_non_in.close();
		    } else
			 cerr << "ERROR:  Cannot open noncoding ORF length file " << length_file << endl;
		    Sequence_Lengths_Non[train_str] = lengths_non;

		    // set prior probability
		    if(gene_count > 0  && nonorf_count > 0)
			 Sequence_Prior[train_str] = log(gene_count/nonorf_count);
		    else
			 Sequence_Prior[train_str] = 0;

		    // now re-write gene lengths using static dist
		    /*
		    length_file = "/fs/szasmg/dakelley/research/gene_prediction/results/2-23/all.lengths.genes.txt";
		    ifstream length_gene_in2(length_file.c_str());
		    vector<double> lengths_gene2;
		    if(length_gene_in2.good()) {
			 Read_Length_Dist(length_gene_in2, lengths_gene2);
			 length_gene_in2.close();
		    } else
			 cerr << "ERROR:  Cannot open gene length file " << length_file << endl;
		    Sequence_Lengths_Gene[train_str] = lengths_gene2;
		    */
	       }
	  }
     }
}


static void Read_Meta_Starts()

// Using classifications, determine all start distributions that we will
// want to score with and save the counts in a hash_map.

{
     string train_str, start_file;
     vector<string> seq_classes, strain_nc;
     start_hash::const_iterator li;

     class_hash::const_iterator ci;
     for(ci = classifications.begin(); ci != classifications.end(); ci++) {
	  seq_classes = ci->second;
	  for(unsigned int i = 0; i < seq_classes.size(); i++) {
	       train_str = seq_classes[i];

	       li = Sequence_Starts_Gene.find(train_str);
	       if(li == Sequence_Starts_Gene.end()) {
		    strain_nc = split(train_str, '|');

		    // gene start dist
		    start_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".starts.genes.txt";
		    ifstream start_gene_in(start_file.c_str());
		    vector<float> starts_gene;
		    if(start_gene_in.good()) {
			 Read_Start_Dist(start_gene_in, starts_gene);
			 start_gene_in.close();
		    } else
			 cerr << "WARNING:  Cannot open gene start codon file " << start_file << endl;
		    Sequence_Starts_Gene[train_str] = starts_gene;

		    // non start dist
		    start_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".starts.non.txt";
		    ifstream start_non_in(start_file.c_str());
		    vector<float> starts_non;
		    if(start_non_in.good()) {
			 Read_Start_Dist(start_non_in, starts_non);
			 start_non_in.close();
		    } else
			 cerr << "WARNING: Cannot open noncoding ORF start codon file " << start_file << endl;
		    Sequence_Starts_Non[train_str] = starts_non;
	       }
	  }
     }    
}


static void Read_Meta_Stops()

// Using classifications, determine all translation table codes
// which tell us the valid stop codons and save in a hash_map.

{
     string train_str, gbk_file, line;
     size_t tt_i;
     int tt_code;
     vector<string> strain_nc;
     transl_hash::const_iterator li;

     class_hash::const_iterator ci;
     for(ci = classifications.begin(); ci != classifications.end(); ci++) {
	  train_str = ci->second[0];

	  li = Sequence_Transl.find(train_str);
	  if(li == Sequence_Transl.end()) {
	       strain_nc = split(train_str, '|');

	       // find GenBank transl_table code
	       gbk_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".gbk";
	       ifstream gbk_in(gbk_file.c_str());
	       tt_i = string::npos;
	       while(getline(gbk_in, line)) {
		    tt_i = line.find("transl_table=");
		    if(tt_i != string::npos) {
			 tt_code = strtol(line.substr(tt_i+13).c_str(), NULL, 10);
			 Sequence_Transl[train_str] = tt_code;
			 break;
		    }
	       }
	       gbk_in.close();
	       
	       // if not found, use default
	       if(tt_i == string::npos)
		    Sequence_Transl[train_str] = 11;
	  }
     }
}


static void Read_Meta_AdjOr()

// Using classifications, determine all adjacent orientation
// distributions that we will want to score with and save the
// counts in a hash_map.

{
     string train_str, adjor_file;
     vector<string> seq_classes, strain_nc;
     adjor_hash::const_iterator li;

     class_hash::const_iterator ci;
     for(ci = classifications.begin(); ci != classifications.end(); ci++) {
	  seq_classes = ci->second;
	  for(unsigned int i = 0; i < seq_classes.size(); i++) {
	       train_str = seq_classes[i];
	       li = Sequence_AdjOr_Gene.find(train_str);
	       if(li == Sequence_AdjOr_Gene.end()) {
		    strain_nc = split(train_str, '|');

		    // gene adjor dist
		    adjor_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".adj_orients.genes.txt";
		    ifstream adjor_gene_in(adjor_file.c_str());
		    vector<float> adjor_gene;
		    if(adjor_gene_in.good()) {
			 Read_Orient_Dist(adjor_gene_in, adjor_gene);
			 adjor_gene_in.close();
		    } else
			 cerr << "WARNING: Cannot open gene adjacent orientation file " << adjor_file << endl;
		    Sequence_AdjOr_Gene[train_str] = adjor_gene;

		    // non adjor dist
		    adjor_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".adj_orients.non.txt";
		    ifstream adjor_non_in(adjor_file.c_str());
		    vector<float> adjor_non;
		    if(adjor_non_in.good()) {
			 Read_Orient_Dist(adjor_non_in, adjor_non);
			 adjor_non_in.close();
		    } else
			 cerr << "WARNING: Cannot open noncoding ORF adjacent orientation file " << adjor_file << endl;
		    Sequence_AdjOr_Non[train_str] = adjor_non;
	       }
	  } 
     }
}


static void Read_Meta_AdjDist()

// Using classifications, determine all adjacent distance
// distributions that we will want to score with and save
// the counts in a hash_map.  Also, process and set the
// distribution max overlap.

{
     string train_str, adjdist_file;
     vector<string> seq_classes, strain_nc;
     adjdist_hash::const_iterator di;

     class_hash::const_iterator ci;
     for(ci = classifications.begin(); ci != classifications.end(); ci++) {
	  seq_classes = ci->second;
	  for(unsigned int i = 0; i < seq_classes.size(); i++) {
	       train_str = seq_classes[i];
	       di = Sequence_AdjDist_ff_Gene.find(train_str);
	       if(di == Sequence_AdjDist_ff_Gene.end()) {
		    strain_nc = split(train_str, '|');

		    // gene adj dist 1,1 dist
		    adjdist_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".adj_dist.1.1.genes.txt";
		    ifstream adjdist_ff_gene_in(adjdist_file.c_str());
		    vector<float> adjdist_ff_gene;
		    if(adjdist_ff_gene_in.good()) {
			 Read_Dist_Dist(adjdist_ff_gene_in, adjdist_ff_gene);
			 adjdist_ff_gene_in.close();
		    } else
			 cerr << "WARNING:  Cannot open gene adjacent distance file " << adjdist_file << endl;
		    Sequence_AdjDist_ff_Gene[train_str] = adjdist_ff_gene;
		    
		    // non adj dist 1,1 dist
		    adjdist_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".adj_dist.1.1.non.txt";
		    ifstream adjdist_ff_non_in(adjdist_file.c_str());
		    vector<float> adjdist_ff_non;
		    if(adjdist_ff_non_in.good()) {
			 Read_Dist_Dist(adjdist_ff_non_in, adjdist_ff_non);
			 adjdist_ff_non_in.close();
		    } else
			 cerr << "WARNING: Cannot open noncoding ORF adjacent distance file " << adjdist_file << endl;
		    Sequence_AdjDist_ff_Non[train_str] = adjdist_ff_non;

		    // gene adj dist 1,-1 dist
		    adjdist_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".adj_dist.1.-1.genes.txt";
		    ifstream adjdist_fr_gene_in(adjdist_file.c_str());
		    vector<float> adjdist_fr_gene;
		    if(adjdist_fr_gene_in.good()) {
			 Read_Dist_Dist(adjdist_fr_gene_in, adjdist_fr_gene);
			 adjdist_fr_gene_in.close();
		    } else
			 cerr << "WARNING: Cannot open gene adjacent distance file " << adjdist_file << endl;
		    Sequence_AdjDist_fr_Gene[train_str] = adjdist_fr_gene;

		    // non adj dist 1,-1 dist
		    adjdist_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".adj_dist.1.-1.non.txt";
		    ifstream adjdist_fr_non_in(adjdist_file.c_str());
		    vector<float> adjdist_fr_non;
		    if(adjdist_fr_non_in.good()) {
			 Read_Dist_Dist(adjdist_fr_non_in, adjdist_fr_non);
			 adjdist_fr_non_in.close();
		    } else
			 cerr << "WARNING: Cannot open noncoding ORF adjacent distance file " << adjdist_file << endl;
		    Sequence_AdjDist_fr_Non[train_str] = adjdist_fr_non;

		    // gene adj dist -1,1 dist
		    adjdist_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".adj_dist.-1.1.genes.txt";
		    vector<float> adjdist_rf_gene;
		    ifstream adjdist_rf_gene_in(adjdist_file.c_str());
		    if(adjdist_rf_gene_in.good()) {
			 Read_Dist_Dist(adjdist_rf_gene_in, adjdist_rf_gene);
			 adjdist_rf_gene_in.close();
		    } else
			 cerr << "WARNING: Cannot open gene adjacent distance file" << adjdist_file << endl;
		    Sequence_AdjDist_rf_Gene[train_str] = adjdist_rf_gene;

		    // non adj dist -1,1 dist
		    adjdist_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".adj_dist.-1.1.non.txt";
		    ifstream adjdist_rf_non_in(adjdist_file.c_str());
		    vector<float> adjdist_rf_non;
		    if(adjdist_rf_non_in.good()) {
			 Read_Dist_Dist(adjdist_rf_non_in, adjdist_rf_non);
			 adjdist_rf_non_in.close();
		    } else
			 cerr << "WARNING: Cannot open noncoding ORF adjacent distance file " << adjdist_file << endl;
		    Sequence_AdjDist_rf_Non[train_str] = adjdist_rf_non;
	       }
	  }
     }

     LogOdds_AdjDist.Set_MaxOverlap(Dist_Max_Overlap);
}


static void Read_Meta_GC()

// Using classifications, determine all GC contents that we will
// want to use to build independent models in a hash_map.

{
     string train_str, gc_file, line;
     vector<string> seq_classes, strain_nc;
     gc_hash::const_iterator di;

     class_hash::const_iterator ci;
     for(ci = classifications.begin(); ci != classifications.end(); ci++) {
	  seq_classes = ci->second;
	  for(unsigned int i = 0; i < seq_classes.size(); i++) {
	       train_str = seq_classes[i];
	       di = Sequence_GC.find(train_str);
	       if(di == Sequence_GC.end()) {
		    strain_nc = split(train_str, '|');

		    gc_file = ICM_dir + "/" + strain_nc[0] + "/" + strain_nc[1] + ".gc.txt";
		    ifstream gc_open(gc_file.c_str());
		    if(gc_open.good()) {
			 getline(gc_open, line);
			 Sequence_GC[train_str] = strtod(line.c_str(), NULL);
		    } else {
			 cerr << "WARNING: GC classification file unavailable " << gc_file << endl;
			 Sequence_GC[train_str] = 0.5;
		    }
		    gc_open.close();
	       }
	  }
     }
}


static void Reverse_Transfer_Qual (vector<int> & buff, int start, int len)

//  Copy to string  buff  the substring of Quality_Values  starting at subscript
//   start  and going to the left for a length of  len .

{    
     buff . resize (len);
     for  (int j = 0;  j < len;  j ++, start --)
	  buff [j] = Quality_Values[start];
}


static int Rev_Next_Stop(int end_point)

// Walk forward from end_point to find the next 
// reverse stop codon, or the end of the sequence.

{
     if(end_point >= 0 && end_point < Sequence_Len)
	  return Rev_Next_Stops[end_point];
     else
	  return end_point;
/*
     Codon_t codon;
     int i = end_point;
     int c, which;

     while(i <= Sequence_Len-3) {
	  //get next codon
	  for(c = 0; c < 3; c++) {
	       codon.Shift_In(Sequence[i]);
	       i++;
	  }
	  codon.Reverse_Complement();

	  if(codon.Must_Be(Fwd_Stop_Pattern, which))
	       return i-3;
     }

     return i;
*/
}


static void  Score_All_Frames ()

// Score the entire sequence in all 6 frames, which helps for
// error prediction because we avoid re-scoring multiple times
//
// I have no idea why, but the codons seem to be scored in frame
// order 0, 2, 1. Which is why when cumulative score is called,
// the frame 1 is given.

{
     string buff;
     vector<double> gene_scores, non_scores;
     
     // reverse sequence
     Reverse_Transfer(buff, Sequence, Sequence_Len-1, Sequence_Len);

     // for each frame
     for(int f = 0; f < 3; f++) {
	  // score sequence
	  Gene_ICM.Frame_Score(buff, gene_scores, f);
	  Indep_Model.Frame_Score(buff, non_scores, f);

	  // un-reverse and combine
	  Frame_Scores[f].resize(Sequence_Len);
	  for(int i = 0; i < Sequence_Len; i++)
	       Frame_Scores[f][i] = gene_scores[Sequence_Len-1-i] - non_scores[Sequence_Len-1-i];
     }

     // reverse sequence
     Complement_Transfer (buff, Sequence, 0, Sequence_Len);
     
     // for each frame
     for(int f = 0; f < 3; f++) {
	  // score sequence
	  Gene_ICM.Frame_Score(buff, gene_scores, f);
	  Indep_Model.Frame_Score(buff, non_scores, f);

	  // un-reverse and combine
	  Frame_Scores[3+f].resize(Sequence_Len);
	  for(int i = 0; i < Sequence_Len; i++)
	       Frame_Scores[3+f][i] = gene_scores[i] - non_scores[i];
     }
}


static void Score_Indels(Orf_t & orf, vector<Start_t> & start_list, vector<Error_t> & errors, double suffix_score, int suffix_j, vector<double> & score, vector<double> & indep_score, int q, int k, int j)

// Branch off into different frames to implicitly predict an insertion and deletion.

{
     int error_end_point;
     int error_suffix_j;
     double error_suffix_score;
     int frame = orf.Get_Frame();
     double prob_err = pow(10.0, -(double)q/10.0);
     double score_penalty = log(prob_err/2.0) - log(1.0-prob_err);

     if(frame > 0)
     {
	  error_suffix_score = suffix_score + score[j] - indep_score[j] + score_penalty;

	  if(error_suffix_score > Indel_Suffix_Score_Threshold)
	  {
	       error_end_point = k + (j % 3); // move to last bp of next valid codon (1-based)
	       error_suffix_j = suffix_j + j + 2 - (j % 3);

	       vector<Error_t> del_indels(errors);
	       Error_t del(k+3, 1); // 1-based coordinate of deleted nt original or following nt in fragment
	       del_indels.push_back(del);
			 
	       if (Dave_Log)
		    cout << Fasta_Header << "\td\t" << (frame>0) << "\t" <<
			 orf.Get_Stop_Position() << "\t" << (k+3) << "\t" << (j%3) << "\n";
			 
	       Score_Orf_Starts(orf, start_list, error_end_point, error_suffix_score, error_suffix_j, del_indels);
	  }


	  error_suffix_score = suffix_score + score[j-1] - indep_score[j-1] + score_penalty;
		    
	  if(error_suffix_score > Indel_Suffix_Score_Threshold)
	  {
	       error_end_point = k - (2 - (j % 3)); // move back one codon
	       error_suffix_j = suffix_j + j + 2 - (j % 3);
			 
	       vector<Error_t> ins_indels(errors);
	       Error_t ins(k+2, 0);  // 1-based coordinate of nt following insertion in original or inserted nt in fragment
	       ins_indels.push_back(ins);
			 
	       if(Dave_Log)
		    cout << Fasta_Header << "\ti\t" << (frame>0) << "\t" <<
			 orf.Get_Stop_Position() << "\t" << (k+2) << "\t" << (j%3) << "\n";

	       Score_Orf_Starts(orf, start_list, error_end_point, error_suffix_score, error_suffix_j, ins_indels);
	  }
     }
     else
     {
	  error_suffix_score = suffix_score + score[j] - indep_score[j] + score_penalty;

	  if(error_suffix_score > Indel_Suffix_Score_Threshold)
	  {
	       error_end_point = k - (j % 3);
	       error_suffix_j = suffix_j + j + 2 - (j % 3);

	       vector<Error_t> del_indels(errors);
	       Error_t del(k-1, 1); // 1-based coordinate of deleted nt original or following nt in fragment
	       del_indels.push_back(del);

	       if(Dave_Log)
		    cout << Fasta_Header << "\td\t" << (frame>0) << "\t" <<
			 orf.Get_Stop_Position() << "\t" << (k-1) << "\t" << (j%3) << "\n";
			 
	       Score_Orf_Starts(orf, start_list, error_end_point, error_suffix_score, error_suffix_j, del_indels);
	  }

	  error_suffix_score = suffix_score + score[j-1] - indep_score[j-1] + score_penalty;

	  if(error_suffix_score > Indel_Suffix_Score_Threshold)
	  {
	       error_end_point = k + 2 - (j % 3);
	       error_suffix_j = suffix_j + j + 2 - (j % 3);

	       vector<Error_t> ins_indels(errors);
	       Error_t ins(k-2, 0); // 1-based coordinate of nt following insertion in original or inserted nt in fragment
	       ins_indels.push_back(ins);

	       if(Dave_Log)
		    cout << Fasta_Header << "\ti\t" << (frame>0) << "\t" <<
			 orf.Get_Stop_Position() << "\t" << (k-2) << "\t" << (j%3) << "\n";

	       Score_Orf_Starts(orf, start_list, error_end_point, error_suffix_score, error_suffix_j, ins_indels);
	  }
     }
}


static void  Score_Orfs_Errors
    (vector <Orf_t> & orf_list, FILE * fp)

//  Compute scores for all orfs in  orf_list  using coding model
//  in global  Gene_ICM , which is assumed to have been built on reverse
//  gene strings.   Indep_Model  is the model of independent,
//  stop-codon-free sequence.  Put orfs that are candidate genes
//  onto  gene_list .  Print log information to  fp .

{
     vector <Start_t>  start_list;
     int  i, n, id = 0;

     // sequence pre-processing
     Score_All_Frames();
     Save_Prev_Stops();

     // reset saved PWM scores
     Meta_PWM_Save.resize(2*Sequence_Len);
     for(unsigned int si = 0; si < 2*Sequence_Len; si++) {
	  pair<double,int> dummy(0.0, 999);
	  Meta_PWM_Save[si] = dummy;
     }
     
     n = orf_list . size ();
     for  (i = 0;  i < n;  i ++)
     {
	  int first_j = 0;
	  int end_point;
	  int frame = orf_list[i].Get_Frame();
	  vector<Error_t> errors;

	  ////////////////////////////////////////////
	  // score starts
	  ////////////////////////////////////////////
	  start_list . clear ();
	  if(frame > 0)
	       end_point = orf_list[i].Get_Stop_Position() - 1;
	  else
	       end_point = orf_list[i].Get_Stop_Position() + 3;

	  Score_Orf_Starts(orf_list[i], start_list, end_point, 0, 0, errors);

	  // boost long ORFs
	  for(unsigned int start_i = 0; start_i < start_list.size(); start_i++)
	       if(start_list[start_i].j > Ignore_Score_Len)
		    start_list[start_i].score = Max(0.0, start_list[start_i].score);

	  ////////////////////////////////////////////
	  // filter 
	  ////////////////////////////////////////////
	  if(!start_list.empty())
	  {
	       // note: reverse starts are sorted in reverse order
	       sort(start_list.begin(), start_list.end(), Start_Cmp);
	  
	       if(frame > 0)
		    first_j = start_list.front().j;
	       else
		    first_j = start_list.back().j;

	       if(first_j+1 >= Min_Gene_Len) {

		    double best_score = -DBL_MAX;
		    int best_j = 0;
		    for(unsigned int s = 0; s < start_list.size(); s++) {
			 if(start_list[s].score > best_score) {
			      best_score = start_list[s].score;
			      best_j = start_list[s].j;
			 }
		    }

		    // if ORF has a good score, add it
		    //if(best_j + 1 >= Min_Gene_Len && best_score > Start_Threshold) {
		    if(best_score > Start_Threshold) {
			 if(frame > 0)
			      Add_Events_Fwd(orf_list[i], start_list, id);
			 else
			      Add_Events_Rev(orf_list[i], start_list, id);
		    }
	       }
	  }
     }
     
     return;
}


static void Score_Orf_Starts(Orf_t & orf, vector<Start_t> & start_list, int end_point, double suffix_score, int suffix_j, vector<Error_t> & errors)
{
     string seq_buff;
     vector<int> qual_buff;
     Start_t  start;
     vector<double> score;
     vector<double> indep_score;
     int first_pos = 0;
     Codon_t  codon;
     bool  orf_is_truncated = false;
     int  which;
     int  frame;
     int  lo, hi, len, lowest_j;
     int  j, k, m;
     int error_end_point;
     int error_suffix_j;
     double error_suffix_score;
     double  next_s;
     double best_score = -DBL_MAX;
     int num_errors = (int)errors.size();
     
     ////////////////////////////////////////////
     // score sequence
     ////////////////////////////////////////////
     frame = orf . Get_Frame ();
     len = orf . Get_Orf_Len ();

     if(frame > 0) {
	  // end_point is the last base of the gene (1-based)
	  // orf.Get_Stop_Position is the first base of the stop codon (1-based)
	  // hi is the last base of the gene (1-based)
	  //   thus, Reverse transfer substracts 1 to make it 0-based
	  // len is the number of nucleotides between stop codons
	  // lo is the last nt of the previous stop codon (1-based)

	  hi = end_point;
	  lo = Fwd_Prev_Stop(end_point-1)+1;
	  len = hi-lo;
	  if(len >= 0)
	  {
	       Reverse_Transfer (seq_buff, Sequence, hi-1, len);
	       if(Allow_Indels)
		    Reverse_Transfer_Qual (qual_buff, hi-1, len);
	  }
	  orf_is_truncated = (lo < 3 && Allow_Truncated_Orfs);
	  k = lo-1;

	  //cout << "ORF: " << "fwd " << Fasta_Header << " " << orf.Get_Frame() << " " << orf.Get_Stop_Position() << " " << orf.Get_Orf_Len() << " " << hi << " " << lo << " " << len << endl;
     }
     else
     {
	  // end_point is the first base of the gene (1-based)
	  // lo is the first base of the gene (1-based)
	  // hi is the first base of the next stop codon (1-based)
	  // len is the number of nucleotides between stop codons
	  lo = end_point;
	  hi = Rev_Next_Stop(end_point-1)+1;
	  len = hi-lo;
	  if(lo-1 < Sequence_Len)
	  {
	       Complement_Transfer (seq_buff, Sequence, lo-1, len);
	       if(Allow_Indels)
		    Complement_Transfer_Qual (qual_buff, lo-1, len);
	  }
	  orf_is_truncated = (Sequence_Len - (hi-1) < 3 && Allow_Truncated_Orfs);
	  k = hi+1;
	  
	  //cout << "ORF: " << "rev " << Fasta_Header << " " << orf.Get_Frame() << " " << orf.Get_Stop_Position() << " " << orf.Get_Orf_Len() << " " << hi << " " << lo << " " << len << endl;
     }
     
     //Gene_ICM . Cumulative_Score (seq_buff, score, 1);
     //Indep_Model . Cumulative_Score (seq_buff, indep_score, 1);
     Cumulative_Frame_Score (frame, lo, hi, score, indep_score);

     
     ////////////////////////////////////////////
     // mutate previous codon
     ////////////////////////////////////////////
     if(Allow_Subs && num_errors < 1) {
	  int error_pos;

	  // set end_point and error point
	  if(frame > 0) {
	       error_end_point = lo - 3;
	       error_pos = lo - 3;
	  } else {
	       error_end_point = hi + 3;
	       error_pos = hi + 1;
	  }

	  // BE VERY CAREFUL HERE BECAUSE IT WILL AFFECT WHETHER WE MUTATE TO GET TO THE END OF THE GENE OR NOT
	  if(error_end_point >= 0 && error_end_point-2 < Sequence_Len) {
	       // set j
	       error_suffix_j = suffix_j + len;

	       // add in error llr
	       error_suffix_score = suffix_score + Pass_Stop_Penalty(frame, lo, hi);
	       if(!score.empty())
		    error_suffix_score += score.back() - indep_score.back();

	       // add new error
	       vector<Error_t> errors_sub(errors);
	       Error_t stop_sub(error_pos, 2);
	       errors_sub.push_back(stop_sub);

	       if(Dave_Log)
		    cout << Fasta_Header << "\t" << frame << "\t" << orf.Get_Stop_Position() << "\t" <<
			 lo << "\t" << hi << "\t" << error_suffix_score << "\n";

	       Score_Orf_Starts(orf, start_list, error_end_point, error_suffix_score, error_suffix_j, errors_sub);
	  }
     }


     ////////////////////////////////////////////
     // find starts
     ////////////////////////////////////////////
     m = score . size ();
     lowest_j = Min (3, Min_Gene_Len - 3);
     for  (j = m - 1;  j >= lowest_j;  j --)
     {
	  if(Allow_Indels && qual_buff[j] <= Indel_Quality_Threshold && num_errors < Indel_Max)
	       Score_Indels(orf, start_list, errors, suffix_score, suffix_j, score, indep_score, qual_buff[j], k, j);

	  codon . Shift_In (seq_buff [j]);

	  if  (j % 3 == 0
	       && (codon . Can_Be (Fwd_Start_Pattern, which)
		   || (first_pos == 0 && orf_is_truncated))
	       && j + 3 + suffix_j >= Min_Gene_Len)
	  {    
	       next_s = score [j - 1] - indep_score [j - 1];
	       // this is the score for the orf without the start
	       // codon--position j is the last base of the start codon
	       start . j = j + 2 + suffix_j;
	       start . pos = k;
	       // k is the 1-based sequence coordinate of the base that
	       // is 2 behind the position represented by j
	       start . score = next_s + suffix_score;
	       if(start.score > best_score)
		    best_score = start.score;
	       start . first = (first_pos == 0);
	       start . errors = errors;
	       
	       if(which >= 0 && first_pos == 0 && orf_is_truncated) {
		    // if first codon is a start, we still need
		    // to consider the truncated gene first
		    start.which = -1;
		    start.truncated = true;
		    start_list.push_back(start);
		    
		    start.first = false; // for the next guy
	       }
		    
	       start . which = which;
	       start . truncated = (which < 0);
	       start_list . push_back (start);
	       
	       if  (first_pos == 0) {
		    first_pos = k;
	       }	       
	  }
	  if  (frame > 0)
	       k ++;
	  else
	       k --;
     }
}


static void Set_Quality_454()

// Set the quality values for a 454 sequence based on
// the homopolymer runs

{
     unsigned int i,q;
     vector<int> Run_Qualities;
     for(q = 0; q < 6; q++)
	  Run_Qualities.push_back(31 - 5*q);

     unsigned int homopolymer_run = 0;
     char last_nt = ' ';

     Quality_Values.resize(Sequence.size());
     for(i = 0; i < Sequence.size(); i++) {
	  if(Sequence[i] != last_nt) {
	       if(i > 0) {
		    // set last nt of run to lower quality
		    if(homopolymer_run < Run_Qualities.size())
			 Quality_Values[i-1] = Run_Qualities[homopolymer_run];
		    else
			 Quality_Values[i-1] = Run_Qualities.back();
	       }

	       homopolymer_run = 1;
	  } else {
	       // set middling nt's to high quality
	       Quality_Values[i-1] = 31;

	       homopolymer_run++;
	  }
	  
	  last_nt = Sequence[i];
     }

     // set last nt of run to lower quality
     if(homopolymer_run < Run_Qualities.size())		    
	  Quality_Values[i-1] = Run_Qualities[homopolymer_run];
     else
	  Quality_Values[i-1] = Run_Qualities.back();
}


static void  Trace_Back
    (FILE * fp, const Event_Node_t & final_event)

//  Trace back through the list of best events starting at
//   final_event . best_pred  and output to  fp  the corresponding
//  set of genes.

{
   Event_Node_t  * p;
   vector <Gene_t>  gene_list;
   Gene_t  gene;
   double  prev_score = 0; // to avoid warning
   int rev_start = 0; // to avoid warning
   int  f, i, j, n;
   vector<Error_t> rev_errors;

   for  (p = final_event . best_pred;  p -> e_type != INITIAL;  p = p -> best_pred)
   {
	switch  (p -> e_type)
        {
	case  FWD_START :
	     j = gene . Get_Stop_Position ();
	     gene . Set_Gene_Len (2 + j - p -> pos);
	     gene . Set_Score (p -> score - p -> best_pred -> score);
	     gene . Set_ID (p -> id);
	     gene . Set_Errors (p -> errors);
	     if  (p -> truncated)
		  gene . Set_Status_Bit (TRUNCATED_START_FLAG);
	     gene_list . push_back (gene);
	     gene . Clear_Status ();
	     break;
	case  FWD_STOP :
	     gene . Set_Stop_Position (p -> pos - 2);
	     gene . Set_Frame (1 + (p -> pos % 3));
	     break;
	case  REV_START :
	     rev_start = p -> pos;
	     prev_score = p -> score;
	     rev_errors = p -> errors;
	     if  (p -> truncated)
		  gene . Set_Status_Bit (TRUNCATED_START_FLAG);
	     break;
	case  REV_STOP :
	     gene . Set_Stop_Position (p -> pos - 2);
	     gene . Set_Frame (- (1 + (p -> pos % 3)));
	     gene . Set_Gene_Len (rev_start - p -> pos);
	     gene . Set_Score (prev_score - p -> score);
	     gene . Set_ID (p -> id);
	     gene . Set_Errors (rev_errors);
	     gene_list . push_back (gene);
	     gene . Clear_Status ();
	     break;
	default :
	     printf ("Bad event type = %d\n", int (p -> e_type));
	     exit (EXIT_FAILURE);
        }
   }

   n = gene_list . size ();

   // Adjust stop positions to be in the range  1 .. Sequence_Len
   // and set the frame accordingly
   for  (i = 0;  i < n;  i ++) {
	j = gene_list [i] . Get_Stop_Position ();
	f = Position_To_Frame (j);
	if  (gene_list [i] . Get_Frame () > 0)
	     gene_list [i] . Set_Frame (f);
        else
	     gene_list [i] . Set_Frame (-1 * f);
   }

   //sort (gene_list . begin (), gene_list . end (), By_ID);

   for  (i = n-1;  i >= 0;  i --) {
	int  start, stop;
	
	if  (gene_list [i] . Get_Frame () > 0) {
	     stop = gene_list [i] . Get_Stop_Position () + 2;
	     start = stop - gene_list [i] . Get_Gene_Len () - 2;
	     if  (gene_list [i] . Get_Status_Bit (TRUNCATED_START_FLAG))
		  start -= 3;
	     // move an artificial start at the beginning of the sequence
	     // off the front to indicate the gene could extend there
	} else {
	     stop = gene_list [i] . Get_Stop_Position ();
	     start = stop + gene_list [i] . Get_Gene_Len () + 2;
	     if  (gene_list [i] . Get_Status_Bit (TRUNCATED_START_FLAG))
		  start += 3;
	     // move an artificial start at the end of the sequence
	     // off the back to indicate the gene could extend there
	}
	
	// separate insertions and deletions
	vector<Error_t> gene_errors = gene_list[i].Get_Errors();
	vector<int> insertions;
	vector<int> deletions;
	vector<int> substitutions;
	for(unsigned int errors_i = 0; errors_i < gene_errors.size(); errors_i++) {
	     if(gene_errors[errors_i].type == 0)
		  insertions.push_back(gene_errors[errors_i].pos);
	     else if(gene_errors[errors_i].type == 1)
		  deletions.push_back(gene_errors[errors_i].pos);
	     else
		  substitutions.push_back(gene_errors[errors_i].pos);
	}
	sort(insertions.begin(), insertions.end());
	sort(deletions.begin(), deletions.end());
	sort(substitutions.begin(), substitutions.end());
	
	// print gene
	fprintf (fp, "orf%05d %8d %8d %+3d %8.2f",
		 gene_list [i] . Get_ID (),  start, stop,
		 gene_list [i] . Get_Frame (),
		 gene_list [i] . Get_Score ());
	
	// print errors
	fprintf(fp, " I:");
	if(!insertions.empty()) {
	     fprintf(fp, "%d", insertions[0]);
	     for(unsigned int ins_i = 1; ins_i < insertions.size(); ins_i++)
		  fprintf(fp, ",%d", insertions[ins_i]);
	}
	fprintf(fp, " D:");
	if(!deletions.empty()) {
	     fprintf(fp, "%d", deletions[0]);
	     for(unsigned int del_i = 1; del_i < deletions.size(); del_i++)
		  fprintf(fp, ",%d", deletions[del_i]);
	}
	fprintf(fp, " S:");
	if(!substitutions.empty()) {
	     fprintf(fp, "%d", substitutions[0]);
	     for(unsigned int sub_i = 1; sub_i < substitutions.size(); sub_i++)
		  fprintf(fp, ",%d", substitutions[sub_i]);
	}
	fprintf(fp, "\n");
   }
   
   return;
}


static void Update_Meta_Null_ICM()

// Update null ICM for the next sequence's classifications

{   
     // get header hash
     string header_prefix = split(string(Fasta_Header))[0];
     vector<string> Seq_Classes = classifications[header_prefix];     

     // null ICM
     float num_classes = (float)Seq_Classes.size();
     Indep_GC_Frac = 0.0;
     for(unsigned int s = 0; s < num_classes; s++)
	  Indep_GC_Frac += Sequence_GC[Seq_Classes[s]];
     Indep_GC_Frac /= num_classes;
     //cout << "GC% " << Indep_GC_Frac << endl;
     Indep_Model . Build_Indep_WO_Stops (Indep_GC_Frac, Stop_Codon);
     Set_Ignore_Score_Len ();
}


static void Update_Meta_RBS()

// Update all RBS PWM(s) for the next sequence's classifications

{
     // get header hash
     string header_prefix = split(string(Fasta_Header))[0];
     vector<string> Seq_Classes = classifications[header_prefix];     

     // RBS PWM
     Meta_Ribosome_PWMs.clear();
     //Meta_Ribosome_PWMs.push_back(*Sequence_RBSs[Seq_Classes[0]]);
     for(unsigned int i = 0; i < Seq_Classes.size(); i++)
	  Meta_Ribosome_PWMs.push_back(*Sequence_RBSs[Seq_Classes[i]]);
}


static void Update_Meta_Length()

// Update length model for the next sequence's classifications.
// Remember Sequence_Lengths vectors hold logarithms.

{
     string sclass;
     unsigned int s, l;

     string header_prefix = split(string(Fasta_Header))[0];
     vector<string> Seq_Classes = classifications[header_prefix];

     vector<double> lengths_gene;
     vector<double> lengths_non;
     LogOdds_Prior = LogOdds_Fudge;

     float num_classes = (float)Seq_Classes.size();

     for(s = 0; s < num_classes; s++)
     {
	  sclass = Seq_Classes[s];

	  LogOdds_Prior += Sequence_Prior[sclass] / num_classes;

	  vector<double> lengths_gene_sc = Sequence_Lengths_Gene[sclass];
	  lengths_gene.resize(lengths_gene_sc.size(), log(0));
	  for(l = 0; l < lengths_gene_sc.size(); l++)
	       lengths_gene[l] = log_add(lengths_gene[l], lengths_gene_sc[l]);

	  vector<double> lengths_non_sc = Sequence_Lengths_Non[sclass];
	  lengths_non.resize(lengths_non_sc.size(), log(0));
	  for(l = 0; l < lengths_non_sc.size(); l++)
	       lengths_non[l] = log_add(lengths_non[l], lengths_non_sc[l]);
     }     
     
     for(l = 0; l < lengths_gene.size(); l++)
	  lengths_gene[l] -= log(num_classes);
     for(l = 0; l < lengths_non.size(); l++)
	  lengths_non[l] -= log(num_classes);

     vector<int> seq_length;
     seq_length.push_back(Sequence_Len / 3);

     LogOdds_Length.Make_Log_Odds(lengths_gene, lengths_non, seq_length, (unsigned int)Min_Gene_Len);

     // print models
     if(Dave_Log) {
	  if(strcmp(Fasta_Header,"read0") == 0)
	       LogOdds_Length.Print(Output_Tag);
     }
}


static void Update_Meta_Start()

// Update start codon model for the next sequence's classifications

{
     string sclass;
     unsigned int s, l;

     string header_prefix = split(string(Fasta_Header))[0];
     vector<string> Seq_Classes = classifications[header_prefix];     

     vector<float> starts_gene;
     vector<float> starts_non;

     float num_classes = (float)Seq_Classes.size();

     for(s = 0; s < num_classes; s++)
     {
	  sclass = Seq_Classes[s];

	  vector<float> starts_gene_sc = Sequence_Starts_Gene[sclass];
	  starts_gene.resize(starts_gene_sc.size());
	  for(l = 0; l < starts_gene_sc.size(); l++)
	       starts_gene[l] += starts_gene_sc[l] / num_classes;

	  vector<float> starts_non_sc = Sequence_Starts_Non[sclass];
	  starts_non.resize(starts_non_sc.size());
	  for(l = 0; l < starts_non_sc.size(); l++)
	       starts_non[l] += starts_non_sc[l] / num_classes;
     }
     
     LogOdds_Start.Make_Log_Odds(starts_gene, starts_non);
}


static void Update_Meta_Stop()

// Update stop codon model for the next sequence's classifications

{
     Codon_t  codon;
     bool errflg;
     string header_prefix = split(string(Fasta_Header))[0];
     vector<string> Seq_Classes = classifications[header_prefix];

     // get code
     Genbank_Xlate_Code = Sequence_Transl[Seq_Classes[0]];

     // set stops
     Set_Stop_Codons_By_Code (Stop_Codon, Genbank_Xlate_Code, errflg);

     // set patterns
     Fwd_Stop_Pattern . clear ();
     Rev_Stop_Pattern . clear ();
     Num_Stop_Codons = Stop_Codon . size ();
     for (int i = 0;  i < Num_Stop_Codons;  i ++) {
	  codon . Set_From (Stop_Codon [i]);
	  Fwd_Stop_Pattern . push_back (codon);
	  codon . Reverse_Complement ();
	  Rev_Stop_Pattern . push_back (codon);
     }

     // if not reset in ICM, must reset here
     if(User_ICM) {
	  Indep_Model . Build_Indep_WO_Stops (Indep_GC_Frac, Stop_Codon);
	  Set_Ignore_Score_Len ();
     }

     //cerr << Fasta_Header << " translation table " << Genbank_Xlate_Code << endl;
}


static void Update_Meta_Adj()

// Update adjacency models for the next sequence's classifications

{
     string sclass;
     unsigned int s, l;

     string header_prefix = split(string(Fasta_Header))[0];
     vector<string> Seq_Classes = classifications[header_prefix];     

     vector<float> adjor_gene;
     vector<float> adjor_non;
     vector<float> adjdist_ff_gene;
     vector<float> adjdist_ff_non;
     vector<float> adjdist_fr_gene;
     vector<float> adjdist_fr_non;
     vector<float> adjdist_rf_gene;
     vector<float> adjdist_rf_non;

     float num_classes = (float)Seq_Classes.size();

     for(s = 0; s < num_classes; s++)
     {
	  sclass = Seq_Classes[s];

	  // adj or
	  vector<float> adjor_gene_sc = Sequence_AdjOr_Gene[sclass];
	  adjor_gene.resize(adjor_gene_sc.size());
	  for(l = 0; l < adjor_gene_sc.size(); l++)
	       adjor_gene[l] += adjor_gene_sc[l] / num_classes;

	  vector<float> adjor_non_sc = Sequence_AdjOr_Non[sclass];
	  adjor_non.resize(adjor_non_sc.size());
	  for(l = 0; l < adjor_non_sc.size(); l++)
	       adjor_non[l] += adjor_non_sc[l] / num_classes;

	  // adj dist
	  vector<float> adjdist_ff_gene_sc = Sequence_AdjDist_ff_Gene[sclass];
	  adjdist_ff_gene.resize(adjdist_ff_gene_sc.size());
	  for(l = 0; l < adjdist_ff_gene_sc.size(); l++)
	       adjdist_ff_gene[l] += adjdist_ff_gene_sc[l] / num_classes;

	  vector<float> adjdist_ff_non_sc = Sequence_AdjDist_ff_Non[sclass];
	  adjdist_ff_non.resize(adjdist_ff_non_sc.size());
	  for(l = 0; l < adjdist_ff_non_sc.size(); l++)
	       adjdist_ff_non[l] += adjdist_ff_non_sc[l] / num_classes;

	  vector<float> adjdist_fr_gene_sc = Sequence_AdjDist_fr_Gene[sclass];
	  adjdist_fr_gene.resize(adjdist_fr_gene_sc.size());
	  for(l = 0; l < adjdist_fr_gene_sc.size(); l++)
	       adjdist_fr_gene[l] += adjdist_fr_gene_sc[l] / num_classes;

	  vector<float> adjdist_fr_non_sc = Sequence_AdjDist_fr_Non[sclass];
	  adjdist_fr_non.resize(adjdist_fr_non_sc.size());
	  for(l = 0; l < adjdist_fr_non_sc.size(); l++)
	       adjdist_fr_non[l] += adjdist_fr_non_sc[l] / num_classes;	  

	  vector<float> adjdist_rf_gene_sc = Sequence_AdjDist_rf_Gene[sclass];
	  adjdist_rf_gene.resize(adjdist_rf_gene_sc.size());
	  for(l = 0; l < adjdist_rf_gene_sc.size(); l++)
	       adjdist_rf_gene[l] += adjdist_rf_gene_sc[l] / num_classes;

	  vector<float> adjdist_rf_non_sc = Sequence_AdjDist_rf_Non[sclass];
	  adjdist_rf_non.resize(adjdist_rf_non_sc.size());
	  for(l = 0; l < adjdist_rf_non_sc.size(); l++)
	       adjdist_rf_non[l] += adjdist_rf_non_sc[l] / num_classes;	  
     }     
     
     LogOdds_AdjOr.Make_Log_Odds(adjor_gene, adjor_non);

     LogOdds_AdjDist.Make_Log_Odds_Fwd_Fwd(adjdist_ff_gene, adjdist_ff_non);
     LogOdds_AdjDist.Make_Log_Odds_Fwd_Rev(adjdist_fr_gene, adjdist_fr_non);
     LogOdds_AdjDist.Make_Log_Odds_Rev_Fwd(adjdist_rf_gene, adjdist_rf_non);
}


static void  Usage
    (void)

//  Print to stderr description of options and command line for
//  this program.

  {
   fprintf (stderr,
       "USAGE:  glimmer-mg [options] <sequence-file> <icm-file> <tag>\n"
       "\n"
       "Read DNA sequences in <sequence-file> and predict genes\n"
       "in them using the Interpolated Context Model in <icm-file>.\n"
       "Output details go to file <tag>.detail and predictions go to\n"
       "file <tag>.predict\n"
       "\n"
       "Options:\n"
       " -b <filename>\n"
       " --rbs_pwm <filename>\n"
       "    Read a position weight matrix (PWM) from <filename> to identify\n"
       "    the ribosome binding site to help choose start sites\n"
       " -c <filename>\n"
       " --class <filename>\n"
       "    Read the sequences classifications from <filename> formatted\n"
       "    as \"fasta_header     genome1 genome2 genome3 ...\n"
       " -f <filename>\n"
       " --features <filename>\n"
       "    Read feature counts for a specific organism from <filename>.\n"
       "    See manual for more information.\n"
       " -g <n>\n"
       " --gene_len <n>\n"
       "    Set minimum gene length to <n>\n"
       " -h\n"
       " --help\n"
       "    Print this message\n"
       " -i\n"
       " --indel\n"
       "    Predict genes in \"indel-mode\" where gene predictions may shift\n"
       "    the coding frame, implicitly predicting an insertion or deletion\n"
       "    in the sequence.\n"
       " -m <filename>\n"
       " --icm <filename>\n"
	    "    Read ICM from <filename> and use to score ORF coding likelihood\n"
       " -o <n>\n"
       " --max_olap <n>\n"
       "    Set maximum overlap length to <n>.  Overlaps this short or shorter\n"
       "    are ignored.\n"
       " -q <n>\n"
       " --ignore_score_len <n>\n"
       "    Do not use the initial score filter on any gene <n> or more\n"
       "    base long\n"
       " -r\n"
       " --circular\n"
       "    Assume circular rather than linear genome, i.e., allow wraparound\n"
       " -s\n"
       " --sub\n"
       "    Predict genes in \"substitution-mode\" where gene predictions may\n"
       "    predict a sequencing error in a stop codon and pass through it.\n"
       " -u\n"
       " --fudge\n"
       "    Value to be added to the log-likelihood ratio score of every ORF.\n"
       " -z <n>\n"
       " --trans_table <n>\n"
       "    Use Genbank translation table number <n> for stop codons\n"
       " -Z <codon-list>\n"
       " --stop_codons <codon-list>\n"
       "    Use comma-separated list of codons as stop codons\n"
       "    Sample format:  -Z tag,tga,taa\n"
       "\n");

   return;
  }
