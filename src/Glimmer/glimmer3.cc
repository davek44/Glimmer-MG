//  A. L. Delcher
//
//  File:  glimmer3.cc
//
//  Last Modified:  Tue May  9 10:25:40 EDT 2006
//
//  This program finds open reading frames in the file named
//  on the command line and scores them using the probability
//  model in the file indicated by the second command-line
//  parameter.
//
//  Copyright (c) 2006 University of Maryland Center for Bioinformatics
//  & Computational Biology


#include  "glimmer_base.hh"
#include  "glimmer3.hh"

static int  For_Edwin = 0;

// Global variables

bool  Allow_Truncated_Orfs = false;
  // If set true by -X option, then score orfs that
  // extend to the end of the sequence
Event_Node_t  * Best_Event [6];
  // Best parse event up to the current point in each reading frame
string  Command_Line;
  // Command, options and parameters that invoked the program
vector <double>  Cumulative_Score [6];
  // Prefix-sum score at each position of the input sequence
  // in each reading frame, plus the independent model
  // Frames are, in order:  +1, +2, +3, -1, -2, -3, ind
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
double  Indep_GC_Frac = -1.0;
  // GC proportion used in simple independent model.
  // Set from counts of input sequences or by -C option
int  Ignore_Score_Len = INT_MAX;
  // Genes at least this long do not count the independent model
  // in their score
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
char  * Orflist_File_Name = NULL;
  // Name of file containing a list of regions (which should be valid
  // orfs) that will be scored separately with no overlap rules
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
bool  Separate_Orf_Input = false;
  // If set true by -M option then input is multifasta file
  // of orfs to be scored separately (like Orflist_Option)
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
int  Threshold_Score = DEFAULT_THRESHOLD_SCORE;
  // Minimum score for an orf to be considered a potential gene
bool  Use_Entropy_Profiles = false;
  // If set true (by the -E option) then show the entropy distance
  // ratio in the output.
bool  Use_First_Start_Codon = DEFAULT_USE_FIRST_START_CODON;
  // If true, automatically use the earliest start codon in a gene;
  // otherwise, try to choose the best start codon

bool Allow_Indels = false;
bool Allow_Subs = false;
  // Artifact of Glimmer-MG
int Dist_Max_Overlap = -1;
  // Overlaps of this many or fewer bases are allowed between adjacent
  // genes (as defined by the adjacent distance distributions
double Event_Threshold = -3;
  // Minimum score for an ORF to be added as an event
char * Feature_File = NULL;
  // Name of file used to read in other models
float LogOdds_Prior = DEFAULT_PRIOR;
  // Log odds ratio between genes and noncoding orf counts  
float LogOdds_Fudge = 1.0;
  // Fudge factor to be added to the LogOdds_Prior
Length_Dist_t LogOdds_Length;
  // Log odds ratio between gene length and noncoding orf length
Start_Dist_t LogOdds_Start(DEFAULT_START_PROB);
  // Log odds ratio between gene start codon and noncoding orf start codon
AdjOr_Dist_t LogOdds_AdjOr;
  // Log odds ratio between gene adjacent orientation and noncoding orf adjacent orientation
AdjDist_Dist_t LogOdds_AdjDist;
  // Log odds ratio between gene adjacent distances and noncoding orf adjacent distances
double Start_Threshold = -6;
  // Minimum score for an ORF to be considered by Add_Events
bool User_RBS = false;
bool User_Start = false;
bool User_Length = false;
bool User_Adj = false;


int  main
(int argc, char * argv [])

{
     FILE  * sequence_fp, * detail_fp, * predict_fp;
     vector <string>  seq_list, hdr_list;
     vector <Orf_t>  orf_list;
     vector <Gene_t>  gene_list;
     string  hdr, filename;
     time_t  now;
     int  i;
     string header, header_prefix;
     int seqs_analyzed = 0;     

     try
     {
	  now = time (NULL);
	  cerr << "Starting at " << ctime (& now) << endl;

	  Verbose = 0;

	  Parse_Command_Line (argc, argv);

	  if  (Ignore_File_Name != NULL)
	       Get_Ignore_Regions ();

	  if  (Orflist_File_Name != NULL)
	       Get_Orf_Pos_List ();

	  // outdated start codon code
	  Set_Start_And_Stop_Codons ();
	  Prob_To_Logs (Start_Prob);

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

	  // get GC
	  if  (! GC_Frac_Set)
	       Set_GC_Fraction ();
	  
	  // make null ICM
	  Indep_Model . Build_Indep_WO_Stops (Indep_GC_Frac, Stop_Codon);
	  Set_Ignore_Score_Len ();

	  // get gene ICM
	  Gene_ICM . Read (ICM_File_Name);

	  // prepare RBS PWM(s)
	  LogOdds_PWM = Ribosome_PWM;
	  LogOdds_PWM.Make_Log_Odds_WRT_GC (Indep_GC_Frac);
	  
	  if  (Separate_Orf_Input)
	  {
	       // separate mode that should probably be it's own program

	       sequence_fp = File_Open (Sequence_File_Name, "r", __FILE__, __LINE__);
	       Read_Sequences (sequence_fp, seq_list, hdr_list, Sequence_Ct);
	       fclose (sequence_fp);

	       if(Detail_Log)
		    Print_Orflist_Headings (detail_fp);
	 	
	       for  (i = 0;  i < Sequence_Ct;  i ++)
		    Score_Separate_Input (seq_list [i], hdr_list [i], i, detail_fp, predict_fp);


	  } else if  (Orflist_File_Name != NULL)
	  {
	       // separate mode that should probably be it's own program

	       sequence_fp = File_Open (Sequence_File_Name, "r", __FILE__, __LINE__);
	       Read_Sequences (sequence_fp, seq_list, hdr_list, Sequence_Ct);
	       fclose (sequence_fp);

	       Sequence = seq_list[0];
	       Fasta_Header = hdr_list[0].c_str();

	       if(Detail_Log) {
		    Print_Orflist_Headings (detail_fp);
	       }
	       Score_Orflist (detail_fp, predict_fp);

	  } else {

	       // get sequences
	       sequence_fp = File_Open (Sequence_File_Name, "r", __FILE__, __LINE__);
	       Read_Sequences (sequence_fp, seq_list, hdr_list, Sequence_Ct);
	       fclose (sequence_fp);	       

	       for(i = 0; i < Sequence_Ct; i++)
	       {	    
		    Fasta_Header = hdr_list[i].c_str();
		    header_prefix = split(string(Fasta_Header))[0];

		    // prepare sequence
		    Sequence = seq_list[i];
		    Sequence_Len = Sequence . length ();
		    for(unsigned int seq_i = 0; seq_i < Sequence_Len; seq_i++)
			 Sequence[seq_i] = tolower(Filter(Sequence[seq_i]));

		    if(Detail_Log){
			 fprintf (detail_fp, "\n\n>%s\n", Fasta_Header);
			 Echo_Specific_Settings (detail_fp, Sequence_Len);
		    }
		    fprintf (predict_fp, ">%s\n", Fasta_Header);
			 
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
		    Score_Orfs (orf_list, gene_list, detail_fp);

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


static void  All_Frame_Score
    (const string & s, int len, int frame, vector <double> & af)

//  Score the first  len  characters of string  s  in all six reading
//  frames using global model  Gene_ICM .   frame  is the
//  frame position in the original genome of the first character of
//   s , where frame positions are numbered  1,2,3,1,2,3  starting
//  with the first character of the genome.  frame also has the
//  direction of the gene in the genome string.
//  **NOTE**  s  is the reverse (but not complemented) of the gene.
//  Store the results in  af  where the order of reading frames
//  is  +1,+2,+3,-1,-2,-3 .   af  is assumed to be large enough
//  to hold the results.

  {
   string  rev_compl;
   const char  * cstr = s . c_str ();

   af [0] = Gene_ICM . Score_String (cstr, len, 1);
   af [1] = Gene_ICM . Score_String (cstr, len, 2);
   af [2] = Gene_ICM . Score_String (cstr, len, 0);

   Reverse_Complement_Transfer (rev_compl, s, 0, len);

   af [3] = Gene_ICM . Score_String (rev_compl . c_str (), len, 1);
   af [4] = Gene_ICM . Score_String (rev_compl . c_str (), len, 0);
   af [5] = Gene_ICM . Score_String (rev_compl . c_str (), len, 2);

   Permute_By_Frame (af, frame);

   return;
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
   fprintf (fp, "Excluded regions file = %s\n",
        Printable (Ignore_File_Name));
   fprintf (fp, "List of orfs file = %s\n",
        Printable (Orflist_File_Name));

   fprintf (fp, "Input %s separate orfs\n",
        Separate_Orf_Input ? "is" : "is NOT");
   fprintf (fp, "Independent (non-coding) scores %s used\n",
        Use_Independent_Score ? "are" : "are NOT");
   if  (! Separate_Orf_Input)
       {
        fprintf (fp, "Circular genome = %s\n", Printable (Genome_Is_Circular));
       }
   if  (! Separate_Orf_Input && Orflist_File_Name == NULL)
       {
        fprintf (fp, "Truncated orfs = %s\n", Printable (Allow_Truncated_Orfs));
        fprintf (fp, "Minimum gene length = %d bp\n", Min_Gene_Len);
        fprintf (fp, "Maximum overlap bases = %d\n", Max_Olap_Bases);
        fprintf (fp, "Threshold score = %d\n", Threshold_Score);
        fprintf (fp, "Use first start codon = %s\n",
             Printable (Use_First_Start_Codon));
        if  (Genbank_Xlate_Code != 0)
            fprintf (fp, "Translation table = %d\n", Genbank_Xlate_Code);
        fprintf (fp, "Start codons = ");
        Print_Comma_Separated_Strings (Start_Codon, fp);
        fputc ('\n', fp);
        fprintf (fp, "Start probs = ");
        n = Start_Prob . size ();
        for  (i = 0;  i < n;  i ++)
          {
           if  (i > 0)
               fputc (',', fp);
           fprintf (fp, "%.3f", Start_Prob [i]);
          }
        fputc ('\n', fp);
        fprintf (fp, "Stop codons = ");
        Print_Comma_Separated_Strings (Stop_Codon, fp);
        fputc ('\n', fp);
       }

   fprintf (fp, "GC percentage = %.1f%%\n", 100.0 * Indep_GC_Frac);
   if  (Use_Independent_Score)
       fprintf (fp, "Ignore score on orfs longer than %s\n",
            Num_Or_Max (Ignore_Score_Len));

   return;
  }


static double  Entropy_Distance_Ratio
    (int start, int len, int fr)

//  Return the distance ratio for the entropy profile for the
//  gene starting at position  start  (in 1-based coordinates)
//  on global  Sequence with length  len  and in reading frame  fr .
//  The ratio is the distance to global  Pos_Entropy_Profile  over
//  the distance to global  Neg_Entropy_Profile .

  {
   string  buff;
   int  count [26] = {0};
   double  ep [20];
   double  pos_dist, neg_dist, ratio;
   char  aa;
   int  i;

   if  (fr > 0)
       Forward_Strand_Transfer (buff, Sequence, On_Seq_0 (start - 1), len);
     else
       Reverse_Strand_Transfer (buff, Sequence, On_Seq_0 (start - 1), len);

   for  (i = 0; i < len;  i += 3)
     {
      aa = Codon_Translation (buff . c_str () + i, Genbank_Xlate_Code);
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


static void  Find_Stops_Reverse
    (const string & s, int len, vector <bool> & has_stop)

//  Set  has_stop [i]  to true iff string  s  has a
//  stop codon in the frame corresponding to  i .
//  The order of frames is  +1,+2,+3,-1,-2,-3 .
//  Use only the first  len  characters of  s .
//   s  is the reverse (but not complemented) of the DNA strand
//  Automatically set  has_stop [6]  to  false, representing the
//  independent model "frame".

  {
   Codon_t  codon;
   int  frame_ss;    // frame subscript
   int  which;
   int  i;

   has_stop . resize (7);
   for  (i = 0;  i < 7;  i ++)
     has_stop [i] = false;

   frame_ss = 1;

   for  (i = len - 1;  i >= 0;  i --)
     {
      codon . Shift_In (s [i]);

      if  (codon . Must_Be (Fwd_Stop_Pattern, which))
          has_stop [frame_ss] = true;
      if  (codon . Must_Be (Rev_Stop_Pattern, which))
          has_stop [frame_ss + 3] = true;

      if  (frame_ss == 2)
          frame_ss = 0;
        else
          frame_ss ++;
     }

   return;
  }


static void  Fix_Wrap
    (int & p, const int n)

//  Change position  p  so that it falls in the interval  1 .. n
//  where it should be assuming a circular coordinate scheme.

  {
   while  (p < 1)
     p += n;

   while  (p > n)
     p -= n;

   return;
  }


static void  Get_Orf_Pos_List
    (void)

//  Read the list of orfs from the file with name in global
//   Orflist_File_Name  and store them in global list
//   Orf_Pos_List .  The format for each
//  line of input is:
//     <tag>  <start>  <stop>  <dir>  <rest of line ignored>  
//  where <start> and <stop> are integer values.  The <stop> position
//  includes the ending stop codon for the orf.  The orf specified
//  is bases <start>..<stop> inclusive, where bases in the input
//  sequence are numbered starting at 1.  <dir> indicates the
//  strand of the gene for cases where it might wraparound the
//  start position of the genome sequence.
//  Blank lines and lines beginning with # are skipped.

  {
   FILE  * fp;
   char  line [MAX_LINE], t [MAX_LINE];
   Orf_Pos_t  orf;
   int  line_ct;

   fp = File_Open (Orflist_File_Name, "r", __FILE__, __LINE__);

   Orf_Pos_List . clear ();
   line_ct = 0;
   while  (fgets (line, MAX_LINE, fp) != NULL)
     {
      char  * p;
      int  a, b, d;

      line_ct ++;

      // set  p  to point to the first non-blank character on the line
      for  (p = line;  * p != '\0' && isspace (* p);  p ++)
        ;
      
      if  (* p == '\0' || * p == '#')
          continue;
      else if  (sscanf (line, "%s %d %d %d", t, & a, & b, & d) == 4)
          {
           orf . tag = strdup (t);
           orf . start = a;
           orf . stop = b;
           orf . dir = d;
           Orf_Pos_List . push_back (orf);
          }
        else
          {
           fprintf (stderr, "ERROR:  Following line %d in file %s is bad--skipped:\n",
                line_ct, Orflist_File_Name);
           fputs (line, stderr);
           fputc ('\n', stderr);
          }
     }

   fclose (fp);

   return;
  }


static void  Integerize_Scores
    (const vector <double> ds, int hi_score, const vector <bool> set_negative,
    vector <int> & is)

//  Convert the scores in  ds  to integers ranging from
//   0 .. hi_score  putting the results into  is .
//  Automatically set to  -1  entries corresponding
//  to values in  set_negative  that are true and ignore them
//  in the calculation.

  {
   vector <double>  v;
   double  min, max, sum;
   int  i, n;

   n = ds . size ();
   is . resize (n);
   v . resize (n);

   min = DBL_MAX;
   max = - DBL_MAX;
   for  (i = 0;  i < n;  i ++)
     if  (! set_negative [i])
         {
          if  (ds [i] > max)
              max = ds [i];
          if  (ds [i] < min)
              min = ds [i];
         }

   if  (min < max + MAX_LOG_DIFF)
       min = max + MAX_LOG_DIFF;
   
   sum = 0.0;
   for  (i = 0;  i < n;  i ++)
     if  (set_negative [i])
         v [i] = -1.0;
     else if  (ds [i] < min)
         v [i] = 0.0;
       else
         {
          v [i] = exp (ds [i] - min);
          sum += v [i];
         }

   for  (i = 0;  i < n;  i ++)
     if  (set_negative [i])
         is [i] = -1;
       else
         {
          is [i] = int (HI_SCORE * (v [i] / sum));
          if  (is [i] >= HI_SCORE)
              is [i] = HI_SCORE - 1;
         }

   return;
  }


static int  On_Seq_0
    (int i)

//  Return the subscript equivalent to  i  on a sequence of
//  length  Sequence_Len  (with subscripts starting at 0)
//  assuming circular wraparounds.

  {
   while  (i < 0)
     i += Sequence_Len;
   while  (Sequence_Len <= i)
     i -= Sequence_Len;

   return  i;
  }


static void  Output_Extra_Start_Info
    (FILE * fp, int i, int lo, int hi, int frame,
     vector <Start_t> & start_list)

//  Print to  fp  additional information about the start sites
//  in  start_list .   i  is the subscript of the orf, and  lo .. hi
//  are its tweeny coordinates.   frame  is the reading frame of the
//  orf.

// DRK: I broke this my meddling with the PWMs. To fix it,
// just follow what I did for Add_Events.

  {
   int  stop_pos;
   int  q, r;

   if  (i == 0)
       printf (">%s\n", Fasta_Header);

   if  (frame > 0)
       stop_pos = hi + 3;
     else
       stop_pos = lo - 2;

   Fix_Wrap (stop_pos, Sequence_Len);
   printf ("# %7d  %+2d\n", stop_pos, frame);

   r = start_list . size ();
   for  (q = 0;  q < r && q < 10;  q ++)
     {
      double  score, combined_score;
      int  j, sep;

      j = start_list [q] . pos;

      if  (frame > 0)
          {
           PWM_Score_Fwd_Start (j, LogOdds_PWM, Ribosome_Window_Size, score, sep);
           combined_score = start_list [q] . score
                + Start_Prob [start_list [q] . which];
           if  (score > 0.0)
                combined_score += score;
           printf ("  %7d %c%c%c %7.3f %7.3f %7.3f %3d\n", j, Sequence [On_Seq_0 (j - 1)],
                Sequence [On_Seq_0 (j)], Sequence [On_Seq_0 (j + 1)],
                start_list [q] . score, score, combined_score, sep);
          }
        else
          {
           PWM_Score_Rev_Start (j, LogOdds_PWM, Ribosome_Window_Size, score, sep);
           combined_score = start_list [q] . score
                + Start_Prob [start_list [q] . which];
           if  (score > 0.0)
                combined_score += score;
           printf ("  %7d %c%c%c %7.3f %7.3f %7.3f %3d\n", j,
                Complement (Sequence [On_Seq_0 (j - 1)]),
                Complement (Sequence [On_Seq_0 (j - 2)]),
                Complement (Sequence [On_Seq_0 (j - 3)]),
                start_list [q] . score, score, combined_score, sep);
          }
     }

   return;
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
        {"start_codons", 1, 0, 'A'},
        {"rbs_pwm", 1, 0, 'b'},
        {"gc_percent", 1, 0, 'C'},
        {"entropy", 1, 0, 'E'},
        {"first_codon", 0, 0, 'F'},
	{"features", 1, 0, 'f'},
        {"gene_len", 1, 0, 'g'},
        {"help", 0, 0, 'h'},
        {"ignore", 1, 0, 'g'},
        {"linear", 0, 0, 'l'},
        {"orf_coords", 1, 0, 'L'},
	{"icm", 1, 0, 'm'},
        {"separate_genes", 1, 0, 'M'},
        {"no_indep", 0, 0, 'n'},
        {"max_olap", 1, 0, 'o'},
        {"start_probs", 1, 0, 'P'},
        {"ignore_score_len", 1, 0, 'q'},
        {"threshold", 1, 0, 't'},
	{"fudge", 1, 0, 'u'},
        {"extend", 0, 0, 'X'},
        {"trans_table", 1, 0, 'z'},
        {"stop_codons", 1, 0, 'Z'},
        {0, 0, 0, 0}
      };

   while  (! errflg && ((ch = getopt_long (argc, argv,
        "A:b:C:E:f:Fg:hi:lL:m:Mno:P:q:t:u:Xz:Z:",
        long_options, & option_index)) != EOF))
#else
   while  (! errflg && ((ch = getopt (argc, argv,
        "A:b:C:E:f:Fg:hi:lL:m:Mno:P:q:t:u:Xz:Z:")) != EOF))
#endif

     switch  (ch)
       {
        case  'A' :
          Command_Line . append (" -A ");
          Command_Line . append (optarg);
          Start_Codon . clear ();
          for  (p = strtok (optarg, ",");  p != NULL;  p = strtok (NULL, ","))
            {
             q = strdup (p);
             Make_Lower_Case (q);
             Start_Codon . push_back (q);
            }
          break;

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

        case  'C' :
          Command_Line . append (" -C ");
          Command_Line . append (optarg);
          Indep_GC_Frac = strtod (optarg, & p) / 100.0;
          if  (p == optarg || Indep_GC_Frac < 0.0 || Indep_GC_Frac > 100.0)
              {
               fprintf (stderr, "ERROR:  Bad independent model GC fraction (-C option)\n"
                    "  value = \"%s\"", optarg);
               errflg = true;
              }
          GC_Frac_Set = true;
          break;

        case  'E' :
          Command_Line . append (" -E ");
          Command_Line . append (optarg);
          if  (strcmp (optarg, "#") != 0)
              Read_Entropy_Profiles (optarg, errflg);
          Use_Entropy_Profiles = true;
          break;
	  
       case  'f' :
	 Command_Line . append (" -f");
	 Use_First_Start_Codon = true;
	 break;

       case 'F' :
	    Command_Line.append(" -F");
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
          Command_Line . append (" -i ");
          Command_Line . append (optarg);
          Ignore_File_Name = optarg;
          break;

        case  'l' :
          Command_Line . append (" -l");
          Genome_Is_Circular = false;
          break;

        case  'L' :
          Command_Line . append (" -L ");
          Command_Line . append (optarg);
          Orflist_File_Name = optarg;
          break;
	 
       case  'm' :
	    Command_Line . append (" -m ");
	    Command_Line . append (optarg);
	    ICM_File_Name = optarg;
	    break;

        case  'M' :
          Command_Line . append (" -M");
          Separate_Orf_Input = true;
          break;

       case  'n' :
          Command_Line . append (" -n");
          Use_Independent_Score = false;
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

        case  'P' :
          Command_Line . append (" -P ");
          Command_Line . append (optarg);
          Start_Prob . clear ();
          for  (p = strtok (optarg, ",");  p != NULL;  p = strtok (NULL, ","))
            Start_Prob . push_back (strtod (p, NULL));
          break;

        case  'q' :
          Command_Line . append (" -q ");
          Command_Line . append (optarg);
          Ignore_Score_Len = strtol (optarg, & p, 10);
          if  (p == optarg || Ignore_Score_Len < 0)
              {
               fprintf (stderr, "ERROR:  Bad ignore independent model length\n"
                    "  (-q option)  value = \"%s\"", optarg);
               errflg = true;
              }
          break;

        case  't' :
          Command_Line . append (" -t ");
          Command_Line . append (optarg);
          Threshold_Score = strtol (optarg, & p, 10);
          if  (p == optarg || Threshold_Score <= 0 || Threshold_Score >= 100)
              {
               fprintf (stderr, "ERROR:  Bad threshold score (-t option)\n"
                    "  value = \"%s\"", optarg);
               errflg = true;
              }
          break;

       case  'u' :
	  Command_Line . append (" -u ");
          Command_Line . append (optarg);
          LogOdds_Fudge = strtod (optarg, & p);
          if  (p == optarg)
              {
               fprintf (stderr, "ERROR:  Bad value for fudge factor (-u option)\n"
                    "  value = \"%s\"", optarg);
               errflg = true;
              }
	  LogOdds_Prior += LogOdds_Fudge;
          break;

        case  'X' :
          Command_Line . append (" -X");
          Allow_Truncated_Orfs = true;
          Genome_Is_Circular = false;
          break;

        case  'z' :
          Command_Line . append (" -z ");
          Command_Line . append (optarg);
          Genbank_Xlate_Code = strtol (optarg, & p, 10);
          Set_Stop_Codons_By_Code (Stop_Codon, Genbank_Xlate_Code, errflg);
          break;

        case  'Z' :
          Command_Line . append (" -Z ");
          Command_Line . append (optarg);
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

   return;
  }


template  <class DT>
static void  Permute_By_Frame
    (vector <DT> & v, int frame)

//  Permute the first 6 entries in  v  so that they
//  represent a reverse sequence of a gene, where the first
//  base of the sequence comes from genome position with
//  frame  frame .  Positions of the genome are numbered 1,2,3,1,2,3...
//  Frame is positive for forward strand genes in the genome and negative
//  for reverse strand genes.  The input values in  v  represent
//  scores for a frame  +3  sequence.

  {
   DT  save;

   switch  (frame)
     {
      case  1 :
        save = v [0];
        v [0] = v [2];
        v [2] = v [1];
        v [1] = save;
        save = v [3];
        v [3] = v [5];
        v [5] = v [4];
        v [4] = save;
        break;
      case  2 :
        save = v [0];
        v [0] = v [1];
        v [1] = v [2];
        v [2] = save;
        save = v [3];
        v [3] = v [4];
        v [4] = v [5];
        v [5] = save;
        break;
      case  3 :
        break;
      case  -1 :
        save = v [0];
        v [0] = v [3];
        v [3] = save;
        save = v [1];
        v [1] = v [5];
        v [5] = save;
        save = v [2];
        v [2] = v [4];
        v [4] = save;
        break;
      case  -2 :
        save = v [0];
        v [0] = v [4];
        v [4] = save;
        save = v [1];
        v [1] = v [3];
        v [3] = save;
        save = v [2];
        v [2] = v [5];
        v [5] = save;
        break;
      case  -3 :
        save = v [0];
        v [0] = v [5];
        v [5] = save;
        save = v [1];
        v [1] = v [4];
        v [4] = save;
        save = v [2];
        v [2] = v [3];
        v [3] = save;
        break;
     }

   return;
  }


static void  Print_Orflist_Headings
    (FILE * fp)

//  Print column headings for separate orf list (-L option) to  fp .

  {
   fputc ('\n', fp);

   fprintf (fp, "%-12s %5s  %8s %8s %8s", "", "", "", "", "");
   if  (Use_Independent_Score)
       fprintf (fp, "  %s\n", "------------- Scores -------------");
     else
       fprintf (fp, "  %s\n", "----------- Scores ------------");
   fprintf (fp, "%-12s %5s  %8s %8s %8s  %7s %5s %s",
        "  ID", "Frame", "Start", "Stop", "Len", "Raw", "InFrm", "F1 F2 F3 R1 R2 R3");
   if  (Use_Independent_Score)
       fprintf (fp, " NC");
   if  (Use_Entropy_Profiles)
       fprintf (fp, " %-4s", "EDR");
   fprintf (fp, "\n");

   return;
  }


static void  Read_Entropy_Profiles
    (const char * fn, bool & errflg)

//  Read positive and negative entropy profiles from the
//  file name  fn .  If not successful, set  errflg  to  true .
//  Save the entropy profiles in globals  Pos_Entropy_Profile
//  and  Neg_Entropy_Profile .

  {
   FILE  * fp;
   char  line [MAX_LINE];
   int  i;

   fp = File_Open (fn, "r");
   fgets (line, MAX_LINE, fp);  // skip header line
   for  (i = 0;  i < 20;  i ++)
     if  (fscanf (fp, "%s %lf %lf\n", line, Pos_Entropy_Profile + i,
             Neg_Entropy_Profile + i) != 3)
         {
          errflg = true;
          return;
         }

   fclose (fp);

   return;
  }


static void  Read_Sequences
    (FILE * fp, vector <string> & seq_list, vector <string> & hdr_list,
     int & seq_ct)

//  Read fasta-format sequences from  fp  (which is already open),
//  convert them to lower-case, and store them in  seq_list .
//  Store the fasta header lines in  hdr_list .  Set  seq_ct  to
//  the number of sequences read.

  {
   string  seq, hdr;
   int  i, len;

   seq_list . clear ();
   hdr_list . clear ();
   seq_ct = 0;

   while  (Fasta_Read (fp, seq, hdr))
     {
      len = seq . length ();
      for  (i = 0;  i < len;  i ++)
        seq [i] = Filter (tolower (seq [i]));

      seq_list . push_back (seq);
      hdr_list . push_back (hdr);
      seq_ct ++;
     }

   return;
  }


static void  Score_Orflist
    (FILE * detail_fp, FILE * summary_fp)

//  Score the entries in global  Orf_Pos_List  using the sequence
//  in global  Sequence  sending detailed results to  detail_fp and
//  summary results to  summary_fp.

  {
   string  buff;
   vector <double>  af, score, indep_score;
   vector <int>  int_score;
   vector <bool>  has_stop;
   int  fr, frame, frame_score;
   int  lo, hi, len;
   int  i, j, m, n;

   if  (Use_Independent_Score)
       af . resize (7);
     else
       af . resize (6);

   n = Orf_Pos_List . size ();
   for  (i = 0;  i < n;  i ++)
     {
      double  gene_score;
      int  start, stop;

      start = Orf_Pos_List [i] . start;
      stop = Orf_Pos_List [i] . stop;

      if  (Orf_Pos_List [i] . dir > 0)
          {
           frame = 1 + (stop % 3);
           fr = 1 + (1 + frame) % 3;
           len = 1 + stop - start - 3;
           if  (len < 0)
               len += Sequence_Len;
           hi = stop - 3;
           if  (hi <= 0)
               hi += Sequence_Len;
           Reverse_Transfer (buff, Sequence, hi - 1, len);
          }
        else
          {
           fr = frame = - ((stop - 1) % 3) - 1;
           len = 1 + start - stop - 3;
           if  (len < 0)
               len += Sequence_Len;
           lo = stop + 2;
           if  (lo >= Sequence_Len)
               lo -= Sequence_Len;
           Complement_Transfer (buff, Sequence, lo, len);
          }

      
      Gene_ICM . Cumulative_Score (buff, score, 1);
      Indep_Model . Cumulative_Score (buff, indep_score, 1);
      m = score . size ();

      if  (Use_Independent_Score)
          af [6] = indep_score [m - 4];   // excludes the start codon
      All_Frame_Score (buff, m - 3, fr, af);
      Find_Stops_Reverse (buff, m - 3, has_stop);
      gene_score = 100.0 * (score [m - 4] - indep_score [m - 4]) / (m - 3);

      Permute_By_Frame (has_stop, fr);
      Integerize_Scores (af, HI_SCORE, has_stop, int_score);
      if  (frame > 0)
          frame_score = int_score [frame - 1];
        else
          frame_score = int_score [2 - frame];

      // print score details
      if(Detail_Log) {
	   fprintf (detail_fp, "%-14s %+3d  %8d %8d %8d  %7.2f %5d",
		    Orf_Pos_List [i] . tag, frame, start, stop, len, gene_score, frame_score);
	   for  (j = 0;  j < 6;  j ++)
		if  (int_score [j] < 0)
		     fprintf (detail_fp, "  -");
		else
		     fprintf (detail_fp, " %2d", int_score [j]);
	   if  (Use_Independent_Score)
		fprintf (detail_fp, " %2d", int_score [6]);
	   if  (Use_Entropy_Profiles)
		fprintf (detail_fp, " %4.2f", Entropy_Distance_Ratio (start, m, frame));
	   fputc ('\n', detail_fp);
      }

      // print overall score
      fprintf (summary_fp, "%-14s %8d %8d %+3d %8.2f\n",
          Orf_Pos_List [i] . tag, start, stop, frame, gene_score);
     }

   return;
  }



static void  Score_Orfs
    (vector <Orf_t> & orf_list, vector <Gene_t> & gene_list, FILE * fp)

//  Compute scores for all orfs in  orf_list  using coding model
//  in global  Gene_ICM , which is assumed to have been built on reverse
//  gene strings.   Indep_Model  is the model of independent,
//  stop-codon-free sequence.  Put orfs that are candidate genes
//  onto  gene_list .  Print log information to  fp .

  {
   string  buff;
   vector <double>  af, score, indep_score;
   vector <bool>  is_start;
   vector <Start_t>  start_list;
   Start_t  start;
   char  tag [MAX_LINE];		
   int  i, n, id = 0;

   if  (Use_Independent_Score)
       af . resize (7);
     else
       af . resize (6);

   gene_list . clear ();

   n = orf_list . size ();
   for  (i = 0;  i < n;  i ++)
     {
      double  first_score = - DBL_MAX; // to avoid warning
      double best_score = - DBL_MAX;
      double  gene_score;
      vector <int>  int_score;
      vector <bool>  has_stop;
      int  first_pos = 0, best_pos = 0;
      int  first_j = 0, best_j = 0;
      double  max, max_rate;
      Codon_t  codon;
      double  s;
      bool  is_tentative_gene, orf_is_truncated = false;
      bool  first_is_truncated = false, best_is_truncated = false;
      int  which;
      int  fr, frame, max_j, orf_start, frame_score;
      int  lo, hi, len, lowest_j;
      int  j, k, m;

      frame = orf_list [i] . Get_Frame ();
      len = orf_list [i] . Get_Orf_Len ();
      if  (frame > 0)
          {
           hi = orf_list [i] . Get_Stop_Position () - 1;
           if  (hi <= 0)
               hi += Sequence_Len;
           lo = hi - len;
           Reverse_Transfer (buff, Sequence, hi - 1, len);
           fr = 1 + (1 + frame) % 3;
           orf_is_truncated = (lo < 3 && Allow_Truncated_Orfs);
           k = orf_list [i] . Get_Stop_Position () - len - 2;
          }
        else
          {
           lo = orf_list [i] . Get_Stop_Position () + 2;
           if  (lo >= Sequence_Len)
               lo -= Sequence_Len;
           hi = lo + len;
           Complement_Transfer (buff, Sequence, lo, len);
           fr = frame;
           orf_is_truncated = (Sequence_Len - hi < 3 && Allow_Truncated_Orfs);
           k = orf_list [i] . Get_Stop_Position () + len + 4;
          }
      // lo .. hi  are the between coordinates of the orf region.

      Gene_ICM . Cumulative_Score (buff, score, 1);
      Indep_Model . Cumulative_Score (buff, indep_score, 1);
      m = score . size ();

      max = 0.0;
      max_j = 0;
      is_start . resize (m, false);
      start_list . clear ();
      lowest_j = Min (3, Min_Gene_Len - 3);
      for  (j = m - 1;  j >= lowest_j;  j --)
        {
         codon . Shift_In (buff [j]);

	 // DRK
	 // not used at all?
         s = score [j] - indep_score [j];
         if  (s > max)
             {
              max = s;
              max_rate = s / (j + 1);
              max_j = j;
             }

         if  (j % 3 == 0
                 && (codon . Can_Be (Fwd_Start_Pattern, which)
                       || (first_pos == 0 && orf_is_truncated))
                 && j + 3 >= Min_Gene_Len)
             {
              double  next_s;

              next_s = score [j - 1] - indep_score [j - 1];
                // this is the score for the orf without the start
                // codon--position j is the last base of the start codon
              is_start [j + 2] = true;
              start . j = j + 2;
              start . pos = k;
                // k is the 1-based sequence coordinate of the base that
                // is 2 behind the position represented by j
              start . score = next_s;
              start . first = (first_pos == 0);

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

              if  (first_pos == 0)
                  {
                   first_score = next_s;
                   first_pos = k;
                   first_j = j + 2;
                   //first_is_truncated = start . truncated;
		   first_is_truncated = (first_pos == 0 && orf_is_truncated);
                  }
              if  (next_s > best_score)
                  {
                   best_score = next_s;
                   best_pos = k;
                   best_j = j + 2;
                   best_is_truncated = start . truncated;
                  }
             }
         if  (frame > 0)
             k ++;
           else
             k --;
        }

      if  (Use_First_Start_Codon)
          {
           best_score = first_score;
           best_pos = first_pos;
           best_j = first_j;
           best_is_truncated = first_is_truncated;
          }

      if  (first_j + 1 < Min_Gene_Len)
          continue;

      if  (frame > 0)
          {
           k = hi + 3;
           orf_start = lo + 1;
          }
        else
          {
           k = lo - 2;
           orf_start = hi;
          }

      if  (Use_Independent_Score)
          af [6] = indep_score [best_j - 3];

//**ALD  Changed  best_j + 1  to  best_j - 2  to omit start codon
//  from score to be consistent with the independent score
      All_Frame_Score (buff, best_j - 2, fr, af);
      Find_Stops_Reverse (buff, best_j - 2, has_stop);

      Permute_By_Frame (has_stop, fr);
      Integerize_Scores (af, HI_SCORE, has_stop, int_score);
      if  (frame > 0)
          frame_score = int_score [frame - 1];
        else
          frame_score = int_score [2 - frame];


      // this filtering works better...

      // boost long ORFs
	  for(unsigned int start_i = 0; start_i < start_list.size(); start_i++)
	       if(start_list[start_i].j > Ignore_Score_Len)
		    start_list[start_i].score = Max(0.0, start_list[start_i].score);

      is_tentative_gene = (first_j + 1 >= Min_Gene_Len && best_score > Start_Threshold);

      /* then this filtering...

      // For now just use the score, will add more later
      is_tentative_gene
           = (best_j + 1 >= Min_Gene_Len && frame_score >= Threshold_Score);


      // If it's long enough to ignore the independent score,
      // rescue it
      if  (! is_tentative_gene && len >= Ignore_Score_Len)
          {
           best_score = first_score;
           best_pos = first_pos;
           best_j = first_j;
           is_tentative_gene = true;
          }
      */

//**ALD  Changed to omit start codon
      gene_score = 100.0 * best_score / (best_j - 2);

      //if  (For_Edwin)
      //    Output_Extra_Start_Info (stdout, i, lo, hi, frame, start_list);

      if  (is_tentative_gene)
          {
           sprintf (tag, "%04d", ++ Gene_ID_Ct);
           Gene_t  gene (orf_list [i]);
           gene . Set_Score (gene_score);
           gene . Set_Gene_Len (best_j + 1);
           gene_list . push_back (gene);
          }
        else
          strcpy (tag, "    ");

      if  (Genome_Is_Circular)
          {
           Fix_Wrap (orf_start, Sequence_Len);
           Fix_Wrap (best_pos, Sequence_Len);
           Fix_Wrap (k, Sequence_Len);
          }
      else if  (orf_is_truncated)
      {
           if  (frame > 0)
	   {
                orf_start -= 3;
                if  (best_is_truncated)
		     best_pos -= 3;
	   }
	   else
	   {
                orf_start += 3;
                if  (best_is_truncated)
		     best_pos += 3;
	   }
      }

      if(Detail_Log) {
	   fprintf (fp, "%4s %+5d %8d %8d %8d  %7d %7d  %7.2f %5d",
		    tag, frame, orf_start, best_pos, k, len, best_j + 1,
		    gene_score, frame_score );
	   for  (j = 0;  j < 6;  j ++)
		if  (int_score [j] < 0)
		     fprintf (fp, "  -");
		else
		     fprintf (fp, " %2d", int_score [j]);
	   if  (Use_Independent_Score)
		fprintf (fp, " %2d", int_score [6]);
	   if  (Use_Entropy_Profiles)
		fprintf (fp, " %4.2f", Entropy_Distance_Ratio (best_pos,
							       best_j + 1, frame));
	   fputc ('\n', fp);
      }

      if  (is_tentative_gene)
	   if(frame > 0)
		Add_Events_Fwd (orf_list [i], start_list, id);
	   else
		Add_Events_Rev (orf_list [i], start_list, id);
     }
   
   return;
  }


static void  Score_Separate_Input
    (const string & seq, const string & hdr, int seq_num, FILE * detail_fp,
     FILE * predict_fp)

//  Score the sequence  seq  with fasta header  hdr  in frame and output
//  the results to  detail_fp  and  predict_fp .

  {
   string  buff;
   vector <double>  af, score, indep_score;
   char  line [MAX_LINE], tag [MAX_LINE], * p;
   vector <int>  int_score;
   vector <bool>  has_stop;
   double  gene_score;
   int  fr, frame, frame_score;
   int  len;
   int  j, m;

   len = seq . length () - 3;  // remove stop codon
   Reverse_Transfer (buff, seq, len - 1, len);
   strcpy (line, hdr . c_str ());
   p = strtok (line, " \t\n");
   if  (p == NULL)
       sprintf (tag, "Seq%04d", seq_num);
     else
       strcpy (tag, p);

   if  (Use_Independent_Score)
       af . resize (7);
     else
       af . resize (6);
   
   frame = 1;  // assume all orfs are in correct reading frame
   fr = 3;     // shifted number for this frame

   Gene_ICM . Cumulative_Score (buff, score, 1);
   Indep_Model . Cumulative_Score (buff, indep_score, 1);
   len = m = score . size ();

   if  (Use_Independent_Score)
       af [6] = indep_score [m - 4];   // excludes the start codon
   All_Frame_Score (buff, m - 3, fr, af);
   Find_Stops_Reverse (buff, m - 3, has_stop);
   gene_score = 100.0 * (score [m - 4] - indep_score [m - 4]) / (m - 3);

   Permute_By_Frame (has_stop, fr);
   Integerize_Scores (af, HI_SCORE, has_stop, int_score);
   if  (frame > 0)
       frame_score = int_score [frame - 1];
     else
       frame_score = int_score [2 - frame];

   // print score details
   if(Detail_Log) {
	fprintf (detail_fp, "%-14s %+3d  %8d %8d %8d  %7.2f %5d",
		 tag, frame, 1, len, len, gene_score, frame_score);
	for  (j = 0;  j < 6;  j ++)
	     if  (int_score [j] < 0)
		  fprintf (detail_fp, "  -");
	     else
		  fprintf (detail_fp, " %2d", int_score [j]);
	if  (Use_Independent_Score)
	     fprintf (detail_fp, " %2d", int_score [6]);
	if  (Use_Entropy_Profiles)
	     fprintf (detail_fp, " %4.2f", Entropy_Distance_Ratio (1, m, frame));
	fputc ('\n', detail_fp);
   }

   // print overall score
   fprintf (predict_fp, "%-14s %8d %8d %+3d %8.2f\n",
       tag, 1, len, frame, gene_score);

   return;
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
   for  (i = 0;  i < n;  i ++)
     {
      if  (Genome_Is_Circular)
          {
           j = On_Seq_1 (gene_list [i] . Get_Stop_Position ());
           gene_list [i] . Set_Stop_Position (j);
          }
        else
          j = gene_list [i] . Get_Stop_Position ();
      f = Position_To_Frame (j);
      if  (gene_list [i] . Get_Frame () > 0)
          gene_list [i] . Set_Frame (f);
        else
          gene_list [i] . Set_Frame (-1 * f);
     }

   //sort (gene_list . begin (), gene_list . end (), By_ID);

   //for  (i = 0;  i < n;  i ++)
   for  (i = n-1;  i >= 0;  i --)
     {
      int  start, stop;

      if  (gene_list [i] . Get_Frame () > 0)
          {
           if  (Genome_Is_Circular)
               {
                stop = On_Seq_1 (gene_list [i] . Get_Stop_Position () + 2);
                start = On_Seq_1 (stop - gene_list [i] . Get_Gene_Len () - 2);
               }
             else
               {
                stop = gene_list [i] . Get_Stop_Position () + 2;
                start = stop - gene_list [i] . Get_Gene_Len () - 2;
                if  (gene_list [i] . Get_Status_Bit (TRUNCATED_START_FLAG))
                    start -= 3;
                  // move an artificial start at the beginning of the sequence
                  // off the front to indicate the gene could extend there
               }
          }
        else
          {
           if  (Genome_Is_Circular)
               {
                stop = On_Seq_1 (gene_list [i] . Get_Stop_Position ());
                start = On_Seq_1 (stop + gene_list [i] . Get_Gene_Len () + 2);
               }
             else
               {
                stop = gene_list [i] . Get_Stop_Position ();
                start = stop + gene_list [i] . Get_Gene_Len () + 2;
                if  (gene_list [i] . Get_Status_Bit (TRUNCATED_START_FLAG))
                    start += 3;
                  // move an artificial start at the end of the sequence
                  // off the back to indicate the gene could extend there
               }
          }

      	// print gene
	fprintf (fp, "orf%05d %8d %8d %+3d %8.2f\n",
		 gene_list [i] . Get_ID (),  start, stop,
		 gene_list [i] . Get_Frame (),
		 gene_list [i] . Get_Score ());
     }

   return;
  }


static void  Usage
    (void)

//  Print to stderr description of options and command line for
//  this program.

  {
   fprintf (stderr,
       "USAGE:  glimmer3 [options] <sequence-file> <icm-file> <tag>\n"
       "\n"
       "Read DNA sequences in <sequence-file> and predict genes\n"
       "in them using the Interpolated Context Model in <icm-file>.\n"
       "Output details go to file <tag>.detail and predictions go to\n"
       "file <tag>.predict\n"
       "\n"
       "Options:\n"
       " -A <codon-list>\n"
       " --start_codons <codon-list>\n"
       "    Use comma-separated list of codons as start codons\n"
       "    Sample format:  -A atg,gtg\n"
       "    Use -P option to specify relative proportions of use.\n"
       "    If -P not used, then proportions will be equal\n"
       " -b <filename>\n"
       " --rbs_pwm <filename>\n"
       "    Read a position weight matrix (PWM) from <filename> to identify\n"
       "    the ribosome binding site to help choose start sites\n"
       " -C <p>\n"
       " --gc_percent <p>\n"
       "    Use <p> as GC percentage of independent model\n"
       "    Note:  <p> should be a percentage, e.g., -C 45.2\n"
       " -E <filename>\n"
       " --entropy <filename>\n"
       "    Read entropy profiles from <filename>.  Format is one header\n"
       "    line, then 20 lines of 3 columns each.  Columns are amino acid,\n"
       "    positive entropy, negative entropy.  Rows must be in order\n"
       "    by amino acid code letter\n"
       " -F\n"
       " --first_codon\n"
       "    Use first codon in orf as start codon\n"
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
       " -i <filename>\n"
       " --ignore <filename>\n"
       "    <filename> specifies regions of bases that are off \n"
       "    limits, so that no bases within that area will be examined\n"
       " -l\n"
       " --linear\n"
       "    Assume linear rather than circular genome, i.e., no wraparound\n"
       " -L <filename>\n"
       " --orf_coords <filename>\n"
       "    Use <filename> to specify a list of orfs that should\n"
       "    be scored separately, with no overlap rules\n"
       " -M\n"
       " --separate_genes\n"
       "    <sequence-file> is a multifasta file of separate genes to\n"
       "    be scored separately, with no overlap rules\n"
       " -m <filename>\n"
       " --icm <filename>\n"
       "    Read ICM from <filename> and use to score ORF coding likelihood\n"
       " -o <n>\n"
       " --max_olap <n>\n"
       "    Set maximum overlap length to <n>.  Overlaps this short or shorter\n"
       "    are ignored.\n"
       " -P <number-list>\n"
       " --start_probs <number-list>\n"
       "    Specify probability of different start codons (same number & order\n"
       "    as in -A option).  If no -A option, then 3 values for atg, gtg and ttg\n"
       "    in that order.  Sample format:  -P 0.6,0.35,0.05\n"
       "    If -A is specified without -P, then starts are equally likely.\n"
       " -q <n>\n"
       " --ignore_score_len <n>\n"
       "    Do not use the initial score filter on any gene <n> or more\n"
       "    base long\n"
       " -t <n>\n"
       " --threshold <n>\n"
       "    Set threshold score for calling as gene to n.  If the in-frame\n"
       "    score >= <n>, then the region is given a number and considered\n"
       "    a potential gene.\n"
       " -u\n"
       " --fudge\n"
       "    Value to be added to the log-likelihood ratio score of every ORF.\n"
       " -X\n"
       " --extend\n"
       "    Allow orfs extending off ends of sequence to be scored\n"
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
