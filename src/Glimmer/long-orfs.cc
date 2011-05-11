//  A. L. Delcher
//
//  File:  long-orfs.cc
//
//  Last Modified:  1 Aug 2005
//
//  This program finds sufficiently long, non-overlapping reading
//  frames in the file named on the command line.  Optionally,
//  entropy distance can be used to filter the selected orfs.
//
//  Copyright (c) 2006 University of Maryland Center for Bioinformatics
//  & Computational Biology



#include  "long-orfs.hh"



// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

static double  Entropy_Cutoff = 1.0;
  // Genes with entropy distance score not below this
  // are eliminated before overlaps are considered
static const char  * Fasta_Header;
  // Header on first line of fasta input file
static bool  Fixed_Min_Len = false;
  // If set true by the -f option, then the specified
  // or default minimum length gene will be used.
static vector <Codon_t>  Fwd_Start_Pattern;
  // Bit patterns representing possible forward start codons
static vector <Codon_t>  Fwd_Stop_Pattern;
  // Bit patterns representing possible forward stop codons
static int  Genbank_Xlate_Code = 0;
  // Holds the Genbank translation table number that determines
  // stop codons and codon translation.
static bool  Genome_Is_Circular = DEFAULT_GENOME_IS_CIRCULAR;
  // If true, input sequences are assumed to be circularly connected
  // so genes will be allowed to wrap around the end
static char  * Ignore_File_Name = NULL;
  // Name of file containing list of regions that cannot be included
  // in gene predictions
static vector <Range_t>  Ignore_Region;
  // List of regions to be skipped
static int  Min_Gene_Len = DEFAULT_MIN_GENE_LEN;
  // Shortest (in nucleotides) gene that will be considered for scoring
static int  Max_Olap_Bases = DEFAULT_MAX_OLAP_BASES;
  // Overlaps of this many or fewer bases are allowed between adjacent
  // genes
static double  Neg_Entropy_Profile [20] = DEFAULT_NEG_ENTROPY_PROF;
  // Entropy distribution of amino-acids in non-genes
static bool  Optimize_Total_Len = false;
  // If set true by the -L option, then minimum gene length
  // that optimizes the total length of genes (instead of their
  // number) will be used.  Not compatible with -f option.
static char  * Output_Filename = NULL;
  // Name of file to which output is sent
static double  Pos_Entropy_Profile [20]  = DEFAULT_POS_ENTROPY_PROF;
  // Entropy distribution of amino-acids in genes
static bool  Print_Output_Header = true;
  // Determines if header information is printed before the
  // list of gene coordinates in the output file
static vector <Codon_t>  Rev_Start_Pattern;
  // Bit patterns representing possible reverse start codons
static vector <Codon_t>  Rev_Stop_Pattern;
  // Bit patterns representing possible reverse stop codons
static string  Sequence;
  // The input sequence to be scored.
static char  * Sequence_File_Name = NULL;
  // Name of the input sequence file
static int  Sequence_Len;
  // Length of genomic sequence string being processed.
static vector <const char *>  Start_Codon;
  // Sequences assumed to be start codons
static vector <const char *>  Stop_Codon;
  // Sequences assumed to be stop codons
static string  Tag;
  // The fasta-header lines of the sequence in  Sequence
static bool  Use_Entropy_Filter = false;
  // If set true by the -t option, then use the
  //  Entropy_Cutoff  value to filter orfs before
  // checking for overlaps
static bool  Without_Stops = false;
  // If set true by the -w option, then output
  // coordinates will *NOT* include the stop codon



int  main
    (int argc, char * argv [])

  {
   try
     {
      FILE  * sequence_fp, * output_fp;
      vector <Orf_t>  orf_list;
      vector <Orf_Interval_t>  interval;
      string  hdr;
      time_t  now;
      int  optimal_len, total_len;
      int  i;

      now = time (NULL);
      cerr << "Starting at " << ctime (& now) << endl;

      Verbose = 0;

      Parse_Command_Line (argc, argv);

      if  (Ignore_File_Name != NULL)
          Get_Ignore_Regions ();

      Set_Start_And_Stop_Codons ();

      if  (strcmp (Output_Filename, "-") == 0)
          output_fp = stdout;
        else
          output_fp = File_Open (Output_Filename, "w", __FILE__, __LINE__);

      Echo_General_Settings (stderr);
      if  (Print_Output_Header)
          Echo_General_Settings (output_fp);

      sequence_fp = File_Open (Sequence_File_Name, "r", __FILE__, __LINE__);

      if  (! Fasta_Read (sequence_fp, Sequence, hdr))
          SIMPLE_THROW ("ERROR:  Failed to read input sequence");

      Sequence_Len = Sequence . length ();
      for  (i = 0;  i < Sequence_Len;  i ++)
        Sequence [i]  = Filter (tolower (Sequence [i]));

      Fasta_Header = hdr . c_str ();

      Find_Orfs (orf_list);

      if  (Use_Entropy_Filter)
          Entropy_Filter (orf_list, Entropy_Cutoff);

      if  (orf_list . size () == 0)
          SIMPLE_THROW ("ERROR:  No valid orfs found below entropy cutoff");

      Get_Intervals (interval, orf_list);

      if  (! Fixed_Min_Len)
          {
           optimal_len = Find_Optimal_Len (interval);
           Remove_Shorter (interval, optimal_len);
           Min_Gene_Len = optimal_len;
          }

      Eliminate_Overlapping (interval, Max_Olap_Bases);

      Echo_Specific_Settings (stderr, Sequence_Len);
      if  (Print_Output_Header)
          Echo_Specific_Settings (output_fp, Sequence_Len);

      Output_Orfs (output_fp, interval, total_len);

      fprintf (stderr, "Number of genes = %d\n", int (interval . size ()));
      fprintf (stderr, "Total bases = %d\n", total_len);

      fclose (sequence_fp);
      fclose (output_fp);
     }
   catch (std :: exception & e)
     {
      cerr << "** Standard Exception **" << endl;
      cerr << e << endl;
      exit (EXIT_FAILURE);
     }

   return  0;
  }



static void  Echo_General_Settings
    (FILE * fp)

//  Output values of global variables and parameter settings
//  to  fp .

  {
   fprintf (fp, "Sequence file = %s\n", Sequence_File_Name);
   fprintf (fp, "Excluded regions file = %s\n",
        Printable (Ignore_File_Name));

   fprintf (fp, "Circular genome = %s\n", Printable (Genome_Is_Circular));
   fprintf (fp, "Initial minimum gene length = %d bp\n", Min_Gene_Len);
   if  (Fixed_Min_Len)
       fprintf (fp, "Fixed minimum gene length\n");
     else
       fprintf (fp, "Determine optimal min gene length to maximize %s\n",
            Optimize_Total_Len ? "total bases" : "number of genes");
   fprintf (fp, "Maximum overlap bases = %d\n", Max_Olap_Bases);
   if  (Genbank_Xlate_Code != 0)
       fprintf (fp, "Translation table = %d\n", Genbank_Xlate_Code);
   fprintf (fp, "Start codons = ");
   Print_Comma_Separated_Strings (Start_Codon, fp);
   fputc ('\n', fp);
   fprintf (fp, "Stop codons = ");
   Print_Comma_Separated_Strings (Stop_Codon, fp);
   fputc ('\n', fp);

   return;
  }



static void  Echo_Specific_Settings
    (FILE * fp, int len)

//  Output values of variables an settings that depend on the
//  current input string, which has length  len .

  {
   fprintf (fp, "Sequence length = %d\n", len);
   fprintf (fp, "Final minimum gene length = %d\n", Min_Gene_Len);

   return;
  }



static void  Eliminate_Overlapping
    (vector <Orf_Interval_t> & interval, int max_olap)

//  Eliminate from  interval  entries that overlap other
//  entries by more than  max_olap .

  {
   vector <int>  highest;
   int  right_wrap;
   int  i, j, n;

   n = interval . size ();
   highest . resize (n);
     // highest [i]  will be the max hi value in interval [0 .. i]

   // set values in  highest and  right_wrap
   right_wrap = 0;
   for  (i = 0;  i < n;  i ++)
     {
      if  (i == 0)
          highest [i] = interval [i] . hi;
        else
          highest [i] = Max (highest [i - 1], interval [i] . hi);
      if  (Genome_Is_Circular)
          right_wrap = Max (right_wrap, interval [i] . hi - Sequence_Len);
     }

   highest [0] = interval [0] . hi;
   for  (i = 1;  i < n;  i ++)
     {
      for  (j = i - 1;  0 <= j;  j --)
        {
         if  (highest [j] <= interval [i] . lo + max_olap)
             break;  // can't sufficiently overlap j or anything before it

         if  (max_olap < Intersect_Size (interval [j] . lo, interval [j] . hi,
                interval [i] . lo, interval [i] . hi))
             interval [j] . deleted = interval [i] . deleted = true;
        }

      // also check any wraparounds
      if  (Genome_Is_Circular
             && interval [i] . lo + max_olap <= right_wrap)
          {
           for  (j = n - 1;  j > i && interval [i] . lo + max_olap
                      <= highest [j] - Sequence_Len;  j --)
             if  (max_olap < Intersect_Size (interval [i] . lo,
                    interval [i] . hi, interval [j] . lo - Sequence_Len,
                    interval [j] . hi - Sequence_Len))
                 interval [j] . deleted = interval [i] . deleted = true;
          }
     }

   // move non-deleted entries to the front of  interval
   for  (i = j = 0;  i < n;  i ++)
     if  (! interval [i] . deleted)
         {
          if  (i != j)
              interval [j] = interval [i];
          j ++;
         }

   interval . resize (j);

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



static void  Entropy_Filter
    (vector <Orf_t> & orf_list, double cutoff)

//  Remove from  orf_list  all entries whose entropy distance
//  is  >= cutoff .

  {
   int  i, j, n;

   n = orf_list . size ();
   for  (i = j = 0;  i < n;  i ++)
     {
      double  score;
      int  frame, len, start, stop;

      stop = orf_list [i] . Get_Stop_Position ();
      len = orf_list [i] . Get_Gene_Len ();
      frame = orf_list [i] . Get_Frame ();
      if  (frame > 0)
          start = On_Seq_1 (stop - len);
        else
          start = On_Seq_1 (stop + len + 2);
      score = Entropy_Distance_Ratio (start, len, frame);
      if  (score < cutoff)
          {
           if  (i != j)
               orf_list [j] = orf_list [i];
           j ++;
          }
     }

   orf_list . resize (j);

   return;
  }



static int  Find_Optimal_Len
    (const vector <Orf_Interval_t> & interval)

//  Find the length  L  such that considering only entries in
//   interval  L  or longer and eliminating from them any that
//  overlap others (by more than  Max_Olap_Bases ) either the sum of
//  the interval lengths or the number of intervals is maximum.
//  Return the optimal value of  L .

  {
   vector <Range_t>  range_list;
   Range_t  new_range;
   int  left_wrap;
     // max positions any interval extends left of zero
   int  right_wrap;
     // max positions any interval extends right of  Sequence_Len
   vector <int>  highest;
   priority_queue <int, vector <int>, greater <int> >  pq;
   int  opt_bases_len, opt_total_bases, total_bases;
     // these will determine the gene length that maximizes the sum
     // of non-overlapping gene lengths
   int  opt_count_len, opt_count, count;
     // these will determine the gene length that maximizes the number
     // of non-overlapping genes
   int  i, j, n;

   n = interval . size ();
   if  (n == 0)
       return  Min_Gene_Len;  // nothing to do; return the existing value

   highest . resize (n);
     // highest [i]  will be the max hi value in interval [0 .. i]

   // set values in  highest ,  left_wrap  and  right_wrap
   // first value determines  left_wrap  since  interval  is sorted
   // ascending by  lo  value
   if  (Genome_Is_Circular && interval [0] . lo < 0)
       left_wrap = -1 * interval [0] . lo;
     else
       left_wrap = 0;
   right_wrap = 0;
   for  (i = 0;  i < n;  i ++)
     {
      if  (i == 0)
          highest [i] = interval [i] . hi;
        else
          highest [i] = Max (highest [i - 1], interval [i] . hi);
      if  (Genome_Is_Circular)
          right_wrap = Max (right_wrap, interval [i] . hi - Sequence_Len);
     }

   for  (i = 0;  i < n;  i ++)
     {
      int  i_len, longest = Min_Gene_Len - 1;

      // first, for every entry, find the longest orf it overlaps
      // entries in  interval  must already be sorted
      // first look for overlaps on the left
      for  (j = i - 1;  0 <= j;  j --)
        {
         if  (highest [j] <= interval [i] . lo + Max_Olap_Bases)
             break;  // can't sufficiently overlap j or anything before it

         if  (Max_Olap_Bases < Intersect_Size (interval [j] . lo, interval [j] . hi,
                interval [i] . lo, interval [i] . hi))
             longest = Max (longest, interval [j] . hi - interval [j] . lo);
        }
      // also check any wraparounds
      if  (Genome_Is_Circular
             && interval [i] . lo + Max_Olap_Bases <= right_wrap)
          {
           for  (j = n - 1;  j > i && interval [i] . lo + Max_Olap_Bases
                      <= highest [j] - Sequence_Len;  j --)
             if  (Max_Olap_Bases < Intersect_Size (interval [i] . lo,
                    interval [i] . hi, interval [j] . lo - Sequence_Len,
                    interval [j] . hi - Sequence_Len))
                 longest = Max (longest, interval [j] . hi - interval [j] . lo);
          }

      // now look for overlaps on the right
      for  (j = i + 1;  j < n;  j ++)
        {
         if  (interval [i] . hi <= interval [j] . lo + Max_Olap_Bases)
             break;  // can't sufficiently overlap j or anything after it

         if  (Max_Olap_Bases < Intersect_Size (interval [j] . lo, interval [j] . hi,
                interval [i] . lo, interval [i] . hi))
             longest = Max (longest, interval [j] . hi - interval [j] . lo);
        }
      // check wraparounds
      if  (Genome_Is_Circular
             && Sequence_Len - interval [i] . hi + Max_Olap_Bases <= left_wrap)
          {
           for  (j = 0;  j < i && Sequence_Len + interval [j] . lo
                      <= interval [i] . hi - Max_Olap_Bases;  j ++)
             if  (Max_Olap_Bases < Intersect_Size (interval [i] . lo,
                    interval [i] . hi, interval [j] . lo + Sequence_Len,
                    interval [j] . hi + Sequence_Len))
                 longest = Max (longest, interval [j] . hi - interval [j] . lo);
          }

      i_len = interval [i] . hi - interval [i] . lo;

      if  (longest < i_len)
          { // in this case setting the min gene length to a value
            // from  longest + 1 .. i_len  inclusive will get the i_len
            // bases included.  Any smaller value and i_len will be
            // excluded because of it's overlap with longest, and any
            // larger value and i_len will be excluded because it's
            // too short itself.  We'll keep track of this in  range_list
           
           new_range . lo = longest + 1;
           new_range . hi = i_len;
           range_list . push_back (new_range);
          }
        else
          ;  // do nothing--can never include orf i's bases
     }

   // we can now get the optimum length from the entries in  range_list
   // first we sort by lo

   sort (range_list . begin (), range_list . end (), Range_Cmp);

   // For a minimum gene length of  m = range_list [i] . lo  the total
   // bases is the sum of entries  range_list [0 .. i] . hi  that
   // are >= m.  We compute these by scanning  range_list  entries in
   // order and using a priority queue to subtract out  hi  entries that
   // are too small.  As we're going, we keep track of the  lo  entry
   // that achieved the maximum total bases.  There is never any point
   // in choosing a minimum gene length that is not equal to a value
   // of  range_list [i] . lo  because any higher value could cause
   //  range_list [] . hi  values to drop out.

   n = range_list . size ();
   if  (n == 0)
       return  Min_Gene_Len;
          // nothing to do; return the existing value
          // can only happen if max overlap gene of every gene is the same
          // size as the gene itself--not very likely

   opt_bases_len = opt_total_bases = total_bases = 0;
   opt_count_len = opt_count = count = 0;
   
   for  (i = 0;  i < n;  i ++)
     {
      total_bases += range_list [i] . hi;
      count ++;
      while  (! pq . empty () && pq . top () < range_list [i] . lo)
        {
         total_bases -= pq . top ();
         count --;
         pq . pop ();
        }
      if  (opt_total_bases < total_bases
             || (opt_total_bases == total_bases && opt_count < count))
          {
           opt_total_bases = total_bases;
           opt_bases_len = range_list [i] . lo;
          }
      if  (opt_count < count
             || (opt_count == count && opt_total_bases < total_bases))
          {
           opt_count = count;
           opt_count_len = range_list [i] . lo;
          }
      pq . push (range_list [i] . hi);
     }

   if  (Optimize_Total_Len)
       return  Max (opt_bases_len, Min_Gene_Len);
     else
       return  Max (opt_count_len, Min_Gene_Len);
  }



static void  Find_Orfs
    (vector <Orf_t> & orf_list)

//  Put in  orf_list  all sufficiently long orfs in global
//  string  Sequence .

  {
   Orf_t  orf;
   Codon_t  codon;

   // Positions stored in these are the first (i.e., lowest-subscript)
   // base of the codon, using positions starting at 1.
   int  first_fwd_start [3] = {INT_MAX, INT_MAX, INT_MAX};
   int  last_rev_start [3] = {0};
   int  prev_fwd_stop [3] = {0}, prev_rev_stop [3] = {0};
   int  first_fwd_stop [3] = {0};
        // Used for wraparound in circular genomes
   int  ignore_start, ignore_stop;
        // indicate next beginning and ending positions of next
        // region to be ignored
   int  ignore_ct;
        // number of ignore regions
   int  ignore_sub;
        // subscript of current ignore region
   bool  hit_ignore = false;
        // indicates if any ignore region has been reached yet
   bool  ignoring = false;
        // indicates current status of ignore region
   int  first_base = 1;
        // position of the first base in the current region being
        // processed
   int  frame, gene_len, orf_len;
        // frame subscripts are 0, 1, 2 for both forward and reverse
        // events.  The frame is based on the *LAST* (i.e., highest-subscript)
        // base of the codon, using positions starting at 0
   int  i, j, n;

   orf_list . clear ();
   n = Sequence_Len;

   if  (n < Min_Gene_Len)
       return;

   if  (Genome_Is_Circular)
       {
        // allow 2-base overhang to catch start and stop codons that
        // span the end of  Sequence
        n += 2;
        Sequence . push_back (Sequence [0]);
        Sequence . push_back (Sequence [1]);
       }
 
   if  (Ignore_Region . size () == 0)
       ignore_start = ignore_stop = INT_MAX;
     else
       {
        ignore_ct = Ignore_Region . size ();
        ignore_start = Ignore_Region [0] . lo;
        ignore_stop = Ignore_Region [0] . hi;
        ignore_sub = 0;
       }

   frame = 0;
   for  (i = 0;  i < n;  i ++)
     {
      if  (i == ignore_start)
          {
           Finish_Orfs (false, prev_rev_stop, last_rev_start, i, orf_list);
           hit_ignore = ignoring = true;
          }
      else if  (i == ignore_stop)
          {
           // reset saved positions to their initial values as if the
           // start of the genome
           for  (j = 0;  j < 3;  j ++)
             {
              first_fwd_start [j] = INT_MAX;
              last_rev_start [j] = 0;
              prev_fwd_stop [j] = 0;
              prev_rev_stop [j] = 0;
             }
           codon . Clear ();
           first_base = i + 1;
           ignoring = false;
           ignore_sub ++;
           if  (ignore_sub >= ignore_ct)
               ignore_start = ignore_stop = INT_MAX;
             else
               {
                ignore_start = Ignore_Region [ignore_sub] . lo;
                ignore_stop = Ignore_Region [ignore_sub] . hi;
               }
          }

      if  (! ignoring)
          {
           int  which;

           codon . Shift_In (Sequence [i]);

           if  (codon . Can_Be (Fwd_Start_Pattern, which)
                   && first_fwd_start [frame] == INT_MAX)
               first_fwd_start [frame] = i - 1;

           if  (codon . Can_Be (Rev_Start_Pattern, which))
               {
                last_rev_start [frame] = i - 1;
               }

           if  (codon . Must_Be (Fwd_Stop_Pattern, which))
               {
                if  (prev_fwd_stop [frame] == 0)
                    {
                     Handle_First_Forward_Stop (frame, i - 1, first_fwd_start [frame],
                          first_base, gene_len, orf_len,
                          Genome_Is_Circular && ! hit_ignore);
                     first_fwd_stop [frame] = i - 1;
                    }
                  else
                    {
                     gene_len = i - first_fwd_start [frame] - 1;
                     orf_len = i - prev_fwd_stop [frame] - 4;
                    }

                if  (gene_len >= Min_Gene_Len)
                    {
                     orf . Set_Stop_Position (i - 1);
                     orf . Set_Frame (1 + (frame + 1) % 3);
                     orf . Set_Gene_Len (gene_len);
                     orf . Set_Orf_Len (orf_len);
                     orf_list . push_back (orf);
                    }

                first_fwd_start [frame] = INT_MAX;
                prev_fwd_stop [frame] = i - 1;
               }

           if  (codon . Must_Be (Rev_Stop_Pattern, which))
               {
                if  (prev_rev_stop [frame] != 0)
                    {
                     gene_len = last_rev_start [frame] - prev_rev_stop [frame];

                     if  (gene_len >= Min_Gene_Len)
                         {
                          orf . Set_Stop_Position (prev_rev_stop [frame]);
                          orf . Set_Frame (-1 - (frame + 1) % 3);
                          orf . Set_Gene_Len (gene_len);
                          orf . Set_Orf_Len (i - prev_rev_stop [frame] - 4);
                          orf_list . push_back (orf);
                         }
                    }
                last_rev_start [frame] = 0;
                prev_rev_stop [frame] = i - 1;
               }
          }

      if  (frame == 2)
          frame = 0;
        else
          frame ++;
     }

   Finish_Orfs (Genome_Is_Circular, prev_rev_stop, last_rev_start,
        Sequence_Len, orf_list);

   if  (Genome_Is_Circular)
       Sequence . resize (Sequence_Len);

   return;
  }



static void  Finish_Orfs
    (bool use_wraparound, const int prev_rev_stop [3],
     const int last_rev_start [3], int last_position,
     vector <Orf_t> & orf_list)

//  Finish reverse-strand orfs because we've hit the end of the
//  genome (or hit an ignore region).  If  use-wraparound  is true,
//  then the orfs can wrap around the end of the (circular) genome;
//  otherwise, not.   prev_rev_stop  has the position of the last-seen
//  reverse stop codons in each frame, and  last_rev_start  has the
//  position of the last-seen reverse start codons in each frame.
//   last_position  is the last available sequence position to use.
//  Add any suitable orfs to  orf_list .

  {
   Orf_t  orf;
   int  frame, gene_len, orf_len;

   for  (frame = 0;  frame < 3;  frame ++)
     {
      Handle_Last_Reverse_Stop (frame, prev_rev_stop, last_rev_start,
           gene_len, orf_len, use_wraparound, last_position);
      if  (gene_len >= Min_Gene_Len)
          {
           orf . Set_Stop_Position (prev_rev_stop [frame]);
           orf . Set_Frame (-1 - (frame + 1) % 3);
           orf . Set_Gene_Len (gene_len);
           orf . Set_Orf_Len (orf_len);
           orf_list . push_back (orf);
          }
     }

   return;
  }



static void  Get_Ignore_Regions
    (void)

//  Read the list of regions from the with name in global
//   Ignore_File_Name .  Sort them and coalesce overlapping regions.
//  Put the results in global  Ignore_Region .  The format for each
//  line of input is:
//     <lo>  <hi>  <rest of line ignored>  
//  where <lo> and <hi> are integer values.  The region specified
//  is bases <lo>..<hi> inclusive, where bases are numbered starting
//  at 1.  If <hi> is less than <lo> the values are silently swapped.
//  There is no provision for circularity.  If more than one sequence
//  is read in to be searched for genes, these regions will be used
//  to screen them *ALL*, which is very likely not at all what is
//  desired.  Blank lines and lines beginning with # are skipped.

  {
   FILE  * fp;
   char  line [MAX_LINE];
   Range_t  range;
   int  i, j, n, line_ct;

   fp = File_Open (Ignore_File_Name, "r", __FILE__, __LINE__);

   line_ct = 0;
   while  (fgets (line, MAX_LINE, fp) != NULL)
     {
      char  * p;
      int  a, b;

      line_ct ++;

      // set  p  to point to the first non-blank character on the line
      for  (p = line;  * p != '\0' && isspace (* p);  p ++)
        ;
      
      if  (* p == '\0' || * p == '#')
          continue;
      else if  (sscanf (line, "%d %d", & a, & b) == 2)
          {
           if  (a < b)
               {
                range . lo = a - 1;
                  // convert to 0-based between coordinates
                range . hi = b;
               }
             else
               {
                range . lo = b - 1;
                range . hi = a;
               }
           Ignore_Region . push_back (range);
          }
        else
          {
           fprintf (stderr, "ERROR:  Following line %d in file %s is bad--skipped:\n",
                line_ct, Ignore_File_Name);
           fputs (line, stderr);
           fputc ('\n', stderr);
          }
     }

   fclose (fp);

   // sort regions by lo value
   sort (Ignore_Region . begin (), Ignore_Region . end (), Range_Cmp);

   // combine overlapping regions and move them to the front of  Ignore_Region
   n = Ignore_Region . size ();

   if  (n <= 1)
       return;

   for  (i = 0, j = 1;  j < n;  j ++)
     if  (Ignore_Region [j] . lo < Ignore_Region [i] . hi)
         {  // overlap
          if  (Ignore_Region [i] . hi < Ignore_Region [j] . hi)
              Ignore_Region [i] . hi = Ignore_Region [j] . hi;
                 // j extends i to the right
         }
       else
         {
          i ++;
          if  (i != j)
              Ignore_Region [i] = Ignore_Region [j];
                // move j region down to front of list
         }

   Ignore_Region . resize (i + 1);

   return;
  }



static void  Get_Intervals
    (vector <Orf_Interval_t> & interval, const vector <Orf_t> & orf_list)

//  Populate  interval  with intervals corresponding to the entries
//  in  orf_list .  Intervals are in 0-based between coordinates and
//  are sorted by lo value then by hi.

  {
   Orf_Interval_t  new_int;
   int  i, n;

   interval . clear ();

   new_int . deleted = false;

   n = orf_list . size ();
   for  (i = 0;  i < n;  i ++)
     {
      int  frame, stop, len;

      frame = orf_list [i] . Get_Frame ();
      stop = orf_list [i] . Get_Stop_Position ();
      len = orf_list [i] . Get_Gene_Len ();
         // does not include the stop codon length
      if  (frame > 0)
          {
           new_int . hi = On_Seq_0 (stop - 1);
           new_int . lo = new_int . hi - len;
          }
        else
          {
           new_int . lo = On_Seq_0 (stop + 2);
           new_int . hi = new_int . lo + len;
          }
        //  new_int . lo  and  ol . hi  are the 0-based between coordinates
        //  of the coding portion of the gene
      new_int . frame = frame;
        // keep track of corresponding entry in  orf_list  so can delete
      interval . push_back (new_int);
     }

   sort (interval . begin (), interval . end (), Orf_Interval_Cmp);

   return;
  }



static void  Handle_First_Forward_Stop
     (int fr, int pos, int start_pos, int first_base, int & gene_len,
      int & orf_len, bool use_wraparound)

//  Handle the case of a forward stop codon, beginning at position
//   pos  in the global  Sequence  (counting starting at 1)  which
//  is in frame subscript  fr  (0, 1 or 2).   start_pos  is the
//  position of the first possible start codon in this frame, or else
//   INT_MAX  if none has been encountered yet.   first_base  is the
//  position of the first base in this region.  Set gene_len
//  to the length of longest possible gene for this orf.  If no gene
//  is possible (e.g., because there is no start codon), then set
//   gene_len  to  0 .  Set  orf_len  to the length of this orf.
//  If  use_wraparound  is true, allow orfs/genes to wrap around
//  through the front of the (circular) sequence.

  {
   if  (use_wraparound)
       {
        Wrap_Through_Front (fr, pos, gene_len, orf_len);
        if  (gene_len == 0 && start_pos != INT_MAX)
            gene_len = pos - start_pos;
       }
     else
       {
        // assume the orf is entirely contained in  Sequence  no
        // matter whether the odd 1 or 2 bases at the front could be
        // a stop or not
        orf_len = pos - first_base;
        orf_len -= orf_len % 3;  // round down
        if  (start_pos == INT_MAX)
            gene_len = 0;
          else
            gene_len = pos - start_pos;
       }

   return;
  }



static void  Handle_Last_Reverse_Stop
     (int fr, const int prev_rev_stop [3], const int last_rev_start [3],
      int & gene_len, int & orf_len, bool use_wraparound, int last_position)

//  Set  orf_len  and  gene_len  to the length of the last orf, and longest
//  gene in it, resp., in reverse reading frame  fr .
//   prev_rev_stop  has the last stop position in  Sequence  in each
//  reverse reading frame, and  last_rev_start  has the corresponding
//  last start locations.    use_wraparound  indicates whether the
//  orfs are allowed to wrap around the end of the (circular) genome.
//   last_position  is the highest-numbered sequence position available

  {
   if  (prev_rev_stop [fr] == 0)
       {
        // no reverse stop in this frame at all
        gene_len = orf_len = 0;
        return;
       }

   if  (use_wraparound)
       {
        int  wrap_fr;
             // the frame at the front of the genome corresponding
             // to  fr
        wrap_fr = (3 + fr - (Sequence_Len % 3)) % 3;

        Wrap_Around_Back (wrap_fr, prev_rev_stop [fr], gene_len, orf_len);

        if  (gene_len == 0 && last_rev_start [fr] > 0)
            gene_len = last_rev_start [fr] - prev_rev_stop [fr];
       }
     else
       {
        orf_len = last_position - prev_rev_stop [fr] - 2;
             // round down to next multiple of 3
        orf_len -= orf_len % 3;

        if  (last_rev_start [fr] == 0)
            gene_len = 0;
          else
            gene_len = last_rev_start [fr] - prev_rev_stop [fr];
       }


   assert (orf_len % 3 == 0);
   assert (gene_len % 3 == 0);

   return;
  }



static int  Intersect_Size
    (int a, int b, int c, int d)

//  Return the number of bases by which region  a .. b  overlaps
//  region  c .. d .  All values are space-based coordinates.

  {
   if  (d <= a || b <= c)
       return  0;

   return  Min (b, d) - Max (a, c);
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



static int  On_Seq_1
    (int i)

//  Return the subscript equivalent to  i  on a sequence of
//  length  Sequence_Len  (with subscripts starting at 1)
//  assuming circular wraparounds.

  {
   while  (i < 1)
     i += Sequence_Len;
   while  (Sequence_Len < i)
     i -= Sequence_Len;

   return  i;
  }



static void  Output_Orfs
    (FILE * fp, const vector <Orf_Interval_t> & interval, int & total_len)

//  Print the regions in  interval  to  fp  and set  total_len  to
//  the sum of their lengths.  Include in the output the frame and
//  entropy distance value of each region.

  {
   double  entropy_ratio;
   int  start, stop, len;
   int  i, n;

   if  (Print_Output_Header)
       fprintf (fp, "\nPutative Genes:\n");

   total_len = 0;

   n = interval . size ();
   for  (i = 0;  i < n;  i ++)
     {
      len = interval [i] . hi - interval [i] . lo;
      total_len += len;
         // does not include the stop codon length
      if  (interval [i] . frame > 0)
          {
           if  (Without_Stops)
               {
                stop = On_Seq_1 (interval [i] . hi);
                start = On_Seq_1 (stop - len + 1);
               }
             else
               {
                stop = On_Seq_1 (interval [i] . hi + 3);
                start = On_Seq_1 (stop - len - 2);
               }
          }
        else
          {
           if  (Without_Stops)
               {
                stop = On_Seq_1 (interval [i] . lo + 1);
                start = On_Seq_1 (stop + len - 1);
               }
             else
               {
                stop = On_Seq_1 (interval [i] . lo - 2);
                start = On_Seq_1 (stop + len + 2);
               }
          }
      entropy_ratio = Entropy_Distance_Ratio (start, len, interval [i] . frame);
      fprintf (fp, "%05d %7d %7d  %+2d  %6.3f\n", i + 1, start, stop,
           interval [i] . frame, entropy_ratio);
     }

   return;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   char  * p, * q;
   bool  errflg = false;
   int  ch;

   optarg = NULL;

#if  ALLOW_LONG_OPTIONS
   int  option_index = 0;
   static struct option  long_options [] = {
        {"start_codons", 1, 0, 'A'},
        {"entropy", 1, 0, 'E'},
        {"fixed", 0, 0, 'f'},
        {"min_len", 1, 0, 'g'},
        {"help", 0, 0, 'h'},
        {"ignore", 1, 0, 'i'},
        {"linear", 0, 0, 'l'},
        {"length_opt", 0, 0, 'L'},
        {"no_header", 0, 0, 'n'},
        {"max_olap", 1, 0, 'o'},
        {"cutoff", 1, 0, 't'},
        {"without_stops", 0, 0, 'w'},
        {"trans_table", 1, 0, 'z'},
        {"stop_codons", 1, 0, 'Z'},
        {0, 0, 0, 0}
      };

   while  (! errflg && ((ch = getopt_long (argc, argv,
        "A:E:fg:hi:lno:t:wz:Z:",
        long_options, & option_index)) != EOF))
#else
   while  (! errflg && ((ch = getopt (argc, argv,
        "A:E:fg:hi:lno:t:wz:Z:")) != EOF))
#endif

     switch  (ch)
       {
        case  'A' :
          Start_Codon . clear ();
          for  (p = strtok (optarg, ",");  p != NULL;  p = strtok (NULL, ","))
            {
             q = strdup (p);
             Make_Lower_Case (q);
             Start_Codon . push_back (q);
            }
          break;

        case  'E' :
          Read_Entropy_Profiles (optarg, errflg);
          break;

        case  'f' :
          Fixed_Min_Len = true;
          break;

        case  'g' :
          Min_Gene_Len = strtol (optarg, & p, 10);
          if  (p == optarg || Min_Gene_Len <= 0)
              {
               fprintf (stderr, "ERROR:  Bad minimum gene length (-g option)\n"
                    "  value = \"%s\"", optarg);
               errflg = true;
              }
          break;

        case  'h' :
          errflg = true;
          break;

        case  'i' :
          Ignore_File_Name = optarg;
          break;

        case  'l' :
          Genome_Is_Circular = false;
          break;

        case  'L' :
          Optimize_Total_Len = true;
          break;

        case  'n' :
          Print_Output_Header = false;
          break;

        case  'o' :
          Max_Olap_Bases = strtol (optarg, & p, 10);
          if  (p == optarg || Max_Olap_Bases < 0)
              {
               fprintf (stderr, "ERROR:  Bad max overlap bases (-o option)\n"
                    "  value = \"%s\"", optarg);
               errflg = true;
              }
          break;

        case  't' :
          Entropy_Cutoff = strtod (optarg, & p);
          Use_Entropy_Filter = true;
          break;

        case  'w' :
          Without_Stops = true;
          break;

        case  'z' :
          Genbank_Xlate_Code = strtol (optarg, & p, 10);
          Set_Stop_Codons_By_Code (Stop_Codon, Genbank_Xlate_Code, errflg);
          break;

        case  'Z' :
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

   if  (optind > argc - 2)
       {
        Usage ();
        exit (EXIT_FAILURE);
       }

   Sequence_File_Name = argv [optind ++];
   Output_Filename = argv [optind ++];

   return;
  }



static void  Print_Comma_Separated_Strings
    (const vector <const char *> & v, FILE * fp)

//  Print the strings in  v  to  fp .  Separate them by
//  commas with no spaces.

  {
   int  i, n;

   n = v . size ();

   if  (n == 0)
       return;

   fprintf (fp, "%s", v [0]);
   for  (i = 1;  i < n;  i ++)
     fprintf (fp, ",%s", v [i]);

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



static void  Remove_Shorter
    (vector <Orf_Interval_t> & interval, int len)

//  Remove from  interval  any entry shorter than  len .

  {
   int  i, j, n;

   n = interval . size ();
   for  (i = j = 0;  i < n;  i ++)
     if  (len <= interval [i] . hi - interval [i] . lo)
         {
          if  (i != j)
              interval [j] = interval [i];
          j ++;
         }

   interval . resize (j);

   return;
  }



static void  Set_Start_And_Stop_Codons
    (void)

//  Set globals  Start_Codon  and  Stop_Codon  to the sequences
//  that are allowed to be start and stop codons for genes.

  {
   Codon_t  codon;
   int  i, n;

   if  (Start_Codon . size () == 0)
       {
        n = sizeof (DEFAULT_START_CODON) / sizeof (char *);
        for  (i = 0;  i < n;  i ++)
          Start_Codon . push_back (DEFAULT_START_CODON [i]);
       }

   if  (Stop_Codon . size () == 0)
       {
        n = sizeof (DEFAULT_STOP_CODON) / sizeof (char *);
        for  (i = 0;  i < n;  i ++)
          Stop_Codon . push_back (DEFAULT_STOP_CODON [i]);
       }

   Fwd_Start_Pattern . clear ();
   Fwd_Stop_Pattern . clear ();
   Rev_Start_Pattern . clear ();
   Rev_Stop_Pattern . clear ();

   n = Start_Codon . size ();
   for  (i = 0;  i < n;  i ++)
     {
      codon . Set_From (Start_Codon [i]);
      Fwd_Start_Pattern . push_back (codon);
      codon . Reverse_Complement ();
      Rev_Start_Pattern . push_back (codon);
     }

   n = Stop_Codon . size ();
   for  (i = 0;  i < n;  i ++)
     {
      codon . Set_From (Stop_Codon [i]);
      Fwd_Stop_Pattern . push_back (codon);
      codon . Reverse_Complement ();
      Rev_Stop_Pattern . push_back (codon);
     }

   return;
  }



static void  Usage
    (void)

//  Print to stderr description of options and command line for
//  this program.

  {
   fprintf (stderr,
       "USAGE:  long-orfs [options] <sequence-file> <output-file>\n"
       "\n"
       "Read DNA sequence in <sequence-file> and find and output the\n"
       "coordinates of long, non-overlapping orfs in it.\n"
       "Output goes to file <output-file> or standard output if <output-file>\n"
       "is \"-\"\n"
       "\n"
       "Options:\n"
       " -A <codon-list>\n"
       " --start_codons <codon-list>\n"
       "    Use comma-separated list of codons as start codons\n"
       "    Sample format:  -A atg,gtg\n"
       " -E <filename>\n"
       " --entropy <filename>\n"
       "    Read entropy profiles from <filename>.  Format is one header\n"
       "    line, then 20 lines of 3 columns each.  Columns are amino acid,\n"
       "    positive entropy, negative entropy.  Rows must be in order\n"
       "    by amino acid code letter\n"
       " -f\n"
       " --fixed\n"
       "    Do *NOT* automatically determine the minimum gene length so as\n"
       "    to maximize the total length of output regions\n"
       " -g <n>\n"
       " --min_len <n>\n"
       "    Only genes with length >= <n> will be considered\n"
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
       " -L\n"
       " --length_opt\n"
       "    Find and use the minimum gene length that maximizes the total\n"
       "    length of non-overlapping genes, instead of maximizing the\n"
       "    number of such genes\n"
       " -n\n"
       " --no_header\n"
       "    Do not include heading information in the output; only output\n"
       "    the orf-coordinate lines\n"
       " -o <n>\n"
       " --max_olap <n>\n"
       "    Set maximum overlap length to <n>.  Overlaps this short or shorter\n"
       "    are ignored.\n"
       " -t <x>\n"
       " --cutoff <x>\n"
       "    Only genes with entropy distance score less than <x> will be considered\n"
       " -w\n"
       " --without_stops\n"
       "    Do *NOT* include the stop codon in the output coordinates.\n"
       "    By default, it is included.\n"
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



static void  Wrap_Around_Back
    (int wfr, int pos, int & gene_len, int & orf_len)

//  Set  orf_len  to the length of the complement-strand orf that
//  wraps around the end of the sequence in global  Sequence .  The
//  stop codon for the orf is at position  pos  (first base of codon
//  numbered starting at 1).   wfr  is the frame subscript of the
//  reading frame to use at the beginning of  Sequence  (i.e., it
//  allows for  Sequence_Len  not being a multiple of 3).  The
//  maximum possible orf length is  Sequence_Len - 3  rounded down
//  to the nearest multiple of 3.  Set  gene_len  to the longest
//  possible gene in that orf, looking only for starts that are completely
//  contained in the start of  Sequence .  If no starts are found,
//  set  gene_len  to  0  (even though there may be starts between
//   pos  and the end of  Sequence ).

  {
   Codon_t  codon;
   int  start_at, check_len, frame, orf_add, which;
   int  i;

   assert (pos > 0);
   check_len = pos - 1;

   start_at = -1;
   orf_add = 0;
     // this is the number of extra bases at the front of the sequence
     // to add to the orf at the back
   frame = 0;
   for  (i = 0;  i < check_len;  i ++)
     {
      codon . Shift_In (Sequence [i]);

      if  (frame == wfr)
          {
           if  (codon . Must_Be (Rev_Stop_Pattern, which))
               {
                orf_add = i - 2;
                break;
               }
             else
               orf_add = i + 1;
          }
      if  (frame == wfr && codon . Can_Be (Rev_Start_Pattern, which))
          start_at = i + 1;

      if  (frame == 2)
          frame = 0;
        else
          frame ++;
     }

   orf_len = orf_add + Sequence_Len - pos - 2;
   orf_len -= orf_len % 3;
   if  (start_at == -1)
       gene_len = 0;
     else
       gene_len = start_at + Sequence_Len - pos - 2;
   
   return;
  }



static void  Wrap_Through_Front
    (int fr, int pos, int & gene_len, int & orf_len)

//  Set  orf_len  to the length of the orf with forward frame subscript
//   fr  with stop codon at position  pos  that wraps around and begins
//  at the end of the sequence in global  Sequence .  Set  gene_len
//  to the longest possible gene in that orf.  Start looking at the
//  beginning of  Sequence  and assume there are no stops between
//  there and  pos .  If no starts are found, set  gene_len  to  0
//  (even though there may be starts between  0  and  pos in  Sequence ).

  {
   Codon_t  codon;
   int  start_at, check_len, which;
   int  i, j, s;

   assert (pos > 0);
   start_at = -1;
   s = (pos - 1) % 3;
   check_len = Sequence_Len + s - pos - 4;

   // Loop back to at most original stop codon.  Do not allow the
   // orf to overlap that stop codon.
   for  (i = 0;  i < check_len;  i += 3)
     {
      for  (j = 0;  j < 3;  j ++)
        {
         s --;
         if  (s < 0)
             s += Sequence_Len;
         codon . Reverse_Shift_In (Sequence [s]);
        }

      if  (codon . Must_Be (Fwd_Stop_Pattern, which))
          break;
      if  (codon . Can_Be (Fwd_Start_Pattern, which))
          start_at = i + 3;

     }

   orf_len = i + 3 * ((pos - 1) / 3);
   if  (start_at == -1)
       gene_len = 0;
     else
       gene_len = start_at + 3 * ((pos - 1) / 3);
   
   return;
  }



