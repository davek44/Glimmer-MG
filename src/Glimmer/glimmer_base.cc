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

bool Dave_Log = false;
bool Length_Log = false;
bool Detail_Log = false;
bool Sequence_Log = false;

// Glimmer3 specific
char  * Ignore_File_Name = NULL;
  // Name of file containing list of regions that cannot be included
  // in gene predictions
vector <Range_t>  Ignore_Region;
vector <double>  Start_Prob;
  // Probability of occurrence of start codons
bool  Use_Independent_Score = DEFAULT_USE_INDEPENDENT_SCORE;
  // If true, let the non-Markov independent model compete with
  // the periodic Markov models to score genes.


// Glimmer-MG specific
vector< pair<double,int> > Meta_PWM_Save;
  // Saved PWM scores so we don't recompute
vector<PWM_t>  Meta_Ribosome_PWMs;
  // Mixture of raw probability RBS PWM models
int Min_Indel_ORF_Len = 15;
  // Minimum ORF size to allow to be scored and considered for indels

void  Add_Events_Fwd
    (const Orf_t & orf, vector <Start_t> & start_list, int & id)

//  Add events for  orf  with possible coding start sites in
//   start_list  to global  Last_Event .   id  is the id number
//  of this orf (which corresponds to the numbers in the detail
//  file.

{
     Event_Node_t  * ne;   // new event
     int  fr, sub;
     int  i, n;
     Event_Node_t rbs_dummy, rbs_dummy_max;
     map< vector<Error_t>, int, vec_error_cmp> error_id_map;
     map< vector<Error_t>, int, vec_error_cmp>::const_iterator error_it;
     map< int, Event_Node_t*> start_event_map;
     map< int, Event_Node_t*>::iterator start_it;

     fr = orf . Get_Frame ();
     n = start_list . size ();
     sub = fr - 1;

     for  (i = 0;  i < n;  i ++)
     {
	  if  (1 + start_list [i] . j >= Min_Gene_Len)
	  {
	       ne = new Event_Node_t;
	       ne -> e_type = FWD_START;
	       ne -> pos = start_list [i] . pos + 2;
	       // event pos is last base of codon; start pos is first
	       ne -> frame = fr;
	       ne -> score = start_list [i] . score + LogOdds_Prior;

	       if(User_RBS)
		    PWM_Score_Fwd_Start (start_list [i] . pos, LogOdds_PWM,
					 Ribosome_Window_Size, ne -> pwm_score, ne -> pwm_sep);
	       else
		    PWM_Meta_Score_Fwd_Start (start_list [i] . pos, ne -> pwm_score, ne -> pwm_sep);
	       Add_PWM_Score (ne);

	       if  (start_list [i] . which >= 0)
		    ne -> score += LogOdds_Start.Score(start_list [i] . which);

	       ne -> score += LogOdds_Length.Score((unsigned int)((1+start_list[i].j)/3),
						   start_list[i].truncated, orf.Get_Stop_Position() > Sequence_Len-2,
						   Sequence_Len/3);
	       ne -> is_first_start = start_list [i] . first;
	       ne -> truncated = start_list [i] . truncated;
	       ne -> best_pred = NULL;
	       ne -> frame_pred = Last_Event [sub];
	       ne -> errors = start_list[i].errors;

	       if(Dave_Log)
		    Print_Start(start_list[i], orf.Get_Stop_Position(), ne);

	       if(ne->score > Event_Threshold) {
		    // hash start if it's the best at it's pos
		    start_it = start_event_map.find(ne->pos);
		    if(start_it == start_event_map.end())
			 start_event_map[ne->pos] = ne;
		    else if(ne->score > start_event_map[ne->pos]->score) {
			 //delete start_event_map[ne->pos];
			 //start_event_map[ne->pos] = ne;
			 delete start_it->second;
			 start_it->second = ne;
		    } else
			 delete ne;
	       } else {
		    delete ne;
	       }
	  }
     }

     // add best starts
     for(start_it = start_event_map.begin(); start_it != start_event_map.end(); start_it++) {
	  ne = start_it->second;

	  ne->frame_pred = Last_Event[sub];
	  Last_Event[sub] = ne;

	  // hash errors
	  error_it = error_id_map.find(ne->errors);
	  if(error_it == error_id_map.end())
	       error_id_map[ne->errors] = ++id;		    
	  ne->id = error_id_map[ne->errors];
     }

     // add stop codons, one for each error set
     if(!start_event_map.empty()) {
	  for(error_it = error_id_map.begin(); error_it != error_id_map.end(); error_it++) {
	       ne = new Event_Node_t;
	       ne -> e_type = FWD_STOP;
	       ne -> id = error_it->second;
	       ne -> pos = orf . Get_Stop_Position () + 2;
	       // event pos is last base of codon; orf pos is first
	       ne -> frame = fr;
	       ne -> is_first_start = false;
	       ne -> truncated = false;
	       ne -> score = 0.0;
	       ne -> best_pred = NULL;
	       ne -> frame_pred = Last_Event [sub];
	       ne -> errors = error_it->first;
	       Last_Event [sub] = ne;
	  }
     }

     return;
}


void  Add_Events_Rev
    (const Orf_t & orf, vector <Start_t> & start_list, int & id)

//  Add events for  orf  with possible coding start sites in
//   start_list  to global  Last_Event .   id  is the id number
//  of this orf (which corresponds to the numbers in the detail
//  file.

{
     Event_Node_t  * ne;   // new event
     int  fr, sub;
     int  i, n;
     Event_Node_t rbs_dummy, rbs_dummy_max;
     map< vector<Error_t>, int, vec_error_cmp> error_id_map;
     map< vector<Error_t>, int, vec_error_cmp>::const_iterator error_it;
     map< int, Event_Node_t*> start_event_map;
     map< int, Event_Node_t*>::iterator start_it;

     fr = orf . Get_Frame ();
     n = start_list . size ();
     sub = 2 - fr;

     for  (i = 0;  i < n;  i ++)
     {
	  if  (1 + start_list [i] . j >= Min_Gene_Len)
	  {
	       ne = new Event_Node_t;
	       ne -> e_type = REV_START;
	       //ne -> id = error_id_map[start_list[i] . errors];
	       ne -> pos = start_list [i] . pos;
	       // both pos's are last base of codon, i.e., highest coord
	       ne -> frame = fr;
	       ne -> score = start_list [i] . score + LogOdds_Prior;
		  
	       if(User_RBS)
		    PWM_Score_Rev_Start (start_list [i] . pos, LogOdds_PWM,
					 Ribosome_Window_Size, ne -> pwm_score, ne -> pwm_sep);
	       else
		    PWM_Meta_Score_Rev_Start (start_list [i] . pos, ne -> pwm_score, ne -> pwm_sep);
	       Add_PWM_Score (ne);
		  
	       if  (start_list [i] . which >= 0)
		    ne -> score += LogOdds_Start.Score(start_list [i] . which);

	       ne -> score += LogOdds_Length.Score((unsigned int)((1+start_list[i].j)/3),
						   start_list[i].truncated, orf.Get_Stop_Position() < 1,
						   Sequence_Len/3);
	       ne -> is_first_start = start_list [i] . first;
	       ne -> truncated = start_list [i] . truncated;
	       ne -> best_pred = NULL;
	       ne -> frame_pred = Last_Event [sub];
	       ne -> errors = start_list[i].errors;
		  
	       if(Dave_Log)
		    Print_Start(start_list[i], orf.Get_Stop_Position(), ne);

	       // if good enough, add event
	       if(ne->score > Event_Threshold) {
		    // hash start if it's the best at it's pos
		    start_it = start_event_map.find(ne->pos);
		    if(start_it == start_event_map.end())
			 start_event_map[ne->pos] = ne;
		    else if(ne->score > start_event_map[ne->pos]->score) {
			 //delete start_event_map[ne->pos];
			 //start_event_map[ne->pos] = ne;
			 delete start_it->second;
			 start_it->second = ne;
		    } else
			 delete ne;
	       } else {
		    delete ne;
	       }
	  }	       
     }

     // hash errors from best starts
     for(start_it = start_event_map.begin(); start_it != start_event_map.end(); start_it++) {
	  ne = start_it->second;
	  error_it = error_id_map.find(ne->errors);
	  if(error_it == error_id_map.end())
	       error_id_map[ne->errors] = ++id;		    
	  ne->id = error_id_map[ne->errors];
     }
     	  
     // make stop events     
     for(error_it = error_id_map.begin(); error_it != error_id_map.end(); error_it++) {
	  ne = new Event_Node_t;
	  ne -> e_type = REV_STOP;
	  ne -> id = error_it->second;
	  ne -> pos = orf . Get_Stop_Position () + 2;
	  // event pos is last base of codon; orf pos is first
	  ne -> frame = fr;
	  ne -> is_first_start = false;
	  ne -> truncated = false;
	  ne -> score = 0.0;
	  ne -> best_pred = NULL;
	  ne -> frame_pred = Last_Event [sub];
	  ne -> errors = error_it->first;
	  Last_Event [sub] = ne;
     }

     // add best starts
     for(start_it = start_event_map.begin(); start_it != start_event_map.end(); start_it++) {
	  ne = start_it->second;

	  ne->frame_pred = Last_Event[sub];
	  Last_Event[sub] = ne;
     }

     return;
}
   


void  Add_PWM_Score
    (Event_Node_t * p)

//  Add all or part of  p -> pwm_score  to  p -> score  depending
//  on the location of the PWM match.

  {
   static const int  LO_SEP = 4, HI_SEP = 10, HI_TAIL = 6;
   double  coeff;

   if  (p -> pwm_score < 0.0)
       return;

   // Use all the pwm_score if the pwm_sep is between LO_SEP and HI_SEP
   // Otherwise, use a fraction of it.
   if  (p -> pwm_sep < LO_SEP)
       coeff = double (p -> pwm_sep) / LO_SEP;
   else if  (p -> pwm_sep <= HI_SEP)
       coeff = 1.0;
   else if  (p -> pwm_sep < HI_SEP + HI_TAIL)
       coeff = double (HI_SEP + HI_TAIL - p -> pwm_sep) / HI_TAIL;
     else
       coeff = 0.0;

   if  (0.0 < coeff)
       p -> score += coeff * p -> pwm_score;

   return;
  }


void AdjDist_Smooth(vector<float> & dists)

// Perform kernel smoothing on an adjacent distance
// distribution in a special way.

{
     const float olap_sigma = 20;
     const float pos_sigma = 30;
     unsigned int d;

     // Regress overlapping counts 3-periodically
     // Note that it doesn't really matter which frame is which
     vector<float> overlap[3];
     for(d = 0; d < (unsigned int)Dist_Max_Overlap-5; d++)
	  overlap[d % 3].push_back(dists[d]);

     for(int i = 0; i < 3; i++) {
	  kernel_smooth(overlap[i], olap_sigma);
	  for(d = i; d < (unsigned int)Dist_Max_Overlap-5; d += 3)
	       dists[d] = overlap[i][(d-i)/3];
     }

     // Don't regress start/stop codon overlaps

     // Regress dists >= 0 as normally
     vector<float> pos;
     for(d = Dist_Max_Overlap; d < dists.size(); d++)
	  pos.push_back(dists[d]);

     kernel_smooth(pos, pos_sigma);

     for(d = Dist_Max_Overlap; d < dists.size(); d++)
	  dists[d] = pos[d-Dist_Max_Overlap];     
}


void Blend_Length(vector<double> & length_dist, vector<double> & par_length_dist, vector<double> & nonpar_length_dist, double par_cumprob)

// Blend the nonparametric and parametric length distributions linearly up to the point
// where 'par_cumprob' of the nonparametric distribution remains.

{
     unsigned int l;
     unsigned int min_aa_len = (unsigned int)ceil((float)Min_Gene_Len/3.0);

     // find the lower blending point
     double tmp_cumprob = 0;
     unsigned int blend_lower = min_aa_len;
     while(blend_lower < nonpar_length_dist.size() && tmp_cumprob < par_cumprob) {
	  tmp_cumprob += exp(nonpar_length_dist[blend_lower]);
	  blend_lower++;
     }

     // find the upper blending point
     tmp_cumprob = 0;
     unsigned int blend_upper = nonpar_length_dist.size()-1;
     while(blend_upper > 0 && tmp_cumprob < par_cumprob) {
	  tmp_cumprob += exp(nonpar_length_dist[blend_upper]);
	  blend_upper--;
     }

     if(blend_lower == nonpar_length_dist.size() || blend_upper == 0) {
	  sprintf (Clean_Exit_Msg_Line,
		   "ERROR:  Could not find quartiles of the nonparametric length distribution\n");
	  Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
     }

     // nonparametric section
     for(l = min_aa_len; l < blend_lower; l++)
	 length_dist[l] = nonpar_length_dist[l]; 

     // blend
     double coeff;
     double blend_dist = (double)blend_upper - (double)blend_lower;
     for(; l <= blend_upper; l++) {
	  //cerr << l << "\t" << par_length_dist[l] << "\t" << nonpar_length_dist[l] << "\t" << coeff << endl;
	  coeff = ((double)l - (double)blend_lower) / blend_dist;
	  length_dist[l] = coeff_log_add(par_length_dist[l], nonpar_length_dist[l], coeff);
     }

     // parametric
     for(; l < length_dist.size(); l++)
	  length_dist[l] = par_length_dist[l];

     // re-normalize
     log_normalize(length_dist, min_aa_len);
}


void  Clear_Events
    (void)

//  Free memory in chains pointed to by  Last_Event .  Note that
//  the initial event is not dynamically allocated (it's the global
//  variable  First_Event ) so it is not cleared.

  {
   Event_Node_t  * p, * q;
   int  i;

   for  (i = 0;  i < 6;  i ++)
     for  (p = Last_Event [i];  p != NULL && p -> e_type != INITIAL;  p = q)
       {
        q = p -> frame_pred;
        delete p;
       }

   return;
  }



void  Complement_Transfer
    (string & buff, const string & s, int start, int len)

//  Copy to string  buff  the substring of  s  starting at subscript
//   start  and going to the right for a length of  len .  Wraparound
//  the end of  s  if necessary.  Convert each character to its
//  Watson-Crick complement as it is copied.

  {
   int  j, n;

   n = s . length ();
   assert (start < n);
   assert (0 <= len);

   buff . resize (len);
   for  (j = 0;  j < len;  j ++, start ++)
     {
      if  (start >= n)
          start -= n;
      buff [j] = Complement (s [start]);
     }

   return;
  }


void  Disqualify
    (Event_Node_t * p, int cutoff)

//  Set the  disqualified  bit true for nodes reachable from
//   p  by  best_pred  pointers that have  pos  values at least
//  as great as  cutoff .

// DK: I think this might be to avoid connecting to anything that
// would already be in the highest scoring parse anyway.

{
     Event_Node_t  * q;
     
     if  (p == NULL)
	  return;
     
     for  (q = p -> best_pred;  q != NULL && cutoff <= q -> pos;  q = q -> best_pred)
	  q -> disqualified = true;
     
     return;
}



void  Do_Fwd_Stop_Codon
    (int i, int frame, int prev_fwd_stop [3], int first_fwd_start [3],
     int first_fwd_stop [3], int first_base, bool hit_ignore,
     vector <Orf_t> & orf_list)

//  Create a new entry for the forward orf ending at sequence subscript  i
//  and add it to  orf_list , if it's sufficiently long.   frame  is
//  the reading frame subscript of this orf.   prev_fwd_stop  indicates
//  the location of the previous forward stop codons.   first_fwd_start
//  has the locations of the first start codon for the current forward
//  reading frames.  Set  first_fwd_stop  to this position if there
//  is no prior stop in this reading frame.   first_base  is the position
//  of the first sequence base after an ignore region, or the start of
//  the sequence if no ignore regions have been encountered, which is
//  indicated by  hit_ignore.

{
     Orf_t  orf;
     int  gene_len, orf_len;
     
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
     
     if  (gene_len >= Min_Gene_Len || ((Allow_Indels || Allow_Subs) && orf_len >= Min_Indel_ORF_Len))
     {
	  orf . Set_Stop_Position (i - 1);
	  orf . Set_Frame (1 + (frame + 1) % 3);
	  orf . Set_Gene_Len (gene_len);
	  orf . Set_Orf_Len (orf_len);
	  orf_list . push_back (orf);
     }
     
     first_fwd_start [frame] = INT_MAX;
     prev_fwd_stop [frame] = i - 1;
     
     return;
}

void  Do_Rev_Stop_Codon
    (int i, int frame, int prev_rev_stop [3], int last_rev_start [3],
     bool hit_ignore, vector <Orf_t> & orf_list)
{
     Orf_t  orf;
     int gene_len, orf_len;
     int orf_stop = 0; // to avoid warning

     if  (prev_rev_stop [frame] == 0)
	  Handle_First_Reverse_Stop (i - 1, last_rev_start [frame],
				     gene_len, orf_stop, hit_ignore);
     else
     {
	  orf_stop = prev_rev_stop [frame];
	  gene_len = last_rev_start [frame] - orf_stop;
     }

     orf_len = i - orf_stop - 4;
		    
     if  (gene_len >= Min_Gene_Len || ((Allow_Indels || Allow_Subs) && orf_len >= Min_Indel_ORF_Len))
     {
	  orf . Set_Stop_Position (orf_stop);
	  orf . Set_Frame (-1 - (frame + 1) % 3);
	  orf . Set_Gene_Len (gene_len);
	  orf . Set_Orf_Len (orf_len);
	  orf_list . push_back (orf);
     }
     last_rev_start [frame] = 0;
     prev_rev_stop [frame] = i - 1;

     return;
}


void  Echo_Specific_Settings
    (FILE * fp, int len)

//  Output values of variables an settings that depend on the
//  current input string, which has length  len .

  {
   fprintf (fp, "Sequence length = %d\n", len);

   return;
  }


int  Find_Uncovered_Position
    (vector <Event_Node_t *> ep)

//  Find a position in  ep  that is not covered by any potential
//  gene, if possible.  If the first gene does not overlap the
//  last gene, then return  0 .  Also return  0  if there is
//  no uncovered position.  The position is regarded as being
//  between bases, and positions are numbered from  0  to  Sequence_Len .

  {
   int  cover_ct, zero_pos;
   int  first_pos, last_pos;
   int  i, n;

   n = ep . size ();

   if  (n <= 1)
       return  0;

   // ep is already sorted ascending by position and the initial
   // event is first in it
   first_pos = ep [1] -> pos - 3;  // between position in front of codon
   last_pos = ep [n - 1] -> pos - Sequence_Len;
     // between position after codon normalized to wrapped front position

   if  (last_pos <= first_pos)
       return  0;  // no overlap between front and back

   cover_ct = 0;
   zero_pos = ep [n - 1] -> pos;
   for  (i = 1;  i < n;  i ++)
     switch  (ep [i] -> e_type)
       {
        case  FWD_START :
          if  (ep [i] -> is_first_start)
              {
               cover_ct ++;
               if  (cover_ct == 1 && 3 <= ep [i] -> pos - zero_pos)
                   {
                    assert (zero_pos >= 1);
                    return  zero_pos;
                   }
              }
          break;

        case  FWD_STOP :
          cover_ct --;
          if  (cover_ct == 0)
              zero_pos = ep [i] -> pos;
          break;

        case  REV_START :
          if  (ep [i] -> is_first_start)
              {
               cover_ct --;
               if  (cover_ct == 0)
                   zero_pos = ep [i] -> pos;
              }
          break;

        case  REV_STOP :
          cover_ct ++;
          if  (cover_ct == 1 && 3 <= ep [i] -> pos - zero_pos)
              {
               assert (zero_pos >= 1);
               return  zero_pos;
              }
          break;

        case  INITIAL :
        case  TERMINAL :
        default :
          sprintf (Clean_Exit_Msg_Line, "ERROR:  Unexpected event type = %s",
               Print_String (ep [i] -> e_type));
          SIMPLE_THROW (Clean_Exit_Msg_Line);
       }

   return  0;
  }



void  Find_Orfs
    (vector <Orf_t> & orf_list)

//  Put in  orf_list  all sufficiently long orfs in global
//  string  Sequence .

{
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
     int  ignore_ct = 0;  // to avoid warning
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
     int  frame;
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
	  // check if this position is the boundary of an ignore region
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
		    last_rev_start [frame] = i - 1;
	       
	       if  (codon . Must_Be (Fwd_Stop_Pattern, which))
		    Do_Fwd_Stop_Codon (i, frame, prev_fwd_stop, first_fwd_start,
				       first_fwd_stop, first_base, hit_ignore, orf_list);
	       
	       if  (codon . Must_Be (Rev_Stop_Pattern, which))
		    Do_Rev_Stop_Codon (i, frame, prev_rev_stop, last_rev_start,
				       hit_ignore, orf_list);
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
     else if  (Allow_Truncated_Orfs)
	  // Treat 3 bp past the end of the sequence as stop codons
	  for  ( ;  i < n + 3;  i ++)
	  {
	       if  (! ignoring)
		    Do_Fwd_Stop_Codon (i, frame, prev_fwd_stop, first_fwd_start,
				       first_fwd_stop, first_base, hit_ignore, orf_list);
	       
	       if  (frame == 2)
		    frame = 0;
	       else
		    frame ++;
	  }
     
     return;
}


void  Finish_Orfs
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
     int  frame, gene_len, orf_len, orf_stop;
     
     for  (frame = 0;  frame < 3;  frame ++)
     {
	  Handle_Last_Reverse_Stop (frame, prev_rev_stop, last_rev_start,
				    gene_len, orf_len, orf_stop, use_wraparound, last_position);
	       
	  if  (gene_len >= Min_Gene_Len || ((Allow_Indels || Allow_Subs) && orf_len >= Min_Indel_ORF_Len))
	  {
	       orf . Set_Stop_Position (orf_stop);
	       orf . Set_Frame (-1 - (frame + 1) % 3);
	       orf . Set_Gene_Len (gene_len);
	       orf . Set_Orf_Len (orf_len);
	       orf_list . push_back (orf);
	  }
     }
     
     return;
  }


int  Frame_To_Sub
    (int f)

//  Return the subscript equivalent of frame  f .

  {
   if  (f > 0)
       return  f - 1;
     else
       return  2 - f;
  }


void  Get_Ignore_Regions
    (void)

//  Read the list of regions from the file with name in global
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


vector<int> Get_Sequence_Lengths()

// Make a vector of all of the sequence lengths (in amino acids) 
// in the file given.

{
     vector<int> lengths;
     string seq, head;
     FILE  * sequence_fp;

     sequence_fp = File_Open (Sequence_File_Name, "r", __FILE__, __LINE__);
     while(Fasta_Read(sequence_fp, seq, head))
	  lengths.push_back(seq.size() / 3);
     fclose (sequence_fp);

     return lengths;
}


void  Handle_First_Forward_Stop
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
	  if  (Allow_Truncated_Orfs && gene_len < Min_Gene_Len)
	       gene_len = orf_len;
     }
     
     return;
  }



void  Handle_First_Reverse_Stop
    (int pos, int last_start, int & gene_len, int & orf_stop, bool hit_ignore)

//  Set  gene_len  to the length of the reverse-strand gene whose start
//  is at  last_start  (left base of start codon, start-at-1) and which
//  extends off the front of the sequence.  Set  orf_stop  to the first,
//  frame-correct position < 1 where the stop codon (left base) could be.
//  It doesn't matter if the 2nd or 3rd base of this stop codon placement
//  overlaps the beginning of the sequence.
//   pos  is the position (start-at-1 coords) of the right bounding stop
//  codon of this gene.   Set  gene_len  to zero and return, however,
//  if either  hit_ignore  is true or  Allow_Truncated_Orfs  is false.

{
     if  (hit_ignore || ! Allow_Truncated_Orfs)
     {
	  gene_len = 0;
	  return;
     }
     
     orf_stop = pos % 3;
     if  (orf_stop > 0)
	  orf_stop -= 3;
     gene_len = last_start - orf_stop;
     
     return;
  }



void  Handle_Last_Reverse_Stop
     (int fr, const int prev_rev_stop [3], const int last_rev_start [3],
      int & gene_len, int & orf_len, int & orf_stop, bool use_wraparound, int last_position)

//  Set  orf_len  and  gene_len  to the length of the last orf, and longest
//  gene in it, resp., in reverse reading frame  fr .
//   prev_rev_stop  has the last stop position in  Sequence  in each
//  reverse reading frame, and  last_rev_start  has the corresponding
//  last start locations.    use_wraparound  indicates whether the
//  orfs are allowed to wrap around the end of the (circular) genome.
//   last_position  is the highest-numbered sequence position available

{
     if  (prev_rev_stop [fr] == 0) {
	  if(fr == 0)
	       orf_stop = -1;
	  else if(fr == 1)
	       orf_stop = 0;
	  else
	       orf_stop = -2;
     }
     else
	  orf_stop = prev_rev_stop [fr];

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
	  orf_len = last_position - orf_stop - 2;
	  // round down to next multiple of 3
	  orf_len -= orf_len % 3;

	  if  (last_rev_start [fr] == 0)
	       gene_len = 0;
          else
	       gene_len = last_rev_start [fr] - orf_stop;
	  if  (Allow_Truncated_Orfs && gene_len < Min_Gene_Len)
	       gene_len = orf_len;
     }

     assert (orf_len % 3 == 0);
     assert (gene_len % 3 == 0);

     return;
}



void  Initialize_Terminal_Events
    (Event_Node_t & first_event, Event_Node_t & final_event,
     Event_Node_t * best_event [6], Event_Node_t * last_event [6])

//  Set up  first_event  and  final_event  and make all
//  entries in  best_event  and  last_event  point to
//   first_event .

  {
   int  i;

   first_event . e_type = INITIAL;
   first_event . pos = 0;
   first_event . score = 0.0;
   first_event . best_pred = NULL;
   first_event . frame_pred = NULL;

   for  (i = 0;  i < 6;  i ++)
     last_event [i] = best_event [i] = & first_event;

   final_event . e_type = TERMINAL;
   final_event . frame_pred = NULL;

   return;
  }


double  Olap_Score_Adjustment
    (int lo, int hi, int f1, int f2)

//  Return the larger of the frame  f1  and  frame  f2  scores
//  on the subsequence from  lo .. hi  of global  Sequence .
//   lo  and  hi  are inclusive, start at 1 coordinates.
//  Because wraparounds may have confused the frames, only the
//  sign of the frames is used.   f1  is assumed to be the
//  frame of the beginnning of the subsequence starting on
//  a codon boundary.   f2  is the corresponding frame at the
//  end of the sequence.

  {
   string  buff;
   double  s1, s2;
   int  len, fs;

   len = 1 + hi - lo;
   if  (len < 1)
       return  0.0;

   if  (lo < 1)
       lo += Sequence_Len;
   if  (lo > Sequence_Len)
       lo -= Sequence_Len;
   if  (hi < 1)
       hi += Sequence_Len;
   if  (hi > Sequence_Len)
       hi -= Sequence_Len;

   lo --;  // convert to subscript
   hi --;

   switch  (len % 3)
     {
      case  0 :
        fs = 1;
        break;
      case  1 :
        fs = 0;
        break;
      case  2 :
      default:
        fs = 2;
        break;
     }
   // fs is the frame subscript to use in scoring in the direction
   // that does not necessarily start on a codon boundary

   if  (f1 > 0)
       {
        Reverse_Transfer (buff, Sequence, hi, len);
        s1 = Gene_ICM . Score_String (buff . c_str (), len, fs)
               - Indep_Model . Score_String (buff . c_str (), len, fs);
       }
     else
       {
        Complement_Transfer (buff, Sequence, lo, len);
        s1 = Gene_ICM . Score_String (buff . c_str (), len, 1)
               - Indep_Model . Score_String (buff . c_str (), len, 1);
       }

   if  (f1 * f2 < 0)
       Reverse_Complement (buff);

   if  (f2 > 0)
       s2 = Gene_ICM . Score_String (buff . c_str (), len, 1)
              - Indep_Model . Score_String (buff . c_str (), len, 1);
     else
       s2 = Gene_ICM . Score_String (buff . c_str (), len, fs)
              - Indep_Model . Score_String (buff . c_str (), len, fs);

   return  Max (s1, s2);
  }


int  On_Seq_1
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


void Parse_Features(const char* feature_file)

// Parse distributions for features length, start codons,
// adjacent gene orientations, and adjacent gene distances.

{
     string line;
     vector<string> lv;
     ifstream ff(feature_file);

     if(!ff.good()) {
	  sprintf (Clean_Exit_Msg_Line, "ERROR:  Cannot open feature file %s\n",feature_file);
	  Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
     }

     float gene_count = 0;
     float nonorf_count = 0;
     vector<double> length_gene;
     vector<double> length_non;
     vector<float> start_gene;
     vector<float> start_non;
     vector<float> adjor_gene;
     vector<float> adjor_non;
     vector<float> adjdist_ff_gene;
     vector<float> adjdist_ff_non;
     vector<float> adjdist_fr_gene;
     vector<float> adjdist_fr_non;
     vector<float> adjdist_rf_gene;
     vector<float> adjdist_rf_non;

     while(getline(ff, line)) {
	  if(line.size() >= 4 && line.substr(0,4) == "DIST") {
	       lv = split(line);
	       if(lv.size() != 3)
		    Parse_Features_Err(line);
	       else {
		    string dist_type = upper(lv[1]);
		    string orf_type = upper(lv[2]);

		    if(dist_type == "START")
			 if(orf_type == "GENE")
			      Read_Start_Dist(ff, start_gene);
			 else if(orf_type == "NON")
			      Read_Start_Dist(ff, start_non);
			 else
			      Parse_Features_Err(line);

		    else if(dist_type == "LENGTH")
			 if(orf_type == "GENE")
			      gene_count = Read_Length_Dist(ff, length_gene);
			 else if(orf_type == "NON")
			      nonorf_count = Read_Length_Dist(ff, length_non);
			 else
			      Parse_Features_Err(line);		   

		    else if(dist_type == "ADJACENT_ORIENTATION")
			 if(orf_type == "GENE")
			      Read_Orient_Dist(ff, adjor_gene);
			 else if(orf_type == "NON")
			      Read_Orient_Dist(ff, adjor_non);
			 else
			      Parse_Features_Err(line);		   

		    else if(dist_type == "ADJACENT_DISTANCE_1_1")
			 if(orf_type == "GENE")
			      Read_Dist_Dist(ff, adjdist_ff_gene);
			 else if(orf_type == "NON")
			      Read_Dist_Dist(ff, adjdist_ff_non);
			 else
			      Parse_Features_Err(line);

		    else if(dist_type == "ADJACENT_DISTANCE_1_-1")
			 if(orf_type == "GENE")
			      Read_Dist_Dist(ff, adjdist_fr_gene);
			 else if(orf_type == "NON")
			      Read_Dist_Dist(ff, adjdist_fr_non);
			 else
			      Parse_Features_Err(line);

		    else if(dist_type == "ADJACENT_DISTANCE_-1_1")
			 if(orf_type == "GENE")
			      Read_Dist_Dist(ff, adjdist_rf_gene);
			 else if(orf_type == "NON")
			      Read_Dist_Dist(ff, adjdist_rf_non);
			 else
			      Parse_Features_Err(line);

		    else
			 Parse_Features_Err(line);
	       }
	  }
     }
     ff.close();

     if(gene_count > 0 && nonorf_count > 0) {
	  vector<int> seq_lengths = Get_Sequence_Lengths();

	  LogOdds_Prior = LogOdds_Fudge + log(gene_count/nonorf_count);
	  LogOdds_Length.Make_Log_Odds(length_gene, length_non, seq_lengths, (unsigned int)Min_Gene_Len);
	  User_Length = true;
     }

     if(!start_gene.empty()) {	  
	  LogOdds_Start.Make_Log_Odds(start_gene, start_non);
	  User_Start = true;
     }

     if(!adjor_gene.empty()) {
	  LogOdds_AdjOr.Make_Log_Odds(adjor_gene, adjor_non);

	  LogOdds_AdjDist.Set_MaxOverlap(Dist_Max_Overlap);
	  LogOdds_AdjDist.Make_Log_Odds_Fwd_Fwd(adjdist_ff_gene, adjdist_ff_non);
	  LogOdds_AdjDist.Make_Log_Odds_Fwd_Rev(adjdist_fr_gene, adjdist_fr_non);
	  LogOdds_AdjDist.Make_Log_Odds_Rev_Fwd(adjdist_rf_gene, adjdist_rf_non);

	  User_Adj = true;
     }

     if(Dave_Log) {
	  LogOdds_Length.Print(Output_Tag);
	  LogOdds_Start.Print(Output_Tag);
	  LogOdds_AdjOr.Print(Output_Tag);
	  LogOdds_AdjDist.Print(Output_Tag);
     }
}

void Parse_Features_Err(string line)
{
     sprintf (Clean_Exit_Msg_Line, "ERROR:  Feature distribution file must have 'DIST <FEATURE> <GENE|NON>\n%s\n",line.c_str());
     Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
}


int  Position_To_Frame
    (int p)

//  Return the reading frame corresponding to a codon beginning in
//  position  p .  Allow  p  to be negative.  For  p = ...,-2,-1,0,1,2,3,4,...
//  frames are, respectively,  ...,1,2,3,1,2,3,1,...

  {
   if  (p >= 0)
       return  1 + ((p + 2) % 3);
     else
       return  3 - ((-1 * p) % 3);
  }



void  Print_Comma_Separated_Strings
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



void  Print_Headings
    (FILE * fp)

//  Print column headings to  fp .

  {
   fputc ('\n', fp);

   fprintf (fp, "%4s %5s %17s %8s  %15s", "", "", "----- Start -----",
        "", "--- Length ----");
   if  (Use_Independent_Score)
       fprintf (fp, "  %s\n", "------------- Scores -------------");
     else
       fprintf (fp, "  %s\n", "----------- Scores ------------");
   fprintf (fp, "%4s %5s %8s %8s %8s  %7s %7s  %7s %5s %s",
        " ID ", "Frame", "of Orf", "of Gene", "Stop", "of Orf", "of Gene",
        "Raw", "InFrm", "F1 F2 F3 R1 R2 R3");
   if  (Use_Independent_Score)
       fprintf (fp, " NC");
   fprintf (fp, "\n");

   return;
  }


void Print_Start(Start_t & start, int stop_pos, Event_Node_t * ne)

// Print out all score components for a start codon.
//
// I screwed up the id numbers by assigning them only after hashing errors,
// but everything should still work fine.

{
     unsigned int error_i;
     int log_start, log_stop;
     int trunc_stop = 0;
     double log_length;

     if(ne->frame > 0) {
	  //log_start = ne->pos-3;
	  log_start = ne->pos-2;
	  if(start.truncated)
	       log_start -= 3;

	  log_stop = stop_pos+2;
	  if(stop_pos > Sequence_Len-2)
	       trunc_stop = 1;

	  log_length = LogOdds_Length.Score((unsigned int)((1+start.j)/3),
					    start.truncated, stop_pos > Sequence_Len-2,
					    Sequence_Len/3);
     } else {
	  log_start = ne->pos;
	  if(start.truncated)
	       log_start += 3;

	  //log_stop = stop_pos-1;
	  log_stop = stop_pos;
	  if(stop_pos < 1)
	       trunc_stop = 1;
	  
	  log_length = LogOdds_Length.Score((unsigned int)((1+start.j)/3),
					    start.truncated, stop_pos < 1,
					    Sequence_Len/3);
     }

     if(ne->pwm_score > -1000)
	  printf("%-10s %8d %8d %d %d %8.2f %8.2f %4d %8.2f %8.2f %8.2f",
		 Fasta_Header, log_start, log_stop, start.truncated, trunc_stop, start.score,
		 ne->pwm_score, ne->pwm_sep, LogOdds_Start.Score(start.which), log_length, ne->score);
     else
	  printf("%-10s %8d %8d %d %d %8.2f %8s %4s %8.2f %8.2f %8.2f",
		 Fasta_Header, log_start, log_stop, start.truncated, trunc_stop, start.score,
		 "-", "-", LogOdds_Start.Score(start.which), log_length, ne->score);

     if(!ne->errors.empty()) {
	  printf(" I:");
	  for(error_i = 0; error_i < ne->errors.size(); error_i++)
	       if(ne->errors[error_i].type == 0)
		    printf("%d,",ne->errors[error_i].pos);
	  printf(" D:");
	  for(error_i = 0; error_i < ne->errors.size(); error_i++)
	       if(ne->errors[error_i].type == 1)
		    printf("%d,",ne->errors[error_i].pos);
	  printf(" S:");
	  for(error_i = 0; error_i < ne->errors.size(); error_i++)
	       if(ne->errors[error_i].type == 2)
		    printf("%d,",ne->errors[error_i].pos);
     }
     printf("\n");
}


const char  * Print_String
    (Event_t e)

//  Return a printable equivalent for  e .

  {
   switch  (e)
     {
      case  INITIAL :
        return  "Initial";
      case  FWD_START :
        return  "F_Start";
      case  FWD_STOP :
        return  "F_Stop";
      case  REV_START :
        return  "R_Start";
      case  REV_STOP :
        return  "R_Stop";
      case  TERMINAL :
        return  "Terminal";
     }
   return  "None";
  }



void  Prob_To_Logs
    (vector <double> & v)

//  Convert the entries in  v  to their natural logarithms.
//  Add psuedo-count value for zero entries.  Normalize all
//  values in case the original values don't sum to 1.0

  {
   double  subtr;
   double  sum = 0.0, sum2 = 0.0;
   int  i, n;

   n = v . size ();
   for  (i = 0;  i < n;  i ++)
     {
      if  (v [i] < 0.0)
          {
           sprintf (Clean_Exit_Msg_Line, "ERROR:  Bad start codon probability %f\n",
                v [i]);
           Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
          }
      sum += v [i];
     }
   if  (sum == 0.0)
       {
        sprintf (Clean_Exit_Msg_Line, "ERROR:  Start codon probabilities all zero\n");
        Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
       }

   for  (i = 0;  i < n;  i ++)
     if  (v [i] == 0.0)
         {
          v [i] = sum * 1e-5;
          sum2 += v [i];
         }
   subtr = log (sum + sum2);

   for  (i = 0;  i < n;  i ++)
     v [i] = log (v [i]) - subtr;

   return;
  }



void  Process_Events
(void)

//  Find the best-scoring collection of genes represented by the
//  sequence of events in the global list of events pointed to by
//   Last_Event .

{
     vector <Event_Node_t *> ep;
     Event_Node_t  * p;
     int  i, n;

     // Make  ep  point to all the events
     // Also make the initial event's position smaller than the
     // position of any other event
     for  (i = 0;  i < 6;  i ++)
     {
	  int  min_pos = 0;

	  for  (p = Last_Event [i];  p != NULL && p -> e_type != INITIAL ;
                p = p -> frame_pred)
	  {
	       ep . push_back (p);
	       min_pos = Min (min_pos, p -> pos - 1);
	  }
	  if  (p == NULL)
          {
	       sprintf (Clean_Exit_Msg_Line, "ERROR:  Missing initial event\n");
	       Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
          }
	  p -> pos = Min (min_pos, p -> pos);
     }
     // Add a single copy of the initial event
     ep . push_back (p);

     n = int (ep . size ());

     // Sort all events into order by their  pos  field
     sort (ep . begin (), ep . end (), Event_Pos_Cmp);

     if  (Genome_Is_Circular)
     {
	  int  reference_pos;

	  reference_pos = Find_Uncovered_Position (ep);
	  if  (reference_pos > 0)
	       Shift_Events (ep, reference_pos);
     }

     // Scan  ep  and by dynamic programming find the best predecessor
     // event for each event.  Save the best event in each frame in
     // global  Best_Event [] .

     for  (i = 0;  i < n;  i ++)
     {
	  switch  (ep [i] -> e_type)
	  {
	  case  INITIAL :
	       Process_Initial_Event (ep [i]);
	       break;
	  case  FWD_START :
	  case  REV_STOP :
	       Process_Fwd_Start_Rev_Stop_Event (ep [i]);
	       break;
	  case  FWD_STOP :
	  case  REV_START :
	       Process_Fwd_Stop_Rev_Start_Event (ep [i]);
	       break;
	  default :
	       sprintf (Clean_Exit_Msg_Line, "ERROR:  Unexpected event type = %d\n",
			int (ep [i] -> e_type));
	       Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
	  }
     }

     return;
}


void  Process_Fwd_Start_Rev_Stop_Event
    (Event_Node_t * ep)

//  Process the forward-start-type or reverse-stop-type event pointed
//  to by  ep  by computing the best score that can be obtained by
//  combining it with prior events.

//  My adjacent gene code made this a lot more complicated.  I'm now
//  first checking all events that follow after the Best_Event (that
//  have score > 0), and then if the Best_Event is a reverse start,
//  I check for other preceding reverse starts.

  {
   int  i, f;
   Event_Node_t * rev_start_ptr, * max_event_ptr, * after_best_ptr;
   float this_score, max_score;
   int distance;

   f = Frame_To_Sub (ep -> frame);

   // Connect  ep  to the highest-scoring prior event and increment
   //  ep -> score  by that score

   // set initial max
   max_event_ptr = Best_Event[0];
   if(max_event_ptr->e_type == INITIAL) {
	max_score = max_event_ptr->score; // i.e. 0
   } else {
	distance = ep->pos - max_event_ptr->pos - 3;
	max_score = max_event_ptr->score +
	     LogOdds_AdjOr.Score(max_event_ptr->e_type, ep->e_type) +
	     LogOdds_AdjDist.Score(max_event_ptr->e_type, ep->e_type, distance);
   }

   for  (i = 0;  i < 6;  i ++) {
	// check all forward stop and reverse start events after best
	for(after_best_ptr = Last_Event[i]; after_best_ptr != Best_Event[i]; after_best_ptr = after_best_ptr->frame_pred) {
	     if((after_best_ptr->e_type == FWD_STOP || after_best_ptr->e_type == REV_START) &&
		after_best_ptr->score > 0) {
		  distance = ep->pos - after_best_ptr->pos - 3;
		  this_score = after_best_ptr->score +
		       LogOdds_AdjOr.Score(after_best_ptr->e_type, ep->e_type) +
		       LogOdds_AdjDist.Score(after_best_ptr->e_type, ep->e_type, distance);
		  
		  if(this_score > max_score) {
		       max_score = this_score;
		       max_event_ptr = after_best_ptr;
		  }
	     }
	}

	// if reverse, check multiple starts
	if(Best_Event[i]->e_type == REV_START) {
	     for(rev_start_ptr = Best_Event[i]; rev_start_ptr->e_type == REV_START; rev_start_ptr = rev_start_ptr->frame_pred) {
		  distance = ep->pos - rev_start_ptr->pos - 3;
		  this_score = rev_start_ptr->score +
		       LogOdds_AdjOr.Score(rev_start_ptr->e_type, ep->e_type) +
		       LogOdds_AdjDist.Score(rev_start_ptr->e_type, ep->e_type, distance);

		  if(this_score > max_score) {
		       max_score = this_score;
		       max_event_ptr = rev_start_ptr;
		  }
	     }

	// if forward, just check stop
	} else if(Best_Event[i]->e_type == FWD_STOP) {
	     distance = ep->pos - Best_Event[i]->pos - 3;
	     this_score = Best_Event[i]->score +
		  LogOdds_AdjOr.Score(Best_Event[i]->e_type, ep->e_type) +
		  LogOdds_AdjDist.Score(Best_Event[i]->e_type, ep->e_type, distance);

	     if(this_score > max_score) {
		  max_score = this_score;
		  max_event_ptr = Best_Event[i];
	     }

	// if initial, ignore adjacency
	} else {
	     assert(Best_Event[i]->e_type == INITIAL);

	     this_score = Best_Event[i]->score;
	     if(this_score > max_score) {
		  max_score = this_score;
		  max_event_ptr = Best_Event[i];
	     }
	}
   }

   ep -> best_pred = max_event_ptr;
   ep -> score += max_score;

   // Make  ep  the last in the chain of events in this reading frame
   ep -> frame_pred = Last_Event [f];
   Last_Event [f] = ep;

   return;
  }


void  Process_Initial_Event
    (Event_Node_t * ep)

//  Process the initial-type event pointed to by  ep  by adding
//  it to the global lists  Best_Event []  and  Last_Event [] .

  {
   int  i;

   for  (i = 0;  i < 6;  i ++)
     Best_Event [i] = Last_Event [i] = ep;

   ep -> pos = 0;
   ep -> score = 0.0;
   ep -> frame_pred = ep -> best_pred = NULL;

   return;
  }



void  Process_Fwd_Stop_Rev_Start_Event
(Event_Node_t * ep)

//  Process the forward-stop-type event pointed to by  ep  by 
//  connecting it to the best previous start codon in the same frame
//  with the same gene id. If that score is better than the best score
//  in the frame, then make  Best_Event  for the frame point to  ep .
//  Also check for allowed overlaps with prior forward starts or reverse stops.

//  Process the reverse-start-type event pointed to by  ep  by computing
//  the best score that can be obtained by combining it with
//  prior events.

{
     Event_Node_t  * p;
     Event_Node_t * best_p;
     double mx;
     int  i, f;
     float old_adj_score, new_adj_score;
     int distance;
     double score_needed;
     Event_Node_t  * q;
     double  adj, diff;
     int  lo;
     unsigned int error_i;
     bool overlap_error;
     const float adj_score_buf = 0.0; // good in theory, slightly bad in practice

     f = Frame_To_Sub (ep -> frame);

     if(ep->e_type == FWD_STOP)
     {
	  // Find the best start codon and make  ep  point back to it
	  mx = -DBL_MAX;
	  best_p = NULL;
	  for  (p = Last_Event [f];  p -> e_type != INITIAL;  p = p -> frame_pred)
	  {
	       //for  (p = Last_Event [f];  p -> e_type == FWD_START;  p = p -> frame_pred) {
	       //cout << ep->pos << " " << p->id << " " << p->pos << " " << p->score << endl;
	       if  (p->id == ep->id && p -> score > mx)
	       {
		    mx = p -> score;
		    best_p = p;
	       }
	  }
	  ep -> best_pred = best_p;
	  ep -> score = mx;
	  //ep->errors = best_p->errors;
     }
     else
     {
	  // Connect  ep  to its corresponding reverse-stop event and increment
	  //  ep -> score  by that score
	  for  (p = Last_Event [f];  p != NULL && (p -> e_type == REV_START || p->id != ep->id);
		p = p -> frame_pred)
	       ;
	  if  (p == NULL || p -> e_type != REV_STOP)
	  {
	       sprintf (Clean_Exit_Msg_Line,
			"ERROR:  No reverse stop for reverse start at pos = %d\n", ep -> pos);
	       Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
	  }
	  ep -> best_pred = p;
	  ep -> score += p -> score;
     }

     // Check any events that represent genes that may overlap this one
     // by less than the allowable overlap threshold and adjust their
     // score and make them point to  ep  if it gives a better score
     if  (Best_Event [f] -> score < ep -> score + adj_score_buf)
     {
	  Disqualify (p, 3 + ep -> pos - Max_Olap_Bases);

	  // only update best event if better w/o adjacency scores
	  if(Best_Event [f] -> score < ep -> score)
	       Best_Event [f] = ep;

	  for  (i = 0;  i < 6;  i ++)
          {
	       for  (p = Last_Event [i];
		     p != NULL && 3 + ep -> pos - p -> pos <= Max_Olap_Bases;
		     p = p -> frame_pred)
	       {
		    // only concerned with left-hand-side of genes
		    if  (! p -> disqualified && (p -> e_type == FWD_START || p -> e_type == REV_STOP))
		    {
			 // our score must be better than the previous connection
			 if  (p -> best_pred == NULL)
			      score_needed = 0.0;
			 else
			      score_needed = p -> best_pred -> score;

			 if  (score_needed < ep -> score + adj_score_buf)
			 {
			      if  (p -> e_type == FWD_START)
				   lo = p -> pos - 2;
			      else
				   lo = p -> pos + 1;

			      // check for error in overlapping region
			      overlap_error = false;
			      for(error_i = 0; error_i < ep->errors.size(); error_i++)
				   if(p->pos-2 <= ep->errors[error_i].pos)
					overlap_error = true;
			      for(error_i = 0; error_i < p->errors.size(); error_i++)
				   if(p->errors[error_i].pos <= ep->pos)
					overlap_error = true;

			      if(!overlap_error)
			      {
				   // INDEL BUG!
				   // Score adjustment in the overlapping region rescores the region
				   // in the frames of the two genes being considered. In the case
				   // of a gene with an indel, even if the indel is not in the overlapping
				   // region, we may get the frame wrong if it has shifted.

				   // adjust ICM score in overlapping region

				   /*
				   if(ep->e_type == FWD_STOP)
					adj = Max(0.0, Olap_Score_Adjustment(lo, ep->pos-3, p->frame, ep->frame));
				   else
					adj = Max(0.0, Olap_Score_Adjustment(lo, ep->pos, p->frame, ep->frame));
				   */
				   adj = 0.0;
				   diff = ep->score - p->best_pred->score - adj;

				   // adjust adjacent gene scores
				   if(p->best_pred == NULL || p->best_pred->e_type == INITIAL)
					old_adj_score = 0;
				   else {
					distance = p->pos - p->best_pred->pos - 3;
					old_adj_score = LogOdds_AdjOr.Score(p->best_pred->e_type, p->e_type) +
					     LogOdds_AdjDist.Score(p->best_pred->e_type, p->e_type, distance);
				   }
				   distance = p->pos - ep->pos - 3;
				   new_adj_score = LogOdds_AdjOr.Score(ep->e_type, p->e_type) +
					LogOdds_AdjDist.Score(ep->e_type, p->e_type, distance);
				   diff += new_adj_score - old_adj_score;

				   // cout << Fasta_Header << ": considering changing (" << p->id << "," << p->pos << ") to (" << ep->id << "," << ep->pos << ") diff=" << diff << endl;

				   if  (diff > 0)
				   {
					p -> score += diff;
					p -> best_pred = ep;
				   
					// change all events in frame to this stop codon?
					// I doubt this happens.
					for  (q = Last_Event [i];  q != p;  q = q -> frame_pred)
					     if  (q -> best_pred == p)
						  q -> score += diff;
				   }
			      }
			 }
		    }
	       }
	  }
	  Requalify (p, 3 + ep -> pos - Max_Olap_Bases);
     }

     // Make  ep  the last in the chain of events in this reading frame
     ep -> frame_pred = Last_Event [f];
     Last_Event [f] = ep;

     return;
}


void  PWM_Meta_Score_Fwd_Start
    (int pos, double & score, int & separation)

//  Find the highest scoring match for Meta_Ribosome_PWMs
//  against the sequence in a window of length Ribosome_Window_Size
//  in front of position  pos  (numbered starting at 1) in the
//  forward direction.  Set  score  to the highest score and
//  set  separation  to the number of positions between the best
//  match and  pos .

//  DRK: I could do a better job here by making the GC% 
//  global, and by updating it on the fly rather than
//  recomputing it for each window. Also, I threw away 
//  some circular genome code.

{
     double  sc;
     int  bottom, lo, sep;
     int  j, n;

     score = 0.0;
     separation = 0;
     if (Meta_Ribosome_PWMs.empty())
	  return;

     int save_pos = pos-1;
     if(Meta_PWM_Save[save_pos].second != 999) {
	  score = Meta_PWM_Save[save_pos].first;
	  separation = Meta_PWM_Save[save_pos].second;

     } else {
	  unsigned int pwm_i;
	  unsigned int pwm_num = Meta_Ribosome_PWMs.size();
	  vector<double> cond_p(pwm_num);
	  double gc_lp;

	  double gc_log = log (0.5 * Indep_GC_Frac);
	  double at_log = log (0.5 * (1.0 - Indep_GC_Frac));
	  double nt_lp[] = {at_log, gc_log, gc_log, at_log};     

	  n = Meta_Ribosome_PWMs[0].Width ();
	  bottom = pos - Ribosome_Window_Size - 1;

	  score = - DBL_MAX;
	  separation = sep = 0;
	  for  (lo = pos - n - 1;  0 <= lo && bottom <= lo;  lo --, sep ++) {
	       // initialize probabilities
	       for(pwm_i = 0; pwm_i < pwm_num; pwm_i++)
		    cond_p[pwm_i] = 1.0;
	       gc_lp = 0.0;

	       // compute probabilities
	       for(j = 0; j < n; j++) {
		    for(pwm_i = 0; pwm_i < pwm_num; pwm_i++) {
			 cond_p[pwm_i] *= Meta_Ribosome_PWMs[pwm_i].Column_Score(Sequence[lo+j],j);
		    }
		    gc_lp += nt_lp[Nucleotide_To_Subscript(Sequence[lo+j])];
	       }

	       // compute log of mixture probability
	       sc = 0.0;
	       for(pwm_i = 0; pwm_i < pwm_num; pwm_i++)
		    sc += cond_p[pwm_i];
	       sc = log(sc / (double)pwm_num) - gc_lp;

	       if  (sc > score) {
		    score = sc;
		    separation = sep;
	       }
	  }

	  // save
	  pair<double,int> this_pwm(score,separation);
	  Meta_PWM_Save[save_pos] = this_pwm;
     }

     return;
}


void  PWM_Meta_Score_Rev_Start
    (int pos, double & score, int & separation)

//  Find the highest scoring match for Meta_Ribosome_PWMs
//  against the sequence in a window of length Ribosome_Window_Size
//  following position  pos  (numbered starting at 1) on the
//  reverse-complement strand.  Set  score  to the highest score and
//  set  separation  to the number of positions between the best
//  match and  pos .

//  DRK: I could do a better job here by making the GC% 
//  global, and by updating it on the fly rather than
//  recomputing it for each window. Also, I threw away 
//  some circular genome code.

{
     double  sc;
     int  top, hi, sep;
     int  j, n;

     score = 0.0;
     separation = 0;
     if (Meta_Ribosome_PWMs.empty())
	  return;

     int rev_pos = Sequence_Len + pos - 1;
     if(Meta_PWM_Save[rev_pos].second != 999) {
	  score = Meta_PWM_Save[rev_pos].first;
	  separation = Meta_PWM_Save[rev_pos].second;

     } else {

	  unsigned int pwm_i;
	  unsigned int pwm_num = Meta_Ribosome_PWMs.size();
	  vector<double> cond_p(pwm_num);
	  double gc_lp;

	  double gc_log = log (0.5 * Indep_GC_Frac);
	  double at_log = log (0.5 * (1.0 - Indep_GC_Frac));
	  double nt_lp[] = {at_log, gc_log, gc_log, at_log};

	  n = Meta_Ribosome_PWMs[0] . Width ();
	  top = pos - 1 + Ribosome_Window_Size;

	  score = - DBL_MAX;
	  separation = sep = 0;
	  for  (hi = pos - 1 + n;  hi < Sequence_Len && hi <= top;  hi ++, sep ++) {
	       // initialize probabilities
	       for(pwm_i = 0; pwm_i < pwm_num; pwm_i++)
		    cond_p[pwm_i] = 1.0;
	       gc_lp = 0.0;

	       // compute probabilities
	       for(j = 0; j < n; j++) {
		    for(pwm_i = 0; pwm_i < pwm_num; pwm_i++) {
			 cond_p[pwm_i] *= Meta_Ribosome_PWMs[pwm_i].Column_Score(Complement(Sequence[hi-j]),j);
		    }
		    gc_lp += nt_lp[Nucleotide_To_Subscript(Complement(Sequence[hi-j]))];
	       }

	       // compute log of mixture probability
	       sc = 0.0;
	       for(pwm_i = 0; pwm_i < pwm_num; pwm_i++)
		    sc += cond_p[pwm_i];
	       sc = log(sc / (double)pwm_num) - gc_lp;

	       if  (sc > score) {
		    score = sc;
		    separation = sep;
	       }
	  }
	  
	  // save
	  pair<double,int> this_pwm(score,separation);
	  Meta_PWM_Save[rev_pos] = this_pwm;
     }

     return;
}


void  PWM_Score_Fwd_Start
    (int pos, const PWM_t & pwm, int window, double & score, int & separation)

//  Find the highest scoring match for  pwm
//  against the sequence in a window of length  window
//  in front of position  pos  (numbered starting at 1) in the
//  forward direction.  Set  score  to the highest score and
//  set  separation  to the number of positions between the best
//  match and  pos .

  {
   double  sc;
   int  bottom, lo, sep;
   int  j, n;

   score = 0.0;
   separation = 0;

   if  (pwm . Is_Empty ())
       return;

   n = pwm . Width ();
   bottom = pos - window - 1;

   score = - DBL_MAX;
   separation = sep = 0;
   for  (lo = pos - n - 1;  0 <= lo && bottom <= lo;  lo --, sep ++)
     {
      sc = 0.0;
      for  (j = 0;  j < n;  j ++)
        sc += pwm . Column_Score (Sequence [lo + j], j);
      if  (sc > score)
          {
           score = sc;
           separation = sep;
          }
     }

   // handle wraparound here
   if  (Genome_Is_Circular)
       for  ( ;  bottom <= lo;  lo --, sep ++)
         {
          sc = 0.0;
          for  (j = 0;  j < n;  j ++)
            {
             int  k;

             k = lo + j;
             if  (k < 0)
                 k += Sequence_Len;
             sc += pwm . Column_Score (Sequence [k], j);
            }
          if  (sc > score)
              {
               score = sc;
               separation = sep;
              }
         }

   return;
  }



void  PWM_Score_Rev_Start
    (int pos, const PWM_t & pwm, int window, double & score, int & separation)

//  Find the highest scoring match for  pwm
//  against the sequence in a window of length  window
//  following position  pos  (numbered starting at 1) on the
//  reverse-complement strand.  Set  score  to the highest score and
//  set  separation  to the number of positions between the best
//  match and  pos .

  {
   double  sc;
   int  top, hi, sep;
   int  j, n;

   if  (pwm . Is_Empty ())
       {
        score = 0.0;
        separation = 0;
        return;
       }

   n = pwm . Width ();
   top = pos - 1 + window;

   score = - DBL_MAX;
   separation = sep = 0;
   for  (hi = pos - 1 + n;  hi < Sequence_Len && hi <= top;  hi ++, sep ++)
     {
      sc = 0.0;
      for  (j = 0;  j < n;  j ++)
        sc += pwm . Column_Score (Complement (Sequence [hi - j]), j);
      if  (sc > score)
          {
           score = sc;
           separation = sep;
          }
     }

   // handle wraparound here
   if(Genome_Is_Circular) {
	for  ( ;  hi <= top;  hi ++, sep ++)
	{
	     sc = 0.0;
	     for  (j = 0;  j < n;  j ++)
	     {
		  int  k;
		  
		  k = hi - j;
		  if  (Sequence_Len <= k)
		       k -= Sequence_Len;
		  sc += pwm . Column_Score (Complement (Sequence [k]), j);
	     }
	     if  (sc > score)
	     {
		  score = sc;
		  separation = sep;
	     }
	}
   }

   return;
  }


void Read_Dist_Dist(ifstream & ff, vector<float> & dist_dist)

// Read distribution from file, but only up to max_length.
// Add pseudocounts.
// Regress.
// Normalize.
// Also, set Dist_Max_Overlap by checking the file.  Compare it to
// Max_Olap_Bases and previous settings.  (-1 means not set yet).

{
     const unsigned int max_dist = 1000;
     const float pseudocount = 0.25;

     string line;
     vector<string> lv;
     int dist;
     float count;
     unsigned int l;

     dist_dist.clear();

     // handle first line separately to parse max overlap
     getline(ff, line);
     lv = split(line);
     dist = strtol(lv[0].c_str(), NULL, 10);
     if(Dist_Max_Overlap == -1) {
	  Dist_Max_Overlap = -1*dist;
	  if(Dist_Max_Overlap != Max_Olap_Bases) {
	       sprintf (Clean_Exit_Msg_Line, "ERROR:  Feature file max overlap %d does not match specified max overlap %d\n",Dist_Max_Overlap,Max_Olap_Bases);
	       Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
	  }
     } else {
	  if(Dist_Max_Overlap != -1*dist) {
	       sprintf (Clean_Exit_Msg_Line, "ERROR:  Max overlap in feature file differs by distribution\n");
	       Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
	  }
     }

     // add first
     count = strtod(lv[1].c_str(), NULL);
     dist_dist.push_back(count);

     // add rest
     while(getline(ff, line))
     {	  
	  lv = split(line);
	  if(lv.size() == 2) {
	       count = strtod(lv[1].c_str(), NULL);
	       // add next
	       dist_dist.push_back(count);
	  } else
	       break;
     }

     dist_dist.resize(Dist_Max_Overlap+max_dist);

     // add pseudocounts
     for(l = 0; l < dist_dist.size(); l++)
	  dist_dist[l] += pseudocount;

     // regress
     AdjDist_Smooth(dist_dist);

     // normalize
     float adjdist_sum = 0;
     for(l = 0; l < dist_dist.size(); l++)
	  adjdist_sum += dist_dist[l];
     for(l = 0; l < dist_dist.size(); l++)
	  dist_dist[l] /= adjdist_sum;
}


float Read_Length_Dist(ifstream & ff, vector<double> & length_dist)

// Read distribution from file, but only up to max_length.
// Add pseudocounts.
// Regress.
// Normalize.
// Return gene/orf count used to set prior probability.

// lenth_dist is filled with log probabilities, not probabilities!

{
     const unsigned int max_length = 2000;
     // nonparametric parameters
     //const double pseudocount = .01;
     const float sigma = 20;
     // parametric parameters
     double par_cumprob = 0.25;
      
     string line;
     vector<string> lv;
     unsigned int len, count, l;
     unsigned int min_aa_len = (unsigned int)ceil((float)Min_Gene_Len/3.0);
     vector<double> nonpar_length_dist(max_length);

     // read length histogram
     while(getline(ff, line))
     {	  
	  lv = split(line);
	  if(lv.size() == 2) {
	       len = strtol(lv[0].c_str(), NULL, 10);
	       count = strtol(lv[1].c_str(), NULL, 10);

	       // assuming they are ordered
	       if(len+1 > max_length)
		    nonpar_length_dist.resize(len+1);
	       nonpar_length_dist[len] = count;
	  } else
	       break;
     }

     // get count
     float total_count = 0;
     for(l = min_aa_len; l < max_length; l++)
	  total_count += nonpar_length_dist[l];

     // compute ML gamma parameters
     double k, theta;
     gamma_ml(k, theta, nonpar_length_dist);

     // compute gamma log probabilities
     vector<double> par_length_dist(max_length);
     double denom = k*log(theta) + gamma(k); // log(exp(gamma(k)));
     for(l = 0; l < max_length; l++)
	  par_length_dist[l] = (k-1)*log((double)l) - (double)l/theta - denom;

     // normalize using min length (distribution is chopped at tail but prob ok)
     log_normalize(par_length_dist, min_aa_len);

	  
     // no pseudocounts!

     // nonparametric regress
     kernel_smooth(nonpar_length_dist, sigma, max_length);
     
     // normalize using min length
     normalize(nonpar_length_dist, min_aa_len);

     // trim back to max_length, and convert to log
     nonpar_length_dist.resize(max_length);
     for(l = min_aa_len; l < max_length; l++)
	  nonpar_length_dist[l] = log(nonpar_length_dist[l]);

     // blend distributions
     length_dist.clear();
     length_dist.resize(max_length);
     Blend_Length(length_dist, par_length_dist, nonpar_length_dist, par_cumprob);

     if(Length_Log) {
	  ofstream len_out("nonpar.txt");
	  for(l = 0; l < max_length; l++)
	       len_out << nonpar_length_dist[l] << endl;
	  len_out.close();
	  len_out.open("par.txt");
	  for(l = 0; l < max_length; l++)
	       len_out << par_length_dist[l] << endl;
	  len_out.close();
	  len_out.open("blend.txt");
	  for(l = 0; l < max_length; l++)
	       len_out << length_dist[l] << endl;
	  len_out.close();
     }

     return total_count;
}


void Read_Orient_Dist(ifstream & ff, vector<float> & orients)

// Read distribution from file
// Add pseudocounts.
// Normalize.

{
     const float pseudocount = 1.0;
     string line;
     vector<string> lv;
     vector<string> ors;
     unsigned int i;

     orients.resize(4);
     while(getline(ff, line)) {
	  lv = split(line);
	  if(lv.size() == 2) {
	       ors = split(lv[0],',');
	       if(ors[0] == "1") {
		    if(ors[1] == "1")
			 orients[0] = strtol(lv[1].c_str(), NULL, 10);
		    else if(ors[1] == "-1")
			 orients[1] = strtol(lv[1].c_str(), NULL, 10);
		    else {
			 sprintf (Clean_Exit_Msg_Line,
				  "ERROR:  Please use 1 and -1 for adjacent gene orientations\n");
			 Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
		    }
	       } else if(ors[0] == "-1") {
		    if(ors[1] == "1")
			 orients[2] = strtol(lv[1].c_str(), NULL, 10);
		    else if(ors[1] == "-1")
			 orients[3] = strtol(lv[1].c_str(), NULL, 10);
		    else {
			 sprintf (Clean_Exit_Msg_Line,
				  "ERROR:  Please use 1 and -1 for adjacent gene orientations\n");
			 Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
		    }
	       } else {
		    sprintf (Clean_Exit_Msg_Line,
			     "ERROR:  Please use 1 and -1 for adjacent gene orientations\n");
		    Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
	       }

	  } else
	       break;
     }

     // add pseudocounts
     for(i = 0; i < orients.size(); i++)
	  orients[i] += pseudocount;

     // normalize
     float orient_sum = 0;
     for(i = 0; i < orients.size(); i++)
	  orient_sum += orients[i];
     for(i = 0; i < orients.size(); i++)
	  orients[i] /= orient_sum;	  
}


void Read_Start_Dist(ifstream & ff, vector<float> & start_dist)
{
     string line, start_codon;
     vector<string> lv;
     int start_code = 0;
     unsigned int s;

     start_dist.clear();
     start_dist.resize(3);

     while(getline(ff, line)) {
	  lv = split(line);
	  if(lv.size() == 2) {
	       start_codon = upper(lv[0]);
	       if(start_codon == "ATG")
		    start_code = 0;
	       else if(start_codon == "GTG")
		    start_code = 1;
	       else if(start_codon == "TTG")
		    start_code = 2;
	       else {
		    sprintf (Clean_Exit_Msg_Line, "ERROR:  Cannot recognize start codons aside from ATG, GTG, and TTG\n");
		    Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
	       }
	       start_dist[start_code] = strtol(lv[1].c_str(), NULL, 10);		    
	  } else
	       break;
     }

     // add pseudocounts
     for(s = 0; s < start_dist.size(); s++)
	  start_dist[s] += 1.0;

     // normalize
     float start_sum = 0;
     for(s = 0; s < start_dist.size(); s++)
	  start_sum += start_dist[s];
     for(s = 0; s < start_dist.size(); s++)
	  start_dist[s] /= start_sum;
}


void  Requalify
    (Event_Node_t * p, int cutoff)

//  Set the  disqualified  bit false for nodes reachable from
//   p  by  best_pred  pointers that have  pos  values at least
//  as great as  cutoff .

  {
   Event_Node_t  * q;

   if  (p == NULL)
       return;

   for  (q = p -> best_pred;  q != NULL && cutoff <= q -> pos;  q = q -> best_pred)
     q -> disqualified = false;

   return;
  }



void  Reverse_Complement_Transfer
    (string & buff, const string & s, int lo, int hi)

//  Copy to string  buff  the reverse complement of the substring
//  of  s  between positions  lo  and  hi  (which are
//  space-based coordinates).

  {
   int  i, j;

   assert (hi <= int (s . length ()));

   buff . resize (hi - lo);
   for  (j = 0, i = hi - 1;  i >= lo;  j ++, i --)
     buff [j] = Complement (s [i]);

   return;
  }



void  Reverse_Transfer
    (string & buff, const string & s, int start, int len)

//  Copy to string  buff  the substring of  s  starting at subscript
//   start  and going to the left for a length of  len .  Wraparound
//  end of  s  if necessary.  Do *NOT* reverse-complement.

{
     int  j, n, f;
     
     n = s . length ();
     assert (start < n);
     assert (0 <= len);
     
     f = 2;
     buff . resize (len);
     for  (j = 0;  j < len;  j ++, start --)
     {
	  //if(Dave_Log)
	  //     cout << "Scoring " << (start+1) << " " << s[start] << " " << (f+1) << "\n";
	  f = (f-1+3)%3;

	  buff [j] = s [start];
	  if  (start <= 0)
	       start += n;	  
     }
     
     return;
  }


void  Set_Final_Event
    (Event_Node_t & fe, Event_Node_t * best_event [6],
     int seq_len)

//  Set final event  fe , representing the end of genome,
//  and make it point back to the best event in  best_event .
//   seq_len  is the length of the entire genome sequence.

  {
   int  i;

   fe . pos = seq_len;
   fe . score = best_event [0] -> score;
   fe . best_pred = best_event [0];

   for  (i = 1;  i < 6;  i ++)
     {
      if  (best_event [i] -> score >= fe . score)
          {
           fe . score = best_event [i] -> score;
           fe . best_pred = best_event [i];
          }
     }

   return;
  }


void  Set_GC_Fraction ()
//    (double & gc, const vector <string> & s)

//  Set  gc  to the fraction of letters in all strings in  s  that are
//  'g' or 'c'.

{
     unsigned int  j, m, ct = 0, total = 0;
     FILE  * sequence_fp;
     string seq, header;
     
     sequence_fp = File_Open (Sequence_File_Name, "r", __FILE__, __LINE__);

     while(Fasta_Read(sequence_fp, seq, header))
     {
	  m = seq.length ();
	  total += m;
	  for  (j = 0;  j < m;  j ++)
	  {
	       seq[j] = Filter(tolower(seq[j]));
	       if  (seq[j] == 'g' || seq[j] == 'c')
		    ct ++;
	  }
     }

     fclose (sequence_fp);

     Indep_GC_Frac = double (ct) / total;
     GC_Frac_Set = true;

     return;
  }

void  Set_Ignore_Score_Len
    (void)

//  Set global  Ignore_Score_Len  to the length of the longest orf
//  that would be expected to occur once at random in a million bases.
//  Assume an over-simplified model with independent stop codons.
     
{
     
//   if  (Ignore_Score_Len == INT_MAX)
//       {
     double  poisson_lambda = 0.0;
     int  i, n;
     
     n = Stop_Codon . size ();
     for  (i = 0;  i < n;  i ++)
     {
	  double  x = 1.0;
	  int  j;
	  
	  for  (j = 0;  j < 3;  j ++)
	       if  (Stop_Codon [i] [j] == 'c' || Stop_Codon [i] [j] == 'g')
		    x *= Indep_GC_Frac / 2.0;
	       else
		    x *= (1.0 - Indep_GC_Frac) / 2.0;
	  
	  poisson_lambda += x;
     }
     
     assert (poisson_lambda != 0.0);
     Ignore_Score_Len
	  = (long int) floor (3.0 * log (2.0 * 1000000 * poisson_lambda)
			      / poisson_lambda);
//  }
     
     return;
}


void  Set_Start_And_Stop_Codons
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
        if  (Start_Prob . size () == 0)
            for  (i = 0;  i < n;  i ++)
              Start_Prob . push_back (DEFAULT_START_PROB [i]);
        else if  (Start_Codon . size () != Start_Prob . size ())
            {
             sprintf (Clean_Exit_Msg_Line,
                  "ERROR:  Different number of start codons & probs (%d & %d, resp.)\n",
                  int (Start_Codon . size ()), int (Start_Prob . size ()));
             Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
            }
       }
   else if  (Start_Prob . size () == 0)
       { // assign equal likelihood
        n = Start_Codon . size ();
        for  (i = 0;  i < n;  i ++)
          Start_Prob . push_back (1.0 / n);
       }
   else if  (Start_Codon . size () != Start_Prob . size ())
       {
        sprintf (Clean_Exit_Msg_Line,
             "ERROR:  Different number of start codons & probs (%d & %d, resp.)\n",
             int (Start_Codon . size ()), int (Start_Prob . size ()));
        Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
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

   n = Num_Start_Codons = Start_Codon . size ();
   for  (i = 0;  i < n;  i ++)
     {
      codon . Set_From (Start_Codon [i]);
      Fwd_Start_Pattern . push_back (codon);
      codon . Reverse_Complement ();
      Rev_Start_Pattern . push_back (codon);
     }

   n = Num_Stop_Codons = Stop_Codon . size ();
   for  (i = 0;  i < n;  i ++)
     {
      codon . Set_From (Stop_Codon [i]);
      Fwd_Stop_Pattern . push_back (codon);
      codon . Reverse_Complement ();
      Rev_Stop_Pattern . push_back (codon);
     }

   return;
  }



void  Shift_Events
    (vector <Event_Node_t *> & ep, int reference_pos)

//  Change the position of all events in  ep  that are before
//   reference_pos  by adding global  Sequence_Len  to them
//  and then sort according to the new positions

  {
   Event_Node_t  * frame_last [6];
   int  f, i, n, q;

   n = ep . size ();
   if  (n <= 1)
       return;

   for  (f = 0;  f < 6;  f ++)
     frame_last [f] = Last_Event [f];

   // Find the lowest-position event in each frame after  reference_pos
   // ep [0] is the initial-state event
   for  (q = n - 1;  q > 0 && reference_pos < ep [q] -> pos;  q --)
     {
      f = Frame_To_Sub (ep [q] -> frame);
      frame_last [f] = ep [q];
     }

   // Break the chain of events in each frame to skip over events
   // before  reference_pos
   for  (f = 0;  f < 6;  f ++)
     if  (reference_pos < frame_last [f] -> pos)
         frame_last [f] -> frame_pred = ep [0];
       else
         Last_Event [f] = ep [0];

   // Add the events before  reference_pos  onto the back of the
   // frame chains after incrementing positions.
   for  (i = 1;  i <= q;  i ++)
     {
      ep [i] -> pos += Sequence_Len;
      ep [i] -> Set_Frame_From_Pos ();
      f = Frame_To_Sub (ep [i] -> frame);
      ep [i] -> frame_pred = Last_Event [f];
      Last_Event [f] = ep [i];
     }

   // Sort all events into order by their  pos  field
   sort (ep . begin (), ep . end (), Event_Pos_Cmp);

   return;
  }



void  Show_Events
    (FILE * fp)

//  Display to  fp  the contents of the global lists of events
//  pointed to by  Last_Event .

  {
   vector <Event_Node_t *> ep;
   Event_Node_t  * p;
   int  i, n;

   for  (i = 0;  i < 6;  i ++)
     for  (p = Last_Event [i];  p != NULL;  p = p -> frame_pred)
       ep . push_back (p);

   n = int (ep . size ());

   // Sort all events into order by their  pos  field
   sort (ep . begin (), ep . end (), Event_Pos_Cmp);

   fprintf (fp, "\n%8s  %-8s  %2s  %10s\n", "Position", "Type", "Fr", "Score");
   for  (i = 0;  i < n;  i ++)
     fprintf (fp, "%8d  %-8s  %+2d  %10.2f\n", ep [i] -> pos,
          Print_String (ep [i] -> e_type), ep [i] -> frame, ep [i] -> score);

   return;
  }


void  Wrap_Around_Back
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



void  Wrap_Through_Front
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



void  Event_Node_t :: Set_Frame_From_Pos
    (void)

// Set the  frame  field of this node to the frame corresponding
// to the value in the  pos  field  but retaining the sign of
// the  frame  field.

  {
   int  f;

   assert (pos > 2);

   f = 1 + (pos % 3);
   if  (frame > 0)
       frame = f;
     else
       frame = -1 * f;

   return;
  }



