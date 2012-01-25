//  Programmers:   Arthur L. Delcher
//                 Doug Harmon
//
//  File:          icm.cc
//
//  Last Updated:  Mon Jun 12 15:34:00 EDT 2006
//
//  Purpose:  Routines for defining and manipulating the
//  Interpolated Context Model (ICM) used by Glimmer2
//
//
//  Copyright (c) 2006 University of Maryland Center for Bioinformatics
//  & Computational Biology


#include  "icm.hh"

using namespace std;

extern int  Verbose;



ICM_t :: ICM_t
    (int w, int d, int p)

//  Constructor for the ICM

  {
   int  i;

   model_len = w;
   model_depth = d;
   periodicity = p;
   num_nodes = (Int_Power (ALPHABET_SIZE, model_depth + 1) - 1) / (ALPHABET_SIZE - 1);
   score = (ICM_Score_Node_t * *)
                  Safe_calloc (periodicity, sizeof (ICM_Score_Node_t *),
                  __FILE__, __LINE__);
   for  (i = 0;  i < periodicity;  i ++)
     score [i] = (ICM_Score_Node_t *)
                   Safe_calloc (num_nodes, sizeof (ICM_Score_Node_t),
                   __FILE__, __LINE__);
   empty = true;
  }



ICM_t :: ~ ICM_t ()

//  Destroy this ICM

  {
   int  i;

   if  (score != NULL)
       {
        for  (i = 0;  i < periodicity;  i ++)
          free (score [i]);
        free (score);
       }
  }



void  ICM_t :: Build_Indep_WO_Stops
    (double gc_frac, const vector <const char *> & stop_codon)

//  Make this model represent generating codons with independent
//  nucleotides with GC-portion of  gc_frac  but without generating
//  any codons in  stop_codon .  The model is built in
//  reverse order of the strings in  stop_codon .

  {
   double  codon_prob [64], base_prob [4];
   double  sum;
   int  pattern [3];
   int  i, j, k, n;

   if  (model_len != 3 || model_depth != 2 || periodicity != 3
           || ALPHABET_SIZE != 4 || num_nodes != 21)
       {
        fprintf (stderr,
             "ERROR:  Incompatible ICM_Training_t for Build_Indep_WO_Stops\n");
        fprintf (stderr,
             "model_len = %d  model_depth = %d  periodicity = %d\n"
             "alphabet_size = %d  num_nodes = %d\n",
             model_len, model_depth, periodicity, ALPHABET_SIZE, num_nodes);
        fprintf (stderr,
             "Should be  %d ,  %d ,  %d ,  %d ,  %d  respectively\n",
             3, 2, 3, 4, 21);
        exit (EXIT_FAILURE);
       }

   // set base_prob to independent probability of a, c, g, t
   base_prob [1] = base_prob [2] = gc_frac / 2.0;        // c, g
   base_prob [0] = base_prob [3] = 0.5 - base_prob [1];  // a, t

   // set codon_prob to independent probabilities of all codons
   pattern [0] = pattern [1] = pattern [2] = 0;
   for  (i = 0;  i < 64;  i ++)
     {
      codon_prob [i] = base_prob [pattern [0]] * base_prob [pattern [1]]
                          * base_prob [pattern [2]];

      // increment pattern
      for  (j = 2;  j >= 0;  j --)
        {
         pattern [j] ++;
         if  (pattern [j] == 4)
             pattern [j] = 0;
           else
             break;
        }
     }

   // set  codon_prob  for all codons in  stop_codon  to near-zero value
   // Note:  Logically reverse  stop_codon  entries since scoring is done
   //  in the reverse direction (i.e., 3' to 5') of orfs
   n = stop_codon . size ();
   for  (i = 0;  i < n;  i ++)
     {
      j = Subscript (stop_codon [i] [0])
            + 4 * Subscript (stop_codon [i] [1])
            + 16 * Subscript (stop_codon [i] [2]);

      codon_prob [j] = 1e-20;
     }

   // normalize probability values
   sum = 0.0;
   for  (i = 0;  i < 64;  i ++)
     sum += codon_prob [i];
   for  (i = 0;  i < 64;  i ++)
     codon_prob [i] /= sum;

   // initialize  score  nodes
   for  (i = 0;  i < periodicity;  i ++)
     for  (j = 0;  j < num_nodes;  j ++)
       {
#if  STORE_MUT_INFO
        score [i] [j] . mut_info = 0.0;
#endif
        for  (k = 0;  k < 4;  k ++)
          score [i] [j] . prob [k] = 0.0;
       }

   // set root values of the ICM tree
   // these are the independent probabilities
   for  (i = 0;  i < periodicity;  i ++)
     {
      ICM_Score_Node_t  * p = score [i];
      int  d1 = Int_Power (4, (3 - i) % 3);

      if  (i == 1)
             // for frame i=1 this is independent of prior base which
             // is in the preceding codon
          p -> mut_info_pos = -1;
        else
          p -> mut_info_pos = 1;
      for  (j = 0;  j < 64;  j ++)
        p -> prob [(j / d1) % 4] += codon_prob [j];
     }

   // set level-1 values of the ICM tree
   for  (i = 0;  i < periodicity;  i ++)
     {
      ICM_Score_Node_t  * p = score [i] + 1;
      int  d1 = Int_Power (4, (3 - i) % 3);
      int  d2 = Int_Power (4, (4 - i) % 3);

      for  (j = 0;  j < 4;  j ++)
        if  (i == 2)
            p [j] . mut_info_pos = -1;
          else
            p [j] . mut_info_pos = 0;

      if  (i != 1)
          for  (j = 0;  j < 64;  j ++)
            p [(j / d2) % 4] . prob [(j / d1) % 4] += codon_prob [j];
     }

   // set level-2 values of the ICM tree
   // only need frame i=0 since other frames are stopped at
   // higher levels
   i = 0;
     {
      ICM_Score_Node_t  * p = score [i] + 5;
      int  d1 = Int_Power (4, (3 - i) % 3);
      int  d2 = Int_Power (4, (4 - i) % 3);
      int  d3 = Int_Power (4, (5 - i) % 3);

      for  (j = 0;  j < 16;  j ++)
        p [j] . mut_info_pos = -1;
      for  (j = 0;  j < 64;  j ++)
        {
         k = 4 * ((j / d2) % 4) + (j / d3) % 4;
         p [k] . prob [(j / d1) % 4] += codon_prob [j];
        }
     }

   // normalize and take logs of all prob values
   for  (i = 0;  i < periodicity;  i ++)
     for  (j = 0;  j < num_nodes;  j ++)
       {
        sum = 0.0;
        for  (k = 0;  k < 4;  k ++)
          sum += score [i] [j] . prob [k];
        for  (k = 0;  k < 4;  k ++)
          score [i] [j] . prob [k]
               = (sum == 0.0 ? 0.0 : log (score [i] [j] . prob [k] / sum));
       }
   
   empty = false;

   return;
  }



void  ICM_t :: Build_Reverse_Codon_WO_Stops
    (double codon_prob [64], const vector <const char *> & stop_codon)

//  Make this model represent generating independent codons
//  with proportions in  codon_prob  but without generating
//  any codons in  stop_codon .  The model is built in
//  reverse order of the strings in  stop_codon  and the entries
//  in  codon_prob  must be in alphabetical order by *REVERSE*
//  (but not reverse-complement) codon string, i.e., forward codons aaa, caa,
//  gaa, taa, ata, ....

  {
   double  sum;
   int  i, j, k, n;

   if  (model_len != 3 || model_depth != 2 || periodicity != 3
           || ALPHABET_SIZE != 4 || num_nodes != 21)
       {
        fprintf (stderr,
             "ERROR:  Incompatible ICM_Training_t for Build_Reverse_Codon_WO_Stops\n");
        fprintf (stderr,
             "model_len = %d  model_depth = %d  periodicity = %d\n"
             "alphabet_size = %d  num_nodes = %d\n",
             model_len, model_depth, periodicity, ALPHABET_SIZE, num_nodes);
        fprintf (stderr,
             "Should be  %d ,  %d ,  %d ,  %d ,  %d  respectively\n",
             3, 2, 3, 4, 21);
        exit (EXIT_FAILURE);
       }

   // set  codon_prob  for all codons in  stop_codon  to near-zero value
   // Note:  Logically reverse  stop_codon  entries since scoring is done
   //  in the reverse direction (i.e., 3' to 5') of orfs
   n = stop_codon . size ();
   for  (i = 0;  i < n;  i ++)
     {
      j = Subscript (stop_codon [i] [0])
            + 4 * Subscript (stop_codon [i] [1])
            + 16 * Subscript (stop_codon [i] [2]);

      codon_prob [j] = 1e-20;
     }

   // normalize probability values
   sum = 0.0;
   for  (i = 0;  i < 64;  i ++)
     sum += codon_prob [i];
   for  (i = 0;  i < 64;  i ++)
     codon_prob [i] /= sum;

   // initialize  score  nodes
   for  (i = 0;  i < periodicity;  i ++)
     for  (j = 0;  j < num_nodes;  j ++)
       {
#if  STORE_MUT_INFO
        score [i] [j] . mut_info = 0.0;
#endif
        for  (k = 0;  k < 4;  k ++)
          score [i] [j] . prob [k] = 0.0;
       }

   // set root values of the ICM tree
   // these are the independent probabilities
   for  (i = 0;  i < periodicity;  i ++)
     {
      ICM_Score_Node_t  * p = score [i];
      int  d1 = Int_Power (4, (3 - i) % 3);

      if  (i == 1)
             // for frame i=1 this is independent of prior base which
             // is in the preceding codon
          p -> mut_info_pos = -1;
        else
          p -> mut_info_pos = 1;
      for  (j = 0;  j < 64;  j ++)
        p -> prob [(j / d1) % 4] += codon_prob [j];
     }

   // set level-1 values of the ICM tree
   for  (i = 0;  i < periodicity;  i ++)
     {
      ICM_Score_Node_t  * p = score [i] + 1;
      int  d1 = Int_Power (4, (3 - i) % 3);
      int  d2 = Int_Power (4, (4 - i) % 3);

      for  (j = 0;  j < 4;  j ++)
        if  (i == 2)
            p [j] . mut_info_pos = -1;
          else
            p [j] . mut_info_pos = 0;

      if  (i != 1)
          for  (j = 0;  j < 64;  j ++)
            p [(j / d2) % 4] . prob [(j / d1) % 4] += codon_prob [j];
     }

   // set level-2 values of the ICM tree
   // only need frame i=0 since other frames are stopped at
   // higher levels
   i = 0;
     {
      ICM_Score_Node_t  * p = score [i] + 5;
      int  d1 = Int_Power (4, (3 - i) % 3);
      int  d2 = Int_Power (4, (4 - i) % 3);
      int  d3 = Int_Power (4, (5 - i) % 3);

      for  (j = 0;  j < 16;  j ++)
        p [j] . mut_info_pos = -1;
      for  (j = 0;  j < 64;  j ++)
        {
         k = 4 * ((j / d2) % 4) + (j / d3) % 4;
         p [k] . prob [(j / d1) % 4] += codon_prob [j];
        }
     }

   // normalize and take logs of all prob values
   for  (i = 0;  i < periodicity;  i ++)
     for  (j = 0;  j < num_nodes;  j ++)
       {
        sum = 0.0;
        for  (k = 0;  k < 4;  k ++)
          sum += score [i] [j] . prob [k];
        for  (k = 0;  k < 4;  k ++)
          score [i] [j] . prob [k]
               = (sum == 0.0 ? 0.0 : log (score [i] [j] . prob [k] / sum));
       }
   
   empty = false;

   return;
  }



void  ICM_t :: Cumulative_Score
    (const string & s, vector <double> & score, int frame)  const

//  Set  score [i]  to be the score of substring  s [0 .. i]
//  for each position in  s .  Use this model with the first base
//  in frame  frame .

  {
   double  result, x;
   const char  * cstr = s . c_str ();
   int  start, stop;
   int  i, n;

   if  (periodicity == 1)
       frame = 0;
   assert (0 <= frame && frame < periodicity);

   n = s . length ();
   score . resize (n);

   result = 0.0;

   stop = Min (model_len - 1, n);
   for  (i = 0;  i < stop;  i ++)
     {
      x = Partial_Window_Prob (i, cstr, frame);
//**ALD
      //printf ("Cumulative_Score:  i = %2d  ch = %c  frame = %d  prob = %6.4f\n", i, s [i], frame, exp (x));
      result += x;
      if  (frame == periodicity - 1)
          frame = 0;
        else
          frame ++;
      score [i] = result;
     }

   for  (start = 0;  i < n;  start ++, i ++)
     {
      x = Full_Window_Prob (cstr + start, frame);
//**ALD
      //printf ("Cumulative_Score:  %-.*s  start = %3d  i = %3d  ch = %c  frame = %d  prob = %6.4f\n",
//	      model_len, cstr + start, start, i, s [i], frame, exp (x));
      result += x;
      if  (frame == periodicity - 1)
          frame = 0;
        else
          frame ++;
      score [i] = result;
     }

   return;
  }



void  ICM_t :: Cumulative_Score_String
    (char * string, int len, int frame, double * cum_score)

//  Set entries in  cum_score  to cumulative score up to each respective base
//  of string  string [0 .. (len - 1)] .
//  Use this model with the first base in frame  frame .
//  Array  cum_score  is assumed to be large enough to hold the results.

  {
   double  result, x;
   int  start, stop = model_len - 1;
   int  i;
  
   if  (periodicity == 1)
       frame = 0;
   assert (0 <= frame && frame < periodicity);

   if  (Verbose > 0)
       printf ("Cumulative_Score_String  len = %d  frame = %d\n", len, frame);

   result = cum_score [0] = 0.0;

   for  (i = 0;  i < model_len - 1;  i ++)
     {
      x = Partial_Window_Prob (i, string, frame);
      if  (Verbose > 0)
          printf ("%7d: %8.3f\n", i, x);
      result += x;
      frame = (frame + 1) % periodicity;
      cum_score [i + 1] = result;
     }

   for  (start = 0;  stop < len;  start ++, stop ++)
     {
      x = Full_Window_Prob (string + start, frame);
      if  (Verbose > 0)
          printf ("%7d: %8.3f\n", start + model_len - 1, x);
      result += x;
      frame = (frame + 1) % periodicity;
      cum_score [stop + 1] = result;
     }

   return;
  }


void  ICM_t :: Display
    (FILE * fp)

//  Print an ASCII, human-readable version of this model to  fp
//  Assume  fp  is already open and that the caller will close it
//  if necessary.

  {
   int  i, period;

   fprintf (fp, "model_len = %d  periodicity = %d  depth = %d  num_nodes = %d\n",
            model_len, periodicity, model_depth, num_nodes);
            
   for  (period = 0;  period < periodicity;  period ++)
     {
      fprintf (fp, "period = %d\n", period);
      for  (i = 0;  i < num_nodes;  i ++)
        {
         int  j;

         fprintf (fp, "%3d:  %2d ", i, score [period] [i] . mut_info_pos);
         for  (j = 0;  j < ALPHABET_SIZE;  j ++)
           fprintf (fp, " %7.4f", exp (score [period] [i] . prob [j]));
         fputc ('\n', fp);
        }
     }
   return;
  }


void  ICM_t :: Frame_Score
    (const string & s, vector <double> & score, int frame)  const

//  Set  score [i]  to be the score of s [i] for each 
//  position in  s using only the frame given.

{
     const char  * cstr = s . c_str ();
     int  start, stop;
     int  i, n;

     assert (0 <= frame && frame < periodicity);

     n = s . length ();
     score . resize (n);

     stop = Min (model_len - 1, n);
     for  (i = 0;  i < stop;  i ++)
	  score [i] = Partial_Window_Prob (i, cstr, frame);

     for  (start = 0;  i < n;  start ++, i ++)
	  score [i] = Full_Window_Prob (cstr + start, frame);

     return;
}


void  ICM_t :: Full_Window_Distrib
    (char * string, int frame, float * dist)

//  Set  dist  to the probabilities of the possible characters
//  using  string  as the context and entries in  score [frame] .

  {
   int  num_node, i, pos, sub;

   num_node = 0;

   for  (i = 0;  i < model_depth;  i ++)
     {
      pos = score [frame] [num_node] . mut_info_pos;

      if  (pos == -1)
          break;

      if  (pos < -1)  // No information here or below in tree, go back up
                      // Shouldn't happen
          {
           num_node = PARENT (num_node);
           pos = score [frame] [num_node] . mut_info_pos;
           break;
          }

      sub = Subscript (string [pos]);

      num_node = (num_node * ALPHABET_SIZE) + sub + 1;
     }

   pos = score [frame] [num_node] . mut_info_pos;
   if  (pos < -1)
       {
        num_node = PARENT (num_node);
        pos = score [frame] [num_node] . mut_info_pos;
       }

   memcpy (dist, score [frame] [num_node] . prob, ALPHABET_SIZE * sizeof (float));

   return;
  }



double  ICM_t :: Full_Window_Prob
    (const char * string, int frame)  const

//  Return the log-probability of the last character in the first
//  model_len  bases of  string  conditioned on the preceding characters
//  using the entries in  score [frame] .

  {
   double  prob;
   int  num_node, i, pos, sub;

   num_node = 0;

   for  (i = 0;  i < model_depth;  i ++)
     {
      pos = score [frame] [num_node] . mut_info_pos;

      if  (pos == -1)
          break;

      if  (pos < -1)  // No information here or below in tree, go back up
                      // Shouldn't happen
          {
           num_node = PARENT (num_node);
           pos = score [frame] [num_node] . mut_info_pos;
           break;
          }

      sub = Subscript (string [pos]);

      num_node = (num_node * ALPHABET_SIZE) + sub + 1;
     }

   pos = score [frame] [num_node] . mut_info_pos;
   if  (pos < -1)
       {
        num_node = PARENT (num_node);
        pos = score [frame] [num_node] . mut_info_pos;
       }

   sub = Subscript (string [model_len - 1]);

   prob = (double) score [frame] [num_node] . prob [sub];

   if  (pos < -1)
       {
        fprintf (stderr, "WARNING:  prob = %.4f  pos = %d in  Full_Window_Prob\n",
                 prob, pos);
        fprintf (stderr, "num_node = %d\n",
                 num_node);
       }

   return  prob;
  }



void  ICM_t :: Input
    (FILE * fp)

//  Input the contents of this model from  fp , which has
//  already been opened.

  {
   char  line [ID_STRING_LEN];
   int  param [NUM_FIXED_LENGTH_PARAMS];
   int  node_id;
   int  prev_node;
   int  period;
   int  i;

   // free memory from previous version
   for  (i = 0;  i < periodicity;  i ++)
     free (score [i]);
   free (score);
   score = NULL;

   // skip the text header line
   if  (fread (line, sizeof (char), ID_STRING_LEN, fp) != unsigned (ID_STRING_LEN))
       {
        fprintf (stderr, "ERROR reading ICM header\n");
        exit (EXIT_FAILURE);
       };    

   if  (fread (param, sizeof (int), NUM_FIXED_LENGTH_PARAMS, fp) != NUM_FIXED_LENGTH_PARAMS)
       {
        fprintf (stderr, "ERROR reading parameters\n");
        exit (EXIT_FAILURE);
       }

   if  (ICM_VERSION_ID != param [0])
       {
        fprintf (stderr, "Bad ICM version = %d  should be %d\n",
                 param [0], ICM_VERSION_ID);
        exit (EXIT_FAILURE);
       }
   if  (ID_STRING_LEN != param [1])
       {
        fprintf (stderr, "Bad ID_STRING_LEN = %d  should be %d\n",
                 param [1], ID_STRING_LEN);
        exit (EXIT_FAILURE);
       }

   model_len = param [2];
   model_depth = param [3];
   periodicity = param [4];
   num_nodes = param [5];

   score = (ICM_Score_Node_t * *) Safe_malloc
             (periodicity * sizeof (ICM_Score_Node_t *), __FILE__, __LINE__);
   for  (i = 0;  i < periodicity;  i ++)
     score [i] = (ICM_Score_Node_t *) Safe_calloc
                   (num_nodes, sizeof (ICM_Score_Node_t), __FILE__, __LINE__);

   period = -1;
   prev_node = 0;
   while  (fread (& node_id, sizeof (int), 1, fp) != 0)
     {
      if  (node_id < 0)
          break;

      if  (node_id == 0)
          period ++;

      // read in the probabilities
      if  (fread (score [period] [node_id] . prob,
                  sizeof (float), ALPHABET_SIZE, fp) != unsigned (ALPHABET_SIZE))
          {
           fprintf (stderr, "ERROR reading icm node = %d  period = %d\n",
                    node_id, period); 
           exit (EXIT_FAILURE);
          }

      // read in the max mutual information position
      if  (fread (& (score [period] [node_id] . mut_info_pos), sizeof (short int), 1, fp)
             != 1)
          {
           fprintf (stderr, "ERROR reading mut_info_pos for node = %d  period = %d\n",
                    node_id, period);
           exit (EXIT_FAILURE);
          }

      // check for cut nodes
      if  (node_id != 0 && prev_node != node_id - 1)
          for  (i = prev_node + 1;  i < node_id;  i ++)
               score [period] [i] . mut_info_pos = -2;

      if  (node_id == 0 && period > 0)
          for  (i = prev_node + 1;  i < num_nodes;  i ++)
            score [period - 1] [i] . mut_info_pos = -2;

      prev_node = node_id;
     }

   if  (period != periodicity - 1)
       {
        fprintf (stderr, "ERROR:  Too few nodes for periodicity = %d\n",
                 periodicity);
        exit (EXIT_FAILURE);
       }

   // check for cut nodes in last period
   if  (prev_node != num_nodes - 1)
       for  (i = prev_node + 1;  i < num_nodes;  i ++)
            score [period] [i] . mut_info_pos = -2;

   empty = false;

   return;
  }


void  ICM_t :: Output
    (FILE * fp, bool binary_form)

//  Output the contents of this model to  fp .
//  If  binary_form  is true, then do it in binary; otherwise,
//  write an ascii text version.

  {
   int  end_marker = -1;
   int  i, frame;

   Write_Header (fp, binary_form);

   for  (frame = 0;  frame < periodicity;  frame ++)
     {
      Output_Node (fp, score [frame] + 0, 0, frame, binary_form);
      for  (i = 1;  i < num_nodes;  i ++)
        if  (score [frame] [i] . mut_info_pos >= -1)
            Output_Node (fp, score [frame] + i, i, frame, binary_form);
     }
   if  (binary_form)
       fwrite (& end_marker, sizeof (int), 1, fp);

   return;
  }



void  ICM_t :: Output_Node
    (FILE * fp, ICM_Score_Node_t * node, int id, int frame, bool binary_form)

//  Output  id  and then contents of  node  to  fp , in binary or ascii text
//  depending on whether  binary_form  is true or not, resp.
//  frame  is the frame within the periodic rotation of this node

  {
   if  (Verbose > 1)
       fprintf (stderr, "output node %d  frame %d\n", id, frame);

   if  (binary_form)
       {
        fwrite (& id, sizeof (int), 1, fp);
        fwrite (node -> prob, sizeof (float), ALPHABET_SIZE, fp);
        fwrite (& node -> mut_info_pos, sizeof (short int), 1, fp);
       }
     else
       {
        char  label [2 * 100];
          // allow extra positions to insert vertical separators
        int  i;

        assert (model_len <= 100);
        for  (i = 0;  i < model_len;  i ++)
          label [i] = '-';
        label [model_len - 1] = '?';
        label [model_len] = '\0';
        
        // put characters in label that represent the restrictions
        // on context positions for this node
        Set_Label_String (label, id, frame);

        if  (Verbose > 1)
            fprintf (stderr, "Label set to %s\n", label);

        fprintf (fp, "%6d  %s", id, label);
#if  STORE_MUT_INFO
        fprintf (fp, " %7.4f", node -> mut_info);
#endif
        for  (i = 0;  i < ALPHABET_SIZE;  i ++)
          fprintf (fp, " %6.3f", exp (node -> prob [i]));
        fputc ('\n', fp);
       }

   return;
  }



double  ICM_t :: Partial_Window_Prob
    (int predict_pos, const char * string, int frame)  const

//  Return the log-probability of character  string [predict_pos]  using
//  only the preceding  (predict_pos - 1)  characters  using the
//  entries in  score [frame] .

  {
   double  prob;
   int  num_node, i, start, pos, sub;

   start = predict_pos - (model_len - 1);
     // Negative position in  string  where this window would start
   num_node = 0;

   for  (i = 0;  i < model_depth;  i ++)
     {
      pos = start + score [frame] [num_node] . mut_info_pos;

      if  (pos < 0)
          break;

      sub = Subscript (string [pos]);

      num_node = (num_node * ALPHABET_SIZE) + sub + 1;
     }

   if  (score [frame] [num_node] . mut_info_pos == -2)
       num_node = PARENT (num_node);

   sub = Subscript (string [predict_pos]);

   prob = (double) score [frame] [num_node] . prob [sub];

   return  prob;
  }



void  ICM_t :: Read
    (char * path)

//  Read this model in from the file specified in  path

  {
   FILE  * fp;

   fp = File_Open (path, "r");     // Should be "rb"?

   Input (fp);

   fclose (fp);

   return;
  }


double  ICM_t :: Score_String
    (const char * string, int len, int frame)  const

//  Return the log-probability score of  string [0 .. (len - 1)]
//  in frame  frame  of this model.

  {
   double  result, x;
   int  start, stop = model_len - 1;
   int  i;
  
   if  (periodicity == 1)
       frame = 0;
   assert (0 <= frame && frame < periodicity);

   if  (Verbose > 0)
       printf ("Score_String  len = %d  frame = %d\n", len, frame);

   result = 0.0;

   for  (i = 0;  i < len && i < model_len - 1;  i ++)
     {
      x = Partial_Window_Prob (i, string, frame);
      if  (Verbose > 0)
          printf ("%7d: %8.3f\n", i, x);
      result += x;
      frame = (frame + 1) % periodicity;
     }

   for  (start = 0;  stop < len;  start ++, stop ++)
     {
      x = Full_Window_Prob (string + start, frame);
      if  (Verbose > 0)
          printf ("%7d: %8.3f\n", start + model_len - 1, x);
      result += x;
      frame = (frame + 1) % periodicity;
     }

   return  result;
  }



void  ICM_t :: Set_Label_String
    (char * label, int id, int frame)

//  Fill in characters of  label  with characters representing the
//  context matches for node  id  based on the ancestors of
//  this node in the tree in frame  frame .

  {
   int  last_separator, separator_ct;
   int  i, parent, mip;

   mip = score [frame] [id]. mut_info_pos;
   if  (mip >= 0)
       label [score [frame] [id]. mut_info_pos] = MAX_MI_CHAR;

   while  (id > 0)
     {
      parent = PARENT (id);
      label [score [frame] [parent] . mut_info_pos]
          = ALPHA_STRING [id - ALPHABET_SIZE * parent - 1];
      id = parent;
     }

   // add separators

   if  (periodicity == 1)
       last_separator = separator_ct = 0;
     else
       {
        if  (frame == 0)
            last_separator = model_len - periodicity;
          else
            last_separator = model_len - frame;
        if  (last_separator < 0)
            last_separator = 0;
        separator_ct = (last_separator + periodicity - 1) / periodicity;
       }

   for  (i = model_len;  i > 0;  i --)
     {
      label [i + separator_ct] = label [i];
      if  (i == last_separator)
          {
           separator_ct --;
           label [i + separator_ct] = SEPARATOR_CHAR;
           last_separator -= periodicity;
          }
     }

   return;
  }



void  ICM_t :: Write_Header
    (FILE * fp, bool binary_form)

//  Send to  fp  the parameter information for this model.
//  binary_form  determines whether the format is binary or
//  Ascii readable (for debugging purposes only)

  {

   if  (! binary_form)
       fprintf (fp, "ver = %.2f  len = %d  depth = %d"
                "  periodicity = %d  nodes = %d\n",
                ICM_VERSION_ID / 100.0, model_len, model_depth,
                periodicity, num_nodes);
     else
       {
        char  line [ID_STRING_LEN] = {'\0'};
        int  param [NUM_FIXED_LENGTH_PARAMS];

        sprintf (line, ">ver = %.2f  len = %d  depth = %d"
                 "  periodicity = %d  nodes = %d\n", 
                 ICM_VERSION_ID / 100.0, model_len, model_depth,
                 periodicity, num_nodes);
        assert (int (strlen (line)) < ID_STRING_LEN);
        fwrite (line, sizeof (char), ID_STRING_LEN, fp);

        param [0] = ICM_VERSION_ID;
        param [1] = ID_STRING_LEN;
        param [2] = model_len;
        param [3] = model_depth;
        param [4] = periodicity;
        param [5] = num_nodes;

        fwrite (param, sizeof (int), NUM_FIXED_LENGTH_PARAMS, fp);
       }

   return;
  }

void ICM_t :: Copy(ICM_t & icm) {
     empty = icm.empty;
     model_len = icm.model_len;
     model_depth = icm.model_depth;
     periodicity = icm.periodicity;
     num_nodes = icm.num_nodes;
     score = icm.score;
}


ICM_Training_t :: ICM_Training_t
    (int w, int d, int p)  :  ICM_t (w, d, p)

//  Constructor for the ICM

  {
   int  (* ptr) [ALPHA_SQUARED];
   int  i, j;

   train = (ICM_Training_Node_t * *)
                  Safe_calloc (periodicity, sizeof (ICM_Training_Node_t *),
                  __FILE__, __LINE__);

   if  (model_depth == 0)
       ptr = NULL;
     else
       ptr = count_memory = (int (*) [ALPHA_SQUARED])
               Safe_calloc (periodicity * num_nodes * (model_len - 1),
                            sizeof (int [ALPHA_SQUARED]), __FILE__, __LINE__);

   for  (i = 0;  i < periodicity;  i ++)
     {
      train [i] = (ICM_Training_Node_t *)
                    Safe_calloc (num_nodes, sizeof (ICM_Training_Node_t),
                    __FILE__, __LINE__);

      for  (j = 0;  j < num_nodes;  j ++)
        {
         train [i] [j] . count = ptr;
         ptr += model_len - 1;
        }
     }
  }



ICM_Training_t :: ~ ICM_Training_t ()

//  Destroy this ICM

  {
   int  i;

   free (count_memory);
   for  (i = 0;  i < periodicity;  i ++)
     free (train [i]);
   free (train);
  }



void  ICM_Training_t :: Complete_Tree
    (const vector <char *> & data)

//  Fill in the remaining nodes below the root in this model.
//  Starting with the root, restrict each character at the
//  mutual-information position to generate the child nodes.
//  For each child count how many times the sequence with the
//  designated characters occurs and from those counts
//  calculate the max mutual-information position for each child.
//  Keep doing this for every node until either  model_depth  is
//  reached or there is not enough data to go any deeper.
//  Use the IMM interpolation scheme to set probabilities
//  in the low-count case.

  {
   int  sub, max_pos, string_ct;
   double  best_info, next_info, used_info;
   int  sum;
   int  symbol;
     // subscript of character in alphabet
   int  first_node, last_node, nodes_on_level;
   int  frame, level;
   int  i, j, k;

   string_ct = int (data . size ());

   first_node = 1;
   nodes_on_level = ALPHABET_SIZE;

   for  (level = 1;  level <= model_depth;  level ++)
     {
      for  (i = 0;  i < string_ct;  i ++)
        Count_Char_Pairs_Restricted (data [i], level);

      last_node = first_node + nodes_on_level - 1;

      for  (frame = 0;  frame < periodicity;  frame ++)
        {
         symbol = 0;
         for  (sub = first_node;  sub <= last_node;
                 sub ++, symbol = (symbol + 1) % ALPHABET_SIZE)
           {
            int  final_char_ct [ALPHABET_SIZE] = {0};
              // number of occurrences of each symbol in the last position

            train [frame] [sub] . mut_info_seq = (short int) symbol;

            if  (score [frame] [PARENT (sub)] . mut_info_pos < 0)
                // Don't process this node; stopped at parent
                {
                 score [frame] [sub] . mut_info_pos = -2;
                 continue;
                }

            // sum over k of  count [i] [k]  is same for any i
            sum = 0;
            for  (i = k = 0;  i < ALPHABET_SIZE;  i ++)
              for  (j = 0;  j < ALPHABET_SIZE;  j ++)
                {
                 sum += train [frame] [sub] . count [0] [k];
                 final_char_ct [j] += train [frame] [sub] . count [0] [k];
                 k ++;
                }

            // find the position pair with the max mutual information
            max_pos = 0;
            best_info = Get_Mutual_Info (train [frame] [sub] . count [0],
                                         ALPHABET_SIZE, sum);
            used_info = best_info;

            for  (i = 1;  i < model_len - 1;  i ++)
              {
               next_info = Get_Mutual_Info (train [frame] [sub] . count [i],
                                            ALPHABET_SIZE, sum);
               if  (next_info >= best_info)
                   {
                    used_info = best_info = next_info;
                    max_pos = i;
                   }
               else if  (next_info >= (best_info / (1.0 + MUT_INFO_BIAS)))
                   {
                    // prefer positions to the right (i.e., closer to the
                    // predicted base) if mutual-information values are
                    // close enough
                    max_pos = i;
                    used_info = next_info;
                   }
              }

            if  (best_info <= MUT_INFO_EPSILON && sum < SAMPLE_SIZE_BOUND)
                // Not enough information gain; don't go down tree any further
                max_pos = -1;

            score [frame] [sub] . mut_info_pos = (short int) max_pos;
#if  STORE_MUT_INFO
            score [frame] [sub] . mut_info = float (used_info);
#endif

            if  (Verbose > 1)
                {
                 fprintf (stderr,
                      "frame = %d  node = %d  mut_info_pos = %d  mut_info = %.3f  cts: ",
                          frame, sub, max_pos, best_info);
                 for  (i = 0;  i < ALPHABET_SIZE;  i ++)
                   fprintf (stderr, " %4d", final_char_ct [i]);
                 fputc ('\n', stderr);
                }

            Interpolate_Probs (frame, sub, final_char_ct);

#if  0
// Should be in a separate method
            print_node (print_string, level, sub, mut_info[max_pos],frame, ending_sum);
#endif
           }
        }

      first_node = last_node + 1;
      nodes_on_level *= ALPHABET_SIZE;

      if  (Verbose > 0)
          fprintf (stderr, "Training done for level %d\n", level);
     }

   return;
  }



void  ICM_Training_t :: Count_Char_Pairs_Restricted
    (const char * string, int level)

//  For each complete window of length  model_len  in  string
//  determine the appropriate frame of the model and the node
//  at level  level  to which it should contribute counts.  Then
//  for each position  j = 0 .. (model_len - 2)  add 1 to
//  ct [j] [p]  where  p  is the index of the character
//  pair at positions  j  and  (model_len - 1)  in the window.

  {
   ICM_Training_Node_t  * node;
   int  start, stop, end, frame, last_char_sub;
   int  i, j;

   start = 0;
   end = int (strlen (string));
   frame = model_len % periodicity;

   for  (stop = model_len - 1;  stop < end;  start ++, stop ++)
     {
      node = Get_Training_Node (string + start, frame, level);
      if  (node != NULL)
          {
           last_char_sub = Subscript (string [stop]);
           for  (i = 0;  i < model_len - 1;  i ++)
             {
              j = ALPHABET_SIZE * Subscript (string [start + i])
                      + last_char_sub;
              node -> count [i] [j] ++;
             }
          }

      frame ++;
      if  (frame == periodicity)
          frame = 0;
     }

   return;
  }



ICM_Training_Node_t *  ICM_Training_t :: Get_Training_Node
    (const char * w, int frame, int level)

//  Find the node at level  level  in the  frame 'th segment
//  of the model that matches the string window  w .
//  Return a pointer to that node if it's valid; otherwise,
//  return  NULL;

  {
   int  i, j, sub;

   sub = 0;

   for  (i = 0;  i < level;  i ++)
     {
      j = score [frame] [sub] . mut_info_pos;
      if  (j < 0)
          return  NULL;

      sub = sub * ALPHABET_SIZE + Subscript (w [j]) + 1;
     }

   return  train [frame] + sub;
  }



void  ICM_Training_t :: Interpolate_Probs
    (int frame, int sub, int ct [])

//  Set the probabilities for node at subscript  sub  in the
//  frame 'th segment of the model using the frequencies in  ct
//  and interpolating with the probabilities in this node's
//  parent if the sum of the counts is sufficiently small

  {
   double  expected, chi2_stat, lambda, total_sum;
   int  parent;
   int  i;

   parent = PARENT (sub);

   total_sum = 0.0;
   for  (i = 0;  i < ALPHABET_SIZE;  i ++)
     total_sum += ct [i];

   // set probabilities directly including small bias from parent probabilities
   // to prevent zero probabilities
   for (i = 0;  i < ALPHABET_SIZE;  i ++)
     score [frame] [sub] . prob [i]
         = (ct [i] + PSEUDO_COUNT * score [frame] [parent] . prob [i])
               / (total_sum + PSEUDO_COUNT);

   // if there are enough samples those probabilities are OK and there
   // is no interpolation
   if  (total_sum >= SAMPLE_SIZE_BOUND)
       return;

   // calculate the chi-squared statistic
   chi2_stat = 0.0;
   for  (i = 0;  i < ALPHABET_SIZE;  i ++)
     {
      expected = total_sum * score [frame] [parent] . prob [i];
      if  (expected > 0.0)
          chi2_stat += pow (ct [i] - expected, 2.0) / expected;
     }
  
   // search for chi2_stat in table to get corresponding significance value
   for  (i = 0;  i < NUM_CHI2_ENTRIES && CHI2_VAL [i] < chi2_stat;  i ++)
     ;

   // determine interpolation coefficient lambda.  The assigned probs will be
   // lambda of this node's plus (1 - lambda) of the parent's
   if  (i == 0)
       lambda = 0.0;
   else if  (i == NUM_CHI2_ENTRIES)
       lambda = 1.0;
     else 
       lambda = CHI2_SIGNIFICANCE [i-1]
                  + ((chi2_stat - CHI2_VAL [i - 1])
                          / (CHI2_VAL [i] - CHI2_VAL [i - 1])) 
                      * (CHI2_SIGNIFICANCE [i] - CHI2_SIGNIFICANCE [i - 1]);

   // further weight lambda by the number of sample windows at this node
   lambda *= total_sum / SAMPLE_SIZE_BOUND;
   if  (lambda > 1.0)
       lambda = 1.0;

   // do the interpolation
   for  (i = 0;  i < ALPHABET_SIZE;  i ++)
     {
      score [frame] [sub] . prob [i] *= lambda;
      score [frame] [sub] . prob [i]
          += (1.0 - lambda) * score [frame] [parent] . prob [i];
     }

   return;
  }



void  ICM_Training_t :: Take_Logs
    (void)

//  Take natural logarithms of all probabilities in this model

  {
   int  i, j, frame;

   for  (frame = 0;  frame < periodicity;  frame ++)
     for  (i = 0;  i < num_nodes;  i ++)
       for  (j = 0;  j < ALPHABET_SIZE;  j ++)
         if  (score [frame] [i] . prob [j] > 0.0)
             score [frame] [i] . prob [j]
                 = float (log (score [frame] [i] . prob [j]));
           else
             score [frame] [i] . prob [j] = - FLT_MAX;

   return;
  }



void  ICM_Training_t :: Train_Model
    (const vector <char *> & data)

//  Calculate the probabilities for this model based on the
//  strings in  data .

  {
   int  frame, string_ct;
   
   string_ct = int (data . size ());

   for  (frame = 0;  frame < periodicity;  frame ++)
     {
      int  offset;
        // where first window should start in the string
      int  final_char_ct [ALPHABET_SIZE] = {0};
        // number of occurrences of each character as last in window
      double  best_info, next_info;
      int  max_pos, sum;
      int  i, j, k;

      offset = frame - (model_len % periodicity);
      if  (offset < 0)
          offset += periodicity;

      if  (model_depth == 0)
          {
           for  (i = 0;  i < string_ct;  i ++)
             Count_Single_Chars (final_char_ct, data [i] + offset,
                                 model_len, periodicity);
           sum = 0;
           for  (i = 0;  i < ALPHABET_SIZE;  i ++)
             sum+= final_char_ct [i];
           for  (i = 0;  i < ALPHABET_SIZE;  i ++)
             score [frame] [0] . prob [i]
               = (final_char_ct [i] + float (PSEUDO_COUNT / ALPHABET_SIZE))
                   / (sum + PSEUDO_COUNT);
           score [frame] [0] . mut_info_pos = -1;
          }
        else
          {
           for  (i = 0;  i < string_ct;  i ++)
             Count_Char_Pairs (train [frame] [0] . count,
                               data [i] + offset, model_len, periodicity);

           // sum over k of  count [i] [k]  is same for any i
           sum = 0;
           for  (i = k = 0;  i < ALPHABET_SIZE;  i ++)
             for  (j = 0;  j < ALPHABET_SIZE;  j ++)
               {
                sum += train [frame] [0] . count [0] [k];
                final_char_ct [j] += train [frame] [0] . count [0] [k];
                k ++;
               }
           for  (j = 0;  j < ALPHABET_SIZE;  j ++)
             score [frame] [0] . prob [j]
               = (final_char_ct [j] + float (PSEUDO_COUNT / ALPHABET_SIZE))
                   / float (sum + PSEUDO_COUNT);

           // find the position pair with the max mutual information
           max_pos = 0;
           best_info = Get_Mutual_Info (train [frame] [0] . count [0],
                                        ALPHABET_SIZE, sum);

           for  (i = 1;  i < model_len - 1;  i ++)
             {
              next_info = Get_Mutual_Info (train [frame] [0] . count [i],
                                           ALPHABET_SIZE, sum);
              if  (next_info >= best_info)
                  {
                   best_info = next_info;
                   max_pos = i;
                  }
              else if  (next_info >= (best_info / (1.0 + MUT_INFO_BIAS)))
                  max_pos = i;
                  // prefer positions to the right (i.e., closer to the
                  // predicted base) if mutual-information values are
                  // close enough
             }

           score [frame] [0] . mut_info_pos = (short int) max_pos;
#if  STORE_MUT_INFO
           score [frame] [0] . mut_info = float (best_info);
#endif
          }
      if  (Verbose > 1)
          {
           fprintf (stderr, "frame = %d  node = %d  mut_info_pos = %d  mut_info = %.3f  cts: ",
                    frame, 0, max_pos, best_info);
           for  (i = 0;  i < ALPHABET_SIZE;  i ++)
             fprintf (stderr, " %4d", final_char_ct [i]);
           fputc ('\n', stderr);
          }

#if  0
// Should be in a separate method
      print_node (print_string, 0, 0, mut_info [Node [frame] [0] . mut_info_pos],
                  frame, final_char_ct);
#endif
     }

   Complete_Tree (data);

   Take_Logs ();

   return;
  }



Fixed_Length_ICM_t :: Fixed_Length_ICM_t
    (int len, int sp, int * perm, ICM_Model_t mt)

//  Construct a  Fixed_Length_ICM_t  of length  len  using
//  perm  for the order of bases

  {
   length = len;
   if  (perm == NULL)
       permutation = NULL;
     else
       {
        int  i;

        permutation = new int [len];
        for  (i = 0;  i < len;  i ++)
          permutation [i] = perm [i];
       }
   special_position = sp;
   model_type = mt;
  }



Fixed_Length_ICM_t :: ~ Fixed_Length_ICM_t
    ()

//  Destroy this  Fixed_Length_ICM_t

  {
   if  (permutation != NULL)
       delete [] permutation;
  }



void  Fixed_Length_ICM_t :: read
    (const char * path)

//  Read the  Fixed_Length_ICM_t  in from the file at  path .

  {
   FILE  * fp;
   char  line [ID_STRING_LEN];
   int  param [NUM_FIXED_LENGTH_PARAMS];
   int  i, n;

   // free memory from previous version
   n = sub_model . size ();
   for  (i = 0;  i < n;  i ++)
     delete sub_model [i];

   fp = File_Open (path, "r");     // Should be "rb"?

   fread (line, sizeof (char), ID_STRING_LEN, fp);    // skip the text header line

   if  (fread (param, sizeof (int), NUM_FIXED_LENGTH_PARAMS, fp)
          != NUM_FIXED_LENGTH_PARAMS)
       {
        fprintf (stderr, "ERROR reading file \"%s\"\n", path);
        exit (EXIT_FAILURE);
       }

   if  (ICM_VERSION_ID != param [0])
       {
        fprintf (stderr, "Bad ICM version = %d  should be %d\n",
                 param [0], ICM_VERSION_ID);
        exit (EXIT_FAILURE);
       }
   if  (ID_STRING_LEN != param [1])
       {
        fprintf (stderr, "Bad ID_STRING_LEN = %d  should be %d\n",
                 param [1], ID_STRING_LEN);
        exit (EXIT_FAILURE);
       }

   length = param [2];
   max_depth = param [3];
   special_position = param [4];
   model_type = ICM_Model_t (param [5]);

   permutation = new int [length];
   fread (permutation, sizeof (int), length, fp);

   for  (i = 0;  i < length;  i ++)
     {
      ICM_t  * p;

      p = new  ICM_t (1, 0, 1);
      p -> Input (fp);

      sub_model . push_back (p);
     }

   return;
  }



double  Fixed_Length_ICM_t :: Score_Window
    (char * w)

//  Return the score of this model on string  w

  {
   static char  * buff = NULL;
   static int  buff_len = 0;
   double  score = 0.0;
   int  i;

   if  (length > buff_len)
       {
        buff = (char *) Safe_realloc (buff, length, __FILE__, __LINE__);
        buff_len = length;
       }

   strncpy (buff, w, length);
   if  (permutation != NULL)
       Permute_String (buff, permutation, length);

   for  (i = 0;  i < length;  i ++)
     {
      if  (buff [i] == '\0')
          {
           fprintf (stderr, "ERROR:  String \"%s\" too short in Score_Window\n",
                    buff);
           exit (EXIT_FAILURE);
          }
      score += sub_model [i] -> Full_Window_Prob (buff, 0);
     }

   return  score;
  }



double  Fixed_Length_ICM_t :: subrange_score
    (char * w, int lo, int hi)

//  Return the score of this model on the portion of window
//  w  between positions  lo  and  hi .   lo  and  hi  are
//  in gap-based coordinates and  w  should point to the beginning
//  of the full-window (i.e., not to where  lo  is).

  {
   static char  * buff = NULL;
   static int  buff_len = 0;
   double  score = 0.0;
   int  i;

   if  (lo < 0 || length < hi || hi < lo)
       {
        fprintf (stderr, "ERROR:  Bad range  lo = %d  hi = %d  in subrange_score\n",
                lo, hi);

        exit (EXIT_FAILURE);
       }

   if  (length > buff_len)
       {
        buff = (char *) Safe_realloc (buff, length, __FILE__, __LINE__);
        buff_len = length;
       }

   strncpy (buff, w, length);
   if  (permutation != NULL)
       Permute_String (buff, permutation, length);

   for  (i = lo;  i < hi;  i ++)
     {
      if  (buff [i] == '\0')
          {
           fprintf (stderr, "ERROR:  String \"%s\" too short in Score_Window\n",
                    buff);
           exit (EXIT_FAILURE);
          }
      score += sub_model [i] -> Full_Window_Prob (buff, 0);
     }

   return  score;
  }



Fixed_Length_ICM_Training_t :: Fixed_Length_ICM_Training_t
    (int len, int md, int sp, int * perm, ICM_Model_t mt)

//  Construct a  Fixed_Length_ICM_Training_t  of length  len  using
//  perm  for the order of bases

  {
   length = len;
   max_depth = md;
   special_position = sp;
   if  (perm == NULL)
       permutation = NULL;
     else
       {
        int  i;

        permutation = new int [len];
        for  (i = 0;  i < len;  i ++)
          permutation [i] = perm [i];
       }
   model_type = mt;
  }



Fixed_Length_ICM_Training_t :: ~ Fixed_Length_ICM_Training_t
       ()

//  Destroy this  Fixed_Length_ICM_Training_t

  {
   int  i, n;

   if  (permutation != NULL)
       delete [] permutation;

   n = int (sub_model . size ());
   for  (i = 0;  i < n;  i ++)
     delete sub_model [i];
  }



void  Fixed_Length_ICM_Training_t :: Output
    (FILE * fp, bool binary_form)

//  Output the contents of this model to  fp .
//  If  binary_form  is true, then do it in binary; otherwise,
//  write an ascii text version.

  {
   int  i;

   Write_Header (fp, binary_form);

   for  (i = 0;  i < length;  i ++)
     sub_model [i] -> Output (fp, binary_form);

   return;
  }



void  Fixed_Length_ICM_Training_t :: Train_Model
    (vector <char *> & data)

//  Calculate the probabilities for this model based on the
//  strings in  data .

  {
   vector <char *>  sub_data;
   ICM_Training_t  * mp;
   int  depth, string_ct;
   int  i, j;
   
   string_ct = int (data . size ());

   if  (permutation != NULL)
       Permute_Data (data, permutation);

   // Make a vector with enough room to hold the substrings
   // of data used to build the sub-models

   for  (j = 0;  j < string_ct;  j ++)
     {
      char  * tmp;

      tmp = (char *) Safe_malloc (length + 1, __FILE__, __LINE__);
      sub_data . push_back (tmp);
     }

   for  (i = 1;  i <= length;  i ++)
     {
      for  (j = 0;  j < string_ct;  j ++)
        {
         strncpy (sub_data [j], data [j], i);
         sub_data [j] [i] = '\0';
        }

      depth = i - 1;
      if  (depth > max_depth)
          depth = max_depth;
      mp = new ICM_Training_t (i, depth, 1);
      mp -> Train_Model (sub_data);
      sub_model . push_back (mp);
     }

   // Free string memory allocated for  sub_data

   for  (j = 0;  j < string_ct;  j ++)
     free (sub_data [j]);

   return;
  }



void  Fixed_Length_ICM_Training_t :: Write_Header
    (FILE * fp, bool binary_form)

//  Send to  fp  the parameter information for this model.
//  binary_form  determines whether the format is binary or
//  Ascii readable (for debugging purposes only)

  {
   int  i;

   if  (! binary_form)
       {
        fprintf (fp, "ver=%.2f  len=%d  depth=%d  special=%d  type=%d",
                 ICM_VERSION_ID / 100.0, length, max_depth, special_position,
                 int (model_type));
        for  (i = 0;  i < length;  i ++)
          {
           if  (i == 0)
               fprintf (fp, "  %d", permutation == NULL ? i : permutation [i]);
             else
               fprintf (fp, ",%d", permutation == NULL ? i : permutation [i]);
          }
        fprintf (fp, "\n");
       }
     else
       {
        char  line [ID_STRING_LEN] = {'\0'};
        char  perm [ID_STRING_LEN] = {'\0'};
        int  param [NUM_FIXED_LENGTH_PARAMS];

        sprintf (line, ">ver=%.2f  len=%d  depth=%d  special=%d  type=%d",
                 ICM_VERSION_ID / 100.0, length, max_depth,
                 special_position, int (model_type));
        for  (i = 0;  i < length;  i ++)
          {
           if  (i == 0)
               sprintf (perm, "  %d", permutation == NULL ? i : permutation [i]);
             else
               sprintf (perm, ",%d", permutation == NULL ? i : permutation [i]);
           strcat (line, perm);
          }
        strcat (line, "\n");

        assert (int (strlen (line)) < ID_STRING_LEN);
        fwrite (line, sizeof (char), ID_STRING_LEN, fp);

        param [0] = ICM_VERSION_ID;
        param [1] = ID_STRING_LEN;
        param [2] = length;
        param [3] = max_depth;
        param [4] = special_position;
        param [5] = int (model_type);

        fwrite (param, sizeof (int), NUM_FIXED_LENGTH_PARAMS, fp);

        if  (permutation != NULL)
             fwrite (permutation, sizeof (int), length, fp);
          else
            {
             int  * tmp;

             tmp = new int [length];
             for  (i = 0;  i < length;  i ++)
               tmp [i] = i;
             fwrite (tmp, sizeof (int), length, fp);
             delete [] tmp;
            }
       }

   return;
  }



void  Count_Char_Pairs
    (int ct [] [ALPHA_SQUARED], char * string, int w, int period)

//  For each complete window of length  w  in  string
//  and for each position  j = 0 .. (w - 2)  add 1 to
//  ct [j] [p]  where  p  is the index of the character
//  pair at positions  j  and  (w - 1)  in the window.
//  The first window starts at the beginning of  string ,
//  and the window advances by  period positions each step
//  as it moves down the string

  {
   int  start, stop, end, i, j, last_char_sub;

   start = 0;
   end = int (strlen (string));

   for  (stop = w - 1;  stop < end;  start += period, stop += period)
     {
      last_char_sub = Subscript (string [stop]);
      for  (i = 0;  i < w - 1;  i ++)
        {
         j = ALPHABET_SIZE * Subscript (string [start + i])
                 + last_char_sub;
         ct [i] [j] ++;
        }
     }

   return;
  }



void  Count_Single_Chars
    (int ct [ALPHABET_SIZE], char * string, int w, int period)

//  For each complete window of length  w  in  string
//  add 1 to  ct [p]  where  p  is the index of the character
//  at the end of the window.
//  The first window starts at the beginning of  string ,
//  and the window advances by  period positions each step
//  as it moves down the string

  {
   int  pos, end, last_char_sub;

   end = int (strlen (string));

   for  (pos = w - 1;  pos < end;  pos += period)
     {
      last_char_sub = Subscript (string [pos]);
      ct [last_char_sub] ++;
     }

   return;
  }



double  Get_Mutual_Info
    (int ct [], int n, int sum)

//  Calculate and return the mutual information for the  n^2
//  counts in  ct  representing frequency of occurrence of
//  pairs of events, where each event has  n  possibilities
//  sum  is the sum of all entries in  ct []

  {
   double  mut_info = 0.0;
   double  * left_prob;
     // probability of each first event in pair
   double  * right_prob;
     // probability of each second event in pair
   int  i, j, k;

   if  (sum == 0)
       return  0.0;

   left_prob = (double * ) Safe_malloc (n * sizeof (double), __FILE__, __LINE__);
   right_prob = (double * ) Safe_malloc (n * sizeof (double), __FILE__, __LINE__);

   for  (i = 0;  i < n;  i ++)
     left_prob [i] = right_prob [i] = 0.0;

   // Calculate  left_prob  and  right_prob
   for  (i = k = 0;  i < n;  i ++)
     for  (j = 0;  j < n;  j ++)
       {
        left_prob [i] += ct [k];
        right_prob [j] += ct [k];
        k ++;
       }
    for  (i = 0;  i < n;  i ++)
      {
       left_prob [i] /= sum;
       right_prob [i] /= sum;
      }

   for  (i = k = 0;  i < n;  i ++)
     for  (j = 0;  j < n;  j ++)
       {
        double  prob = double (ct [k]) / sum;

        if  (prob != 0.0 && left_prob [i] != 0.0 && right_prob [j] != 0.0)
            mut_info += prob * log (prob / (left_prob [i] * right_prob [j]));

        k ++;
       }

   free (left_prob);
   free (right_prob);

   return mut_info;
  }



void  Permute_Data
    (vector <char *> & data, int * perm)

//  Rearrange the characters in each string in  data  according
//  to the permutation in  perm .

  {
   int  len;
   int  i, n;

   n = data . size ();
   if  (n == 0)
       return;

   len = strlen (data [0]);

   for  (i = 0;  i < n;  i ++)
     Permute_String (data [i], perm, len);

   return;
  }



void  Permute_String
    (char * s, int * perm, int n)

//  Rearrange the characters in  s  according
//  to the permutation in  perm .

  {
   static char  * buff = NULL;
   static int  buff_len = 0;
   int  i;

   if  (n > buff_len)
       {
        buff = (char *) Safe_realloc (buff, n, __FILE__, __LINE__);
        buff_len = n;
       }

   for  (i = 0;  i < n;  i ++)
     buff [i] = s [perm [i]];
   strncpy (s, buff, n);

   return;
  }



int  Subscript
    (char ch)

//  Return the subscript equivalent (used in offsets of the
//  model) for character  ch .

  {
   const char  * p;

   p = strchr (ALPHA_STRING, tolower (Filter (ch)));
   //p = strchr (ALPHA_STRING, ch);
   if  (p == NULL)
       {
        fprintf (stderr, "ERROR:  Bad character %c in subscript conversion",
                 ch);
        exit (EXIT_FAILURE);
       }

   return  int (p - ALPHA_STRING);
  }



