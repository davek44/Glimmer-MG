//    Programmer:  Arthur L. Delcher
//          File:  score-fixed.cc
//  Last Updated:  Mon Jun  7 10:35:06 EDT 2004
//
//  Compute scores for a set of fixed length strings using
//  the model in the file on the command line.


#include  "score-fixed.hh"


static char  * Pos_Model_Path;
  // Name of file containing the positive model
static char  * Neg_Model_Path;
  // Name of file containing the negative model
static bool  Simple_Output = false;
  // If true, output is just string number (starting with 0)
  // and either 1 or -1 indicating which model scores higher
static bool  Use_Neg_ICM_Model = false;
  // If true then negative model is usual (streaming) ICM
  // instead of fixed-length ICM
static bool  Use_Null_Neg_Model = false;
  // If true then negative model is ignored and value
  // from it is automatically zero


//**ALD  Gets rid of make undefined reference error
int  Unused = Filter ('a');



int  main
    (int argc, char * argv [])

  {
   Fixed_Length_ICM_t  pos_model;
   ICM_t  neg_icm_model;
   Fixed_Length_ICM_t  neg_fixed_model;
   char  * string = NULL, * tag = NULL;
   long int  string_size = 0, tag_size = 0;
   int  string_num = 0;

   Parse_Command_Line (argc, argv);

   pos_model . read (Pos_Model_Path);
   fprintf (stderr, "pos model  len = %d  special = %d  type = %d\n",
            pos_model . getModelLength (),
            pos_model . getSpecialPosition (),
            pos_model . getModelType ());
   if  (Use_Null_Neg_Model)
       fprintf (stderr, "Using null negative model\n");
   else if  (Use_Neg_ICM_Model)
       neg_icm_model . Read (Neg_Model_Path);
     else
       {
        neg_fixed_model . read (Neg_Model_Path);
        fprintf (stderr, "neg model  len = %d  special = %d  type = %d\n",
                 neg_fixed_model . getModelLength (),
                 neg_fixed_model . getSpecialPosition (),
                 neg_fixed_model . getModelType ());
       }

   while  (Read_String (stdin, string, string_size, tag, tag_size))
     {
      double  pos_score, neg_score;
      double  avg_pos_score, avg_neg_score;
      int  len;

      string_num ++;
      len = strlen (string);

      pos_score = pos_model . score (string);
      if  (Use_Null_Neg_Model)
          neg_score = 0.0;
      else if  (Use_Neg_ICM_Model)
          neg_score = neg_icm_model . Score_String (string, strlen (string), 1);
        else
          neg_score = neg_fixed_model . score (string);

      avg_pos_score = pos_score / len;
      avg_neg_score = neg_score / len;

      if  (Simple_Output)
          printf ("%6d %3d\n", string_num - 1,
                  pos_score >= neg_score ? 1 : -1);
        else
          printf ("%5d:  %10.4f %9.5f   %10.4f %9.5f   %9.5f\n",
                  string_num, pos_score, avg_pos_score,
                  neg_score, avg_neg_score, avg_pos_score - avg_neg_score );
     }

   return  0;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   bool  errflg = false;
   int  ch;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "hINs")) != EOF))
     switch  (ch)
       {
        case  'h' :
          errflg = true;
          break;

        case  'I' :
          Use_Neg_ICM_Model = true;
          break;

        case  'N' :
          Use_Null_Neg_Model = true;
          break;

        case  's' :
          Simple_Output = true;
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = TRUE;
       }

   if  (errflg || (Use_Null_Neg_Model && optind > argc - 1)
          || (! Use_Null_Neg_Model && optind != argc - 2))
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   Pos_Model_Path = argv [optind ++];
   if  (! Use_Null_Neg_Model)
       Neg_Model_Path = argv [optind ++];

   return;
  }



int  Read_String
    (FILE * fp, char * & s, long int & s_size, char * & tag,
     long int & tag_size)

//  Read next string from  fp  (assuming FASTA format) into  s [0 .. ]
//  which has  s_size  characters.  Allocate extra memory if needed
//  and adjust  s_size  accordingly.  Return  TRUE  if successful,  FALSE
//  otherwise (e.g., EOF).  Put FASTA header line into  tag [0 .. ]
//  (and adjust  tag_size  if needed).

  {
   int  ch, ct;

   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     ;

   if  (ch == EOF)
       return  FALSE;

   ct = 0;
   while  ((ch = fgetc (fp)) != EOF && ch != '\n' && isspace (ch))
     ;
   if  (ch == EOF)
       return  FALSE;
   if  (ch != '\n' && ! isspace (ch))
       ungetc (ch, fp);
   while  ((ch = fgetc (fp)) != EOF && ch != '\n')
     {
      if  (ct >= tag_size - 1)
          {
           tag_size += INCR_SIZE;
           tag = (char *) Safe_realloc (tag, tag_size);
          }
      tag [ct ++] = char (ch);
     }
   tag [ct ++] = '\0';

   ct = 0;
   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     {
      if  (isspace (ch))
          continue;

      if  (ct >= s_size - 1)
          {
           s_size += INCR_SIZE;
           s = (char *) Safe_realloc (s, s_size);
          }
      s [ct ++] = char (ch);
     }
   s [ct ++] = '\0';

   if  (ch == '>')
       ungetc (ch, fp);

   return  TRUE;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s [options]  <pos-model>  <neg-model>  < input-file\n"
           "\n"
           "Read sequences from  stdin  and score them using fixed-length\n"
           "model in file  <model> .  Output scores to  stdout.\n"
           "Output columns are:  sequence number, positive total score,\n"
           "  positive score per base, negative total score,\n"
           "  negative score per base, delta positive/negative per-base scores.\n"
           "\n"
           "Options:\n"
           " -h        Print this message\n"
           " -I        Negative model is regular ICM, not fixed-length ICM\n"
           " -N        Use NULL negative model, i.e., constant zero\n"
           " -s        Output simple format of string num and 1 or -1\n"
           "\n",
           command);

   return;
  }



