//    Programmer:  Arthur L. Delcher
//          File:  build-fixed.cc
//  Last Updated:  Fri Jun  4 16:31:05 EDT 2004
//                
//  This program reads (from  stdin ) a set of fixed_length strings in
//  multi-fasta format.  It then builds and outputs to  stdout
//  a fixed-length interpolated context model (ICM) that matches the input.
//
//  Copyright (c) 2006 University of Maryland Center for Bioinformatics
//  & Computational Biology


#include  "build-fixed.hh"


static FILE  * Index_File_fp = NULL;
  // File containing a list of subscripts of strings to train model
static int  Model_Depth = DEFAULT_MODEL_DEPTH;
  // Maximum number of positions to use in Markov context
static int  Model_Len = DEFAULT_MODEL_LEN;
  // Width of Markov context and character to be predicted
static ICM_Model_t  Model_Type = UNKNOWN_TYPE;
  // Type of model
static int  * Permutation = NULL;
  // Describes how to re-order the characters before building the model
static int  Permutation_Len = 0;
  // Length of above permutation; must match length of input strings
static bool  Print_Binary = true;
  // Print model as a binary file iff this is true; otherwise print
  // as text file
static int  Special_Position = -1;
  // Designated position in model, e.g., for splice junction
static vector <char *>  Training_Data;
  // Holds the training strings


//**ALD  Gets rid of make undefined reference error
int  Unused = Filter ('a');



int  main
    (int argc, char * argv [])
  {
   int  string_ct;
     // Number of strings read from training file
   int  i;


   Parse_Command_Line (argc, argv);

   Read_Training_Data (stdin);
   string_ct = Training_Data . size ();

   if  (string_ct <= 0)
       {
        fprintf (stderr, "ERROR:  No strings read to train model\n");
        exit (EXIT_FAILURE);
       }

   if  (Index_File_fp != NULL)
       {
        // Read the file of subscripts, make a list of the strings
        // they refer to and use that for training.

        vector <char *>  list;
        int  sub;

        while  (fscanf (Index_File_fp, "%d", & sub) == 1)
          list . push_back (Training_Data [sub]);

        Training_Data = list;
        string_ct = Training_Data . size ();
       }

   Model_Len = strlen (Training_Data [0]);
   for  (i = 1;  i < string_ct;  i ++)
     if  (int (strlen (Training_Data [i])) != Model_Len)
         {
          fprintf (stderr, "ERROR:  String #%d has length = %d\n",
                   i, int (strlen (Training_Data [i])));
          fprintf (stderr, "        different from string #0 length = %d\n",
                   Model_Len);
          exit (EXIT_FAILURE);
         }
   if  (Permutation != NULL && Permutation_Len != Model_Len)
       {
        fprintf (stderr, "ERROR:  Permutation len = %d  string_len = %d\n",
                 Permutation_Len, Model_Len);
        exit (EXIT_FAILURE);
       }

   // create the model
   if  (Special_Position > Model_Len)
       {
        fprintf (stderr, "ERROR:  Bad special position = %d\n",
                 Special_Position);
       }
   Fixed_Length_ICM_Training_t  model (Model_Len, Model_Depth, Special_Position,
                                       Permutation, Model_Type);

   model . Train_Model (Training_Data);

   model . Output (stdout, Print_Binary);

   return 0;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   char  * p;
   int  ch, errflg = FALSE;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "bd:hi:p:s:tv:")) != EOF))
     switch  (ch)
       {
        case  'b' :
          Print_Binary = true;
          break;
          
        case  'd' :
          Model_Depth = int (strtol (optarg, & p, 10));
          if  (p == optarg || Model_Depth <= 0)
              {
               fprintf (stderr, "Bad model depth value \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
          break;
          
        case  'h' :
          errflg = TRUE;
          break;

        case  'i' :
          Index_File_fp = File_Open (optarg, "r");
          break;

        case  'p' :
          {
           vector <int>  perm;
           int  i, j, n;

           for  (p = strtok (optarg, ", ");  p != NULL;  p = strtok (NULL, ", "))
             perm . push_back (atoi (p));
           n = perm . size ();
           Permutation = (int *) Safe_calloc (n, sizeof (int), __FILE__,
                                      __LINE__);
           for  (i = 0;  i < n;  i ++)
             if  (Permutation [perm [i]] == 0)
                 Permutation [perm [i]] = 1;
               else
                 {
                  fprintf (stderr, "ERROR:  Illegal permutation\n");
                  for  (j = 0;  j <= i;  j ++)
                    fprintf (stderr, " %d", perm [j]);
                  fprintf (stderr, " <-- duplicate\n");
                  exit (EXIT_FAILURE);
                 }
           for  (i = 0;  i < n;  i ++)
             if  (Permutation [i] == 0)
                 {
                  fprintf (stderr, "ERROR:  Illegal permutation--missing %d\n", i);
                  exit (EXIT_FAILURE);
                 }
           for  (i = 0;  i < n;  i ++)
             Permutation [i] = perm [i];
           Permutation_Len = n;
          }
          break;

        case  's' :
          Special_Position = strtol (optarg, NULL, 10);
          break;

#if  0    // ALD removed on 22 May 2006
        case  'T' :
          Model_Type = ICM_Model_t (strtol (optarg, NULL, 10));
          break;
#endif

        case  't' :
          Print_Binary = false;
          break;
          
        case  'v' :
          Verbose = int (strtol (optarg, & p, 10));
          if  (p == optarg)
              {
               fprintf (stderr, "Bad verbose value \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
          break;
          
        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc - 0)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   return;
  }



static int  Read_String
    (FILE * fp, char * & s, long int & s_size, char * & tag, long int & tag_size)

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



static void  Read_Training_Data
    (FILE  * fp)

// Read in training strings from  fp .  Format is multifasta, i.e., for
// each string a header line (starting with '>') followed by arbitrarily
// many data lines.  Save strings in global  Training_Data

  {
   char  * string = NULL, * tag = NULL;
   char  * p;
   long int  string_size = 0, tag_size = 0;

   while  (Read_String (fp, string, string_size, tag, tag_size))
     {
      p = strdup (string);
      Training_Data . push_back (p);
     }

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s [<options>]  < <input-file>  > <output-file>\n"
           "\n"
           "Read sequences from  stdin  and output to  stdout \n"
           "the fixed-length interpolated context model built from them\n"
           "\n"
           "Options:\n"
           " -d <num>  Set depth of model to <num>\n"
           " -h        Print this message\n"
           " -i <fn>   Train using strings specified by subscripts in file\n"
           "           named <fn>\n"
           " -p n1,n2,...,nk  Permutation describing re-ordering of\n"
           "           character positions of input string to build model\n"
           " -s <num>  Specify special position in model\n"
           " -t        Output model as text (for debugging only)\n"
           " -v <num>  Set verbose level; higher is more diagnostic printouts\n"
           "\n",
           command);

   return;
  }



