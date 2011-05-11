//  Programmers:   Arthur L. Delcher
//                 Doug Harmon
//  File:          build-icm.cc
//  Last Updated:  Mon Jun 12 15:34:00 EDT 2006
//                
//  This program reads (from  stdin ) a set of strings in
//  multi-fasta format.  It then builds and outputs to  stdout
//  an interpolated context model (ICM) that matches the input.
//
//  Copyright (c) 2006 University of Maryland Center for Bioinformatics
//  & Computational Biology


#include  "build-icm.hh"


static int  Genbank_Xlate_Code = 0;
  // Holds the Genbank translation table number that determines
  // stop codons and codon translation.
static int  Model_Depth = DEFAULT_MODEL_DEPTH;
  // Maximum number of positions to use in Markov context
static int  Model_Len = DEFAULT_MODEL_LEN;
  // Width of Markov context and character to be predicted
static int  Model_Periodicity = DEFAULT_PERIODICITY;
  // Number of different models to cycle through
static char  * Output_Filename = NULL;
  // Name of file to which the ICM is written
static bool  Print_Binary = true;
  // Print model as a binary file iff this is true; otherwise print
  // as text file
static bool  Reverse_Strings = false;
  // If true, then use the reverse (not the reverse complement)
  // of input strings to train the model.
static bool  Skip_In_Frame_Stop_Strings = false;
  // If true then ignore any input string with an in-frame stop codon
static vector <const char *>  Stop_Codon;
  // Sequences assumed to be stop codons
static vector <char *>  Training_Data;
  // Holds the training strings




//**ALD  Gets rid of make undefined reference error
int  Unused = Filter ('a');



int  main
    (int argc, char **argv)
  {
   FILE  * output_fp;
   int  string_ct;
     // Number of strings read from training file


   Parse_Command_Line (argc, argv);

   if  (strcmp (Output_Filename, "-") == 0)
       output_fp = stdout;
   else if  (Print_Binary)
       output_fp = File_Open (Output_Filename, "wb");
     else
       output_fp = File_Open (Output_Filename, "w");

   // create the model
   ICM_Training_t  model (Model_Len, Model_Depth, Model_Periodicity);

   Read_Training_Data (stdin);
   string_ct = Training_Data . size ();
   if  (string_ct == 0)
       {
        fprintf (stderr, "ERROR:  Cannot create model--no input data\n");
        fclose (output_fp);
        exit (EXIT_FAILURE);
       }

   if  (Skip_In_Frame_Stop_Strings)
       {
        bool  skip;
        int  i, j, k, s, len, ct = 0;

        Set_Stop_Codons ();

        int  num_stops = Stop_Codon . size ();

        for  (i = k = 0;  i < string_ct;  i ++)
          {
           skip = false;

           // Assuming data has been converted to lower-case if needed

           len = strlen (Training_Data [i]);

           for  (j = 0;  j < len - 2 && ! skip;  j += 3)
             for  (s = 0;  s < num_stops && ! skip;  s ++)
               skip = (strncmp (Training_Data [i] + j, Stop_Codon [s], 3) == 0);

           if  (skip)
               ct ++;
             else
               Training_Data [k ++] = Training_Data [i];
          }

        fprintf (stderr,
                 "Skipped %d strings with in-frame stops of %d total strings\n",
                 ct, string_ct);
        Training_Data . resize (k);
       }

   if  (Reverse_Strings)
       {
        int  i, n;

        n = Training_Data . size ();
        for  (i = 0;  i < n;  i ++)
          Reverse_String (Training_Data [i]);
       }

   model . Train_Model (Training_Data);

   model . Output (output_fp, Print_Binary);

   fclose (output_fp);

   return 0;
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
        {"depth", 1, 0, 'd'},
        {"no_stops", 0, 0, 'F'},
        {"help", 0, 0, 'h'},
        {"period", 1, 0, 'p'},
        {"reverse", 0, 0, 'r'},
        {"text", 0, 0, 't'},
        {"verbose", 1, 0, 'v'},
        {"width", 1, 0, 'w'},
        {"trans_table", 1, 0, 'z'},
        {"stop_codons", 1, 0, 'Z'},
        {0, 0, 0, 0}
      };

   while  (! errflg && ((ch = getopt_long (argc, argv,
        "d:Fhp:rtv:w:z:Z:",
        long_options, & option_index)) != EOF))
#else
   while  (! errflg && ((ch = getopt (argc, argv,
        "d:Fhp:rtv:w:z:Z:")) != EOF))
#endif

     switch  (ch)
       {
        case  'd' :
          Model_Depth = int (strtol (optarg, & p, 10));
          if  (p == optarg || Model_Depth <= 0)
              {
               fprintf (stderr, "Bad model depth value \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
          break;
          
        case  'F' :
          Skip_In_Frame_Stop_Strings = true;
          break;

        case  'h' :
          errflg = TRUE;
          break;

        case  'p' :
          Model_Periodicity = int (strtol (optarg, & p, 10));
          if  (p == optarg || Model_Periodicity <= 0)
              {
               fprintf (stderr, "Bad model period value \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
          break;
          
        case  'r' :
          Reverse_Strings = true;
          break;
          
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
          
        case  'w' :
          Model_Len = int (strtol (optarg, & p, 10));
          if  (p == optarg || Model_Len <= 0)
              {
               fprintf (stderr, "Bad model length value \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
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
          errflg = TRUE;
       }

   if  (errflg || optind != argc - 1)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   Output_Filename = argv [optind ++];

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
      Make_Lower_Case (p);
      Training_Data . push_back (p);
     }

   return;
  }



static void  Set_Stop_Codons
    (void)

//  Set global  Stop_Codon  to the sequences
//  that are allowed to be stop codons for genes, if
//  not already set.

  {
   int  i, n;

   if  (Stop_Codon . size () == 0)
       {
        n = sizeof (DEFAULT_STOP_CODON) / sizeof (char *);
        for  (i = 0;  i < n;  i ++)
          Stop_Codon . push_back (DEFAULT_STOP_CODON [i]);
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
           "USAGE:  build-icm [options] output_file < input-file\n"
           "\n"
           "Read sequences from standard input and output to  output-file\n"
           "the interpolated context model built from them.\n"
           "Input also can be piped into the program, e.g.,\n"
           "  cat abc.in | build-icm xyz.icm\n"
           "If <output-file> is \"-\", then output goes to standard output\n"
           "\n"
           "Options:\n"
           " -d <num>\n"
           "    Set depth of model to <num>\n"
           " -F\n"
           "    Ignore input strings with in-frame stop codons\n"
           " -h\n"
           "    Print this message\n"
           " -p <num>\n"
           "    Set period of model to <num>\n"
           " -r\n"
           "    Use the reverse of input strings to build the model\n"
           " -t\n"
           "    Output model as text (for debugging only)\n"
           " -v <num>\n"
           "    Set verbose level; higher is more diagnostic printouts\n"
           " -w <num>\n"
           "    Set length of model window to <num>\n"
           "\n");

   return;
  }



