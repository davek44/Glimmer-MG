//  A. L. Delcher
//
//  File:  anomaly.cc
//
//  Last Modified:  Fri May 19 17:10:05 EDT 2006
//
//  This program reads the sequence in the first command-line
//  file and then takes the start and end positions specified in
//  the second command-line file and checks for anomalous
//  start/stop codons and frame shifts.


#include  "anomaly.hh"


// Global variables

static bool  Check_Previous_Stop = false;
  // Determines whether to check if codon before start coordinate
  // is a valid stop codon
static bool  Check_Start_Codon = true;
  // Determines whether to check if first codon is a valid start
static char  * Coord_File_Name = NULL;
  // From the second command-line argument
static int  Num_Start_Codons;
  // Number of different start codon patterns
static int  Num_Stop_Codons;
  // Number of different stop codon patterns
static char  * Sequence_File_Name = NULL;
  // From the first command-line argument
static vector <const char *>  Start_Codon;
  // Sequences assumed to be start codons
static vector <const char *>  Stop_Codon;
  // Sequences assumed to be stop codons


int  main
    (int argc, char * argv [])

  {
   FILE  * fp;
   string  Data, hdr;
   char  * Buffer, Line [MAX_LINE], Name [MAX_LINE];
   char  Codon [4] = "tag";
   int  Direction, Frame_Shift;
   long int  Buffer_Len, Gene_Len;
   long int  i, j, Begin, End, Len, Start, Stop;
   int  problem_ct = 0, ok_ct = 0;
   
   try
     {
      Parse_Command_Line (argc, argv);

      Set_Start_And_Stop_Codons ();

      fp = File_Open (Sequence_File_Name, "r");

      Buffer = (char *) Safe_malloc (INIT_SIZE);
      Buffer_Len = INIT_SIZE;

      Fasta_Read (fp, Data, hdr);

      fclose (fp);

      Len = Data . length ();
      Data . insert (Data . begin (), 'x');
        // Put a dummy character at the front of Data so subscripts
        // will start at 1

      fp = File_Open (Coord_File_Name, "r");

      while  (fgets (Line, MAX_LINE, fp) != NULL)
        {
         bool  problem = false;

         if  (sscanf (Line, "%s %ld %ld", Name, & Start, & End) != 3)
             {
              printf ("Bad line:  %s\n...Skipping\n", Line);
              continue;
             }

         if  (Start < End && End - Start <= Len / 2 || Start - End > Len / 2)
             {
              Direction = +1;
              Gene_Len = 1 + End - Start;
              if  (Gene_Len < 0)
                  Gene_Len += Len;

              if  (Buffer_Len < Gene_Len + 1)
                  Buffer = (char *) Safe_realloc (Buffer, 1 + Gene_Len);
              Buffer_Len = 1 + Gene_Len;
              for  (i = 0;  i < Gene_Len;  i ++)
                {
                 if  (Start + i <= Len)
                     j = Start + i;
                   else
                     j = Start + i - Len;
                 Buffer [i] = tolower (Data [j]);
                }
              Buffer [i] = '\0';
             }
           else
             {
              Direction = -1;
              Gene_Len = 1 + Start - End;
              if  (Gene_Len < 0)
                  Gene_Len += Len;

              if  (Buffer_Len < Gene_Len + 1)
                  Buffer = (char *) Safe_realloc (Buffer, 1 + Gene_Len);
              Buffer_Len = 1 + Gene_Len;
              for  (i = 0;  i < Gene_Len;  i ++)
                {
                 if  (Start - i >= 1)
                     j = Start - i;
                   else
                     j = Start - i + Len;
                 Buffer [i] = Complement (tolower (Data [j]));
                }
              Buffer [i] = '\0';
             }

         if  (Check_Previous_Stop)
             {
              if  (Direction == +1)
                  {
                   for  (i = 3;  i > 0;  i --)
                     if  (Start - i < 1)
                         Codon [i] = tolower (Data [Start - i + Len]);
                       else
                         Codon [i] = tolower (Data [Start - i]);
                  }
                else
                  {
                   for  (i = 3;  i > 0;  i --)
                     if  (Start + i > Len)
                         Codon [i] = Complement (tolower (Data [Start + i - Len]));
                       else
                         Codon [i] = Complement (tolower (Data [Start + i]));
                  }
              if  (! Is_Stop_Codon (Codon))
                  {
                   printf ("%-10s %8ld %8ld no stop before start\n",
                                  Name, Start, End);
                   problem = true;
                  }
             }
         if  (Check_Start_Codon && ! Is_Start_Codon (Buffer))
             {
              printf ("%-10s has bad start codon = %.3s\n", Name, Buffer);
              problem = true;
             }
         if  (! Is_Stop_Codon (Buffer + Gene_Len - 3))
             {
              printf ("%-10s has bad stop codon = %s\n", Name, Buffer + Gene_Len - 3);
              problem = true;
              for  (j = Gene_Len;  j < Len;  j += 3)
                {
                 for  (i = 0;  i < 3;  i ++)
                   if  (Direction == +1)
                       {
                        if  (Start + i + j > Len)
                            Codon [i] = tolower (Data [Start + i + j - Len]);
                          else
                            Codon [i] = tolower (Data [Start + i + j]);
                       }
                     else
                       {
                        if  (Start - i - j < 1)
                            Codon [i] = Complement (tolower (Data [Start - i - j + Len]));
                          else
                            Codon [i] = Complement (tolower (Data [Start - i - j]));
                       }
                 if  (Is_Stop_Codon (Codon))
                     break;
                }
              assert (j < Len);
              printf ("           next stop occurs at offset %ld"
                   "  Gene_Len = %ld  diff = %+ld\n",
                   j, Gene_Len, j - Gene_Len + 3);
             }

         Frame_Shift = (Gene_Len % 3);
         if  (Frame_Shift)
             {
              printf ("%-10s %8ld %8ld has %+d frame shift\n",
                              Name, Start, End, Frame_Shift);
              problem = true;

              for  (i = 0;  i < Gene_Len - 3;  i += 3)
                if  (Is_Stop_Codon (Buffer + i))
                    break;
              if  (i < Gene_Len - 3)
                  {
                   Stop = Start + Direction * (i - 1);
                   if  (Stop < 1)
                       Stop += Len;
                   else if  (Stop > Len)
                       Stop -= Len;
                   printf ("   Best prefix is %8ld %8ld   Len = %ld\n",
                                Start, Stop, i);
                  }
                else
                  {
                   printf ("   No stop found in start frame\n");
                   continue;
                  }

              for  (i = Gene_Len - 6;  i >= 0;  i -= 3)
                if  (Is_Stop_Codon (Buffer + i))
                    break;
              i += 3;
              Begin = Start + Direction * i;
              if  (Begin < 1)
                  Begin += Len;
              else if  (Stop > Len)
                  Begin -= Len;
              printf ("   Best suffix is %8ld %8ld   Len = %ld\n",
                           Begin, End, Gene_Len - i - 3);

             }
           else
             {
              for  (i = 0;  i < Gene_Len - 3;  i += 3)
                if  (Is_Stop_Codon (Buffer + i))
                    {
                     printf ("%-10s has stop codon %.3s at offset %ld"
                          "  Gene_Len = %ld  diff = %+ld\n",
                          Name, Buffer + i, i, Gene_Len, Gene_Len - 3 - i);
                     problem = true;
                    }
             }
         if  (problem)
             problem_ct ++;
           else
             ok_ct ++;
        }

      fprintf (stderr, "     OK orfs = %7d\n", ok_ct);
      fprintf (stderr, "Problem orfs = %7d\n", problem_ct);
     }
   catch (std :: exception & e)
     {
      cerr << "** Standard Exception **" << endl;
      cerr << e << endl;
      exit (EXIT_FAILURE);
     }

   return  0;
  }



static bool  Is_Start_Codon
    (const char * p)

//  Return true iff the first 3 characters of p match a
//  string in global  Start_Codon .

  {
   int  i;

   for  (i = 0;  i < Num_Start_Codons;  i ++)
     if  (strncmp (p, Start_Codon [i], 3) == 0)
         return  true;

   return  false;
  }



static bool  Is_Stop_Codon
    (const char * p)

//  Return true iff the first 3 characters of p match a
//  string in global  Stop_Codon .

  {
   int  i;

   for  (i = 0;  i < Num_Stop_Codons;  i ++)
     if  (strncmp (p, Stop_Codon [i], 3) == 0)
         return  true;

   return  false;
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

   while  (! errflg && ((ch = getopt (argc, argv, "A:stZ:")) != EOF))
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

        case  's' :
          Check_Start_Codon = FALSE;
          break;

        case  't' :
          Check_Previous_Stop = TRUE;
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
   Coord_File_Name = argv [optind ++];

   return;
  }



static void  Set_Start_And_Stop_Codons
    (void)

//  Set globals  Start_Codon  and  Stop_Codon  to the sequences
//  that are allowed to be start and stop codons for genes.

  {
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

   Num_Start_Codons = Start_Codon . size ();
   Num_Stop_Codons = Stop_Codon . size ();

   return;
  }



static void  Usage
    (void)

//  Print to stderr description of options and command line for
//  this program.

  {
   fprintf (stderr,
       "USAGE:  anomaly [options] <sequence-file> <coord-file>\n"
       "\n"
       "Read DNA sequence in <sequence-file> and for each region specified\n"
       "by the coordinates in <coord-file>, check whether the region\n"
       "represents a normal gene, i.e., it begins with a start codon, ends\n"
       "with a stop codon, and has no frame shifts.\n"
       "Output goes to standard output.\n"
       "\n"
       "Options:\n"
       " -A <codon-list>\n"
       "    Use comma-separated list of codons as start codons\n"
       "    Sample format:  -A atg,gtg\n"
       "    Default start codons are atg,gtg,ttg\n"
       " -s\n"
       "    Omit the check that the first codon is a start codon.\n"
       " -t\n"
       "    Check whether the codon preceding the start coordinate position\n"
       "    is a stop codon.  This is useful if the coordinates represent\n"
       "    the entire region between stop codons.\n"
       " -Z <codon-list>\n"
       "    Use comma-separated list of codons as stop codons\n"
       "    Sample format:  -Z tag,tga,taa\n"
       "\n");

   return;
  }



