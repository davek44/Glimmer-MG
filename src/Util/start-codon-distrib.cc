//  A. L. Delcher
//
//  File:  start-codon-distrib.cc
//
//  Last Modified:  Tue May 23 10:30:35 EDT 2006
//
//  This program reads the sequence in the file named as the
//  first command-line argument and then the list of gene coordinates
//  in the second command-line file.  It then counts the number of
//  different start codons for each of the genes in the second file.
//  If the genome is circular, the direction of each gene is assumed
//  to be the shorter direction around the circle; otherwise the
//  order of the coordinates determines the direction.  Alternately,
//  the 4th field of the coordinate file can specify the direction
//  of the gene.


#include  "start-codon-distrib.hh"


// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

static bool  Comma3_Output = false;
  // If set true by the -3 option, then output only a comma separated
  // list (no spaces) of atg, gtg, ttg start proportions, in that order
static char  * Coord_Info = NULL;
  // Name of file with list of coordinates (or a single coordinate
  // specification).
static bool  Is_Circular = true;
  // Determines whether the input sequence is regarded as a circular
  // genome.
static char  * Sequence_File_Name = NULL;
  // Name of file with input sequence
static bool  Use_Direction = false;
  // If set true (by -d option) then use 4th field of coords
  // to determine direction in a circular genome.



int  main
    (int argc, char * argv [])

  {
   FILE  * fp;
   vector <Count_Entry_t>  entry;
   string  sequence, hdr;
   char  line [MAX_LINE], tag [MAX_LINE];
   char  codon [4];
   long int  seq_len, start, end;
   int  i, n, total_entries = 0;

   Verbose = 0;

   Parse_Command_Line (argc, argv);

   fp = File_Open (Sequence_File_Name, "r");

   if  (! Fasta_Read (fp, sequence, hdr))
       {
        sprintf (Clean_Exit_Msg_Line, "ERROR:  Failed to open file %s",
             Sequence_File_Name);
        Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
       }
   fclose (fp);

   seq_len = sequence . length ();

   if  (strcmp (Coord_Info, "-") == 0)
       fp = stdin;
     else
       fp = File_Open (Coord_Info, "r");

   codon [3] = '\0';
   while  (fgets (line, MAX_LINE, fp) != NULL)
     {
      int  dir;

      if  (Use_Direction)
          {
           if  (sscanf (line, "%s %ld %ld %d", tag, & start, & end, & dir) != 4)
               {
                fprintf (stderr, "ERROR:  Skipped following coord line\n");
                fputs (line, stderr);
                continue;
               }
          }
        else
          {
           if  (sscanf (line, "%s %ld %ld", tag, & start, & end) != 3)
               {
                fprintf (stderr, "ERROR:  Skipped following coord line\n");
                fputs (line, stderr);
                continue;
               }
           if  ((start < end && (! Is_Circular || end - start <= seq_len / 2))
                    || (Is_Circular && start - end > seq_len / 2))
               dir = +1;
             else
               dir = -1;
          }

      if  (dir > 0)
          {
           for  (i = 0;  i < 3;  i ++)
             codon [i] = tolower (sequence [Seq_Sub (start + i, seq_len)]);
           Incr (entry, codon);
          }
        else
          {
           for  (i = 0;  i < 3;  i ++)
             codon [i] = Complement (tolower (sequence [Seq_Sub (start - i, seq_len)]));
           Incr (entry, codon);
          }
      total_entries ++;
     }

   if  (Comma3_Output)
       {
        if  (total_entries == 0)
            total_entries = 1;
        n = entry . size ();
        for  (i = 0;  i < n;  i ++)
          if  (strcmp (entry [i] . codon, "atg") == 0)
              {
               printf ("%.3f", double (entry [i] . count) / total_entries);
               break;
              }
        if  (i == n)
            printf ("%.3f", 0.0);
        for  (i = 0;  i < n;  i ++)
          if  (strcmp (entry [i] . codon, "gtg") == 0)
              {
               printf (",%.3f", double (entry [i] . count) / total_entries);
               break;
              }
        if  (i == n)
            printf (",%.3f", 0.0);
        for  (i = 0;  i < n;  i ++)
          if  (strcmp (entry [i] . codon, "ttg") == 0)
              {
               printf (",%.3f\n", double (entry [i] . count) / total_entries);
               break;
              }
        if  (i == n)
            printf (",%.3f\n", 0.0);
       }
     else
       {
        sort (entry . begin (), entry . end (), Count_Entry_Cmp);
        n = entry . size ();
        for  (i = 0;  i < n;  i ++)
          {
           printf (" %s   %6d  %5.1f%%\n", entry [i] . codon, entry [i] . count,
                Percent (entry [i] . count, total_entries));
          }
        printf ("Total: %6d\n", total_entries);
       }

   return  0;
  }



static bool  Count_Entry_Cmp
    (Count_Entry_t const & a, Count_Entry_t const & b)

//  Compare  a  and  b  first by  count  (for descending order sort),
//  or if equal, by  codon  alphabetically.
    
  {
   if  (a . count > b . count)
       return  true;

   return  (a . count == b . count && strcmp (a . codon, b . codon) < 0);
  }



static void  Incr
    (vector <Count_Entry_t> & entry, const char s [4])

//  Search for string  s  in  entry .  If found, increment the
//  corresponding count; otherwise, create a new entry and set
//  its count to 1.

  {
   Count_Entry_t  e;
   int  i, n;

   n = entry . size ();
   for  (i = 0;  i < n;  i ++)
     if  (strcmp (entry [i] . codon, s) == 0)
         {
          entry [i] . count ++;
          return;
         }

   strcpy (e . codon, s);
   e . count = 1;
   entry . push_back (e);

   return;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   bool  errflg = false;
   int  ch;

   optarg = NULL;

#if  ALLOW_LONG_OPTIONS
   int  option_index = 0;
   static struct option  long_options [] = {
        {"dir", 0, 0, 'd'},
        {"help", 0, 0, 'h'},
        {"minlen", 1, 0, 'l'},
        {"nostart", 0, 0, 's'},
        {"nostop", 0, 0, 't'},
        {"nowrap", 0, 0, 'w'},
        {"3comma", 0, 0, '3'},
        {0, 0, 0, 0}
      };

   while  (! errflg
        && ((ch = getopt_long (argc, argv, "dhw3",
                     long_options, & option_index)) != EOF))
#else
   while  (! errflg
        && ((ch = getopt (argc, argv, "dhw3")) != EOF))
#endif

     switch  (ch)
       {
        case  'd' :
          Use_Direction = true;
          break;

        case  'h' :
          Usage ();
          exit (EXIT_SUCCESS);

        case  'w' :
          Is_Circular = false;
          break;

        case  '3' :
          Comma3_Output = true;
          break;

        default :
          errflg = true;
       }

   if  (errflg || optind > argc - 2)
       {
        Usage ();
        exit (EXIT_FAILURE);
       }

   Sequence_File_Name = argv [optind ++];
   Coord_Info = argv [optind ++];

   return;
  }



static long int  Seq_Sub
    (long int sub, long int seq_len)

//  Return the subscript (in zero-based coordinates) of position
//   sub  (which is in 1-based coordinates) allowing for circular
//  wraparounds (off either end) of a sequence of length  seq_len .

  {
   sub --;

   while  (sub < 0)
     sub += seq_len;

   while  (seq_len - 1 <= sub)
     sub -= seq_len;

   return  sub;
  }



static void  Usage
    (void)

//  Print to stderr description of options and command line for
//  this program.

  {
   fprintf (stderr,
       "USAGE:  start-codon-distrib [options] <sequence-file> <coords>\n"
       "\n"
       "Read fasta-format <sequence-file> and count the number of\n"
       "different start codons for the genes specified in <coords>.\n"
       "By default, <coords> is the name of a file containing lines of the form\n"
       "  <tag>  <start>  <stop>  [<frame>] ...\n"
       "Coordinates are inclusive counting from 1, e.g., \"1 3\"\n"
       "represents the 1st 3 characters of the sequence.\n"
       "Output goes to stdout.\n"
       "\n"
       "Options:\n"
       " -d\n"
       " --dir\n"
       "    Use the 4th column of each input line to specify the direction\n"
       "    of the sequence.  Positive is forward, negative is reverse\n"
       "    The input sequence is assumed to be circular\n"
       " -h\n"
       "    Print this message\n"
       " -w\n"
       " --nowrap\n"
       "    Use the actual input coordinates without any wraparound\n"
       "    that would be needed by a circular genome.  Without this\n"
       "    option, the output sequence is the shorter of the two ways\n"
       "    around the circle.  Use the -d option to specify direction\n"
       "    explicitly.\n"
       " -3\n"
       " --3comma\n"
       "    output only a comma separated list (no spaces) of atg, gtg, ttg\n"
       "start proportions, in that order\n"
       "\n");

   return;
  }



