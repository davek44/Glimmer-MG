//  A. L. Delcher
//
//  File:  fasta.cc
//
//  Last Modified:  23 October 2003
//
//  Routines to manipulate FASTA format files


#include  "fasta.hh"
#include  "kelley.hh"


void  Fasta_Print
    (FILE * fp, const char * s, const char * hdr, int fasta_width)

//  Print string  s  in fasta format to  fp .  Put string  hdr
//  on header line, unless it's  NULL  in which case do not print
//  a header line at all.  Print at most  fasta_width  characters per
//  line.

  {
   int  ct = 0;

   if  (hdr != NULL)
       fprintf (fp, ">%s\n", hdr);

   while  (* s != '\0')
     {
      if  (ct == fasta_width)
          {
           fputc ('\n', fp);
           ct = 0;
          }
      fputc (* s, fp);
      s ++;
      ct ++;
     }

   fputc ('\n', fp);

   return;
  }



void  Fasta_Print_N
    (FILE * fp, const char * s, int n, const char * hdr, int fasta_width)

//  Print first  n  bytes of  string  s  in fasta format to  fp .
//  Put string  hdr  on header line.  Print at most  fasta_width
//  characters per line.

  {
   int  ct = 0, i;

   if  (hdr != NULL)
       fprintf (fp, ">%s\n", hdr);

   for  (i = 0;  i < n;  i ++)
     {
      if  (ct == fasta_width)
          {
           fputc ('\n', fp);
           ct = 0;
          }
      fputc (s [i], fp);
      ct ++;
     }

   fputc ('\n', fp);

   return;
  }



void  Fasta_Print_Skip
    (FILE * fp, const char * s, const char * skip, const char * hdr,
     int fasta_width)

//  Print string  s  in fasta format to  fp  but omit any characters
//  that occur in string  skip .  Put string  hdr
//  on header line, unless it's  NULL  in which case do not print
//  a header line at all.  Print at most  fasta_width  characters per
//  line.

  {
   int  ct = 0;

   if  (hdr != NULL)
       fprintf (fp, ">%s\n", hdr);

   while  (* s != '\0')
     {
      if  (strchr (skip, * s) == NULL)
          {
           if  (ct == fasta_width)
               {
                fputc ('\n', fp);
                ct = 0;
               }
           fputc (* s, fp);
           ct ++;
          }
      s ++;
     }

   fputc ('\n', fp);

   return;
  }


bool  Fasta_Qual_Vec_Read
(FILE * fp, vector<int> & q, string & hdr)

//  Read next fasta-like-format quality value sequence from
//  file  fp  (which must already be open) into vector  q.
//  Put the faster header line (without the '>' and trailing spaces) into
//  string  hdr .  Return true if successfully read; false, otherwise.

  {
   bool  have_value;
   int  ch, val;

   q . clear();
   hdr . erase ();

   // skip till next '>' if necessary
   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     ;

   if  (ch == EOF)
       return  false;

   // skip spaces if any
   while  ((ch = fgetc (fp)) != EOF && ch == ' ')
     ;
   if  (ch == EOF)
       return  false;
   ungetc (ch, fp);

   // put rest of line into  hdr
   while  ((ch = fgetc (fp)) != EOF && ch != '\n')
     hdr . push_back (char (ch));

   // put all numbers up till next '>' into  q
   have_value = false;
   val = 0;
   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     {
      if  (isspace (ch))
          {
           if  (have_value)
		q . push_back (val);
           have_value = false;
           val = 0;
          }
      else if  (isdigit (ch))
          {
           have_value = true;
           val = 10 * val + ch - '0';
          }
     }

   if  (ch == '>')
       ungetc (ch, fp);

   return  true;
  }


bool  Fasta_Qual_Read
    (FILE * fp, string & q, string & hdr)

//  Read next fasta-like-format quality value sequence from
//  file  fp  (which must already be open) into string  q 
//  (encoded by adding the quality value to the  QUALITY_OFFSET  value).
//  Put the faster header line (without the '>' and trailing spaces) into
//  string  hdr .  Return  true  if a string is successfully,
//  read; false, otherwise.

  {
   bool  have_value;
   int  ch, val;

   q . erase ();
   hdr . erase ();

   // skip till next '>' if necessary
   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     ;

   if  (ch == EOF)
       return  false;

   // skip spaces if any
   while  ((ch = fgetc (fp)) != EOF && ch == ' ')
     ;
   if  (ch == EOF)
       return  false;
   ungetc (ch, fp);

   // put rest of line into  hdr
   while  ((ch = fgetc (fp)) != EOF && ch != '\n')
     hdr . push_back (char (ch));

   // put all numbers up till next '>' into  q
   have_value = false;
   val = 0;
   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     {
      if  (isspace (ch))
          {
           if  (have_value)
               q . push_back (char (val + QUALITY_OFFSET));
           have_value = false;
           val = 0;
          }
      else if  (isdigit (ch))
          {
           have_value = true;
           val = 10 * val + ch - '0';
          }
     }

   if  (ch == '>')
       ungetc (ch, fp);

   return  true;
  }



bool  Fasta_Read
    (FILE * fp, string & s, string & hdr)

//  Read next fasta-format string from file  fp  (which must
//  already be open) into string  s .  Put the faster
//  header line (without the '>' and trailing spaces) into
//  string  hdr .  Return  true  if a string is successfully,
//  read; false, otherwise.

  {
   int  ch;

   s . erase ();
   hdr . erase ();

   // skip till next '>' if necessary
   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     ;

   if  (ch == EOF)
       return  false;

   // skip spaces if any
   while  ((ch = fgetc (fp)) != EOF && ch == ' ')
     ;
   if  (ch == EOF)
       return  false;
   ungetc (ch, fp);

   // put rest of line into  hdr
   while  ((ch = fgetc (fp)) != EOF && ch != '\n')
     hdr . push_back (char (ch));
   
   // trim up to first space
   //hdr = split(hdr)[0];

   // put everything up till next '>' into  s
   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     {
      if  (! isspace (ch))
          s . push_back (char (ch));
     }

   if  (ch == '>')
       ungetc (ch, fp);

   return  true;
  }



