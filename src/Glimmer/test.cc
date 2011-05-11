//  A. Delcher
//
//  Test file for syntatic/system things


#include  "delcher.hh"
// #include  "fasta.hh"
// #include  "gene.hh"

// Test comment

int  main
    (int argc, char * argv [])

  {
   string  s, hdr;
   time_t  now;
   int  a, b;

   now = time (NULL);
   cerr << "Starting at " << ctime (& now) << endl;

   while  (scanf ("%d %d", & a, & b) == 2)
     {
      printf ("%3d  %3d  %9d\n", a, b, Int_Power (a, b));
     }

   return  0;
  }

