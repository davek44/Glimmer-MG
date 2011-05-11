//  A. L. Delcher
//
//  File:  fasta.hh
//
//  Last Modified:  23 October 2003
//
//  Routines to manipulate FASTA format files


#ifndef  __FASTA_H_INCLUDED
#define  __FASTA_H_INCLUDED


#include  "delcher.hh"
#include  <string>
#include  <vector>


const int  DEFAULT_FASTA_WIDTH = 60;
  // Max number of characters to print on a FASTA data line
const char  QUALITY_OFFSET = '0';
  // Value added to qualities to create a printable character


void  Fasta_Print
    (FILE * fp, const char * s, const char * hdr = NULL,
     int fasta_width = DEFAULT_FASTA_WIDTH);
void  Fasta_Print_N
    (FILE * fp, const char * s, int n, const char * hdr = NULL,
     int fasta_width = DEFAULT_FASTA_WIDTH);
void  Fasta_Print_Skip
    (FILE * fp, const char * s, const char * skip, const char * hdr = NULL,
     int fasta_width = DEFAULT_FASTA_WIDTH);
bool  Fasta_Qual_Vec_Read
    (FILE * fp, vector<int> & q, string & hdr);
bool  Fasta_Qual_Read
    (FILE * fp, string & q, string & hdr);
bool  Fasta_Read
    (FILE * fp, string & s, string & hdr);


#endif
