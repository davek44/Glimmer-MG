//  A. L. Delcher
//
//  File:  anomaly.hh
//
//  Last Modified:  Fri May 19 17:10:05 EDT 2006
//
//  Declarations for  anomaly.cc



#ifndef  __ANOMALY_HH_INCLUDED
#define  __ANOMALY_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "gene.hh"


static bool  Is_Start_Codon
    (const char * p);
static bool  Is_Stop_Codon
    (const char * p);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Set_Start_And_Stop_Codons
    (void);
static void  Usage
    (void);


#endif
