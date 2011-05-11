//  A. L. Delcher
//
//  File:  xlate_tables.hh
//
//  Last Modified:  Tue May  9 10:25:40 EDT 2006
//
//  DNA to amino-acid translation tables



#ifndef  __XLATE_TABLES_HH_INCLUDED
#define  __XLATE_TABLES_HH_INCLUDED

const bool  IS_AMINO [26] = {
   true,   // A
   false,  // B
   true,   // C
   true,   // D
   true,   // E
   true,   // F
   true,   // G
   true,   // H
   true,   // I
   false,  // J
   true,   // K
   true,   // L
   true,   // M
   true,   // N
   false,  // O
   true,   // P
   true,   // Q
   true,   // R
   true,   // S
   true,   // T
   false,  // U
   true,   // V
   true,   // W
   false,  // X
   true,   // Y
   false   // Z
  };

// The Standard Code
const char  CODON_XLATE_TABLE_1 [] = 
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// The Vertebrate Mitochondrial Code
const char  CODON_XLATE_TABLE_2 [] =
  "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// The Yeast Mitochondrial Code
const char  CODON_XLATE_TABLE_3 [] = 
  "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// The Mold, Protozoan, and Coelenterate Mitochondrial Code
//   and the Mycoplasma/Spiroplasma Code
const char  CODON_XLATE_TABLE_4 [] = 
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// The Invertebrate Mitochondrial Code 
const char  CODON_XLATE_TABLE_5 [] = 
  "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// The Ciliate, Dasycladacean and Hexamita Nuclear Code
const char  CODON_XLATE_TABLE_6 [] = 
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// The Echinoderm and Flatworm Mitochondrial Code
const char  CODON_XLATE_TABLE_9 [] = 
  "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// The Euplotid Nuclear Code
const char  CODON_XLATE_TABLE_10 [] = 
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// The Alternative Yeast Nuclear Code
const char  CODON_XLATE_TABLE_12 [] = 
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// The Ascidian Mitochondrial Code
const char  CODON_XLATE_TABLE_13 [] = 
  "KNKNTTTTGGRSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// The Alternative Flatworm Mitochondrial Code
const char  CODON_XLATE_TABLE_14 [] = 
  "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// Blepharisma Nuclear Code
const char  CODON_XLATE_TABLE_15 [] = 
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// Chlorophycean Mitochondrial Code
const char  CODON_XLATE_TABLE_16 [] = 
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// Trematode Mitochondrial Code
const char  CODON_XLATE_TABLE_21 [] = 
  "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// Scenedesmus obliquus mitochondrial Code
const char  CODON_XLATE_TABLE_22[] = 
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVLY*Y*SSS*CWCLFLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

// Thraustochytrium Mitochondrial Code
const char  CODON_XLATE_TABLE_23 [] = 
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF";
// aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
// aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
// acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt

#endif
