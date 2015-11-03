#!/usr/bin/perl

use strict;

$| = 1;

system('wget -O .taxonomyData/.1_fromNCBI/bacteria_taxonomy.html "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=2&lvl=50&srchmode=1&keep=1&unlock" >& /dev/null');

system('wget -O .taxonomyData/.1_fromNCBI/archaea_taxonomy.html "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&name=Archaea&lvl=50&srchmode=1&keep=1&unlock" >& /dev/null');
