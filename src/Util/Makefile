# Makefile for Glimmer3/src/Util

LOCAL_WORK = $(shell cd ../..; pwd)

UTIL_SRCS = entropy-profile.cc entropy-score.cc entropy-fasta.cc extract.cc multi-extract.cc \
  start-codon-distrib.cc uncovered.cc window-acgt.cc

SOURCES = $(UTIL_SRCS)
OBJECTS = $(UTIL_OBJS)

PROGS = entropy-profile entropy-score entropy-fasta extract multi-extract start-codon-distrib \
  uncovered window-acgt

LIBRARIES = 

include  $(LOCAL_WORK)/src/c_make.glm

all:    $(OBJECTS) $(LIBRARIES) $(PROGS)

entropy-profile:  entropy-profile.o libGLMcommon.a

entropy-score:  entropy-score.o libGLMcommon.a

entropy-fasta: entropy-fasta.o libGLMcommon.a

extract:  extract.o libGLMcommon.a

multi-extract:  multi-extract.o libGLMcommon.a

start-codon-distrib:  start-codon-distrib.o libGLMcommon.a

uncovered:  uncovered.o libGLMcommon.a

window-acgt:  window-acgt.o libGLMcommon.a

