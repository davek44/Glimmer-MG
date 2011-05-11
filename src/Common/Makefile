# Makefile for Glimmer3/src/Common

LOCAL_WORK = $(shell cd ../..; pwd)

COMMON_SRCS = delcher.cc fasta.cc gene.cc kelley.cc
COMMON_OBJS = $(COMMON_SRCS:.cc=.o)

SOURCES = $(COMMON_SRCS)
OBJECTS = $(COMMON_OBJS)

PROGS = 

LIBRARIES = libGLMcommon.a

include  $(LOCAL_WORK)/src/c_make.glm

all:    $(OBJECTS) $(LIBRARIES) $(PROGS)

libGLMcommon.a:  $(COMMON_OBJS)
