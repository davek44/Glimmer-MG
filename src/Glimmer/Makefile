# Makefile for Glimmer3/src/Glimmer

LOCAL_WORK = $(shell cd ../..; pwd)

GLIMMER_SRCS = anomaly.cc glimmer_base.cc glimmer3.cc glimmer-mg.cc long-orfs.cc test.cc

SOURCES = $(GLIMMER_SRCS)
OBJECTS = $(GLIMMER_OBJS)

PROGS = anomaly glimmer3 glimmer-mg long-orfs test

LIBRARIES = 

include  $(LOCAL_WORK)/src/c_make.glm

all:    $(OBJECTS) $(LIBRARIES) $(PROGS)

anomaly:  anomaly.o libGLMcommon.a

glimmer3:  glimmer3.o glimmer_base.o libGLMcommon.a libGLMicm.a

glimmer-mg: glimmer-mg.o glimmer_base.o libGLMcommon.a libGLMicm.a

long-orfs:  long-orfs.o libGLMcommon.a

test:  test.o libGLMcommon.a
