# Makefile for Glimmer3/src/ICM

LOCAL_WORK = $(shell cd ../..; pwd)

ICM_SRCS = icm.cc build-icm.cc build-fixed.cc score-fixed.cc
ICM_OBJS = $(ICM_SRCS:.cc=.o)

SOURCES = $(ICM_SRCS)
OBJECTS = $(ICM_OBJS)

PROGS = build-icm build-fixed score-fixed

LIBRARIES = libGLMicm.a

include  $(LOCAL_WORK)/src/c_make.glm

all:    $(OBJECTS) $(LIBRARIES) $(PROGS)

build-icm:  build-icm.o libGLMicm.a libGLMcommon.a

build-fixed:  build-fixed.o libGLMicm.a libGLMcommon.a

score-fixed:  score-fixed.o libGLMicm.a libGLMcommon.a

libGLMicm.a:  $(ICM_OBJS)
