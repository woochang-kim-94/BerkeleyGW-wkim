PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

SOBJS = sutil.o psutil.o \
	sssmv.o sssr2.o sssmm.o sssr2k.o \
	ssstd2.o ssstrd.o slastd.o \
	psssmv.o psssr2k.o \
	pssstd2.o pssstrd.o pslastd.o
DOBJS = dutil.o pdutil.o \
	dssmv.o dssr2.o dssmm.o dssr2k.o \
	dsstd2.o dsstrd.o dlastd.o \
	pdssmv.o pdssr2k.o \
	pdsstd2.o pdsstrd.o pdlastd.o
OBJS = $(SOBJS) $(DOBJS)
LIB = libsseig.a

.PHONY: sreal dreal lib

default: lib
all: lib sreal dreal

lib: $(LIB)
$(LIB): $(OBJS)
	$(AR) ru $(LIB) $(OBJS)

sreal: stests
stests_srcs = $(wildcard stest*.f90 pstest*.f90)
stests: $(stests_srcs:.f90=.x)
stest%.x: stest%.o $(LIB)
	$(LINK) $(FOPTS) -o $@ $^
pstest%.x: pstest%.o $(LIB)
	$(LINK) $(FOPTS) -o $@ $^

dreal: dtests
dtests_srcs = $(wildcard dtest*.f90 pdtest*.f90)
dtests: $(dtests_srcs:.f90=.x)
dtest%.x: dtest%.o $(LIB)
	$(LINK) $(FOPTS) -o $@ $^
pdtest%.x: pdtest%.o $(LIB)
	$(LINK) $(FOPTS) -o $@ $^
