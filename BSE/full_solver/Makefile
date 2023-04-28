PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

SSEIG_ROOT = $(PREFIX)/BSE/SSEIG
SSEIG_LIB = $(SSEIG_ROOT)/libsseig.a
$(SSEIG_LIB):
	cd $(SSEIG_ROOT) && $(MAKE) lib

OBJS = pdbseig.o pdembed1.o pdembed2.o pdbssolver1.o pdbssolver1_svd.o \
	pzbseig.o pzembed1.o pzbssolver1.o pzbsaux1.o pdsorteig.o \
	pzembed2.o \
	zbssolver1.o zbsaux1.o dsorteig.o \
	pdbsabsp.o pdbslanczos.o \
	pzbsabsp.o pzbslanczos.o pzlaconj.o dbsabsp_aux.o approx_delta.o \
	pdlanczos.o pzlanczos.o \
	pddot_all.o pzdotc_all.o \
	pdstein.o pzstein.o \
	pdlaprnt.o pzlaprnt.o

FULL_SOLVER_OBJS := pdbseig.o pdembed1.o pdembed2.o pdbssolver1.o pdbssolver1_svd.o \
	pzbseig.o pzembed1.o pzbssolver1.o pzbsaux1.o pdsorteig.o pzembed2.o \
	zbssolver1.o zbsaux1.o dsorteig.o pdbsabsp.o pdbslanczos.o \
	pzbsabsp.o pzbslanczos.o pzlaconj.o dbsabsp_aux.o approx_delta.o \
	pdlanczos.o pzlanczos.o pddot_all.o pzdotc_all.o \
	pdstein.o pzstein.o pdlaprnt.o pzlaprnt.o \
	zbsabsp2.o pzsplanczos.o zembed3.o pzdotu_all.o

FULL_SOLVER_LIB = libfull_solver.a
lib: $(FULL_SOLVER_LIB)
$(FULL_SOLVER_LIB): $(FULL_SOLVER_OBJS)
	$(AR) ru $(FULL_SOLVER_LIB) $(FULL_SOLVER_OBJS)

default: lib
all: lib dmain.x dabsp.x zmain.x zabsp.x

dmain.x: dmain.o $(OBJS) $(SSEIG_LIB)
	$(LINK) $(FOPTS) -o $@ $^

dabsp.x: dabsp.o $(OBJS) $(SSEIG_LIB)
	$(LINK) $(FOPTS) -o $@ $^

zabsp.x: zabsp.o $(OBJS) $(SSEIG_LIB)
	$(LINK) $(FOPTS) -o $@ $^

zmain.x: zmain.o $(OBJS) $(SSEIG_LIB)
	$(LINK) $(FOPTS) -o $@ $^
