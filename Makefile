CC = gcc

FFTWDIR = /usr/local
FFTWINC = -I$(FFTWDIR)/include
FFTWLIB = -L$(FFTWDIR)/lib -lfftw3

# define WALLCLOCK on the CFLAGS line for walltime instead of cpu time (default)
# -U__STRICT_ANSI__ - M_PI for GCC
# TODO -fPIC?
CFLAGS = -O3 ${FFTWINC} -std=c11 -m64 -U__STRICT_ANSI__

LDFLAGS = -lm -m64

# naive
NAIVESRC = primitive.c pmls.c naive_synthesis.c \
	makeweights.c csecond.c

NAIVEOBJ = primitive.o pmls.o naive_synthesis.o \
	makeweights.o csecond.o

# semi-naive
SEMISRC = pmls.c cospmls.c seminaive.c \
	csecond.c primitive.c makeweights.c

SEMIOBJ = pmls.o cospmls.o seminaive.o \
	csecond.o primitive.o makeweights.o

# semi-naive spherical transform and convolution
FSTSEMISRC = $(SEMISRC) naive_synthesis.c

FSTSEMIOBJ = $(SEMIOBJ) naive_synthesis.o

CONVSEMISRC = $(FSTSEMISRC)

CONVSEMIOBJ = $(FSTSEMIOBJ)

SRC_DIR = src
OBJ_DIR = bin
TEST_DIR = tests

LIB_NAME = s2kit

SRC_LIB_FILES = $(wildcard $(SRC_DIR)/*.c)
SRC_TEST_FILES = $(wildcard $(TEST_DIR)/*.c)
ALL_SRC_FILES = $(SRC_LIB_FILES) $(SRC_TEST_FILES)


configure:
	mkdir -p $(OBJ_DIR)

# TODO: not working
lib:
	make \
	clean \
	configure

	$(CC) -c $(CFLAGS) ${FFTWLIB} $(LDFLAGS) $(SRC_LIB_FILES)
	ar rcs $(OBJ_DIR)/$(LIB_NAME).a $(wildcard *.o)
	rm *.o

test:
	make \
	lib \
	test_naive \
	test_semi \
	test_conv_semi_fly \
	test_conv_semi_memo \
	test_s2_semi_fly \
	test_s2_semi_memo_fwd \
	test_s2_semi_memo_inv \
	test_s2_semi_memo 

all:
	make \
	configure \
	legendre \
	sphere

clean:
	rm -vf *.o
	rm -vrf ./bin

legendre:
	make \
	configure \
	test_naive \
	test_semi

sphere:
	make \
	configure \
	test_s2_semi_memo \
	test_s2_semi_memo_fwd \
	test_s2_semi_memo_inv \
	test_s2_semi_fly \
	test_conv_semi_memo \
	test_conv_semi_fly

depend:
	makedepend ${FFTWINC} ${ALL_SRC_FILES}


# make definitions for the individual executables

test_naive: $(NAIVEOBJ) test_naive.o
	$(CC) $(CFLAGS) $(NAIVEOBJ) test_naive.o \
	$(LDFLAGS) -o bin/test_naive

test_semi: $(SEMIOBJ) test_semi.o
	$(CC) $(CFLAGS) $(SEMIOBJ) test_semi.o \
	${FFTWLIB} $(LDFLAGS) -o bin/test_semi

test_s2_semi_memo: $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo.o \
	${FFTWLIB} $(LDFLAGS) -o bin/test_s2_semi_memo

test_s2_semi_memo_fwd: $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_fwd.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_fwd.o \
	${FFTWLIB} $(LDFLAGS) -o bin/test_s2_semi_memo_fwd

test_s2_semi_memo_inv: $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_inv.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_inv.o \
	${FFTWLIB} $(LDFLAGS) -o bin/test_s2_semi_memo_inv

test_s2_semi_fly: $(FSTSEMIOBJ) FST_semi_fly.o test_s2_semi_fly.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_fly.o test_s2_semi_fly.o \
	${FFTWLIB} $(LDFLAGS) -o bin/test_s2_semi_fly

test_conv_semi_memo: $(CONVSEMIOBJ) FST_semi_memo.o test_conv_semi_memo.o
	$(CC) $(CFLAGS) $(CONVSEMIOBJ) FST_semi_memo.o test_conv_semi_memo.o \
	${FFTWLIB} $(LDFLAGS) -o bin/test_conv_semi_memo

test_conv_semi_fly: $(CONVSEMIOBJ) FST_semi_fly.o test_conv_semi_fly.o
	$(CC) $(CFLAGS) $(CONVSEMIOBJ) FST_semi_fly.o test_conv_semi_fly.o \
	${FFTWLIB} $(LDFLAGS) -o bin/test_conv_semi_fly


# dependencies for make depend

FST_semi_fly.o: primitive.h makeweights.h pmls.h cospmls.h naive_synthesis.h
FST_semi_fly.o: seminaive.h FST_semi_fly.h

FST_semi_memo.o: makeweights.h FST_semi_memo.h
FST_semi_memo.o: cospmls.h primitive.h naive_synthesis.h seminaive.h

cospmls.o: primitive.h pmls.h

pmls.o: primitive.h

seminaive.o: cospmls.h

test_conv_semi_fly.o: FST_semi_fly.h

test_conv_semi_memo.o: FST_semi_memo.h cospmls.h

test_naive.o: pmls.h makeweights.h naive_synthesis.h csecond.h

test_s2_semi_fly.o: makeweights.h FST_semi_fly.h csecond.h

test_s2_semi_memo.o: makeweights.h cospmls.h FST_semi_memo.h csecond.h

test_s2_semi_memo_fwd.o: makeweights.h cospmls.h FST_semi_memo.h csecond.h

test_s2_semi_memo_inv.o: makeweights.h cospmls.h FST_semi_memo.h csecond.h

test_semi.o: makeweights.h cospmls.h primitive.h seminaive.h csecond.h
