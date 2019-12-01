CC = gcc

FFTWDIR = /usr/local
FFTWINC = -I$(FFTWDIR)/include
FFTWLIB = -L$(FFTWDIR)/lib -lfftw3

# define WALLCLOCK on the CFLAGS line for walltime instead of cpu time (default)
# -U__STRICT_ANSI__ - M_PI for GCC
CFLAGS = -O3 ${FFTWINC} -std=c11 -m64 -fPIC -U__STRICT_ANSI__

LDFLAGS = -lm -m64 -fPIC

# naive
NAIVESRC = primitive.c pmls.c naive_synthesis.c \
	makeweights.c csecond.c

NAIVEOBJ = primitive.c pmls.o naive_synthesis.o \
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

ALLSRC = FST_semi_fly.c FST_semi_memo.c cospmls.c csecond.c \
	makeweights.c naive_synthesis.c pmls.c primitive.c \
	seminaive.c test_conv_semi_fly.c test_conv_semi_memo.c \
	test_naive.c test_s2_semi_fly.c test_s2_semi_memo.c \
	test_s2_semi_memo_for.c test_s2_semi_memo_inv.c \
	test_semi.c


configure:
	mkdir -p bin

all:
	make \
	configure \
	legendre \
	sphere

legendre:
	make \
	configure \
	test_naive \
	test_semi

sphere:
	make \
	configure \
	test_s2_semi_memo \
	test_s2_semi_memo_for \
	test_s2_semi_memo_inv \
	test_s2_semi_fly \
	test_conv_semi_memo \
	test_conv_semi_fly

depend:
	makedepend ${FFTWINC} ${ALLSRC}

clean:
	rm -vf *.o
	rm -vrf ./bin


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

test_s2_semi_memo_for: $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_for.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_for.o \
	${FFTWLIB} $(LDFLAGS) -o bin/test_s2_semi_memo_for

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

test_s2_semi_memo_for.o: makeweights.h cospmls.h FST_semi_memo.h csecond.h

test_s2_semi_memo_inv.o: makeweights.h cospmls.h FST_semi_memo.h csecond.h

test_semi.o: makeweights.h cospmls.h primitive.h seminaive.h csecond.h
