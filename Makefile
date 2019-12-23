CC = gcc

FFTWDIR = /usr/local
FFTWINC = -I$(FFTWDIR)/include
FFTWLIB = -L$(FFTWDIR)/lib -lfftw3

SRC_DIR = src
SRC_INC = -I$(SRC_DIR)/
TEST_DIR = test
BUILD_DIR = build
BIN_DIR = bin

# define WALLCLOCK on the CFLAGS line for walltime instead of cpu time (default)
# -U__STRICT_ANSI__ - math constants (e.g. M_PI) for GCC
CFLAGS = -O3 ${FFTWINC} ${SRC_INC} -std=c11 -m64 -fPIC -U__STRICT_ANSI__

LDFLAGS = -lm -m64 -fPIC

# naive
NAIVESRC = primitive.c pmls.c naive_synthesis.c makeweights.c csecond.c
NAIVEOBJ = $(addsuffix .o,$(basename $(NAIVESRC)))

# semi-naive
SEMISRC = pmls.c cospmls.c seminaive.c csecond.c primitive.c makeweights.c
SEMIOBJ = $(addsuffix .o,$(basename $(SEMISRC)))

# semi-naive spherical transform and convolution
FSTSEMISRC = $(SEMISRC) naive_synthesis.c
FSTSEMIOBJ = $(SEMIOBJ) naive_synthesis.o
FSTSEMIOBJFLY  = $(FSTSEMIOBJ) FST_semi_fly.o
FSTSEMIOBJMEMO = $(FSTSEMIOBJ) FST_semi_memo.o

CONVSEMISRC = $(FSTSEMISRC)
CONVSEMIOBJ = $(FSTSEMIOBJ)
CONVSEMIOBJFLY  = $(CONVSEMIOBJ) FST_semi_fly.o
CONVSEMIOBJMEMO = $(CONVSEMIOBJ) FST_semi_memo.o


configure:
	mkdir -p $(BIN_DIR)
	mkdir -p $(BUILD_DIR)

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
	test_s2_semi_memo_fwd \
	test_s2_semi_memo_inv \
	test_s2_semi_fly \
	test_conv_semi_memo \
	test_conv_semi_fly

clean:
	rm -rf $(BIN_DIR)
	rm -rf $(BUILD_DIR)


# make definitions for the individual executables

test_naive: $(NAIVEOBJ) test_naive.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(NAIVEOBJ)) $(BUILD_DIR)/test_naive.o \
	$(LDFLAGS) -o $(BIN_DIR)/test_naive

test_semi: $(SEMIOBJ) test_semi.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(SEMIOBJ)) $(BUILD_DIR)/test_semi.o \
	${FFTWLIB} $(LDFLAGS) -o $(BIN_DIR)/test_semi

test_s2_semi_memo: $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(FSTSEMIOBJMEMO)) $(BUILD_DIR)/test_s2_semi_memo.o \
	${FFTWLIB} $(LDFLAGS) -o $(BIN_DIR)/test_s2_semi_memo

test_s2_semi_memo_fwd: $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_fwd.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(FSTSEMIOBJMEMO)) $(BUILD_DIR)/test_s2_semi_memo_fwd.o \
	${FFTWLIB} $(LDFLAGS) -o $(BIN_DIR)/test_s2_semi_memo_fwd

test_s2_semi_memo_inv: $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_inv.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(FSTSEMIOBJMEMO)) $(BUILD_DIR)/test_s2_semi_memo_inv.o \
	${FFTWLIB} $(LDFLAGS) -o $(BIN_DIR)/test_s2_semi_memo_inv

test_s2_semi_fly: $(FSTSEMIOBJ) FST_semi_fly.o test_s2_semi_fly.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(FSTSEMIOBJFLY)) $(BUILD_DIR)/test_s2_semi_fly.o \
	${FFTWLIB} $(LDFLAGS) -o $(BIN_DIR)/test_s2_semi_fly

test_conv_semi_memo: $(CONVSEMIOBJ) FST_semi_memo.o test_conv_semi_memo.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(CONVSEMIOBJMEMO)) $(BUILD_DIR)/test_conv_semi_memo.o \
	${FFTWLIB} $(LDFLAGS) -o $(BIN_DIR)/test_conv_semi_memo

test_conv_semi_fly: $(CONVSEMIOBJ) FST_semi_fly.o test_conv_semi_fly.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(CONVSEMIOBJFLY)) $(BUILD_DIR)/test_conv_semi_fly.o \
	${FFTWLIB} $(LDFLAGS) -o $(BIN_DIR)/test_conv_semi_fly


# explicit rule for every file to be compiled in BUILD_DIR

FST_semi_fly.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/FST_semi_fly.o $(SRC_DIR)/FST_semi_fly.c

FST_semi_memo.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/FST_semi_memo.o $(SRC_DIR)/FST_semi_memo.c

cospmls.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/cospmls.o $(SRC_DIR)/cospmls.c

makeweights.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/makeweights.o $(SRC_DIR)/makeweights.c

naive_synthesis.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/naive_synthesis.o $(SRC_DIR)/naive_synthesis.c

pmls.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/pmls.o $(SRC_DIR)/pmls.c

primitive.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/primitive.o $(SRC_DIR)/util/primitive.c

seminaive.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/seminaive.o $(SRC_DIR)/seminaive.c

csecond.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/csecond.o $(TEST_DIR)/util/csecond.c

test_conv_semi_fly.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_conv_semi_fly.o $(TEST_DIR)/test_conv_semi_fly.c

test_conv_semi_memo.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_conv_semi_memo.o $(TEST_DIR)/test_conv_semi_memo.c

test_naive.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_naive.o $(TEST_DIR)/test_naive.c

test_s2_semi_fly.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_s2_semi_fly.o $(TEST_DIR)/test_s2_semi_fly.c

test_s2_semi_memo.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_s2_semi_memo.o $(TEST_DIR)/test_s2_semi_memo.c

test_s2_semi_memo_fwd.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_s2_semi_memo_fwd.o $(TEST_DIR)/test_s2_semi_memo_fwd.c

test_s2_semi_memo_inv.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_s2_semi_memo_inv.o $(TEST_DIR)/test_s2_semi_memo_inv.c

test_semi.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_semi.o $(TEST_DIR)/test_semi.c
