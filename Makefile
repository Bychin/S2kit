CC = gcc

FFTW_DIR = /usr/local
FFTW_INC = -I$(FFTW_DIR)/include
FFTW_LIB = -L$(FFTW_DIR)/lib -lfftw3

SRC_DIR = src
SRC_INC = -I$(SRC_DIR)/

TEST_DIR = test
BUILD_DIR = build
BIN_DIR = bin

# define WALLCLOCK on the CFLAGS line for walltime instead of cpu time (default)
# -U__STRICT_ANSI__ - math constants (e.g. M_PI) for GCC
CFLAGS = -O3 $(FFTW_INC) $(SRC_INC) -std=c11 -m64 -fPIC -U__STRICT_ANSI__

LDFLAGS = -lm -m64 -fPIC

# naive Legendre transform
# TODO check src
DLT_NAIVE_SRC = chebyshev_nodes.c pmm.c l2_norms.c vector_funcs.c pml.c naive.c weights.c csecond.c
DLT_NAIVE_OBJ = $(addsuffix .o,$(basename $(DLT_NAIVE_SRC)))

# semi-naive Legendre transform
DLT_SEMI_SRC = chebyshev_nodes.c pmm.c pml.c cospml.c seminaive.c csecond.c l2_norms.c vector_funcs.c weights.c
DLT_SEMI_OBJ = $(addsuffix .o,$(basename $(DLT_SEMI_SRC)))

# semi-naive spherical transform and convolution
FST_SEMI_SRC = $(DLT_SEMI_SRC) naive.c util.c
FST_SEMI_OBJ = $(DLT_SEMI_OBJ) naive.o util.o
FST_SEMI_OBJ_FLY  = $(FST_SEMI_OBJ) FST_semi_fly.o
FST_SEMI_OBJ_MEMO = $(FST_SEMI_OBJ) FST_semi_memo.o

CONV_SEMI_SRC = $(FST_SEMI_SRC)
CONV_SEMI_OBJ = $(FST_SEMI_OBJ)
CONV_SEMI_OBJ_FLY  = $(CONV_SEMI_OBJ) FST_semi_fly.o
CONV_SEMI_OBJ_MEMO = $(CONV_SEMI_OBJ) FST_semi_memo.o


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
	test_DLT_naive \
	test_DLT_semi

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

test_DLT_naive: $(DLT_NAIVE_OBJ) test_DLT_naive.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(DLT_NAIVE_OBJ)) $(BUILD_DIR)/test_DLT_naive.o \
	$(LDFLAGS) -o $(BIN_DIR)/test_DLT_naive

test_DLT_semi: $(DLT_SEMI_OBJ) test_DLT_semi.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(DLT_SEMI_OBJ)) $(BUILD_DIR)/test_DLT_semi.o \
	$(FFTW_LIB) $(LDFLAGS) -o $(BIN_DIR)/test_DLT_semi

test_s2_semi_memo: $(FST_SEMI_OBJ) FST_semi_memo.o test_s2_semi_memo.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(FST_SEMI_OBJ_MEMO)) $(BUILD_DIR)/test_s2_semi_memo.o \
	$(FFTW_LIB) $(LDFLAGS) -o $(BIN_DIR)/test_s2_semi_memo

test_s2_semi_memo_fwd: $(FST_SEMI_OBJ) FST_semi_memo.o test_s2_semi_memo_fwd.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(FST_SEMI_OBJ_MEMO)) $(BUILD_DIR)/test_s2_semi_memo_fwd.o \
	$(FFTW_LIB) $(LDFLAGS) -o $(BIN_DIR)/test_s2_semi_memo_fwd

test_s2_semi_memo_inv: $(FST_SEMI_OBJ) FST_semi_memo.o test_s2_semi_memo_inv.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(FST_SEMI_OBJ_MEMO)) $(BUILD_DIR)/test_s2_semi_memo_inv.o \
	$(FFTW_LIB) $(LDFLAGS) -o $(BIN_DIR)/test_s2_semi_memo_inv

test_s2_semi_fly: $(FST_SEMI_OBJ) FST_semi_fly.o test_s2_semi_fly.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(FST_SEMI_OBJ_FLY)) $(BUILD_DIR)/test_s2_semi_fly.o \
	$(FFTW_LIB) $(LDFLAGS) -o $(BIN_DIR)/test_s2_semi_fly

test_conv_semi_memo: $(CONV_SEMI_OBJ) FST_semi_memo.o test_conv_semi_memo.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(CONV_SEMI_OBJ_MEMO)) $(BUILD_DIR)/test_conv_semi_memo.o \
	$(FFTW_LIB) $(LDFLAGS) -o $(BIN_DIR)/test_conv_semi_memo

test_conv_semi_fly: $(CONV_SEMI_OBJ) FST_semi_fly.o test_conv_semi_fly.o
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(CONV_SEMI_OBJ_FLY)) $(BUILD_DIR)/test_conv_semi_fly.o \
	$(FFTW_LIB) $(LDFLAGS) -o $(BIN_DIR)/test_conv_semi_fly


# explicit rule for every file to be compiled in BUILD_DIR

FST_semi_fly.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/FST_semi_fly.o $(SRC_DIR)/FST_semi_fly.c

FST_semi_memo.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/FST_semi_memo.o $(SRC_DIR)/FST_semi_memo.c

cospml.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/cospml.o $(SRC_DIR)/legendre_polynomials/cospml.c

pml.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/pml.o $(SRC_DIR)/legendre_polynomials/pml.c

pmm.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/pmm.o $(SRC_DIR)/legendre_polynomials/pmm.c

l2_norms.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/l2_norms.o $(SRC_DIR)/legendre_polynomials/util/l2_norms.c

vector_funcs.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/vector_funcs.o $(SRC_DIR)/legendre_polynomials/util/vector_funcs.c

naive.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/naive.o $(SRC_DIR)/legendre_transform/naive.c

seminaive.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/seminaive.o $(SRC_DIR)/legendre_transform/seminaive.c

weights.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/weights.o $(SRC_DIR)/legendre_transform/weights.c

chebyshev_nodes.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/chebyshev_nodes.o $(SRC_DIR)/util/chebyshev_nodes.c

util.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/util.o $(SRC_DIR)/util/util.c

csecond.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/csecond.o $(TEST_DIR)/util/csecond.c

test_conv_semi_fly.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_conv_semi_fly.o $(TEST_DIR)/test_conv_semi_fly.c

test_conv_semi_memo.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_conv_semi_memo.o $(TEST_DIR)/test_conv_semi_memo.c

test_DLT_naive.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_DLT_naive.o $(TEST_DIR)/test_DLT_naive.c

test_DLT_semi.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_DLT_semi.o $(TEST_DIR)/test_DLT_semi.c

test_s2_semi_fly.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_s2_semi_fly.o $(TEST_DIR)/test_s2_semi_fly.c

test_s2_semi_memo.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_s2_semi_memo.o $(TEST_DIR)/test_s2_semi_memo.c

test_s2_semi_memo_fwd.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_s2_semi_memo_fwd.o $(TEST_DIR)/test_s2_semi_memo_fwd.c

test_s2_semi_memo_inv.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/test_s2_semi_memo_inv.o $(TEST_DIR)/test_s2_semi_memo_inv.c
