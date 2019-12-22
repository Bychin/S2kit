#!/bin/bash

BIN="./bin"
DATA="./data"

echo "--- MAKE ---"
make clean
make all
echo "--- END OF MAKE ---"

echo "--- EXECUTE TESTS ---"

echo "--- EXECUTE test_naive 0 1024 100 ---"
$BIN/test_naive 0 1024 100

echo "--- EXECUTE test_semi 0 1024 100 ---"
$BIN/test_semi 0 1024 100

echo "--- EXECUTE test_s2_semi_fly 512 1 ---"
$BIN/test_s2_semi_fly 512 1

echo "--- EXECUTE test_s2_semi_memo 512 1 ---"
$BIN/test_s2_semi_memo 512 1

echo "--- EXECUTE test_conv_semi_fly 64 ---"
$BIN/test_conv_semi_fly $DATA/s64.dat $DATA/f64.dat $DATA/o64_conv_semi_fly.dat 64
echo "--- DIFF for test_conv_semi_fly 64 ---"
diff $DATA/o64_conv_semi_fly_original.dat $DATA/o64_conv_semi_fly.dat
echo "--- END OF DIFF for test_conv_semi_fly 64 ---"

echo "--- EXECUTE test_conv_semi_fly 128 ---"
$BIN/test_conv_semi_fly $DATA/s128.dat $DATA/f128.dat $DATA/o128_conv_semi_fly.dat 128
echo "--- DIFF for test_conv_semi_fly 128 ---"
diff $DATA/o128_conv_semi_fly_original.dat $DATA/o128_conv_semi_fly.dat
echo "--- END OF DIFF for test_conv_semi_fly 128 ---"

echo "--- EXECUTE test_conv_semi_memo 64 ---"
$BIN/test_conv_semi_memo $DATA/s64.dat $DATA/f64.dat $DATA/o64_conv_semi_memo.dat 64
echo "--- DIFF for test_conv_semi_memo 64 ---"
diff $DATA/o64_conv_semi_memo_original.dat $DATA/o64_conv_semi_memo.dat
echo "--- END OF DIFF for test_conv_semi_memo 64 ---"

echo "--- EXECUTE test_conv_semi_memo 128 ---"
$BIN/test_conv_semi_memo $DATA/s128.dat $DATA/f128.dat $DATA/o128_conv_semi_memo.dat 128
echo "--- DIFF for test_conv_semi_memo 128 ---"
diff $DATA/o128_conv_semi_memo_original.dat $DATA/o128_conv_semi_memo.dat
echo "--- END OF DIFF for test_conv_semi_memo 128 ---"

echo "--- END OF TESTS ---"

make clean
