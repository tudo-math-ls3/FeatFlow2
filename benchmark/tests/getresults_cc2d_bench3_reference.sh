#! /bin/sh

#####################################################################
# Script to extract the reference results from the calculated tests.
# To be called from inside of the benchmark directory.
# Creates a subdirectors "results_cc2d_bench3_reference" and places
# the reference results of all levels there.
#####################################################################

mkdir ./results_cc2d_bench3_reference

# Level 1+2+3 can directly be obtained.

cp ./logs/cc2d_bench3_reference_q2-lv1-DTMIN1/bdforces_lv1 ./results_cc2d_bench3_reference/bdforces_lv1
cp ./logs/cc2d_bench3_reference_q2-lv2-DTMIN1/bdforces_lv2 ./results_cc2d_bench3_reference/bdforces_lv2
cp ./logs/cc2d_bench3_reference_q2-lv3-DTMIN1/bdforces_lv3 ./results_cc2d_bench3_reference/bdforces_lv3

cp ./logs/cc2d_bench3_reference_q2-lv1-DTMIN1/pointvalues ./results_cc2d_bench3_reference/pointvalues_lv1
cp ./logs/cc2d_bench3_reference_q2-lv2-DTMIN1/pointvalues ./results_cc2d_bench3_reference/pointvalues_lv2
cp ./logs/cc2d_bench3_reference_q2-lv3-DTMIN1/pointvalues ./results_cc2d_bench3_reference/pointvalues_lv3

# Build level 4 from the subtests
head --lines=+2 ./logs/cc2d_bench3_reference_q2-lv4-DTMIN1/bdforces_lv4 > ./results_cc2d_bench3_reference/bdforces_lv4
head --lines=+2 ./logs/cc2d_bench3_reference_q2-lv4-DTMIN1/pointvalues > ./results_cc2d_bench3_reference/pointvalues_lv4
for ((  i = 1 ;  i <= 8;  i++  ))
do
  tail --lines=+3 ./logs/cc2d_bench3_reference_q2-lv4-DTMIN$i/bdforces_lv4 >> ./results_cc2d_bench3_reference/bdforces_lv4
  tail --lines=+3 ./logs/cc2d_bench3_reference_q2-lv4-DTMIN$i/pointvalues >> ./results_cc2d_bench3_reference/pointvalues_lv4
done

# Build level 5 from the subtests
head --lines=+2 ./logs/cc2d_bench3_reference_q2-lv5-DTMIN$i/bdforces_lv5 > ./results_cc2d_bench3_reference/bdforces_lv5
head --lines=+2 ./logs/cc2d_bench3_reference_q2-lv5-DTMIN$i/pointvalues > ./results_cc2d_bench3_reference/pointvalues_lv5
for ((  i = 1 ;  i <= 40;  i++  ))
do
  tail --lines=+3 ./logs/cc2d_bench3_reference_q2-lv5-DTMIN$i/bdforces_lv5 >> ./results_cc2d_bench3_reference/bdforces_lv5
  tail --lines=+3 ./logs/cc2d_bench3_reference_q2-lv5-DTMIN$i/pointvalues >> ./results_cc2d_bench3_reference/pointvalues_lv5
done

# Build level 6 from the subtests
head --lines=+2 ./logs/cc2d_bench3_reference_q2-lv6-DTMIN$i/bdforces_lv6 > ./results_cc2d_bench3_reference/bdforces_lv6
head --lines=+2 ./logs/cc2d_bench3_reference_q2-lv6-DTMIN$i/pointvalues > ./results_cc2d_bench3_reference/pointvalues_lv6
for ((  i = 1 ;  i <= 160;  i++  ))
do
  tail --lines=+3 ./logs/cc2d_bench3_reference_q2-lv6-DTMIN$i/bdforces_lv6 >> ./results_cc2d_bench3_reference/bdforces_lv6
  tail --lines=+3 ./logs/cc2d_bench3_reference_q2-lv6-DTMIN$i/pointvalues >> ./results_cc2d_bench3_reference/pointvalues_lv6
done
