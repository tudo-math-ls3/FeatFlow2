#! /bin/sh

collect_results() {
  nresultcount="$1"
  directorypattern="$2"
  sourcefileforces="$3"
  sourcefilepointvalues="$4"
  resultfileforces="$5"
  resultfilepointvalues="$6"
  
  # Collect the results and write to the files.
  head --lines=+2 "$directorypattern"1/$sourcefileforces > $resultfileforces
  head --lines=+2 "$directorypattern"1/$sourcefilepointvalues > $resultfilepointvalues
  for ((  i = 1 ;  i <= $nresultcount;  i++  ))
  do
    tail --lines=+3 "$directorypattern"$i/$sourcefileforces >> $resultfileforces
    tail --lines=+3 "$directorypattern"$i/$sourcefilepointvalues >> $resultfilepointvalues
  done
  
}


#####################################################################
# Script to extract the reference results from the calculated tests.
# To be called from inside of the benchmark directory.
# Creates a subdirectors "results_cc2d_bench2_reference" and places
# the reference results of all levels there.
#####################################################################

mkdir ./results_cc2d_bench2_reference

collect_results 4   ./logs/cc2d_bench2_reference_q1t_org_step4_dt1-lv3-DTMIN bdforces_lv3 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv3_dt1 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv3_dt1
collect_results 10  ./logs/cc2d_bench2_reference_q1t_org_step4_dt1-lv4-DTMIN bdforces_lv4 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv4_dt1 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv4_dt1
collect_results 40  ./logs/cc2d_bench2_reference_q1t_org_step4_dt1-lv5-DTMIN bdforces_lv5 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv5_dt1 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv5_dt1
collect_results 100 ./logs/cc2d_bench2_reference_q1t_org_step4_dt1-lv6-DTMIN bdforces_lv6 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv6_dt1 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv6_dt1

collect_results 2   ./logs/cc2d_bench2_reference_q2_step4_dt1-lv2-DTMIN bdforces_lv2 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv2_dt1 ./results_cc2d_bench2_reference/pointvalues_q2_lv2_dt1
collect_results 4   ./logs/cc2d_bench2_reference_q2_step4_dt1-lv3-DTMIN bdforces_lv3 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv3_dt1 ./results_cc2d_bench2_reference/pointvalues_q2_lv3_dt1
collect_results 10  ./logs/cc2d_bench2_reference_q2_step4_dt1-lv4-DTMIN bdforces_lv4 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv4_dt1 ./results_cc2d_bench2_reference/pointvalues_q2_lv4_dt1
collect_results 40  ./logs/cc2d_bench2_reference_q2_step4_dt1-lv5-DTMIN bdforces_lv5 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv5_dt1 ./results_cc2d_bench2_reference/pointvalues_q2_lv5_dt1
# collect_results 100 ./logs/cc2d_bench2_reference_q2_step4_dt1-lv6-DTMIN bdforces_lv6 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv6_dt1 ./results_cc2d_bench2_reference/pointvalues_q2_lv6_dt1


collect_results 4   ./logs/cc2d_bench2_reference_q1t_org_step4_dt2-lv3-DTMIN bdforces_lv3 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv3_dt2 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv3_dt2
collect_results 10  ./logs/cc2d_bench2_reference_q1t_org_step4_dt2-lv4-DTMIN bdforces_lv4 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv4_dt2 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv4_dt2
collect_results 40  ./logs/cc2d_bench2_reference_q1t_org_step4_dt2-lv5-DTMIN bdforces_lv5 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv5_dt2 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv5_dt2
collect_results 100 ./logs/cc2d_bench2_reference_q1t_org_step4_dt2-lv6-DTMIN bdforces_lv6 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv6_dt2 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv6_dt2

collect_results 2   ./logs/cc2d_bench2_reference_q2_step4_dt2-lv3-DTMIN bdforces_lv3 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv3_dt2 ./results_cc2d_bench2_reference/pointvalues_q2_lv3_dt2
collect_results 4   ./logs/cc2d_bench2_reference_q2_step4_dt2-lv3-DTMIN bdforces_lv3 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv3_dt2 ./results_cc2d_bench2_reference/pointvalues_q2_lv3_dt2
collect_results 10  ./logs/cc2d_bench2_reference_q2_step4_dt2-lv4-DTMIN bdforces_lv4 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv4_dt2 ./results_cc2d_bench2_reference/pointvalues_q2_lv4_dt2
collect_results 40  ./logs/cc2d_bench2_reference_q2_step4_dt2-lv5-DTMIN bdforces_lv5 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv5_dt2 ./results_cc2d_bench2_reference/pointvalues_q2_lv5_dt2
# collect_results 100 ./logs/cc2d_bench2_reference_q2_step4_dt2-lv6-DTMIN bdforces_lv6 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv6_dt2 ./results_cc2d_bench2_reference/pointvalues_q2_lv6_dt2


collect_results 4   ./logs/cc2d_bench2_reference_q1t_org_step4_dt3-lv3-DTMIN bdforces_lv3 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv3_dt3 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv3_dt3
collect_results 10  ./logs/cc2d_bench2_reference_q1t_org_step4_dt3-lv4-DTMIN bdforces_lv4 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv4_dt3 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv4_dt3
collect_results 40  ./logs/cc2d_bench2_reference_q1t_org_step4_dt3-lv5-DTMIN bdforces_lv5 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv5_dt3 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv5_dt3
collect_results 100 ./logs/cc2d_bench2_reference_q1t_org_step4_dt3-lv6-DTMIN bdforces_lv6 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv6_dt3 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv6_dt3

collect_results 2   ./logs/cc2d_bench2_reference_q2_step4_dt3-lv2-DTMIN bdforces_lv2 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv2_dt3 ./results_cc2d_bench2_reference/pointvalues_q2_lv2_dt3
collect_results 4   ./logs/cc2d_bench2_reference_q2_step4_dt3-lv3-DTMIN bdforces_lv3 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv3_dt3 ./results_cc2d_bench2_reference/pointvalues_q2_lv3_dt3
collect_results 10  ./logs/cc2d_bench2_reference_q2_step4_dt3-lv4-DTMIN bdforces_lv4 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv4_dt3 ./results_cc2d_bench2_reference/pointvalues_q2_lv4_dt3
collect_results 40  ./logs/cc2d_bench2_reference_q2_step4_dt3-lv5-DTMIN bdforces_lv5 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv5_dt3 ./results_cc2d_bench2_reference/pointvalues_q2_lv5_dt3


collect_results 4   ./logs/cc2d_bench2_reference_q1t_org_step4_dt4-lv3-DTMIN bdforces_lv3 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv3_dt4 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv3_dt4
collect_results 10  ./logs/cc2d_bench2_reference_q1t_org_step4_dt4-lv4-DTMIN bdforces_lv4 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv4_dt4 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv4_dt4
collect_results 40  ./logs/cc2d_bench2_reference_q1t_org_step4_dt4-lv5-DTMIN bdforces_lv5 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv5_dt4 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv5_dt4
collect_results 100 ./logs/cc2d_bench2_reference_q1t_org_step4_dt4-lv6-DTMIN bdforces_lv6 pointvalues ./results_cc2d_bench2_reference/bdforces_q1t_org_lv6_dt4 ./results_cc2d_bench2_reference/pointvalues_q1t_org_lv6_dt4

collect_results 2   ./logs/cc2d_bench2_reference_q2_step4_dt4-lv2-DTMIN bdforces_lv2 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv2_dt4 ./results_cc2d_bench2_reference/pointvalues_q2_lv2_dt4
collect_results 4   ./logs/cc2d_bench2_reference_q2_step4_dt4-lv3-DTMIN bdforces_lv3 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv3_dt4 ./results_cc2d_bench2_reference/pointvalues_q2_lv3_dt4
collect_results 10  ./logs/cc2d_bench2_reference_q2_step4_dt4-lv4-DTMIN bdforces_lv4 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv4_dt4 ./results_cc2d_bench2_reference/pointvalues_q2_lv4_dt4
collect_results 40  ./logs/cc2d_bench2_reference_q2_step4_dt4-lv5-DTMIN bdforces_lv5 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv5_dt4 ./results_cc2d_bench2_reference/pointvalues_q2_lv5_dt4
collect_results 400 ./logs/cc2d_bench2_reference_q2_step4_dt4-lv6-DTMIN bdforces_lv6 pointvalues ./results_cc2d_bench2_reference/bdforces_q2_lv6_dt4 ./results_cc2d_bench2_reference/pointvalues_q2_lv6_dt4
