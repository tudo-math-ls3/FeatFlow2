#!/bin/bash

logfile=../log/flagship_20110329_132433.474.log

# Generate data file for inner convergence behavior
echo "Generating data file for inner convergence behavior"

sed -e '/^Number of nonlinear iterations/ { N; N; N; N; N; N; N; N; N; N; N; N; s/Improvement of residual[^\n]*\n//; s/Convergence rate[^\n]*\n//; s/####[^\n]*\n//g; s/Nonlinear solution[^\n]*\n//; s/Number of nonlinear iterations://g; s/Norm of final residual://g; s/Norm of initial residual://g; s/\n/\t/g; p; }; d;' /data/peppone/mmoelle1/Featflow2/application/flagship/log/flagship_20110329_132433.474.log > inner.dat

# Generate data file for outer convergence behavior
echo "Generating data file for outer convergence behavior"

sed -e '/^Two-level theta-scheme/ { s/^.*Time = \([0-9.E\-]*\) .*$/\1/; p; }; /^Coupled solution/ { N; N; N; s/Coupled solution[^\n]*\n//; s/Number of outer iterations://; s/Norm of final residual://; s/Norm of initial residual://; s/\n/\t/g; p; }; d;' /data/peppone/mmoelle1/Featflow2/application/flagship/log/flagship_20110329_132433.474.log | sed 'N; s/\n/ /;' > outer.dat
