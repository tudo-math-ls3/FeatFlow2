# -*- mode: sh -*-
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default settings for *all* tests defined in this directory
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description of the file format:
# This file is supposed to be a simple ASCII text file. Settings are
# defined per line via
#   keyword  = value
# The scope of a keyword/value setting is from its definition till
# being overwritten or till the end of file.
# A test definition starts with
#   testid       = ARBITRARY_TESTNAME
# All keyword/value settings defined earlier - as well as those defined in
# the file defaults.fbdef - are automatically assigned to this test definition.
# The value of any keyword can, of course, be subsequently overwritten
# via a line "keyword = value". This updated keyword/value setting will then
# be the default for any subsequently defined test.
# Anything behind a hash sign (#) is regarded as comment.
#
# Every keyword found is exported as environment variable (in uppercase).
# The case of the corresponding values is left untouched.
#
# Example for keyword/value settings inheritance: The following block
#
#	datfile = data/master_app1.dat
#	mglevels = 1,2,3,4
#
#	testid   = FOO1
#	appl = baz
#
#	testid   = FOO2
#	datfile = data/master_app2.dat
#
#	testid   = FOO3
#	appl = foobar
#	mglevels = 5,6
#
# is equivalent to
#
#	testid   = FOO1
#	appl = baz
#	datfile = data/master_app1.dat
#	mglevels = 1,2,3,4
#
#	testid   = FOO2
#	appl = baz
#	datfile = data/master_app2.dat
#	mglevels = 1,2,3,4
#
#	testid   = FOO3
#	appl = foobar
#	datfile = data/master_app2.dat
#	mglevels = 5,6
#


# FEAT2 only supports serial execution, don't try to run any application
# with MPI.
execmode = serial


# comma-separated list of names of hosts a benchmark test should
# be run on.
# To run a test on any host, specify 'all' or 'any' (without the
# single ticks). To disable execution in general, specify 'none'.
#
# Let's just hope that nobody will ever have access to a machine
# baptised any of those keywords. :-)
#
# Example:
#	allowhosts       = cassini,venera,viking
allowhosts       = all


# comma-separated list of names of hosts a benchmark test should
# *never* be run on.
# Special values: 'all', 'any', 'none'
# Complementary to the 'allowhosts' keyword in the following sense:
# A host added via 'allowhosts' is automatically removed from the
# whitelist when it occurs in a subsequent 'denyhosts' statement.
# Vice versa holds for the blacklist.
denyhosts        = dantooine,tatooine,jerusalem,jericho,peppone


# comma-separated lists of build IDs (standard FEAT2 five dash-separated tokens
# to identify architecture, cpu, operating system, compiler suite, BLAS
# implementation) a benchmark test should be (never) be run on.
# The mechanism works exactly like 'allowhosts' and supports the same special
# values 'all', 'any' and 'none'. The remark on the complementary nature of
# 'allowhosts'/'denyhosts' holds for 'allowbuildID'/'denybuildID', too.
#
# Example:
#	allowbuildIDs = all
#	denybuildIDs  = pc-athlon64-linux-gcc-goto,pc64-opteron-linux-intel-mkl
# will run everywhere but where the two forbidden build IDs match.
allowbuildIDs = all
denybuildIDs  = none
