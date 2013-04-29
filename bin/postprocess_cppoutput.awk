# This awk script is used in f90cpp and by all Makefiles generated by
# the FEAT2 configure script.
#
# It parses the output cpp generates and looks for lines starting with
# a hash character ('#'). These lines are usually interpreted by the
# C/C++ compiler to ensure that error message of the C preprocessed
# source code match those of the original code. Fortran compiler are
# not set up like this. If we debug the code or the runtime environments
# reports an error in line X, this awk code ensures that the error really
# occurred in line X of the code we are editing.
#
# Known drawbacks:
# Script only works if lines have been omitted, not if some have been
# added.
#
# Developer notice:
#    Evaluate the control commands from the
#    parsed file to the appropriate number of lines. Without this the
#    parsed file has a different number of lines which could cause that
#    error messages of the compiler point to the wrong line (source file
#    and compile file are different without translation).
#
#    Comments on the awk statement:
#    cpp will add lines of the format
#      # <line number within original source file where to continue> <name of file parsed>
#    e.g.
#      # 1 "../../sbblas/sourcemvline.f"
#      # 1 "<built-in>"
#      # 1 "<command line>"
#      # 1 "../../sbblas/sourcemvline.f"
#    Whenever including a file via #include there is an additional fourth
#    argument (which is '1' or '2') in the output of cpp. For our purpose
#    only lines with empty fourth argument or a fourth argument equal to
#    2 are of interest. There is one problem because we do not filter for
#    argument 3 being the current file name: recursive includes via
#    #include are not possible as they would ratten the line number
#    recovery tried here.

# BEGIN {
#     # initialise variable
#     corr = 1;
# }
# {
#     if ($1 == "#") {
# 	if ($4 == "" || $4 == "2") {
# 	    j = corr;
# 	    k = $2;
# 	    for (i = j + 1; i <= k; i++)
# 		print "! QQX";
# 	    corr = corr + k - j;
# 	}
#     } else {
# 	print;
# 	corr = corr + 1;
#     }
# }

BEGIN {
    # Initialise variables
    line = 1;
    filename = "";
}
{
    # The first line states the file name 
    # Store name of the file which is currently processed
    # and stringify it so that it can be compared with 
    # filenames generated by cpp.
#    if (filename != FILENAME) {
#	filename = "\"" FILENAME "\"";
#    }
    if (filename == "") {
	filename = $3;
    }

    # Read through file line by line and perform special
    # treatment of lines starting by #-symbol
    if ($1 == "#") {

	# Check of the fourth field is either empty or 2
	# which indicates returning from a file which has
	# been included in one of the previous lines.
	if ($4 == "" || $4 == "2") {

	    # To allow for recursive includes we need to
	    # check of we are back in the main file or
	    # just in one of the recursively included ones
	    if ($3 == filename) {
		# Store linenumber where to continue in
		# the main file generated by cpp
		l = $2;
		for (i = line+1; i <= l; i++)
		    print "! QQX ";
		line = l;
	    }
	}
	# Ignore all other lines starting ith #-symbol
    } else {
	# Just print current line
     	print;
     	line = line + 1;
    }
}
