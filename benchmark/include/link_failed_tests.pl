#!/usr/bin/env perl
#
# Script which turns output of a FEAST benchmark run into HTML format
# adding hyperreferences to and among every failed test (i.e. crashed, 
# deviant or without reference solution).
#
# It is not mandatory that the tests were run (or the output is) in 
# alphabetical order. When traversing the links, however, the "next"
# and "previous" links will respect alphabetical ordering.

use strict;
use warnings;
use Env qw(HOME LOGNAME PATH SHELL);

# Provide some constants defined in C fcntl.h
# (in particular SEEK_SET)
use Fcntl 'SEEK_SET';

my $debugScript = 0;

# Expect input:
# arg1 = from-address 
# arg2 = subject-prefix
# arg3 = subject-postfix
# arg4 = recipient
# arg5 = benchmark log file
# arg6 = list of conspicuous tests
my ($from_address, $prefix, $postfix, $recipient, $logfile) = @ARGV;

die "$0: No from statement given.\n"   if ($from_address =~ m/^\s*$/);
die "$0: No subject prefix given.\n"   if ($prefix =~ m/^\s*$/);
die "$0: No subject postfix given.\n"   if ($postfix =~ m/^\s*$/);
die "$0: No recipient given.\n" if ($recipient =~ m/^\s*$/);
die "$0: No logfile given.\n"   if ($logfile =~ m/^\s*$/);

# Extract conspicuous tests from log file
open (INFILE, "<", $logfile) or die "$0: Cannot open file <$logfile>: $!\n";
my @allconspicoustests = ();
my ($testsTotal, $testsSuccessfull, $testsConspicuous, $testsCrashed) = (0, 0, 0, 0);
my $skip = 1;
while (<INFILE>) {
    $skip = 0 if (m/^SUMMARY:/);
    next if ($skip);

    if (m/^tests in total\s*:\s*(\d+)/) {
	$testsTotal = $1;
    } elsif (m/^tests successfully passed\s*:\s*(\d+)/) {
	$testsSuccessfull = $1;
    } elsif (m/^tests (deviant|unverifiable)\s*:\s*(\d+)/) {
	$testsConspicuous += $2;
    } elsif (m/^tests failing execution\s*:\s*(\d+)/) {
	$testsCrashed += $1;
	last;
    }
}

if ($debugScript) {
    print "Parsing <$logfile> had the following result:\n" .
	  "* Total number of tests               : $testsTotal\n" .
	  "* Number of successfull tests         : $testsSuccessfull\n" .
	  "* Number of deviant/unverifiable tests: $testsConspicuous\n" .
	  "* Number of failed tests              : $testsCrashed\n\n";
}

# If any test failed, try to extract the test IDs that
# follow now in the log file
if ($testsConspicuous > 0 || $testsCrashed > 0) {
    FILE: while (<INFILE>) {
	if (m/^The tests with deviating results have the following IDs:/) {
	    # Read till we find first empty line indicating
	    # end of tests that deviated
	    $_ = <INFILE>;
	    last FILE if (eof INFILE);
	    chomp($_);
	    while (! m/^\s*$/) {
		push @allconspicoustests, $_;
		$_ = <INFILE>;
	    	last FILE if (eof INFILE);
		chomp($_);
	    }
	}

	if (m/^The tests that (could not be verified due to a missing|crashed during execution have the following IDs)/) {
	    # Skip next line as it's still text.
	    $_ = <INFILE>;

	    # Read till we find first empty line indicating
	    # end of tests that crashed / did not have reference solutions
	    $_ = <INFILE>;
	    last FILE if (eof INFILE);
	    chomp($_);
	    while (! m/^\s*$/) {
		push @allconspicoustests, $_;
		$_ = <INFILE>;
	    	last FILE if (eof INFILE);
		chomp($_);
	    }
	}
    }
}

# rewind file
seek(INFILE, 0, SEEK_SET);


# Store all conspicuous tests in a hash to be later on able to
# simply check if a hash entry with key "test id" is defined. If so,
# we will know that the test currently treated is on the list of
# conspicuous ones.
# (Very nifty: The corresponding value of each key is the next
#  conspicuous test ID. To be used for creating "next" links.)
my %nextFailedTest;
my $nextConspTest = "";
foreach my $entry (reverse sort @allconspicoustests) {
   $nextFailedTest{$entry} = $nextConspTest;
   $nextConspTest = $entry;
}
undef $nextConspTest;
# Same thing for links to previous failed test
my %prevFailedTest;
my $prevConspTest = "";
my $firstFailedTest = "";
foreach my $entry (sort @allconspicoustests) {
   $prevFailedTest{$entry} = $prevConspTest;
   $prevConspTest = $entry;

   # Determine first failed test
   $firstFailedTest = $entry if ($firstFailedTest eq "");
}
undef $prevConspTest;


my $message =
qq{From: $from_address
To: $recipient
Subject: $prefix $testsTotal|$testsConspicuous|$testsCrashed [tot|dev|crash] $postfix
Content-type: text/html; charset="us-ascii"

<html>
<head>
    <style type="text/css">
	div.link, div.anchor {
		margin: 0;
		font-size: 11pt;
		font-family: "Courier New", Courier, monospace, sans-serif;
	}

	pre {
		margin: 0;
	}
    </style>
</head>

<body>
<pre>
};

my $isConspicuous = 0;
my $testID = "";
my $firstTest = 1;
LINE: while (<INFILE>) {
    $_ =~ s/</&lt;/g;
    $_ =~ s/>/&gt;/g;

    # Start of output of first test?
    # Add link to first failed test (if any)
    # Add it directly before first
    #   # =====================
    #   # Test description:
    if (1 == $firstTest && $_ =~ m/^# [=]+$/) {
	$firstTest = 0;

	# Add link now
	$message .= qq{</pre>\n<div class="anchor">};
	if ($firstFailedTest ne "") {
	    $message .= qq{<a href="#$firstFailedTest">(first failed test)</a>\n};
	}
	# Always add link to summary at end of the e-mail.
	$message .= qq{<a href="#summary">(summary)</a>\n};
	$message .= qq{</div>\n<pre>\n\n};
	$message .= $_;
	next LINE;
    }

    # Start of output of a new test? Is test ID one of the conspicuous ones?
    if ($_ =~ m/^TESTID\s+=\s+(\S+)\s*$/) {
	$testID = $1;
        if (exists($nextFailedTest{$testID})) {
	    $message .= qq{</pre><a name="$testID"></a><pre>\n};
	}
    }

    # Insert link to next failed test
    if ($_ =~ m/^TEST (DEVIANT|FAILED|WITHOUT REFERENCE SOLUTION)/  &&  defined($nextFailedTest{$testID})) {
	# Case 1: Test failed
	#
	# typical fragment:
	#----
	#forrtl: error (78): process killed (SIGTERM)
	#TEST FAILED
	#
        #No information on FEAST run available due to execution failure.
	#
	#
	#FINISH  : 01:50:26
	#----
	#
	# => Add link immediately after "TEST FAILED"
	$message .= $_;
	my $link = "";
	if ($prevFailedTest{$testID} ne "") {
	    $link .= qq{<a href="#$prevFailedTest{$testID}">(previous failed test)</a>\n};
	} else {
	    $link .= qq{<a href="#">(beginning of e-mail)</a>\n};
	}
	if ($nextFailedTest{$testID} ne "") {
	    $link .= qq{<a href="#$nextFailedTest{$testID}">(next failed test)</a>\n};
	}
	# Always add link to summary at end of the e-mail.
	$link .= qq{<a href="#summary">(summary)</a>\n};
	$message .= qq{</pre>\n<div class="anchor">$link</div>\n<pre>\n};

	# Case 2: Test with deviating results
	#
	# typical fragment:
	#TEST DEVIANT
	#
	#LEVEL #NEQ   #NIT   abs. res     hmin        AR         L2err     true err    I_eff     L_oo err
	#1c1
	#< 5      1024   12    7.830E-13  2.085E-04  1.605E+02   1.316E-04  6.032E-03  1.046E+00  4.815E-04
	#---
	#> 5      1024   12    7.829E-13  2.085E-04  1.605E+02   1.316E-04  6.032E-03  1.046E+00  4.815E-04
	#
	#Difference:
	#LEVEL  #NEQ   #NIT    abs. res        hmin        AR       L2 err     Grad err        I_eff     L_oo err
	#    -    -     -    1.28e-04           -         -            -            -            -            -
	#    -    -     -           -           -         -            -            -            -            -
	#    -    -     -           -           -         -            -            -            -            -
	#Note: Values for abs. res, hmin, AR, L2 err, Grad err, I_eff and L_oo err are relative, not absolute!
	#
	#
	#FINISH  : 06:34:55
	#
	# => Add additional link before "FINISH"
	#
	#
	# Case 3: Test without reference solution
	#
	# typical fragment:
	#----
	#TEST WITHOUT REFERENCE SOLUTION
	#
	#LEVEL #NEQ   #NIT   abs. res     hmin        AR         L2err     true err    I_eff     L_oo err
	#2      256    11    8.808E-13  1.562E-02  1.000E+00   5.619E-05  -2.859E-05  1.676E+00   -
	#
	#
	#FINISH  : 06:34:55
	#----
	#
	# => Add additional link before "FINISH"
	if ($_ =~ m/^TEST (DEVIANT|WITHOUT REFERENCE SOLUTION)/) {
	    while (<INFILE>) {
		$_ =~ s/</&lt;/g;
		$_ =~ s/>/&gt;/g;
		# In case the input file is corrupted, exit and continue with outer loop
		redo LINE if ($_ =~ m/^TESTID\s+=\s+(\S+)\s*$/);

		# As soon as we reach FINISH, add link and continue with outer loop
		if ($_ =~ m/^FINISH\s+:/) {
		    $message .= qq{</pre>\n<div class="anchor">$link</div>\n<pre>\n\n};
		    $message .= $_;
		    $_ = <INFILE>;
		    redo LINE;
		} else {
		    $message .= $_;
		}
	    } # loop over input file
	}
	# Clear line as it has already been processed.
	$_ = "";
    }

    # Create an anchor at the end
    if ($_ =~ m/^SUMMARY:/) {
        chomp($_);
	$message .= qq{</pre><a name="summary"></a><pre>\n};
    }

    # Create links to conspicuous tests at the end of the email
    if ($_ =~ m/^(\S+)\s*$/  &&  exists($nextFailedTest{$1})) {
        chomp($_);
	$message .= qq{</pre>\n<div class="link"><a href="#$1">} . $_ . qq{</a></div>\n<pre>\n};
        $_ = "";
    }

    $message .= $_;
}
close(INFILE);
$message .= qq{
</pre>
</body>
</html>};

# Send message
open (MAILPRG, "| /usr/lib/sendmail $recipient") or die "$0: Cannot talk to sendmail program: $!\n";
# add extra header fields to mail:
# it is quite convenient to see PATH and other environment variable settings
# as part of the mail, convenient for debugging.
print MAILPRG "X-Cron-Env: <HOME=$ENV{HOME}>\n";
print MAILPRG "X-Cron-Env: <LOGNAME=$ENV{LOGNAME}>\n";
print MAILPRG "X-Cron-Env: <PATH=$ENV{PATH}>\n";
print MAILPRG "X-Cron-Env: <SHELL=$ENV{SHELL}>\n";
print MAILPRG $message;
close MAILPRG;

1;
