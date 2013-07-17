#!/usr/bin/env perl
#
# Script which collects the results of a FEAST benchmark executed
# via a queueing system and creates a logfile containing individual
# results and overall statistics at the end.
# Assumption: 1 test per file! See help screen for implications otherwise.


# Allow proper programming only
# (a bit like Fortran90's "implicit none")
use strict;
use warnings;

# handling of command line options
use Getopt::Long 2.33 qw(:config);

# Portably parse a pathname into directory, basename and extensions
use File::Basename;


# Some global variables
(my $progname=$0) =~ s/^.*\/(.*)/$1/;
my $debugScript = 0;

# String that marks the end of a FEAST benchmark test's log output
my $fieldDelimiter = qr{^RUNTIME :};

my $logcontent;
my (@listOfDeviantTests, @listOfUnveriTests, @listOfCrashedTests);
my ($testall, $testpass, $testdeviant, $testunverified, $testexecfailed) = (0,0,0,0,0);

# variables for command line options
my $cl_help    = "";
my $cl_version = "";


# Parse command line options
GetOptions (
    "help"    => \$cl_help,
    "version" => \$cl_version,
    );

# Show version if requested
&show_version() if ($cl_version);

# Show help if requested
&show_help() if ($cl_help);

# No argument given
unless (@ARGV) {
    die "$0: The basename of the logfiles to gather needs to be given,\n" .
	"e.g. '$0 output'. Try $0 --help for more information.\n";
}

my $basename = &basename($ARGV[0]);

# Get names of all files that match 'output.*'
my $directory = &dirname($ARGV[0]);
opendir LOCALDIR, $directory or
    die "$progname: Cannot open directory <$directory>: $!\n";
my @files = grep /^$basename.*/, readdir LOCALDIR;
closedir LOCALDIR;


# Parse those files
unless (@files) {
    die "$0: No files found with basename <$basename>.\n";
}
foreach my $file (sort @files) {
    $file = $directory . "/" . $file;
    open(INFILE, "<", $file) or
	warn "$progname: Cannot open file <$file>: $!\n";
    print STDERR "Parsing $file... " if ($debugScript);
    my $show = 1;
    my $testID = "";
    my $testfinished = 0;
    while (<INFILE>) {
	# Find out test ID
	if (m/^TESTID\s+=\s+(\S+)\s*$/) {
	    $testall++;
	    $testID = $1;
	    $show = 1;
	    $testfinished = 0;
	}

	# Copy content to standard output till delimiter is found
	print if ($show);
	if (m/$fieldDelimiter/) {
	    $show = 0;

	    # test finished somehow, either killed by queueing system
	    # because of wallclock time limit or ended by itself.
	    $testfinished = 1;
	}

	# As soon as we reached statistics section (show == 0),
	# try to find out if the test succeeded, deviated, failed or is new
	if (! $show) {
	    if (m/^tests successfully passed\s*:\s*(\d+)/) {
		$testpass += $1;
	    } elsif (m/^tests deviant\s*:\s*(\d+)/) {
		$testdeviant += $1;
		push @listOfDeviantTests, $testID if ($1 > 0);
	    } elsif (m/^tests unverifiable\s*:\s*(\d+)/) {
		$testunverified += $1;
		push @listOfUnveriTests, $testID if ($1 > 0);
	    } elsif (m/^tests failing execution\s*:\s*(\d+)/) {
		$testexecfailed += $1;
		push @listOfCrashedTests, $testID if ($1 > 0);
	    }
	}
    }
    if ($testfinished == 0) {
        # If test did not finish (= logfile truncated), add test ID
        # to list of crashed tests
	$testexecfailed += 1;
	push @listOfCrashedTests, $testID;
    }

    close(INFILE);
    print STDERR "done\n" if ($debugScript);
    print "\n";
}

# Final statistics
print "\nSUMMARY:\n\n";

print "tests in total           : $testall\n";
print "tests successfully passed: $testpass\n";
print "tests deviant            : $testdeviant\n";
print "tests unverifiable       : $testunverified\n";
print "tests failing execution  : $testexecfailed\n\n";

# Print list of failed tests
if ($testdeviant > 0) {
  print "\n";
  print "The tests with deviating results are coded in the following files:\n";
  print join("\n", @listOfDeviantTests) . "\n";
}

# Print list of tests that could not be verified due to
# missing reference solution
if ($testunverified > 0) {
  print "\n";
  print "The tests that could not be verified due to a missing\n";
  print "reference solution are coded in the following files\n";
  print join("\n", @listOfUnveriTests) . "\n";
}

# Print list of tests that crashed
if ($testexecfailed > 0) {
  print "\n";
  print "The tests that crashed during execution are coded in\n";
  print "the following files\n";
  print join("\n", @listOfCrashedTests) . "\n";
}
print "\n";


# Function which shows valid options and IDs to this script
sub show_help {
    my $counter = 0;
    print qq{        Usage: $progname [options] <basename>

A Perl script used by 
* fbenchmark2/bin/runregressiontest
to collect the results of a FEAST benchmark executed via a 
queueing system. It creates a logfile containing the non-redundant
data of all files matching <basename>.* and adds statistics at the
end about successful, deviant and crashed tests.

The current implementation implicitly assumes that every logfile
only contains the output of a single test ID run: If more than one
test ID output is contained, every test that has hit the wallclock time
limit (and hence has an incomplete log output) and that is not the
last one in a logfile is not detected as crashed test!

Available options:
--help              Print this message
--version           Show version information
};
    exit 0;
}



# Function which shows this script's version
sub show_version {
    print $progname . " v" . &get_version() . "\n";
    print "Written by Sven H.M. Buijssen.\n";

    exit 0;
}



# Determine version of this script
# (from the CVS ID in the header of this script or
#  from the hard-coded VERSION constant)
sub get_version {
    my $version = "";

    # Open this script for reading
    open(FILE, "<", $0);
    if (! eof(FILE)) {
        while (<FILE>) {
            if (m/^# \$Id: gather_logfiles_to_mail.pl,v ([\d\.]+) /) {
                $version = $1;
                last;
            }
        }
    }
    close(FILE);

    # Fall back to hard-coded version number if version number unset
    $version = "[not available]" if ($version eq "");

    return $version;
}
