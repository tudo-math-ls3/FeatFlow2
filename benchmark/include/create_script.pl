#!/usr/bin/env perl
#
# A Perl script used by the Makefile in FEAT2 Benchmark & Regression
# directory to create the benchmark control script (intuitively called runtests).
#
# This script parses an fbconf file (or more generally a given ASCII
# file containing a list of test IDs, one or multiple per line), looks
# up the settings associated with these test IDs (stored in one of the
# fbdef files in subdirectory 'tests') and creates instruction blocks
# (environment variable declarations in uppercase in sh syntax) for
# all the tests requested. These instructions are integrated into
# 'runtests'. If a test ID is found for which there is no definition
# in any of the fbdef files, it will report a warning and ignore the
# test ID.
#
# A number of special keywords in *.fbdef files are supported that are not
# stored in the 'runtests' script, but instead influence whether and how
# test IDs are included in 'runtests'. These keywords are:
#   allowhosts:
#      comma-separated list of names of hosts a benchmark test should
#      be run on.
#      To run a test on any host, specify 'all' or 'any' (without the
#      single ticks). To disable execution in general, specify 'none'.
#
#      Examples:
#         allowhosts       = all
#         allowhosts       = ashenvale,nightwish,kittyhawk,yggdrasill,oldman
#
#
#   allowbuildIDs:
#      comma-separated list of build IDs (standard FEAT2 five dash-separated tokens
#      to identify architecture, cpu, operating system, compiler suite and BLAS
#      implementation) a benchmark test should be run on. Like for 'allowhosts'
#      special values 'all', 'any' and 'none' are supported. Additionally, simple
#      Perl regular expressions are supported.
#
#      Examples:
#         allowbuildIDs   = all
#         allowbuildIDs   = pc-opteron-linux64-intel-goto
#         allowbuildIDs   = pc-.*-linux64-intel-.*
#
#      Note that regexp support is not 100% foolproof. A setting of
#         allowbuildIDs   = pc-opteron-linux64-intel-goto, pc-.*-linux64-intel-.*
#      is not undone by a single
#         denybuildIDs   = pc-.*-linux64-intel-.*
#      but only by
#         denbybuildIDs  = pc-opteron-linux64-intel-goto, pc-.*-linux64-intel-.*
#
#      More pathologic cases:
#         allowbuildIDs   = all
#         denybuildIDs    = pc-.*-linux32-.*-.*
#         allowbuildIDs   = pc-athlon64-linux32-gcc-.*
#      => should allow pc-athlon64-linux32-gcc-.* and does. Ok.
#
#         allowbuildIDs   = all
#         allowbuildIDs   = pc-.*-linux32-.*-.*
#         denybuildIDs    = pc-athlon64-linux32-gcc-.*
#      => should allow pc-athlon64-linux32-gcc-.* and does. Ok.
#
#
#   denyhosts:
#      comma-separated list of names of hosts a benchmark test should
#      never be run on. Complementary to 'allowhosts' keyword, see below.
#      Special values: 'all', 'any', 'none'
#
#      Examples:
#         denyhosts   = none
#         denyhosts   = oldman,nightwish
#
#      Remark on allow/denyhosts:
#      The following equivalences hold:
#         denyhosts   = none    <=>       allowhosts = all
#         denyhosts   = any     <=>       allowhosts = none
#      A host added via 'allowhosts' is automatically removed from the whitelist
#      when it occurs in a subsequent 'denyhosts' statement and vice versa holds
#      for the blacklist.
#
#
#   denybuildIDs:
#      comma-separated list of build IDs (standard FEAT2 five dash-separated tokens
#      to identify architecture, cpu, operating system, compiler suite, BLAS
#      implementation) a benchmark test should never be run on.
#      The mechanism works exactly like 'allowbuildIDs' and supports the same special
#      values 'all', 'any' and 'none' as well as simple Perl regular expressions.
#      The remark on the complementary nature of 'allowhosts'/'denyhosts' holds for
#      'allowbuildIDs'/'denybuildIDs', too.
#
#      Examples:
#         denybuildIDs   = pc-opteron-linux64-intel-goto
#         denybuildIDs   = pc-opteron-linux64-intel-goto, sun4u-sparcv9-sunos-.*-.*
#
#
#   execmode:
#      comma-separated list of execution modes. Valid values are 'mpi' and
#      'serial'. The use of 'parallel' is equivalent to 'mpi', but deprecated.
#      We might add 'openmp' in future.
#



# Allow proper programming only
# (a bit like Fortran90's "implicit none")
use strict;
use warnings;

# Try to import execmode from environment (set by calling Makefile)
use Env;

# handling of command line options
use Getopt::Long 2.33 qw(:config prefix_pattern=--);

# Portably create directories recursively
use File::Path qw(mkpath);

# Portably rename/move directories
use File::Copy qw(move);

# Portably try every conceivable way to determine
# the current hostname
use Sys::Hostname;


# Some constants
(my $progname=$0) =~ s/^.*\/(.*)/$1/;
my $debugScript = 0;

# Name of current host
my $host = &hostname();

# Name of directory containing the files that code
# our tests
my $testsdir = "tests";

# Extension of files encoding FEAT2 benchmark tests
my $testfileext = "fbdef";

# Immediately print messages to screen, don't buffer them,
# especially helpful when debugging.
use IO::Handle;
STDOUT->autoflush(1);

use constant VERSION => 2.12;



# Parse command line
#
# If GetOptions() returns false, the function detected one or
# more errors during option parsing. So, this script should die.
#
# Pass the variables by reference (syntax: prefix them with a backslash),
# when passing an array explicitly
my $unknownoption=0;
# variables for command line options
my %cl = ();
GetOptions (
            "append-to-files=s"             => \$cl{'append-to-files'},
            "log-base-directory=s"          => \$cl{'log-base-directory'},
            "overwrite-log-directory"       => \$cl{'overwrite-log-directory'},
            "help"                          => \$cl{'help'},
            "version"                       => \$cl{'version'},
            ) || ($unknownoption=1);

# Show version if requested
&show_version() if ($cl{'version'});

# Show help if requested
&show_help() if ($cl{'help'});

# At least one file has to be given on the command line.
# (Note: perl starts to count from zero!)
if ($#ARGV < 0) {
    &show_help();
    die "\n";
}


# Set default execution mode
$ENV{"MPI"} = "NO" unless ($ENV{"MPI"});

# --log-base-directory can also be given via env variable LOG_BASE_DIRECTORY=<string>
if ($ENV{"LOG_BASE_DIRECTORY"} && ! defined($cl{'log-base-directory'})) {
    $cl{'log-base-directory'} = $ENV{"LOG_BASE_DIRECTORY"};
}
$cl{'log-base-directory'} ||= "logs";
# --overwrite-log-directory can also be given via env variable OVERWRITE_LOG_DIRECTORY=1
if ($ENV{"OVERWRITE_LOG_DIRECTORY"} && ! $cl{'overwrite-log-directory'}) {
    $cl{'overwrite-log-directory'} = $ENV{"OVERWRITE_LOG_DIRECTORY"};
}


# Set build ID
my $buildID = $ENV{"ID"} || "";
$buildID =~ s/^(\w+-\w+-\w+-\w+-\w+)(-\w+|)$/$1/;


# Step 1:
# Find all files that code tests,
# i.e. all files that match *.$testfileext,
# but exclude defaults.$testfileext as it should be the first on the list
opendir TESTSDIR, $testsdir
    or die "\n$progname: Cannot open directory <$testsdir>: $!\n\n";
my @testfiles = sort grep { /.+\.$testfileext$/ && ! /^defaults.$testfileext$/  && ! /^\./ } readdir TESTSDIR;
closedir TESTSDIR;

# Define a hash for all the information that code a test
# Access:
#   $test{$testid}{'class'}
#   $test{$testid}{'descr'}
#   $test{$testid}{'grid'}
# etc.
my %test;

# Step 2:
# Read in default settings from defaults.$testfileext, if exists
my %default;
my @defaultSettingsFiles;
push @defaultSettingsFiles, "defaults." . $testfileext;
foreach my $testfile (@defaultSettingsFiles) {
    $testfile = $testsdir . "/" . $testfile;
    if ( -e $testfile) {
	print STDERR "Parsing test coding file <$testfile>...\n" if ($debugScript);
	open(TESTFILE, $testfile)
	    or warn "$progname: WARNING: Cannot read defaults from file <$testfile>: $!\n";
	my $lineno = 0;
	LINE: while (<TESTFILE>) {
	    $lineno++;
	    my $origline = $_;

	    # Remove any inlined comments
	    $_ =~ s/\#.+$//;

	    # Remove trailing white space
	    $_ =~ s/\s+$//;

	    # Ignore comments (i.e. lines starting with a '#')
	    next LINE if ($_ =~ m/^\s*\#/);
	    # Ignore empty lines
	    next LINE if ($_ =~ m/^\s*$/);

	    # Match a line of format
	    # keyword  =  entry
	    m%^\s*(\S+)\s*=\s*(.*)$%;

	    # Catch case where a line in a FEAT2 benchmark test definition file
	    # is invalidly formatted.
	    if (! (defined($1) && defined($2))) {
		die "\n$progname: ERROR:\n" .
		    "  Line $lineno in file <$testfile> has unknown format:\n" .
		    "  genuine line   : <$origline>\n" .
		    "  normalised line: <$_>\n\n";
	    }

	    # Store new settings as defaults for next entry.
	    #
	    # A value, however, may contain multiple keyword-value definitions
	    # following the scheme:
	    #     appl      = sbblas,MATCONST:NO,AUX1:MV
	    # We need to split it up and store them separately as keyword-value pairs:
	    #     APPL     = sbblas
	    #     MATCONST = NO
	    #     AUX1     = MV
	    # otherwise we can override discrete settings later on, e.g. override
	    # only the value of AUX1.
	    my %hash = &splitup($1, $2);
	    foreach my $keyword (keys %hash) {
		my $value = $hash{$keyword};
		# Remove leading white space
		# (any trailing white space has been remove on the complete line already)
		$value =~ s/^\s*//;

		# Store it (keyword in uppercase).
		if ($value ne "") {
		    # Special handling of ALLOW*, DENY* because they are realised as
		    # hashes, not scalar values and may require updating the
		    # complimentary setting (i.e. allowhosts = jerusalem should cancel
		    # out a previous 'denyhosts = jerusalem' setting.)
		    if ($keyword =~ m/^\b(allow|deny)(hosts|buildIDs)\b$/i) {

			# 'allowhosts'/'denyhosts' and 'allowbuildIDs'/'denybuildIDs' are
			# complimentary and hash structures, while all other values are
			# scalars.  Take additional measures to ensure correct input and
			# consistency of complimentary settings.
			&update_allowdeny_settings($keyword, $value, \%default, "Storing as new default value:");
		    } else {
			printf STDERR
			    "Storing as new default value: '" .
			    uc($keyword) . "' => '" .
			    $value. "'\n" if ($debugScript);
			$default{uc($keyword)} = $value;
		    }
		} else {
		    if (exists($default{uc($keyword)})) {
			# Unset keyword
			printf STDERR
			    "Unsetting default value for: '" .
			    uc($keyword) . "'\n" if ($debugScript);
			delete $default{uc($keyword)};
		    }
		}
	    }
	    undef %hash;
	}
	close(TESTFILE);
	print STDERR "Finished parsing test coding file <$testfile>.\n" if ($debugScript);
    } else {
	warn "$progname: WARNING: Cannot open file <$testfile>: $!\n";
    }
}
undef @defaultSettingsFiles;


# Step 3:
# Read in all test coding files and
# store the information in hash %test
my $testid = "";
my $numTestfiles = 0;
foreach my $testfile (@testfiles) {
    $numTestfiles++;
    $testfile = $testsdir . "/" . $testfile;
    print STDERR "Parsing test coding file <$testfile>...\n" if ($debugScript);
    open(TESTFILE, $testfile)
	or die "\n$progname: ERROR: Cannot open file <$testfile>: $!\n\n";

    # Hash that contains for every keyword defined througout a test coding file
    # the latest settings. Used to supplement incomplete definitions.
    my %inherited;

    # For backwards compatibility, all tests are supposed to be suitable for
    # MPI execution.
    $inherited{EXECMODE} = "serial";

    # Copy defaults to inherited
    foreach my $keyword (keys %default) {
	# Scalars can simply be copied, but using an assignment operator for
	# arrays and hashes, too, would only copy references to these arrays and
	# hashes and not create independent structures! Those can not be altered
	# independently later on. All test IDs would point to the very same,
	# single array and hash objects.
	if (ref($default{$keyword}) eq "") {
	    # Dealing with a scalar value
	    print STDERR "Copying inherited scalar value to inherited{$keyword}: $default{$keyword}\n" if ($debugScript);
	    $inherited{$keyword} = $default{$keyword};
	} elsif (ref($default{$keyword}) eq "ARRAY") {
	    # Dealing with an array
	    print STDERR "Copying inherited array value to inherited{$keyword}: $default{$keyword}\n" if ($debugScript);
	    @{ $inherited{$keyword} } = ();
	    foreach my $token (@{ $default{$keyword} }) {
		push @{ $inherited{$keyword} }, $token;
	    }
	} elsif (ref($default{$keyword}) eq "HASH") {
	    # Dealing with a hash
	    print STDERR "Copying inherited hash value to inherited{$keyword}: $default{$keyword}\n" if ($debugScript);
	    %{ $inherited{$keyword} } = ();
	    foreach my $key (keys %{ $default{$keyword} }) {
		$inherited{$keyword}{$key} = $default{$keyword}{$key};
	    }
	}

    }

    my $lineno = 0;
    LINE: while (<TESTFILE>) {
	$lineno++;
	my $origline = $_;

	# Remove any inlined comments
	$_ =~ s/\#.+$//;

	# Remove trailing white space
	$_ =~ s/\s+$//;

	# Ignore comments (i.e. lines starting with a '#')
	next LINE if ($_ =~ m/^\s*\#/);
	# Ignore empty lines
	next LINE if ($_ =~ m/^\s*$/);

	# Store test ID
	if ($_ =~ m%^\s*testid\s*=\s*(.+)$%i) {
	    $testid = $1;
	    # Init settings with defaults (= latest settings)
	    foreach my $entry (sort keys %inherited) {
		# Scalars can simply be copied, but using an assignment operator for
		# arrays and hashes, too, would only copy references to these arrays and
		# hashes and not create independent structures! Those can not be altered
		# independently later on. All test IDs would point to the very same,
		# single array and hash objects.
		if (ref($inherited{$entry}) eq "") {
		    # Dealing with a scalar value
		    print STDERR "Copying inherited scalar value to test{$testid}{$entry}: $inherited{$entry}\n" if ($debugScript);
		    $test{$testid}{$entry} = $inherited{$entry};
		} elsif (ref($inherited{$entry}) eq "ARRAY") {
		    # Dealing with an array
		    print STDERR "Copying inherited array value to test{$testid}{$entry}: $inherited{$entry}\n" if ($debugScript);
		    @{ $test{$testid}{$entry} } = ();
		    foreach my $token (@{ $inherited{$entry} }) {
			push @{ $test{$testid}{$entry} }, $token;
		    }
		} elsif (ref($inherited{$entry}) eq "HASH") {
		    # Dealing with a hash
		    print STDERR "Copying inherited hash value to test{$testid}{$entry}: $inherited{$entry}\n" if ($debugScript);
		    %{ $test{$testid}{$entry} } = ();
		    foreach my $key (keys %{ $inherited{$entry} }) {
			$test{$testid}{$entry}{$key} = $inherited{$entry}{$key};
		    }
		}
	    }
	} else {
	    # Match a line of format
	    # keyword  =  entry
	    m%^\s*(\S+)\s*=\s*(.*)$%;

	    # Catch case where a line in a FEAT2 benchmark test definition file
	    # is invalidly formatted.
	    if (! (defined($1) && defined($2))) {
		die "\n$progname: ERROR:\n" .
		    "  Line $lineno in file <$testfile> has unknown format:\n" .
		    "  genuine line   : <$origline>\n" .
		    "  normalised line: <$_>\n\n";
	    }

	    # Store new settings as defaults for next entry.
	    #
	    # A value, however, may contain multiple keyword-value definitions
	    # following the scheme:
	    #     appl      = sbblas,MATCONST:NO,AUX1:MV
	    # We need to split it up and store them separately as keyword-value pairs:
	    #     APPL     = sbblas
	    #     MATCONST = NO
	    #     AUX1     = MV
	    # otherwise we can override discrete settings later on, e.g. override
	    # only the value of AUX1.
	    my %hash = &splitup($1, $2);
	    foreach my $keyword (keys %hash) {
		my $value = $hash{$keyword};
		# Remove leading white space
		# (any trailing white space has been remove on the complete line already)
		$value =~ s/^\s*//;

		if ($testid ne "") {
		    if ($value ne "") {
			# Special handling of ALLOW*, DENY* because they are realised as
			# hashes, not scalar values and may require updating the
			# complimentary setting (i.e. allowhosts = jerusalem should cancel
			# out a previous 'denyhosts = jerusalem' setting.)
			if ($keyword =~ m/^\b(allow|deny)(hosts|buildIDs)\b$/i) {
			    # 'allowhosts'/'denyhosts' and 'allowbuildIDs'/'denybuildIDs' are
			    # complimentary and hash structures, while all other values are
			    # scalars.  Take additional measures to ensure correct input and
			    # consistency of complimentary settings.
			    &update_allowdeny_settings($keyword, $value, \%{ $test{$testid} }, "Storing in test{$testid}{$keyword}:");
			} else {
			    # Store it for current ID (keyword in uppercase).
			    print STDERR "Storing in test{$testid}{$keyword}: $value\n" if ($debugScript);
			    $test{$testid}{uc($keyword)} = $value;
			}
		    } else {
			if (exists($test{$testid}{uc($keyword)})) {
			    # Unset keyword
			    printf STDERR
				"Unsetting test{$testid}{$keyword}\n" if ($debugScript);
			    delete $test{$testid}{uc($keyword)};
			}
		    }
		}

		if ($value ne "") {
		    # Store new settings as defaults for next entry
		    printf STDERR "Storing as new inherited value: '$keyword' => $value\n" if ($debugScript);
		    if ($keyword =~ m/^\b(allow|deny)(hosts|buildIDs)\b$/i) {
			# 'allowhosts'/'denyhosts' and 'allowbuildIDs'/'denybuildIDs' are
			# complimentary. Take additional measures to ensure consistency
			# of settings.
			&update_allowdeny_settings($keyword, $value, \%inherited, "Storing as new inherited value:");
		    } else {
			$inherited{uc($keyword)} = $value;
		    }
		} else {
		    if (exists($inherited{uc($keyword)})) {
			# Unset keyword
			delete $inherited{uc($keyword)};
		    }
		}
	    }
	    undef %hash;
	}
    }
    # Don't inherit defaults from other *.$testfileext files
    undef %inherited;
    $testid = "";

    close(TESTFILE);
    print STDERR "Finished parsing test coding file <$testfile>.\n" if ($debugScript);
}

if ($numTestfiles == 0) {
    die "\n$progname: ERROR:\n" .
	"  Not a single file matching *.$testfileext has been found in directory\n" .
	"  <$testsdir/>. Without such FEAT2 benchmark test definition\n" .
	"  files it is not possible to create a FEAT2 benchmark run script.\n\n";
}


if ($debugScript) {
    print STDERR "\nWhat has been stored in hash %test:\n";
    # Print all entries in hash %test
    foreach my $testid (sort keys %test) {
	foreach my $entry (sort keys %{ $test{$testid} }) {
	    print STDERR "test{$testid}{$entry} = " . $test{$testid}{$entry} . "\n";
	}
    }
}


# Step 4:
# Read in the IDs in the file given as argument
# and create an executable script
my @idsToCode = ();
foreach my $file (@ARGV) {
    print STDERR "Parsing file <$file>...\n" if ($debugScript);
    open(FILE, $file)
	or warn "$progname: WARNING: Cannot open file <$file>: $!\n";

   LINE2: while (<FILE>) {
	# Ignore comments (i.e. lines starting with a '#')
	next LINE2 if ($_ =~ m/^\#/);
	# Ignore empty lines
	next LINE2 if ($_ =~ m/^\s*$/);

	# Remove any inlined comments
	$_ =~ s/\#.+$//;

	# Remove white spaces (trim)
	chomp($_);
	$_ =~ s/^\s*(\S+)\s*$/$1/;

	print STDERR "Found ID <" . $_ . ">.\n" if ($debugScript);
	push @idsToCode, split('\s+', $_);
    }
    close(FILE);
    print STDERR "Finished parsing file <$file>.\n" if ($debugScript);
}


# Step 5:
my @subentry    = ();
my @subsubentry = ();
my $testidsFound = 0;
my @vars2export = ();
ID: foreach my $testid (@idsToCode) {
    # Check whether anything at all has been defined in any of
    # the *.$testfileext files for the requested ID.
    if (scalar(keys %{ $test{$testid} }) == 0) {
	    warn "$progname: WARNING:\n" .
		 "  Test case with ID <$testid> is undefined!\n";
	    next ID;
    }

    # Check whether a minimal set of keywords has been defined
    # in any of the *.$testfileext files for the requested ID.
    foreach my $check ( 'CLASS', 'DESCR', 'APPL', 'MGLEVELS', 'DATFILE' ) {
	unless (defined($test{$testid}{$check})) {
	    warn "$progname: WARNING:\n" .
		 "  Test case with ID <$testid> is incompletely defined!\n" .
		 "  Keyword <$check> not found.\n";
	    next ID;
	}
    }


    if (exists($test{$testid}{'HOSTPERMISSIONS'})) {
	my $allowed = 0;
	foreach my $entry ( @{ $test{$testid}{'HOSTPERMISSIONS'} } ) {
	    $entry =~ m/^(.)(.*)$/;
	    my ($perm, $regexp) = ($1, $2);

	    if ($perm eq "+") {
		$perm = 1;
	    } else {
		$perm = 0;
	    }
	    print STDERR "Splitting up $entry into $perm and $regexp.\n" if ($debugScript);

	    if ($regexp eq "all" || $host =~ m/$regexp/) {
		$allowed = $perm;
	    }
	}

	if (! $allowed) {
	    my ($allowedHosts, $deniedHosts) = ("", "");
	    foreach my $entry ( @{ $test{$testid}{'HOSTPERMISSIONS'} } ) {
		$entry =~ m/^(.)(.*)$/;
		my ($perm, $regexp) = ($1, $2);

		if ($perm eq "+") {
		    $allowedHosts .= ", " . $regexp;
		} else {
		    $deniedHosts .= ", " . $regexp;
		}
	    }
	    $allowedHosts =~ s/^, //;
	    $deniedHosts  =~ s/^, //;

	    warn "$progname: WARNING:\n" .
	    	"  Test case with ID <$testid> is coded not to be run on current host, <$host>:\n" .
	    	"  list of allowed hosts: " . $allowedHosts . "\n" .
	    	"  list of denied hosts : " . $deniedHosts . "\n" .
	    	"  Test will be skipped.\n";
	    next ID;
	}
    }


    if (exists($test{$testid}{'BUILDIDPERMISSIONS'})) {
	my $allowed = 0;
	foreach my $entry ( @{ $test{$testid}{'BUILDIDPERMISSIONS'} } ) {
	    $entry =~ m/^(.)(.*)$/;
	    my ($perm, $regexp) = ($1, $2);

	    if ($perm eq "+") {
		$perm = 1;
	    } else {
		$perm = 0;
	    }
	    print STDERR "Splitting up $entry into $perm and $regexp.\n" if ($debugScript);

	    if ($regexp eq "all" || $buildID =~ m/$regexp/) {
		$allowed = $perm;
	    }
	}

	if (! $allowed) {
	    my ($allowedBuildIDs, $deniedBuildIDs) = ("", "");
	    foreach my $entry ( @{ $test{$testid}{'BUILDIDPERMISSIONS'} } ) {
		$entry =~ m/^(.)(.*)$/;
		my ($perm, $regexp) = ($1, $2);

		if ($perm eq "+") {
		    $allowedBuildIDs .= ", " . $regexp;
		} else {
		    $deniedBuildIDs .= ", " . $regexp;
		}
	    }
	    $allowedBuildIDs =~ s/^, //;
	    $deniedBuildIDs  =~ s/^, //;

	    warn "$progname: WARNING:\n" .
		"  Test case with ID <$testid> is coded not to be run for current build ID, <$buildID>:\n" .
		"  list of allowed build IDs: " . $allowedBuildIDs . "\n" .
		"  list of denied build IDs : " . $deniedBuildIDs . "\n" .
	    	"  Test will be skipped.\n";
	    next ID;
	}
    }

    # For backwards compatibility, still support EXECMODE=PARALLEL
    # but translate it silently to EXECMODE=MPI
    $test{$testid}{EXECMODE} =~ s/\bparallel\b/mpi/gi;


    # Check whether definition of requested ID has an execmode
    # that matches current one.
    my $requestedExecmode = "";
    if ("YES" =~ m/\b$ENV{"MPI"}\b/i) {
	$requestedExecmode = "mpi";
    }
    elsif ("NO" =~ m/\b$ENV{"MPI"}\b/i) {
	$requestedExecmode = "serial";
    } else {
	warn "$progname: WARNING:\n" .
	    "  Environment variable MPI is set to <" . $ENV{"MPI"} . ">. Only YES and NO are permitted..\n" .
	next ID;
    }
    if ($test{$testid}{EXECMODE} !~ m/\b$requestedExecmode\b/i) {
	warn "$progname: WARNING:\n" .
	    "  Test case with ID <$testid> is coded to work with execmode <$test{$testid}{EXECMODE}>.\n" .
	    "  Requested, however, is for execmode <$requestedExecmode>. Test will be skipped.\n";
	next ID;
    } else {
	# A single runtests script can only do either MPI or serial tests
	# So, if multiple execmodes are given, override value with the currently active one.
	$test{$testid}{EXECMODE} = $requestedExecmode;
    }


    # Write to screen or append to file?
    if ($cl{'append-to-files'}) {
	my $filename = $cl{'append-to-files'} . $testid;
	open(STDOUT, '>>', $filename)
	    or warn "$progname: WARNING: Cannot append to file <$filename>: $!\n";
    }

    $testidsFound++;
    # Re-initialise array that stores all environment variables set
    @vars2export = ();

    print STDOUT "# ================================================================\n\n";

    # Description of individual test case
    print STDOUT "echo '# ====================='\n";
    print STDOUT "echo '# Test description:'\n";
    print STDOUT "echo '# " . $test{$testid}{DESCR} . "'\n";
    print STDOUT "echo\n";

    # Set log directory if not already set
    $test{$testid}{LOGDIR} ||= $cl{'log-base-directory'} . "/" . $testid;

    # Create the log directory FEAT2 will use for this test ID -
    # if the directory does not already exist
    if ( -d $test{$testid}{LOGDIR} ) {
	# Directory does already exist.
	# Rename or overwrite?
	unless ($cl{'overwrite-log-directory'}) {
	    if ( -d "$test{$testid}{LOGDIR}.prev" ) {
		unlink <$test{$testid}{LOGDIR}.prev/*>;
		rmdir $test{$testid}{LOGDIR} . ".prev" ||
		    die "\n$progname: ERROR:\n" .
		        "  Cannot remove directory <$test{$testid}{LOGDIR}.prev>: $!\n";
	    }
	    # rename existing (clobbering any existing *.prev directory) and create a new one.
	    File::Copy::move("$test{$testid}{LOGDIR}", "$test{$testid}{LOGDIR}.prev") ||
		die "\n$progname: ERROR:\n" .
		    "  This script tried to rename a directory named <" . $test{$testid}{LOGDIR} . ">,\n" .
		    "  but an error occured: $?\n\n";
	    File::Path::mkpath($test{$testid}{LOGDIR}) ||
		die "\n$progname: ERROR:\n" .
		    "  This script tried to create a directory named <" . $test{$testid}{LOGDIR} . ">,\n" .
		    "  but an error occured: $?\n\n";
	}
    } else {
	# Directory did not yet exist.

	if ( -e $test{$testid}{LOGDIR} ) {
	    # some file system object with the desired name
	    # already exists and it is not a directory. Issue an error.
	    die "\n$progname: ERROR:\n" .
		"  This script was about to create a directory named <" . $test{$testid}{LOGDIR} . ">.\n" .
		"  But there exists already such a file system object and it is no directory. Please check.\n\n";
	} else {
	    File::Path::mkpath($test{$testid}{LOGDIR}) ||
		die "\n$progname: ERROR:\n" .
		    "  This script tried to create a directory named <" . $test{$testid}{LOGDIR} . ">,\n" .
		    "  but an error occured: $?\n\n";
	}
    }


    # Pure cosmetics: Determine length of longest environment variable
    my $fieldlength = &get_max_col_length( keys %{ $test{$testid} } ) + 1;

    # Special treatment for test id
    # (as we want it to appear as first item)
    &formatted_print_env_variable("TESTID", $testid, $fieldlength);
    push @vars2export, "TESTID";


    # Export all remaining items but the special keywords
    # ('buildidpermissions' 'hostpermissions', 'execmode')
    foreach my $entry (sort keys %{ $test{$testid} }) {
	# Skip entries we already handled
	next if ("DESCR TESTID EXECMODE BUILDIDPERMISSIONS HOSTPERMISSIONS" =~ m/\b$entry\b/i);

	&formatted_print_env_variable($entry, $test{$testid}{$entry}, $fieldlength);
	push @vars2export, $entry;
    }


    # Finally, export execmode. This has to be done at the end as it might be necessary
    # to override previously set variables because of the execmode.

    # WARNING:
    # The scripts used to set up and start batch job scripts (bin/*schedule_*tests)
    # use 'sed' regular expression substitutions to dynamically set/override the
    # EXECMODE value when creating the batch script.
    # Whenever changing the syntax here, do not forget to adapt them!
    if ($test{$testid}{EXECMODE} =~ m/\bmpi\b/i) {
	# MPI case
	&formatted_print_env_variable("EXECMODE", $test{$testid}{EXECMODE}, $fieldlength);
    } elsif ($test{$testid}{EXECMODE} =~ m/\bserial\b/i) {
	# serial case
	&formatted_print_env_variable("EXECMODE", $test{$testid}{EXECMODE}, $fieldlength);

    } else {
	# unknown case
	&formatted_print_env_variable("EXECMODE", $test{$testid}{EXECMODE}, $fieldlength);
	warn "$progname: WARNING:\n" .
	    "  Execution mode <$test{$testid}{EXECMODE}> is requested for test case with ID <$testid>.\n" .
	    "  This is an unhandled case, so it might happen that things don't work as intended!\n";
    }

    my $length = length("export");
    print STDOUT "export";
    foreach my $entry (sort keys %{ $test{$testid} }) {
	unless ($entry =~ m/\b(HOST|BUILDID)PERMISSIONS\b/i) {
	    $length += length($entry) + 1;

	    if ($length < 80) {
		print STDOUT " $entry";
	    } else {
		print STDOUT "\nexport $entry";
		$length = length("export");
	    }
	}
    }
    print STDOUT "\n";
    print STDOUT "export TESTID\n";

    # All non-interactive MPI jobs need to have the environment variables
    # FEAT2 uses in its master.dat's explicitly set with the 'mpirun' command.
    print STDOUT qq{\nvars2export="\$vars2exportAlways } . join(" ", @vars2export) . qq{"\n};

    # Applications are named according the scheme:
    #   <basename>-<application>
    # where <application> is both found in src_<application> and as value of the
    # keyword 'appl' in a *.$testfileext file in subdirectory 'tests'.
    print STDOUT "\nappname=\${APPL_BASENAME}-\${APPL}\n";

    # Include the file that actually starts and evaluate a test
    print STDOUT "\nfb_runAndEvalBenchmarkTest\n\n";

    # Now, explicitly unset all previously set environment variables to
    # provide a clean environment for subsequent tests.
    print STDOUT "# Now, explicitly unset all previously set environment variables to\n" .
	         "# provide a clean environment for subsequent tests.\n";
    $length = length("unset");
    print STDOUT "unset";
    foreach my $entry (sort keys %{ $test{$testid} }) {
	unless ($entry =~ m/\b(HOST|BUILDID)PERMISSIONS\b/i) {
	    $length += length($entry) + 1;

	    if ($length < 80) {
		print STDOUT " $entry";
	    } else {
		print STDOUT "\nunset $entry";
		$length = length("unset");
	    }
	}
    }
    print STDOUT "\nunset TESTID\n\n";

    if ($cl{'append-to-files'}) {
	close(STDOUT);
    }
}

if ($testidsFound == 0) {
    die "\n\n$progname: ERROR:\n" .
	"  Not a single test case found in " . join(' ', @ARGV) . "\n" .
	"  is found to be valid. Either no definition in one\n" .
	"  of the following files:\n  * " .
	join("\n  * ", sort @testfiles) . "\n" .
	"  could be found or the current hostname is not listed as valid host\n" .
	"  for any of the test IDs.\n\n";
}


### Internals functions ###

# Formatted print of an environment variable and its value
# WARNING:
# The scripts used to set up and start batch job scripts (bin/*schedule_*tests)
# use 'sed' regular expression substitutions to dynamically set/override
# the EXECMODE exported to runtests with this function.
# Whenever changing the syntax here, do not forget to adapt them!
sub formatted_print_env_variable {
    my ($variablename, $value, $fieldlength) = (@_);

    print STDOUT
	# Print environment variable definition: keyword=value
	$variablename . "=" . $value . "\n";
    return;
}


# Determine length of longest item in given list
sub get_max_col_length {
    my @cols = @_;
    my $maxlength = 0;

    foreach my $entry (@cols) {
	$maxlength = &max($maxlength, length($entry));
    }

    return $maxlength;
}


# Determine maximum of a list of numbers
sub max {
    my $max = shift(@_);
    foreach my $foo (@_) {
	$max = $foo if $max < $foo;
    }
    return $max;
}


# Take a keyword/value pair, detect any keyword/value pairs contained
# within and return a hash of keyword/value pairs.
# Note: Value must not contain any whitespaces to distinguish
#       keyword/value pair from free text fields like 'description'
#       which also may contain commas and colons.
#
# Intention: Handle cases like the following:
#     appl      = sbblas,MATCONST:NO,AUX1:MV
# We need to split it up and store them separately as keyword-value pairs:
#     APPL     = sbblas
#     MATCONST = NO
#     AUX1     = MV
#
# On the other hand, values that contain commas, but no colons should not
# be split up. Think about cases like:
#     descr         = fictitious boundary tests, NCC (4 processes, 4 macros), low mg level
#     execmode      = mpi, serial    # deprecated: parallel, serial
sub splitup {
    my ($keyword, $value) = (@_);
    my %hash = ();

    # Handle
    if ($value =~ m/.+,\S+:\S+/) {
	# Split up entry
	my @subentry = split(',', $value);
	# Handle special case of first item where keyword is already known.
	# Store it (keyword in uppercase).
	$hash{$keyword} = $subentry[0];

	# If the value does have subitems, export them as well
	# (again keyword in uppercase).
	for (my $i = 1; $i <= $#subentry; $i++) {
	    my @subsubentry = split(':', $subentry[$i]);
	    $hash{$subsubentry[0]} = $subsubentry[1];
	}
    } else {
	$hash{uc($keyword)} = $value;
    }

    return %hash;
}


sub findInArray () {
    my ( $arrayRef, $string ) = @_;

    for (my $i = 0; $i <= $#{ $arrayRef }; $i++) {
	if (${ $arrayRef }[$i] eq $string) {
	    return $i;
	}
    }
    return -1;
}


sub addToArray () {
    my ( $arrayRef, $string ) = @_;

    # Add only if not already on the list
    return if ( grep { $_ eq $string} @{ $arrayRef } );

    push @{ $arrayRef }, $string;
    return;
}


sub clearArray () {
    my ( $arrayRef ) = @_;

    @{ $arrayRef } = ();
    return;
}


sub removeIthTokenFromArray () {
    my ( $arrayRef, $string, $i ) = @_;

    if ($i >= 0) {
	splice @{ $arrayRef }, $i, 1;
    } else {
	warn "$0: removeFromArray(): Could not locate $i-th string <$string> and remove from array.\n";
    }
    return;
}


sub removeFromArray () {
    my ( $arrayRef, $string ) = @_;

    my $i = &findInArray($arrayRef, $string);
    &removeIthTokenFromArray($arrayRef, $string, $i);
    return;
}


sub update_allowdeny_settings {
    my ($keyword, $value, $hashref, $message) = @_;

    # 'ALLOWHOSTS=...' given.
    if ($keyword =~ m/^\ballowhosts\b$/i) {
	# allowhosts = all|any => 1) discard any previous allowhosts setting,
	#                            replace with "all"
	#                         2) set denyhosts to "none"
	# allowhosts = none    => 1) discard any previous allowhosts setting,
	#                            replace with "none"
	#                         2) set denyhosts to "all"
	# allowhosts = tesla, previous setting neither one of "all|any|none":
	#                      => 1) append
	#                      => 2) remove tesla from denyhosts list, if present

	foreach my $entry (split(/\s*,\s*/, lc($value))) {
	    # Deal with special values 'all' / 'any'
	    if ($entry eq "all" || $entry eq "any") {
		# AllowHosts: all
		&clearArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} });
		&addToArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} }, "+all");

	    # Deal with special values 'none'
	    } elsif ($entry eq "none") {
		# AllowHosts: none <=> DenyHosts: all
		&clearArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} });
		&addToArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} }, "-all");
	    }

	    # A host is given, append/replace/don't append as appropriate.
	    else {
		# Remove host from denylist (if present), add it to allowlist
		my $index = &findInArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} }, "-" . $entry);
		&removeIthTokenFromArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} }, "-" . $entry, $index) if ($index >= 0);
		&addToArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} }, "+" . $entry);
	    }
	}
    }

    # 'DENYHOSTS=...' given.
    if ($keyword =~ m/^\bdenyhosts\b$/i) {
        # denyhosts = all|any => 1) discard any previous denyhosts setting,
        #                           replace with "all"
        #                        2) set allowhosts to "none"
        # denyhosts = none    => 1) discard any previous denyhosts setting,
        #                           replace with "none"
        #                        2) set allowhosts to "all"
        # denyhosts = tesla, previous setting neither one of "all|any|none":
        #                      => 1) append
        #                      => 2) remove tesla from allowhosts list, if present

	foreach my $entry (split(/\s*,\s*/, lc($value))) {
	    # Deal with special values 'all' / 'any'
	    if ($entry eq "all" || $entry eq "any") {
		# AllowHosts: all
		&clearArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} });
		&addToArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} }, "-all");

	    # Deal with special values 'none'
	    } elsif ($entry eq "none") {
		# AllowHosts: none <=> DenyHosts: all
		&clearArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} });
		&addToArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} }, "+all");
	    }

	    # A host is given, append/replace/don't append as appropriate.
	    else {
		# Remove host from allowlist (if present), add it to denylist
		my $index = &findInArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} }, "+" . $entry);
		&removeIthTokenFromArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} }, "+" . $entry, $index) if ($index >= 0);
		&addToArray(\@{ ${ $hashref }{'HOSTPERMISSIONS'} }, "-" . $entry);
	    }
	}
    }


    # 'ALLOWBUILDIDS=...' given.
    if ($keyword =~ m/^\ballowBuildIds\b$/i) {
	# allowBuildIds = all|any => 1) discard any previous allowBuildIds setting,
	#                               replace with "all"
	#                            2) set denyBuildIds to "none"
	# allowbuildids = none    => 1) discard any previous allowBuildIds setting,
	#                               replace with "none"
	#                            2) set denyBuildIds to "all"
	# # allowBuildIds = pc-opteron-linux64-.*.-.*, previous setting "all|any":
	# #                         => 1) append to list to support cases like
	# #                                denyBuildIDs = all
	# #                                allowBuildIDs  = pc-.*-linux64-intel-.*
	# #                                denyBuildIDs = pc-athlon64-linux64-intel-.*
	# #                               which forbids all Linux64 build IDs using Intel
	# #                               compiler except for on Athlon64 cpus.
	# #                            2) remove pc-opteron-linux64-.*.-.* from denyBuildIds list, if present
	# # allowBuildIds = pc-opteron-linux64-.*.-.*, previous setting "none":
	# #                         => 1) replace
	# #                            2) remove pc-opteron-linux64-.*.-.* from denyBuildIds
	# #                               list, if present
	# allowBuildIds = pc-opteron-linux64-.*.-.*, previous setting neither one of "all|any|none":
	#                         => 1) append
	#                         => 2) remove pc-opteron-linux64-.*.-.* from denyBuildIds
	#                               list, if present

	foreach my $entry (split(/\s*,\s*/, lc($value))) {
	    # Deal with special values 'all' / 'any'
	    if ($entry eq "all" || $entry eq "any") {
		# AllowBuildIDs: all
		&clearArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} });
		&addToArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} }, "+all");

	    # Deal with special values 'none'
	    } elsif ($entry eq "none") {
		# AllowBuildIDs: none <=> DenyBuildIDs: all
		&clearArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} });
		&addToArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} }, "-all");
	    }

	    # A host is given, append/replace/don't append as appropriate.
	    else {
		# Remove host from denylist (if present), add it to allowlist
		my $index = &findInArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} }, "-" . $entry);
		&removeIthTokenFromArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} }, "-" . $entry, $index) if ($index >= 0);
		&addToArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} }, "+" . $entry);
	    }
	}
    }

    # 'DENYBUILDIDS=...' given.
    if ($keyword =~ m/^\bdenyBuildIds\b$/i) {
	# denyBuildIds = all|any => 1) discard any previous denyBuildIds setting,
	#                              replace with "all"
	#                           2) set allowBuildIds to "none"
	# denyBuildIds = none    => 1) discard any previous denyBuildIds setting,
	#                             replace with "none"
	#                           2) set allowBuildIds to "all"
	# # denyBuildIds = pc-opteron-linux64-.*.-.*, previous setting "all|any":
	# #                        => 1) append to list to support cases like
	# #                               allowBuildIDs = all
	# #                               denyBuildIDs  = pc-.*-linux64-intel-.*
	# #                               allowBuildIDs = pc-athlon64-linux64-intel-.*
	# #                              which forbids all Linux64 build IDs using Intel
	# #                              compiler except for on Athlon64 cpus.
	# #                           2) remove pc-opteron-linux64-.*.-.* from allowBuildIds list, if present
	# # denyBuildIds = pc-opteron-linux64-.*.-.*, previous setting "none":
	# #                        => 1) replace
	# #                           2) remove pc-opteron-linux64-.*.-.* from allowBuildIds
	# #                              list, if present
	# denyBuildIds = pc-opteron-linux64-.*.-.*, previous setting neither one of "all|any|none":
	#                        => 1) append
	#                        => 2) remove pc-opteron-linux64-.*.-.* from allowBuildIds
	#                              list, if present

	foreach my $entry (split(/\s*,\s*/, lc($value))) {
	    # Deal with special values 'all' / 'any'
	    if ($entry eq "all" || $entry eq "any") {
		# AllowBuildIDs: all
		&clearArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} });
		&addToArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} }, "-all");

	    # Deal with special values 'none'
	    } elsif ($entry eq "none") {
		# AllowBuildIDs: none <=> DenyBuildIDs: all
		&clearArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} });
		&addToArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} }, "+all");
	    }

	    # A host is given, append/replace/don't append as appropriate.
	    else {
		# Remove host from allowlist (if present), add it to denylist
		my $index = &findInArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} }, "+" . $entry);
		&removeIthTokenFromArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} }, "+" . $entry, $index) if ($index >= 0);
		&addToArray(\@{ ${ $hashref }{'BUILDIDPERMISSIONS'} }, "-" . $entry);
	    }
	}
    }

    return;
}


# Determine version of this script
# (from the CVS ID in the header of this script or
#  from the hard-coded VERSION constant)
sub get_version {
    my $version = "";
    my $additional = "";

    # Open this script for reading
    open(FILE, "<", $0);
    if (! eof(FILE)) {
        while (<FILE>) {
            if (m/^# \$Id: create_script.pl ([\d\.]+) /) {
                $version = $1;
                last;
            }
        }
    }
    close(FILE);

    # Fall back to hard-coded version number if version number unset
    $version = VERSION if ($version eq "");

    return $version;
}


# Function which shows this script's version
sub show_version {
    print $progname . " revision " . &get_version() . "\n";
    print "Written by Sven H.M. Buijssen.\n";

    exit 0;
}


# Function which shows valid options and IDs to this script
sub show_help {
    print
        "Usage: $progname [options] <file>\n" .
	"\n" .
	"where <file> is an ASCII file containing a list of test IDs.\n" .
	"\n" .
	"This script parses an fbconf file (or more generally a given ASCII\n" .
	"file containing a list of test IDs, one or multiple per line), looks\n" .
	"up the settings associated with these test IDs (stored in one of the\n" .
	"fbdef files in subdirectory 'tests') and creates instruction blocks\n" .
	"(environment variable declarations in uppercase in sh syntax) for\n" .
	"all the tests requested. These instructions are integrated into\n" .
	"'runtests'. If a test ID is found for which there is no definition\n" .
	"in any of the fbdef files, it will report a warning and ignore the\n" .
	"test ID\n" .
	"\n" .
        "Command line options:\n" .
        "---------------------\n" .
        "--append-to-files <string>\n" .
        "                    Instead of writing to screen, append output to a\n" .
	"                    file, one per test ID. File names are constructed\n" .
	"                    according to\n" .
	"                        <string>.<test ID>\n" .
        "--log-base-directory <directory>\n" .
        "                    Name of base directory for all log files.\n" .
	"                    Default: LOG_BASE_DIR variable from Makefile\n" .
        "--overwrite-log-directory\n" .
        "                    Every test ID writes its logs to separate directory.\n" .
	"                    By default, if such a directory does already exist,\n" .
	"                    the existing one is renamed to \"<test ID>.prev\".\n" .
	"                    Specifying this option this behaviour can be turned off\n" .
	"                    Alternatively, the flag can be set by means of the\n" .
	"                    environment variable OVERWRITE_LOG_DIRECTORY=1.\n" .
        "--help              Print this message\n" .
        "--version           Show version information\n" .
    exit 0;
}
