#!/usr/bin/env perl
#
# Perl script to generate one or multiple tests.
#
# mage_test.pl: Command line script to generate test(s) from .FBDEF files.
#
# Call:
#
#     mage_test.pl [test_id] [test_id] [test_list.fbconf] ... OPTION=value
#
# In:
#
#     test_id          = one or multiple test id's to generate
#     test_list.fbconf = one or multiple .fbconf files with test is's.
#     OPTION           = one or multiple options to set in the tests.
#                        The options are specified in the form "OPTION=value".

# Global settings for the scripts to generate
my %ENVIRONMENT = ( "DEBUGGER" => 0,
                    "KEEPGOING" => 0,
                    "MPIENV" => "",
                    "OVERWRITE_LOG_DIRECTORY" => "",
                    "SUPPRESSOUTPUT" => "1",
                    
                    # Filename of the output file.
                    "SCRIPT" => "runtests",
                    
                    # Create one file for each test id.
                    # If set to 0, all tests are written ton one script.
                    "MULTIPLEFILES" => "0",
                    
                    # Temporary filename with test id's
                    "TEMPIDFILE" => "_tmp_.fbconf"
                    
                  );
                  
# List of directories where .fbconf and .fbdef file are searched.
my @DIRECTORIES = qw ( . ./tests );

# Basename of all FEAT2 benchmark applications
$APPBASE = "feat2benchmark";

# Read the template scrtips and do some variable replacements.
# All variables in the script marked with "<SET_VIA_MAKEFILE>" are
# replaced with the values specified in %ENVIRONMENT.
#
# Call:
#    @template = read_template_script($filename)
#
# with
#    $filename = name of the template script.
#
# Returns
#    @template = The template script with the environment variables
#                from %ENVIRONMENT expanded.
sub read_template_script ($) {
  my ($filename) = @_;
  
  # Read the file -- quick and dirty.
  local(*FHANDLE);
  open(FHANDLE, $filename) 
  or die("$0: ERROR. Could not open file <$filename>!\n");
  @raw_data=<FHANDLE>; 
  close(FHANDLE);
  
  # Loop through all lines. If a line "ABC=<SET_VIA_MAKEFILE>" is found,
  # "<SET_VIA_MAKEFILE>" is replaced by the corresponding value in
  # %ENVIDONMENT.
  foreach (@raw_data) {
    
    # Remove line break.
    chomp ($_);
    
    # Check if this is an option we have to set via %ENVIRONMENT.
    if ($_ =~ m/\s*(\S*)\w*=\w*\<SET_VIA_MAKEFILE\>/) {
      
      # Replace by a value from %ENVIRONMENT
      if (exists($ENVIRONMENT{$1})) {
        my $value = $ENVIRONMENT{$1};
        $_ =~ s/\<SET_VIA_MAKEFILE\>/$value/;
      }
      else {
        die "$0: ERROR. makefile option <$1> not found in \<$_\>.\n";
      }
      
    }
  }
  
  return @raw_data;
}

# Searches the current environment for the variables in
# %ENVIRONMENT, reads the variables and overwrites default values
# in %ENVIRONMENT.
sub get_environment () {
  foreach my $key (%ENVIRONMENT) {
    # Does $key exist in the current environment?`
    if (exists($ENV{$key})) {
      # Overwrite the default value
      $ENVIRONMENT{$key} = $ENV{$key};
    }
  }
}

# Searches for a file in the folders specified by @DIRECTORIES and th
# current directory.
#
# Call:
#    $path = find_file ($filename)
#
# with
#    $filename = file to search
#    $path     = fully qualified path to the file.
#                "" is returned if the file is not found.
sub find_file ($) {
  my ($filename) = @_;
  
  # Is the file in the current directory?
  if (-e $filename) {
    return $filename;
  }
  
  # Search for the file.
  foreach (@DIRECTORIES) {
    if (-e $_ . "/" . $filename) {
      return $_ . "/" . $filename;
    }
  }
  
  # Not found.
  return "";
}

# Checks if an ID refers to a file.
#
# Call:
#    $filename = check_id_fbconf($testid)
#
# with:
#    $testid    = a test id to check
#    $filename  = returns the filename of a .fbconf file identified
#                 by $testid or "" if $testid is not a file.
sub check_id_fbconf ($) {
  my ($filename) = @_;
  
  # Is there a file named $filename
  my $fname = find_file($filename);
  
  if ($fname eq "") {
    # No, is there a file .fbconf?
    $fname = find_file($filename . ".fbconf");
  }

  if ($fname eq "") {
    # No, is there a file .fbconf2?
    $fname = find_file($filename . ".fbconf2");
  }
  
  return $fname;
}

# Reads and evaluates an .fbconf file. Each line contains either a test id
# or the name of another .fbconf file; in the latter case, the .fbconf files
# have to be recursively evaluated.
#
# Call:
#    @template = read_and_eval_fbconf($filename)
#
# with
#    $filename = name of the .fbconf file.
#
# Returns
#    @template = List of test ids from the file and its subfiles (if necessary).
sub read_and_eval_fbconf ($) {
  my ($filename) = @_;

  my @idlist;
  
  if ($filename ne "") {
    # Open and read the file -- quick and dirty.
    local(*FHANDLE);
    open(FHANDLE, $filename) 
    or die("$0: ERROR. Could not open file <$filename>!\n");
    @raw_data=<FHANDLE>; 
    close(FHANDLE);
    
    # Check all lines
    foreach (@raw_data) {
      chomp ($_);
      
      # A .fbconf file?
      if ($_ =~ m/.fbconf[2]?$/) {
        # Read recursively.
        @newids = read_and_eval_fbconf ($_);
        push (@idlist,@newids);
      }
      else {
        # Check if this refers to a file.
        my $newfilename = check_id_fbconf ($_);
        if ($newfilename ne "") {
          # Read recursively.
          @newids = read_and_eval_fbconf ($newfilename);
          push (@idlist,@newids);
        }
        else {
          # Push the id.
          push (@idlist,$_);
        }
      }
    }
  }
  
  # return the list.
  return @idlist;
}

# Walks through all directories specified in @DIRECTORIES and parses the
# .fbdef/.fbdef2 files there. Returns a hash containing all test ids in all
# these files.
#
# Call:
#    get_all_testids (\@idfiles, \%idlist)
#
# with
#    @idfiles   = Reference to a list that receives all .fbdef files.
#    %idlist    = List of all id's in all the files.
sub get_all_testids (\@\%) {
  my ($idfiles, $idlist) = @_;
  
  my @files;
  foreach my $dir (@DIRECTORIES) {
    # Find the .fbdef/.fbdef2 files there.
    my @files1 = glob("$dir/*.fbdef");
    my @files2 = glob("$dir/*.fbdef2");
    
    # Add the dir to each filename.
    foreach (@files1) {
      push (@$idfiles,$dir . "/" .$_);
    }
    foreach (@files2) {
      push (@$idfiles,$dir . "/" .$_);
    }
  }
  
  # Read all the files. Find the strings "testid = ..." and save
  # all id's to a hash. Throw a warning if a hash is multiply defined.
  # Save the filename in the hash.
  foreach my $filename (@$idfiles) {
    local(*FHANDLE);
    open(FHANDLE, $filename) 
    or die("$0: ERROR. Could not open file <$filename>!\n");
    
    # Parse the file.
    my $line = <FHANDLE>;
    if ($line =~ m/\s*testid\s*=\s*(\S*)/) {
      if (exists($idlist->{$1})) {
        print "$0: WARNING. Test id <$1> multiply defined.\n";
        print "$0: Original definition in   $idlist->{$1}.\n";
        print "$0: Additional definition in $filename.\n";
      }
      else {
        # Store the id in the hash. Value is the filename.
        $idlist->{$1} = $filename;
      }
    }
    
    close (FHANDLE)    
  }
}

# Generates scripts for a set of test id's and sets up the corresponding
# directories
#
# Call:
#    generate_tests (@testids,@testheader)
#
# With
#    @testids    = list of test id's to generate scripts for.
#    @testheader = template script header for all scripts.
sub generate_tests (\@\@) {
  my ($testids,$testheader) = @_;

  # Get the script filename
  my $basefilename = $ENVIRONMENT{"SCRIPT"};
  
  # Create one or multiple files?
  my $multiplefiles = $ENVIRONMENT{"MULTIPLEFILES"};
  
  # Filename that holds the ID's while caööing create_script.pl.
  my $tempfilename = $ENVIRONMENT{"TEMPIDFILE"};
  
  # Create a shell-invoke command for calling the create_script.pl script.
  # This command sets up an environment with environment variables
  # from %ENVIRONMENT and invoked create_script.pl with this modified
  # environment.
  my $screatescriptinvoke = "";
  foreach (keys(%ENVIRONMENT)) {
    $screatescriptinvoke .= $_ . "=" . $ENVIRONMENT{$_} . " ";
  }
  $screatescriptinvoke .= "./include/create_script.pl";
  
  # Loop through the test id's.
  if ($multiplefiles == 0) {

    # Get a filename
    $filename = $basefilename;
    
    # Open a new file and write the header
    open (FHANDLE,">$filename");

    # Write the header
    foreach (@$testheader) {
      print FHANDLE "$_\n";
    }
    
    # Applications are named according the scheme:
    #   <APPL_BASENAME>-<application>
    # where <application> is both found 
    # * in {kernel,apps,area51}_<application> and 
    # * as value of the keyword 'appl' in a *.fbdef file in subdirectory 'tests'.
    #
    # Define <APPL_BASENAME>
    print FHANDLE "\nAPPL_BASENAME=$APPBASE\n\n";

    # Close the file 
    close (FHANDLE);
    
    # Print all test id's to a temporary .fbconf file.
    # Open a new file and write the header
    open (FHANDLE,">$tempfilename");
    foreach (@$testids) {
      print FHANDLE "$_\n";
    }
    close (FHANDLE);
    
    # Invoke create_script.pl, write all configurations into that file.
    my $output=`$screatescriptinvoke $tempfilename >> $filename`;
    
    unlink ($tempfilename);
    
    # Write a footer.
    open (FHANDLE,">>$filename");
    print FHANDLE "\nfb_footer\n";
    close (FHANDLE);
    
    # Make the script executable
    chmod 0755,$filename;
    
    print "# Script named <$filename> has been created.\n"
   
  }
  else {
    
    # Loop through the test ids
    foreach my $testid (@$testids) {
      
      # Get a unique filename
      $filename = $basefilename . "." . $testid . ".sh";
      
      # Open a new file and write the header
      open (FHANDLE,">$filename");
  
      # Write the header
      foreach (@$testheader) {
        print FHANDLE "$_\n";
      }
      
      # Applications are named according the scheme:
      #   <APPL_BASENAME>-<application>
      # where <application> is both found 
      # * in {kernel,apps,area51}_<application> and 
      # * as value of the keyword 'appl' in a *.fbdef file in subdirectory 'tests'.
      #
      # Define <APPL_BASENAME>
      print FHANDLE "\nAPPL_BASENAME=$APPBASE\n\n";
  
      # Close the file 
      close (FHANDLE);
      
      # Print the test id's to a temporary .fbconf file.
      # Open a new file and write the header
      open (FHANDLE,">$tempfilename");
      print FHANDLE "$testid\n";
      close (FHANDLE);
      
      # Invoke create_script.pl, write all configurations into that file.
      my $output=`$screatescriptinvoke $tempfilename >> $filename`;
      
      unlink ($tempfilename);
      
      # Write a footer.
      open (FHANDLE,">>$filename");
      print FHANDLE "\nfb_footer\n";
      close (FHANDLE);
      
      # Make the script executable
      chmod 0755,$filename;
      
      print "# Script named <$filename> has been created.\n"    
    }
  }
}

################################################################################
# Main program.
#
# Check command line options
if ($#ARGV < 0) {
  die "\nmake_test.pl: Command line script to generate test(s) from .FBDEF files.\n\n";
}

# Get the default environment variables.
get_environment ();

# Classify the command line options.
my @testids;
foreach (@ARGV) {
  if ($_ =~ m/(\S*)=(\S*)/) {
    # This is an option, save in the environment.
    $ENVIRONMENT{$1} = $2;
    next;
  }
  
  if ($_ =~ m/.fbconf[2]?$/) {
    # Filename of an FBCONF file. Read the ID's in the FBCONF file and append
    # them to the testids list.
    my @newids = read_and_eval_fbconf($_);
    push (@testids,@newids);
  }
  else {
    # That is (hopefully) a test id.
    # Check if there is a file named "filename.fbconf". If yes, it's another
    # test suite.
    my $filename = check_id_fbconf ($_);
    if ($filename ne "") {
      # That's a file.
      my @newids = read_and_eval_fbconf($filename);
      push (@testids,@newids);
    }
    else {
      # This is a test id.
      push (@testids,$_);
    }
  }
}

my $testidcount = @testids;
if ($testidcount == 0) {
  die "$0: ERROR. No test id's specified!\n";
}

# Read the template
my @template = read_template_script("include/runtests.template");

# Get all test id's and corresponding test files.
#my @testids;
#my @testfiles;

# Generate the script
generate_tests (@testids,@template);
