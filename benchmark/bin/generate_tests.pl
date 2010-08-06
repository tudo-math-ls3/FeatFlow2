#!/usr/bin/env perl
#
# Perl script for the generation of parameter tests.
#
# generate_tests.pl: Command line script to generate .FBDEF files.
#
# Call:
#
#     generate_tests.pl filename [list]
#
# In:
#
#     filename = name of a .fbgen file
#     list     = optional command to generate a list of configurations.
# 
# The utility parses the .fbgen file 'filename' and creates a
# .fbdef compatible configuration set from it which is printed to the
# terminal. If the optional command 'list' is specified, a list of all
# configuration id's in the .fbgen file is printed to the terminal.
#
# The file "demotestset.fbgen" defines a demo file which can be
# parsed with this test generator.

# Allow proper programming only
# (a bit like Fortran90's "implicit none")
use strict;
use warnings;

use Text::ParseWords;

use constant DEBUG => 0;

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
# Left trim function to remove leading whitespace
sub ltrim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($)
{
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}

# Expands a parameter list.
# Call:
#      expand_paramlist(params)
#
# In:
#   count  = expectzed length of the list or =0 if list can have arbitrary length.
#   params = List of parameters to be parsed.
#                
# Returns:
#   Fully expanded list with all parameters.
#
sub expand_paramlist ($$) {
  
  # Get the parameters as a list, trim them.
  my $count = shift;
  my $params = shift;
  my @paramlist = quotewords(",", 0, $params);
  
  my $paramcount = 0;
  
  print STDERR "Analysing list: @paramlist\n"
    if (DEBUG==1);

  # Generate a list from start to end. Note that the sign of $step defines
  # whether we go forward or backward.
  my @elements = ();
    
  foreach my $par (@paramlist) {
    
    # Stop if we reached the number of arguments
    if (($count != 0) && ($paramcount >= $count)) {
      last;
    }

    # Get the trimmed argument and save it.
    $par = trim($par);

    # Detect if there are quotation marks or non-number characters
    # in the string. If yes, this is a list of strings, so we are done.
    if ( $par !~ m/^[\d:\s\.,]*$/ ) {
      push (@elements,$par);
      ++$paramcount;
      next;
    }

    # Check all arguments. If they match the form "xx.xx:yy.yy:zz.zz", we have
    # a sequence which must be evaluated to single numbers.  
    #if ( $par =~ m/\s*([\d\.]*)?\s*:?\s*([\d\.]*)?\s*:?\s*([\d\.]*)?\s*/ ) {
    if ( $par =~ m/\s*([\d\.]*)\s*/ ) {

      # Try to match more to fill $1, $2, $3 as far as possible.
      $par =~ m/\s*([\d\.]*)\s*:\s*([\d\.]*)\s*/;
      $par =~ m/\s*([\d\.]*)\s*:\s*([\d\.]*)\s*:\s*([\d\.]*)\s*/;

      # Is any of the arguments empty?
      my $start;
      my $step;
      my $end;
      
      # Empty step -> step=1
      $step = 1;
      if (defined($3)) {
        if ($3 ne "") {
          $step = $3;
        }
      }
      
      # Empty start -> start is calculated
      $start = "";
      if (defined ($1)) {
        $start = $1;
      }
      
      # Empty end-> end is calculated.
      $end = "";
      if (defined($2)) {
        $end = $2;
      }
      
      if ($start eq "") {
        die "Cannot determine end of sequence!\nParams: $par\n"
          if (($2 eq "") || ($count == 0) );
        $start = $2 - $step*($count-1);
      }
      
      # Empty end -> Calculate
      if ($end eq "") {
        if (!defined ($3)) {
          # No step given -> start=end.
          $end = $start;
        }
        else {
          die "Cannot determine end of sequence!\nParams: $par\n"
            if ($count == 0);
          $end = $start + $step*$count;
        }
      }
      
      print STDERR "Element '$par' evaluated to: $start:$end:$step\n"
        if (DEBUG==1);
      
      if ( $step > 0) {
        while ($start <= ($end + 1E-14)) {
          push (@elements,$start);
          $start = $start + $step;
          ++$paramcount;
          if ( ($count != 0) && ($paramcount >= $count) ) {
            last;
          }
        }
      }
      elsif ( $step < 0) {
        while ($start >= ($end - 1E-14)) {
          push (@elements,$start);
          $start = $start + $step;
          ++$paramcount;
          if ( ($count != 0) && ($paramcount >= $count) ) {
            last;
          }
        }
      }
      else {
        # Stepsize=0. That only works if start=end, otherwise 
        if ( ($count != 0) || ($start == $end) ) {
          push (@elements,$start);
          ++$paramcount;
          if ( ($count != 0) && ($paramcount >= $count) ) {
            last;
          }
        }
        else {
          die "Stepsize of 0 is not allowed!\nParams: $par\n";
        }
      }
    }
  }
  
  # Return the list, finish.
  return @elements; 
}

# Reads an FBGEN file.
# Call:
#    readfbgen (filename)
#
# In: 
#    filename = name/path of the file to read.
#   
# Returns
#    A hash which encapsules the configuration.
sub readfbgenfile ($) {
  
  # Get the filename of the file to read.
  my $filename = $_[0];

  # Open the file  
  local(*FBGENfile);
  open(FBGENfile, $filename)
    or die "\nERROR: Cannot open file <$filename>: $!\n\n";

  print STDERR "File opened: $filename\n"
    if (DEBUG==1);

  # A test is described by the following variables:
  # 1.) Name of the test, application, test class
  my $testname = "";
  my $application = "";
  my $class = "";
  my $descr = "";
  my $comment = "";

  # 2.) A list of include files whcih to include prior to the test definition.
  my @testinclude = ();
  
  # 3.) A list of all variable names.
  my @varnames = ();
  
  # 4.) A hash of all aliases. For each variable in $varnames, one alias
  # is associated. The keys are the variable names.
  my %varalias;
  
  # 5.) A hash of all values for all variables. The keys are the variable names.
  my %varvalues;
  
  # 6.) A hash that contains modifiers vor a variable.
  # Modifiers are: "*": Increase the index of a variable together with it's
  # parent.
  my %varmodifiers;

  # Read the file, line by line.
  my $currentline = "";
  
  LINE: while (defined($currentline = <FBGENfile>)) {
    
    # Remove trailing line break;
    chomp($currentline);
    
    print STDERR "Line: $currentline\n"
      if (DEBUG==1);
      
    # Remove any inlined comments
    $currentline =~ s/\#.+$//;
    
    # Remove trailing white space
    $currentline =~ s/\s+$//;
    
    # Ignore comments (i.e. lines starting with a '#')
    next LINE if ($currentline =~ m/^\s*\#/);
    
    # Ignore empty lines
    next LINE if ($currentline =~ m/^\s*$/);

    # Does the line end with a "\" character? If yes, we have to append
    # the next line to it.
    while ($currentline =~ m/.*\\$/) {
      
      # Remove the last character, probably the last two.
      chop($currentline);

      # If there is double-backslash remove the newline.
      if ($currentline =~ m/.*\\$/) {
        chop($currentline);
        $currentline .= "\n";
      }
      
      # Read the next line
      if (defined (my $nextline = <FBGENfile>)) {
        chomp ($nextline);
        $currentline .= $nextline;
      }
    }
      
    # If the test name is not defined yet, match the test name.
    if ($testname eq "") {
      if ($currentline =~ m/^\s*testname\s*=\s*(.+)\s*$/) {
        $testname = $1;
        print STDERR "Name of the test: $testname\n"
          if (DEBUG==1);
        next LINE;
      }
    }

    # Application name
    if ($application eq "") {
      if ($currentline =~ m/^\s*appl\s*=\s*(.+)\s*$/) {
        $application = $1;
        print STDERR "Name of the application: $application\n"
          if (DEBUG==1);
        next LINE;
      }
    }

    # Test class
    if ($class eq "") {
      if ($currentline =~ m/^\s*class\s*=\s*(.+)\s*$/) {
        $class = $1;
        print STDERR "Name of the test: $class\n"
          if (DEBUG==1);
        next LINE;
      }
    }

    # Test description
    if ($descr eq "") {
      if ($currentline =~ m/^\s*descr\s*=\s*(.+)\s*$/) {
        $descr = $1;
        print STDERR "Description of the test: $descr\n"
          if (DEBUG==1);
        next LINE;
      }
    }

    # Test comment
    if ($comment eq "") {
      if ($currentline =~ m/^\s*comment\s*=(.*)\s*/s) {
        $comment = $1;
        print STDERR "Comment: $comment\n"
          if (DEBUG==1);
        next LINE;
      }
    }

    # Append include files.
    if ($currentline =~ m/^\s*testinclude\s*=\s*(.+)\s*$/) {
      push (@testinclude,$1);
      print STDERR "File to include: $1\n"
        if (DEBUG==1);
      next LINE;
    }
    
    # Append variable names.
    if ($currentline =~ m/^\s*(\w+)([\*\.]?)\s*:?\s*(\w*)\s*=\s*(.+)\s*$/) {
      
      # Get the variable name and a possible modifier.
      my $varname = $1;
      my $varmod = $2;
      
      # Get the alias if defined. If not, the name is the alias.
      my $alias = $3;
      if ($alias eq "") {
        # The alias is the name itself -- without a possibly trailing "*".
        $alias = $varname;
      }
      
      my $varvalue = $4;

      print STDERR "Variable found: $varname$varmod:$alias: '$varvalue'\n"
        if (DEBUG==1);

      # Push the variable name.
      push (@varnames,$varname);

      # Remember the alias.
      $varalias{$1} = $alias;
      
      # Put the values into the hash.
      $varvalues{$varname} = $varvalue;

      # Put the modifiers into the hash.
      $varmodifiers{$varname} = $varmod;
            
      next LINE;
    }
      
  }

  close (FBGENfile);
  print STDERR "File closed: $filename\n"
    if (DEBUG==1);

  # Generate the return value.
  my %configuration =
    ( "testname" => $testname,
      "appl" => $application,
      "class" => $class,
      "descr" => $descr,
      "comment" => $comment,
      "testinclude" => \@testinclude,
      "varnames" => \@varnames,
      "varalias" => \%varalias,
      "varmodifiers" => \%varmodifiers,
      "varvalues" => \%varvalues );
  
  return %configuration;
}

# Expands all values in the configuration.
# Call:
#    expand_values(\%configuration)
#
# In:
#    configuration = a configuration hash.
#
# Changes the values in the configuration from a string to an (expanded) list
# of values.
sub expand_values($) {
  
  my $configuration = shift;
  
  # Loop through all variable names
  my $varnames = $configuration->{"varnames"};
  
  my $lastcount = 0;
  
  foreach my $var (@$varnames) {
    
    # Get the line;
    my $unexpanded_line = $configuration->{"varvalues"}->{$var};
    my @expanded_line = ();
    print STDERR "Expanding <$var> = $unexpanded_line\n"
      if (DEBUG==1);

    # If the variable name has a trailing "*", it is a dependent variable
    # and has as many entries as its parent -- the last non-"*" parameter.
    if ($var =~ m/.*\*$/) {
      # Dependent variable.
      @expanded_line = expand_paramlist($lastcount,$unexpanded_line);
    } else {
      # Standard variable. Expand the line and remember the number of
      # entries for a possible 'dependend' variable.
     
      @expanded_line = expand_paramlist(0,$unexpanded_line);
      $lastcount = @expanded_line;
    }
    print STDERR "Expanded line: @expanded_line\n"
      if (DEBUG==1);
    

    # Replace the line with a reference to the expanded line.
    $configuration->{"varvalues"}->{$var} = \@expanded_line;
    
  }
}

# Generates all test-id's from a fully qualified configuration
# Call:
#    generate_tests(\%configuration)
#
# In:
#    onlyids       = Whether to generate only an ID-list or
#       the whole test configuration. =1: Generate only list of ID's.
#    configuration = a configuration hash.
#
# The configuration is printed to STDOUT.
sub generate_tests($$) {
  
  # Get the parameters
  my $onlyids = shift;
  my $configuration = shift;
  
  # Get a sorted list of all variables etc.
  my $varnames = $configuration->{"varnames"};
  my @varnamessorted = sort @$varnames;
  print STDERR "Variables to process: @$varnames\n"
    if (DEBUG==1);

  my $varalias = $configuration->{"varalias"};
  my $varvalues = $configuration->{"varvalues"};
  my $varmod = $configuration->{"varmodifiers"};
  my $testinclude = $configuration->{"testinclude"};
  
  # Create a hash for indexing of the current configuration.
  my %varindex;
  foreach (@$varnames) {
    $varindex{$_} = 0;
  }
  
  # Get the application name and test class.
  my $testname = $configuration->{"testname"};
  my $application = $configuration->{"appl"};
  my $class = $configuration->{"class"};
  my $descr = $configuration->{"descr"};
  my $comment = $configuration->{"comment"};
  
  # Print a header
  if ($onlyids == 0) {
    print "\n";
    print "# ============================================================\n";
    print "# $testname\n";
    print "# ============================================================\n";
    print "\n";
    if ($comment ne "") {
      # Print the comment line by line
      my @lines = split("\n",$comment);
      foreach(@lines) {
        chomp;
        print "#$_\n";
      }
    }
    print "\n";
    print "appl = $application\n" if ($application ne "");
    print "class = $class\n" if ($class ne "");
    print "\n";
  }
  
  # Loop until the index of the 1st element gets out-of-bounds.
  # If that's the case, we are done. 
  LOOP: while (1==1) {
   
    # Get the actual test id.
    my $testid = $testname;
    foreach my $varname (@varnamessorted) {
      # Get the corresponding modifier
      my $varmod = $varmod->{$varname};
      
      # Append the current index and variable name to the test id.
      # Variables with a "." at the end are not included.
      if ($varmod !~ m/\./) {
        $testid = $testid . "-" . $varname . ($varindex{$varname}+1);
      }
    }
    
    # Remove any "*" or ".".
  $testid =~ s/[\*\.]//;
    
    print STDERR "Generated test-id: $testid\n"
      if (DEBUG==1);
    
    # Print the test.
    if ($onlyids == 0) {
      print "\n";
      print "testid = $testid\n";
      print "descr = $descr\n";
      
      # Print the includes.
      foreach my $incl (@$testinclude) {
        print "include $incl\n";
      }
      
      # Print all parameters (or so to say, their alias values)
      # with the corresponding value.
      foreach my $varname (@$varnames) {
        
        # Get the alias
        my $alias = $varalias->{$varname};
        
        # Get the value.
        my $idx = $varindex{$varname};
        my $parlist = $varvalues->{$varname};
        my $value = $$parlist[$idx];
        
        # Check if there is the term "$(...)" in the value.
        # If yes, try to replace it with the current value of
        # the corresponding varible.
        while ( $value =~ m/\$\((\%?)([^\)]+)[\*\.]?\)/ ) {
          # Try to get the referring value.
          my $varname2 = $2;
          my $idx2 = $varindex{$varname2};
          die "ERROR: Could not determine value of referring ".
              "variable <$varname2>\nOriginal value: $value\n"
            if ( !defined($idx2) );

          print STDERR "Detected replacement variable '$varname2' in '$value'."
            if (DEBUG==1);

          my $parlist2;
          my $value2;
          
          if ($1 eq "%") {
            $value2 = $idx2;
          }
          else {
            $parlist2 = $varvalues->{"NLMAX"};
            $value2 = $$parlist2[$idx2];
          }
          print STDERR " Replacing '\$($1$varname2)' by '$value2'.\n"
            if (DEBUG==1);
          $value =~ s/\$\(\%?([^\)]+)[\*\.]?\)/$value2/;
        }
        
        # Probably invoke the shell to execute any subcommand.
        if ( $value =~ m/`([^`]+)`/ ) {
          my $value2 = `$1`;
          chomp($value2);
          
          print STDERR "Replacing command string '$1' by '$value2'.\n"
            if (DEBUG==1);
          $value =~ s/`([^`]+)`/$value2/;
        }
        
        print "$alias = $value\n";
        
      }
      
    }
    else {
      print "$testid\n";
    }
    
    # Increase from the rear to the front
    for (my $idx=@$varnames-1; $idx >= 0; --$idx) {
      
      # Corresponding variable name...
      my $varname = $$varnames[$idx];
      
      # Corresponding modifiers...
      my $varmod = $varmod->{$varname};
      
      # Get the corresponding parameter list of parameter $idx.
      my $parlist = $varvalues->{$varname};
      
      # Increase the counter from the rear. If the number gets
      # too large, put it back to 1 and increase the
      # previous position. Otherwise, we can continue
      # the loop directly, as the previous positions
      # stay unchanged.
      #
      # If the variable modifier contains a '.', the variable
      # is not increased.
      #
      # RULE DISABLED: The . below is more general.
      # if ( $varmod =~ m/\./ ) {
      #   next;
      # }
      
      if (++($varindex{$varname}) >= @$parlist) {
        $varindex{$varname} = 0;
      }
      else {
        # We had no overflow.
        # If this is a 'master' parameter, we can stop here, as the previous
        # parameters do not change.
        # If this is a 'child' parameter, also increase the previous
        # parameter until we increased the 'master' of the group.
        if ( $varmod !~ m/[\*.]/ ) {
          next LOOP;
        }
      }
    }
    
    # FOR-Loop finished, i.e. the index at position 0 had an overflow.
    # That means, all configurations are done.
    last;
  }
  
}

# Main program.
#
# Check command line options
if ($#ARGV < 0) {
  die "\ngenerate_tests.pl: Command line script to generate .FBDEF files.\n\n" .
      "Call:\n\n" .
      "    generate_tests.pl filename [list]\n\n" .
      "In:\n\n" .
      "    filename = name of a .fbgen file\n" .
      "    list     = optional command to generate a list of configurations.\n" .
      "\n" .
      "The utility parses the .fbgen file 'filename' and creates a\n" .
      ".fbdef compatible configuration set from it which is printed to the\n" .
      "terminal. If the optional command 'list' is specified, a list of all\n" .
      "configuration id's in the .fbgen file is printed to the terminal.\n\n";
}

my %configuration = readfbgenfile($ARGV[0]);
expand_values(\%configuration);

# Get the 2nd parameter.
my $command = "";
if ($#ARGV >= 1) {
  $command = $ARGV[1];
}

# Generate the tests.
generate_tests($command eq "list",\%configuration);
