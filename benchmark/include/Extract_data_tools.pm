#!/usr/bin/env perl
#
# This is a parser library that contains various parser routines to extract
# data from output files. The routines here allow to extract data or blocks
# of data from different files, rearrange them and produce .CSV files from them
# which can be used by 3rd party tools to generate LATEX tables.
#
# Quick start:
#
# For extracting data, the following basic routines are provided:
#
# 1.) grep_datablocks
#  Allows to parse a datafile. The file is searched for a string $anchor.
#  Whereever the $anchor is found, the routine extracts that line plus
#  optionally some lines above/below the line.
#
# 2.) split_to_pairs
#  Splits text data read by grep_datablocks into name=value pairs.
#  All names are stored in a @names array, all values in a @vals array.
#  Both arrays have the same length.
#
# 3.) group_results
#  Creates a hash %resulthash from a (@names,@vals) pair and a list @columns.
#  @columns is a list of column headlines which are searched in @names.
#  All values whose name match the same column name are collected into
#  a list and stored in a hash with the column name defining the key.
#  The hash indexes the values in the form
#
#     $resulthash{$columnname) = \@list-of-values
#
#  Afterwards, one can fetch all results of the column 
#  $columnname by accessing the hash %resulthash.
#
#  Each call to group_results appends more values to the %resulthash hash,
#  so one can collect data of multiple files to the columns. Each call then
#  gives a couple of new rows.
#
# 4.) regroup_results
#  Allows to extract a column from multiple hashes and form a new
#  %resulthash hash. The routine accepts a list @resulthashlist
#  which is a collection of pointers to result-hashes, so
#
#    @resulthashlist = ( \%hash_1, \%hash_2, \%hash_3, ... ).
#
#  The parameter $column is the name of the column to be extracted from
#  %hash_1, %hash_2,... . The parameter $newcolumns defines the new names
#  of the column $column in %hash_1, %hash_2,... . The result is the hash:
#
#     $resulthash{$newcolumns[$i]} = \@ ($hash_$i->{$column})
#
#  Each call to group_results appends more values to the %resulthash hash,
#  so one can collect data of multiple files to the columns. Each call then
#  gives a couple of new rows.
#
# 5.) resulthash_to_list
#  Extracts some or all keys from a resulthash-hash and creates a list of
#  columns with the corresponding values. This defines a table.
#  One specifies a hash %resulthash with the results and an optional
#  list of keys @columns. The routine generates a list containing
#  lists of values according to the columns in the form
#
#    $table[$i] = \@ ($resulthash{$columns[$i]})
#
#  If @columns is not specified, the result is a table containing all
#  keys, in the order of the keys of %resulthash.
#
# 6.) extract_column_from_hashlist
#  Combination of regroup_results and resulthash_to_list. This routine does not
#  create a hash from a list of hashes and a column name, but it creates a list
#  of columns with values. This defines a table.
#  One specifies a column name $column and a list of hashes with
#  results @resulthashlist,
#
#    @resulthashlist = ( \%hash_1, \%hash_2, \%hash_3, ... ).
#
#  and one obtains a list @table with
#
#    $table[$i] = \@ ($hash_$i->{$column}).
#
#  The advantage is that the table does not need unique keys.
#
# 7.) list_to_resultgroups
#  The inverse operation of resulthash_to_list. Takes a list of columns
#  and a list of column names and creates a resulthash hash with the column
#  names as keys and the columns as values.
#
# 8.) list_to_CSV
#  Allows to export a list of columns @table into CSV data.
#  @table contains a list of pointers to other lists with values, so
#  
#    $table[$i] = \@values-of-column-$i
#
#  Such a list can be generated e.g. by extract_column_from_hashlist or
#  resulthash_to_list. The routine produces a list of strings containing
#  CSV output.


package include::Extract_data_tools;

use strict;
use warnings;

# List of the exported subroutines
use Exporter;
our @ISA = ('Exporter');
our @EXPORT = qw(grep_datablocks 
                 filter_datablocks 
                 strip_timestamps
                 split_to_pairs
                 get_all_logfiles
                 read_fbconf
                 group_results
                 regroup_results
                 fillup_lists
                 extract_column_from_hashlist
                 resulthash_to_list
                 list_to_resultgroups
                 merge_lists
                 getline_CSV
                 list_to_CSV
                 empty_CSV_line);

# Directorya of the logfiles.
my $defaultlogdir = "./logs";
my $defaultlog = "runlog";

use constant DEBUG => 0;

# Expands logfile paths.
#
# Call:
#    get_all_logfiles (\@testids,[$filename],[$path])
#
# Returns a list of fully qualified file paths to output log files.
# $filename is the filename of a log file which must be the same
# for all tests; if not specified, $defaultlog is used.
#
# $path is the base directory of all tests; if not present, $defaultlogdir
# is used.
#
# \@testids is a reference to a list of test ids whose log directories
# should be expanded to logfiles.
#
# The routine returns a list with each element having the form
#    "$path/$testid/$filename$.
sub get_all_logfiles (\@;$$) {
  my ($testids,$filename,$path) = @_;
  
  if (!defined ($filename)) {
    $filename = $defaultlog;
  }
  if (!defined ($path)) {
    $path = $defaultlogdir;
  }
  
  my @logfiles = ();
  
  foreach (@$testids) {
    push (@logfiles,$path . "/" . $_ . "/" . $filename);
  }
  
  return @logfiles;
}

# Read an FBCONF file.
#
# Call:
#    read_fbconf ($filename)
#
# Reads the .fbconf file $filename and returns a list of all test ids
# from the file as list.
sub read_fbconf ($) {
  my $filename = shift;
  local (*FILE);
  open (FILE,$filename)
  or die "Cannot open file <$filename>\n";
  
  my @testids = ();
  
  while (<FILE>) {
    chomp;
    push (@testids,$_)
  }
  
  close (FILE);
  
  return @testids;
}

# Greps for blocks of data in a data file.
#
# Call:
#    grep_datablocks ($filename,$anchor,$prevlines,$nextlines)
#
# In:
#    filename      : Name/Path of the file
#    anchor        : Text to grep for.
#    prevlines     : Number of lines to return before the anchor.
#    nextlines     : Number of lines to return after the anchor.
#
# The routine greps for the string [anchor] in the file [filename]
# it returns a list of all lines that are found, plus [prevlines] lines
# before and [nextlines] after the found item.

sub grep_datablocks($$$$) {
  my ($filename,$anchor,$prevlines,$nextlines) = @_;
  
  local (*FILE);
  
  # Open the file
  open (FILE,$filename)
  or die "Cannot open file <$filename>\n";
  
  # Create a list of lines that are read.
  my @outbuffer = ();
  my @tempbuffer = ();
  
  # Read the lines.
  while (<FILE>) {
    chomp;
    print "Line: $_\n" if (DEBUG==1);
    if (@tempbuffer > $prevlines) {
      for (my $i=0; $i < (@tempbuffer-1); ++$i) {
        $tempbuffer[$i] = $tempbuffer[$i+1];
      }
      $tempbuffer[@tempbuffer-1] = $_;
    }
    else {
      $tempbuffer[@tempbuffer] = $_;
    }
    
    # Does the item match?
    if (index($_,$anchor) ne -1) {
      print "Line matched\n" if (DEBUG==1);
      
      # Get up to nextlines more lines.
      for (my $i=0; $i < $nextlines; ++$i) {
        my $line = <FILE>;
        if (!defined($line)) {
          last;
        }
        chomp ($line);
        print "Additional line: $line\n" if (DEBUG==1);
        $tempbuffer[@tempbuffer] = $line;
      }
      
      # Append the temp buffer to the buffer.
      push (@outbuffer,@tempbuffer);
      
      # Create a new temp buffer.
      @tempbuffer = ();
    }
  }
  
  # Close the file.
  close (FILE);
  
  # return empty buffer, not found.
  return @outbuffer;
}

# Filter for blocks of data in a data file.
#
# Call:
#    filter_datablocks (@output,$anchor,$prevlines,$nextlines)
#
# In:
#    output        : text data to be filtered.
#    anchor        : Text to grep for.
#    prevlines     : Number of lines to return before the anchor.
#    nextlines     : Number of lines to return after the anchor.
#
# The routine greps for the string [anchor] in the lines [output] and
# returns a list of all lines that are found, plus [prevlines] lines
# before and [nextlines] after the found item.
#
# prevlines/nextlines are allowed to be negative. In this case, the
# routine filters the corresponding lines after/before the anchor.
#
# The routine works similar to grep_datablocks but filters a list
# instead of a file.

sub filter_datablocks(\@$$$) {
  my ($output,$anchor,$prevlines,$nextlines) = $_;
  
  # Create a list of lines that are read.
  my @outbuffer = ();
  
  # Read the lines.
  my $iline;
  for ($iline = 0; $iline < @$output; ++$iline) {
    my $line = $output->[$iline];
    print "Line: $line\n" if (DEBUG==1);
    
    # Does the item match?
    if (index($line,$anchor) ne -1) {
      print "Line matched\n" if (DEBUG==1);
      
      # Get the lines.
      for (my $i = ($iline-$prevlines >= 0 ? $iline-$prevlines : 0);
              $i <= ($iline+$nextlines < @$output ? $iline+$nextlines : @$output-1);
            ++$i) {
        $outbuffer[@outbuffer] = $output->[$i];
      }
    }
  }
  
  # return empty buffer, not found.
  return @outbuffer;
}

# Remove all timestamps.
#
# Call:
#    strip_timestamps \@output
#
# Loops through all the lines in [output] and removes possible timestamps
# at the beginning of the line.
sub strip_timestamps(\@) {
  my $list = shift;
  foreach (@$list) {
    $_ =~ s/\d\d\.\d\d\.\d\d\d\d\s\d\d\:\d\d\:\d\d\:\s//;
  }
}

# Split a list into pairs
#
# Call:
#    split_to_pairs (\@output,\@name,\@vals)
#
# Splits an output list @output into pairs containing a name and a value.
# Each line of the output is analysed for the pattern
#     "name: value"  or "name=value"
# The name is stored in @names, the value in @vals.

sub split_to_pairs(\@\@\@)
{
  my ($output,$names,$vals) = @_;
  my $icount = 0;
  foreach (@$output) {
    if ($_ =~ m/(.*)[\:=](.*)/) {
      my $name = $1;
      my $val = $2;
      $name =~ s/^\s*//;
      $name =~ s/\s*$//;
      $val =~ s/^\s*//;
      $val =~ s/\s*$//;
      print "Match: $_ = <$name> <$val>\n"
      if (DEBUG==1);
      @$names[$icount] = $name;
      @$vals[$icount] = $val;
      $icount++;
    }
  }
}

# table the results of an output file.
#
# Call:
#    group_results (\@names,\@values,\@columns,\%resulthash)
#
# Analyses the list @names=@values and table results into columns.
# @columns is a list of all items to be searches as name in @names.
# The routine searches for these items and incorporates the values
# into the hash %resulthash. Old values in the hash are not deleted,
# new values are appended.
sub group_results (\@\@\@\%) {
  my ($names,$vals,$columns,$resulthash) = @_;
  
  # Loop through the lines and analyse. Find the columns.
  for (my $iline=0; $iline < @$names; ++$iline) {
    my $name = $names->[$iline];
    my $value = $vals->[$iline];
    if (exists $resulthash->{$name}) {
      # Attach the value.
      my $vallist = $resulthash->{$name};
      push (@$vallist,$value);
    }
    else {
      # Create a new value.
      my @vallist;
      push (@vallist,$value);
      $resulthash->{$name} = \@vallist;
    }
  }
}

# Regroup column data.
#
# Call:
#    regroup_results ($column,@newcolumns,\@resulthashlist,\%resulthash))
#
# Extracts the column $column from all resultgroup-hashes in the list
# @resultgrouplist. Appends the data to the hast %resulthash
# using the keys defined in @newcolumns. For every has in the
# @resulthashlist list, the @newcolumns list must define a new name
# for the column $column when being incorporated into %resulthash.
sub regroup_results ($\@\@) {
  my ($column, $newcolumns, $resulthashlist, $resulthash) = $_;
  
  # Loop through all hashes in the hash-list.
  for (my $ihash=0; $ihash < @$resulthashlist; ++$ihash) {
    # Get the new column name and the corresponding has from the lists.
    my $columnname = $newcolumns->[$ihash];
    my $resultgroup = $resulthashlist->[$ihash];
    
    # Append the data to %resulthash; if the key does not exist, create it.
    if (!exists $resulthash->{$columnname}) {
      my @newlist = ();
      $resulthash->{$columnname} = \@newlist;
    }
    
    my $datalist = $resulthash->{$columnname};
    push (@$datalist,@$resultgroup);
  }
}

# Fill lists to maximum length.
#
# Call:
#   fillup_lists (\%resulthash, [$rowcount], [$defaultval])
#
# Accepts a hash %resulthash which is assumed to be a hash of
# lists. Each list is filled up with the value $defaultval to have
# the same number of rows. If $rowcount is present, the list must have
# at least $rowcount rows. if $defaultval is not present, "" is assumed
# as default value.
sub fillup_lists (\%;$$) {
  my ($resulthash,$rowcount,$defaultval) = @_;
  
  if (!defined ($rowcount)) {
    $rowcount = 0;
  }
  
  # Determine the maximum number of rows in each column.
  foreach my $key (keys %$resulthash) {
    my $vallist = $resulthash->{$key};
    if (@$vallist > $rowcount) {
      $rowcount = @$vallist;
    }
  }

  if (!defined ($defaultval)) {
    $defaultval = "";
  }
  
  # Append $defaultval to the lists if they are not long enough.
  foreach my $key (keys %$resulthash) {
    my $vallist = $resulthash->{$key};
    # Copy the list in memory and fill up.
    my @newlist;
    push (@newlist,@$vallist);
    for (my $i = @$vallist; $i < $rowcount; ++$i) {
      push (@newlist,$defaultval);
    }
    $resulthash->{$key} = \@newlist;
  }
}

# Extract a column from a hash list.
#
# Call:
#    extract_column_from_hashlist ($column, \@resulthashlist, \@table)
#
# @resulthashlist is assumed to be a list of resulthash-hashes, each
# hash corresponding to one file and containing all extracted values
# from a file.
# extract_column_from_hashlist loops now through all hashes in this list
# and extracts the column $column from all the hashes. It appends
# a list of all value-lists of the columns in the order specified
# by @grouplist to the columns in @table.
# The list data is copied in memory.
sub extract_column_from_hashlist ($\@\@) {
  my ($column,$resulthashlist,$table) = @_;
  
  # Loop through all table.
  for (my $ihash=0; $ihash < @$resulthashlist; $ihash++) {
    # Get the hash with the results.
    my $resulthash = $resulthashlist->[$ihash];
  
    # $resulthash is a pointer to a hash with the values.
    # Get the list of the values associated to the column.
    my $valuelist = $resulthash->{$column};
    
    if (!defined($valuelist)) {
      die "Column <$column> undefined!\n" .
      "Existing keys: <" . (keys %$resulthash) . ">\n";
    }
    
    # Append the data to the columns; create columns if necessary.
    if (!defined($table->[$ihash])) {
      my @newlist = ();
      $table->[$ihash] = \@newlist;
    }
    my $targetlist = $table->[$ihash];
    push (@$targetlist,@$valuelist);
  }
}

# Extract columns from resulthash.
#
# Call:
#    resulthash_to_list (\%resulthash,[\@columns])
#
# Extracts some columns from a resulthash list and stores them
# to a list of value-lists. @columns is a list of all column headlines
# that should be extracted. If not present, all table are extracted.
# The list data is copied in memory.
sub resulthash_to_list (\%;\@) {
  my ($resulthash,$columns) = @_;
  
  my @columnkeys;
  
  if (!defined ($columns)) {
    @columnkeys = keys(%$resulthash);
    $columns = \@columnkeys;
  }
  
  my @table = ();
  
  # Loop through all table.
  foreach my $column (@$columns) {
    # $resulthash is a pointer to a hash with the values.
    # Get the list of the values associated to the column.
    my $valuelist = $resulthash->{$column};
    
    # Attach the pointer to the result columns.
    my @newlist;
    push (@newlist,$valuelist);
    $table[@table] = \@newlist;
  }
  return @table;
}

# Convert a list of value-lists to a resulthash hash.
#
# Call:
#    list_to_resultgroups (\@table,\%resulthash,\@columns)
#
# Interprets all lists in @table as value-lists and appends the data
# to the hash %resulthash. @columns defines a list of
# headlines for the columns, which are used as keys for the hash.
# The list data is copied in memory.
sub list_to_resultgroups (\@\%\@) {
  my ($table,$resulthash,$columns) = @_;
  
  # Loop through all data table.
  for (my $igroup = 0; $igroup < @$columns; ++$igroup) {
    my $groupname = $columns->[$igroup];
    my $valuelist = $table->[$igroup];
    
    # Does the group already exist?
    if (exists $resulthash->{$groupname}) {
      # Append the data.
      my $olddata = $resulthash->{$groupname};
      push (@$olddata,@$valuelist);
    }
    else {
      # Copy each list and form a hash from it.
      my @newlist;
      push (@newlist,copy_list($valuelist));
      $resulthash->{$groupname} = \@newlist;
    }
  }
}

# Incorporate the columns of a table into another table.
#
# Call:
#    merge_lists (\@source_table,\@dest_table,$ipos,$iblocksize)
#
# The columns of the table are incorporated into the destination table @table.
# They are inserted in the columns
#    $ipos, $ipos+$iblocksize, $ipos+2*$iblocksize, $ipos+3*$iblocksize, ...
# If the destination table does not contain enough columns, the remaining
# columns of @source_table are appended.
# The return value is the new table.
sub merge_lists (\@\@$$) {
  my ($sourcetable,$desttable,$ipos,$iblocksize) = @_;
  
  # Create a new list.
  my @newlist = ();
  
  # Shuffle the entries.
  my $icol = 0;
  my $relpos = 0;
  foreach my $element (@$desttable) {
    # If $relpos reaches position $ipos, insert the column and decrease
    # $relpos by $blocksize.
    if ($relpos++ == $ipos) {
      push (@newlist,$sourcetable->[$icol++]);
      $relpos -= $iblocksize;
    }
    push (@newlist,$element);
  }
  
  # Append remaining elements
  for (my $i=$icol; $i < @$sourcetable; ++$i) {
    push (@newlist,$sourcetable->[$icol]);
  }
  
  # return the list
  return @newlist;
}

# Generates one line of CSV output.
#
# Call:
#    getline_CSV (\@listofvalues)
#
# The values in the list @listofvalues are concatenated and divided by
# marker characters. The formed line is returned as string.
sub getline_CSV(\@) {
  my ($columns) = @_;
  my $line = "";
  my $colcount = 0;
  foreach (@$columns) {
    $line .= (($colcount++ > 0) ? ";" : "") . "\"" . $_ . "\"";
  }
  return $line;
}

# Convert result-list to CSV data.
#
# Call:
#    list_to_CSV (\@table)
#
# Converts a data list to CSV output. @table is a list of value-lists
# that contains columns of data. The routine generates a list of strings in CSV
# format that can be written to a CSV file.
sub list_to_CSV(\@;\@) {
  my ($table) = @_;
  
  my @lines = ();
  
  # Add the headlines.
  my $line = "";
  my $colcount = 0;
  
  # Determine number of rows from the first column.
  my $col = $table->[0];
  my $rowcount = @$col;
  
  # Loop through the rows.
  for (my $irow=0; $irow < $rowcount; ++$irow) {
    # New line.
    $line = "";

    # Loop through the columns.
    for (my $icol=0; $icol < @$table; ++$icol) {
      # Create the line.
      my $col = $table->[$icol];
      my $data = $col->[$irow];
      $line .= (($icol > 0) ? ";" : "") . "\"" . $data . "\"";
    }
    
    # Save the line.
    $lines[@lines] = $line;
    
    print "CSV: $line\n"
    if (DEBUG==1);
  }
  
  return @lines;
}

# Creates an empty line with as many columns as defined in @columns.
#
# Call:
#    empty_CSV_line (\@columns)
#
# Returns:
#    A string representing an empty CSV line.
sub empty_CSV_line (\@) {
  my $columns = shift;
  
  # Add the data.
  my $line = "";
  my $colcount = 0;
  foreach (@$columns) {
    $line .= (($colcount++ > 0) ? ";" : "") . "\"\"";
  }
  return $line;
}

# default module initialisation; return TRUE to mark the module as
# successfully initialised.
1;
