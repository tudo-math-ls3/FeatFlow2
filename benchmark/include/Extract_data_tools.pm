#!/usr/bin/env perl
#
# This is a parser library that contains various parser routines to extract
# data from output files. The routines here allow to extract data or blocks
# of data from different files, rearrange them and produce .CSV files from them
# which can be used by 3rd party tools to generate LATEX tables.
#

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
                 getline_CSV
                 list_to_CSV
                 empty_CSV_line
                 getline_GNUPLOT_comment
                 list_to_GNUPLOT
                 merge_table_columns
                 merge_table_entries
                 table_push_row
                 table_unshift_row
                 create_empty_table
                 direct_grep_to_table
                 format_column_width
                 format_columns_width
                 insert_column
                 delete_column
                 delete_empty_columns
                 column_to_table
                 print_table);

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

# Greps for blocks of data in a data file.
# Returns a list with the blocks found.
#
# Call:
#    grep_separate_datablocks ($filename,$anchor,$prevlines,$nextlines)
#
# In:
#    filename      : Name/Path of the file
#    anchor        : Text to grep for.
#    prevlines     : Number of lines to return before the anchor.
#    nextlines     : Number of lines to return after the anchor.
#
# The routine greps for the string [anchor] in the file [filename]
# it returns a list of columns with all lines that are found, plus 
# [prevlines] lines before and [nextlines] after the found item.
#
# grep_datablocks returns just a list of all lines found.
# grep_separate_datablocks returns a list of columns, each column a list
# of strings fond in one block. That way, one can figure out
# which strings belong to which block.
sub grep_separate_datablocks($$$$) {
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
      
      # Append the temp buffer to the buffer as new column.
      my @col = ();
      push (@col,@tempbuffer);
      push (@outbuffer,\@col);
      
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
#    regroup_results ($col,@newcolumns,\@resulthashlist,\%resulthash))
#
# Extracts the column $col from all resultgroup-hashes in the list
# @resultgrouplist. Appends the data to the hast %resulthash
# using the keys defined in @newcolumns. For every has in the
# @resulthashlist list, the @newcolumns list must define a new name
# for the column $col when being incorporated into %resulthash.
sub regroup_results ($\@\@) {
  my ($col, $newcolumns, $resulthashlist, $resulthash) = $_;
  
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
#    extract_column_from_hashlist ($col, \@resulthashlist, \@table)
#
# @resulthashlist is assumed to be a list of resulthash-hashes, each
# hash corresponding to one file and containing all extracted values
# from a file.
# extract_column_from_hashlist loops now through all hashes in this list
# and extracts the column $col from all the hashes. It appends
# a list of all value-lists of the columns in the order specified
# by @grouplist to the columns in @table.
# The list data is copied in memory.
sub extract_column_from_hashlist ($\@\@) {
  my ($col,$resulthashlist,$table) = @_;
  
  # Loop through all table.
  for (my $ihash=0; $ihash < @$resulthashlist; $ihash++) {
    # Get the hash with the results.
    my $resulthash = $resulthashlist->[$ihash];
  
    # $resulthash is a pointer to a hash with the values.
    # Get the list of the values associated to the column.
    my $valuelist = $resulthash->{$col};
    
    if (!defined($valuelist)) {
      die "Column <$col> undefined!\n" .
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
#    @datatable = resulthash_to_list (\%resulthash,[\@columns])
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
  foreach my $col (@$columns) {
    # $resulthash is a pointer to a hash with the values.
    # Get the list of the values associated to the column.
    my $valuelist = $resulthash->{$col};
    
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


# Creates an empty table with $colnum columns.
#
# Call:
#    create_empty_table ($colnum,@datatable)
#
# In:
#    colnum    : number of columns.
#    datatable : reference to a list that receives the table.
sub create_empty_table ($\@) {
  
  my ($colnum,$datatable) = @_;
  
  # Each column is a list, so just create a list with $colnum sublists.
  # Column 0 is the row number.
  @$datatable = ();
  for (my $i=0; $i < $colnum; ++$i) {
    my @col = ( );
    push (@$datatable,\@col);
  }
}

# Appends a row to a table.
#
# Call:
#    table_push_row (@row,@datatable)
#
# In:
#    row          : A list of values, one value for each column.
#    datatable    : A list of columns representing the table.
sub table_push_row (\@\@) {
  
  my ($row,$datatable) = @_;
  
  # Create a new row. First column is the row number.
  # The remaining columns are filled with $emptystring.
  my $firstcol = $datatable->[0];
  my $rowcount = @$firstcol;
  my @myrow = ( );
  
  # Copy the values
  push (@myrow,@$row);
  
  # Fill with "" if too short.
  my $columncount = @$datatable;
  while (@myrow < $columncount) {
    push (@myrow,"")
  }
  
  my $size=@$datatable;
  print "table_push_row: Size=$size, Row: <@myrow>\n" if (DEBUG == 1);
  
  # Append the column.
  for (my $i=0; $i < @$datatable; ++$i) {
    my $col = $datatable->[$i];
    print "table_push_row: Appending <$myrow[$i]>\n" if (DEBUG == 1);
    push (@$col,$myrow[$i]);
    print "table_push_row: Col $i: @$col\n" if (DEBUG == 1);
  }
      
}

# Inserts a row to the beginning of a table.
#
# Call:
#    table_unshift_row (@row,@datatable)
#
# In:
#    row          : A list of values, one value for each column.
#    datatable    : A list of columns representing the table. The first column
#                   is the row number (maintained automatically).
sub table_unshift_row (\@\@) {
  
  my ($row,$datatable) = @_;
  
  # Create a new row. First column is the row number.
  # The remaining columns are filled with $emptystring.
  my $firstcol = $datatable->[0];
  my $rowcount = @$firstcol;
  my @myrow = ( );
  
  # Copy the values
  push (@myrow,@$row);
  
  # Fill with "" if too short.
  my $columncount = @$datatable;
  while (@myrow < $columncount) {
    push (@myrow,"")
  }
  
  # Append the column.
  for (my $i=0; $i < @$datatable; ++$i) {
    my $col = $datatable->[$i];
    unshift (@$col,$myrow[$i]);
  }
}

# Inserts a column @column to row $col in table $datatable
#
# Call:
#    insert_column (@column,$col,@datatable)
#
#    column    : List of entries; column to insert.
#    col       : Position in the datatable where to insert.
#    datatable : The table to modify.
sub insert_column (\@$\@) {
  my ( $column, $col, $datatable ) = @_;
  
  splice (@$datatable,$col,0,$column);
}

# Deletes column $col from a table @datatable.
#
# Call:
#    delete_column ($col,@datatable)
#
#    col       : Column number to delete.
#    datatable : The table to modify.
sub delete_column ($\@) {
  my ( $col, $datatable ) = @_;
  my $size2 = @$datatable;
  
  if ($col > @$datatable-1) {
    my $size = @$datatable;
    print "delete_column: Warning: Try to delete column $col of $size.\n";
    return;
  }
  
  splice (@$datatable,$col,1);
}

# Splits each entry on the whitespaces and creates a table.
#
# Call:
#    @datatable = column_to_table
sub column_to_table (\@) {
  my ( $column ) = @_;
  
  # 1st loop through the new columns: determine table width.
  my $maxsize = 0;
  foreach (@$column) {
    my @row = split (/\s+/,$_);
    if (@row > $maxsize) {
      $maxsize = @row;
    }
  }
  
  # 2nd loop: Add the rows.
  my @datatable;
  create_empty_table ($maxsize,@datatable);
  foreach (@$column) {
    my @row = split (/\s+/,$_);
    table_push_row (@row,@datatable);
  }
  return @datatable;
}

# Directly grep results from a file and put them into a table.
#
# Call:
#    direct_grep_to_table ($filename,$anchor,$prevlines,$nextlines,@columnids,
#                          $minrows,@datatable)
#
# In:
#    filename      : Name of the datafile.
#    anchor        : Text to grep for.
#    prevlines     : Number of lines to return before the anchor.
#    nextlines     : Number of lines to return after the anchor.
#    columnids     : List of all ids that should be captured ("id=...")
#    minrows       : Minimum number of rows that must be found.
#                    If less rows are found, the remaining rows
#                    are willed with $emptystring.
#    datatable     : Data table. A list of columns. Each column is a list
#                    of strings. The first column contains a row id which is
#                    automatically increased by the number of rows found.
#                    The remaining rows contain the data in the order specified
#                    by @columnids.
#
# Greps for the string [anchor] in the file [filename]
# plus [prevlines] lines before and [nextlines] after the anchor.
# Searches for all id in @columnids and puts their data to the table @table.
sub direct_grep_to_table ($$$$\@$\@) {

  my ($filename,$anchor,$prevlines,$nextlines,$columnids,
      $minrows,$datatable) = @_;

  # Grep for all datablocks in the file.
  my @datalines = ();
  if (-e $filename) {
    @datalines = grep_separate_datablocks ($filename,$anchor,$prevlines,$nextlines);
  }
  else
  {
    print "$0: Warning. File not found: <$filename>\n";
  }
  
  # Create a hash that contains for every column id the column number.
  my %colnum;
  for (my $i = 0; $i < @$columnids; ++$i) {
    
    # Cut away any spaces and a possible ending character.
    my $name = $columnids->[$i];
    $name =~ s/^\s*//;
    $name =~ s/\s*[\:=]?$//;

    $colnum{$name} = $i;
  }

  # Size of the table, number of columns.
  my $firstcol = $datatable->[0];
  
  if (!defined($firstcol)) {
    print "direct_grep_to_table: Warning. Table empty.\n";
    return;
  }
  
  my $initrowcount = @$firstcol;
  my $rowcount = @$firstcol;
  
  # Loop through all entries found.
  foreach my $block (@datalines) {
    
    # Strip timestamps and divide into pairs: name/value
    strip_timestamps (@$block);
    
    my @names = ();
    my @vals = ();
    split_to_pairs (@$block,@names,@vals);
    print "Names: @names\nValues: @vals\n"
    if (DEBUG==1);
    
    # Create a new row, filled with "".
    my @row = ( );
    foreach (@$columnids) {
      push (@row,"");
    }
    
    # Fill the columns.
    for (my $i = 0; $i < @names; ++$i) {
      # Get the column id from the hash -- if it is a column.
      my $name = $names[$i];
      my $id = $colnum{$name};
      if (defined($id)) {
        print "$name = $id = <$vals[$i]>\n"
        if (DEBUG == 1);
        $row[$id] = $vals[$i];
      }
    }
    
    # Add the row to the table.
    table_push_row(@row,@$datatable);
    ++$rowcount;
  }
  
  while ($rowcount < $initrowcount + $minrows) {
    # Append the missing rows as empty rows.
    my @row;
    foreach (@$columnids) {
      push (@row,"");
    }
    table_push_row(@row,@$datatable);
    ++$rowcount;
  }

}

# Formats a column such that all entries have the same width.
# Fills up with spaces.
#
# Call:
#     format_column_width (@column)
#
# In:
#    column    : list of entries.
sub format_column_width (\@) {
  
  my ( $column ) = @_;
  
  # Determine max. length.
  my $len = 0;
  
  foreach (@$column) {
    my $slen = length($_);
    if ($slen > $len) {
      $len = $slen;
    }
  }
  
  # Fill the columns, right adjusted.
  foreach (@$column) {
    my $slen = length($_);
    $_ = " "x($len-$slen) . $_;
  }
}

# Formats all column such that all entries in a column have the same width.
# Fills up with spaces.
#
# Call:
#     format_column_width (@table)
#
# In:
#    column    : list of entries.
sub format_columns_width (\@) {
  
  my ( $datatable ) = @_;
  
  # Format the columns.
  foreach (@$datatable) {
    format_column_width (@$_);
  }
}

# Deletes all empty columns in a table.
#
# Call:
#    delete_empty_columns (@datatable,$headcount)
#
# In:
#    datatable : Table with data.
#    headcount : OPTIONAL: Number of lines in the header which do not contain 
#                data. These lines are ignored upon checking for empty columns.
#                =0: check all rows.
sub delete_empty_columns (\@;$) {
  my ( $datatable, $headcount ) = @_;
  if (!defined $headcount) {
    $headcount = 0;
  }

  # Check all columns.
  for (my $colid=@$datatable-1; $colid >= 0; --$colid) {
    my $column = $datatable->[$colid];
  
    # Check if the column is empty.
    my $empty = 1;
    for (my $rowid=$headcount; $rowid < @$column; ++$rowid) {
      if ($column->[$rowid] ne "") {
        $empty = 0;
        last;
      }
    }
    
    if ($empty == 1) {
      print "Column $colid empty.\n"
      if (DEBUG==1);

      # Delete the column
      delete_column ($colid,@$datatable);
    }
  }
  my $size = @$datatable;
  print "Remaining columns: $size\n"
  if (DEBUG==1);
}

# Incorporate the columns of a table into another table.
#
# Call:
#    @datatable = merge_table_columns (\@source_table,\@dest_table,$ipos,$iblocksize)
#
# The columns of the table are incorporated into the destination table @table.
# They are inserted in the columns
#    $ipos, $ipos+$iblocksize, $ipos+2*$iblocksize, $ipos+3*$iblocksize, ...
# If the destination table does not contain enough columns, the remaining
# columns of @source_table are appended.
# If $iblocksize=0, @source_table is completely inserted at position $ipos.
# The return value is the new table.
sub merge_table_columns (\@\@$$) {
  my ($sourcetable,$desttable,$ipos,$iblocksize) = @_;
  
  # Create a new list.
  my @newlist = ();
  
  # Shuffle the entries.
  my $icol = 0;
  my $relpos = 0;
  foreach my $element (@$desttable) {
    # If $relpos reaches position $ipos, insert the column and decrease
    # $relpos by $blocksize.
    while (($relpos == $ipos) && ($icol < @$sourcetable)) {
      push (@newlist,$sourcetable->[$icol++]);
      $relpos -= $iblocksize;
    }
    push (@newlist,$element);
    ++$relpos;
  }
  
  # Append remaining elements
  for (my $i=$icol; $i < @$sourcetable; ++$i) {
    push (@newlist,$sourcetable->[$i]);
  }
  
  # return the list
  return @newlist;
}

# Incorporate a table in another one.
#
# Call:
#    @datatable = merge_table_entries (\@source_table,\@dest_table)
#
# In:
#   @source_table   : table which is to be incorporated to @dest_table
#   @dest_table     : destination table
#
# Whereever there are nonempty entries (!= "") in @source_table, the routine
# overwrites the entries in @dest_table with @source_table. Old values in
# @dest_table are overwritten. If the entry in @source_table is empty (""),
# the value in @dest_table is kept.
sub merge_table_entries (\@\@) {
  my ($sourcetable,$desttable) = @_;
  
  # Create a new list.
  my @newlist = ();
  
  # Loop through all the columns.
  my $columns = @$sourcetable;
  if (@$desttable > $columns) {
    $columns = @$desttable;
  }
  
  for (my $icol = 0; $icol < $columns; ++$icol) {
    # Get the column in the source and destination table.
    my $sourcecolumn = $sourcetable->[$icol];
    my $destcolumn = $desttable->[$icol];
    
    # A column may be undefined; in that case, create a new one.
    if (!defined ($sourcecolumn)) {
      my $rowcount = @$destcolumn;
      my @sourcecol = ("") x $rowcount;
      $sourcecolumn = \@sourcecol;
    }

    if (!defined ($destcolumn)) {
      my $rowcount = @$destcolumn;
      my @destcol = ("") x $rowcount;
      $destcolumn = \@destcol;
    }
   
    # Number of rows
    my $rowcount = @$sourcecolumn;
    
    # Loop through the rows, create a new column.
    my @newcolumn = ();
    for (my $i = 0; $i < $rowcount; ++$i) {
      if ($sourcecolumn->[$i] ne "") {
        push (@newcolumn,$sourcecolumn->[$i])
      }
      else {
        push (@newcolumn,$destcolumn->[$i])
      }
    }
    
    # Save the column.
    push (@newlist,\@newcolumn)
  }
  
  # Return the new table.
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

# Generates one line of GNUPLOT comment output.
#
# Call:
#    getline_GNUPLOT (\@listofvalues)
#
# The values in the list @listofvalues are concatenated and divided by
# marker characters. The formed line is returned as string.
sub getline_GNUPLOT_comment(\@) {
  my ($columns) = @_;
  my $line = "";
  my $colcount = 0;
  foreach (@$columns) {
    $line .= (($colcount++ > 0) ? "  " : "# ") . $_;
  }
  return $line;
}

# Convert result-list to GNUPLOT data.
#
# Call:
#    list_to_GNUPLOT (\@table)
#
# Converts a data list to GNUPLOT output. @table is a list of value-lists
# that contains columns of data. The routine generates a list of strings in GNUPLOT
# format that can be written to a GNUPLOT file.
sub list_to_GNUPLOT(\@) {
  my ($datatable) = @_;
  
  my @lines = ();
  
  # Add the headlines.
  my $line = "";
  my $colcount = @$datatable;
  if ($colcount <= 0) {
    # Table empty.
    print "list_to_GNUPLOT: Warning. Table empty.\n";
    return @lines;
  }  
  
  # Determine number of rows from the first column.
  my $col = $datatable->[0];
  my $rowcount = @$col;
  
  # Loop through the rows.
  for (my $irow=0; $irow < $rowcount; ++$irow) {
    # New line.
    $line = "";

    # Loop through the columns.
    for (my $icol=0; $icol < @$datatable; ++$icol) {
      # Create the line.
      my $col = $datatable->[$icol];
      my $data = $col->[$irow];
      $line .= "  " . $data;
    }
    
    # Save the line.
    $lines[@lines] = $line;
    
    print "GNUPLOT: $line\n"
    if (DEBUG==1);
  }
  
  return @lines;
}

# Prints a table to screen.
#
# Call:
#    print_table (@datatable)
sub print_table (\@) {
  my ($datatable) = @_;

  my @lines = list_to_GNUPLOT(@$datatable);
  foreach (@lines) {
    print "$_\n";
  }
}

# default module initialisation; return TRUE to mark the module as
# successfully initialised.
1;
