#!/usr/bin/env perl
#
# Nonblockingly read a single key from STDIN,
# adapted from http://search.cpan.org/dist/TermReadKey/ReadKey.pm
#
# Advantages:
# * returns the key that got actually pressed (removing it from the buffer
#   such that subsequent calls of this script do not always return true)
#
# Disdvantages:
# * requires non-default Perl module Term::ReadKey,
#   hence the script might not run
#   (i.e. non-default package libterm-readkey-perl on Ubuntu 10.04.3,
#         not installed on SOLARIS Sparc)

use warnings;
use strict;

use Term::ReadKey;
ReadMode 4; # Turn off controls keys
my $key = "";
$key = ReadKey(-1);
ReadMode 0; # Reset tty mode before exiting

exit ($key ? 1 : 0);
