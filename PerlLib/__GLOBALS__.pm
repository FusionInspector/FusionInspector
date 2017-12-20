package __GLOBALS__;

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);

our @EXPORT = qw($TRINITY_HOME);

our ($TRINITY_HOME);



BEGIN {
# Removed by Cicada Dennis on 2017/12/20.
# Remove this because TRINITY_HOME is not always needed.
# it is only needed when the Trinity option is on the command line.
#    unless ($TRINITY_HOME = $ENV{TRINITY_HOME}) {
#        confess "Error, need TRINITY_HOME env var set to Trinity base installation directory";
#    }
}


1;

#EOM

