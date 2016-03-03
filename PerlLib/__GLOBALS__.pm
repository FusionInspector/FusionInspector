package __GLOBALS__;

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);

our @EXPORT = qw($TRINITY_HOME);

our ($TRINITY_HOME);



BEGIN {
    unless ($TRINITY_HOME = $ENV{TRINITY_HOME}) {
        confess "Error, need TRINITY_HOME env var set to Trinity base installation directory";
    }
}


1;

#EOM

