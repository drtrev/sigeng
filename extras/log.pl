#!/usr/bin/perl -w
# Format logfile
# NOTE: assumes log file ends in .txt
# TODO: only-sig option doesn't work when ANOVA is printed
# TODO: use pdflatex!

use strict;
use warnings;
use File::Basename;

use Getopt::Long;

sub fixhere {
  my $str = shift;
  $str =~ s/^[^\S\n]+//gm;
  return $str;
}

sub help {
  # ridiculous thing so that here document can be indented!
  print fixhere <<"    END";
    Usage: log.pl [--only-sig] <logfile>
    (Will run pandoc and produce tex and pdf with approx same basename)
    May clobber <logfile>.*
    END
}

if ($#ARGV < 0) {
  print "Error, no log file specified.\n";
  exit;
}

my $logfile;
my $help=0;
my $nolatex=0;
my $onlysig=0;

GetOptions('help' => \$help, 'nolatex' => \$nolatex, 'onlysig' => \$onlysig);

if ($help) {
  help;
  exit;
}

$logfile = shift @ARGV;
#print "Logfile: $logfile\n";
#print "Nolatex: $nolatex\n";
#print "Onlysig: $onlysig\n";
#exit;

if (!-e $logfile) {
  print "Error, log file not found.\n";
  exit;
}

sub removenonsig {
  my $filebase = shift @_;
  my $filepath = shift @_;
  my $filesuff = shift @_;

  my $filename = "$filepath$filebase$filesuff";
  my $outname = "$filepath$filebase-onlysig.txt";

  open(FILE, "<$filename") || die("log.pl: removenonsig: Could not open file \"$filename\" - $!");

  my @store;

  while(<FILE>) {
    push @store, $_;
    
    if ($_ =~ /Found no trends/) {
      if ($store[$#store-1] =~ /Found no sig/) {
        # we've got no trends and no sig, so remove last 6 lines
        for my $i (1..6) { pop @store; }
      }
    }
  }

  close(FILE);

  open(FILE, ">$outname") || die("Could not open file \"$outname\" - $!");
  for my $i (@store) {
    print FILE $i;
  }
  close(FILE);


  return $outname;
}

my ($logname, $logpath, $logsuffix) = fileparse($logfile, "\\.txt");

if ($onlysig) {
  $logfile = removenonsig($logname, $logpath, $logsuffix);
  if (!-e $logfile) {
    print "Error, new log file not found.\n";
    exit;
  }
  ($logname, $logpath, $logsuffix) = fileparse($logfile, "\\.txt");
}

print "$logname, $logpath, $logsuffix\n";

my $logtexfile = "$logpath$logname.tex";
my $logtexfileoutbase = "$logpath${logname}2";
my $logtexfileout = "$logtexfileoutbase.tex";

# make html version
my $command = "pandoc -s $logfile -o $logpath$logname.html";
print "$command\n";
system($command);

if (!$nolatex) {
  $command = "pandoc -s $logfile -o $logtexfile";
  print "$command\n";
  system($command);

  if (!-e "$logtexfile") {
    print "Error, tex file not created by pandoc.\n";
    exit;
  }

  print "Reading file \"$logtexfile\", writing to \"$logtexfileout\"\n";
  open(LOGFILE, "<$logtexfile") || die("Could not open tex logfile \"$logtexfile\" - $!");
  open(OUTFILE, ">$logtexfileout") || die("Could not open outfile \"$logtexfileout\" - $!");

  while(<LOGFILE>) {
    s#\\documentclass\[\]\{.*\}#\\documentclass\[a4paper\]\{book\}#;

    # hack - add hspace of 0 and it seems happier to have paragraphs/headings alone in quote environment
    s#\\begin\{quote\}#\\begin\{quote\}\\hspace\{0cm\}\\vspace\{-1cm\}#g;
    s#\\end\{quote\}#\\hspace\{0cm\}\\end\{quote\}#g;
    print OUTFILE;

    if ($_ =~ /documentclass/) {
      print OUTFILE "\\usepackage[cm]{fullpage}\n";
    }
  }

  close(OUTFILE);
  close(LOGFILE);

  # let the latex files appear in /tmp

  chdir("/tmp");

  $command = "latex $logtexfileout > /dev/null 2>&1";
  print "$command\n";
  system($command);

  $command = "dvips $logtexfileoutbase.dvi > /dev/null 2>&1";
  print "$command\n";
  system($command);

  $command = "ps2pdf $logtexfileoutbase.ps > /dev/null 2>&1";
  print "$command\n";
  system($command);
}

