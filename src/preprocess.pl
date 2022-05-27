use strict;

die 'Sintassi: adaptBiogridWithScoredSeedGene Biogrid seedGene outLink outGene\n' if (@ARGV < 4);

my $biogrid = shift @ARGV;
my $seedGene = shift @ARGV;
my $outLink = shift @ARGV;
my $outGene = shift @ARGV;

my @mat;
my $id1;
my $id2;
my %gene;
my $ngene = 0;
open FLINK, "<$biogrid";
open FoutLINK, ">$outLink";
while (<FLINK>) {
  chomp;
  my ($node1, $node2) = split /\s+/;
  unless ($node1 eq $node2) {
    if (!defined($gene{$node1})) {
      $gene{$node1} = $ngene;
      $ngene++;
    }
    if (!defined($gene{$node2})) {
      $gene{$node2} = $ngene;
      $ngene++;
    }
    $id1 = $gene{$node1};
    $id2 = $gene{$node2};
    if ($mat[$id1][$id2] == 0) {
      $mat[$id1][$id2] = 1;
      $mat[$id2][$id1] = 1;
      print FoutLINK "$id1 $id2\n";
    }
  }
}
close FLINK;
close FoutLINK;

open FIN, "<$seedGene";
my ($name_gene, $score);
my %scoreSeedGene;
my $maxScore = 0;
my $nseedgene = 0;
my $notfoundseedgene = 0;
while (<FIN>) {
  chomp;
  ($name_gene, $score) = split /\s+/;
  if (defined($gene{$name_gene})) {
    $scoreSeedGene{$name_gene} = $score;
    $nseedgene++;
  }
  else {
    print "Error, not found seed gene $name_gene\n";
    $notfoundseedgene++;
  }
  $maxScore = $score if ($score > $maxScore);
}
close FIN;
print "$notfoundseedgene seed genes not found\n";
print "$nseedgene seed genes present\n";
print "Maximum score $maxScore\n";

open FoutGene, ">$outGene";
for (keys %gene) {
  if (defined($scoreSeedGene{$_})) {
    my $adaptScore = $maxScore - $scoreSeedGene{$_};
    print FoutGene "$gene{$_} $_ $scoreSeedGene{$_}\n";
  }
  else {
    print FoutGene "$gene{$_} $_ 0.0\n";
  }
}
close FoutGene;


