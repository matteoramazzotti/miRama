#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;

Getopt::Long::Configure ("bundling","no_ignore_case");

GetOptions( 
	'help|h' => \my $help,
	'verbose|v+' => \my $verbose,
	'targets|t=s' => \my $targets, # "human_3_UTR.fasta"
	'mirnas|m=s' => \my $mirdb, # "mirbase/mature.fa"
	'filter|f=s' => \my $mirsource,
	'mode|M=s' => \my $strand,
	'slice|s=i' => \my $slice,
	'range|r=s' => \my $range,
	'out-folder|o=s' => \my $outFolder
);

if ($help || !$targets || !$mirdb) {
	print STDERR "USAGE: mirSeed.pl\n",
	"\t-t, --targets\t\ttarget_genes.fasta\n",
	"\t-m, --mirnas\t\tmicroRNA_database.fasta\n",
	"\t-f, --filter\t\tmiRNA_source\n",
	"\t-M, --mode\t\tmode selected: seed,plus or minus, Default: \'seed\'\n",
	"\t-s, --slice\t\tslice size for plus and minus modes, Default: \'6\'\n",
	"\t-r, --range\t\trange for seed mode, Default: \'2,8\'\n",
	"\t-o, --out-folder\toutput directory, Default: \'./\'\n",
	"\t-v, --verbose\t\tenable verbose mode, Default: \'False\'\n",
	"\n";
	print STDERR "     example: perl mirSeed.pl --targets human_3_UTR.fasta --mirnas mature.fa --filter bta-miR --mode seed --range 2,7 \n\n";
	print STDERR "     example: perl mirSeed.pl --targets human_3_UTR.fasta --mirnas mature.fa --filter DE_miRNA.names --mode plus --slice 8\n\n";	
	exit;
}

$mirsource = $mirsource ? $mirsource : "hsa"; #a file containing newline spearated list of microRNA names or a three-letter name organism (e.g. hsa, bta etc)
$strand = $strand ? $strand : "seed"; #plus/minus/seed
$range = $range ? $range : "2,8"; #the range of mature mirna or a range e.g. 2,8 (the seed)
$slice = $slice ? $slice : "6"; #the slicing size of mature mirna e.g. 6 (the seed)
$outFolder = $outFolder ? $outFolder : "./";


print STDERR "\nLoading gene targets: ";
open(IN,"$targets") or die "No gene target found\n";
my @targets;
my %targets;
my $name;
while (<IN>) {
	chomp;
	if ($_ =~ ">") {
		$name = $_;
		$name =~ s/>//; 
		$name =~ s/ \[.+//;
		push(@targets,$name);
	} else {
		$targets{$name} .= $_;
	}
}
close(IN);

print STDERR "\n\t",scalar keys %targets, " sequences loaded\n";
print STDERR "\nLoading miRNA seq:\n";

my $cnt = 0;
my %mirseq; 
open(IN,$mirdb) or die "No mirbase db file!\n";
while(<IN>) {
	chomp;
	$cnt++;
	if ($cnt == 1) {
		$name = $_;	
		$name =~ s/>//;
		$name =~ s/ .+//g;
	}
	if ($cnt == 2) {
		$mirseq{$name} = $_;
		$mirseq{$name} =~ s/U/T/g;
		$cnt = 0;
	}
}
close IN;
print STDERR "\t",scalar keys %mirseq, " microRNAs loaded\n\n";

# load names of microRNA of interest if $ARGV[0] = file
# else load only three letters organism name (e.g. bta, hsa etc.)

my @in_mir;
my $tot = 0;
my $valid = 0;
if (-f $mirsource) {
	open(IN,$mirsource);
	while(<IN>) {
		chomp;
		next if (!$_ || $_ =~ /^#/);
		$tot++;
		if ($mirseq{$_}) {
			push(@in_mir,$_);
			$valid++;
		}
	}
	close IN;
	print STDERR "miRNAs from file: $valid / $tot valid \n";
} else {
	#an organism name (eg hsa,bta,mmu etc.) to filter by org
	foreach my $x (keys(%mirseq)) {
		$tot++;
		if ($x =~ /$mirsource/) {
			push(@in_mir,$x);
			$valid++;
		}
	}
	print STDERR "miRNAs from regex $ARGV[0]: $valid / $tot valid\n";
}
print STDERR "\t",scalar @in_mir," microRNAs loaded\n\n";

my $ind = 0;
my %genes;
my %mir_targets;
my %gene_targeted; # for transposed %mir_targets;
my %n_gene_targeted; # for number of mirnas interacting with each gene (for later sorting)
my %n_mir_targets; # for number of genes interacting with each mirna (for later sorting)
foreach my $mir (@in_mir) {
	$ind++;
	my $cnt = 0;
	print STDERR "$ind: scanning $mir: ";
	#if "seed" only the mirna seed is extracted
	#the mir is sliced in chunks of size $slice to test seed-like perfect matches with target
	if ($strand eq "seed") {
		print STDERR "---------------MIR $mir\n";
		my @coord = split(",",$slice);
		print STDERR "$coord[0],$coord[1]\t",$mirseq{$mir},"\n";
		my $start = $coord[0];
		my $end = $coord[1];
		my $sub = substr($mirseq{$mir}, $start - 1, $end - 1);
		print STDERR $sub,"\n";
		my @res = analyze($sub);
		next if (scalar(@res) == 0);
		foreach my $i (0..$#res) {
			my @parts = split("\\|",$res[$i]);
			next if (!$parts[2]);
			push @{$genes{$mir}},$parts[2];
		}
		summarize($mir);
	} else {
		#if "strand" the mirna sequence is progressively scanned in chunks of $slice size
		foreach my $start (1..length($mirseq{$mir})) {
			if ($strand eq "minus") {
				$mirseq{$mir} =~ tr/ATGC/TACG/;
				$mirseq{$mir} = reverse($mirseq{$mir});
			}

			my $sub = substr($mirseq{$mir}, $start, $slice);
			print "sub $sub, size",length($sub),"\n";

			if (length($sub) >= $slice) {
				my @res = analyze($sub);
				next if (scalar @res == 0);
				foreach my $i (0..$#res) {
					my @parts = split("\\|",$res[$i]);
					push @{$genes{$mir}},$parts[2];
				}
				summarize($mir)
			}
		}
	}
}



# print tables
open(OUT,">","$outFolder/miRNAs_vs_genes_$slice.tsv");
print OUT "miRNA ID\tNumber of Target Genes\tTarget Genes\n"; 
foreach my $mir (sort ({$n_mir_targets{$a} <=> $n_mir_targets{$b}} keys(%n_mir_targets))) {
	my @genes = split(" ",$mir_targets{$mir});
	print OUT "$mir\t",scalar(@genes),"\t",$mir_targets{$mir},"\n";
}
close(OUT);

open(OUT,">","$outFolder/genes_vs_miRNAs_$slice.tsv");
print OUT "gene ID\tNumber of miRNAs\tmiRNAs\n"; 
foreach my $gene (sort ({$n_gene_targeted{$a} <=> $n_gene_targeted{$b}} keys(%n_gene_targeted))) {
	my @mirnas = keys(%{$gene_targeted{$gene}});
	print OUT "$gene\t",scalar(@mirnas),"\t",join(" ",@mirnas),"\n";
}
close(OUT);

sub summarize {
	#makes the list of genes unique
	my $mir = shift;
	my %h = map { $_ => 1 } @{$genes{$mir}};
	my @array_out = sort keys %h;
	#stores the list of matched genes
	$mir_targets{$mir} = join " ",@array_out;
	$n_mir_targets{$mir} = scalar(@array_out);
	print STDERR scalar(@array_out)," targets\n";

	#stores the transposed version (different format)
	foreach my $gene (@array_out) {
		${$gene_targeted{$gene}}{$mir} = 1;
		$n_gene_targeted{$gene} = scalar(keys(%{$gene_targeted{$gene}}));
	}

}
sub analyze {
	my $in=shift;
	my @res;
	#each target matching (perfect match)
	#between mirna and mRNA is collected
	foreach my $t (sort keys %targets) {
		if ($targets{$t} =~ $in) {
			push(@res,$t)
		}
	}
	return(@res);
}