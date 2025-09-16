use strict;
use warnings;

use Getopt::Long;

my @library;
GetOptions(
	"main-bam|b=s" => \my $mainBam,
	"library-bam|B=s@" => \@library,
	"help|h" => \my $help,
	'verbose|v+' => \my $verbose,
) or die "Error\n";

@library = map { split /,/ } @library;

if (!@library || !$mainBam || $help) {
	print STDERR "\n USAGE: mirScreen.pl --main-bam mirBam --file otherBam1:type,otherBam2:type\n\n"; 
	print STDERR "\t-b, --main-bam\t\tbam file containing reads to check\n",
	"\t-B, --library-bam\tbam files to match the reads against separated by \':\' with a type flag to record matches\n",
	"\t-v, --verbose\t\tenable verbose mode, Default: \'False\'\n",
	"\t-h, --help\t\tdisplay help\n",
	"\n";
	print STDERR "     example: perl mirScreen.pl -v -b bam1 -B bam2:trna -B bam3:mature\n\n";
	print STDERR "     example: perl mirScreen.pl -v -b bam1 -B bam2:trna,bam3:mature\n\n";	
	exit;
}

my @files = ($mainBam,@library);

#the rationale here is to map reads against  different small rna datasets and provide evidence that a read is actually that a of a miRNA
#by excluding from results all reads that don't map to mirna only.
my $umicheck = 1;
my $cnt = 0;
my @types;
my %umi_dup = ();
my %umi_good = ();
my %cnts = ();
my %targets = ();
my %cnt_tot = ();
foreach my $file (@files) {
	$cnt++;
	my ($file,$type) = split (/:/,$file);
	$type = 'mir' if ($cnt == 1);
	push(@types,$type);
	print STDERR "Loading $file ($type)..." if ($verbose);
	open(IN,"samtools view $file |");
	my %seen = ();

	my $tot1 = 0;
	my $tot2 = 0;
	while (my $line = <IN>) {
		my @tmp = split (/\t/,$line);
		#the last : split if tmp[0] is the UMI sequence
		my @umi = split(/[:_]/,$tmp[0]);
		my $umi = pop @umi;
		#tmp[2] is the target
		my $target = $tmp[2];
#		print STDERR "$umi -> $type -> $tmp[2]\n" if ($type ne 'mir');
		#print STDERR "UMI $umi is a PCR duplicate" if ($umicheck && $umi_dup{$umi});
		#$umi_seen{$umi} = 1;
		$tot1++;
		next if ($tmp[2] eq '*');
		#in case umi are used for deduplication
		#the umi/target pair are recorded and marked for discard if they occur twice or more in mir
		$umi_dup{$umi} = 1 if ($umicheck && $type eq 'mir' && $cnts{$umi."@".$type."@".$target});
		#the umi will not be further (in types != mir) if dupli
		next if ($umi_dup{$umi});
		$umi_good{$umi} = 1; #this will be used for output
		$tot2++;
		#counts are only recorded for mir targets, other matches are simply indicated as counts per type
		#$cnt{$umi."@".$type."@".$target}++ if ($type eq 'mir');
		$cnts{$umi."@".$type}++ if ($type ne 'mir');
		$targets{$umi} = $target if ($type eq 'mir');
		$cnt_tot{$type}++;
		#$target{$tmp[2]."@".$type} = 1;
		#$targets_seen{$tmp[2]} = 1;
		#print STDERR $line,"\n";
		#print STDERR $tmp[0],"\t",$tmp[2],"\n";
	}
	close IN;
	print STDERR "$type: $tot1 total reads, ",scalar keys %umi_good," non duplicated\n" if ($type eq 'mir' && $verbose);
	print STDERR "$type: ",$cnt_tot{$type},"/$tot1 reads have matches\n" if ($verbose);
}
#shift @types;
print "UMI\tTARGETMIR\t",join "\t",@types,"\n";
foreach my $umi (sort keys %umi_good) {
	next if !$targets{$umi};
#	print $umi;
	foreach my $type (@types) {
		if ($type eq 'mir') {
			print $umi,"\t",$targets{$umi},"\t1";
		} else {
			print "\t",$cnts{$umi."@".$type} if $cnts{$umi."@".$type};
			print "\t0" if !$cnts{$umi."@".$type};
		}
	}
	print "\n";
}