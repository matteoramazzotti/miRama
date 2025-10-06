#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

Getopt::Long::Configure ("bundling","no_ignore_case");

#see https://resources.qiagenbioinformatics.com/manuals/biomedicalgenomicsanalysis/current/index.php?manual=Create_UMI_Reads_miRNA.html#sec:createumiformirna

GetOptions( 
	'help|h' => \my $help,
	'verbose|v+' => \my $verbose,
	'fastq|f=s' => \my $file,
	'targetdb|t=s' => \my $targetdb,
	'mode|m=s' => \my $mode,
	'out|o=s' => \my $out,
	'outmode|O=s' => \my $outmode,
	'min=i' => \my $min,
	'max=i' => \my $max,
	'umilen=i' => \my $umilen,
	'guess' => \my $guess,
	'common=s' => \my $common,
	'mis=i' => \my $mis
);


#defaults
my $direct = 0;
my $logcounter = 100;
my $buffer = 100;
$common = $common ? $common : 19;
$min = $min ? $min : 15;
$max = $max ? $max : 55;
$umilen = $umilen ? $umilen : 12;
$out = $out ? $out : 'sample';
$outmode = $outmode ? $outmode : 'fastq';
$mis = $mis ? $mis : 0;
$mode = $mode ? $mode : 'full';


if ((!$file && (-t STDIN)) || $help) {
	print STDERR "\n USAGE: miRama.pl [file[.gz]] [-targetdb file.fasta] -mode [full|prep|map] [-mis N] -out sample -outmode fastq -min 15 -max 55 -umilen 12 -common 19 -buffer 100\n\n"; 
	exit;
}



#reads can be in a file or piped in
my $in;
if (-e $file) {
	open($in,$file) if ($file !~ /\.gz$/);
	open($in,"zcat $file|") if ($file =~ /\.gz$/);
} else {
	$in = *STDIN;
}

my %db = load_db(db => $targetdb) if (-e $targetdb);

my $seqout = $out.".$outmode"; # output file for valid reads in fast[aq] format 
my $mapout = $out.".mapped.all" if ($targetdb); # output file with mapping of valid reads to $tarfetdb
my $cntout1 = $out.".counts.all" if ($targetdb); # output file with counts for all targets in targetdb
my $cntout2 = $out.".counts.nozero" if ($targetdb); # output file with counts for mapped targets in targetdb (i.e. zero counts are not written)

# $mapout $cntout1 $cntout2 unused

#INITIALIZATION

#Illumina reads for microRNA analysis from Qiagen should be
#GACTCNTAGCGGTGGA|AACTGTAGGCACCATCAAT|AAAACTCCAGCT|AGATCGGAAGAGCACACGTCTGAACTC
#RRRRRRRRRRRRRRRR|CCCCCCCCCCCCCCCCCCC|UUUUUUUUUUUU|TTTTTTTTTTTTTTTTTTTTTTTTTTT
#R: actual read
#C: common sequence
#U: UMI sequenmce
#T: trash at end

#the common sequence can be unknown, so I devised a strategy
#for autodetect the common sequence using kmer (of a given length) 
#frequencies on a number (buffer) of reads


# my %cnt = ();
# my %com = ();
# my @buf = ();
my $nreads = 0;
my $outseq = 0;
my $totmapped = 0;
my $valid = 0; 
#my $cnt = 0;



# goto MAP if ($mode eq "map");
if ($mode eq "map") {
	my %mapped = reads_mapper(input=> $in,mis=> $mis);
	# open (OUT,">$out.log");
	# print OUT "TOT\t$nreads\nVALID\t$valid\nMAPPED\t$totmapped\nUMIs\t",scalar keys %umis,"\nTARGETS\t",scalar keys %mapped,"\n" if ($targetdb);
	# print OUT "TOT\t$nreads\nVALID\t$valid\nUMIs\t",scalar keys %umis,"\n" if (!$targetdb);
	# close OUT;
	exit;
}


if ($common !~ /[ATGC]/) {
	$guess = $common;
	$common = "auto";
	print STDERR "common set to $common with length $guess, scanning $buffer sequences\n";
} else {
	print STDERR "common set to $common\n";
}
open(SEQOUT,">",$seqout);
my (@bufname,@bufseq,@bufqual,%mapped,%umis) = reads_scanner(
	input => $in,
	target_db => $targetdb
);
#post processing of buffered sequences

reads_post_processor(
	names => \@bufname,
	sequences => \@bufseq,
	quality => \@bufqual,
	min => $min,
	max => $max,
	mis => $mis,
	common => $common,
	umilen => $umilen,
	mapped => \%mapped,
	umis => \%umis
);
close(SEQOUT);




# MAP:
# #in this mode correct reads are taken from a file  
# if ($mode eq 'map') {
# 	my $choice = 0; #?
# 	my $readmode = 'fastq';
# 	my $linenum = 0;
# 	while(my $line = <$in>) {
# 		$linenum++;
# 		chomp $line;
# 		$readmode = 'fasta' if (!$readmode && $linenum == 1 && $line =~ />/);
# 		if ($readmode eq 'fasta') {
# 			if ($line =~ />/) {
# 				$name = $line;
# 				$name =~ s/^>//;
# 				$cnt++;
# 			} else {
# 				my $map;
# 				$map = scan_db_eq($line) if ($mis == 0);
# 				$map = scan_db_mis($line) if ($mis > 0);
# 				$mapped{$map}++ if ($map);
# 				$totmapped++ if ($map);
# 				print STDERR "\rSUMMARY: $totmapped / $cnt mapped";
# 			}
# 		} else {
# 			if ($linenum == 1) {
# 				$name = $line;
# 				$name =~ s/^>//;
# 				$cnt++;
# 			} 
# 			if ($linenum == 2) {
# 				my $map;
# 				$map = scan_db_eq($line) if ($mis == 0);
# 				$map = scan_db_mis($line) if ($mis > 0);
# 				$mapped{$map}++ if ($map);
# 				$totmapped++ if ($map);
# 				print STDERR "\rSUMMARY: $totmapped / $cnt mapped";
# 			}
# 			$linenum = 0 if ($linenum == 4);	
# 		}
# 	}
# }


open (OUT,">$out.log");
print OUT "TOT\t$nreads\nVALID\t$valid\nMAPPED\t$totmapped\nUMIs\t",scalar keys %umis,"\nTARGETS\t",scalar keys %mapped,"\n" if ($targetdb);
print OUT "TOT\t$nreads\nVALID\t$valid\nUMIs\t",scalar keys %umis,"\n" if (!$targetdb);
close OUT;

if ($targetdb) {
	print "\nSUMMARY: $totmapped / $nreads mapped\n\n";
	open(OUT,">$out.counts.txt");
	foreach my $name (keys %db) {
		$mapped{$name} = 0 if (!$mapped{$name});
		print OUT $name,"\t",$mapped{$name},"\n";
	}
	close OUT;
	open(OUT,">$out.counts.nozero.txt");
	foreach my $name (keys %db) {
			$mapped{$name} = 0 if (!$mapped{$name});
			print OUT $name,"\t",$mapped{$name},"\n" if ($mapped{$name} > 0);
	}
	close OUT;
}
# map mode outputs this? 
#SUMMARY: 30809 / 4363692 mapped
#SUMMARY: 30809 / 0 mapped
sub load_db {
	my %args = (
		db => "",
		@_
	);
	my $targetdb = $args{db};
	print STDERR "Loading $targetdb: ";
	open(DB,$targetdb);
	my $name;
	my %db;
	while (my $line = <DB>) {
		chomp $line;
		if ($line =~ />/) {
			$name = $line;
			$name =~ s/>//;
		} else {
			$line =~ s/U/T/g;
			$db{$name} .= $line;
		}
	}
	close DB;
	print STDERR scalar keys %db," targets in database\n";
	return(%db);
}

sub reads_scanner {
	my %args = (
		input => undef,
		target_db => undef,
		@_
	);

	my $in = $args{input};
	my $targetdb = $args{target_db};
	my $name;
	my $toprint;
	my @bufname;
	my @bufseq;
	my @bufqual;
	my $validread;
	my $umi;
	my %mapped = ();
	my %umis = ();
	my %com = ();
	my $cnt = 0;

	while(my $line = <$in>) {
		chomp $line;
		$cnt++;
		if ($cnt == 1) {
			$line =~ s/^@/>/ if ($outmode eq 'fasta');
			$name = $line;
			$name =~ s/ /:/;
			$toprint = '';
		}
		if ($cnt == 2) {
			$nreads++;
			if ($common eq 'auto' && $nreads < $buffer) {
				#these sequences are only used for guessing common and should be recomputed at the end, so they are put in a buffer
				push(@bufname,$name);
				push(@bufseq,$line);
				#print STDERR $line,"\n";
				while($line =~ /(.{$guess})/g) {
					pos($line) = pos($line)-($guess-1);
					$com{$1}++;
					#print $1,"-$guess\n";
				}
				#@top = sort {$com{$b} <=> $com{$a}} keys %com;
				#print STDERR "$nreads / $buffer: ",scalar keys %com," $top[0] -> $com{$top[0]}\n";
				#<STDIN>;
				next;
			}
			if ($common eq 'auto' && $nreads == $buffer) {
				my @top = sort {$com{$b} <=> $com{$a}} keys %com;
				$common = $top[0];
					print STDERR "\ncommon is $common with count $com{$common} / $buffer\n\n";
				next;
			}
			#the $common is not auto if sufficient reads have been evaluated and the true common sequence has been identified
			next if ($common eq 'auto');
			#the common sequence is established (or given, anyway not 'auto')
			$toprint = '';
			$validread = '';
			#$valid = 0;
			#the actual read (the part before the common sequence) is filterd by size
			if ($line =~ /(.+?)$common(.{$umilen})/) {
				if (length($1)>=$min && length($1)<=$max) {
					$validread = $1;
					$umi = $2;
					$umis{$umi}++;
					$valid++;
					if ($targetdb) {
						#scan_db returns zero if unmatched
						my $map = scan_db_eq($validread);
						$mapped{$map}++ if ($map);
						$totmapped++ if ($map);
						$toprint = "$name:UMI:$umi\n$validread";
					}
				}
			}
		}
		if ($cnt == 4) {
			if ($toprint) {
				print SEQOUT $toprint,"\n+\n",substr($line,0,length($validread)),"\n" if ($outmode eq 'fastq');
				print SEQOUT $toprint,"\n" if ($outmode eq 'fasta');
			}
			push(@bufqual,$line);
			$cnt = 0;
		}
		if ($targetdb) {
			print STDERR "\r$nreads: $valid valid, $totmapped mapped in ",scalar keys %umis," UMIs and ",scalar keys %mapped, " targets" if ($nreads > $buffer && $nreads % $logcounter == 0);
		} else {
			print STDERR "\r$nreads: $valid valid in ", scalar keys %umis," UMIs" if ($nreads > $buffer && $nreads % $logcounter == 0);
		}
	}
	return(@bufname,@bufseq,@bufqual,%mapped,%umis);
}
sub reads_post_processor {
	my %args = (
		names => [],
		sequences => [], 
		quality => [],
		min => 15,
		max => 55,
		umilen => 12,
		common => 19,
		mapped => (),
		umis => (),
		mis => 0,
		@_
	);
	my @bufname = @{$args{names}};
	my @bufqual = @{$args{quality}};
	my @bufseq = @{$args{sequences}};
	my %mapped = %{$args{mapped}};
	my %umis = %{$args{umis}};
	my $min = $args{min};
	my $max = $args{max};
	my $mis = $args{mis};
	my $common = $args{common};
	my $umilen = $args{umilen};

	my $validread;
	my $toprint;
	my $umi;


	foreach my $n (0..$#bufname) {
		if ($bufseq[$n] =~ /(.+?)$common(.{$umilen})/) {
			$toprint = '';
			if (length($1)>=$min && length($1)<=$max) {
				$validread = $1;
				$umi = $2;
				$umis{$umi}++;
				$valid++;
				if ($targetdb) {
					#scan_db returns zero if unmatched
					my $map;
					$map = scan_db_eq($validread) if ($mis == 0);
					$map = scan_db_mis($validread) if ($mis > 0);
					$mapped{$map}++ if ($map);
					$totmapped++ if ($map);
					$toprint = "$bufname[$n]|umi:$umi\n$validread";
				}
			}
			if ($toprint) {
				print SEQOUT $toprint,"\n+\n",substr($bufqual[$n],0,length($validread)),"\n" if ($outmode eq 'fastq');
				print SEQOUT $toprint,"\n" if ($outmode eq 'fasta');
			}
		}
	}
	#the last time	
	if ($targetdb) {
		print STDERR "\r$nreads: $valid valid, $totmapped mapped in ",scalar keys %umis," UMIs and ",scalar keys %mapped, " targets";
	} else {
		print STDERR "\r$nreads: $valid valid in ", scalar keys %umis," UMIs";
	}
	return(%mapped);
}
sub reads_mapper {
	my %args = (
		input => undef,
		mis => 0,
		@_
	);
	my $in = $args{input};
	my $mis = $args{mis};
	my %mapped = ();

	my $name;
	my $cnt = 0;
	my $choice = 0; #?
	my $readmode = 'fastq';
	my $linenum = 0;
	while(my $line = <$in>) {
		$linenum++;
		chomp $line;
		$readmode = 'fasta' if (!$readmode && $linenum == 1 && $line =~ />/);
		if ($readmode eq 'fasta') {
			if ($line =~ />/) {
				$name = $line;
				$name =~ s/^>//;
				$cnt++;
			} else {
				my $map;
				$map = scan_db_eq($line) if ($mis == 0);
				$map = scan_db_mis($line) if ($mis > 0);
				$mapped{$map}++ if ($map);
				$totmapped++ if ($map);
				print STDERR "\rSUMMARY: $totmapped / $cnt mapped";
			}
		} else {
			if ($linenum == 1) {
				$name = $line;
				$name =~ s/^>//;
				$cnt++;
			} 
			if ($linenum == 2) {
				my $map;
				$map = scan_db_eq($line) if ($mis == 0);
				$map = scan_db_mis($line) if ($mis > 0);
				$mapped{$map}++ if ($map);
				$totmapped++ if ($map);
				print STDERR "\rSUMMARY: $totmapped / $cnt mapped";
			}
			$linenum = 0 if ($linenum == 4);	
		}
	}
	return(%mapped);
}

sub scan_db_eq {
	my $in = shift;
	foreach my $k (keys %db) {
		return $k if ($db{$k} =~ /$in/i);
	}
	return(0);
}

sub scan_db_mis {
	my $seq = shift;
	my $pos = 0;
	my $lenf = length($seq);
	foreach my $k (keys %db) {
		while (1) {
			my $part = substr($db{$k},$pos,$lenf);
			#hamming distance betwen two strings
			my $edit = ($part ^ $seq) =~ tr/\001-\255//;
			return $k if ($edit <= $mis);
			$pos++;
			last if ($pos == length($seq));
		}
	}
}
