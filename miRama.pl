#!/usr/bin/perl
#see https://resources.qiagenbioinformatics.com/manuals/biomedicalgenomicsanalysis/current/index.php?manual=Create_UMI_Reads_miRNA.html#sec:createumiformirna
if (!$ARGV[0] || $ARGV[0] eq "-h") {
	print STDERR "\n USAGE: miRama.pl [file[.gz]] [-targetdb file.fasta] -mode [full|prep|map] [-mis N] -out sample -outmode fastq -min 15 -max 55 -umilen 12 -common 19 -buffer 100\n\n"; 
	exit;
}

#defaults
$logcounter = 100;
$common=19;
$min=15;
$max=55;
$umilen=12;
$buffer = 100;
$out = 'sample';
$direct = 0;
$outmode = 'fastq';
$mis = 0;
$mode = 'full';

#reads can be in a file or piped in
if (-e $ARGV[0]) {
	$file = shift @ARGV;
	open($in,$file) if ($file !~ /\.gz$/);
	open($in,"zcat $file|") if ($file =~ /\.gz$/);
} else {
	$in = *STDIN;
}

#options processing
for($i=0;$i<=$#ARGV;$i++) {
	if ($ARGV[$i] eq '-targetdb') {
		$targetdb = $ARGV[$i+1];
		die "ERROR: $targetdb not found" if (!-e $targetdb);
		$i++;
	}
		if ($ARGV[$i] eq '-out') {
				$out = $ARGV[$i+1];
				$i++;
		}
	   if ($ARGV[$i] eq '-outmode') {
				$outmode = $ARGV[$i+1];
				$i++;
		}
		if ($ARGV[$i] eq '-min') {
				$min = $ARGV[$i+1];
				$i++;
		}
		if ($ARGV[$i] eq '-max') {
				$max = $ARGV[$i+1];
				$i++;
		}
		if ($ARGV[$i] eq '-umilen') {
				$umilen = $ARGV[$i+1];
				$i++;
		}
		if ($ARGV[$i] eq '-guess') {
				$guess = $ARGV[$i+1];
				$i++;
		}
		if ($ARGV[$i] eq '-common') {
				$common = $ARGV[$i+1];
				$i++;
		}
	    if ($ARGV[$i] eq '-mode') {
				$mode = $ARGV[$i+1];
		$i++
		}
	    if ($ARGV[$i] eq '-mis') {
				$mode = $ARGV[$i+1];
		$i++
		}

#	print STDERR "$ARGV[$i] -> $ARGV[$i+1]\n";
}
#exit;

&load_db if (-e $targetdb);

#output file for valid reads in fast[aq] format 
$seqout = $out.".$outmode" if ($outmode);
#output file with mapping of valid reads to $tarfetdb
$mapout = $out.".mapped.all" if ($targetdb);
#output file with counts for all targets in targetdb
$cntout1 = $out.".counts.all" if ($targetdb);
#output file with counts for mapped targets in targetdb (i.e. zero counts are not written)
$cntout2 = $out.".counts.nozero" if ($targetdb);

goto MAP if ($mode eq "map");

open(SEQOUT,">$seqout") if ($seqout);

if ($common !~ /[ATGC]/) {
	$guess = $common;
	$common = "auto";
	print STDERR "common set to $common with length $guess, scanning $buffer sequences\n";
} else {
	print STDERR "common set to $common\n";
}

#INITIALIZATION
%umis = ();
%cnt = ();
%com = ();
%mapped = ();
@buf = ();
$nreads = 0;
$outseq = 0;
$totmapped = 0;
$valid = 0;
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

while($line = <$in>) {
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
#			print STDERR $line,"\n";
			while($line =~ /(.{$guess})/g) {
				pos($line) = pos($line)-($guess-1);
				$com{$1}++;
#				print $1,"-$guess\n";
			}
#			@top = sort {$com{$b} <=> $com{$a}} keys %com;
#			print STDERR "$nreads / $buffer: ",scalar keys %com," $top[0] -> $com{$top[0]}\n";
#			<STDIN>;
			next;
		}
		if ($common eq 'auto' && $nreads == $buffer) {
			@top = sort {$com{$b} <=> $com{$a}} keys %com;
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
		$line =~ /(.+?)$common(.{$umilen})/;
		#the actual read (the part before the common sequence) is filterd by size
		if (length($1)>=$min && length($1)<=$max) {
			$validread = $1;
			$umi = $2;
			$umis{$umi}++;
			$valid++;
			if ($targetdb) {
			#scan_db returns zero if unmatched
				$map = scan_db_eq($validread);
				$mapped{$map}++ if ($map);
				$totmapped++ if ($map);
				$toprint = "$name:UMI:$umi\n$validread";
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
#post processing of buffered sequences
foreach $n (0..$#bufname) {
	$bufseq[$n] =~ /(.+?)$common(.{$umilen})/;
	$toprint = '';
		if (length($1)>=$min && length($1)<=$max) {
			$validread = $1;
				$umi = $2;
				$umis{$umi}++;
				$valid++;
				if ($targetdb) {
				#scan_db returns zero if unmatched
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
close SEQOUT;
#the last time	
if ($targetdb) {
	print STDERR "\r$nreads: $valid valid, $totmapped mapped in ",scalar keys %umis," UMIs and ",scalar keys %mapped, " targets";
} else {
	print STDERR "\r$nreads: $valid valid in ", scalar keys %umis," UMIs";
}

MAP:
#in this mode correct reads are taken from a file  
if ($mode eq 'map') {
	$choice = 0;
	$readmode = 'fastq';
	while($line = <$in>) {
		$linenum++;
		chomp $line;
		$readmode = 'fasta' if (!$readmode && $linenum == 1 && $line =~ />/);
		if ($readmode eq 'fasta') {
				if ($line =~ />/) {
						$name = $line;
				$name =~ s/^>//;
				$cnt++;
				} else {
					$map = scan_db_eq($line) if ($mis == 0);
					$map = scan_db_mis($line) if ($mis > 0);
					$mapped{$map}++ if ($map);
					$totmapped++ if ($map);
					print STDERR "\rSUMMARY: $totmapped / $cnt mapped";
				}
		}
		else {
			if ($linenum == 1) {
				$name = $line;
				$name =~ s/^>//;
				$cnt++;
			} 
			if ($linenum == 2) {
				$map = scan_db_eq($line) if ($mis == 0);
				$map = scan_db_mis($line) if ($mis > 0);
				$mapped{$map}++ if ($map);
				$totmapped++ if ($map);
				print STDERR "\rSUMMARY: $totmapped / $cnt mapped";
			}
			$linenum = 0 if ($linenum == 4);	
		}
	}
}

open (OUT,">$out.log");
print OUT "TOT\t$nreads\nVALID\t$valid\nMAPPED\t$totmapped\nUMIs\t",scalar keys %umis,"\nTARGETS\t",scalar keys %mapped,"\n" if ($targetdb);
print OUT "TOT\t$nreads\nVALID\t$valid\nUMIs\t",scalar keys %umis,"\n" if (!$targetdb);
close OUT;

if ($targetdb) {
	print "\nSUMMARY: $totmapped / $nreads mapped\n\n";
	open(OUT,">$out.counts.txt");
	foreach $name (keys %db) {
		$mapped{$name} = 0 if (!$mapped{$name});
		print OUT $name,"\t",$mapped{$name},"\n";
	}
	close OUT;
	open(OUT,">$out.counts.nozero.txt");
	foreach $name (keys %db) {
			$mapped{$name} = 0 if (!$mapped{$name});
			print OUT $name,"\t",$mapped{$name},"\n" if ($mapped{$name} > 0);
	}
	close OUT;
}

sub load_db {
	print STDERR "Loading $targetdb: ";
	open(DB,$targetdb);
	while ($line = <DB>) {
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
	print STDERR scalar keys %db," targets in databse\n";
}

sub scan_db_eq {
	my $in = shift;
	foreach $k (keys %db) {
		return $k if ($db{$k} =~ /$in/i);
	}
	return(0);
}

sub scan_db_mis {
	my $seq = shift;
	my $pos = 0;
	$lenf = length($seq);
	foreach $k (keys %db) {
		while (1) {
			$part = substr($db{$k},$pos,$lenf);
			#hamming distance betwen two strings
			$edit = ($part ^ $seq) =~ tr/\001-\255//;
			return $k if ($edit <= $mis);
			$pos++;
			last if ($pos == length($seq));
		}
	}
}

