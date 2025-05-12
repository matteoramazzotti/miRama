if (!$ARGV[0] || $ARGV[0] eq "-h") {
	print STDERR "\n USAGE: mirCollect.pl screen_folder [filter1[,N]] .. [filterN[,N]] outfile\n\n"; 
	print STDERR "screeen folder: the folder containing .screen file produced by miRscreen.pl\n";
	print STDERR "filters: columns of screen files (from 4th) to consider in filtering\n\n";
	print STDERR "N: threshold value. Default is 0l i.e. if > 0 the miR is discarded\n\n";
	print STDERR "Example: with header \"UMI|TARGETMIR|mir|mature|trna|other\"\n\n";
	print STDERR "- mirCollect.pl screen_folder trna other\n";
	print STDERR "  will dicard all miRs that match trnas and other ncRNA\n";
	print STDERR "- mirCollect.pl screen_folder trna,5 other\n";
	print STDERR "  will dicard all miRs that match more than 5 trnas and any other ncRNA\n\n";
	exit;
}
#the rationale here is to map reads against  different small rna datasets and provide evidence that a read is actually that a of a miRNA
#by excluding from results all reads that don't map to mirna only.
$debug = 0;
$log = 1;
$outfile = pop @ARGV;
open(OUT,">$outfile.out.txt");
open(LOG,">$outfile.log.txt");
print LOG "file\tmiRs\tmiRreads\ttotReads\ttotMirs\n";
$dir = shift @ARGV;
opendir DIR,$dir;
@files = readdir DIR;
foreach $p (@ARGV) {
	($f,$v) = split (/,/,$p);
	$filt{$f} = 0;
	$filt{$f} = $v if ($v);
#	print STDERR "$f: $v\n";
}
$filt_tot = scalar keys %filt;
%cnt = ();
closedir DIR;
foreach $file (sort @files) {
	next if ($file =~ /^\./);
	print LOG "$file\t";
	print STDERR "$file: " if ($debug);
	open(IN,"$dir/$file");
	$tot = 0;
	$kept = 0;
	%local = ();
	while ($line = <IN>) {
		next if ($line !~ /\w/);
		print $line if ($debug);
		chomp $line;
		$tot++;
		@tmp = split (/\t/,$line);
		@head = @tmp if ($tot == 1); #the header is skipped
		next if ($tot == 1); #the header is skipped
		$keep = -1;
		foreach $f (3..$#head) {
			$keep++ if (!exists($filt{$head[$f]}) || (exists($filt{$head[$f]}) && ($tmp[$f] <= $filt{$head[$f]})));
			print STDERR "pos$f $head[$f] - value $tmp[$f], has filter at $filt{$head[$f]}, score = $keep, expected $filt_tot\n" if exists($filt{$head[$f]}) && $debug;
			print STDERR "pos$f $head[$f] - value $tmp[$f], no filter, score = $keep, expected $filt_tot\n" if !exists($filt{$head[$f]}) && $debug;
		}
		if ($keep >= $filt_tot) {
			$cnt{$tmp[1]."@".$file}++;
			$all{$tmp[1]} = 1;
			$local{$tmp[1]} = 1;
			$kept++;
		}
		print STDERR "$tmp[1]: ",$cnt{$tmp[1]."@".$file},"\nKEPT: $kept\n" if ($debug);
		<STDIN> if ($debug);
	}
	close IN;
	print STDERR "$file: ",scalar keys %local," unique miRs loaded in $kept/$tot reads. ",scalar keys %all, " miRs available\n";
	print LOG "$file\t",scalar keys %local,"\t$kept\t$tot\t",scalar keys %all,"\n";
	<STDIN> if ($debug);
}

#OUTPUT
print OUT "\t",join "\t", sort @files,"\n";
foreach $m (sort keys %all) {
	print OUT $m;
	foreach $f (sort @files) {
		next if ($f =~ /^\./);
		print OUT "\t",$cnt{$m."@".$f} if ($cnt{$m."@".$f});
		print OUT "\t0" if (!$cnt{$m."@".$f});
	}
	print OUT "\n";
}
