if (!$ARGV[0] || $ARGV[0] eq "-h") {
	print STDERR "\n USAGE: mirScreen.pl mirBam otherBam1:type otherBam2:type"; 
	exit;
}
#the rationale here is to map reads against  different small rna datasets and provide evidence that a read is actually that a of a miRNA
#by excluding from results all reads that don't map to mirna only.
$umicheck = 1;
$cnt = 0;
@type = ();
foreach $file (@ARGV) {
	$cnt++;
	($file,$type) = split (/:/,$file);
	$type = 'mir' if ($cnt == 1);
	push(@types,$type);
	print STDERR "Loading $file ($type)...";
	open(IN,"samtools view $file |");
	%seen = ();
	$tot1 = 0;
	$tot2 = 0;
	while ($line = <IN>) {
		@tmp = split (/\t/,$line);
		#the last : split if tmp[0] is the UMI sequence
		@umi = split(/[:_]/,$tmp[0]);
		$umi = pop @umi;
		#tmp[2] is the target
		$target = $tmp[2];
#		print STDERR "$umi -> $type -> $tmp[2]\n" if ($type ne 'mir');
		#print STDERR "UMI $umi is a PCR duplicate" if ($umicheck && $umi_dup{$umi});
		#$umi_seen{$umi} = 1;
		$tot1++;
		next if ($tmp[2] eq '*');
		#in case umi are used for deduplication
		#the umi/target pair are recorded and marked for discard if they occur twice or more in mir
		$umi_dup{$umi} = 1 if ($umicheck && $type eq 'mir' && $cnt{$umi."@".$type."@".$target});
		#the umi will not be further (in types != mir) if dupli
		next if ($umi_dup);
		$umi_good{$umi} = 1; #this will be used for output
		$tot2++;
		#counts are only recorded for mir targets, other matches are simply indicated as counts per type
		#$cnt{$umi."@".$type."@".$target}++ if ($type eq 'mir');
		$cnt{$umi."@".$type}++ if ($type ne 'mir');
		$target{$umi} = $target if ($type eq 'mir');
		$cnt_tot{$type}++;
		#$target{$tmp[2]."@".$type} = 1;
		#$targets_seen{$tmp[2]} = 1;
#		print STDERR $line,"\n";
#		print STDERR $tmp[0],"\t",$tmp[2],"\n";
	}
	close IN;
	print STDERR "$type: $tot1 total reads, ",scalar keys %umi_good," non duplicated\n" if ($type eq 'mir');
	print STDERR "$type: ",$cnt_tot{$type},"/$tot1 reads have matches\n";
}
#shift @types;
print "UMI\tTARGETMIR\t",join "\t",@types,"\n";
foreach $umi (sort keys %umi_good) {
	next if !$target{$umi};
#	print $umi;
	foreach $type (@types) {
		if ($type eq 'mir') {
			print $umi,"\t",$target{$umi},"\t1";
		} else {
			print "\t",$cnt{$umi."@".$type} if $cnt{$umi."@".$type};
			print "\t0" if !$cnt{$umi."@".$type};
		}
	}
	print "\n";
}
