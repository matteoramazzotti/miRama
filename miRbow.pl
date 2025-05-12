#!/usr/bin/perl
if (!$ARGV[0] || $ARGV[0] eq "-h") {
	print STDERR "\n USAGE: mirBow.pl dir [umi]"; 
	exit;
}
$umicheck = 1 if ($ARGV[1] && $ARGV[1] eq 'umi');

print STDERR "Using UMIs for dedup\n" if ($umicheck);
@files = split(/\n/,`ls -1 $ARGV[0] | head -n3`);
%cnt = ();
foreach $file (@files) {
	print STDERR "Loading data from $file...\n";
	open(IN,"$ARGV[0]/$file") if ($file = ~ /\.sam$/);
	open(IN,"samtools view $ARGV[0]/$file|") if ($file = ~ /\.bam$/);
	%seen = ();
	while ($line = <IN>) {
		next if $line =~ /^@/;
		next if $line !~ /miR/;
		@tmp = split (/\t/,$line);
		#the last : split if tmp[0] is the UMI sequence
		@umi = split(/[:_]/,$tmp[0]);
		$umi = pop @umi;
		#tmp[2] is the target miR
		#print STDERR "$umi -> $tmp[2]\n";
		$cnt{$tmp[2]."@".$file}++ if (!$seen{$tmp[2]."@",$umi});
		print "$tmp[2] in $umi already seen" if ($seen{$tmp[2]."@".$umi});
        $seen{$tmp[2]."@".$umi}++ if ($umicheck);
		$mirs{$tmp[2]} = 1;
#		print STDERR $line,"\n";
#		print STDERR $tmp[0],"\t",$tmp[2],"\n";	
	}
}
foreach $k (keys %seen) {
	print $k,"\t",$seen{$k},"\n" if ($seen{$k} > 1);
}
print "\t",join "\t",@files,"\n";

foreach $m (sort keys %mirs) {
	print $m;
	foreach $s (@files) {
		print "\t",$cnt{$m."@".$s} if ($cnt{$m."@".$s});
		print "\t0" if (!$cnt{$m."@".$s});
	}
	print "\n";
}
