#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use bytes;
use Encode;
use FindBin;
use Getopt::Long;
#use File::HomeDir;
use Sys::Hostname;

my $reffasta;
my $basevcf;
my $list;
my $rprefix;
my $workingdir;
my $tsid;
my $cpu = 1;
my $help;
my $debug;
my $policy = "keep_larger2";		# if multiple SVs, keep larger 2 SVs and dicard other

GetOptions('--refseq=s' => \$reffasta, '--base=s' => \$basevcf, '--output=s' => \$rprefix, '--wdir=s' => \$workingdir, '--seqid=s' => \$tsid, '--list=s' => \$list, '--help' => \$help);

my $err = 0;
if(! $reffasta || ! -e $reffasta){
	print "! missing --refseq [reference.fasta]...\n";
	$err++;
}
if(! $basevcf || ! -e $basevcf){
	print "! missing --base [VCF file], will skip integration with base VCF...\n";
}
if(! $list || ! -e $list){
	print "! missing --list [list (.tsv)] that describe the path of VCF(s)...\n";
	$err++;
}
if(! $tsid){
	print "! missing --seqid [target seq ID]...\n";
	$err++;
}
if(! $cpu || $cpu =~ /\D/){
	$cpu = 1;
}
if($err > 0){
	Help();
	goto END;
}

print "\n";
unless(-e "./$workingdir/log"){
	system("mkdir ./$workingdir/log");
}

my ($hdbseq, $SID) = Open_fasta_as_hash($reffasta, $tsid);
my ($htarget, $TF, $errv) = Read_list($list);
if($errv > 0){
	print "! abort script due to missing file...\n";
	goto END;
}

print "\n! collecting data from target VCF...\n";
my $htdata = {};
foreach my $sample (@{$TF}){
	$htdata = Read_target($sample, $htarget->{$sample}, $htdata, $tsid, $hdbseq->{$tsid}{seq});
}

print "\n! integraing variants between targets...\n";
my $hvarmap = $htdata->{varmap};
#my @TMPK = keys(%{$hvarmap});
#my $tmpk = join(" ", @TMPK);
#print "$tmpk\n";

my $num_SID = @{$SID};
if(! $cpu || $cpu =~ /\D/){
	$cpu = $num_SID;
}
if(! $rprefix){
	$rprefix = "each_multisample_".$tsid;
}
my $rfile0 = "./$workingdir/log/".$rprefix.".vcf";

my $d = int($num_SID / $cpu) + 1;
my $multivcf =<<"EOS";
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=CONFLICT,Number=.,Type=String,Description="Sample names for which there are multiple paths in the graph with conflicting alleles">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
EOS
$multivcf .= $htdata->{header};
$multivcf .= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".join("\t", @{$TF})."\n";
for(my $i = 0; $i < $d; $i++){
	my $p0 = $i * $cpu;
	my $p1 = ($i + 1) * $cpu;
	
	my $thrs = [];
	for(my $j = $p0; $j < $p1; $j++){
		if($SID->[$j]){
			if(! $hvarmap->{DEL}{$SID->[$j]} && ! $hvarmap->{INS}{$SID->[$j]}){
				print "! no data for [$SID->[$j], skip...\n";
			}
			else{
				if($hvarmap->{DEL}{$SID->[$j]}){
					my $thr = threads->new(\&MergeSVblock, $SID->[$j], $TF, $hvarmap->{DEL}{$SID->[$j]}, $htdata->{data}{$SID->[$j]}, $hdbseq->{$SID->[$j]}{seq}, $j, "DEL", $policy);
					push(@{$thrs}, $thr);
				}
				if($hvarmap->{INS}{$SID->[$j]}){
					my $thr = threads->new(\&MergeSVblock, $SID->[$j], $TF, $hvarmap->{INS}{$SID->[$j]}, $htdata->{data}{$SID->[$j]}, $hdbseq->{$SID->[$j]}{seq}, $j, "INS", $policy);
					push(@{$thrs}, $thr);
				}
			}
		}
	}
	if(@{$thrs}){
		foreach my $thr (@{$thrs}){
			$multivcf .= $thr->join();
		}
	}
}
if(! $multivcf){
	print "! missing multisample VCF, abort script...\n";
	goto END;
}
else{
	SAVE($rfile0, $multivcf);
}

if($basevcf && -e $basevcf){
	print "\n! merge [$rfile0] to [$basevcf]...\n";
	my $prefixb = Path2pref($basevcf);
	my $rfile1 = "./$workingdir/log/".$prefixb."_plus_".$rprefix.".vcf";
	NrIntegrateVCF($basevcf, $multivcf, $TF, $SID, $rfile1);
}

END:{
	print "\n! End of script.\n\n";
}


################################################################################
#-------------------------------------------------------------------------------
sub Help{
my $mode = shift;

my $help_command;
if(! $mode || $mode eq '1'){
	$help_command =<<"EOS";

Basic usage: 
  /path/to/Vcfintegrate_cactusMG_asm2sv.pl -r [refseq.fasta] -b [base VCF] -l [list of VCF (.tsv)]

--refseq or -r  : reference fasta
--base or -b    : base VCF (assume it originates from cactus-minigraph, optional)
--list or -l    : list of target VCF that will be merged or integrated to base VCF (mandatory)
--help or -h    : display usage

[note] If variants or SV are overlapping in the same genomic region in both base and target VCF(s), data from base will be removed then those of target VCF is described.

[example of list file (tab-delimited tsv)]
sample1	/path/to/sample1.vcf
sample2	/path/to/sample2.vcf
sample3	/path/to/sample3.vcf

column1: sample name (this must be same with those described in base VCF)
column2: file path
EOS
}

print $help_command;

}


#-------------------------------------------------------------------------------
sub NrIntegrateVCF{
my ($basevcf, $multivcf, $TF, $SID, $rfile1) = @_;

my $err = 0;
my $hbp = {};
my $nrh = {};
my $header = "";
my $AoVR = [];
if($multivcf){
	print "! loading data...\n";
	my $cnt = 0;
	my @L = split(/\n/, $multivcf);
	foreach my $line (@L){
		if(! $line){
			next;
		}
		my @A = split(/\t/, $line);
		my $numA = @A;
		if($numA < 10){
			if($line !~ /CHROM/ && ! $nrh->{$line}){
				$header .= $line."\n";
				$nrh->{$line} = 1;
			}
			next;
		}
		if($line =~ /CHROM/){
			next;
		}
		
		my $pstart = $A[1];
		my $pend = $A[1] + length($A[3]);
		for(my $i = $pstart; $i < $pend; $i++){
			$hbp->{$A[0]}{$i} = 1;
		}
		push(@{$AoVR}, \@A);
		$cnt++;
	}
	print " - [$cnt] variants\n";
}
if(-e $basevcf){
	my $fh;
	if($basevcf =~ /\.gz/){
		open($fh, "gzip -dc $basevcf |") or die;
	}
	else{
		open($fh, "<", $basevcf) or die;
	}
	
	print "! reading [$basevcf]...\n";
	my $hcol = {};
	my $cnt_selected = 0;
	my $cnt_overlap = 0;
	my $cnt_nogenotype = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if(! $line){
			next;
		}
		my @A = split(/\t/, $line);
		my $numA = @A;
		if($numA < 10){
			if($line !~ /CHROM/ && ! $nrh->{$line}){
				$header .= $line."\n";
				$nrh->{$line} = 1;
			}
			next;
		}
		if($line =~ /CHROM/){
			for(my $i = 9; $i < $numA; $i++){
				$hcol->{$A[$i]} = $i;
			}
			my $num_sample = 0;
			foreach my $sample (@{$TF}){
				if(defined $hcol->{$sample}){
					$num_sample++;
				}
			}
			if($num_sample == 0){
				$err = 1;
				last;
			}
			next;
		}
		
		my $judge_overlap = "false";
		my $pstart = $A[1];
		my $pend = $A[1] + length($A[3]);
		for(my $i = $pstart; $i < $pend; $i++){
			if($hbp->{$A[0]}{$i}){
				$judge_overlap = "true";
			}
		}
		if($judge_overlap eq 'true'){
			$cnt_overlap++;
			next;
		}
		
		my @RA;
		for(my $i = 0; $i < 9; $i++){
			push(@RA, $A[$i]);
		}
		my $num_genotyped = 0;
		foreach my $sample (@{$TF}){
			if(defined $hcol->{$sample}){
				my $j = $hcol->{$sample};
				if(defined $A[$j] && $A[$j] ne '.'){
					push(@RA, $A[$j]);
					$num_genotyped++;
				}
				else{
					push(@RA, ".");
				}
			}
		}
		if($num_genotyped == 0){
			$cnt_nogenotype++;
		}
		else{
			push(@{$AoVR}, \@RA);
			$cnt_selected++;
		}
	}
	close $fh;
	
	print " - [$cnt_overlap] overlapping (skip)\n";
	print " - [$cnt_nogenotype] without genotype (skip)\n";
	print " - [$cnt_selected] selected\n";
}
if(@{$AoVR}){
	print "! sorting combined VCF...\n";
	@{$AoVR} = sort {$a->[1] <=> $b->[1]} @{$AoVR};
	@{$AoVR} = sort {$a->[0] cmp $b->[0]} @{$AoVR};
	
	my $hdata = {};
	my $cnt_total = 0;
	foreach my $A (@{$AoVR}){
		$hdata->{$A->[0]} .= join("\t", @{$A})."\n";
		$cnt_total++;
	}
	
	my $rvcf = $header;
	$rvcf .= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".join("\t", @{$TF})."\n";
	my $num_data = 0;
	foreach my $sid (@{$SID}){
		if($hdata->{$sid}){
			$rvcf .= $hdata->{$sid};
			$hdata->{$sid} = "";
			$num_data++;
		}
	}
	if($num_data > 0){
		print " - [$cnt_total] nr variants\n";
		open(my $rfh, ">", $rfile1);
		print $rfh $rvcf;
		close $rfh;
		print "! output [$rfile1]\n";
	}
	else{
		print "! missing data (no vcf)...\n";
	}
}

}


#-------------------------------------------------------------------------------
sub MergeSVblock{
my ($sid, $TF, $hposmap, $hnativevcf, $dbseq, $j, $indel, $policy) = @_;

print " thread [$j] : [$sid] [$indel]...\n";
my @POS = keys(%{$hposmap});
@POS = sort {$a <=> $b} @POS;

my $debug;
my $hrdata = {};
my $nrvar = {};
my $nrvpid = {};
my $nrbp = {};
my @VPID;
my $cnt_revpos = 0;
my $cnt_shrink = 0;
foreach my $bp (@POS){
	if($nrbp->{$bp} || ! defined $hposmap->{$bp}{pos0} || ! defined $hposmap->{$bp}{pos1}){
		next;
	}
	
	my $pstart = $hposmap->{$bp}{pos0};
	my $pend = $hposmap->{$bp}{pos1};
	my $subPSE = {};
	my $subvaridL = "";
	my $subnrv = {};
	$subPSE->{$pstart} = 1;
	$subPSE->{$pend} = 1;
	
	my $maxelem = 0;
	my $judge_enough = "false";
	for(my $t = 0; $t < 100; $t++){
		for(my $i = $pstart; $i <= $pend; $i++){
			if(! $nrbp->{$i} && $hposmap->{$i}){
				if($hposmap->{$i}{pos0} < $pstart){
					$pstart = $hposmap->{$i}{pos0};
				}
				if($hposmap->{$i}{pos1} > $pend){
					$pend = $hposmap->{$i}{pos1};
				}
				$nrbp->{$i} = 1;
			}
		}
		if($pstart != $hposmap->{$bp}{pos0} || $pend != $hposmap->{$bp}{pos1}){
			$subPSE->{$pstart} = 1;
			$subPSE->{$pend} = 1;
			$hposmap->{$bp}{pos0} = $pstart;
			$hposmap->{$bp}{pos1} = $pend;
			$cnt_revpos++;
			next;
		}
		else{
			$subPSE->{$pstart} = 1;
			$subPSE->{$pend} = 1;
			
			for(my $i = $pstart; $i <= $pend; $i++){
				$hposmap->{$i}{pos0} = $pstart;
				$hposmap->{$i}{pos1} = $pend;
				
				if($hposmap->{$i}{varid}){
					my $tmpvarinfo = $hposmap->{$i}{varid};
					if(! $subnrv->{$tmpvarinfo}){
						$subvaridL .= $tmpvarinfo;
						$subnrv->{$tmpvarinfo} = 1;
					}
				}
				$nrbp->{$i} = 1;
			}
			$judge_enough = "true";
			last;
		}
	}
	if($judge_enough eq 'false'){
		print "! round search not enough...\n";
		die;
	}
	elsif($subvaridL){
		my @PSE = keys(%{$subPSE});
		@PSE = sort {$a <=> $b} @PSE;
		my $vpid = join(">", @PSE);
		
		my @VI = split(/\n/, $subvaridL);
		foreach my $varstr (@VI){
			unless($varstr){
				next;
			}
			my @Info = split(/\t/, $varstr);
			unless($Info[1]){
	#			print "! missing varid for $Info[0]\n$varstr\n";
				next;
			}
			my $sample = $Info[0];
			my $varid = $Info[1];
			if(! $nrvar->{$vpid}{$varstr}){
				$hrdata->{$vpid}{pos0} = $pstart;
				$hrdata->{$vpid}{pos1} = $pend;
				$hrdata->{$vpid}{cnt}{$sample} += 1;
				$hrdata->{$vpid}{vcf}{$sample} .= $hnativevcf->{$sample}{$varid}."\n";
				$nrvar->{$vpid}{$varstr} = 1;
				
				if(! $nrvpid->{$vpid}){
					push(@VPID, $vpid);
					$nrvpid->{$vpid} = 1;
				}
			}
		}
	}
}

my $num_VPID = 0;
my $rvcf = "";
foreach my $vpid (@VPID){
	my $pstart = $hrdata->{$vpid}{pos0};
	my $pend = $hrdata->{$vpid}{pos1};
	my $seqREF = substr($dbseq, $pstart - 1, $pend - $pstart + 1);
	my $lenR = length($seqREF);
	
	my $hsubvars = {};
	my $numvar = 0;
	my $num_allalt = 0;
	my @GTs;
	foreach my $sample (@{$TF}){
		if($hrdata->{$vpid}{vcf}{$sample}){
			my @eVCF = split(/\n/, $hrdata->{$vpid}{vcf}{$sample});
			my $AoV = [];
			foreach my $line (@eVCF){
				my @A = split(/\t/, $line);
				my $sublenA = length($A[4]);
				my @tmp;
				push(@tmp, $A[1]);
				push(@tmp, $A[3]);
				push(@tmp, $A[4]);
				push(@tmp, $line);
				push(@{$AoV}, \@tmp);
			}
			@{$AoV} = sort {$a->[0] <=> $b->[0]} @{$AoV};
			
			my $seqALT = "";
			my $prev_refpos = $pstart;
			my $i = 0;
			my $probe;
			foreach my $A (@{$AoV}){
				my $test_seqALT = substr($dbseq, $prev_refpos - 1, $A->[0] - $prev_refpos);
				my $lentest = length($test_seqALT);
				if($A->[0] < $prev_refpos){
#					print "! contradicting position: $sample | $vpid | $hrdata->{$vpid}{pos0} - $hrdata->{$vpid}{pos1} | $sid | $A->[0] < $prev_refpos\n";
#					foreach my $E (@{$AoV}){
#						my $errlenR = length($E->[1]);
#						print " $E->[0] $errlenR\n";
#					}
					next;
				}
				my $tmplenA2 = length($A->[2]);
				$probe .= " $vpid | $sample | $prev_refpos - 1, $A->[0] - $prev_refpos\n";
				$probe .= " $vpid | $sample | len(ALT) = $tmplenA2\n";
				
				$seqALT .= substr($dbseq, $prev_refpos - 1, $A->[0] - $prev_refpos);
				$seqALT .= $A->[2];
				my $sublenR = length($A->[1]);
				$prev_refpos = $A->[0] + $sublenR;
				$i++;
			}
			$seqALT .= substr($dbseq, $prev_refpos - 1, $pend - $prev_refpos + 1);
			$probe .= " $vpid | $sample | $prev_refpos - 1, $pend - $prev_refpos + 1\n";
			
			my $sublenA = length($seqALT);
			if($sublenA > 1000000){
#				print "\n$sample | $sublenA\n";
#				print "$probe\n";
#				die;
				next;
			}
			
			my $gtnum = ".";
			if(! $hsubvars->{$seqALT}){
				$numvar++;
				$hsubvars->{$seqALT}{n} = $numvar;
				$gtnum = $numvar;
			}
			else{
				$gtnum = $hsubvars->{$seqALT}{n};
			}
			$hsubvars->{$seqALT}{cnt} += 1;
			$num_allalt++;
			push(@GTs, $gtnum);
		}
		else{
			push(@GTs, ".");
		}
	}
	
	my @ALTs = keys(%{$hsubvars});
	my $AoALTs = [];
	foreach my $seqALT (@ALTs){
		my $sublenA = length($seqALT);
		my $lendfabs = abs($lenR - $sublenA);
		my @tmp;
		push(@tmp, $hsubvars->{$seqALT}{n});
		push(@tmp, $hsubvars->{$seqALT}{cnt});
		push(@tmp, sprintf("%.1f", $hsubvars->{$seqALT}{cnt} / $num_allalt) );
		push(@tmp, $seqALT);
		push(@tmp, $sublenA);
		push(@tmp, $lendfabs);
		push(@{$AoALTs}, \@tmp);
	}
	
	if($policy eq 'keep_larger2'){
		@{$AoALTs} = sort {$a->[5] <=> $b->[5]} @{$AoALTs};
		
		my $hkeepn = {};
		my $revn = 1;
		my $AoALTs2 = [];
		for(my $i = 0; $i < 2; $i++){
			if($AoALTs->[$i]){
				my $oldn = $AoALTs->[$i][0];
				$hkeepn->{$oldn} = $revn;
				
				$AoALTs->[$i][0] = $revn;
				push(@{$AoALTs2}, $AoALTs->[$i]);
			}
			$revn++;
		}
		
		my @GTs2;
		foreach my $gtnum (@GTs){
			if($hkeepn->{$gtnum}){
				push(@GTs2, $hkeepn->{$gtnum});
			}
			else{
				push(@GTs2, ".");
			}
		}
		@GTs = @GTs2;
		
		@{$AoALTs2} = sort {$a->[0] <=> $b->[0]} @{$AoALTs2};
		
		my @sortedALTs;
		my @ACs;
		my @AFs;
		my @LAs;
		foreach my $tmp (@{$AoALTs2}){
			push(@ACs, $tmp->[1]);
			push(@AFs, $tmp->[2]);
			push(@sortedALTs, $tmp->[3]);
			push(@LAs, $tmp->[4]);
		}
		my $vcfALT = join(",", @sortedALTs);
	#	my $vcfALT = join(",", @LAs);
		my $vcfinfo = "AC=".join(",", @ACs)."\;"."AF=".join(",", @AFs)."\;AN=".$numvar;
		$rvcf .= $sid."\t".$pstart."\t".$vpid."\t".$seqREF."\t".$vcfALT."\t60\t.\t".$vcfinfo."\tGT\t".join("\t", @GTs)."\n";
		$num_VPID++;
		$cnt_shrink++;
	}
	else{
		@{$AoALTs} = sort {$a->[0] <=> $b->[0]} @{$AoALTs};
		
		my @sortedALTs;
		my @ACs;
		my @AFs;
		my @LAs;
		foreach my $tmp (@{$AoALTs}){
			push(@ACs, $tmp->[1]);
			push(@AFs, $tmp->[2]);
			push(@sortedALTs, $tmp->[3]);
			push(@LAs, $tmp->[4]);
		}
		my $vcfALT = join(",", @sortedALTs);
	#	my $vcfALT = join(",", @LAs);
		my $vcfinfo = "AC=".join(",", @ACs)."\;"."AF=".join(",", @AFs)."\;AN=".$numvar;
		$rvcf .= $sid."\t".$pstart."\t".$vpid."\t".$seqREF."\t".$vcfALT."\t60\t.\t".$vcfinfo."\tGT\t".join("\t", @GTs)."\n";
		$num_VPID++;
	}
}

if($cnt_shrink > 0){
	print "! thread [$j] : [$sid] -> $num_VPID SVs | $cnt_revpos revised, $cnt_shrink pruned\n";
}
else{
	print "! thread [$j] : [$sid] -> $num_VPID SVs | $cnt_revpos revised\n";
}

#my $rtmpfile = "_".$sid."_tmp.$indel.vcf";
#open(my $rfh, ">", $rtmpfile);
#print $rfh $rvcf;
#close $rfh;

return $rvcf;
}


#-------------------------------------------------------------------------------
sub Read_target{
my $sample = shift;
my $file = shift;
my $htdata = shift;
my $tsid = shift;
my $dbseq = shift;

my $salt = "ACGTacgt";
my $cnt = 0;
if(-e $file){
	open(my $fh, "<", $file) or die;
	my $header = "";
	my $nrh = {};
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line){
			my @A = split(/\t/, $line);
			my $numA = @A;
			if($numA < 10 && $line !~ /CHROM/ && $line =~ /contig\=/ && $line =~ /\#/){
				my @elem = split(/contig\=\<ID\=/, $line);
				@elem = split(/,/, $elem[1]);
				if(! $nrh->{$elem[0]}){
					$header .= $line."\n";
					$nrh->{$elem[0]} = 1;
				}
			}
			elsif($numA >= 10 && $line !~ /CHROM/ && $A[4] !~ /\<IN/){
				if($A[0] ne $tsid){
					next;
				}
				
				my $RefAlt = $A[3]."_".$A[4];
				my $lenR = length($A[3]);
				my $lenA = length($A[4]);
				my $lendf = $lenR - $lenA;
				
				my $tmpREF = $A[3];
				my $refseq_check = substr($dbseq, $A[1] - 1, $lenR);
				if($refseq_check ne $A[3]){
					print "! unmatching refseq at $A[0] $A[1], skip...\n";
					next;
				}
				
				if($lenR > 200){
					$tmpREF = substr($A[3], 0, 200);
				}
				my $tmpALT = $A[4];
				if($lenR > 200){
					$tmpALT = substr($A[4], 0, 200);
				}
				$tmpREF = crypt($tmpREF, '$1$' . $salt);
				$tmpALT = crypt($tmpALT, '$1$' . $salt);
				$tmpREF =~ s/$salt//;
				$tmpALT =~ s/$salt//;
				$tmpREF =~ s/\t//g;
				$tmpALT =~ s/\t//g;
				$tmpREF =~ s/\n//g;
				$tmpALT =~ s/\n//g;
				
				my $varid = $A[0]."_".$A[1]."_".$lenR.$tmpREF."_".$lenA.$tmpALT;
				$htdata->{data}{$A[0]}{$sample}{$varid} = $line;
				
				my $indel = "null";
				if($lenR <= 1 && $lendf < 0){
					$indel = "INS";
				}
				elsif($lendf > 0){
					$indel = "DEL";
				}
				
				if($indel eq 'DEL'){
					my $pstart = $A[1];
					my $pend = $A[1] + $lenR - 1;
					for(my $i = $pstart; $i < $pend; $i++){
						$htdata->{varmap}{$indel}{$A[0]}{$i}{varid} .= $sample."\t".$varid."\n";
						
						if(! defined $htdata->{varmap}{$indel}{$A[0]}{$i}{pos0} && ! defined $htdata->{varmap}{$indel}{$A[0]}{$i}{pos1}){
							$htdata->{varmap}{$indel}{$A[0]}{$i}{pos0} = $pstart;
							$htdata->{varmap}{$indel}{$A[0]}{$i}{pos1} = $pend;
						}
						
						if(defined $htdata->{varmap}{$indel}{$A[0]}{$i}{pos0} && $pstart < $htdata->{varmap}{$indel}{$A[0]}{$i}{pos0}){
							$htdata->{varmap}{$indel}{$A[0]}{$i}{pos0} = $pstart;
						}
						if(defined $htdata->{varmap}{$indel}{$A[0]}{$i}{pos1} && $pend > $htdata->{varmap}{$indel}{$A[0]}{$i}{pos1}){
							$htdata->{varmap}{$indel}{$A[0]}{$i}{pos1} = $pend;
						}
					}
					$cnt++;
				}
				elsif($indel eq 'INS'){
					$htdata->{varmap}{$indel}{$A[0]}{$A[1]}{varid} .= $sample."\t".$varid."\n";
					$htdata->{varmap}{$indel}{$A[0]}{$A[1]}{pos0} = $A[1];
					$htdata->{varmap}{$indel}{$A[0]}{$A[1]}{pos1} = $A[1];
					$cnt++;
				}
			}
		}
	}
	close $fh;
	print " - [$sample] -> [$cnt] variants\n";
	
	$htdata->{header} = $header;
}

return $htdata;
}


#-------------------------------------------------------------------------------
sub Read_list{
my $file = shift;

my $hash = {};
my @ID;
my $err = 0;
if(-e $file){
	print "\n! reading list [$file]...\n";
	open(my $fh, "<", $file) or die;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line){
			my @A = split(/\t/, $line);
			if($A[0] && ! $hash->{$A[0]}){
				if($A[1] && -e $A[1]){
					print " [$A[0]] = [$A[1]] <OK>\n";
					$hash->{$A[0]} = $A[1];
					push(@ID, $A[0]);
				}
				elsif(! $A[1]){
					print " [$A[0]] = missing path for [$A[1]] (skip)\n";
					$err++;
				}
				elsif(! -e $A[1]){
					print " [$A[0]] = [$A[1]] missing file (skip)\n";
					$err++;
				}
			}
		}
	}
	close $fh;
}

return ($hash, \@ID, $err);
}


#-------------------------------------------------------------------------------
sub Open_fasta_as_hash{
my $file = shift;
my $tsid = shift;

#print "! open [$file] as hash ...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
my $cnt = 0;
my $ID;
my $seq;
my $total_len = 0;
my @SID;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	my @A = split(/\t/, $line);
	if($line =~ /\>/){
		unless($cnt == 0){
			if($hash->{$ID}){
				print "! duplicated seq ID : [$ID]\n";
			}
			
			if($ID eq $tsid){
				my $len = length($seq);
				$hash->{$ID}{seq} = $seq;
				$hash->{$ID}{len} = $len;
				$total_len += $len;
				push(@SID, $ID);
			}
		}
		
		$ID = $line;
		$ID =~ s/\>//;
		my @tmp = split(/\s/, $ID);
		$ID = $tmp[0];
		$seq = "";
		$cnt++;
	}
	else{
		$line =~ s/\.//g;
		if($line){
			$seq .= $line;
		}
	}
}
close $fh;

if($seq){
	if($hash->{$ID}){
		print "! duplicated seq ID : [$ID]\n";
	}
	
	if($ID eq $tsid){
		my $len = length($seq);
		$hash->{$ID}{seq} = $seq;
		$hash->{$ID}{len} = $len;
		$total_len += $len;
		push(@SID, $ID);
	}
}

my $numID = @SID;

print "! reference = [$file]\n";
print " - [$numID] sequence ID, [$total_len] bp\n";

return ($hash, \@SID);
}


#----------------------------------------------------------
sub ADD{
my $file = shift;
my $str = shift;

if($str){
	open(my $fh, ">>", $file) or die;
	print $fh $str;
	close $fh;
	
	#print "! output [$file]\n";
}

}


#----------------------------------------------------------
sub SAVE{
my $file = shift;
my $str = shift;

if($str){
	open(my $fh, ">", $file) or die;
	print $fh $str;
	close $fh;
	print "! output [$file]\n";
}

}


#-------------------------------------------------------------------------------
sub Path2pref{
my $file = shift;

if($file && $file ne 'null'){
	my @tmp = split(/\//, $file);
	my $n = @tmp;
	$tmp[$n - 1] =~ s/\.gvcf\.gz//;
	$tmp[$n - 1] =~ s/\.gvcf//;
	$tmp[$n - 1] =~ s/\.vcf\.gz//;
	$tmp[$n - 1] =~ s/\.vcf//;
	return $tmp[$n - 1];
}

}





