#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use FindBin;

my $gid = shift;
my $dbfasta = shift;
my $rsid = shift;
my $rpos0 = shift;
my $rpos1 = shift;
my $qfasta = shift;
my $qsid = shift;
my $qpos0 = shift;
my $qpos1 = shift;
my $rfile_combined = shift;
my $rfile_tmpheader = shift;
my $resume_from_native = shift;

my $logfile_combined = "log_analysis_hist.txt";
my $native_combined = "native_sv.vcf";
$resume_from_native = 0;

if(! $gid || ! $rfile_combined){
	goto END;
}

my $dbpref = Path2pref($dbfasta);
my $qpref = Path2pref($qfasta);

my $sub_dbfasta = "_".$dbpref."_".$rsid."_".$rpos0."-".$rpos1.".fasta";
my $sub_qfasta = "_".$qpref."_".$qsid."_".$qpos0."-".$qpos1.".fasta";
my $rev_dbfasta = "_rev_".$dbpref."_".$rsid."_".$rpos0."-".$rpos1.".fasta";
my $rev_dbprefx = "_rev_".$dbpref."_".$rsid."_".$rpos0."-".$rpos1;

my $sample = "_db-".$dbpref."_".$rsid."_".$rpos0."-".$rpos1."_q-".$qpref."_".$qsid."_".$qpos0."-".$qpos1;
my $maf = $sample.".maf";
my $psl = $sample.".psl";
my $sam = $sample.".sam";
my $bam = $sample.".bam";
my $bai = $sample.".bam.bai";
my $vcf = $sample.".vcf";
my $tmpvcf = $sample.".tmpvcf";

my $tmpdir = "/media/nvmer/analysis/_tmp_Bt2sub_FindCDSmut";
if(-e $tmpdir){
	$sub_dbfasta = "$tmpdir/"."_".$dbpref."_".$rsid."_".$rpos0."-".$rpos1.".fasta";
	$sub_qfasta = "$tmpdir/"."_".$qpref."_".$qsid."_".$qpos0."-".$qpos1.".fasta";
	$rev_dbfasta = "$tmpdir/"."_rev_".$dbpref."_".$rsid."_".$rpos0."-".$rpos1.".fasta";
	$maf = "$tmpdir/".$sample.".maf";
	$psl = "$tmpdir/".$sample.".psl";
	$sam = "$tmpdir/".$sample.".sam";
	$bam = "$tmpdir/".$sample.".bam";
	$bai = "$tmpdir/".$sample.".bam.bai";
	$vcf = "$tmpdir/".$sample.".vcf";
	$tmpvcf = "$tmpdir/".$sample.".tmpvcf";
}

print "\n---------------------------------------------------------------------------------\n";
print "! gene = [$gid]\n";
print "! db = [$dbpref] [$rsid] [$rpos0 - $rpos1]\n";
print "! q  = [$qpref] [$qsid] [$qpos0 - $qpos1]\n";

my $prev_analysis = SearchLog($logfile_combined, $gid, $dbfasta, $rsid, $rpos0, $rpos1, $qfasta, $qsid, $qpos0, $qpos1);
if($prev_analysis eq 'true'){
	if($resume_from_native && -e $native_combined){
		print "! resume mode\n";
		my $hdbseq_tmp = Open_fasta_as_hash($dbfasta, $rsid);
		AddVcf_fromnative($native_combined, $rfile_combined, $rfile_tmpheader, $hdbseq_tmp->{$rsid}{seq}, $rsid, $rpos0, $rpos1);
	}
	else{
		print "! analysis already done, skip...\n";
	}
	goto END;
}

if($resume_from_native && -e $native_combined){
	goto Skip;
}

my $hdbseq = Open_fasta_as_hash($dbfasta, $rsid);
my $hqseq = Open_fasta_as_hash($qfasta, $qsid);

my $rseq = Prepfasta($dbpref, $hdbseq->{$rsid}{seq}, $rsid, $rpos0, $rpos1, "db", 1);
my $qseq = Prepfasta($qpref, $hqseq->{$qsid}{seq}, $qsid, $qpos0, $qpos1, "q", 100);

SAVE($sub_dbfasta, $rseq);
SAVE($sub_qfasta, $qseq);

if(-e $sam){
	system("rm $sam");
}
if(-e $vcf){
	system("rm $vcf");
}

my $program = "minimap2";
my $judge_skip = 0;

if($program eq 'both' || $program eq 'last'){
	my $bin_lastdb = "lastdb";
	my $bin_lastal = "lastal";
	my $bin_mafconvert = "maf-convert";
	
	LASTINDEX($sub_dbfasta, $bin_lastdb);
	
	my $lastal_e = 25;
	my $lastal_q = 3;
	my $lastal_j = 4;
	my $lastal_a = 1;
	my $lastal_b = 1;
	my $lastsp_s = 35;
	
#	my $cmd3 = "$bin_lastal -e $lastal_e -q $lastal_q -j $lastal_j -P 1 -a $lastal_a -b $lastal_b $sub_dbfasta $qfasta | last-split -s $lastsp_s > $maf";
	my $cmd3 = "$bin_lastal -P 1 $sub_dbfasta $sub_qfasta | last-split -s $lastsp_s > $maf";
	print "! cmd = [$cmd3]\n";
	if(system("timeout 300 $cmd3") != 0){
		print "! lastal failed, skip analysis...\n";
		$judge_skip = 1;
		goto Skip;
	}
	
	my $cntfail_mafconv = 0;
	RedoMafconv:{
		my $RedoMafconv = 1;
	}
	
#	my $cmd4 = "$bin_mafconvert psl $maf > $psl";
#	print "! cmd = [$cmd4]\n";
#	if(system("$cmd4") != 0){
#		$cntfail_mafconv++;
#		if($cntfail_mafconv <= 3){
#			goto RedoMafconv;
#		}
#	}
	
	if($program eq 'both'){
		system("$bin_mafconvert psl $maf > $psl");
		Revfasta($psl, $sub_dbfasta, $rseq, $rev_dbfasta, $rev_dbprefx);
	}
	else{
		my $cmd5 = "$bin_mafconvert sam $maf >> $sam";
		print "! cmd = [$cmd5]\n";
		if(system("$cmd5") != 0){
			$cntfail_mafconv++;
			if($cntfail_mafconv <= 3){
				goto RedoMafconv;
			}
		}
		AddHeaderSAM($sam, $rsid, $hdbseq->{$rsid}{len});
	}
}
if($program eq 'both' || $program eq 'minimap2'){
	my $bin_minimap2 = "minimap2";
	my $min2ndratio_th = 0.99;
	
	my $cmd3 = "";
	if($program eq 'both'){
		$cmd3 = "$bin_minimap2 -ax map-hifi --MD -p $min2ndratio_th -t 2 $rev_dbfasta $sub_qfasta > $sam";
		if(! -e $rev_dbfasta || ! -e $sub_qfasta){
			print "! missing required file(s)...\n";
			goto Skip;
		}
	}
	else{
		$cmd3 = "$bin_minimap2 -ax map-hifi --MD -p $min2ndratio_th -t 2 $sub_dbfasta $sub_qfasta > $sam";
		if(! -e $sub_dbfasta || ! -e $sub_qfasta){
			print "! missing required file(s)...\n";
			goto Skip;
		}
	}
	
	print "! cmd = [$cmd3]\n";
	if(system("timeout 300 $cmd3") != 0){
		AddLog($logfile_combined, $gid, $dbfasta, $rsid, $rpos0, $rpos1, $qfasta, $qsid, $qpos0, $qpos1, "timeup");
		print "! minimap2 timeup, skip analysis...\n";
		$judge_skip = 1;
		goto Skip;
	}
}

if(-e $sam){
	my $bin_samtools = "samtools";
	my $cmd6 = "$bin_samtools view -bS $sam > $bam -@ 1";
	my $cmd7 = "$bin_samtools sort -T $sample -@ 1 -o $bam $bam";
	my $cmd8 = "$bin_samtools index $bam";
	
	my $bin_sniffles = "sniffles";
	my $cmd9 = "$bin_sniffles --input $bam --vcf $vcf";
#	my $cmd9 = "$bin_sniffles --input $bam --vcf $tmpvcf";
#	my $cmd10 = "$bin_sniffles --input $bam --genotype-vcf $tmpvcf --vcf $vcf";
	
	system($cmd6);
	system($cmd7);
	system($cmd8);
	system($cmd9);
#	system($cmd10);
}
if(-e $vcf){
	AddLog($logfile_combined, $gid, $dbfasta, $rsid, $rpos0, $rpos1, $qfasta, $qsid, $qpos0, $qpos1, "analyzed");
	AddVcf($native_combined, $rfile_combined, $rfile_tmpheader, $vcf, $hdbseq->{$rsid}{seq}, $rpos0, $rpos1);
}

Skip:{
	if($judge_skip){
		print "\n! [$gid] | process abondoned due to timeup (this entry should not be used for pangenome construction)...\n";
	}
}

Delete("$sub_dbfasta");
Delete("$sub_qfasta");
Delete("$sub_dbfasta.bck");
Delete("$sub_dbfasta.des");
Delete("$sub_dbfasta.prj");
Delete("$sub_dbfasta.sds");
Delete("$sub_dbfasta.ssp");
Delete("$sub_dbfasta.suf");
Delete("$sub_dbfasta.tis");
Delete("$sam");
Delete("$bam");
Delete("$bai");
Delete("$vcf");
Delete("$tmpvcf");
Delete("$psl");
Delete("$maf");

END:{
	print "\n";
}


################################################################################
#-------------------------------------------------------------------------------
sub Revfasta{
my ($psl, $sub_dbfasta, $dbseq, $rev_dbfasta, $rev_dbprefx) = @_;

my @SL = split(/\n/, $dbseq);
my $rev_dbseq = "";
foreach my $line (@SL){
	if($line =~ /\>/){
		next;
	}
	$rev_dbseq .= $line;
}
$dbseq = $rev_dbseq;

print "! reading [$psl]...\n";
open(my $fh, "<", $psl) or die;
my $hash = {};
my $nrh = {};
my $alenF = 0;
my $alenR = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line && ! $nrh->{$line}){
		$nrh->{$line} = 1;
		
		my @A = split(/\t/, $line);
		my $numA = @A;
		if($A[0] >= 1000){
			if($A[8] eq '+'){
				$alenF += $A[0];
				$hash->{F}{pos0} .= $A[15]."\n";
				$hash->{F}{pos1} .= $A[16]."\n";
				$hash->{F}{cnt} += 1;
			}
			elsif($A[8] eq '-'){
				$alenR += $A[0];
				$hash->{R}{pos0} .= $A[15]."\n";
				$hash->{R}{pos1} .= $A[16]."\n";
				$hash->{R}{cnt} += 1;
			}
		}
	}
}
close $fh;

my @P0;
my @P1;
if($alenF > $alenR){
	@P0 = split(/\n/, $hash->{F}{pos0});
	@P1 = split(/\n/, $hash->{F}{pos1});
	print "! plus strand is predominant\n";
}
else{
	@P0 = split(/\n/, $hash->{R}{pos0});
	@P1 = split(/\n/, $hash->{R}{pos1});
	print "! minus strand is predominant\n";
}

my $numP = @P0;
my $AoP = [];
for(my $i = 0; $i < $numP; $i++){
	my @tmp;
	push(@tmp, $P0[$i]);
	push(@tmp, $P1[$i]);
	push(@{$AoP}, \@tmp);
}
@{$AoP} = sort {$a->[0] <=> $b->[0]} @{$AoP};

my $revseq = ">$rev_dbprefx\n";
my $prev0 = -1;
for(my $i = 0; $i < $numP; $i++){
	if($i == 0){
		$revseq .= substr($dbseq, 0, $AoP->[$i][0] - 1);
		$revseq .= substr($dbseq, $AoP->[$i][0] - 1, $AoP->[$i][1] - $AoP->[$i][0] + 1);
		$prev0 = $AoP->[$i][1];
	}
	else{
		my $numN = $AoP->[$i][0] - $prev0 - 1;
		for(my $j = 0; $j < $numN; $j++){
			$revseq .= "N";
		}
		$revseq .= substr($dbseq, $AoP->[$i][0] - 1, $AoP->[$i][1] - $AoP->[$i][0] + 1);
		$prev0 = $AoP->[$i][1];
	}
}
$revseq .= substr($dbseq, $prev0);

open(my $rfh, ">", $rev_dbfasta);
print $rfh $revseq;
close $rfh;
print "! output [$rev_dbfasta]\n";

}


#-------------------------------------------------------------------------------
sub AddVcf_fromnative{
my ($rfile_native, $rfile, $tmpheader, $dbseq, $rsid, $rpos0, $rpos1) = @_;

print "! reading [$rfile_native] to resume combine...\n";
open(my $fh, "<", $rfile_native) or die;
my $r = "";
my $header = "";
my $cnt_fix = 0;
my $cnt_unknown = 0;
my $cnt_borderaln = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line){
		my @A = split(/\t/, $line);
		my $numA = @A;
		if($numA < 10 || $A[0] =~ /\#CHROM/){
			$header .= $line."\n";
		}
		else{
			#	chr01	1206003	Sniffles2.INS.1S0	N	GGTATTGAACCTCATCATAATGGCATTCAACAATACTATGTAAGAACATACAATGTTTCGCCACTAAGAAGCAATTATCATATAAAATAATTATTATTAATTATAGGACTTATTCCAAACCAAACATGAAACTTGGAAGAAGAGTTTTACTATTTATAATGATAGTGTCAGCCAC	60	PASS	PRECISE;SVTYPE=INS;SVLEN=175;END=1206003;SUPPORT=40;COVERAGE=40,40,40,40,40;STRAND=+;AF=1.000;STDEV_LEN=0.000;STDEV_POS=0.000;SUPPORT_LONG=0	GT:GQ:DR:DV	1/1:60:0:40
			#	chr01	1213302	Sniffles2.DEL.3S0	N	<DEL>	60	PASS	PRECISE;SVTYPE=DEL;SVLEN=-3187;END=1216489;SUPPORT=40;COVERAGE=40,40,40,40,40;STRAND=+;AF=1.000;STDEV_LEN=0.000;STDEV_POS=0.000	GT:GQ:DR:DV	1/1:60:0:40
			#	chr01	1712919	Sniffles2.DUP.2S0	N	<DUP>	60	PASS	PRECISE;SVTYPE=DUP;SVLEN=198;END=1713117;SUPPORT=40;COVERAGE=40,40,40,40,40;STRAND=+;AF=0.571;STDEV_LEN=0.000;STDEV_POS=0.000	GT:GQ:DR:DV	0/1:60:30:40
			#	chr01	8230421	Sniffles2.INV.1S0	N	<INV>	60	PASS	PRECISE;SVTYPE=INV;SVLEN=18170;END=8248591;SUPPORT=40;COVERAGE=0,0,40,40,80;STRAND=+;AF=0.400;STDEV_LEN=0.000;STDEV_POS=0.000	GT:GQ:DR:DV	0/1:60:60:40
			
			if($rsid eq $A[0] && $A[3] eq 'N'){
				my @tag = split(/\;/, $A[7]);
				my $svtype = "null";
				my $pend = "null";
				foreach my $str (@tag){
					if($str =~ /SVTYPE\=/){
						$str =~ s/SVTYPE\=//;
						$svtype = $str;
					}
					elsif($str =~ /END\=/){
						$str =~ s/END\=//;
						$pend = $str;
					}
					if($svtype ne 'null' && $pend ne 'null'){
						last;
					}
				}
				
				my $pstart = $A[1];
				my $pfront = $A[1] - 1;
				
				if($rpos0 <= $pstart && $pstart <= $rpos1 && $rpos0 <= $pend && $pend <= $rpos1){
					my $judge_border_aln = "false";
					my $neighbor = 100;
					if($rpos0 - $neighbor <= $pstart && $pstart <= $rpos0 + $neighbor){
						$judge_border_aln = "true";
					}
					if($rpos1 - $neighbor <= $pstart && $pstart <= $rpos1 + $neighbor){
						$judge_border_aln = "true";
					}
					if($rpos0 - $neighbor <= $pend && $pend <= $rpos0 + $neighbor){
						$judge_border_aln = "true";
					}
					if($rpos1 - $neighbor <= $pend && $pend <= $rpos1 + $neighbor){
						$judge_border_aln = "true";
					}
					if($judge_border_aln eq 'true'){
						$cnt_borderaln++;
						next;
					}
					
					if($svtype eq 'DEL' && $pend ne 'null'){
						$A[1] -= 1;
						$A[3] = substr($dbseq, $pfront - 1, 1).substr($dbseq, $pstart - 1, $pend - $pstart + 1);
						$A[4] = substr($dbseq, $pfront - 1, 1);
						$r .= join("\t", @A)."\n";
						$cnt_fix++;
					}
					elsif($svtype eq 'INS' && $pend ne 'null'){
						$A[1] -= 1;
						$A[3] = substr($dbseq, $pfront - 1, 1);
						$A[4] = substr($dbseq, $pfront - 1, 1).$A[4];
						$r .= join("\t", @A)."\n";
						$cnt_fix++;
					}
					elsif($svtype eq 'INV' && $pend ne 'null'){
						$A[3] = substr($dbseq, $pstart - 1, $pend - $pstart + 1);
						$A[4] = $A[3];
						$A[4] =~ tr/acgtACGT/tgcaTGCA/;
						$A[4] = reverse($A[4]);
						$r .= join("\t", @A)."\n";
						$cnt_fix++;
					}
					elsif($svtype eq 'DUP' && $pend ne 'null'){
						$A[3] = substr($dbseq, $pstart - 1, $pend - $pstart + 1);
						$A[4] = $A[3].$A[3];
						$r .= join("\t", @A)."\n";
						$cnt_fix++;
					}
				}
			}
		}
	}
}
close $fh;

print "! [$cnt_borderaln] alignment at target border (skip)\n";
print "! [$cnt_fix] <DEL> or <INS>\n";
print "! [$cnt_unknown] other category\n";

if($r){
	open(my $rfh, ">>", $rfile);
	print $rfh $r;
	close $rfh;
}

}


#-------------------------------------------------------------------------------
sub AddVcf{
my ($rfile_native, $rfile, $tmpheader, $vcf, $dbseq, $rpos0, $rpos1) = @_;

print "! reading [$vcf] to combine...\n";
open(my $fh, "<", $vcf) or die;
my $r = "";
my $rn = "";
my $header = "";
my $cnt_fix = 0;
my $cnt_unknown = 0;
my $cnt_borderaln = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line){
		my @A = split(/\t/, $line);
		my $numA = @A;
		if($numA < 10 || $A[0] =~ /\#CHROM/){
			$header .= $line."\n";
		}
		else{
			#	chr01	1206003	Sniffles2.INS.1S0	N	GGTATTGAACCTCATCATAATGGCATTCAACAATACTATGTAAGAACATACAATGTTTCGCCACTAAGAAGCAATTATCATATAAAATAATTATTATTAATTATAGGACTTATTCCAAACCAAACATGAAACTTGGAAGAAGAGTTTTACTATTTATAATGATAGTGTCAGCCAC	60	PASS	PRECISE;SVTYPE=INS;SVLEN=175;END=1206003;SUPPORT=40;COVERAGE=40,40,40,40,40;STRAND=+;AF=1.000;STDEV_LEN=0.000;STDEV_POS=0.000;SUPPORT_LONG=0	GT:GQ:DR:DV	1/1:60:0:40
			#	chr01	1213302	Sniffles2.DEL.3S0	N	<DEL>	60	PASS	PRECISE;SVTYPE=DEL;SVLEN=-3187;END=1216489;SUPPORT=40;COVERAGE=40,40,40,40,40;STRAND=+;AF=1.000;STDEV_LEN=0.000;STDEV_POS=0.000	GT:GQ:DR:DV	1/1:60:0:40
			#	chr01	1712919	Sniffles2.DUP.2S0	N	<DUP>	60	PASS	PRECISE;SVTYPE=DUP;SVLEN=198;END=1713117;SUPPORT=40;COVERAGE=40,40,40,40,40;STRAND=+;AF=0.571;STDEV_LEN=0.000;STDEV_POS=0.000	GT:GQ:DR:DV	0/1:60:30:40
			#	chr01	8230421	Sniffles2.INV.1S0	N	<INV>	60	PASS	PRECISE;SVTYPE=INV;SVLEN=18170;END=8248591;SUPPORT=40;COVERAGE=0,0,40,40,80;STRAND=+;AF=0.400;STDEV_LEN=0.000;STDEV_POS=0.000	GT:GQ:DR:DV	0/1:60:60:40
			
			if($A[3] eq 'N'){
				my @tag = split(/\;/, $A[7]);
				my $svtype = "null";
				my $pend = "null";
				foreach my $str (@tag){
					if($str =~ /SVTYPE\=/){
						$str =~ s/SVTYPE\=//;
						$svtype = $str;
					}
					elsif($str =~ /END\=/){
						$str =~ s/END\=//;
						$pend = $str;
					}
					if($svtype ne 'null' && $pend ne 'null'){
						last;
					}
				}
				
				my $pstart = $A[1];
				my $pfront = $A[1] - 1;
				
				my $judge_border_aln = "false";
				my $neighbor = 100;
				if($rpos0 - $neighbor <= $pstart && $pstart <= $rpos0 + $neighbor){
					$judge_border_aln = "true";
				}
				if($rpos1 - $neighbor <= $pstart && $pstart <= $rpos1 + $neighbor){
					$judge_border_aln = "true";
				}
				if($rpos0 - $neighbor <= $pend && $pend <= $rpos0 + $neighbor){
					$judge_border_aln = "true";
				}
				if($rpos1 - $neighbor <= $pend && $pend <= $rpos1 + $neighbor){
					$judge_border_aln = "true";
				}
				if($judge_border_aln eq 'true'){
					$cnt_borderaln++;
					next;
				}
				
				if($svtype eq 'DEL' && $pend ne 'null'){
					$A[1] -= 1;
					$A[3] = substr($dbseq, $pfront - 1, 1).substr($dbseq, $pstart - 1, $pend - $pstart + 1);
					$A[4] = substr($dbseq, $pfront - 1, 1);
					$r .= join("\t", @A)."\n";
					$cnt_fix++;
				}
				elsif($svtype eq 'INS' && $pend ne 'null'){
					$A[1] -= 1;
					$A[3] = substr($dbseq, $pfront - 1, 1);
					$A[4] = substr($dbseq, $pfront - 1, 1).$A[4];
					$r .= join("\t", @A)."\n";
					$cnt_fix++;
				}
				elsif($svtype eq 'INV' && $pend ne 'null'){
					$A[3] = substr($dbseq, $pstart - 1, $pend - $pstart + 1);
					$A[4] = $A[3];
					$A[4] =~ tr/acgtACGT/tgcaTGCA/;
					$A[4] = reverse($A[4]);
					$r .= join("\t", @A)."\n";
					$cnt_fix++;
				}
				elsif($svtype eq 'DUP' && $pend ne 'null'){
					$A[3] = substr($dbseq, $pstart - 1, $pend - $pstart + 1);
					$A[4] = $A[3].$A[3];
					$r .= join("\t", @A)."\n";
					$cnt_fix++;
				}
	#			else{
	#				$r .= $line."\n";
	#				$cnt_unknown++;
	#			}
			}
			$rn .= $line."\n";
		}
	}
}
close $fh;

print "! [$cnt_borderaln] alignment at target border (skip)\n";
print "! [$cnt_fix] <DEL> or <INS>\n";
print "! [$cnt_unknown] other category\n";

if($header){
	open(my $rfh, ">>", $tmpheader);
	print $rfh $header;
	close $rfh;
}
if($r){
	open(my $rfh, ">>", $rfile);
	print $rfh $r;
	close $rfh;
}
if($rn){
	open(my $rfh, ">>", $rfile_native);
	print $rfh $rn;
	close $rfh;
}

}


#-------------------------------------------------------------------------------
sub SearchLog{
my ($file, $gid, $dbfasta, $rsid, $rpos0, $rpos1, $qfasta, $qsid, $qpos0, $qpos1) = @_;

my $judge = "false";
if(-e $file){
	open(my $fh, "<", $file) or die;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line){
			my @A = split(/\t/, $line);
			if($gid eq $A[0] && $dbfasta eq $A[1] && $rsid eq $A[2] && $rpos0 eq $A[3] && $rpos1 eq $A[4] && $qfasta eq $A[5] && $qsid eq $A[6] && $qpos0 eq $A[7] && $qpos1 eq $A[8]){
				$judge = "true";
				last;
			}
		}
	}
	close $fh;
}

return $judge;
}


#-------------------------------------------------------------------------------
sub AddLog{
my ($file, $gid, $dbfasta, $rsid, $rpos0, $rpos1, $qfasta, $qsid, $qpos0, $qpos1, $info) = @_;

my $str = "$gid\t$dbfasta\t$rsid\t$rpos0\t$rpos1\t$qfasta\t$qsid\t$qpos0\t$qpos1\t$info\n";
open(my $fh, ">>", $file);
print $fh $str;
close $fh;

}


#-------------------------------------------------------------------------------
sub AddHeaderSAM{
my ($file, $sid, $len) = @_;

open(my $fh, "<", $file) or die;
my $header = "";
my $r = "";
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line){
		my @A = split(/\t/, $line);
		my $numA = @A;
		if($numA < 4 && $A[0] =~ /\@/){
			$header .= $line."\n";
		}
		else{
			$r .= $line."\n";
		}
	}
}
close $fh;

$header .= "\@SQ\tSN:".$sid."\tLN:".$len."\n";
$r = $header.$r;

open(my $rfh, ">", $file);
print $rfh $r;
close $rfh;

}


#-------------------------------------------------------------------------------
sub Prepfasta{
my ($pref, $seq, $sid, $pos0, $pos1, $str, $multiply) = @_;

#my $fasta = ">".$str."-".$pref."_".$sid."_".$pos0."-".$pos1."\n";
my $len = length($seq);
print " - $pref -> $len bp\n";

my $fasta = "";
if($multiply == 1){
	my $revseq5 = "";
	if($pos0 > 1){
		for(my $i = 0; $i < $pos0; $i++){
			$revseq5 .= "N";
		}
	}
	
	my $revseq3 = "";
	if($pos1 < $len){
		for(my $i = $pos1; $i < $len; $i++){
			$revseq3 .= "N";
		}
	}
	
	$fasta = ">".$sid."\n".$revseq5.substr($seq , $pos0, $pos1 - $pos0 + 1).$revseq3."\n";
}
else{
	my $tmpseq = substr($seq , $pos0, $pos1 - $pos0 + 1);
	for(my $i = 0; $i < $multiply; $i++){
		$fasta .= ">".$sid."_".$i."\n".$tmpseq."\n";
	}
}

return $fasta;
}


#-------------------------------------------------------------------------------
sub Delete{
my $file = shift;

if($file && -e $file){
#	print "rm $file\n";
	system("rm $file");
}

}


#-------------------------------------------------------------------------------
sub LASTINDEX{
my $fasta = shift;
my $bin_lastdb = shift;

my $reference = $fasta;
#$reference =~ s/\.fasta//;
#$reference =~ s/\.fa//;

print "! creating lastal index for [$reference]\n";
system("$bin_lastdb $reference $fasta");

}


#-------------------------------------------------------------------------------
sub Open_fasta_as_hash{
my $file = shift;
my $target = shift;

#print "! reading [$file] as hash ...\n";
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
			if($ID eq $target){
				my $len = length($seq);
				$hash->{$ID}{seq} = $seq;
				$hash->{$ID}{len} = $len;
				$total_len += $len;
				push(@SID, $ID);
				$seq = "";
				last;
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
#		$line =~ s/\.//g;
#		if($line){
			$seq .= $line;
#		}
	}
}
close $fh;

if($seq){
	if($hash->{$ID}){
		print "! duplicated seq ID : [$ID]\n";
	}
	if($ID eq $target){
		my $len = length($seq);
		$hash->{$ID}{seq} = $seq;
		$hash->{$ID}{len} = $len;
		$total_len += $len;
		push(@SID, $ID);
	}
}

my $numID = @SID;

print " [$file] -> [$target] [$total_len] bp\n";

return $hash;
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
	
	#print "! output [$file]\n";
}

}


#-------------------------------------------------------------------------------
sub Path2pref{
my $file = shift;

if($file && $file ne 'null'){
	my @tmp = split(/\//, $file);
	my $n = @tmp;
	$tmp[$n - 1] =~ s/\.fasta//;
	$tmp[$n - 1] =~ s/\.fa//;
	return $tmp[$n - 1];
}

}



