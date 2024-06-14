#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use FindBin;

my $dbfasta = shift;
my $qfasta = shift;
my $asm2sv_tsv = shift;
my $cpu = shift;
my $rdir = shift;
my $rvcf_candidate = shift;
my $rvcf_final = shift;
my $sample_uniqueID = shift;
my $chrinfo = shift;

my $neighbor_bp = 5000;
#$cpu = 24;
my $resume_from_native = 0;		# to resume process in an emergency case

my $script_path = $FindBin::Bin;
my $wpath = getcwd();
#my $HOME = File::HomeDir->my_home;
my $bin_aln2var = $script_path."/aln2var2.pl";
my $bin_sub = $script_path."/rev2sv_sub.pl";

my $err = 0;
$err += FileCheck2($dbfasta, "dbfasta");
$err += FileCheck2($qfasta, "qfasta");
$err += FileCheck2($asm2sv_tsv, "Asm2sv rfile_final");
$err += FileCheck2($chrinfo, "seqname alias between genomes");
if(! $rdir){
	print "! missing output directory...\n";
	$err++;
}
if(! $bin_aln2var || ! -e $bin_aln2var){
	print "! missing $bin_aln2var ...\n";
	$err++;
}
if($err > 0){
	print "! abort script...\n";
	goto END;
}

unless(-e $rdir){
	system("mkdir $rdir");
}
chdir $rdir;

LnFile($dbfasta);
LnFile($qfasta);
LnFile($asm2sv_tsv);
LnFile($chrinfo);

my $dbprefix = Path2pref($dbfasta);
my $qprefix = Path2pref($qfasta);
if(! $sample_uniqueID || $sample_uniqueID eq 'null'){
	$sample_uniqueID = $dbprefix;
}

my $log_dir = "log";
my $log_ess0 = "./log/log_rev2sv.txt";
#if(-d $log_dir){
#	print "! once delete [$log_dir]\n";
#	system("rm -R $log_dir");
#}
unless(-d $log_dir){
	system("mkdir $log_dir");
}

my $htmp_chrinfo = Read_seqinfo($chrinfo);
if(! $htmp_chrinfo->{seqID}{$dbprefix}){
	print "! error: no seqname alias for [$dbprefix]...\n";
	die;
}
if(! $htmp_chrinfo->{seqID}{$qprefix}){
	print "! error: no seqname alias for [$qprefix]...\n";
	die;
}

my $logfile_combined = "log_analysis_hist.txt";
my $combined_tmpheader = "tmp_combined_header.vcf";
my ($list, $CMD, $cnt_entry) = Read_Asm2sv_summary($asm2sv_tsv, $rvcf_candidate, $combined_tmpheader, $bin_aln2var, $resume_from_native, $logfile_combined, $htmp_chrinfo->{seqID}, $dbprefix, $qprefix);
if($cnt_entry == 0){
	print "! no entry found, skip...\n";
	goto Sort;
}

if($resume_from_native){
#	if(-e $rvcf_candidate){
#		print "! once delete [$rvcf_candidate]\n";
#		system("rm $rvcf_candidate");
#	}
	if(-e $rvcf_final){
		print "! once delete [$rvcf_final]\n";
		system("rm $rvcf_final");
	}
	goto Sort;
}

my $d = int($cnt_entry / $cpu) + 1;
my $thrs = [];
for(my $i = 0; $i < $cpu; $i++){
	my $p0 = $i * $d;
	my $p1 = ($i + 1) * $d;
	
	my $each_cmd = "";
	for(my $j = $p0; $j < $p1; $j++){
		if($CMD->[$j]){
			$each_cmd .= $CMD->[$j]."\n";
		}
	}
	
	if($each_cmd){
		my $sh = "./log/_batch_aln2sv_".$i.".sh";
		my $logfile = "./log/_log_aln2sv_".$i.".txt";
		SAVE($sh, $each_cmd);
		my $thr = threads->new(\&Batch_exec, $sh, $logfile, $i);
		push(@{$thrs}, $thr);
	}
}
if(@{$thrs}){
	foreach my $thr (@{$thrs}){
		$thr->join();
	}
}

Sort:{
	my $sortvcf = 1;
}

print "! sorting VCF with [$cpu] threads...\n";
print " - neighbor span = [$neighbor_bp] bp\n";

my ($hdbseq, $SID) = Open_fasta_as_hash($dbfasta);
my $num_SID = @{$SID};

my @RFs;
$d = int($num_SID / $cpu) + 1;
for(my $i = 0; $i < $d; $i++){
	my $p0 = $i * $cpu;
	my $p1 = ($i + 1) * $cpu;
	
	my $thrs = [];
	for(my $j = $p0; $j < $p1; $j++){
		if($SID->[$j]){
			my $each_rvcf = "./log/_each_integrated_".$SID->[$j].".vcf";
			my $logfile = "./log/_log_rev2sv_".$j.".txt";
			my $each_cmd = "$bin_sub $each_rvcf $rvcf_candidate $rvcf_final $combined_tmpheader $dbfasta $SID->[$j] $sample_uniqueID $asm2sv_tsv $resume_from_native $neighbor_bp $chrinfo $dbprefix $qprefix $j > $logfile 2>&1";
			
			if($resume_from_native && -e $each_rvcf){
				system("rm $each_rvcf");
			}
#			if(! -e $each_rvcf){
				my $sh = "./log/_batch_rev2sv_".$j.".sh";
				SAVE($sh, $each_cmd);
				push(@RFs, $each_rvcf);
				my $thr = threads->new(\&Batch_exec2, $sh, $logfile, $j, $SID->[$j]);
				push(@{$thrs}, $thr);
#			}
#			else{
#				print " thread [$j] : [$SID->[$j] already done, skip...\n";
#				push(@RFs, $each_rvcf);
#			}
		}
	}
	if(@{$thrs}){
		foreach my $thr (@{$thrs}){
			$thr->join();
		}
	}
}

if(-e $rvcf_final){
	system("rm $rvcf_final");
}
foreach my $rvcf (@RFs){
	CollectSAVE($rvcf_final, $rvcf);
}
print "! final result = [$rvcf_final]\n";

END:{
	my $end = "";
}


################################################################################
#-------------------------------------------------------------------------------
sub Read_seqinfo{
my $file = shift;

print "! reading seq ID alias info [$file]...\n";
open(my $fh, "<", $file);
my $hash = {};
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	$line =~ s/\"//g;
	
	if($line =~ /\#/){
		next;
	}
	
	# example of chrinfo (tsv)
	#	prefix		convert_ID	original_ID		alias_prefix		norder
	#	Gmax_189	chr01		Gm01			Williams82_189		10
	#	Gmax_189	chr02		Gm02			Williams82_189		10
	#	Gmax_189	chr03		Gm03			Williams82_189		10
	
	my @A = split(/\t/, $line);
	my @B = split(/,/, $A[2]);
	$hash->{seqID}{$A[0]}{$A[2]} = $A[1];
	$hash->{seqID}{$A[0]}{$B[0]} = $A[1];
	
	if(defined $A[3] && ! $hash->{alias_name}{$A[0]} && $A[3] ne '-' && $A[3] ne 'NA' && $A[3] ne 'na' && $A[3] ne '0'){
		$hash->{alias_name}{$A[0]} = $A[3];
	}
	if(defined $A[4] && ! $hash->{norder}{$A[0]} && $A[4] ne '-' && $A[4] ne 'NA' && $A[4] ne 'na' && $A[4] ne '0'){
		$hash->{norder}{$A[0]} = $A[4];
		$hash->{norder_status} = 1;
	}
#	print "$A[0] $B[0] $A[1]\n";
	$cnt++;
}
close $fh;

print "! [$cnt] lines\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub CollectSAVE{
my ($rfile, $file) = @_;

if(-e $file){
	open(my $fh, "<", $file) or die;
	my $r = "";
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line){
			$r .= $line."\n";
		}
	}
	close $fh;
	
	open(my $rfh, ">>", $rfile);
	print $rfh $r;
	close $rfh;
}

}


#-------------------------------------------------------------------------------
sub Batch_exec{
my ($sh, $logfile, $i) = @_;

print " thread [$i] : start...\n";
if(system("bash $sh > $logfile 2>&1") != 0){
	print "! thread [$i] : failed.\n";
}
else{
	print "! thread [$i] : completed.\n";
}

}


#-------------------------------------------------------------------------------
sub Batch_exec2{
my ($sh, $logfile, $i, $sid) = @_;

print " thread [$i] : [$sid] start...\n";
if(system("bash $sh > $logfile 2>&1") != 0){
	print "! thread [$i] : [$sid] failed.\n";
}
else{
	print "! thread [$i] : [$sid] completed.\n";
}

}


#-------------------------------------------------------------------------------
sub Read_Asm2sv_summary{
my $file = shift;
my $rvcf_final = shift;
my $tmpheader_vcf = shift;
my $bin_aln2var = shift;
my $resume_from_native = shift;
my $logfile_combined = shift;
my $hseqinfo = shift;
my $dbprefix = shift;
my $qprefix = shift;

my $hprev = {};
if($logfile_combined && -e $logfile_combined){
	open(my $lfh, "<", $logfile_combined);
	while(my $line = <$lfh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line){
			my @A = split(/\t/, $line);
			$hprev->{$A[0]} = 1;
		}
	}
	close $lfh;
}

print "! reading [$file]...\n";
open(my $fh, "<", $file) or die;
my $list = "";
my $nrcmds = {};
my $cnt_all = 0;
my $cnt_entry = 0;
my $cnt_skip = 0;
my $cnt_prev = 0;
my $cnt_nonnumeric = 0;
my $cnt_missingseq = 0;
my $cnt_unmatchseq = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ /directory/ && $line =~ /dbfasta/){
		next;
	}
	else{
		my @A = split(/\t/, $line);
		my $gid = $A[23];
		my $rsid = $A[2];
		my $rpos0 = $A[3];
		my $rpos1 = $A[4];
		my $qsid = $A[8];
		my $qpos0 = $A[9];
		my $qpos1 = $A[10];
		my $normality_score = $A[21];
		my $asm2sv_judge = $A[22];
		$cnt_all++;
		
		if($rpos0 =~ /\D/ || $rpos1 =~ /\D/ || $qpos0 =~ /\D/ || $qpos1 =~ /\D/){
			$cnt_nonnumeric++;
			next;
		}
		
		if(! $hseqinfo->{$dbprefix}{$rsid}){
			print "! error: missing seqname alias for DB [$dbprefix] [$rsid]...\n";
			print "$file\n$line\n";
			die;
		}
		if(! $hseqinfo->{$qprefix}{$qsid}){
			print "! error: missing seqname alias for q [$qprefix] [$qsid]...\n";
			print "$file\n$line\n";
			die;
		}
		my $rsid_alias = $hseqinfo->{$dbprefix}{$rsid};
		my $qsid_alias = $hseqinfo->{$qprefix}{$qsid};
		
		if($rsid_alias ne $qsid_alias){
			$cnt_unmatchseq++;
			next;
		}
		if($normality_score && $normality_score =~ /\d/ && $normality_score > 0.99 && $normality_score < 1.01 && $asm2sv_judge && $asm2sv_judge eq 'present'){
			$cnt_skip++;
			next;
		}
		if($hprev->{$gid}){
			$cnt_prev++;
			next;
		}
		
		if($rpos1 < $rpos0){
			$rpos0 = $A[4];
			$rpos1 = $A[3];
		}
		if($qpos1 < $qpos0){
			$qpos0 = $A[4];
			$qpos1 = $A[3];
		}
		
		$list .= $gid."\t".$rsid."\t".$rpos0."\t".$rpos1."\t".$qsid."\t".$qpos0."\t".$qpos1."\n";
		my $each_cmd .= "$bin_aln2var $gid $dbfasta $rsid $rpos0 $rpos1 $qfasta $qsid $qpos0 $qpos1 $rvcf_final $tmpheader_vcf";
		if($resume_from_native){
			$each_cmd .= " resume_from_native";
		}
		$nrcmds->{$each_cmd} = 1;
		$cnt_entry++;
	}
}
close $fh;

my @CMD = keys(%{$nrcmds});

print " - [$cnt_all] entries, of these\n";
print " - [$cnt_prev] previously analyzed (skip)\n";
print " - [$cnt_unmatchseq] not identical Chr group (skip)\n";
print " - [$cnt_skip] almost same with reference genome according to Asm2sv result (skip)\n";
print " - [$cnt_entry] are analysis targets here\n";

if($cnt_nonnumeric > 0){
	print " - [$cnt_nonnumeric] with non-numeric position (skip)\n";
}
if($cnt_missingseq > 0){
	print " - [$cnt_missingseq] with missing seq (skip)\n";
}

return ($list, \@CMD, $cnt_entry);
}


#-------------------------------------------------------------------------------
sub Open_fasta_as_hash{
my $file = shift;

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
#				print "! duplicated seq ID : [$ID]\n";
			}
			
			my $len = length($seq);
			$hash->{$ID}{seq} = $seq;
			$hash->{$ID}{len} = $len;
			$hash->{$ID}{n} = $cnt;
			$total_len += $len;
			push(@SID, $ID);
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
#		print "! duplicated seq ID : [$ID]\n";
	}
	
	my $len = length($seq);
	$hash->{$ID}{seq} = $seq;
	$hash->{$ID}{len} = $len;
	$hash->{$ID}{n} = $cnt;
	$total_len += $len;
	push(@SID, $ID);
}

my $numID = @SID;
#print "! [$numID] sequences with [$total_len] bp\n";

return ($hash, \@SID);
}


#-------------------------------------------------------------------------------
sub Rmfiles{
my $dir = shift;

my @GB0 = glob("$dir/*");
foreach my $file (@GB0){
	if($file =~ /BDconfidenthit/ || $file =~ /bN_db-/ || $file =~ /bP_db-/){
		next;
	}
	system("rm $file");
}

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


#-------------------------------------------------------------------------------
sub LnFile{
my $file = shift;

if($file && $file ne 'null' && ! -e $file){
	system("ln -s ../$file ./");
}

}


#-------------------------------------------------------------------------------
sub FileCheck1{
my $file = shift;
my $str = shift;

my $err = 0;
if(! $file){
	print "! missing $str ...\n";
	$err++;
}
elsif(! -e $file){
	print "! missing [$file] ...\n";
	$err++;
}

return $err;
}


#-------------------------------------------------------------------------------
sub FileCheck2{
my $file = shift;
my $str = shift;

my $err = 0;
if(! $file){
	print "! missing $str ...\n";
	$err++;
}
elsif(! -e $file){
	print "! missing [$file] ...\n";
	$err++;
}
elsif($file =~ /\//){
	print "! [$file] must not be PATH...\n";
	$err++;
}

return $err;
}

