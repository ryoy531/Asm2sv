#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "Blast_BDH.pl version 1.01\n";
$version .= "last update: [2019\/2\/10]\n";
$version .= "copyright: ryoichi yano [ryoichiy104\@gmail.com]\n";

#-------------------------------------------------------------------------------

#bidirectional blastを実行

my $rdir = shift;
my $qgenome1 = shift;
my $qgenome2 = shift;
my $qgff1 = shift;
my $qgff2 = shift;
my $qfasta1 = shift;
my $qfasta2 = shift;
my $program = shift;
my $maxhit = shift;
my $pcnt_th = shift;
my $pcov_th = shift;
my $eval_th = shift;
my $core = shift;
my $q = shift;
my $nohup = shift;
my $eaddr = shift;
my $out = shift;
my $timestamp = shift;

$nohup = "pipe";
$eaddr = "n";

unless(defined $rdir){
	goto END;
}
unless($nohup && $nohup eq 'pipe'){
	print "$version";
}

if($qfasta1){
	if($qfasta1 eq 't'){
		$program = 1;
		$maxhit = 20;
		$eval_th = 1;
		$core = 1;
		$nohup = "n";
		$eaddr = "n";
		$qfasta1 = "";
	}
	elsif($qfasta1 eq 'p'){
		$program = 2;
		$maxhit = 20;
		$eval_th = 1;
		$core = 1;
		$nohup = "n";
		$eaddr = "n";
		$qfasta1 = "";
	}
}

#--------------------------------------------------------------//
my $wpath = getcwd();
my $script_log = "";
unless($script_log){
	$script_log = "/media/hdd0/script_logs";
	unless(-e $script_log){
		$script_log = getcwd();
	}
}
#--------------------------------------------------------------//

unless($qfasta1){
	$qfasta1 = File_select(".fasta", "fasta1");
}

unless($qfasta2){
	$qfasta2 = File_select(".fasta", "fasta2");
}

unless($program){
	RE:{
		my $RE = 1;
	}
	
	print "\n---------------------------------------------------------------------\n";
	print "[1] blastn\n";
	print "[2] blastp\n";
	print "---------------------------------------------------------------------\n";
	print "Select program: ";
	$program = <STDIN>;
	$program =~ s/\n//;
	$program =~ s/\r//;
	unless($program){
		goto RE;
	}
	unless($program eq '1' || $program eq '2'){
		goto RE;
	}
	print "! [$program]\n";
}

unless($maxhit){
	print "\nEnter max hit for each query seq (default 99 = infinite): ";
	$maxhit = <STDIN>;
	$maxhit =~ s/\n//;
	$maxhit =~ s/\r//;
	unless($maxhit){
		$maxhit = 99;
	}
	unless($maxhit =~ /[0-9]/){
		if($maxhit =~ /[a-z]/i){
			print "! error: not numerous input [$maxhit]\n";
			goto END;
		}
	}
	if($maxhit < 1){
		$maxhit = 1;
	}
	print "! [$maxhit]\n";
}

my $default_pcnt_th;
my $default_pcov_th;
my $default_eval_th;
if($program eq '1'){
	$default_pcnt_th = 95;
	$default_pcov_th = 50;
	$default_eval_th = 1e-150;
}
elsif($program eq '2'){
	$default_pcnt_th = 50;
	$default_pcov_th = 50;
	$default_eval_th = 1e-150;
}

unless($pcnt_th){
	print "\nEnter threshold for \%identity (default $default_pcnt_th): ";
	$pcnt_th = <STDIN>;
	$pcnt_th =~ s/\n//;
	$pcnt_th =~ s/\r//;
	unless($pcnt_th){
		$pcnt_th = $default_pcnt_th;
	}
	unless($pcnt_th =~ /[0-9]/){
		if($pcnt_th =~ /[a-z]/i){
			print "! error: not numerous input [$pcnt_th]\n";
			goto END;
		}
	}
	if($pcnt_th < 1){
		$pcnt_th = 1;
	}
	print "! [$pcnt_th]\n";
}

unless($pcov_th){
	print "\nEnter threshold for \%alignment coverage (default $default_pcov_th): ";
	$pcov_th = <STDIN>;
	$pcov_th =~ s/\n//;
	$pcov_th =~ s/\r//;
	unless($pcov_th){
		$pcov_th = $default_pcov_th;
	}
	unless($pcov_th =~ /[0-9]/){
		if($pcov_th =~ /[a-z]/i){
			print "! error: not numerous input [$pcov_th]\n";
			goto END;
		}
	}
	if($pcov_th < 1){
		$pcov_th = 1;
	}
	print "! [$pcov_th]\n";
}

unless($eval_th){
	print "\nEnter threshold for E-value (default $default_eval_th): ";
	$eval_th = <STDIN>;
	$eval_th =~ s/\n//;
	$eval_th =~ s/\r//;
	unless($eval_th){
		$eval_th = $default_eval_th;
	}
	
	my $tmp_th = $eval_th;
	$tmp_th =~ s/e-//i;
	unless($tmp_th =~ /[0-9]/){
		if($tmp_th =~ /[a-z]/i){
			print "! error: not numerous input [$eval_th]\n";
			goto END;
		}
	}
	print "! [$eval_th]\n";
}

unless($core){
	print "\nEnter CPU threads used for each analysis (default 8): ";
	$core = <STDIN>;
	$core =~ s/\n//;
	$core =~ s/\r//;
	unless($core){
		$core = 8;
	}
	unless($core =~ /[0-9]/){
		if($core =~ /[a-z]/i){
			print "! error: not numerous input [$core]\n";
			goto END;
		}
	}
	print "! [$core]\n";
}

my $bin;
my $str;
if($program eq '1'){
	$bin = "blastn";
	$str = "bN";
}
elsif($program eq '2'){
	$bin = "blastp";
	$str = "bP";
}

my $analysis_log = "";
$analysis_log .= "\n------------------------------------------------------------------------------------\n";
$analysis_log .= "fasta1                   [$qfasta1]\n";
$analysis_log .= "fasta2                   [$qfasta2]\n";
$analysis_log .= "program                  [$bin]\n";
#$analysis_log .= "max hit per query       [$maxhit]\n";
$analysis_log .= "\%identity threshold      [$pcnt_th]\n";
$analysis_log .= "\%coverage threshold      [$pcov_th]\n";
$analysis_log .= "E-value threshold        [$eval_th]\n";
$analysis_log .= "CPU thread               [$core]\n";
$analysis_log .= "nohup                    [$nohup]\n";
$analysis_log .= "e-mail                   [$eaddr]\n";
$analysis_log .= "------------------------------------------------------------------------------------\n";

unless($nohup && $nohup eq 'pipe'){
	print $analysis_log;
}

unless($q){
	print "\nOK? (Y/n): ";
	$q = <STDIN>;
	$q =~ s/\n//;
	$q =~ s/\r//;
	unless($q){
		$q = "y";
	}
}

if($q =~ /y/i){
	my $gprefix1 = $qgenome1;
	if($gprefix1 =~ /\.fasta/){
		$gprefix1 =~ s/\.fasta//;
	}
	elsif($gprefix1 =~ /\.fa/){
		$gprefix1 =~ s/\.fa//;
	}
	
	my $prefix1 = $qfasta1;
	if($prefix1 =~ /\.fasta/){
		$prefix1 =~ s/\.fasta//;
	}
	elsif($prefix1 =~ /\.fa/){
		$prefix1 =~ s/\.fa//;
	}
	
	my $prefix2 = $qfasta2;
	if($prefix2 =~ /\.fasta/){
		$prefix2 =~ s/\.fasta//;
	}
	elsif($prefix2 =~ /\.fa/){
		$prefix2 =~ s/\.fa//;
	}
	
	unless(-e $rdir){
		system("mkdir $rdir");
	}
	
	unless($timestamp){
		my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime;
		$year += 1900;
		$mon = Num2Month($mon);
		if($hour < 10){
			$hour = "0".$hour;
		}
		if($min < 10){
			$min = "0".$min;
		}
		if($sec < 10){
			$sec = "0".$sec;
		}
		unless($timestamp){
			$timestamp = $year.$mon.$mday."-".$hour.$min.$sec;
		}
	}
	
	#---------------------------------------------------------------------------//
	chdir $rdir;
	unless(-e $qgff1){
		system("ln -s ../$qgff1 ./");
	}
	unless(-e $qgff2){
		system("ln -s ../$qgff2 ./");
	}
	unless(-e $qfasta1){
		system("ln -s ../$qfasta1 ./");
	}
	unless(-e $qfasta2){
		system("ln -s ../$qfasta2 ./");
	}
	
	unless($nohup && $nohup eq 'pipe'){
		print "\n";
	}
	
	MakeBlastDB($qfasta1, $program);
	if($qfasta1 ne $qfasta2){
		MakeBlastDB($qfasta2, $program);
	}
	
	my $lnh1 = Getseqlength($qfasta1);
	my $lnh2 = Getseqlength($qfasta2);
	
	my $result1 = $str."_db-".$prefix2."_q-".$prefix1.".txt";
	my $result2 = $str."_db-".$prefix1."_q-".$prefix2.".txt";
	
	unless($nohup && $nohup eq 'pipe'){
		print "\n";
	}
	
	print "! starting blast q=[$qfasta1] db=[$qfasta2]...\n";
	unless(-e $result1){
		my $splits = Split_query($qfasta1, $core, $prefix1, $nohup);
		my $j = @{$splits};
		
		my $thrs = [];
		my @eachRs;
		my @SHs;
		for(my $k = 0; $k < $j; $k++){
			my $sh_file = $str."_db-".$prefix2."_q-".$prefix1.".$k.sh";
			my $each_result = $str."_db-".$prefix2."_q-".$prefix1.".$k.txt";
			my $sh = "$bin -db $qfasta2 -query $splits->[$k] -num_threads 1 -out $each_result -outfmt 6 -evalue $eval_th -max_target_seqs $maxhit\n";
			
			unless(-e $each_result){
				open(my $shfh, ">", $sh_file) or die;
				print $shfh $sh;
				close $shfh;
				
				unless($nohup && $nohup eq 'pipe'){
					print "! thread [$k] starting ...\n";
				}
				
				my $thr = threads->new(\&Bash_exe, $sh_file, $k, $nohup);
				push(@{$thrs}, $thr);
				push(@eachRs, $each_result);
				push(@SHs, $sh_file);
			}
			else{
				print "! [$each_result] already exists. skipping...\n";
			}
		}
		
		foreach my $thr (@{$thrs}){
			$thr->join();
		}
		
		my $cmd = "cat ".join(" ", @eachRs)." > $result1";
		system("$cmd");
		
		foreach my $each_result (@eachRs){
			system("rm $each_result");
		}
		foreach my $sh_file (@SHs){
			system("rm $sh_file");
		}
		foreach my $each_fasta (@{$splits}){
			system("rm $each_fasta");
		}
	}
	else{
		print "! [$result1] already exists. skipping blast...\n";
	}
	
	unless($nohup && $nohup eq 'pipe'){
		print "\n";
	}
	
	if($qfasta1 eq $qfasta2){
		print "! [$qfasta1] = [$qfasta2], skip bidirectional hit analysis...\n";
		goto JUMP;
	}
	
	print "! starting blast q=[$qfasta2] db=[$qfasta1]...\n";
	unless(-e $result2){
		my $splits = Split_query($qfasta2, $core, $prefix2, $nohup);
		my $j = @{$splits};
		
		my $thrs = [];
		my @eachRs;
		my @SHs;
		for(my $k = 0; $k < $j; $k++){
			my $sh_file = $str."_db-".$prefix1."_q-".$prefix2.".$k.sh";
			my $each_result = $str."_db-".$prefix1."_q-".$prefix2.".$k.txt";
			my $sh = "$bin -db $qfasta1 -query $splits->[$k] -num_threads 1 -out $each_result -outfmt 6 -evalue $eval_th -max_target_seqs $maxhit\n";
			
			unless(-e $each_result){
				open(my $shfh, ">", $sh_file) or die;
				print $shfh $sh;
				close $shfh;
				
				unless($nohup && $nohup eq 'pipe'){
					print "! thread [$k] starting ...\n";
				}
				
				my $thr = threads->new(\&Bash_exe, $sh_file, $k, $nohup);
				push(@{$thrs}, $thr);
				push(@eachRs, $each_result);
				push(@SHs, $sh_file);
			}
			else{
				print "! [$each_result] already exists. skipping...\n";
			}
		}
		
		foreach my $thr (@{$thrs}){
			$thr->join();
		}
		
		my $cmd = "cat ".join(" ", @eachRs)." > $result2";
		system("$cmd");
		
		foreach my $each_result (@eachRs){
			system("rm $each_result");
		}
		foreach my $sh_file (@SHs){
			system("rm $sh_file");
		}
		foreach my $each_fasta (@{$splits}){
			system("rm $each_fasta");
		}
	}
	else{
		print "! [$result2] already exists. skipping blast...\n";
	}
	
	my $info1 = IntegrateLines($result1, $lnh1, $pcnt_th, $pcov_th, $eval_th, $prefix1);		#combine results when >=2 lines are present in blast result for each query
	my $info2 = IntegrateLines($result2, $lnh2, $pcnt_th, $pcov_th, $eval_th, $prefix2);
	
	my @ID1 = keys(%{$lnh1});
	my @ID2 = keys(%{$lnh2});
	@ID1 = sort {$a cmp $b} @ID1;
	@ID2 = sort {$a cmp $b} @ID2;
	my $num_ID1 = @ID1;
	my $num_ID2 = @ID2;
	
	my $hgff1 = Read_gff($qgff1);
	my $hgff2 = Read_gff($qgff2);
	my $g2t1 = $hgff1->{g2t};
	my $t2g1 = $hgff1->{t2g};
	my $g2t2 = $hgff2->{g2t};
	my $t2g2 = $hgff2->{t2g};
	
	my $cnt_bdh1 = 0;
	my $cnt_bdh2 = 0;
	my $cnt_cfh = 0;
	my $r1 = "ID 1=[".$gprefix1."],gene_id [1],Asm2sv_partner [1],ID 2=[".$prefix2."],gene_id [2],Asm2sv_BDH_parter,[1] length,[2] length,[1] alignment length,[2] alignment length,[1] \%identity,[2] \%identity,[1] \%coverage,[2] \%coverage,[1] score,[2] score,[1] Eval,[2] Eval,k,l\n";
	
	my @Header = split(/,/, $r1);
	my $numR = @Header;
	my @emptyR;
	for(my $i = 4; $i < $numR; $i++){
		push(@emptyR, "-");
	}
	
	foreach my $id1 (@ID1){
		unless($t2g1->{$id1}){
			print "! missing gene ID for [$id1]...\n";
			next;
		}
		$emptyR[2] = $lnh1->{$id1};
		
		my $gid1 = $t2g1->{$id1}; 			# e.g. Gmax_Enrei_pseudomol_v3.31_Glyma.05G137100.v4.JE1_miniprot
		my $asm2sv_gid1 = $gid1;
		my $tmp_qprefix1 = $gprefix1."_";
		$asm2sv_gid1 =~ s/$tmp_qprefix1//;
		my @tmpID1 = split(/_/, $asm2sv_gid1);
		my $num_tmpID1 = @tmpID1;
		if($num_tmpID1 <= 2){
			$asm2sv_gid1 = $tmpID1[0];
		}
		else{
			my @tmpID1r;
			for(my $i = 0; $i < $num_tmpID1 - 1; $i++){
				push(@tmpID1r, $tmpID1[$i]);
			}
			$asm2sv_gid1 = join("_", @tmpID1r);
		}
		
		unless($info1->{$id1}{num_hit}){
			$r1 .= $id1.",".$gid1.",".$asm2sv_gid1.",NA,".join(",", @emptyR)."\n";
			next;
		}
		
		my $hBDBH_partner = {};
		my $num_hit = $info1->{$id1}{num_hit};
		my $sw = 0;
		for(my $k = 0; $k < $num_hit; $k++){
			if($info1->{$id1}{$k}{hit}){
				my $hit1 = $info1->{$id1}{$k}{hit};
				
				unless($t2g2->{$hit1}){
					print "! missing gene ID for [$hit1]...\n";
					next;
				}
				my $gid_hit1 = $t2g2->{$hit1};
				
				if($info2->{$hit1}{num_hit}){
					my $num_hit1 = $info2->{$hit1}{num_hit};
					
					for(my $l = 0; $l < $num_hit1; $l++){
						if($info2->{$hit1}{$l}{hit}){
							my $hit2 = $info2->{$hit1}{$l}{hit};
							
							unless($t2g1->{$hit2}){
								print "! missing gene ID for [$hit2]...\n";
								next;
							}
							my $gid_hit2 = $t2g1->{$hit2};
							
							if($id1 eq $hit2 || $gid1 eq $gid_hit2){
								my $judge_BDH_partner = "false";
								if($asm2sv_gid1 eq $gid_hit1){
									if($k == 0 && $l == 0){
										$judge_BDH_partner = "BDBH";
										$cnt_bdh1++;
									}
									else{
										$judge_BDH_partner = "BDH";
										$cnt_cfh++;
									}
								}
								
								$r1 .= $id1.",".$gid1.",".$asm2sv_gid1.",".$hit1.",".$gid_hit1.",".$judge_BDH_partner.",".$lnh1->{$id1}.",".$lnh2->{$hit1}.",".$info1->{$id1}{$k}{alen}.",".$info2->{$hit1}{$l}{alen}.",".$info1->{$id1}{$k}{pcnt}.",".$info2->{$hit1}{$l}{pcnt}.",".$info1->{$id1}{$k}{pcov}.",".$info2->{$hit1}{$l}{pcov}.",".$info1->{$id1}{$k}{score}.",".$info2->{$hit1}{$l}{score}.",".$info1->{$id1}{$k}{Eval}.",".$info2->{$hit1}{$l}{Eval}.",".$k.",".$l."\n";
									
								$sw = 1;
								last;
							}
						}
					}
				}
			}
		}
		
		if($sw eq '0'){
			$r1 .= $id1.",".$gid1.",".$asm2sv_gid1.",NA,".join(",", @emptyR)."\n";
		}
	}
	
	print "! [$num_ID1] in [$qfasta1]\n";
	print "! [$num_ID2] in [$qfasta2]\n";
	print "! [$cnt_bdh1] bidirectional best hits\n";
#	print "! [$cnt_bdh2] bidirectional best hits (transcript variant)\n";
	print "! [$cnt_cfh] bidirectional confident hits\n";
	my $rfile1 = "BDconfidenthit_pcnt-".$pcnt_th."_pcov-".$pcov_th."_Eval-".$eval_th.".csv";
	SAVE($rfile1, $r1);
	#---------------------------------------------------------------------------//
	
	JUMP:{
		my $jump = 1;
	}
}

END:{
	unless($nohup && $nohup eq 'pipe'){
		print "\n! End of script.\n\n";
	}
}
END2:{
	my $end2 = "";
}


################################################################################
#-------------------------------------------------------------------------------
sub Read_gff{
my $gff3 = shift;

print "! reading [$gff3]...\n";
open(my $fh, "<", $gff3) or die;
my $g2t = {};
my $t2g = {};
my $hash = {};
my @GID;
my $gcnt = 0;
my @L;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /\#/){
		next;
	}
	
	my @A = split(/\t/, $line);
	unless($line){
#		print "! missing line\n";
		next;
	}
	
	my @tag = split(/\;/, $A[8]);
	
	if($A[2] eq 'gene'){
		my $gid;
		foreach my $val (@tag){
			if($val =~ /ID\=/){
				$gid = $val;
				$gid =~ s/ID\=//;
				last;
			}
		}
		$gcnt++;
		
		push(@GID, $gid);
		$hash->{$gid}{gff} = $line."\n";
	}
	elsif($A[2] eq 'transcript' || $A[2] eq 'mRNA' || $A[2] eq 'tRNA' || $A[2] eq 'snoRNA' || $A[2] eq 'rRNA' || $A[2] eq 'SRP_RNA' || $A[2] eq 'snRNA' || $A[2] eq 'microRNA'){
		my $tid;
		my $gid;
		foreach my $val (@tag){
			if($val =~ /ID\=/){
				$tid = $val;
				$tid =~ s/ID\=//;
			}
			elsif($val =~ /Parent\=/){
				$gid = $val;
				$gid =~ s/Parent\=//;
			}
		}
		
		$g2t->{$gid} .= $tid."\n";
		$t2g->{$tid} = $gid;
		$hash->{$tid}{gff} .= $line."\n";
		$hash->{$gid}{num_variant} += 1;
	}
	
	push(@L, $line);
}
close $fh;

foreach my $line (@L){
	my @A = split(/\t/, $line);
	my @tag = split(/\;/, $A[8]);
	
#	elsif($A[2] eq 'exon' || $A[2] eq 'CDS' || $A[2] eq 'start_codon' || $A[2] eq 'stop_codon'){
	if($A[2] eq 'CDS'){
		my $tid;
		foreach my $val (@tag){
			if($val =~ /Parent\=/){
				$tid = $val;
				$tid =~ s/Parent\=//;
			}
		}
		
		my @tmp = split(/,/, $tid);
		if($tid =~ /AT/){
			$tid = $tmp[0];
		}
		
		unless($t2g->{$tid}){
			next;
		}
		
		my $gid = $t2g->{$tid};
		$hash->{$gid}{num_CDS} += 1;
		$hash->{$tid}{gff} .= $line."\n";
	}
}

print " [$gcnt] genes\n";

my $rh = {};
$rh->{hash} = $hash;
$rh->{g2t} = $g2t;
$rh->{t2g} = $t2g;

return $rh;
}


#-------------------------------------------------------------------------------
sub IntegrateLines{
my $file = shift;
my $lnh = shift;
my $pcnt_th = shift;
my $pcov_th = shift;
my $eval_th = shift;
my $prefix = shift;

open(my $fh, "<", $file);
my $hash = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	my @A = split(/\t/, $line);
	$hash->{$A[0]} .= $line."\n";
}
close $fh;

my @ID = keys(%{$hash});
@ID = sort {$a cmp $b} @ID;
my $numID = @ID;

my $info = {};
my $cnt_hit = 0;
foreach my $id (@ID){
#QUERY	HIT	PERCENT_IDENTITY	ALIGNMENT_LENGTH	MISMATCH_COUNT	GAP_OPEN_COUNT	QUERY_START	QUERY_END	HIT_START	HIT_END	E-VALUE	SCORE
	unless($lnh->{$id}){
		print "! error : missing seq length info for [$id] [$file]\n";
		die;
	}
	
	my @L = split(/\n/, $hash->{$id});
	my $ihash = {};
	foreach my $l (@L){
		my @A = split(/\t/, $l);
		my $hit = $A[1];
		$ihash->{$hit} .= $l."\n";
	}
	
	my @H = keys(%{$ihash});
	@H = sort {$a cmp $b} @H;
	
	my $AoA = [];
	foreach my $hit (@H){
		my $sum_alen = 0;
		my $sum_mis = 0;
		my $sum_gap = 0;
		my $sum_score = 0;
		my $sum_pcnt = 0;
		my $max_pcnt = 0;
		my $min_eval = 10000;
		my $potential_candidate = 0;
		
		my @eL = split(/\n/, $ihash->{$hit});
		
		my $cnt = 0;
		foreach my $l (@eL){
			my @A = split(/\t/, $l);
			$sum_pcnt += $A[2];
			$sum_alen += $A[3];
			$sum_mis += $A[4];
			$sum_gap += $A[5];
			$sum_score += $A[11];
			
			if($A[2] > $max_pcnt){
				$max_pcnt = $A[2];
			}
			
			if($A[10] < $min_eval){
				$min_eval = $A[10];
			}
			
			$cnt++;
		}
		
		if($id =~ /$hit/ || $hit =~ /$id/){
			$potential_candidate = 1;
		}
		
		my $pcnt_avr = sprintf("%.2f", $sum_pcnt / $cnt);
#		my $pcnt = sprintf("%.2f", $sum_pcnt / $cnt);
		my $pcnt = sprintf("%.2f", $max_pcnt);
		my $pcov = sprintf("%.2f", $sum_alen / $lnh->{$id} * 100);
		my $Eval = $min_eval;
		my $sum_misgap = $sum_gap + $sum_mis;
		
		my @tmp;
		push(@tmp, $hit);
		push(@tmp, $pcnt);
		push(@tmp, $pcov);
		push(@tmp, $Eval);
		push(@tmp, $sum_score);
		push(@tmp, $sum_alen);
		push(@tmp, $sum_misgap);
		push(@tmp, $pcnt_avr);
		push(@tmp, $potential_candidate);
		push(@{$AoA}, \@tmp);
	}
	
	@{$AoA} = sort {$a->[3] <=> $b->[3]} @{$AoA};
	@{$AoA} = sort {$b->[1] <=> $a->[1]} @{$AoA};
	@{$AoA} = sort {$b->[2] <=> $a->[2]} @{$AoA};
	@{$AoA} = sort {$b->[8] <=> $a->[8]} @{$AoA};
	@{$AoA} = sort {$b->[4] <=> $a->[4]} @{$AoA};
	
	$info->{$id}{num_hit} = @{$AoA};
	
	my $k = 0;
	foreach my $A (@{$AoA}){
		if($A->[1] >= $pcnt_th && $A->[2] >= $pcov_th && $A->[3] <= $eval_th){
			$info->{$id}{$k}{hit} = $A->[0];
			$info->{$id}{$k}{pcnt} = $A->[1];
			$info->{$id}{$k}{pcov} = $A->[2];
			$info->{$id}{$k}{Eval} = $A->[3];
			$info->{$id}{$k}{score} = $A->[4];
			$info->{$id}{$k}{alen} = $A->[5];
			$info->{$id}{$k}{misgap} = $A->[6];
			
			$info->{$id}{hitinfo}{$A->[0]}{pcnt} = $A->[1];
			$info->{$id}{hitinfo}{$A->[0]}{pcov} = $A->[2];
			$info->{$id}{hitinfo}{$A->[0]}{Eval} = $A->[3];
			$info->{$id}{hitinfo}{$A->[0]}{score} = $A->[4];
			$info->{$id}{hitinfo}{$A->[0]}{alen} = $A->[5];
			$info->{$id}{hitinfo}{$A->[0]}{misgap} = $A->[6];
			$info->{$id}{hitinfo}{$A->[0]}{pcnt_avr} = $A->[7];
			
			$k++;
			$cnt_hit++;
		}
	}
}

print "! [$prefix] : numID = [$numID] redundant_hit_count = [$cnt_hit] (\%identity >= $pcnt_th, \%coverage >= $pcov_th, E-value <= $eval_th)\n";

return $info;
}


#-------------------------------------------------------------------------------
sub Getseqlength{
my $query = shift;

print "! reading [$query]...\n";

my $fh;
if($query =~ /\.gz/){
	open($fh, "gzip -dc $query |") or die;
}
else{
	open($fh, "<", $query) or die;
}
my $cnt = 0;
my $ID;
my $seq;
my $hash = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ />/){
		if($seq){
			$hash->{$ID} = length($seq);
		}
		$ID = $line;
		$ID =~ s/\>//;
		my @tmp = split(/\s/, $ID);
		$ID = $tmp[0];
		$seq = "";
		$cnt++;
	}
	else{
		$seq .= $line;
	}
}
close $fh;

if($seq){
	$hash->{$ID} = length($seq);
}

return $hash;
}


#-------------------------------------------------------------------------------
sub Bash_exe{
my $sh_file = shift;
my $num = shift;
my $nohup = shift;

system("bash $sh_file\n");

unless($nohup && $nohup eq 'pipe'){
	print "! thread [$num] completed.\n";
}

return $num;
}


#-------------------------------------------------------------------------------
sub Split_query{
my $query = shift;
my $core = shift;
my $qstr = shift;
my $nohup = shift;

print "reading [$query]...\n";
my $fh;
if($query =~ /\.gz/){
	open($fh, "gzip -dc $query |") or die;
}
else{
	open($fh, "<", $query) or die;
}
my $cnt = 0;
my $ID;
my $seq;
my @IDs;
my $hash = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ />/){
		unless($cnt == 0){
			push(@IDs, $ID);
			$hash->{$ID} = $seq;
		}
		$ID = $line;
		$ID =~ s/>//;
		$seq = "";
		$cnt++;
	}
	else{
		$seq .= $line;
	}
}
close $fh;

if($seq){
	push(@IDs, $ID);
	$hash->{$ID} = $seq;
}

print "! [$cnt] sequences\n";

my $d = int($cnt / $core) + 1;

unless($nohup && $nohup eq 'pipe'){
	print "splitting query sequence into [$core] files...\n";
	print "each file will contain [$d] sequences...\n";
}
else{
	print "blast with [$core] partitions...\n";
}

my $seg = $d;
my $tmpfa;
my @rfiles;
my $n = 0;
for(my $i = 0; $i < $cnt; $i++){
	if($i == $seg){
		my $fasta = "$qstr.$n.fa";
		open(my $rfh, ">", $fasta) or die;
		print $rfh $tmpfa;
		close $rfh;
#		print "! [$fasta]\n";
		
		push(@rfiles, $fasta);
		
		$n++;
		$tmpfa = "";
		$seg += $d;
	}
	$tmpfa .= ">$IDs[$i]\n$hash->{$IDs[$i]}\n";
}

if($tmpfa){
	my $fasta = "$qstr.$n.fa";
	open(my $rfh, ">", $fasta) or die;
	print $rfh $tmpfa;
	close $rfh;
#	print "! [$fasta]\n";
	
	push(@rfiles, $fasta);
}
#print "! done\n\n";

return \@rfiles;
}


#-------------------------------------------------------------------------------
sub MakeBlastDB{
my $fasta = shift;
my $program = shift;

my $blastdb = $fasta;

if($program eq '2'){
	my $db1 = $blastdb."\.phd";
	my $db2 = $blastdb."\.phi";
	my $db3 = $blastdb."\.phr";
	my $db4 = $blastdb."\.pin";
	my $db5 = $blastdb."\.pog";
	my $db6 = $blastdb."\.psd";
	my $db7 = $blastdb."\.psi";
	my $db8 = $blastdb."\.psq";
	unless(-e $db1 && -e $db2 && -e $db3 && -e $db4 && -e $db5 && -e $db6 && -e $db7 && -e $db8){
		print "! making blast db for [$blastdb]...\n";
		`makeblastdb -in $blastdb -dbtype prot -hash_index`;
	}
	else{
		print "! blast db already exists for [$blastdb]\n";
	}
}
elsif($program eq '1'){
	my $db1 = $blastdb."\.nhd";
	my $db2 = $blastdb."\.nhi";
	my $db3 = $blastdb."\.nhr";
	my $db4 = $blastdb."\.nin";
	my $db5 = $blastdb."\.nog";
	my $db6 = $blastdb."\.nsd";
	my $db7 = $blastdb."\.nsi";
	my $db8 = $blastdb."\.nsq";
	unless(-e $db1 && -e $db2 && -e $db3 && -e $db4 && -e $db5 && -e $db6 && -e $db7 && -e $db8){
		print "! making blast db for [$blastdb]...\n";
		`makeblastdb -in $blastdb -dbtype nucl -hash_index`;
	}
	else{
		print "! blast db already exists for [$blastdb]\n";
	}
}

}


#-------------------------------------------------------------------------------
sub Num2Month{
my $mon = shift;

if($mon == 0){$mon = "Jan";}
elsif($mon == 1){$mon = "Feb";}
elsif($mon == 2){$mon = "Mar";}
elsif($mon == 3){$mon = "Apr";}
elsif($mon == 4){$mon = "May";}
elsif($mon == 5){$mon = "Jun";}
elsif($mon == 6){$mon = "Jul";}
elsif($mon == 7){$mon = "Aug";}
elsif($mon == 8){$mon = "Sep";}
elsif($mon == 9){$mon = "Oct";}
elsif($mon == 10){$mon = "Nov";}
elsif($mon == 11){$mon = "Dec";}

return $mon;
}


#----------------------------------------------------------
sub SAVE{
my $file = shift;
my $str = shift;

open(my $fh, ">", $file) or die;
print $fh $str;
close $fh;

print "! output [$file]\n";

}


#----------------------------------------------------------
sub APPEND{
my $file = shift;
my $header = shift;
my $str = shift;

if(-e $file){
	open(my $fh, ">>", $file) or die;
	print $fh $str;
	close $fh;
}
else{
	$str = $header.$str;
	
	open(my $fh, ">", $file) or die;
	print $fh $str;
	close $fh;
}

}


#----------------------------------------------------------
sub OPEN{
my $file = shift;

open(my $fh, "<", $file) or die;
my $str;
while(my $line = <$fh>){
	$str .= $line;
}
close $fh;

return $str;
}


#----------------------------------------------------------
sub File_select{
my $keyword = shift;
my $str = shift;

my @files2;
my $q1;

unless($keyword){
	print "\nEnter keyword for file search: ";
	$keyword = <STDIN>;
	$keyword =~ s/\n//;
	$keyword =~ s/\r//;
}

my @files = glob("*");

foreach my $file (@files){
	if($file =~ /$keyword/){
		unless($file =~ /\.fai/ || $file =~ /\.fa\.n../ || $file =~ /\.fa\.p../ || $file =~ /\.fasta\.n../ || $file =~ /\.fasta\.p../ || $file =~ /\.ffa\.n../ || $file =~ /\.csv/){
			push(@files2, $file);
		}
	}
}

print "\n---------------------------------------------------------------------\n";
my $cnt0 = 0;
foreach my $file (@files2){
	print "[$cnt0] $file\n";
	$cnt0++;
}
print "---------------------------------------------------------------------\n";

unless($q1){
	print "Select $str: ";
	$q1 = <STDIN>;
	$q1 =~ s/\n//;
	$q1 =~ s/\r//;
}
unless($q1){
	$q1 = 0;
}

if($files2[$q1] && -e $files2[$q1]){
	print "[$files2[$q1]]\n";
}

return $files2[$q1];
}

