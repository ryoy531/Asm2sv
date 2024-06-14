#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "Find_orthologue_partner.pl version 1.01\n";
$version .= "last update: [2020\/1\/4]\n";
$version .= "copyright: ryoichi yano [ryoichiy104\@gmail.com]\n";

#-------------------------------------------------------------------------------

my $rdir = shift;
my $bin_gffread = shift;
my $bin_Blat2hint_mRNA = shift;
my $bin_blat = shift;
my $bin_blat2hint = shift;
my $bin_Blast_BDH = shift;
my $genome1 = shift;
my $genome2 = shift;
my $gff1 = shift;
my $gff2 = shift;
my $maxhit = shift;
my $th_condition = shift;
my $core = shift;
my $q = "y";

my $err = 0;
unless($rdir){
	print "! missing directory name...\n";
	$err++;
}
$err += FileCheck1($bin_gffread, "gffread");
$err += FileCheck1($bin_Blat2hint_mRNA, "Blat2hint_mRNA.pl");
$err += FileCheck1($bin_blat, "blat");
$err += FileCheck1($bin_blat2hint, "blat2hints.pl");
$err += FileCheck1($bin_Blast_BDH, "Blast_BDH.pl");
$err += FileCheck1($genome1, "genome1");
$err += FileCheck1($genome2, "genome2");
$err += FileCheck1($gff1, "gff1");
$err += FileCheck1($gff2, "gff2");

if($err > 0){
	print "! abort script due to error...\n";
	die;
}

$maxhit = 100;
my $pcnt_th = 90;
my $pcov_th = 1;
my $eval_th = 1e-30;

my $pep_pcnt_th = 90;
my $pep_pcov_th = 1;
my $pep_eval_th = 1e-50;

#my $pcnt_th = 50;
#my $pcov_th = 1;
#my $eval_th = 1e-3;

#my $pep_pcnt_th = 20;
#my $pep_pcov_th = 1;
#my $pep_eval_th = 1e-3;

unless($th_condition){
	print "\n-------------------------------------------------------------------------------\n";
	print "[1] n (90\%, 90\%, 1e-50) | p (90\%, 90\%, 1e-50)\n";
	print "[2] n (90\%, 1\%, 1e-50) | p (90\%, 1\%, 1e-50)\n";
	print "[3] n (70\%, 90\%, 1e-30) | p (70\%, 90\%, 1e-30)\n";
	print "[4] n (70\%, 1\%, 1e-30) | p (70\%, 1\%, 1e-30)\n";
	print "\n";
	print "*blastn/p thresholds (\%similarity, \%alignment coverage, E-value)\n";
	print "-------------------------------------------------------------------------------\n";
	print "Select blast threshold condition (default 4): ";
	$th_condition = <STDIN>;
	$th_condition =~ s/\n//;
	$th_condition =~ s/\r//;
	unless($th_condition){
		$th_condition = 4;
	}
	unless($th_condition eq '1' || $th_condition eq '2' || $th_condition eq '3' || $th_condition eq '4'){
		$th_condition = 4;
	}
	print "! [$th_condition]\n";
}

if($th_condition eq '1'){
	$pcnt_th = 90;
	$pcov_th = 90;
	$eval_th = 1e-50;
	
	$pep_pcnt_th = 90;
	$pep_pcov_th = 90;
	$pep_eval_th = 1e-50;
}
elsif($th_condition eq '2'){
	$pcnt_th = 90;
	$pcov_th = 1;
	$eval_th = 1e-50;
	
	$pep_pcnt_th = 90;
	$pep_pcov_th = 1;
	$pep_eval_th = 1e-50;
}
elsif($th_condition eq '3'){
	$pcnt_th = 70;
	$pcov_th = 90;
	$eval_th = 1e-30;
	
	$pep_pcnt_th = 70;
	$pep_pcov_th = 90;
	$pep_eval_th = 1e-30;
}
elsif($th_condition eq '4'){
	$pcnt_th = 70;
	$pcov_th = 1;
	$eval_th = 1e-30;
	
	$pep_pcnt_th = 70;
	$pep_pcov_th = 1;
	$pep_eval_th = 1e-30;
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
			die;
		}
	}
	print "! [$core]\n";
}

my $analysis_log = "";
$analysis_log .= "\n------------------------------------------------------------------------------------\n";
$analysis_log .= "genome fasta (query)     [$genome1]\n";
$analysis_log .= "genome fasta (hint)      [$genome2]\n";
$analysis_log .= "gff1 (query)             [$gff1]\n";
$analysis_log .= "gff2 (hint)              [$gff2]\n";
$analysis_log .= "max hit per query        [$maxhit]\n";
$analysis_log .= "threshold condition      [$th_condition]\n";
$analysis_log .= "CPU thread               [$core]\n";
$analysis_log .= "------------------------------------------------------------------------------------\n";

print "$analysis_log\n";

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
	my $prefix1 = Gff2prefix($gff1);
	my $prefix2 = Gff2prefix($gff2);
	
	my $prefix3 = $genome1;
	if($prefix3 =~ /\.fasta/){
		$prefix3 =~ s/\.fasta//;
	}
	elsif($prefix3 =~ /\.fa/){
		$prefix3 =~ s/\.fa//;
	}
	
	unless(-e $rdir){
		system("mkdir $rdir");
	}
	
	#---------------------------------------------------------------------------//
	chdir $rdir;
	unless(-e $genome1){
		system("ln -s ../$genome1 ./");
	}
	unless(-e $genome2){
		system("ln -s ../$genome2 ./");
	}
	unless(-e $gff1){
		system("ln -s ../$gff1 ./");
	}
	unless(-e $gff2){
		system("ln -s ../$gff2 ./");
	}
	
	#------------------------------------------------------//
	print "! searching for candidate genomic positions...\n";
	my $ori_qfasta2 = $prefix2."_exon_ori.fasta";
	my $pfasta2 = $prefix2."_protein.fasta";
	
	unless(-e $ori_qfasta2 && -e $pfasta2){
		print "! converting gff to sequence [$gff2]...\n";
		my $cmd = "$bin_gffread $gff2 -g $genome2 -w $ori_qfasta2 -y $pfasta2 > /dev/null 2>&1";
		CMD($cmd);
	}
	else{
		print "! [$ori_qfasta2] [$pfasta2] already exist. skip gffread...\n";
	}
	
	Fasta_nameclean($pfasta2);
	Fasta_nameclean($ori_qfasta2);
	
	my $b2h_gff = "Blat2hint_mRNA_db-".$prefix3."_q-".$prefix2."_exon_ori/hint.E.gff";
	my $hint_gff2 = $prefix2."_hint.E.gff";
	unless(-e $b2h_gff){
		print "! making hint by Blat2hint_mRNA.pl ...\n";
		my $cmd_blat2hint = "perl $bin_Blat2hint_mRNA $bin_blat $bin_blat2hint $genome1 $ori_qfasta2 90 $core y pipe n";
		CMD($cmd_blat2hint);
	}
	else{
		print "! [$b2h_gff] already exists. skip Blat2hint_mRNA.pl ...\n";
	}
	
	if(-e $b2h_gff){
		if(-e $hint_gff2){
			system("rm $hint_gff2");
		}
		system("cp $b2h_gff $hint_gff2");
	}
	else{
		print "! missing [$b2h_gff], abort script...\n";
		die;
	}
	
	#------------------------------------------------------//
	my $bin = "";
	my $str = "";
	my $hgff2ori = Gff_read_as_hash($gff2);
	my $htid2gid2ori = $hgff2ori->{tid2gid};
	my $hgid2tid2ori = $hgff2ori->{gid2tid};
	my $hgene2ori = $hgff2ori->{gene};
	
	my $update_hint_blast = 0;
	
	if($update_hint_blast eq '1'){
		my $ori_gfasta2 = $prefix2."_genomic_ori.fasta";
		my $genomic_flanking_span = 20 * 1000;		# 20 kb
		unless(-e $ori_gfasta2){
			Gff3togenomicDNA($genome2, $gff2, $ori_gfasta2, $genomic_flanking_span);
		}
		else{
			print "! [$ori_gfasta2] already exists...\n";
		}
		
		$bin = "blastn";
		$str = "bgN";
		my $gresult2 = $str."_db-".$prefix1."_q-".$prefix2.".txt";
		
		unless(-e $gresult2){
			MakeBlastDB($genome1, 1);
			print "! starting genomie blast q=[$ori_gfasta2] db=[$genome1]...\n";
			my $splits = Split_query($ori_gfasta2, $core, $prefix2, "n");
			my $j = @{$splits};
			
			my $thrs = [];
			my @eachRs;
			my @SHs;
			for(my $k = 0; $k < $j; $k++){
				my $sh_file = $str."_db-".$prefix1."_q-".$prefix2.".$k.sh";
				my $each_result = $str."_db-".$prefix1."_q-".$prefix2.".$k.txt";
				my $sh = "$bin -db $genome1 -query $splits->[$k] -num_threads 1 -out $each_result -outfmt 6 -evalue 1e-80 -max_target_seqs 5\n";
				
				unless(-e $each_result){
					open(my $shfh, ">", $sh_file) or die;
					print $shfh $sh;
					close $shfh;
					
					my $thr = threads->new(\&Bash_exe, $sh_file, $k, "n");
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
			
			my $cmd = "cat ".join(" ", @eachRs)." > $gresult2";
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
			print "! [$gresult2] already exists. skipping blast...\n";
		}
		
		my $hint_gff2_bkup = $prefix2."_hint.E.gff.original";
		if(! -e $hint_gff2_bkup){
			system("cp $hint_gff2 $hint_gff2_bkup");
		}
		
		Keep_confident_hintGff($hint_gff2, $gresult2, $htid2gid2ori, 97);		# discard BLAT alignment record if it positions outside of BLAST hit region
	}
	
	#------------------------------------------------------//
	my $hgff1 = Gff_read_as_hash($gff1);
	my $hgff2 = Gff_read_as_hash($hint_gff2);
	
	my $hCDS1 = $hgff1->{hCDS};
	my $hexon1 = $hgff1->{hexon};
	my $htidall1 = $hgff1->{htidall};
	my $hgene1 = $hgff1->{gene};
	my $htid2gid1 = $hgff1->{tid2gid};
	my $hgid2tid1 = $hgff1->{gid2tid};
	
	my $hCDS2 = $hgff2->{hCDS};
	my $hexon2 = $hgff2->{hexon};
	my $htidall2 = $hgff2->{htidall};
	my $hgene2 = $hgff2->{gene};
	my $htid2gid2 = $hgff2->{tid2gid};
	my $hgid2tid2 = $hgff2->{gid2tid};
	
	my @TID1 = keys(%{$htidall1});
	@TID1 = sort {$a cmp $b} @TID1;
	
	my @TID2 = keys(%{$htidall2});
	@TID2 = sort {$a cmp $b} @TID2;
	
#	my $g2t1 = Tid2hash(\@TID1);		#hash of gid -> tid
#	my $g2t2 = Tid2hash(\@TID2);		#hash of gid -> tid
	
	#------------------------------------------------------//
	my $qfasta1 = $prefix1."_exon.fasta";
	my $qfasta2 = $prefix2."_exon.fasta";
	my $pfasta1 = $prefix1."_protein.fasta";
	
	unless(-e $qfasta1 && -e $pfasta1){
		print "! converting gff to sequence [$gff1]...\n";
		my $cmd = "$bin_gffread $gff1 -g $genome1 -w $qfasta1 -y $pfasta1 > /dev/null 2>&1";
		CMD($cmd);
	}
	else{
		print "! [$qfasta1] [$pfasta1] already exist. skip gffread...\n";
	}
	
	unless(-e $qfasta2){
		print "\n! converting gff to sequence [$hint_gff2]...\n";
		my $cmd = "$bin_gffread $hint_gff2 -g $genome1 -w $qfasta2 > /dev/null 2>&1";
		CMD($cmd);
	}
	else{
		print "! [$qfasta2] already exist. skip gffread...\n";
	}
	
	Fasta_nameclean($qfasta1);
	Fasta_nameclean($qfasta2);
	Fasta_nameclean($pfasta1);
	my $hseq_qexon = Open_fasta_as_hash($qfasta1);
	my $hseq_qprot = Open_fasta_as_hash($pfasta1);
	
	my $CDScheck_info = Check_protein_coding($hseq_qexon, $hseq_qprot);		#analyze whether transcript is protein-coding or non-coding, return hash
	
	#------------------------------------------------------//
	my $cpun = $core;
	my $cpup = $core;
	if($cpun > 16){
		$cpun = 16;
	}
	
	my $q1_prefix = Fasta2prefix($qfasta1);
	my $q2_prefix = Fasta2prefix($qfasta2);
	my $p1_prefix = Fasta2prefix($pfasta1);
	my $p2_prefix = Fasta2prefix($pfasta2);
	my $rdir_bdhn = "_blastBDHn";
	my $rdir_bdhp = "_blastBDHp";
	my $log_bdhn = "log_Blast_BDHn.txt";
	my $log_bdhp = "log_Blast_BDHp.txt";
	my $result1 = "./$rdir_bdhn/bN_db-".$q2_prefix."_q-".$q1_prefix.".txt";
	my $result2 = "./$rdir_bdhn/bN_db-".$q1_prefix."_q-".$q2_prefix.".txt";
	my $presult1 = "./$rdir_bdhp/bP_db-".$p2_prefix."_q-".$p1_prefix.".txt";
	my $presult2 = "./$rdir_bdhp/bP_db-".$p1_prefix."_q-".$p2_prefix.".txt";
	my $rfile_bdhn = "./$rdir_bdhn/BDconfidenthit_pcnt-".$pcnt_th."_pcov-".$pcov_th."_Eval-".$eval_th.".csv";
	my $rfile_bdhp = "./$rdir_bdhp/BDconfidenthit_pcnt-".$pep_pcnt_th."_pcov-".$pep_pcov_th."_Eval-".$pep_eval_th.".csv";
	my $cmdn = "$bin_Blast_BDH $rdir_bdhn $genome1 $genome1 $gff1 $hint_gff2 $qfasta1 $qfasta2 1 20 $pcnt_th $pcov_th $eval_th $cpun y n n > $log_bdhn 2>&1";
	my $cmdp = "$bin_Blast_BDH $rdir_bdhp $genome1 $genome2 $gff1 $gff2 $pfasta1 $pfasta2 2 20 $pep_pcnt_th $pep_pcov_th $pep_eval_th $cpup y n n > $log_bdhp 2>&1";

	if(! -e $result1 && ! -e $result2){
		print "! reciprocal BLASTn with [$cpun] threads...\n";
		if(system($cmdn) != 0){
			print "! failed.\n";
			die;
		}
	}
	else{
		print "! already done, skip BLASTn...\n";
	}
	Rmfiles("$rdir_bdhn");
	
	if(! -e $presult1 && ! -e $presult2){
		print "! reciprocal BLASTp with [$cpup] threads...\n";
		if(system($cmdp) != 0){
			print "! failed.\n";
			die;
		}
	}
	else{
		print "! already done, skip BLASTp...\n";
	}
	Rmfiles("$rdir_bdhp");	
	
	$err += FileCheck1($result1, "BLASTn result1");
	$err += FileCheck1($result2, "BLASTn result2");
	$err += FileCheck1($presult1, "BLASTp result1");
	$err += FileCheck1($presult2, "BLASTp result2");
	
	if($err > 0){
		print "! missing some required files...\n";
		die;
	}
	
	#------------------------------------------------------//
	print "! integrating results [blastn]...\n";
	my $lnh1 = Getseqlength($qfasta1);
	my $lnh2 = Getseqlength($qfasta2);
	my $info1 = IntegrateLines($result1, $lnh1, $pcnt_th, $pcov_th, $eval_th, $prefix1);		#combine results when >=2 lines are present in blast result for each query
	my $info2 = IntegrateLines($result2, $lnh2, $pcnt_th, $pcov_th, $eval_th, $prefix2);
	
	my @ID1 = keys(%{$lnh1});
	my @ID2 = keys(%{$lnh2});
	@ID1 = sort {$a cmp $b} @ID1;
	@ID2 = sort {$a cmp $b} @ID2;
	my $num_ID1 = @ID1;
	my $num_ID2 = @ID2;
	
	my $BDhitsAll = {};
	my $BDhits = {};
	my $nrhits = {};
	my $cnt_cfh = 0;
	my $r1 = "1=[".$prefix1."],2=[".$prefix2."],[1] length,[2] length,[1] alignment length,[2] alignment length,[1] \%identity,[2] \%identity,[1] \%coverage,[2] \%coverage,[1] score,[2] score,[1] Eval,[2] Eval,k,l\n";
	foreach my $id1 (@ID1){
		if(! $info1->{$id1}{num_hit}){
			$r1 .= $id1.",NA,-,-,-,-,-,-,-,-,-,-\n";
			next;
		}
		
		my $num_hit = $info1->{$id1}{num_hit};
		my $sw = 0;
		
		for(my $k = 0; $k < $num_hit; $k++){
			if($info1->{$id1}{$k}{hit}){
				my $hit1 = $info1->{$id1}{$k}{hit};
				
				if($info2->{$hit1}{num_hit}){
					my $num_hit1 = $info2->{$hit1}{num_hit};
					
					for(my $l = 0; $l < $num_hit1; $l++){
						if($info2->{$hit1}{$l}{hit}){
							my $hit2 = $info2->{$hit1}{$l}{hit};
							
							if($id1 eq $hit2){
								$r1 .= $id1.",".$hit1.",".$lnh1->{$id1}.",".$lnh2->{$hit1}.",".$info1->{$id1}{$k}{alen}.",".$info2->{$hit1}{$l}{alen}.",".$info1->{$id1}{$k}{pcnt}.",".$info2->{$hit1}{$l}{pcnt}.",".$info1->{$id1}{$k}{pcov}.",".$info2->{$hit1}{$l}{pcov}.",".$info1->{$id1}{$k}{score}.",".$info2->{$hit1}{$l}{score}.",".$info1->{$id1}{$k}{Eval}.",".$info2->{$hit1}{$l}{Eval}.",".$k.",".$l."\n";
								my $revname1 = $hit1;
								$revname1 =~ s/grp\=//;
								
								$BDhitsAll->{$id1}{$revname1} += $info1->{$id1}{$k}{score} + $info2->{$hit1}{$l}{score};
								
								if($BDhits->{$id1}{$revname1}){
									if($BDhits->{$id1}{$revname1} < ($info1->{$id1}{$k}{score} + $info2->{$hit1}{$l}{score})){
										$BDhits->{$id1}{$revname1} = $info1->{$id1}{$k}{score} + $info2->{$hit1}{$l}{score};
									}
								}
								else{
									$BDhits->{$id1}{$revname1} = $info1->{$id1}{$k}{score} + $info2->{$hit1}{$l}{score};
								}
								
								$nrhits->{$id1} .= $revname1.",";
								
								$sw = 1;
								$cnt_cfh++;
							}
						}
					}
				}
			}
		}
		
		if($sw eq '0'){
			$r1 .= $id1.",NA,-,-,-,-,-,-,-,-,-,-\n";
		}
	}
	
	print "! [$num_ID1] in [$qfasta1]\n";
	print "! [$num_ID2] in [$qfasta2]\n";
	print "! [$cnt_cfh] bidirectional confident hits\n";
	my $rfile1 = "BDconfidenthit_transcript_".$prefix1."_pcnt-".$pcnt_th."_pcov-".$pcov_th."_Eval-".$eval_th.".csv";
	SAVE($rfile1, $r1);
	
	#------------------------------------------------------//
	print "! integrating presults [blastp]...\n";
	my $plnh1 = Getseqlength($pfasta1);
	my $plnh2 = Getseqlength($pfasta2);
	my $pinfo1 = IntegrateLines($presult1, $plnh1, $pep_pcnt_th, $pep_pcov_th, $pep_eval_th, $prefix1);		#combine presults when >=2 lines are present in blast presult for each query
	my $pinfo2 = IntegrateLines($presult2, $plnh2, $pep_pcnt_th, $pep_pcov_th, $pep_eval_th, $prefix2);
	
	my @pID1 = keys(%{$plnh1});
	my @pID2 = keys(%{$plnh2});
	@pID1 = sort {$a cmp $b} @pID1;
	@pID2 = sort {$a cmp $b} @pID2;
	my $pnum_ID1 = @pID1;
	my $pnum_ID2 = @pID2;
	
	my $pBDhits = {};
	my $nrphits = {};
	my $cnt_pcfh = 0;
	my $pr1 = "1=[".$prefix1."],2=[".$prefix2."],[1] length,[2] length,[1] alignment length,[2] alignment length,[1] \%identity,[2] \%identity,[1] \%coverage,[2] \%coverage,[1] score,[2] score,[1] Eval,[2] Eval,k,l\n";
	foreach my $id1 (@pID1){
		if(! $pinfo1->{$id1}{num_hit}){
			$pr1 .= $id1.",NA,-,-,-,-,-,-,-,-,-,-\n";
			next;
		}
		
		my $pnum_hit = $pinfo1->{$id1}{num_hit};
		my $sw = 0;
		
		for(my $k = 0; $k < $pnum_hit; $k++){
			if($pinfo1->{$id1}{$k}{hit}){
				my $hit1 = $pinfo1->{$id1}{$k}{hit};
				
				if($pinfo2->{$hit1}{num_hit}){
					my $pnum_hit1 = $pinfo2->{$hit1}{num_hit};
					
					for(my $l = 0; $l < $pnum_hit1; $l++){
						if($pinfo2->{$hit1}{$l}{hit}){
							my $hit2 = $pinfo2->{$hit1}{$l}{hit};
							
							if($id1 eq $hit2){
								$pr1 .= $id1.",".$hit1.",".$plnh1->{$id1}.",".$plnh2->{$hit1}.",".$pinfo1->{$id1}{$k}{alen}.",".$pinfo2->{$hit1}{$l}{alen}.",".$pinfo1->{$id1}{$k}{pcnt}.",".$pinfo2->{$hit1}{$l}{pcnt}.",".$pinfo1->{$id1}{$k}{pcov}.",".$pinfo2->{$hit1}{$l}{pcov}.",".$pinfo1->{$id1}{$k}{score}.",".$pinfo2->{$hit1}{$l}{score}.",".$pinfo1->{$id1}{$k}{Eval}.",".$pinfo2->{$hit1}{$l}{Eval}.",".$k.",".$l."\n";
								$BDhitsAll->{$id1}{$hit1} += $pinfo1->{$id1}{$k}{score} + $pinfo2->{$hit1}{$l}{score};
								
								if($pBDhits->{$id1}{$hit1}){
									if($pBDhits->{$id1}{$hit1} < ($pinfo1->{$id1}{$k}{score} + $pinfo2->{$hit1}{$l}{score})){
										$pBDhits->{$id1}{$hit1} = $pinfo1->{$id1}{$k}{score} + $pinfo2->{$hit1}{$l}{score};
									}
								}
								else{
									$pBDhits->{$id1}{$hit1} = $pinfo1->{$id1}{$k}{score} + $pinfo2->{$hit1}{$l}{score};
								}
								
								$nrphits->{$id1} .= $hit1.",";
								
								$sw = 1;
								$cnt_pcfh++;
							}
						}
					}
				}
			}
		}
		
		if($sw eq '0'){
			$pr1 .= $id1.",NA,-,-,-,-,-,-,-,-,-,-\n";
		}
	}
	
	print "! [$pnum_ID1] in [$pfasta1]\n";
	print "! [$pnum_ID2] in [$pfasta2]\n";
	print "! [$cnt_pcfh] bidirectional confident hits\n";
	my $rfile2 = "BDconfidenthit_protein_".$prefix1."_pcnt-".$pep_pcnt_th."_pcov-".$pep_pcov_th."_Eval-".$pep_eval_th.".csv";
	SAVE($rfile2, $pr1);
	
	#------------------------------------------------------//
	print "\n! Searching for partner that has the same genome position...\n";
	my @ID1h = keys(%{$BDhits});
	@ID1h = sort {$a cmp $b} @ID1h;
	
	my $summary_partner = "query\torthologue\tanother gene as a partner of query\tnr-blastn-hit\tnr-blastp-hit\tremark\tscore\n";
	my $cnt_partner1 = 0;
	my $cnt_partner2 = 0;
	my $cnt_partner3 = 0;
	my $cnt_partner4 = 0;
	my $cnt_nopartner1 = 0;
	my $cnt_nopartner2 = 0;
	my $cnt_nopartner3 = 0;
	my $cnt_nopartner4 = 0;
	my $hscore_tid_compare = {};
	my $hexonmatch_tid_compare = {};
	
	foreach my $hid1 (@ID1){
		if($BDhitsAll->{$hid1}){					#$BDhitsAll -> both blastn and p
			my $subh = $BDhitsAll->{$hid1};
			my @ID2h = keys(%{$subh});
			@ID2h = sort {$a cmp $b} @ID2h;
			
			my $Aopt = [];
			for(my $i = 0; $i < 9; $i++){
				my $rth = 1 - (0.1 * $i);
				
				if($i == 0){
					$rth = 0.99;
				}
				
				foreach my $hid2 (@ID2h){
		#			$B[0] = $A[3];
		#			$B[1] = $A[4];
		#			$B[2] = $tid;
		#			$B[3] = $A[6];		#strand
		#			$B[4] = $A[0];		#chromosome
					
					my $AoLE1 = $hexon1->{$hid1};
					my $AoLE2 = $hexon2->{$hid2};
					
					unless($AoLE2->[0][4]){			#some gene may not be aligned to reference by blat, skip such case
						next;
					}
					
					my $score = 0;
					if($BDhits->{$hid1}{$hid2} && $pBDhits->{$hid1}{$hid2}){
						$score = $BDhits->{$hid1}{$hid2} + $pBDhits->{$hid1}{$hid2} * 10;
					}
					elsif($BDhits->{$hid1}{$hid2} && ! $pBDhits->{$hid1}{$hid2}){
						$score = $BDhits->{$hid1}{$hid2};
					}
					elsif(! $BDhits->{$hid1}{$hid2} && $pBDhits->{$hid1}{$hid2}){
						$score = $pBDhits->{$hid1}{$hid2};
					}
					
					$hscore_tid_compare->{$hid1}{$hid2} = $score;
					
					if($AoLE1->[0][4] eq $AoLE2->[0][4] && $AoLE1->[0][3] eq $AoLE2->[0][3]){		#if both chromosome and strand are identical...
						my $compdata = Compare_exon_structure($hid1, $hid2, $AoLE1, $AoLE2, $AoLE1->[0][3], $AoLE2->[0][3], 400, 500, $rth);
						
						if($compdata->{judge} eq 'true'){
							if($compdata->{score}){
								$score = $score * $compdata->{score} * 10;
							}
							
							$hscore_tid_compare->{$hid1}{$hid2} = $score;
							$hexonmatch_tid_compare->{$hid1}{$hid2}{end_exon_match} = $compdata->{end_exon_match};
							$hexonmatch_tid_compare->{$hid1}{$hid2}{scoret} = $compdata->{scoret};
							$hexonmatch_tid_compare->{$hid1}{$hid2}{scoreh} = $compdata->{scoreh};
							$hexonmatch_tid_compare->{$hid1}{$hid2}{num_hit} = $compdata->{num_hit};
							$hexonmatch_tid_compare->{$hid1}{$hid2}{numt} = $compdata->{numt};
							$hexonmatch_tid_compare->{$hid1}{$hid2}{numh} = $compdata->{numh};
							
							my @pt;
							push(@pt, $hid2);
							push(@pt, $score);
							push(@{$Aopt}, \@pt);
						}
					}
				}
			}
			@{$Aopt} = sort {$b->[1] <=> $a->[1]} @{$Aopt};
			
			my $partner;
			if(@{$Aopt}){
				$partner = $Aopt->[0][0];
			}
			
			unless($nrhits->{$hid1}){
				$nrhits->{$hid1} = "-";
			}
			unless($nrphits->{$hid1}){
				$nrphits->{$hid1} = "-";
			}
			
			$nrhits->{$hid1} =~ s/grp\=//g;
			
			if($partner){
				$partner =~ s/grp\=//;
				if($BDhits->{$hid1}{$partner} && $pBDhits->{$hid1}{$partner}){
					my $score = $BDhits->{$hid1}{$partner} + $pBDhits->{$hid1}{$partner} * 10;
					$summary_partner .= $hid1."\t".$partner."\t-\t".$nrhits->{$hid1}."\t".$nrphits->{$hid1}."\talso_blastnp_BDH\t".$score."\n";
					$cnt_partner1++;
				}
				elsif($BDhits->{$hid1}{$partner} && ! $pBDhits->{$hid1}{$partner}){
					my $score = $BDhits->{$hid1}{$partner};
					$summary_partner .= $hid1."\t".$partner."\t-\t".$nrhits->{$hid1}."\t".$nrphits->{$hid1}."\tpartner_is_blastn_BDH_but_not_blastp_BDH\t".$score."\n";
					$cnt_partner2++;
				}
				elsif(! $BDhits->{$hid1}{$partner} && $pBDhits->{$hid1}{$partner}){
					my $score = $pBDhits->{$hid1}{$partner};
					$summary_partner .= $hid1."\t".$partner."\t-\t".$nrhits->{$hid1}."\t".$nrphits->{$hid1}."\tpartner_is_blastp_BDH_but_not_blastn_BDH\t".$score."\n";
					$cnt_partner3++;
				}
				else{
					$cnt_partner4++;
				}
			}
			else{
				my $same_nphit = Search_same_hit_np($nrhits->{$hid1}, $nrphits->{$hid1});
				
				if($same_nphit ne '-'){
					$summary_partner .= $hid1."\tNA\t".$same_nphit."\t".$nrhits->{$hid1}."\t".$nrphits->{$hid1}."\tno_partner_but_BDH_is_same_in_np\t0\n";
#					$summary_partner .= $hid1."\t".$same_nphit."\t-\t".$nrhits->{$hid1}."\t".$nrphits->{$hid1}."\tno_partner_but_BDH_is_same_in_np\t0\n";
					$cnt_nopartner1++;
				}
				else{
					$summary_partner .= $hid1."\tNA\t-\t".$nrhits->{$hid1}."\t".$nrphits->{$hid1}."\tno_partner_and_BDH_is_distinct_in_np\t0\n";
					$cnt_nopartner2++;
				}
			}
		}
		else{
			$summary_partner .= $hid1."\tNA\t-\t-\t-\t-\t0\n";
			$cnt_nopartner4++;
		}
	}
	
	print "! [$cnt_partner1] have blastn/p BDH with same genomic position\n";
	print "! [$cnt_partner2] have blastn BDH with same genomic position but no blastp BDH\n";
	print "! [$cnt_partner3] have blastp BDH with same genomic position but no blastn BDH\n";
	print "! [$cnt_partner4] (other)\n";
	print "! [$cnt_nopartner1] have no partner but have identical BDH in blastn/p\n";
	print "! [$cnt_nopartner2] have no partner and have distinct blastn/p BDH\n";
	print "! [$cnt_nopartner4] have no partner\n";
	
#	my $rfile3 = "list_partner_ID.tsv";
	my $rfile3 = "list_partner_ID"."_c".$th_condition.".tsv";
	SAVE($rfile3, $summary_partner);
	
	#------------------------------------------------------//
	print "\n! Gene-based search for partner...\n";
	my @SPT = split(/\n/, $summary_partner);
	my $hSPT = {};
	my $rSPT = {};
	foreach my $line (@SPT){
		if($line =~ /query\torthologue\t/){
			next;
		}
		my @A = split(/\t/, $line);
		
#		my $gid = Tid2gid($A[0]);
#		my $hgid = Tid2gid($A[1]);
		
		my $gid;
		if($A[0] ne 'NA' && ! defined $htid2gid1->{$A[0]}){
			print "! missing gid for [$A[0]] (a1)\n";
			print "! line=[$line]\n";
			die;
		}
		else{
			$gid = $htid2gid1->{$A[0]};
		}
		
		my $hgid;
		if($A[1] eq 'NA'){
			$hgid = "NA";
		}
		else{
			if(! defined $htid2gid2ori->{$A[1]}){
				print "! missing gid for [$A[1]] (a2)\n";
				print "! line=[$line]\n";
				die;
			}
			else{
				$hgid = $htid2gid2ori->{$A[1]};
			}
		}
		
		$hSPT->{$gid}{str} .= $line."\n";
		
		unless($hSPT->{$gid}{partner}){
			if($hgid ne 'NA'){
				$hSPT->{$gid}{partner} = $hgid;
				$hSPT->{$gid}{$hgid}{info} = $A[5];
				$hSPT->{$gid}{$hgid}{score} = $A[6];
				$hSPT->{$gid}{best_info} = $A[5];
				$hSPT->{$gid}{best_score} = $A[6];
				$hSPT->{$gid}{cnt} += 1;
				
				if($A[5] eq 'also_blastnp_BDH'){
					if($rSPT->{$hgid}{partner}){
						unless($rSPT->{$hgid}{partner} =~ /$gid/){
							$rSPT->{$hgid}{partner} .= $gid.",";
							$rSPT->{$hgid}{cnt} += 1;
						}
					}
					else{
						$rSPT->{$hgid}{partner} .= $gid.",";
						$rSPT->{$hgid}{cnt} += 1;
					}
				}
			}
			else{
				$hSPT->{$gid}{partner} = $hgid;
			}
		}
		else{
			if($hgid ne 'NA'){
				if($hSPT->{$gid}{partner} eq 'NA'){
					$hSPT->{$gid}{partner} = $hgid;
					$hSPT->{$gid}{$hgid}{info} = $A[5];
					$hSPT->{$gid}{$hgid}{score} = $A[6];
					$hSPT->{$gid}{best_info} = $A[5];
					$hSPT->{$gid}{best_score} = $A[6];
					$hSPT->{$gid}{cnt} += 1;
					
#					$rSPT->{$hgid}{partner} .= $gid.",";
#					$rSPT->{$hgid}{cnt} += 1;
					
					if($A[5] eq 'also_blastnp_BDH'){
						if($rSPT->{$hgid}{partner}){
							unless($rSPT->{$hgid}{partner} =~ /$gid/){
								$rSPT->{$hgid}{partner} .= $gid.",";
								$rSPT->{$hgid}{cnt} += 1;
							}
						}
						else{
							$rSPT->{$hgid}{partner} .= $gid.",";
							$rSPT->{$hgid}{cnt} += 1;
						}
					}
				}
				elsif($hSPT->{$gid}{partner} !~ $hgid){
					if($hSPT->{$gid}{best_info} eq 'also_blastnp_BDH'){
						if($A[5] eq 'also_blastnp_BDH'){
							if($hSPT->{$gid}{partner}){
								unless($hSPT->{$gid}{partner} =~ /$hgid/){
									$hSPT->{$gid}{partner} .= ",".$hgid;
									$hSPT->{$gid}{cnt} += 1;
								}
							}
							else{
								$hSPT->{$gid}{partner} .= ",".$hgid;
								$hSPT->{$gid}{cnt} += 1;
							}
							
							if($rSPT->{$hgid}{partner}){
								unless($rSPT->{$hgid}{partner} =~ /$gid/){
									$rSPT->{$hgid}{partner} .= $gid.",";
									$rSPT->{$hgid}{cnt} += 1;
								}
							}
							else{
								$rSPT->{$hgid}{partner} .= $gid.",";
								$rSPT->{$hgid}{cnt} += 1;
							}
							
							$hSPT->{$gid}{$hgid}{info} .= ",".$A[5];
							$hSPT->{$gid}{$hgid}{score} += $A[6];
							
							if($hSPT->{$gid}{best_score} < $A[6]){
								$hSPT->{$gid}{best_info} = $A[5];
								$hSPT->{$gid}{best_score} = $A[6];
							}
						}
					}
					elsif($hSPT->{$gid}{best_info} eq 'partner_is_blastn_BDH_but_not_blastp_BDH'){
						if($A[5] eq 'also_blastnp_BDH'){
							$hSPT->{$gid}{partner} = $hgid;
							$hSPT->{$gid}{$hgid}{info} = $A[5];
							$hSPT->{$gid}{$hgid}{score} = $A[6];
							$hSPT->{$gid}{cnt} += 1;
							
							$hSPT->{$gid}{best_info} = $A[5];
							$hSPT->{$gid}{best_score} = $A[6];
							
							if($rSPT->{$hgid}{partner}){
								$rSPT = Erase_entry($rSPT, $hSPT->{$gid}{partner}, $gid);		#$rSPT, previous_hgid that should edit, gid
								$rSPT->{$hgid}{partner} .= $gid.",";
								$rSPT->{$hgid}{cnt} += 1;
							}
							else{
								$rSPT->{$hgid}{partner} .= $gid.",";
								$rSPT->{$hgid}{cnt} += 1;
							}
						}
						elsif($A[5] eq 'partner_is_blastn_BDH_but_not_blastp_BDH'){
							if($hSPT->{$gid}{partner}){
								unless($hSPT->{$gid}{partner} =~ /$hgid/){
									$hSPT->{$gid}{partner} .= ",".$hgid;
									$hSPT->{$gid}{cnt} += 1;
								}
							}
							else{
								$hSPT->{$gid}{partner} .= ",".$hgid;
								$hSPT->{$gid}{cnt} += 1;
							}
							
#							if($rSPT->{$hgid}{partner}){
#								unless($rSPT->{$hgid}{partner} =~ /$gid/){
#									$rSPT->{$hgid}{partner} .= $gid.",";
#									$rSPT->{$hgid}{cnt} += 1;
#								}
#							}
#							else{
#								$rSPT->{$hgid}{partner} .= $gid.",";
#								$rSPT->{$hgid}{cnt} += 1;
#							}
							
							$hSPT->{$gid}{$hgid}{info} .= ",".$A[5];
							$hSPT->{$gid}{$hgid}{score} += $A[6];
							
							if($hSPT->{$gid}{best_score} < $A[6]){
								$hSPT->{$gid}{best_info} = $A[5];
								$hSPT->{$gid}{best_score} = $A[6];
							}
						}
					}
					elsif($hSPT->{$gid}{best_info} eq 'partner_is_blastp_BDH_but_not_blastn_BDH'){
#						if($A[5] eq 'also_blastnp_BDH' || $A[5] eq 'partner_is_blastn_BDH_but_not_blastp_BDH'){
						if($A[5] eq 'also_blastnp_BDH'){
							$hSPT->{$gid}{partner} = $hgid;
							$hSPT->{$gid}{$hgid}{info} = $A[5];
							$hSPT->{$gid}{$hgid}{score} = $A[6];
							$hSPT->{$gid}{cnt} += 1;
							
							$hSPT->{$gid}{best_info} = $A[5];
							$hSPT->{$gid}{best_score} = $A[6];
							
							if($rSPT->{$hgid}{partner}){
								$rSPT = Erase_entry($rSPT, $hSPT->{$gid}{partner}, $gid);		#$rSPT, previous_hgid that should edit, gid
								$rSPT->{$hgid}{partner} .= $gid.",";
								$rSPT->{$hgid}{cnt} += 1;
							}
							else{
								$rSPT->{$hgid}{partner} .= $gid.",";
								$rSPT->{$hgid}{cnt} += 1;
							}
						}
						elsif($A[5] eq 'partner_is_blastp_BDH_but_not_blastn_BDH'){
							if($hSPT->{$gid}{partner}){
								unless($hSPT->{$gid}{partner} =~ /$hgid/){
									$hSPT->{$gid}{partner} .= ",".$hgid;
									$hSPT->{$gid}{cnt} += 1;
								}
							}
							else{
								$hSPT->{$gid}{partner} .= ",".$hgid;
								$hSPT->{$gid}{cnt} += 1;
							}
							
#							if($rSPT->{$hgid}{partner}){
#								unless($rSPT->{$hgid}{partner} =~ /$gid/){
#									$rSPT->{$hgid}{partner} .= $gid.",";
#									$rSPT->{$hgid}{cnt} += 1;
#								}
#							}
#							else{
#								$rSPT->{$hgid}{partner} .= $gid.",";
#								$rSPT->{$hgid}{cnt} += 1;
#							}
							
							$hSPT->{$gid}{$hgid}{info} .= ",".$A[5];
							$hSPT->{$gid}{$hgid}{score} += $A[6];
							
							if($hSPT->{$gid}{best_score} < $A[6]){
								$hSPT->{$gid}{best_info} = $A[5];
								$hSPT->{$gid}{best_score} = $A[6];
							}
						}
					}
				}
			}
		}
	}
	
	my $cntg_partner1 = 0;
	my $cntg_partner11 = 0;
	my $cntg_partner12 = 0;
	my $cntg_partner13 = 0;
	my $cntg_partner14 = 0;
	my $cntg_partner15 = 0;
	my $cntg_partner1b = 0;
	my $cntg_partner11b = 0;
	my $cntg_partner12b = 0;
	my $cntg_partner13b = 0;
	my $cntg_partner14b = 0;
	my $cntg_partner15b = 0;
	my $cntg_partner1c = 0;
	my $cntg_partner11c = 0;
	my $cntg_partner12c = 0;
	my $cntg_partner13c = 0;
	my $cntg_partner14c = 0;
	my $cntg_partner15c = 0;
	my $cntg_partner2c = 0;
	my $cntg_partner21c = 0;
	my $cntg_partner22c = 0;
	my $cntg_partner23c = 0;
	my $cntg_partner24c = 0;
	my $cntg_partner25c = 0;
	my $cntg_partner3 = 0;
	my $cntg_partner4 = 0;
	my $cntg_nopartner = 0;
	my $cntg_total = 0;
	
	my $unmapped_qexon;
	my $unmapped_qprot;
	
	my $nr_2gid = {};
	my $selected_gid = {};
	my $seqID_cnv;
	my $seqID_nohit;
	my $max_nss1 = 1;
	
	my $genebase_summary = "gene ID\tprotein coding\tcandidate partner\tpartner (genomic position + blast n/p)\tpartner (genomic position + blastn or p)\tblast n/p BDH but different genomic position\tBDH has another partner\tstatus summary1\tstatus summary2\n";
	my $genebase_posinfo;
	
	foreach my $line (@SPT){
		if($line =~ /query\torthologue\t/){
			next;
		}
		my @A = split(/\t/, $line);
		
#		my $gid = Tid2gid($A[0]);
		my $gid = $htid2gid1->{$A[0]};
		
		if($selected_gid->{$gid}){
			next;
		}
		
		my @TID1 = split(/\n/, $hgid2tid1->{$gid});
		my $protein_coding = "FALSE";
		foreach my $tid1 (@TID1){
			if($CDScheck_info->{$tid1}){
				if($CDScheck_info->{$tid1} eq 'protein coding'){
					$protein_coding = "TRUE";
					last;
				}
			}
			else{
				print "! error : missing protein coding info for [$tid1]\n";
				sleep(1);
			}
		}
		
		if($hSPT->{$gid}{partner} ne 'NA'){
			my $partner = $hSPT->{$gid}{partner};
			
			my $nss1 = 0;
			my $nss2 = 0;
			my $nss3 = 0;
			my @SS = split(/,/, $partner);
			my $AoP = [];
			my $Aost = [];
			
			my $AoSS1 = [];
			foreach my $ssid (@SS){
				if($ssid){
					if(defined $hgff2ori->{gid2ch}{$ssid}){
						my @tss;
						unless($hgff2ori->{gid2ch}{$ssid} =~ /chr00/){
							unless(defined $hgff2ori->{gid2ch}{$ssid}){
								print "! missing chr info for [$ssid] (1)\n";
								die;
							}
							unless(defined $hgff1->{gid2ch}{$gid}){
								print "! missing chr info for [$gid] (2)\n";
								die;
							}
							
							if($hgff2ori->{gid2ch}{$ssid} eq $hgff1->{gid2ch}{$gid}){
								push(@tss, $ssid);
								push(@tss, 0);
							}
							else{
								push(@tss, $ssid);
								push(@tss, 1);
							}
						}
						else{
							push(@tss, $ssid);
							push(@tss, 2);
						}
						push(@{$AoSS1}, \@tss);
					}
					else{
						print "! missing chr ID info : [$ssid], abort script...\n";
						die;
					}
				}
			}
			@{$AoSS1} = sort {$a->[1] <=> $b->[1]} @{$AoSS1};
			
			my @SS1r;
			foreach my $tss (@{$AoSS1}){
				if($tss->[0]){
					my $ssid = $tss->[0];
					
					if($nss1 >= 1){
						if($hgff2ori->{gid2ch}{$ssid}){
							if($hgff2ori->{gid2ch}{$ssid} =~ /chr00/ || $hgff2ori->{gid2ch}{$ssid} eq 'chr0'){
								next;
							}
						}
					}
					
					push(@SS1r, $ssid);
					$nss1++;
					
					if($rSPT->{$ssid}{partner}){
						my @SS2 = split(/,/, $rSPT->{$ssid}{partner});
						foreach my $ssid2 (@SS2){
							if($ssid2){
								if($hSPT->{$ssid2}{$ssid}{score}){
									my $sub_jscore = Judge_score($hSPT, $ssid2, $ssid);
									
									if($nss2 == 0){
										my @stmp;
										push(@stmp, $ssid2);
										push(@stmp, $sub_jscore);
										push(@{$Aost}, \@stmp);
										$nss3++;
									}
									elsif($Aost->[0][1] == $sub_jscore){
										my @stmp;
										push(@stmp, $ssid2);
										push(@stmp, $sub_jscore);
										push(@{$Aost}, \@stmp);
										$nss3++;
									}
									
									$nss2++;
								}
							}
						}
					}
					
					my $judge_score = Judge_score($hSPT, $gid, $ssid);
#					if($hSPT->{$gid}{$ssid}{info} =~ /also_blastnp_BDH/){
#						$judge_score = 10;
#					}
#					elsif($hSPT->{$gid}{$ssid}{info} =~ /partner_is_blastn_BDH_but_not_blastp_BDH/ && $hSPT->{$gid}{$ssid}{info} =~ /partner_is_blastp_BDH_but_not_blastn_BDH/){
#						$judge_score = 5;
#					}
#					elsif($hSPT->{$gid}{$ssid}{info} =~ /partner_is_blastp_BDH_but_not_blastn_BDH/){
#						$judge_score = 2;
#					}
#					elsif($hSPT->{$gid}{$ssid}{info} =~ /partner_is_blastn_BDH_but_not_blastp_BDH/){
#						$judge_score = 1;
#					}
					
					my @P;
					push(@P, $ssid);
					push(@P, $judge_score);
					push(@{$AoP}, \@P);
				}
			}
			
			if($nss1 > $max_nss1){
				$max_nss1 = $nss1;
			}
			
			my $npartner = join(",", @SS1r);
			@{$AoP} = sort {$b->[1] <=> $a->[1]} @{$AoP};
			
			if($nss1 > 1){
				if($AoP->[0][1] == 10){
					$genebase_summary .= $gid."\t".$protein_coding."\t-\t".$partner."\t-\t-\t-\tblastn/p + same genomic position\t";
					$genebase_summary .= "gid has >=2 partners";
					$genebase_summary .= "\n";
					$cntg_partner11b++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\tgid has >=2 partners";
					foreach my $subpartner (@SS1r){
						$genebase_posinfo .= "\t".$subpartner."\t".$hgff2ori->{gid2ch}{$subpartner}."\t".$hgff2ori->{gid2pos0}{$subpartner}."\t".$hgff2ori->{gid2pos1}{$subpartner}."\t".$hgff2ori->{gid2strand}{$subpartner}."\t".$hgff2ori->{gid2gnum}{$subpartner};
						$genebase_posinfo .= "\t".SameChr_or_Not($hgff1->{gid2ch}{$gid}, $hgff2ori->{gid2ch}{$subpartner});		#return string : 'same_chr' or 'diff_chr'
					}
					$genebase_posinfo .= "\n";
				}
				elsif($AoP->[0][1] == 5){
					$genebase_summary .= $gid."\t".$protein_coding."\t-\t".$partner."\t-\t-\t-\tblastn/p in distinct variants + same genomic position\t";
					$genebase_summary .= "gid has >=2 partners";
					$genebase_summary .= "\n";
					$cntg_partner12b++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\tgid has >=2 partners";
					foreach my $subpartner (@SS1r){
						$genebase_posinfo .= "\t".$subpartner."\t".$hgff2ori->{gid2ch}{$subpartner}."\t".$hgff2ori->{gid2pos0}{$subpartner}."\t".$hgff2ori->{gid2pos1}{$subpartner}."\t".$hgff2ori->{gid2strand}{$subpartner}."\t".$hgff2ori->{gid2gnum}{$subpartner};
						$genebase_posinfo .= "\t".SameChr_or_Not($hgff1->{gid2ch}{$gid}, $hgff2ori->{gid2ch}{$subpartner});		#return string : 'same_chr' or 'diff_chr'
					}
					$genebase_posinfo .= "\n";
				}
				elsif($AoP->[0][1] == 2){
					$genebase_summary .= $gid."\t".$protein_coding."\t-\t-\t".$partner."\t-\t-\tblastp + same genomic position\t";
					$genebase_summary .= "gid has >=2 partners";
					$genebase_summary .= "\n";
					$cntg_partner13b++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\t-";
					$genebase_posinfo .= "\t-\t-\t-\t-\tNA\t-\t-\n";
				}
				elsif($AoP->[0][1] == 1){
					$genebase_summary .= $gid."\t".$protein_coding."\t-\t-\t".$partner."\t-\t-\tblastn + same genomic position\t";
					$genebase_summary .= "gid has >=2 partners";
					$genebase_summary .= "\n";
					$cntg_partner14b++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\t-";
					$genebase_posinfo .= "\t-\t-\t-\t-\tNA\t-\t-\n";
				}
				else{
					$cntg_partner15b++;
				}
				
				$cntg_partner1b++;
			}
			elsif($nss2 > 1){
				if($Aost->[0][1] == 10){
					$genebase_summary .= $gid."\t".$protein_coding."\t-\t".$partner."\t-\t-\t".$rSPT->{$partner}{partner}."\tblastn/p + same genomic position\t";
					$genebase_summary .= "partner has >=2 gids";
					$genebase_summary .= "\n";
					$cntg_partner21c++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\tpartner has >=2 gids";
					foreach my $subpartner (@SS1r){
						$genebase_posinfo .= "\t".$subpartner."\t".$hgff2ori->{gid2ch}{$subpartner}."\t".$hgff2ori->{gid2pos0}{$subpartner}."\t".$hgff2ori->{gid2pos1}{$subpartner}."\t".$hgff2ori->{gid2strand}{$subpartner}."\t".$hgff2ori->{gid2gnum}{$subpartner};
						$genebase_posinfo .= "\t".SameChr_or_Not($hgff1->{gid2ch}{$gid}, $hgff2ori->{gid2ch}{$subpartner});		#return string : 'same_chr' or 'diff_chr'
					}
					$genebase_posinfo .= "\n";
				}
				elsif($Aost->[0][1] == 5){
					$genebase_summary .= $gid."\t".$protein_coding."\t-\t".$partner."\t-\t-\t".$rSPT->{$partner}{partner}."\tblastn/p in distinct variants + same genomic position\t";
					$genebase_summary .= "partner has >=2 gids";
					$genebase_summary .= "\n";
					$cntg_partner22c++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\tpartner has >=2 gids";
					foreach my $subpartner (@SS1r){
						$genebase_posinfo .= "\t".$subpartner."\t".$hgff2ori->{gid2ch}{$subpartner}."\t".$hgff2ori->{gid2pos0}{$subpartner}."\t".$hgff2ori->{gid2pos1}{$subpartner}."\t".$hgff2ori->{gid2strand}{$subpartner}."\t".$hgff2ori->{gid2gnum}{$subpartner};
						$genebase_posinfo .= "\t".SameChr_or_Not($hgff1->{gid2ch}{$gid}, $hgff2ori->{gid2ch}{$subpartner});		#return string : 'same_chr' or 'diff_chr'
					}
					$genebase_posinfo .= "\n";
				}
				elsif($Aost->[0][1] == 2){
					$genebase_summary .= $gid."\t".$protein_coding."\t-\t".$partner."\t-\t-\t".$rSPT->{$partner}{partner}."\tblastp + same genomic position\t";
					$genebase_summary .= "partner has >=2 gids";
					$genebase_summary .= "\n";
					$cntg_partner23c++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\t-";
					$genebase_posinfo .= "\t-\t-\t-\t-\tNA\t-\t-\n";
				}
				elsif($Aost->[0][1] == 1){
					$genebase_summary .= $gid."\t".$protein_coding."\t-\t".$partner."\t-\t-\t".$rSPT->{$partner}{partner}."\tblastn + same genomic position\t";
					$genebase_summary .= "partner has >=2 gids";
					$genebase_summary .= "\n";
					$cntg_partner24c++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\t-";
					$genebase_posinfo .= "\t-\t-\t-\t-\tNA\t-\t-\n";
				}
				else{
					$cntg_partner25c++;
				}
				
				$cntg_partner2c++;
				$nr_2gid->{$partner} = 1;
			}
			else{
				if($AoP->[0][1] == 10){
					$genebase_summary .= $gid."\t".$protein_coding."\t".$npartner."\t".$partner."\t-\t-\t-\tblastn/p + same genomic position\t";
					$genebase_summary .= "1-to-1 ortholog";
					$genebase_summary .= "\n";
					$cntg_partner11++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\tprimary";
					foreach my $subpartner (@SS1r){
						$genebase_posinfo .= "\t".$subpartner."\t".$hgff2ori->{gid2ch}{$subpartner}."\t".$hgff2ori->{gid2pos0}{$subpartner}."\t".$hgff2ori->{gid2pos1}{$subpartner}."\t".$hgff2ori->{gid2strand}{$subpartner}."\t".$hgff2ori->{gid2gnum}{$subpartner};
						$genebase_posinfo .= "\t".SameChr_or_Not($hgff1->{gid2ch}{$gid}, $hgff2ori->{gid2ch}{$subpartner});		#return string : 'same_chr' or 'diff_chr'
					}
					$genebase_posinfo .= "\n";
				}
				elsif($AoP->[0][1] == 5){
					$genebase_summary .= $gid."\t".$protein_coding."\t".$npartner."\t".$partner."\t-\t-\t-\tblastn/p in distinct variants + same genomic position\t";
					$genebase_summary .= "-";
					$genebase_summary .= "\n";
					$cntg_partner12++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\tprimary";
					foreach my $subpartner (@SS1r){
						$genebase_posinfo .= "\t".$subpartner."\t".$hgff2ori->{gid2ch}{$subpartner}."\t".$hgff2ori->{gid2pos0}{$subpartner}."\t".$hgff2ori->{gid2pos1}{$subpartner}."\t".$hgff2ori->{gid2strand}{$subpartner}."\t".$hgff2ori->{gid2gnum}{$subpartner};
						$genebase_posinfo .= "\t".SameChr_or_Not($hgff1->{gid2ch}{$gid}, $hgff2ori->{gid2ch}{$subpartner});		#return string : 'same_chr' or 'diff_chr'
					}
					$genebase_posinfo .= "\n";
				}
				elsif($AoP->[0][1] == 2){
					$genebase_summary .= $gid."\t".$protein_coding."\t-\t-\t".$partner."\t-\t-\tblastp + same genomic position\t";
					$genebase_summary .= "-";
					$genebase_summary .= "\n";
					$cntg_partner13++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\t-";
					$genebase_posinfo .= "\t-\t-\t-\t-\tNA\t-\t-\n";
				}
				elsif($AoP->[0][1] == 1){
					$genebase_summary .= $gid."\t".$protein_coding."\t-\t-\t".$partner."\t-\t-\tblastn + same genomic position\t";
					$genebase_summary .= "-";
					$genebase_summary .= "\n";
					$cntg_partner14++;
					
					$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\t-";
					$genebase_posinfo .= "\t-\t-\t-\t-\tNA\t-\t-\n";
				}
				else{
					$cntg_partner15++;
				}
				
				$cntg_partner1++;
			}
			
			$cntg_total++;
		}
		elsif($A[2] ne '-'){
			my @htmp = split(/\./, $A[2]);
			my $num_htmp = @htmp;
			
			my $hgid2;
			if($num_htmp == 1){
				$hgid2 = $A[2];
			}
			else{
				for(my $i = 1; $i <= $num_htmp; $i++){
					my @tmpid;
					for(my $j = 0; $j < $i; $j++){
						push(@tmpid, $htmp[$j]);
					}
					$hgid2 = join("\.", @tmpid);
					
					if($rSPT->{$hgid2}{partner}){
						last;
					}
				}
			}
			
			my $subsw = 0;
			if($rSPT->{$hgid2}{partner}){
				my $cnt_partner = Count_rSPT_partner($rSPT->{$hgid2}{partner});
				
				if($cnt_partner == 1){
					$subsw = 1;
				}
			}
			
			if($subsw eq '1'){
				$genebase_summary .= $gid."\t".$protein_coding."\t-\t-\t-\t".$hgid2."\t".$rSPT->{$hgid2}{partner}."\tblastn/p but different genomic position\t";
				$genebase_summary .= "another gid as BDH pair";
				$genebase_summary .= "\n";
				$cntg_partner3++;
				$cntg_total++;
				
				$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\t-";
				$genebase_posinfo .= "\t-\t-\t-\t-\tNA\t-\t-\n";
			}
			else{
				$genebase_summary .= $gid."\t".$protein_coding."\t-\t-\t-\t".$hgid2."\t-\tblastn/p but different genomic position\t";
				$genebase_summary .= "other >2 gids as BDH pair";
				$genebase_summary .= "\n";
				$cntg_partner4++;
				$cntg_total++;
				
				$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\t-";
				$genebase_posinfo .= "\t-\t-\t-\t-\tNA\t-\t-\n";
			}
		}
		else{
			$genebase_summary .= $gid."\t".$protein_coding."\t-\t-\t-\t-\t-\t-\t";
			$genebase_summary .= "-";
			$genebase_summary .= "\n";
			$cntg_nopartner++;
			$cntg_total++;
			
			$genebase_posinfo .= $gid."\t".$hgff1->{gid2ch}{$gid}."\t".$hgff1->{gid2pos0}{$gid}."\t".$hgff1->{gid2pos1}{$gid}."\t".$hgff1->{gid2strand}{$gid}."\t".$hgff1->{gid2gnum}{$gid}."\t-";
			$genebase_posinfo .= "\t-\t-\t-\t-\tNA\t-\t-\n";
			
			if($hgid2tid1->{$gid}){
				my @eachTID = split(/\n/, $hgid2tid1->{$gid});
				foreach my $tid (@eachTID){
					$unmapped_qexon .= $hseq_qexon->{$tid}{seq};
					$unmapped_qprot .= $hseq_qprot->{$tid}{seq};
				}
			}
		}
		
		$selected_gid->{$gid} = 1;
	}
	
	#-----------------------------------------------------//
	my $Chars = Gen_AoChars();
#	my $RGS = Revise_genebase_summary($genebase_summary, $Chars);
	
#	$genebase_summary = $RGS->{r};
#	my $horth = $RGS->{horth};
#	my $hanot = $RGS->{hanot};
#	my $hp2g = $RGS->{hp2g};
#	my $hg2p = $RGS->{hg2p};
#	my $hNA = $RGS->{hNA};
	
#	$rh->{hCDS} = $hCDS;
#	$rh->{hexon} = $hexon;
#	$rh->{htidall} = $htidall;
#	$rh->{gene} = $gene;
#	$rh->{gid2ch} = $gid2ch;
#	$rh->{gid2pos0} = $gid2pos0;
#	$rh->{gid2pos1} = $gid2pos1;
#	$rh->{gid2strand} = $gid2strand;
#	$rh->{gid2gnum} = $gid2gnum;
#	my $hgnum1 = $hgff1->{gid2gnum};
#	my $hgnum2 = $hgff2ori->{gid2gnum};
	
	#In case of judgement = "gid has >=2 partners", analyze whether partners are paralogs for gid or merely splitted (or partial) genes
	my @GSL = split(/\n/, $genebase_summary);
	my $n_genebase_summary = shift(@GSL)."\tparalogs\ttandem paralogs\tCNV\tCNV (target vs hint)\tCNV num\tcandidate for consensus ID\tAlt-ID in replace.csv";
	$n_genebase_summary .= "\thomolog remark";
	$n_genebase_summary .= "\texon match status (11=1st-last, 10=1st-internal, 1=internal-last, 0=internal)\tnum exon\tnum matched exon\tratio matched exon";
	$n_genebase_summary .= "\t\%identity (t)\t\%alignment coverage (t)\tE-value (t)\tsum_score (t)\t\%identity (p)\t\%alignment coverage (p)\tE-value (p)\tsum_score (p)\n";
	my $hcnt_hid = Assign_and_Check_redundancy(\@GSL);
	
	foreach my $l (@GSL){
		my @A = split(/\t/, $l);
		if($A[8] eq 'gid has >=2 partners'){					#one query gene has >=2 hits in the genomic position
			unless($A[3] eq '-'){
				my @H = split(/,/, $A[3]);
				
				my $PRLG = {};
				my $PRLGnear = {};
				foreach my $h1 (@H){
					my @TID1 = split(/\n/, $hgid2tid2ori->{$h1});
					
					foreach my $h2 (@H){
						if($h1 eq $h2){
							next;
						}
						
						my $h1chr = $hgff2ori->{gid2ch}{$h1};
						my $h2chr = $hgff2ori->{gid2ch}{$h2};
						my $h1pos0 = $hgff2ori->{gid2pos0}{$h1};
						my $h2pos0 = $hgff2ori->{gid2pos0}{$h2};
						my $h1pos1 = $hgff2ori->{gid2pos1}{$h1};
						my $h2pos1 = $hgff2ori->{gid2pos1}{$h2};
						my $h1gnum = $hgff2ori->{gid2gnum}{$h1};
						my $h2gnum = $hgff2ori->{gid2gnum}{$h2};
						
						my $dist_gnum = abs($h1gnum - $h2gnum);
						
						my $dist_bp = 10000000000000;
						if($h1chr eq $h2chr){
							if($h1pos1 < $h2pos0){
								$dist_bp = abs($h2pos0 - $h1pos1);
							}
							elsif($h1pos0 > $h2pos1){
								$dist_bp = abs($h1pos0 - $h2pos1);
							}
							elsif($h1pos0 <= $h2pos0 && $h2pos0 <= $h1pos1){
								$dist_bp = 0;
							}
							elsif($h2pos0 <= $h1pos0 && $h1pos0 <= $h2pos1){
								$dist_bp = 0;
							}
							elsif($h1pos0 <= $h2pos1 && $h2pos1 <= $h1pos1){
								$dist_bp = 0;
							}
							elsif($h2pos0 <= $h1pos1 && $h1pos1 <= $h2pos1){
								$dist_bp = 0;
							}
						}
						
						my @TID2 = split(/\n/, $hgid2tid2ori->{$h2});
						foreach my $ht1 (@TID1){
							my $AoLE1 = $hexon2->{$ht1};
							
							foreach my $ht2 (@TID2){
								my $AoLE2 = $hexon2->{$ht2};
								
								if($AoLE1->[0][3] eq $AoLE2->[0][3]){
									my $compdata = Compare_exon_structure($ht1, $ht2, $AoLE1, $AoLE2, $AoLE1->[0][3], $AoLE2->[0][3], 400, 500, 0.1);
									
									if($compdata->{judge} eq 'true'){
										my $tmp;
										if($h2 gt $h1){
											$tmp = $h1.",".$h2;
										}
										elsif($h1 gt $h2){
											$tmp = $h2.",".$h1;
										}
										
										$PRLG->{$tmp} = 1;
										
										if($dist_bp < 10000 || $dist_gnum < 10){
											$PRLGnear->{$tmp} = 1;
										}
									}
								}
							}
						}
					}
				}
				
				my @kPRLG = keys(%{$PRLG});
				@kPRLG = sort {$a cmp $b} @kPRLG;
				
				my @kPRLGnear = keys(%{$PRLGnear});
				@kPRLGnear = sort {$a cmp $b} @kPRLGnear;
				
				if(@kPRLG){
					my $CGR = Concatenate_graph(\@kPRLG);
					my $tophit = Select_gid_maxscore($hscore_tid_compare, $hgid2tid1, $hgid2tid2ori, $A[0], $CGR);
					$A[2] = $tophit;
					
					my $binfo = Select_transcript_maxscore($hscore_tid_compare, $hexonmatch_tid_compare, $info1, $pinfo1, $hgid2tid1, $hgid2tid2ori, $A[0], $A[2]);
					my $rev_hid = RevID($tophit);
					$rev_hid .= ".refversion";
					
					$n_genebase_summary .= join("\t", @A)."\t".join("\|", @{$CGR})."\t";
					my $num_CGR = Count_CNV_CGR($CGR);
					
					if(@kPRLGnear){		#in case of tandem reapeat
						my $CGRnear = Concatenate_graph(\@kPRLGnear);
						$CGRnear = Sort_gid_by_position($CGRnear, $hgff2ori);
						
						$n_genebase_summary .= join("\|", @{$CGRnear})."\tCNV (tandem)\ttarget < hint\t".$num_CGR."\t".$rev_hid."\t".$A[0];
						$n_genebase_summary .= "\t-";
						$n_genebase_summary .= "\t".$binfo->{end_exon_match}."\t".$binfo->{numt}."\t".$binfo->{num_hit}."\t".$binfo->{scoret};
						$n_genebase_summary .= "\t".$binfo->{pcnt_avr}."\t".$binfo->{pcov}."\t".$binfo->{Eval}."\t".$binfo->{score};
						$n_genebase_summary .= "\t".$binfo->{pep_pcnt_avr}."\t".$binfo->{pep_pcov}."\t".$binfo->{pep_Eval}."\t".$binfo->{pep_score}."\n";
					}
					else{
						$n_genebase_summary .= "-\tCNV\ttarget < hint\t".$num_CGR."\t".$rev_hid."\t".$A[0];
						$n_genebase_summary .= "\t-";
						$n_genebase_summary .= "\t".$binfo->{end_exon_match}."\t".$binfo->{numt}."\t".$binfo->{num_hit}."\t".$binfo->{scoret};
						$n_genebase_summary .= "\t".$binfo->{pcnt_avr}."\t".$binfo->{pcov}."\t".$binfo->{Eval}."\t".$binfo->{score};
						$n_genebase_summary .= "\t".$binfo->{pep_pcnt_avr}."\t".$binfo->{pep_pcov}."\t".$binfo->{pep_Eval}."\t".$binfo->{pep_score}."\n";
					}
				}
				else{
					my $homolog_str = $A[3];
					$homolog_str =~ s/,/\;/;
					
					$n_genebase_summary .= join("\t", @A)."\t-\t-\t-\t-\t-\t".$A[0]."\t-";
					$n_genebase_summary .= "\thomolog of ".$homolog_str;
					$n_genebase_summary .= "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n";
				}
			}
			else{
				$n_genebase_summary .= join("\t", @A)."\t-\t-\t-\t-\t-\t".$A[0]."\t-";
				$n_genebase_summary .= "\t-";
				$n_genebase_summary .= "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n";
			}
		}
		elsif($A[8] eq 'another gid as BDH pair'){
			unless($hcnt_hid->{$A[5]}{$A[0]}){
				print "! error : missing [hcnt_hid->{$A[5]}{$A[0]}]\n";
			}
			
			if($hcnt_hid->{$A[5]}{orth}){
				my @G = split(/,/, $hcnt_hid->{$A[5]}{nr});
				my $PRLGnear = {};
				foreach my $g1 (@G){
					foreach my $g2 (@G){
						if($g1 eq $g2){
							next;
						}
						
						my $g1chr = $hgff1->{gid2ch}{$g1};
						my $g2chr = $hgff1->{gid2ch}{$g2};
						my $g1pos0 = $hgff1->{gid2pos0}{$g1};
						my $g2pos0 = $hgff1->{gid2pos0}{$g2};
						my $g1pos1 = $hgff1->{gid2pos1}{$g1};
						my $g2pos1 = $hgff1->{gid2pos1}{$g2};
						my $g1gnum = $hgff1->{gid2gnum}{$g1};
						my $g2gnum = $hgff1->{gid2gnum}{$g2};
						
						my $dist_gnum = abs($g1gnum - $g2gnum);
						
						my $dist_bp = 10000000000000;
						if($g1chr eq $g2chr){
							if($g1pos1 < $g2pos0){
								$dist_bp = abs($g2pos0 - $g1pos1);
							}
							elsif($g1pos0 > $g2pos1){
								$dist_bp = abs($g1pos0 - $g2pos1);
							}
							elsif($g1pos0 <= $g2pos0 && $g2pos0 <= $g1pos1){
								$dist_bp = 0;
							}
							elsif($g2pos0 <= $g1pos0 && $g1pos0 <= $g2pos1){
								$dist_bp = 0;
							}
							elsif($g1pos0 <= $g2pos1 && $g2pos1 <= $g1pos1){
								$dist_bp = 0;
							}
							elsif($g2pos0 <= $g1pos1 && $g1pos1 <= $g2pos1){
								$dist_bp = 0;
							}
						}
						
						if($dist_bp < 10000 || $dist_gnum < 10){
							my $tmp;
							if($g2 gt $g1){
								$tmp = $g1.",".$g2;
							}
							elsif($g1 gt $g2){
								$tmp = $g2.",".$g1;
							}
							
							$PRLGnear->{$tmp} = 1;
						}
					}
				}
				
				my @kPRLGnear = keys(%{$PRLGnear});
				@kPRLGnear = sort {$a cmp $b} @kPRLGnear;
				
				my $tandem_repeat = "-";
				my $CNV_info = "CNV";
				if(@kPRLGnear){		#in case of tandem reapeat
					my $CGRnear = Concatenate_graph(\@kPRLGnear);
					$CGRnear = Sort_gid_by_position($CGRnear, $hgff1);
					$tandem_repeat = join("\|", @{$CGRnear});
					$CNV_info = "CNV (tandem)";
				}
				
				my $nc = $hcnt_hid->{$A[5]}{$A[0]} - 1;
				unless($Chars->[$nc]){
					$Chars->[$nc] = "other";
				}
				my $distinguisher = ".prlg-".$Chars->[$nc];
				
				my $binfo = Select_transcript_maxscore($hscore_tid_compare, $hexonmatch_tid_compare, $info1, $pinfo1, $hgid2tid1, $hgid2tid2ori, $A[0], $A[5]);
				my $rev_hid = RevID($A[5]);
				$rev_hid .= $distinguisher.".refversion";
#				$rev_hid .= $distinguisher;
				
#				$n_genebase_summary .= join("\t", @A)."\t".$hcnt_hid->{$A[5]}{nr}."\t".$tandem_repeat."\t".$CNV_info."\ttarget > hint\t".$hcnt_hid->{$A[5]}{cnt}."\t".$rev_hid."\t".$A[0];
				#disable this line because it is unable distinguish two IDs by 'REGEXP' of mysql
				
				$n_genebase_summary .= join("\t", @A)."\t".$hcnt_hid->{$A[5]}{nr}."\t".$tandem_repeat."\t".$CNV_info."\ttarget > hint\t".$hcnt_hid->{$A[5]}{cnt}."\t".$A[0]."\t-";
				$n_genebase_summary .= "\tparalog of ".$A[5];
				$n_genebase_summary .= "\t".$binfo->{end_exon_match}."\t".$binfo->{numt}."\t".$binfo->{num_hit}."\t".$binfo->{scoret};
				$n_genebase_summary .= "\t".$binfo->{pcnt_avr}."\t".$binfo->{pcov}."\t".$binfo->{Eval}."\t".$binfo->{score};
				$n_genebase_summary .= "\t".$binfo->{pep_pcnt_avr}."\t".$binfo->{pep_pcov}."\t".$binfo->{pep_Eval}."\t".$binfo->{pep_score}."\n";
			}
			else{
				$A[8] = 'another gid as BDH pair candidate';
				
				my $binfo = Select_transcript_maxscore($hscore_tid_compare, $hexonmatch_tid_compare, $info1, $pinfo1, $hgid2tid1, $hgid2tid2ori, $A[0], $A[5]);
				my $homolog_str = $A[5];
				$homolog_str =~ s/,/\;/;
				
				$n_genebase_summary .= join("\t", @A)."\t-\t-\t-\t-\t-\t".$A[0]."\t-";
				$n_genebase_summary .= "\thomolog of ".$homolog_str;
				$n_genebase_summary .= "\t".$binfo->{end_exon_match}."\t".$binfo->{numt}."\t".$binfo->{num_hit}."\t".$binfo->{scoret};
				$n_genebase_summary .= "\t".$binfo->{pcnt_avr}."\t".$binfo->{pcov}."\t".$binfo->{Eval}."\t".$binfo->{score};
				$n_genebase_summary .= "\t".$binfo->{pep_pcnt_avr}."\t".$binfo->{pep_pcov}."\t".$binfo->{pep_Eval}."\t".$binfo->{pep_score}."\n";
			}
		}
		elsif($A[8] eq 'other >2 gids as BDH pair'){
			my $homolog_str = $A[5];
			$homolog_str =~ s/,/\;/;
			
			$n_genebase_summary .= join("\t", @A)."\t-\t-\t-\t-\t-\t".$A[0]."\t-";
			$n_genebase_summary .= "\thomolog of ".$homolog_str;
			$n_genebase_summary .= "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n";
		}
		elsif($A[8] eq 'partner has >=2 gids'){
			my $homolog_str = $A[3];
			$homolog_str =~ s/,/\;/;
			
			$n_genebase_summary .= join("\t", @A)."\t-\t-\t-\t-\t-\t".$A[0]."\t-";
			$n_genebase_summary .= "\thomolog of ".$homolog_str;
			$n_genebase_summary .= "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n";
		}
		elsif($A[8] eq '1-to-1 ortholog'){
			my $binfo = Select_transcript_maxscore($hscore_tid_compare, $hexonmatch_tid_compare, $info1, $pinfo1, $hgid2tid1, $hgid2tid2ori, $A[0], $A[2]);
			my $rev_hid = RevID($A[2]);
			$rev_hid .= ".refversion";
			
			$n_genebase_summary .= join("\t", @A)."\t-\t-\t-\t-\t-\t".$rev_hid."\t".$A[0];
			$n_genebase_summary .= "\t-";
			$n_genebase_summary .= "\t".$binfo->{end_exon_match}."\t".$binfo->{numt}."\t".$binfo->{num_hit}."\t".$binfo->{scoret};
			$n_genebase_summary .= "\t".$binfo->{pcnt_avr}."\t".$binfo->{pcov}."\t".$binfo->{Eval}."\t".$binfo->{score};
			$n_genebase_summary .= "\t".$binfo->{pep_pcnt_avr}."\t".$binfo->{pep_pcov}."\t".$binfo->{pep_Eval}."\t".$binfo->{pep_score}."\n";
		}
		else{
			$n_genebase_summary .= join("\t", @A)."\t-\t-\t-\t-\t-\t".$A[0]."\t-";
			$n_genebase_summary .= "\t-";
			$n_genebase_summary .= "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n";
		}
	}
	#-----------------------------------------------------//
	
	my $genebase_posinfo_header = "gene ID\tchr\tpos_s\tpos_e\tstrand\tgnum\tstatus";
	for(my $i = 0; $i < $max_nss1; $i++){
		$genebase_posinfo_header .= "\tpartner$i\tchr$i\tpos_s$i\tpos_e$i\tstrand$i\tgnum$i\tremark$i";
	}
	$genebase_posinfo_header .= "\n";
	$genebase_posinfo = $genebase_posinfo_header.$genebase_posinfo;
	
	print "! [$cntg_total] total gene\n";
	print "! [$cntg_partner1] have some blast n or p BDH with same genomic position, Of these...\n";
	print "     [$cntg_partner11] have blastn/p BDH at least in one transcript\n";
	print "     [$cntg_partner12] have common blastn/p BDH in distinct transcripts\n";
	print "     [$cntg_partner13] have blastp BDH at least in one transcript\n";
	print "     [$cntg_partner14] have blastn BDH at least in one transcript\n";
	print "     [$cntg_partner15] (other)\n";
#	print "! [$cntg_partner1c] have partner which is shared with other gene, but primary gid can form unique pair\n";
#	print "     [$cntg_partner11c] have blastn/p BDH at least in one transcript\n";
#	print "     [$cntg_partner12c] have common blastn/p BDH in distinct transcripts\n";
#	print "     [$cntg_partner13c] have blastp BDH at least in one transcript\n";
#	print "     [$cntg_partner14c] have blastn BDH at least in one transcript\n";
#	print "     [$cntg_partner15c] (other)\n";
	print "! [$cntg_partner2c] have partner which is shared with other gene\n";
	print "     [$cntg_partner21c] have blastn/p BDH at least in one transcript\n";
	print "     [$cntg_partner22c] have common blastn/p BDH in distinct transcripts\n";
	print "     [$cntg_partner23c] have blastp BDH at least in one transcript\n";
	print "     [$cntg_partner24c] have blastn BDH at least in one transcript\n";
	print "     [$cntg_partner25c] (other)\n";
	print "! [$cntg_partner1b] have >=2 partners\n";
	print "     [$cntg_partner11b] have blastn/p BDH at least in one transcript\n";
	print "     [$cntg_partner12b] have common blastn/p BDH in distinct transcripts\n";
	print "     [$cntg_partner13b] have blastp BDH at least in one transcript\n";
	print "     [$cntg_partner14b] have blastn BDH at least in one transcript\n";
	print "     [$cntg_partner15b] (other)\n";
	print "! [$cntg_partner3] have identical BDH in blastn/p which is the partner of another gene\n";
	print "! [$cntg_partner4] have identical BDH in blastn/p\n";
	print "! [$cntg_nopartner] have no partner\n";
	print "\n";
	
#	my $rfile4 = "list_partner_ID_genebase.tsv";
	my $rfile4 = "list_partner_ID_genebase"."_c".$th_condition.".tsv";
	my $rfile4z = "list_partner_ID_genebase"."_c".$th_condition.".zip";
	SAVE($rfile4, $n_genebase_summary);
#	ZIP($rfile4z, $rfile4);
	
#	my $rfile5 = "summary_posinfo_genebase.tsv";
	my $rfile5 = "summary_posinfo_genebase"."_c".$th_condition.".tsv";
	SAVE($rfile5, $genebase_posinfo);
	
	my $file_unmapped_qexon = "unmapped_".$prefix1."_transcript"."_c".$th_condition.".fasta";
	my $file_unmapped_qprot = "unmapped_".$prefix1."_protein"."_c".$th_condition.".fasta";
	SAVE($file_unmapped_qexon, $unmapped_qexon);
	SAVE($file_unmapped_qprot, $unmapped_qprot);
	
	my $list_not_consensus_gff2 = Search_hintgene_not_consensus($hgene2ori, $n_genebase_summary);
	my $rfile6 = "list_NOTpartner_hintGff_c".$th_condition.".tsv";
	SAVE($rfile6, $list_not_consensus_gff2);
	
	#------------------------------------------------------//
	my @NR2gid = keys(%{$nr_2gid});
	@NR2gid = sort {$a cmp $b} @NR2gid;
	
	foreach my $pid (@NR2gid){
		my @G = split(/,/, $rSPT->{$pid}{partner});
		my $numG = @G;
		my $AoT = [];
		foreach my $gid (@G){
			if($gid && $hSPT->{$gid}{partner} =~ /$pid/){
				my @tmp;
				push(@tmp, $gid);
				push(@tmp, $hSPT->{$gid}{$pid}{score});
				push(@{$AoT}, \@tmp);
			}
			else{
#				print "$gid $hSPT->{$gid}{partner} $hSPT->{$gid}{info} $pid xx\n";
			}
		}
		@{$AoT} = sort {$b->[1] <=> $a->[1]} @{$AoT};
		
#		print "! [$pid] > [$AoT->[0][0]] [$AoT->[0][1]] | [$AoT->[1][0]] [$AoT->[1][1]] | [$numG]\n";
#		sleep(1);
	}
	
	#------------------------------------------------------//
	my $hcumpos1 = Fasta2cumpos($genome1);
	my $hcumpos2 = Fasta2cumpos($genome2);
	print "\n";
	
	my $genebase_cumposinfo = "gene ID\tchr\tcumpos_s\tcumpos_e\tstrand\tgnum\tstatus\tpartner0\tchr0\tpos_s0\tpos_e0\tstrand0\tgnum0\tremark0\n";
	my @GL = split(/\n/, $genebase_posinfo);
	my $hcnt_chr = {};
	my $nrseqID1 = {};
	my $nrseqID2 = {};
	my $nrchrID2 = {};
	foreach my $line (@GL){
		if($line =~ /gene/ && $line =~ /partner/){
			next;
		}
		
		my @A = split(/\t/, $line);
		unless($A[1] =~ /mitochon/ || $A[1] =~ /chloroplast/ || $A[8] =~ /mitochon/ || $A[8] =~ /chloroplast/){
			unless($A[8] eq '-'){
				my $cumpos1s = $A[2] + $hcumpos1->{$A[1]}{cumpos};
				my $cumpos1e = $A[3] + $hcumpos1->{$A[1]}{cumpos};
				my $cumpos2s = $A[9] + $hcumpos2->{$A[8]}{cumpos};
				my $cumpos2e = $A[10] + $hcumpos2->{$A[8]}{cumpos};
				
				$genebase_cumposinfo .= $A[0]."\t".$A[1]."\t".$cumpos1s."\t".$cumpos1e."\t".$A[4]."\t".$A[5]."\t".$A[6]."\t".$A[7]."\t".$A[8]."\t".$cumpos2s."\t".$cumpos2e."\t".$A[11]."\t".$A[12]."\t".$A[13]."\n";
			}
			
			if($A[8] =~ /chr00/){
				$A[8] = "unanchored_scaffolds";
			}
			elsif($A[8] ne '-'){
				$nrchrID2->{$A[8]} = 1;
			}
			
			$nrseqID1->{$A[1]} += 1;
			$nrseqID2->{$A[8]} += 1;
			$hcnt_chr->{nr}{$A[1]}{$A[8]} += 1;
			
			if($A[8] ne '-'){
				if($A[1] eq $A[8]){
					$hcnt_chr->{same}{$A[1]} += 1;
				}
				else{
					$hcnt_chr->{diff}{$A[1]} += 1;
				}
			}
		}
	}
	
	my $rfile7 = "summary_cumposinfo_genebase"."_c".$th_condition.".tsv";
	SAVE($rfile7, $genebase_cumposinfo);
	
	my @Nrc1 = keys(%{$nrseqID1});
	my @Nrc2 = keys(%{$nrseqID2});
	my @Nrch2 = keys(%{$nrchrID2});
	@Nrc1 = sort {$a cmp $b} @Nrc1;
	@Nrc2 = sort {$a cmp $b} @Nrc2;
	@Nrch2 = sort {$a cmp $b} @Nrch2;
	my $num_chr2 = @Nrch2;
	
	my $stats_tr = "-\tall\tsame\tdifferent\t".join("\t", @Nrc2)."\n";
	foreach my $ch1 (@Nrc1){
		unless($nrseqID1->{$ch1}){
			$nrseqID1->{$ch1} = 0;
		}
		unless($hcnt_chr->{same}{$ch1}){
			$hcnt_chr->{same}{$ch1} = 0;
		}
		unless($hcnt_chr->{diff}{$ch1}){
			$hcnt_chr->{diff}{$ch1} = 0;
		}
		
		$stats_tr .= $ch1."\t".$nrseqID1->{$ch1}."\t".$hcnt_chr->{same}{$ch1}."\t".$hcnt_chr->{diff}{$ch1};
		foreach my $ch2 (@Nrc2){
			unless($hcnt_chr->{nr}{$ch1}{$ch2}){
				$hcnt_chr->{nr}{$ch1}{$ch2} = 0;
			}
			$stats_tr .= "\t".$hcnt_chr->{nr}{$ch1}{$ch2};
		}
		$stats_tr .= "\n";
	}
	
	my $rfile8 = "stats_genecount"."_c".$th_condition.".tsv";
	SAVE($rfile8, $stats_tr);
	
	#---------------------------------------------------------------------------//
}

END:{
	my $end1 = "";
}


################################################################################
#-------------------------------------------------------------------------------
sub Keep_confident_hintGff{
my ($hint_gff, $gresult, $t2g, $pcnt_th) = @_;

print "! searching for confident BLAT records based on BLAST alignment...\n";
print "! reading [$hint_gff]...\n";
open(my $fh, "<", $hint_gff) or die;
my $hash = {};
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
#	Enrei331_chr01	b2h	tts	65083	65083	0	-	.	grp=Glyma.01G000100.1;pri=4;src=E
#	Enrei331_chr01	b2h	exon	65083	65553	0	-	.	grp=Glyma.01G000100.1;pri=4;src=E
#	Enrei331_chr01	b2h	exon	65868	66149	0	-	.	grp=Glyma.01G000100.1;pri=4;src=E
#	Enrei331_chr01	b2h	exon	65655	65720	0	-	.	grp=Glyma.01G000100.1;pri=4;src=E
	
	unless($line){
		next;
	}
	
	my @A = split(/\t/, $line);
	if($A[8] && $A[8] =~ /grp\=/){
		my @tag = split(/\;/, $A[8]);
		my $tid = $tag[0];
		$tid =~ s/grp\=//;
		
		if($t2g->{$tid}){
			my $gid = $t2g->{$tid};
			
			if(! $hash->{$gid}{max}){
				$hash->{$gid}{max} = $A[4];
			}
			elsif($hash->{$gid}{max} < $A[4]){
				$hash->{$gid}{max} = $A[4];
			}
			
			if(! $hash->{$gid}{min}){
				$hash->{$gid}{min} = $A[3];
			}
			elsif($hash->{$gid}{min} > $A[3]){
				$hash->{$gid}{min} = $A[3];
			}
			
			$hash->{$gid}{sid} = $A[0];
			$hash->{$gid}{data} .= $line."\n";
			$cnt++;
		}
	}
}
close $fh;

print "! [$cnt] lines\n";
$cnt = 0;

print "! reading [$gresult]...\n";
open(my $bfh, "<", $gresult) or die;
my $bhash = {};
my $rblast = "";
while(my $line = <$bfh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
#	Glyma.02G192700 Enrei331_chr02  99.123  10947   47      26      22263   33187   38383119        38372200        0.0     19638
#	Glyma.02G192700 Enrei331_chr02  98.621  5945    43      20      16510   22438   38388809        38382888        0.0     10488
#	Glyma.02G192700 Enrei331_chr02  99.134  5657    26      8       33196   38838   38372221        38366574        0.0     10154
#	Glyma.02G192700 Enrei331_chr02  98.923  5016    16      9       1       4978    38405146        38400131        0.0     8929
	
	unless($line){
		next;
	}
	
	my @A = split(/\t/, $line);
	if($A[2] > $pcnt_th && $A[11]){
		$bhash->{$A[0]}{$A[1]}{score} += $A[11];
		$bhash->{$A[0]}{$A[1]}{Pos} .= $A[8]."\n";
		$bhash->{$A[0]}{$A[1]}{Pos} .= $A[9]."\n";
		$rblast .= $line."\n";
		$cnt++;
	}
}
close $bfh;

print "! [$cnt] lines\n";
$cnt = 0;

if($rblast){
	system("rm $gresult");
	open(my $rfh, ">", $gresult);
	print $rfh $rblast;
	close $rfh;
}

print "! analyzing entries...\n";
my @GID = keys(%{$bhash});
my $hh = {};
foreach my $gid (@GID){
	my $htmp = $bhash->{$gid};
	my @SID = keys(%{$htmp});
	my $AoA = [];
	foreach my $sid (@SID){
		my @tmp;
		push(@tmp, $sid);
		push(@tmp, $htmp->{$sid}{score});
		push(@{$AoA}, \@tmp);
	}
	@{$AoA} = sort {$b->[1] <=> $a->[1]} @{$AoA};
	
	my $tophit = $AoA->[0][0];
	my @Pos = split(/\n/, $htmp->{$tophit}{Pos});
	@Pos = sort {$a <=> $b} @Pos;
	my $np = @Pos;
	
	$hh->{$gid}{sid} = $tophit;
	$hh->{$gid}{pos0} = $Pos[0];
	$hh->{$gid}{pos1} = $Pos[$np - 1];
}

@GID = keys(%{$hash});
@GID = sort {$a cmp $b} @GID;
my $revgff = "";
my $rcnt = 0;
my $dcnt = 0;
my $gid_discard = "";
foreach my $gid (@GID){
	if($hh->{$gid} && $hh->{$gid}{sid} eq $hash->{$gid}{sid}){
		if($hh->{$gid}{pos0} < $hash->{$gid}{min} && $hash->{$gid}{min} < $hh->{$gid}{pos1}){
			$revgff .= $hash->{$gid}{data};
			$rcnt++;
		}
		elsif($hh->{$gid}{pos0} < $hash->{$gid}{max} && $hash->{$gid}{max} < $hh->{$gid}{pos1}){
			$revgff .= $hash->{$gid}{data};
			$rcnt++;
		}
		else{
			$gid_discard .= $gid."\tunmatched position\t$hash->{$gid}{min} - $hash->{$gid}{max}\t$hh->{$gid}{pos0} - $hh->{$gid}{pos1}\n";
			$dcnt++;
		}
	}
	else{
		$gid_discard .= $gid."\tunmatched seqID\n";
		$dcnt++;
	}
}

if($revgff){
	print "! keep [$rcnt], discard [$dcnt] records\n";
	system("rm $hint_gff");
	open(my $rfh, ">", $hint_gff);
	print $rfh $revgff;
	close $rfh;
	print "! output [$hint_gff]\n";
}
if($gid_discard){
	open(my $rfh, ">discarded_gid.tsv");
	print $rfh $gid_discard;
	close $rfh;
}

}


#-------------------------------------------------------------------------------
sub Gff3togenomicDNA{
my ($genome, $gff, $ori_gfasta, $d) = @_;

print "! preparing genomic seq from [$genome] + [$gff]...\n";
open(my $fh, "<", $genome) or die;
my $hash = {};
my $ID;
my $seq;
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /\>/){
		if($seq){
			$seq =~ tr/acgtn/ACGTN/;
			$hash->{$ID} = $seq;
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
	$seq =~ tr/acgtn/ACGTN/;
	$hash->{$ID} = $seq;
}
#print "[$cnt] chrs or scaffolds\n";

unless(defined $d){
	$d = 10 * 1000;
}

open(my $gfh, "<", $gff) or die;
my $gfasta;
while(my $line = <$gfh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	my @A = split(/\t/, $line);
	if($A[2] && $A[2] eq 'gene'){
		if(defined $hash->{$A[0]}){
			my @T = split(/\;/, $A[8]);
			my $gid = $T[0];
			$gid =~ s/ID\=//;
			
			my $seqlen = length($hash->{$A[0]});
			my $pos0 = $A[3] - 1 - $d;
			my $pos1 = $A[4] + $d;
			
			my $seqg = substr($hash->{$A[0]}, $A[3] - 1, $A[4] - $A[3] + 1);
			my $seq5 = "";
			my $seq3 = "";
			
			if($pos0 >= 0){
				$seq5 = substr($hash->{$A[0]}, $pos0, $d);
			}
			else{
				$seq5 = substr($hash->{$A[0]}, 0, $A[3] - 1);
			}
			if($pos1 <= $seqlen){
				$seq3 = substr($hash->{$A[0]}, $A[4], $d);
			}
			else{
				$seq3 = substr($hash->{$A[0]}, $A[4]);
			}
			
			$seqg =~ tr/acgtn/ACGTN/;
			$seq5 =~ tr/ACGTN/acgtn/;
			$seq3 =~ tr/ACGTN/acgtn/;
			my $seq = $seq5.$seqg.$seq3;
			
			if($A[6] eq '-'){
				$seq =~ tr/ACGT/TGCA/;
				$seq =~ tr/acgt/tgca/;
				$seq = reverse($seq);
			}
			
			$gfasta .= ">".$gid."\n".$seq."\n";
		}
		else{
			print "! missing seq info for [$A[0]] ...\n";
		}
	}
}
close $gfh;

if($gfasta){
	open(my $rfh, ">", $ori_gfasta);
	print $rfh $gfasta;
	close $rfh;
	print "! output [$ori_gfasta]\n";
}

}


#-------------------------------------------------------------------------------
sub Search_hintgene_not_consensus{
my $hgene2ori = shift;
my $n_genebase_summary = shift;

my @GID = keys(%{$hgene2ori});
@GID = sort {$a cmp $b} @GID;

my @L = split(/\n/, $n_genebase_summary);
my $hpartner = {};
foreach my $line (@L){
	# gene ID	protein coding	candidate partner	...
	# GlymaJER.01G000001.1	TRUE	-	-	-	...
	# GlymaJER.01G000002.1	TRUE	Glyma.01G000285	Glyma.01G000285	...
	
	my @A = split(/\t/, $line);
	$hpartner->{$A[2]} = 1;
}

my $list = "";
foreach my $gid (@GID){
	unless($hpartner->{$gid}){
		$list .= $gid."\n";
	}
}

return $list;
}


#-------------------------------------------------------------------------------
sub Check_protein_coding{
my $hseq_qexon = shift;
my $hseq_qprot = shift;

my @ID = keys(%{$hseq_qexon});
@ID = sort {$a cmp $b} @ID;

my $hash = {};
foreach my $id (@ID){
	if($hseq_qprot->{$id}){
		$hash->{$id} = "protein coding";
	}
	else{
		$hash->{$id} = "non-coding";
	}
}

return $hash;
}


#-------------------------------------------------------------------------------
sub Count_rSPT_partner{
my $str = shift;

my @tmp = split(/,/, $str);
my @G;
foreach my $g (@tmp){
	if($g){
		push(@G, $g);
	}
}
my $n = @G;

return $n;
}


#-------------------------------------------------------------------------------
sub Count_CNV_CGR{
my $CGR = shift;

my @Num;
foreach my $cgr (@{$CGR}){
	my @G = split(/,/, $cgr);
	my $n = @G;
	push(@Num, $n);
}
my $str = join("\|", @Num);

return $str;
}


#-------------------------------------------------------------------------------
sub Select_transcript_maxscore{
my $hscore_tid_compare = shift;
my $hexonmatch_tid_compare = shift;
my $info1 = shift;
my $pinfo1 = shift;
my $g2t1 = shift;
my $g2t2 = shift;
my $gid1 = shift;
my $gid2 = shift;

my @TID1 = split(/\n/, $g2t1->{$gid1});
my @TID2 = split(/\n/, $g2t2->{$gid2});

my $AoA = [];
foreach my $tid1 (@TID1){
	foreach my $tid2 (@TID2){
		if($hscore_tid_compare->{$tid1}{$tid2}){
			my @A;
			push(@A, $tid1);
			push(@A, $tid2);
			push(@A, $hscore_tid_compare->{$tid1}{$tid2});
			push(@{$AoA}, \@A);
		}
	}
}
@{$AoA} = sort {$b->[2] <=> $a->[2]} @{$AoA};

my $top1 = $AoA->[0][0];
my $top2 = $AoA->[0][1];

unless(defined $top1 || defined $top2){
	print "! missing topID info for [$gid1] [$gid2]\n";
	sleep(1);
}

my $r = {};
if($info1->{$top1}{hitinfo}{$top2}){
	$r->{pcnt} = $info1->{$top1}{hitinfo}{$top2}{pcnt};
	$r->{pcov} = $info1->{$top1}{hitinfo}{$top2}{pcov};
	$r->{Eval} = $info1->{$top1}{hitinfo}{$top2}{Eval};
	$r->{score} = $info1->{$top1}{hitinfo}{$top2}{score};
	$r->{pcnt_avr} = $info1->{$top1}{hitinfo}{$top2}{pcnt_avr};
}
else{
	$r->{pcnt} = "-";
	$r->{pcov} = "-";
	$r->{Eval} = "-";
	$r->{score} = "-";
	$r->{pcnt_avr} = "-";
}

if($pinfo1->{$top1}{hitinfo}{$top2}){
	$r->{pep_pcnt} = $pinfo1->{$top1}{hitinfo}{$top2}{pcnt};
	$r->{pep_pcov} = $pinfo1->{$top1}{hitinfo}{$top2}{pcov};
	$r->{pep_Eval} = $pinfo1->{$top1}{hitinfo}{$top2}{Eval};
	$r->{pep_score} = $pinfo1->{$top1}{hitinfo}{$top2}{score};
	$r->{pep_pcnt_avr} = $pinfo1->{$top1}{hitinfo}{$top2}{pcnt_avr};
}
else{
	$r->{pep_pcnt} = "-";
	$r->{pep_pcov} = "-";
	$r->{pep_Eval} = "-";
	$r->{pep_score} = "-";
	$r->{pep_pcnt_avr} = "-";
}

if($hexonmatch_tid_compare->{$top1}{$top2}){
	$r->{end_exon_match} = $hexonmatch_tid_compare->{$top1}{$top2}{end_exon_match};
	$r->{scoret} = $hexonmatch_tid_compare->{$top1}{$top2}{scoret};
	$r->{scoreh} = $hexonmatch_tid_compare->{$top1}{$top2}{scoreh};
	$r->{num_hit} = $hexonmatch_tid_compare->{$top1}{$top2}{num_hit};
	$r->{numt} = $hexonmatch_tid_compare->{$top1}{$top2}{numt};
	$r->{numh} = $hexonmatch_tid_compare->{$top1}{$top2}{numh};
}
else{
	$r->{end_exon_match} = "-";
	$r->{scoret} = "-";
	$r->{scoreh} = "-";
	$r->{num_hit} = "-";
	$r->{numt} = "-";
	$r->{numh} = "-";
}

return $r;
}


#-------------------------------------------------------------------------------
sub RevID{
my $id = shift;

#if($id =~ /\./){
#	my @tmp = split(/\./, $id);
#	$id = $tmp[0];
#}

return $id;
}


#-------------------------------------------------------------------------------
sub Assign_and_Check_redundancy{
my $GSL = shift;

my $hash = {};
foreach my $l (@{$GSL}){
	my @A = split(/\t/, $l);
	if($A[8] eq '1-to-1 ortholog' && $A[2] ne '-'){
		$hash->{$A[2]}{orth} = $A[0];
		$hash->{$A[2]}{nr} = $A[0];
		$hash->{$A[2]}{cnt} += 1;
	}
}

foreach my $l (@{$GSL}){
	my @A = split(/\t/, $l);
	if($A[8] eq '-'){
		next;
	}
	elsif($A[8] eq 'gid has >=2 partners'){					#one query gene has >=2 hits in the genomic position
		unless($A[3] eq '-'){
			my @H = split(/,/, $A[3]);
			foreach my $h1 (@H){
				if($hash->{$h1}{orth}){
					print "! [$h1] has been selected as ortholog of [$hash->{$h1}{orth}]\n";
				}
				$hash->{$h1}{cnt} += 1;
				$hash->{$h1}{$A[0]} = $hash->{$h1}{cnt};
			}
		}
	}
	elsif($A[8] eq 'another gid as BDH pair'){
#			if($hash->{$A[5]}{orth}){
#				print "! [$A[5]] has been selected as ortholog of [$hash->{$A[5]}{orth}]\n";
#			}
		
		$hash->{$A[5]}{nr} .= ",".$A[0];
		$hash->{$A[5]}{cnt} += 1;
		$hash->{$A[5]}{$A[0]} = $hash->{$A[5]}{cnt};
	}
	elsif($A[8] eq 'partner has >=2 gids'){
		if($hash->{$A[3]}{orth}){
			print "! [$A[3]] has been selected as ortholog of [$hash->{$A[5]}{orth}]\n";
		}
		$hash->{$A[3]}{cnt} += 1;
		$hash->{$A[3]}{$A[0]} = $hash->{$A[3]}{cnt};
	}
}

return $hash;
}


#-------------------------------------------------------------------------------
sub Sort_gid_by_position{
my $A = shift;
my $hgff2ori = shift;

my @R;
foreach my $str (@{$A}){
	my @G = split(/,/, $str);
	
	my $AoA = [];
	foreach my $gid (@G){
		my $gnum = $hgff2ori->{gid2gnum}{$gid};
		my @tmp;
		push(@tmp, $gid);
		push(@tmp, $gnum);
		push(@{$AoA}, \@tmp);
	}
	@{$AoA} = sort {$a->[1] <=> $b->[1]} @{$AoA};
	
	my @Sorted;
	foreach my $tmp (@{$AoA}){
		push(@Sorted, $tmp->[0]);
	}
	my $nstr = join(",", @Sorted);
	push(@R, $nstr);
}

return \@R;
}


#-------------------------------------------------------------------------------
sub Concatenate_graph{
my $A = shift;

my $group = {};
foreach my $str (@{$A}){
	my @B = split(/,/, $str);
	$group->{$B[0]}{$B[1]} = 1;
	$group->{$B[1]}{$B[0]} = 1;
}

my @nrID = keys(%{$group});
@nrID = sort {$a cmp $b} @nrID;

my $selected = {};
my $G = {};
my $g = 1;
foreach my $id (@nrID){
	if($selected->{$id}){
		next;
	}
	
	if($group->{$id}){
		my @QID;
		push(@QID, $id);
		my $tested = {};
		
		for(my $k = 0; $k < 100; $k++){
			my $nums = 0;
			foreach my $qid (@QID){
				if($tested->{$qid}){
					$nums++;
					next;
				}
				
				$selected->{$qid} = $g;
				$G->{$g}{$qid} = 1;
				
				my $tmp1 = $group->{$qid};
				my @T2 = keys(%{$tmp1});
				@T2 = sort {$a cmp $b} @T2;
				
				foreach my $tid2 (@T2){
					$selected->{$tid2} = $g;
					$G->{$g}{$tid2} = 1;
					
					unless($tested->{$tid2}){
						push(@QID, $tid2);
					}
				}
				
				$tested->{$qid} = 1;
			}
			
			my $numQ = @QID;
			if($numQ == $nums){
				last;
			}
		}
		$g++;
	}
}

my @nrG = keys(%{$G});
@nrG = sort {$a <=> $b} @nrG;
my $numG = @nrG;

my @R;
for(my $i = 1; $i <= $numG; $i++){
	my $tmp = $G->{$i};
	my @snrID = keys(%{$tmp});
	my $snrID_str = join(",", @snrID);
	push(@R, $snrID_str);
}

return \@R;
}


#-------------------------------------------------------------------------------
sub Select_gid_maxscore{
my $hscore_tid_compare = shift;
my $g2t1 = shift;
my $g2t2 = shift;
my $gid1 = shift;
my $GID2 = shift;

my @TID1 = split(/\n/, $g2t1->{$gid1});
my $AoA = [];
foreach my $str2 (@{$GID2}){
	my @subGID2 = split(/,/, $str2);
	foreach my $subgid2 (@subGID2){
		my @TID2 = split(/\n/, $g2t2->{$subgid2});
		
		foreach my $tid1 (@TID1){
			foreach my $tid2 (@TID2){
				if($hscore_tid_compare->{$tid1}{$tid2}){
					my @A;
					push(@A, $subgid2);
					push(@A, $tid2);
					push(@A, $hscore_tid_compare->{$tid1}{$tid2});
					push(@{$AoA}, \@A);
				}
			}
		}
	}
}
@{$AoA} = sort {$b->[2] <=> $a->[2]} @{$AoA};

return $AoA->[0][0];
}


#-------------------------------------------------------------------------------
sub Revise_genebase_summary{
my $lines = shift;
my $Char2 = shift;

my @L = split(/\n/, $lines);
my $horth = {};
my $hanot = {};
my $hp2g = {};
my $hg2p = {};
my $hNA = {};
my $cnt_orth = 0;
my $cnt_anot = 0;
my $cnt_p2g = 0;
my $cnt_g2p = 0;
my $cnt_NA = 0;

#collecting data as hash
foreach my $l (@L){
	my @A = split(/\t/, $l);
	if($A[1] ne '-'){										#1-to-1 ortholog partner
		$horth->{Fwd}{$A[0]} = $A[1];
		$horth->{Rev}{$A[1]} = $A[0];
		$cnt_orth++;
	}
	elsif($A[7] eq 'BDH forms pair with another gid'){		#another query paralog is present for hit gene as BDH, this may be CNV candidate
		unless($A[4] eq '-'){
			$hanot->{Fwd}{$A[0]} = $A[4];
			$hanot->{Rev}{$A[4]} .= $A[0].",";
			$cnt_anot++;
		}
	}
	elsif($A[7] eq 'partner has >=2 gids'){					#one hit gene has >=2 query genes in the genomic position
		unless($A[2] eq '-'){
			$hp2g->{Fwd}{$A[0]} = $A[2];
			$hp2g->{Rev}{$A[2]} .= $A[0].",";
			$cnt_p2g++;
		}
	}
	elsif($A[7] eq 'gid has >=2 partners'){					#one query gene has >=2 hits in the genomic position
		unless($A[2] eq '-'){
			$hg2p->{Fwd}{$A[0]} = $A[2];
			
			my @H = split(/,/, $A[2]);
			foreach my $h (@H){
				$hg2p->{Rev}{$h} .= $A[0].",";
			}
			$cnt_g2p++;
		}
	}
	elsif($A[7] eq '-'){
		$hNA->{$A[0]} = 1;
		$cnt_NA++;
	}
}

#attaching distinguisher (e.g. paraA, paraB ...etc)
my $r = shift(@L)."\tparalog_distinguisher\n";
foreach my $l (@L){
	my @A = split(/\t/, $l);
	if($A[1] ne '-'){										#1-to-1 ortholog partner
		$r .= $l."\t-\n";
	}
	elsif($A[7] eq 'BDH forms pair with another gid'){		#another query paralog is present for hit gene as BDH, this may be CNV candidate
		my $dc;
		unless($A[4] eq '-'){
			my @H = split(/,/, $hanot->{Rev}{$A[4]});
			my $n = 0;
			foreach my $h (@H){
				if($A[0] eq $h){
					$dc = $Char2->[$n];
					$n++;
				}
			}
		}
		
		if($dc){
			$r .= $l."\tprlg-".$dc."\n";
		}
	}
	elsif($A[7] eq 'partner has >=2 gids'){					#one hit gene has >=2 query genes in the genomic position
		$r .= $l."\t-\n";
	}
	elsif($A[7] eq 'gid has >=2 partners'){					#one query gene has >=2 hits in the genomic position
		$r .= $l."\t-\n";
	}
	elsif($A[7] eq '-'){
		$r .= $l."\t-\n";
	}
}

my $rh = {};
$rh->{horth} = $horth;
$rh->{hanot} = $hanot;
$rh->{hp2g} = $hp2g;
$rh->{hg2p} = $hg2p;
$rh->{hNA} = $hNA;
$rh->{r} = $r;

return $rh;
}


#-------------------------------------------------------------------------------
sub Gen_AoChars{

my @Char;
for(my $i = 0; $i < 26; $i++){
	if($i == 0){push(@Char, "A");}
	elsif($i == 1){push(@Char, "B");}
	elsif($i == 2){push(@Char, "C");}
	elsif($i == 3){push(@Char, "D");}
	elsif($i == 4){push(@Char, "E");}
	elsif($i == 5){push(@Char, "F");}
	elsif($i == 6){push(@Char, "G");}
	elsif($i == 7){push(@Char, "H");}
	elsif($i == 8){push(@Char, "I");}
	elsif($i == 9){push(@Char, "J");}
	elsif($i == 10){push(@Char, "K");}
	elsif($i == 11){push(@Char, "L");}
	elsif($i == 12){push(@Char, "M");}
	elsif($i == 13){push(@Char, "N");}
	elsif($i == 14){push(@Char, "O");}
	elsif($i == 15){push(@Char, "P");}
	elsif($i == 16){push(@Char, "Q");}
	elsif($i == 17){push(@Char, "R");}
	elsif($i == 18){push(@Char, "S");}
	elsif($i == 19){push(@Char, "T");}
	elsif($i == 20){push(@Char, "U");}
	elsif($i == 21){push(@Char, "V");}
	elsif($i == 22){push(@Char, "W");}
	elsif($i == 23){push(@Char, "X");}
	elsif($i == 24){push(@Char, "Y");}
	elsif($i == 25){push(@Char, "Z");}
}

my @Char2 = @Char;
for(my $i = 0; $i < 26; $i++){
	my $str = $Char[$i];
	
	for(my $i = 0; $i < 26; $i++){
		if($i == 0){$str .= "A";}
		elsif($i == 1){$str .= "B";}
		elsif($i == 2){$str .= "C";}
		elsif($i == 3){$str .= "D";}
		elsif($i == 4){$str .= "E";}
		elsif($i == 5){$str .= "F";}
		elsif($i == 6){$str .= "G";}
		elsif($i == 7){$str .= "H";}
		elsif($i == 8){$str .= "I";}
		elsif($i == 9){$str .= "J";}
		elsif($i == 10){$str .= "K";}
		elsif($i == 11){$str .= "L";}
		elsif($i == 12){$str .= "M";}
		elsif($i == 13){$str .= "N";}
		elsif($i == 14){$str .= "O";}
		elsif($i == 15){$str .= "P";}
		elsif($i == 16){$str .= "Q";}
		elsif($i == 17){$str .= "R";}
		elsif($i == 18){$str .= "S";}
		elsif($i == 19){$str .= "T";}
		elsif($i == 20){$str .= "U";}
		elsif($i == 21){$str .= "V";}
		elsif($i == 22){$str .= "W";}
		elsif($i == 23){$str .= "X";}
		elsif($i == 24){$str .= "Y";}
		elsif($i == 25){$str .= "Z";}
	}
	
	push(@Char2, $str);
}

return \@Char2;
}


#-------------------------------------------------------------------------------
sub Make_table{
my $num_chr = shift;
my $table_pt = shift;

my @P;
for(my $i = 0; $i < $num_chr + 2; $i++){
	if(($i + 1) == $table_pt){
		$P[$i] = $table_pt;
	}
	else{
		$P[$i] = "n";
	}
}

my $str = join("\t", @P);
$str =~ s/n//g;

return $str;
}


#-------------------------------------------------------------------------------
sub CMD{
my $cmd = shift;

print "! cmd=[$cmd]\n";
if(system($cmd) != 0){
	print "! failed.\n";
	die;
}

}


#-------------------------------------------------------------------------------
sub SameChr_or_Not{
my $ID1 = shift;
my $ID2 = shift;

my $str = "same_chr";
if($ID1 ne $ID2){
	$str = "diff_chr";
}

return $str;
}


#-------------------------------------------------------------------------------
sub Fasta2cumpos{
my $file = shift;

print "\n! calculating cumlative seq length [$file]...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
my $ID;
my $seq;
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ />/){
		unless($cnt == 0){
			$hash->{$ID}{len} = length($seq);
		}
		
		$cnt++;
		
		$ID = $line;
		$ID =~ s/\>//;
		$seq = "";
		
		my @tmp = split(/\s/, $ID);
		$ID = $tmp[0];
	}
	else{
		$line =~ s/\.//;
		$seq .= $line;
	}
}
close $fh;

if($seq){
	$hash->{$ID}{len} = length($seq);
}

my @IDs0 = keys(%{$hash});
@IDs0 = sort {$a cmp $b} @IDs0;

my @IDs;
foreach my $id (@IDs0){
	unless($id =~ /mitochon/ || $id =~ /chloroplast/){
		unless($id =~ /chr00/ || $id =~ /scaffold/){
			push(@IDs, $id);
		}
	}
}
foreach my $id (@IDs0){
	unless($id =~ /mitochon/ || $id =~ /chloroplast/){
		if($id =~ /chr00/ || $id =~ /scaffold/){
			push(@IDs, $id);
		}
	}
}

my $cumpos = 0;
foreach my $id (@IDs){
	print " [$id] : [$cumpos]\n";
	$hash->{$id}{cumpos} = $cumpos;
	$cumpos += $hash->{$id}{len};
}

return $hash;
}


#-------------------------------------------------------------------------------
sub Open_fasta_as_hash{
my $file = shift;

print "! open fasta as hash [$file]...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
my $ID;
my $seq;
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ />/){
		unless($cnt == 0){
			$hash->{$ID}{seq} = ">".$ID."\n".$seq."\n";
			$hash->{$ID}{len} = length($seq);
		}
		
		$cnt++;
		
		$ID = $line;
		$ID =~ s/\>//;
		$seq = "";
		
		my @tmp = split(/\s/, $ID);
		$ID = $tmp[0];
	}
	else{
		$line =~ s/\.//;
		$seq .= $line;
	}
}
close $fh;

if($seq){
	$hash->{$ID}{seq} = ">".$ID."\n".$seq."\n";
	$hash->{$ID}{len} = length($seq);
}

print "! [$cnt] entries\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Tid2gid{
my $str = shift;

my $ver = "null";
if($str =~ /MELO3C/){
	if($str =~ /T/){
		$str =~ s/T/\./;
		$ver = "CM3";
	}
	elsif($str =~ /P/){
		$str =~ s/P/\./;
		$ver = "CM3";
	}
	elsif($str =~ /\.2/){
		$ver = "CM4";
	}
}
elsif($str =~ /MELO\.jh/){
	$ver = "dotID2";
}
elsif($str =~ /Glyma\./){
	$ver = "dotID";
}

my @tmp = split(/\./, $str);
my $num_tmp = @tmp;

my $gid;
if($ver eq 'CM4'){
	$gid = $tmp[0].".2";
}
elsif($ver eq 'CM3'){
	$gid = $tmp[0];
}
elsif($ver eq 'dotID'){
	$gid = $tmp[0].".".$tmp[1];
}
elsif($ver eq 'dotID2'){
	$gid = $tmp[0].".".$tmp[1].".".$tmp[2];
}
else{
	$gid = $tmp[0];
}

return $gid;
}


#-------------------------------------------------------------------------------
sub Judge_score{
my $hSPT = shift;
my $gid = shift;
my $ssid = shift;

my $judge_score = 0;
if($hSPT->{$gid}{$ssid}{info} =~ /also_blastnp_BDH/){
	$judge_score = 10;
}
elsif($hSPT->{$gid}{$ssid}{info} =~ /partner_is_blastn_BDH_but_not_blastp_BDH/ && $hSPT->{$gid}{$ssid}{info} =~ /partner_is_blastp_BDH_but_not_blastn_BDH/){
	$judge_score = 5;
}
elsif($hSPT->{$gid}{$ssid}{info} =~ /partner_is_blastp_BDH_but_not_blastn_BDH/){
	$judge_score = 2;
}
elsif($hSPT->{$gid}{$ssid}{info} =~ /partner_is_blastn_BDH_but_not_blastp_BDH/){
	$judge_score = 1;
}

return $judge_score;
}


#-------------------------------------------------------------------------------
sub Erase_entry{
my $rSPT = shift;
my $prev_hgid = shift;
my $gid = shift;

my @tmp = split(/,/, $rSPT->{$prev_hgid}{partner});
my @ntmp;
foreach my $hid (@tmp){
	if($hid && $hid ne $gid){
		push(@ntmp, $hid);
	}
}
$rSPT->{$prev_hgid}{partner} = join(",", @ntmp);

return $rSPT;
}


#-------------------------------------------------------------------------------
sub Tid2hash{
my $TID = shift;

my $hash = {};
foreach my $tid (@{$TID}){
	my $ver = "null";
	if($tid =~ /MELO3C/){
		if($tid =~ /T/){
			$tid =~ s/T/\./;
			$ver = "CM3";
		}
		elsif($tid =~ /P/){
			$tid =~ s/P/\./;
			$ver = "CM3";
		}
		elsif($tid =~ /\.2/){
			$ver = "CM4";
		}
	}
	elsif($tid =~ /MELO\.jh/){
		$ver = "dotID2";
	}
	elsif($tid =~ /Glyma\./){
		$ver = "dotID";
	}
	
	my @tmp = split(/\./, $tid);
	my $num_tmp = @tmp;
	
	my $gid;
	if($ver eq 'CM4'){
		$gid = $tmp[0].".2";
	}
	elsif($ver eq 'CM3'){
		$gid = $tmp[0];
	}
	elsif($ver eq 'dotID'){
		$gid = $tmp[0].".".$tmp[1];
	}
	elsif($ver eq 'dotID2'){
		$gid = $tmp[0].".".$tmp[1].".".$tmp[2];
	}
	else{
		$gid = $tmp[0];
	}
	
	$hash->{$gid} .= $tid."\n";
}

return $hash;
}


#-------------------------------------------------------------------------------
sub Search_same_hit_np{
my $str1 = shift;
my $str2 = shift;

my @A = split(/,/, $str1);
my @B = split(/,/, $str2);

my $hit = "-";
foreach my $id1 (@A){
	foreach my $id2 (@B){
		if($id1 eq $id2){
			$hit = $id1;
			last;
		}
	}
	
	if($hit ne '-'){
		last;
	}
}

return $hit;
}


#-------------------------------------------------------------------------------
sub Compare_exon_structure{
my $tid = shift;
my $grp = shift;
my $Aot = shift;
my $Aoh = shift;
my $tid_strand = shift;
my $grp_strand = shift;
my $ew = shift;
my $uw = shift;
my $rth = shift;
my $prnt = shift;

my $numt = @{$Aot};
my $numh = @{$Aoh};

if($tid_strand eq '-'){
	@{$Aot} = sort {$b->[0] <=> $a->[0]} @{$Aot};
}
if($grp_strand eq '-'){
	@{$Aoh} = sort {$b->[0] <=> $a->[0]} @{$Aoh};
}

my $r;
my $num_hit = 0;
my $point = 0;
my $mt = {};
my $mh = {};
my @Hpos;
for(my $i = 0; $i < $numt; $i++){
	my $val_t0 = $Aot->[$i][0];
	my $val_t1 = $Aot->[$i][1];
	
	if($mt->{$i}){
		next;
	}
	
	for(my $j = 0; $j < $numh; $j++){
		unless(defined $Aoh->[$j][0] || defined $Aoh->[$j][1]){
			print "! missing value for [$grp] [$j] [$numh] when analyzing [$tid]\n";
			sleep(1);
			next;
		}
		unless($tid_strand eq $grp_strand){
			next;
		}
		
		my $val_h0 = $Aoh->[$j][0];
		my $val_h1 = $Aoh->[$j][1];
		
		my $d0 = abs($val_t0 - $val_h0);
		my $d1 = abs($val_t1 - $val_h1);
		
		my $th0;
		my $th1;
		if($i == 0 || $j == 0){
			if($numt == 1){
				if($numh == 1){
					$th0 = $uw;
					$th1 = $uw;
				}
				else{
					if($j == 0){
						$th0 = $uw;
						$th1 = $ew;
					}
					elsif($j == $numh - 1){
						$th0 = $ew;
						$th1 = $uw;
					}
					else{
						$th0 = $ew;
						$th1 = $ew;
					}
				}
			}
			else{
				$th0 = $uw;
				$th1 = $ew;
			}
		}
		elsif($i == $numt - 1 || $j == $numh - 1){
			$th0 = $ew;
			$th1 = $uw;
		}
		else{
			$th0 = $ew;
			$th1 = $ew;
		}
		
#		if($d0 <= $th0 && $d1 <= $th1){
		if($d0 <= $th0 || $d1 <= $th1){
			unless($mh->{$j}){
				$r .= "match   | T0=[$val_t0] H0=[$val_h0] d0=[$d0] th=[$th0] | T1=[$val_t1] H1=[$val_h1] d1=[$d1] th=[$th1]\n";
				
				my $subpt0 = 0;
				if($d0 == 0){
					$subpt0 = 1.5;
				}
				else{
					$subpt0 = 1 / $d0;
				}
				
				my $subpt1 = 0;
				if($d1 == 0){
					$subpt1 = 1.5;
				}
				else{
					$subpt1 = 1 / $d1;
				}
				
				$point += $subpt0 + $subpt1;
				$num_hit++;
				$mt->{$i} = 1;
				$mh->{$j} = 1;
				
				push(@Hpos, $val_t0);
				push(@Hpos, $val_t1);
				
				last;
			}
		}
		else{
			unless($mh->{$j}){
				$r .= "unmatch | T0=[$val_t0] H0=[$val_h0] d0=[$d0] th=[$th0] | T1=[$val_t1] H1=[$val_h1] d1=[$d1] th=[$th1]\n";
			}
		}
	}
}

my $judge = "false";
my $scoret = ($num_hit / $numt);
my $scoreh = ($num_hit / $numh);

if($scoret >= $rth && $scoreh >= $rth){
	$judge = "true";
}
elsif($scoreh >= 0.95){
	$judge = "true";
}

my $hpos_max = 0;
my $hpos_min = 0;
foreach my $hpos (@Hpos){
	unless($hpos_max){
		$hpos_max = $hpos;
	}
	elsif($hpos_max < $hpos){
		$hpos_max = $hpos;
	}
	
	unless($hpos_min){
		$hpos_min = $hpos;
	}
	elsif($hpos_min > $hpos){
		$hpos_min = $hpos;
	}
}

$r .= "! num_exon : T=[$numt] H=[$numh] match=[$num_hit] | hpos0=[$hpos_min] hpos1=[$hpos_max] | judge=[$judge] | [$tid] vs [$grp]\n";

#if($tid =~ /MELO.efh000269/){
#	print "$r";
#}

my $rhash = {};
$rhash->{judge} = $judge;
$rhash->{score} = sprintf("%.2f", $scoret + $scoreh);
$rhash->{point} = $point;
$rhash->{scoret} = sprintf("%.2f", $scoret);
$rhash->{scoreh} = sprintf("%.2f", $scoreh);
$rhash->{rth} = $rth;
$rhash->{str} = $tid."\t".$grp."\t".$numt."\t".$numh."\t".$num_hit."\t".$judge."\n";
$rhash->{min} = $hpos_min;
$rhash->{max} = $hpos_max;
$rhash->{num_hit} = $num_hit;
$rhash->{numt} = $numt;
$rhash->{numh} = $numh;

my $f = 0;
my $e = $numt - 1;
if($mt->{$f} && $mt->{$e}){
	$rhash->{end_exon_match} = 11;
}
elsif($mt->{$f}){
	$rhash->{end_exon_match} = 10;
}
elsif($mt->{$e}){
	$rhash->{end_exon_match} = 1;
}
else{
	$rhash->{end_exon_match} = 0;
}

return $rhash;
}


#----------------------------------------------------------
sub Fasta_nameclean{
my $query = shift;

#print "! parsing query ...\n";
open(my $fh, "<", $query) or die;
my $nfasta = "";
my $IDinfo = "";
my $ID;
my $ID0;
my $seq;
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ />/){
		unless($cnt == 0){
			$seq =~ tr/acgtn/ACGTN/;
			$nfasta .= ">".$ID."\n".$seq."\n";
			$IDinfo .= $ID.",".$ID0."\n";
		}
		
		$ID0 = $line;
		$ID0 =~ s/\>//;
		$seq = "";
		
		$ID = $ID0;
		my @tmp = split(/\s/, $ID);
		$tmp[0] =~ s/grp\=//;
		$ID = $tmp[0];
		
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
	$seq =~ tr/acgtn/ACGTN/;
	$nfasta .= ">".$ID."\n".$seq."\n";
	$IDinfo .= $ID.",".$ID0."\n";
}

#print "! [$cnt0] entries\n";

if($nfasta){
	system("rm $query");
	
	open(my $rfh, ">", $query) or die;
	print $rfh $nfasta;
	close $rfh;
}
#print "! output [$rfile0]\n";

}


#-------------------------------------------------------------------------------
sub Integrate_gffs{
my $list_integrate = shift;
my $sortedGID = shift;
my $g2t1 = shift;
my $g2t2 = shift;
my $htidall1 = shift;
my $htidall2 = shift;
my $hgene1 = shift;
my $hgene2 = shift;
my $igff1 = shift;
my $using_ESThint_mode = shift;

my @L = split(/\n/, $list_integrate);
my $hcand = {};
foreach my $line (@L){
	my @A = split(/\t/, $line);
	my @C = split(/,/, $A[0]);
	foreach my $cgid (@C){
		$hcand->{gid}{$cgid} = $line;
	}
}

my $ngff = "##gff-version 3\n";
my $selected = {};
foreach my $gid (@{$sortedGID}){
	if($selected->{$gid}){
		next;
	}
	
	if($hcand->{gid}{$gid}){
		my @A = split(/\t/, $hcand->{gid}{$gid});
		my @C = split(/,/, $A[0]);
		my $hintgid = $A[1];
		my $altname = $A[3];
		
		my $i = 1;
		my $subgff;
		my @eachTID2 = split(/\n/, $g2t2->{$hintgid});
		my @Pos;
		
		if($using_ESThint_mode eq 'false'){
			foreach my $tid (@eachTID2){
				my $altname_tid = $altname.".t".$i;
				$i++;
				
				$htidall2->{$tid} =~ s/$tid/$altname_tid/g;
				$htidall2->{$tid} =~ s/$hintgid/$altname/g;
				$subgff .= $htidall2->{$tid};
				
				my @E = split(/\n/, $htidall2->{$tid});
				foreach my $e (@E){
					my @subA = split(/\t/, $e);
					push(@Pos, $subA[3]);
					push(@Pos, $subA[4]);
				}
			}
		}
		
		foreach my $cgid (@C){
			my @eachTID1 = split(/\n/, $g2t1->{$cgid});
			foreach my $tid (@eachTID1){
				my $altname_tid = $altname.".t".$i;
				$i++;
				
				unless(defined $htidall1->{$tid}){
					print "! missing data for [$cgid] [$tid]\n";
					die;
				}
				
				$htidall1->{$tid} =~ s/$tid/$altname_tid/g;
				$htidall1->{$tid} =~ s/$cgid/$altname/g;
				$subgff .= $htidall1->{$tid};
				
				my @E = split(/\n/, $htidall1->{$tid});
				foreach my $e (@E){
					my @subA = split(/\t/, $e);
					push(@Pos, $subA[3]);
					push(@Pos, $subA[4]);
				}
			}
			$selected->{$cgid} = 1;
		}
		
		@Pos = sort {$a <=> $b} @Pos;
		my $npos = @Pos;
		my $minp = $Pos[0];
		my $maxp = $Pos[$npos - 1];
		
		$hgene1->{$altname} =~ s/\n//;
		my @G = split(/\t/, $hgene1->{$altname});
		$G[3] = $minp;
		$G[4] = $maxp;
		
		$ngff .= join("\t", @G)."\n".$subgff;
	}
	else{
		$ngff .= $hgene1->{$gid};
		
		my @eachTID1 = split(/\n/, $g2t1->{$gid});
		foreach my $tid (@eachTID1){
			$ngff .= $htidall1->{$tid};
		}
		
		$selected->{$gid} = 1;
	}
}

if($ngff){
	open(my $rfh, ">", $igff1);
	print $rfh $ngff;
	close $rfh;
	print "! output [$igff1]\n";
}

}


#-------------------------------------------------------------------------------
sub Search_partner{
my $file = shift;
my $gff1 = shift;
my $gff2 = shift;
my $hexon1 = shift;
my $hexon2 = shift;

open(my $fh, "<", $file);
my $hhit = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	my @A = split(/,/, $line);
	my $id1 = $A[0];
	my $hit1 = $A[1];
	
	$hhit->{id1}{cnt} += 1;
	$hhit->{id1}{hit} .= $hit1.",";
}

my @TID = keys(%{$hhit});
@TID = sort {$a cmp $b} @TID;

my $r;
my $cntm = 0;
my $sum = 0;
foreach my $hit1 (@TID){
#	if($hhit->{$hit1}{cnt} >= 1 && $hhit->{$hit1}{cnt} <= 4){
	if($hhit->{$hit1}{cnt} >= 2 && $hhit->{$hit1}{cnt} <= 4){		#changed '1' to '2'
		my @ID1 = split(/,/, $hhit->{$hit1}{id1});
		@ID1 = sort {$a cmp $b} @ID1;
		
		my $subexon2 = $hexon2->{$hit1};
		my @subTID2 = keys(%{$subexon2});
		@subTID2 = sort {$a cmp $b} @subTID2;
		
		my $hjudge = {};
		foreach my $tid2 (@subTID2){
			my $Aoexon2 = $hexon2->{$hit1}{$tid2};
			
			foreach my $id1 (@ID1){					#each gene ID1
				my $subexon1 = $hexon1->{$id1};
				my @subTID1 = keys(%{$subexon1});
				@subTID1 = sort {$a cmp $b} @subTID1;
				
				foreach my $tid1 (@subTID1){		#each transcript ID for each ID1
					my $Aoexon1 = $hexon1->{$id1}{$tid1};
					
					if($Aoexon1->[0][3] eq $Aoexon2->[0][3] && $Aoexon1->[0][4] eq $Aoexon2->[0][4]){		#both strand and chromosome are identical
						foreach my $exon1 (@{$Aoexon1}){
							foreach my $exon2 (@{$Aoexon2}){
								my $exon1_p0 = $exon1->[0];
								my $exon1_p1 = $exon1->[1];
								my $exon2_p0 = $exon2->[0];
								my $exon2_p1 = $exon2->[1];
								
								my $dp0 = abs($exon1_p0 - $exon2_p0);
								my $dp1 = abs($exon1_p1 - $exon2_p1);
								
								if($dp0 <= 400 && $dp1 <= 400){
									$hjudge->{$id1} = 1;
									last;
								}
							}
							
							if($hjudge->{$id1}){
								last;
							}
						}
					}
				}
			}
		}
		
		my $num_true = 0;
		foreach my $id1 (@ID1){
			if($hjudge->{$id1}){
				$num_true += 1;
			}
		}
		
		if($num_true == @ID1){
			$r .= $hhit->{$hit1}{id1}."\t".$hit1."\t".$hhit->{$hit1}{cnt}."\t".$ID1[0]."\n";
			#[0] list of genes to be integrated in gff1
			#[1] hint gene from gff2, it will be integrated to gff1
			#[2] num of genes in [0]
			#[3] altname
			
			$cntm++;
			$sum += $hhit->{$hit1}{cnt};
		}
	}
}

print "! [$sum] can be merged to [$cntm] in [$gff1] based on [$gff2]\n";

return $r;
}


#-------------------------------------------------------------------------------
sub Gff_read_as_hash{
my $gff3 = shift;

print "! reading information from [$gff3]...\n";
open(my $fh, "<", $gff3) or die;
my $hash = {};
my $htidall = {};
my $gene = {};
my $gid2ch = {};
my $gid2pos0 = {};
my $gid2pos1 = {};
my $gid2strand = {};
my $gid2gnum = {};
my $hgnum = {};
my $tid2gid = {};
my $gid2tid = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	$line =~ s/\"//g;
	
	unless($line){
		next;
	}
	elsif($line =~ /\#/){
		next;
	}
	
	my @A = split(/\t/, $line);
	my @tag = split(/\;/, $A[8]);
	
	if($A[2] eq 'gene'){
		$hgnum->{$A[0]} += 1;
		
		my $gid;
		foreach my $val (@tag){
			if($val =~ /ID\=/){
				$gid = $val;
				$gid =~ s/ID\=//;
				last;
			}
		}
		
		$gene->{$gid} = join("\t", @A)."\n";
		$gid2ch->{$gid} = $A[0];
		$gid2pos0->{$gid} = $A[3];
		$gid2pos1->{$gid} = $A[4];
		$gid2strand->{$gid} = $A[6];
		$gid2gnum->{$gid} = $hgnum->{$A[0]};
	}
	elsif($A[2] eq 'transcript' || $A[2] eq 'mRNA'){
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
		
		$tid2gid->{$tid} = $gid;
		$gid2tid->{$gid} .= $tid."\n";
		$htidall->{$tid} .= join("\t", @A)."\n";
	}
	elsif($A[2] eq 'exon'){
		my $tid;
		foreach my $val (@tag){
			if($val =~ /Parent\=/){
				$tid = $val;
				$tid =~ s/Parent\=//;
			}
		}
		
		unless($tid){
			foreach my $val (@tag){
				if($val =~ /grp\=/){
					$tid = $val;
					$tid =~ s/grp\=//;
				}
			}
		}
		
		$hash->{$tid}{exon} .= $line."\n";
		$htidall->{$tid} .= join("\t", @A)."\n";
	}
	elsif($A[2] eq 'CDS'){
		my $tid;
		foreach my $val (@tag){
			if($val =~ /Parent\=/){
				$tid = $val;
				$tid =~ s/Parent\=//;
			}
		}
		
		$hash->{$tid}{CDS} .= $line."\n";
		$htidall->{$tid} .= join("\t", @A)."\n";
	}
	elsif($A[2] eq 'start_codon'){
		my $tid;
		foreach my $val (@tag){
			if($val =~ /Parent\=/){
				$tid = $val;
				$tid =~ s/Parent\=//;
			}
		}
		
		$htidall->{$tid} .= join("\t", @A)."\n";
	}
	elsif($A[2] eq 'stop_codon'){
		my $tid;
		foreach my $val (@tag){
			if($val =~ /Parent\=/){
				$tid = $val;
				$tid =~ s/Parent\=//;
			}
		}
		
		$htidall->{$tid} .= join("\t", @A)."\n";
	}
}
close $fh;

my @TID = keys(%{$hash});
my $numTID = @TID;
my $hCDS = {};
my $hexon = {};
foreach my $tid (@TID){
	my @tmp = split(/\./, $tid);
	my $gid = $tmp[0];
	
	unless($gff3 =~ /\.E\.gff/){
		if($hash->{$tid}{CDS}){
			my @CDS = split(/\n/, $hash->{$tid}{CDS});
			
			my $AoCDS = [];
			foreach my $line (@CDS){
				my @A = split(/\t/, $line);
				my @B;
				$B[0] = $A[3];
				$B[1] = $A[4];
				$B[2] = $tid;
				$B[3] = $A[6];		#strand
				$B[4] = $A[0];		#chromosome
				push(@{$AoCDS}, \@B);
			}
			@{$AoCDS} = sort {$a->[0] <=> $b->[0]} @{$AoCDS};
	#		$hCDS->{$gid}{$tid} = $AoCDS;
			$hCDS->{$tid} = $AoCDS;
		}
	}
	
	my @LE = split(/\n/, $hash->{$tid}{exon});
	
	my $AoLE = [];
	foreach my $line (@LE){
		my @A = split(/\t/, $line);
		my @B;
		$B[0] = $A[3];
		$B[1] = $A[4];
		$B[2] = $tid;
		$B[3] = $A[6];		#strand
		$B[4] = $A[0];		#chromosome
		push(@{$AoLE}, \@B);
	}
	@{$AoLE} = sort {$a->[0] <=> $b->[0]} @{$AoLE};
#	$hexon->{$gid}{$tid} = $AoLE;
	$hexon->{$tid} = $AoLE;
}

print "! [$numTID] transcripts with CDS/exon\n";

my $rh = {};
$rh->{hCDS} = $hCDS;
$rh->{hexon} = $hexon;
$rh->{htidall} = $htidall;
$rh->{gene} = $gene;
$rh->{gid2ch} = $gid2ch;
$rh->{gid2pos0} = $gid2pos0;
$rh->{gid2pos1} = $gid2pos1;
$rh->{gid2strand} = $gid2strand;
$rh->{gid2gnum} = $gid2gnum;
$rh->{tid2gid} = $tid2gid;
$rh->{gid2tid} = $gid2tid;

return $rh;
}


#-------------------------------------------------------------------------------
sub Convert_tid2gid{
my $hash = shift;

my @ID = keys(%{$hash});
@ID = sort {$a cmp $b} @ID;

my $nr = {};
foreach my $tid (@ID){
	my $gid = $tid;
	if($tid =~ /Glyma\./){
		my $pos0 = index($tid, "\.", 0);
		my $pos1 = index($tid, "\.", $pos0 + 1);
		$gid = substr($tid, 0, $pos1);
	}
	elsif($tid =~ /Glyma/){
		my $pos0 = index($tid, "\.", 0);
		$gid = substr($tid, 0, $pos0);
	}
	elsif($tid =~ /AT/){
		my $pos0 = index($tid, "\.", 0);
		$gid = substr($tid, 0, $pos0);
	}
	elsif($tid =~ /MELO/ && $tid =~ /P/){
		my $pos0 = index($tid, "P", 0);
		$gid = substr($tid, 0, $pos0);
	}
	elsif($tid =~ /MELO/ && $tid =~ /T/){
		my $pos0 = index($tid, "T", 0);
		$gid = substr($tid, 0, $pos0);
	}
	elsif($tid =~ /MELO/ && $tid =~ /\./){
		my $pos0 = index($tid, "\.", 0);
		$gid = substr($tid, 0, $pos0);
		$gid .= ".2";
	}
	elsif($tid =~ /\./){
		my $pos0 = index($tid, "\.", 0);
		$gid = substr($tid, 0, $pos0);
	}
	$nr->{$gid} = 1;
}

return $nr;
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
		
		my $pcnt_avr = sprintf("%.2f", $sum_pcnt / $cnt);
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
		push(@{$AoA}, \@tmp);
	}
	
	@{$AoA} = sort {$a->[3] <=> $b->[3]} @{$AoA};
	@{$AoA} = sort {$b->[1] <=> $a->[1]} @{$AoA};
	@{$AoA} = sort {$b->[2] <=> $a->[2]} @{$AoA};
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
my $verbose = shift;

system("bash $sh_file\n");

unless($verbose && $verbose eq 'true'){
	print "! thread [$num] completed.\n";
}

return $num;
}


#-------------------------------------------------------------------------------
sub Split_query{
my $query = shift;
my $core = shift;
my $qstr = shift;
my $verbose = shift;

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

unless($verbose && $verbose eq 'true'){
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


#----------------------------------------------------------
sub ZIP{
my $zip = shift;
my $file = shift;

if($zip && -e $zip){
	system("rm $zip");
}
if($file && -e $file){
	system("zip $zip $file");
}

}


#-------------------------------------------------------------------------------
sub Fasta2prefix{
my $prefix = shift;

if($prefix =~ /\.fasta/){
	$prefix =~ s/\.fasta//;
}
elsif($prefix =~ /\.fa/){
	$prefix =~ s/\.fa//;
}

return $prefix;
}


#-------------------------------------------------------------------------------
sub Gff2prefix{
my $prefix = shift;

if($prefix =~ /\.gff3/){
	$prefix =~ s/\.gff3//;
}
elsif($prefix =~ /\.gff/){
	$prefix =~ s/\.gff//;
}

return $prefix;
}


#-------------------------------------------------------------------------------
sub Rmfiles{
my $dir = shift;

my @GB0 = glob("$dir/*");
foreach my $file (@GB0){
	if($file =~ /bN_db-/ || $file =~ /bP_db-/){
		next;
	}
	system("rm $file");
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

