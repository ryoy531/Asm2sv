#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use FindBin;

my $bin_gffread = shift;
my $script_path = shift;
my $dbfasta = shift;
my $dbgff = shift;
my $qfasta = shift;
my $qgff = shift;
my $asm2sv_tsv = shift;
my $asm2sv_gff = shift;
my $cpu = shift;
my $rdir = shift;
my $final_rfile = shift;
my $final_rfile_IDalias = shift;
my $strict_chr_identical = shift;
my $chrinfo = shift;
my $gid_alias = shift;

unless($script_path){
	goto END;
}

my $wpath = getcwd();
my $bin_Blast_BDH = $script_path."/scripts/Blast_BDH.pl";
my $bin_Blat2hint_mRNA = $script_path."/scripts/Blat2hint_mRNA.pl";
my $bin_Find_orthologue_partner = $script_path."/scripts/Find_orthologue_partner.pl";
my $bin_blat = $script_path."/bin/blat";
my $bin_blat2hint = $script_path."/bin/blat2hints.pl";

my $err = 0;
$err += FileCheck1($bin_Blast_BDH, "Blast_BDH.pl");
$err += FileCheck1($bin_Blat2hint_mRNA, "Blat2hint_mRNA.pl");
$err += FileCheck1($bin_Find_orthologue_partner, "Find_orthologue_partner.pl");
$err += FileCheck1($bin_blat, "blat");
$err += FileCheck1($bin_blat2hint, "blat2hints.pl");
$err += TestBlat($bin_blat);
$err += FileCheck2($dbfasta, "dbfasta");
$err += FileCheck2($dbgff, "dbgff");
$err += FileCheck2($qfasta, "qfasta");
$err += FileCheck2($qgff, "qgff");
$err += FileCheck2($asm2sv_tsv, "Asm2sv rfile_final");
$err += FileCheck2($asm2sv_gff, "Asm2sv rgff_final");
if(! $rdir){
	print "! missing output directory...\n";
}
if(! $bin_gffread || ! -e $bin_gffread){
	print "! missing gffread...\n";
}
if(! $bin_Blast_BDH || ! -e $bin_Blast_BDH){
	print "! missing Blast_BDH.pl ...\n";
}
if($err > 0){
	print "! abort script...\n";
	die;
}
if(! $chrinfo || ! -e $chrinfo){
	$chrinfo = "null";
}
if(! $gid_alias || ! -e $gid_alias){
	$gid_alias = "null";
}

unless(-e $rdir){
	system("mkdir $rdir");
}
chdir $rdir;

LnFile($dbfasta);
LnFile($dbgff);
LnFile($qfasta);
LnFile($qgff);
LnFile($asm2sv_tsv);
LnFile($asm2sv_gff);
LnFile($chrinfo);
LnFile($gid_alias);

my $dbpref = Path2pref($dbfasta);
my $qpref = Path2pref($qfasta);

my $log_dir = "log";
my $log_ess0 = "./log/log_find_partner_round1.txt";
unless(-d $log_dir){
	system("mkdir $log_dir");
}

my $rdir_find0 = "_Find_orthologue_partner";
my $rfile_find0 = "./$rdir_find0/list_partner_ID_genebase_c4.tsv";
print "! step.1 analyze gene partnership with [$bin_Find_orthologue_partner]... (may take a while)\n";
if(! -e $rfile_find0){
	print " - log file = [$log_ess0]\n";
	my $cmd_ess0 = "$bin_Find_orthologue_partner $rdir_find0 $bin_gffread $bin_Blat2hint_mRNA $bin_blat $bin_blat2hint $bin_Blast_BDH $dbfasta $qfasta $dbgff $qgff 20 4 $cpu > $log_ess0 2>&1";
	#my $cmd_ess0 = "$bin_Find_orthologue_partner $rdir_find0 $bin_gffread $bin_Blat2hint_mRNA $bin_blat $bin_blat2hint $bin_Blast_BDH $dbfasta $qfasta $dbgff $qgff 20 4 $cpu";
	if(system($cmd_ess0) != 0){
		print "! failed.\n";
	}
	else{
		print "! done.\n";
	}
}
else{
	print "! [$rfile_find0] already exists, skip...\n";
}
print "! step.1 completed.\n";

my $db_transcript = $dbpref."_transcript.fasta";
my $db_cds = $dbpref."_cds.fasta";
my $db_protein = $dbpref."_protein.fasta";
my $q_transcript = $qpref."_transcript.fasta";
my $q_cds = $qpref."_cds.fasta";
my $q_protein = $qpref."_protein.fasta";

print "! step.2 compare [$dbfasta] and [$qfasta] ...\n";
$err += Gffread2seq($bin_gffread, $dbfasta, $dbgff, $db_transcript, $db_cds, $db_protein);
$err += Gffread2seq($bin_gffread, $qfasta, $qgff, $q_transcript, $q_cds, $q_protein);

if($err > 0){
	print "! abort script due to error...\n";
	die;
}

my $cpun = $cpu;
if($cpun > 16){
	$cpun = 16;
}

my $log_bdhn = "./$log_dir/log_Blast_BDHn.txt";
my $log_bdhp = "./$log_dir/log_Blast_BDHp.txt";
my $rdir_bdhn = "_blastBDHn";
my $rdir_bdhp = "_blastBDHp";
my $bdhn_th_pcnt = 70;
my $bdhn_th_pcov = 1;
my $bdhn_th_eval = "1e-30";
my $bdhp_th_pcnt = 70;
my $bdhp_th_pcov = 1;
my $bdhp_th_eval = "1e-30";
my $rfile_bdhn = "./$rdir_bdhn/BDconfidenthit_pcnt-".$bdhn_th_pcnt."_pcov-".$bdhn_th_pcov."_Eval-".$bdhn_th_eval.".csv";
my $rfile_bdhp = "./$rdir_bdhp/BDconfidenthit_pcnt-".$bdhp_th_pcnt."_pcov-".$bdhp_th_pcov."_Eval-".$bdhp_th_eval.".csv";
my $cmdn = "$bin_Blast_BDH $rdir_bdhn $dbfasta $qfasta $dbgff $qgff $db_transcript $q_transcript 1 20 $bdhn_th_pcnt $bdhn_th_pcov $bdhn_th_eval $cpun y n n > $log_bdhn 2>&1";
my $cmdp = "$bin_Blast_BDH $rdir_bdhp $dbfasta $qfasta $dbgff $qgff $db_protein $q_protein 2 20 $bdhp_th_pcnt $bdhp_th_pcov $bdhp_th_eval $cpu y n n > $log_bdhp 2>&1";

unless(-e $rfile_bdhn){
	print "! reciprocal BLASTn with [$cpun] threads...\n";
	if(system($cmdn) != 0){
		print "! failed.\n";
		die;
	}
}
else{
	print "! [$rfile_bdhn] already exists, skip BLASTn...\n";
}
Rmfiles("$rdir_bdhn");

unless(-e $rfile_bdhp){
	print "! reciprocal BLASTp with [$cpu] threads...\n";
	if(system($cmdp) != 0){
		print "! failed.\n";
		die;
	}
}
else{
	print "! [$rfile_bdhp] already exists, skip BLASTp...\n";
}
Rmfiles("$rdir_bdhp");

print "! collecting information...\n";
my $HGFFdb = Gff_to_hash($dbgff);
my $HGFFq = Gff_to_hash($qgff);
my $hgffdb = $HGFFdb->{hgff};
my $hgffq = $HGFFq->{hgff};

my $hBDH = {};
$hBDH = Read_BDconfident($hBDH, $rfile_bdhn, "nuc");
$hBDH = Read_BDconfident($hBDH, $rfile_bdhp, "prot");

my $cpu_np = $cpu;
if($cpu_np > 8){
	$cpu_np = 8;
}

print "! searching for pairs based on [$asm2sv_tsv] with [$cpu_np] threads\n";
my $hseqinfo = {};
if($chrinfo ne 'null' && -e $chrinfo){
	print " - [$chrinfo] is specified as seq ID alias information.\n";
	my $htmp_chrinfo = Read_seqinfo($chrinfo);
	$hseqinfo = $htmp_chrinfo->{seqID};
}

my $rscfile = "fitpos_info.tsv";
my $DBG = $HGFFdb->{GID};
my $num_DBG = @{$DBG};
if(! -e $rscfile){
	my $screening = "dbgene\tqgene\tfit\tsum_score (BLASTn)\tsum_score (BLASTp)\tAsm2sv_judge\n";
	my $th_flanking = 5000;
	my $npart = int($num_DBG / $cpu_np) + 1;
	my $thrs = [];
	for(my $i = 0; $i < $cpu_np; $i++){
		my $p0 = $i * $npart;
		my $p1 = ($i + 1) * $npart;
		my $thr = threads->new(\&FitPos0, $DBG, $p0, $p1, $rfile_bdhn, $rfile_bdhp, $asm2sv_tsv, $hgffdb, $hgffq, $hseqinfo, $dbpref, $qpref, $th_flanking, 1, $strict_chr_identical, $i);
		push(@{$thrs}, $thr);
	}
	if(@{$thrs}){
		foreach my $thr (@{$thrs}){
			$screening .= $thr->join();
		}
	}
	SAVE($rscfile, $screening);
}
else{
	print "! [$rscfile] already exists, skip...\n";
}

my ($hfitdb, $hfitq, $record) = Fitinfo2hash($rscfile);
my $rec_rfile = "fitpos_sorted.tsv";
SAVE($rec_rfile, $record);
print "! output [$rec_rfile]\n";

print "! integrating results...\n";
my $hpartnerR2 = {};
foreach my $dbgene (@{$DBG}){
	if($hfitdb->{$dbgene}){
		my $qcand = $hfitdb->{$dbgene};
		if($hfitq->{$qcand} && $hfitq->{$qcand} eq $dbgene){
			$hpartnerR2->{$dbgene}{$qcand} = 1;
		}
	}
}

my ($rev_result, $rev_alias, $cnt_alias, $rev_contradict) = Integrate2R($rfile_find0, $hpartnerR2, $hfitdb, $hfitq, $hBDH, $hgffdb, $hgffq, $hseqinfo, $dbpref, $qpref, $gid_alias);

if(! $final_rfile || $final_rfile eq 'null'){
	$final_rfile = $asm2sv_tsv;
	$final_rfile =~ s/\.tsv//;
	$final_rfile .= "_partner.tsv";
}
if(! $final_rfile_IDalias || $final_rfile_IDalias eq 'null'){
	$final_rfile_IDalias = $asm2sv_tsv;
	$final_rfile_IDalias =~ s/\.tsv//;
	$final_rfile_IDalias .= "_partner_aliasID.tsv";
}
if($rev_contradict){
	SAVE("contradicting_seqID.tsv", $rev_contradict);
}

chdir $wpath;

SAVE($final_rfile, $rev_result);
print "! output [$final_rfile]\n";

if($cnt_alias > 0){
	SAVE($final_rfile_IDalias, $rev_alias);
	print "! output [$final_rfile_IDalias]\n";
}

END:{
	my $end = "";
}


################################################################################
#-------------------------------------------------------------------------------
sub Read_seqinfo{
my $file = shift;

#print "! reading seq ID alias info [$file]...\n";
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

#print "! [$cnt] lines\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Integrate2R{
my ($rfile_find0, $hpartnerR2, $hfitdb, $hfitq, $hBDH, $hgffdb, $hgffq, $hseqinfo, $dbpref, $qpref, $gid_alias) = @_;

my $halias = {};
if($gid_alias && $gid_alias ne 'null' && -e $gid_alias){
	print " - [$gid_alias] is specified as gene ID alias\n";
	open(my $fh, "<", $gid_alias);
	my $icns = -1;
	my $iori = -1;
	my $cnt_skip = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		unless($line){
			next;
		}
		my @A = split(/,/, $line);
		my $numA = @A;
		if($line =~ /qID_consensus/ && $line =~ /qID_original/){
			for(my $i = 0; $i < $numA; $i++){
				if($A[$i] eq 'qID_consensus'){
					$icns = $i;
				}
				elsif($A[$i] eq 'qID_original'){
					$iori = $i;
				}
			}
		}
		else{
			if($icns < 0 || ! $A[$icns]){
				print " - missing column for 'qID_consensus', skip...\n";
				return;
			}
			if($iori < 0 || ! $A[$iori]){
				print " - missing column for 'qID_original', skip...\n";
				return;
			}
			if(! $A[$icns] || $A[$icns] eq '-' || $A[$icns] eq 'NA' || $A[$icns] eq 'na' || $A[$icns] eq 'N/A'){
				$A[$icns] = $A[$iori];
			}
			if(! $halias->{$A[$icns]}){
				$halias->{$A[$icns]} = $A[$iori];		# convert qID_consensus to qID_original
			}
			if($A[$iori] && $A[$iori] =~ /GlymaJER/ && $A[0] =~ /GlymaJER/ && $A[0] ne $A[$iori]){
#				print " $A[0] -> $A[$iori]\n";
				$halias->{$A[0]} = $A[$iori];			# convert "old qID_original (A[0])" to qID_original
			}
		}
	}
	close $fh;
}

my $rev = "gene ID\tprotein coding\tcandidate partner\tpartner (genomic position + blast n/p)\tpartner (genomic position + blastn or p)\tblast n/p BDH but different genomic position\tBDH has another partner\tstatus summary1\tstatus summary2\tparalogs\ttandem paralogs\tCNV\tCNV (target vs hint)\tCNV num\tcandidate for consensus ID\tAlt-ID in replace.csv\thomolog remark\t%identity (t)\t%alignment coverage (t)\tE-value (t)\tscore (t)\tBLASTn_status\t%identity (p)\t%alignment coverage (p)\tE-value (p)\tscore (p)\tBLASTp_status\tAsm2sv_judge\tgene1 Chr\tgene1 alias-Chr\tgene1 pos0\tgene1 pos1\tgene1 strand\tgene2 Chr\tgene2 alias-Chr\tgene2 pos0\tgene2 pos1\tgene2 strand\n";
my $rev_alias = $rev;
my $rev_contradict = "";
my $cnt_alias = 0;
my $cnt_seqalias = 0;
my $cnt_contrachr = 0;
if(-e $rfile_find0){
	#print "! reading [$rfile_find0]...\n";
	open(my $fh, "<", $rfile_find0) or die;
	my @L;
	my $hprevq = {};
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		unless($line){
			next;
		}
		elsif($line =~ /gene ID/ && $line =~ /genomic position/){
			next;
		}
		my @A = split(/\t/, $line);
		if($A[8] eq '1-to-1 ortholog' && $A[2] ne '-'){
			$hprevq->{ortho}{$A[2]} = $A[0];
		}
		push(@L, $line);
	}
	close $fh;

	my $cnt_ortho_confirmed = 0;
	my $cnt_replace_ortho1 = 0;
	my $cnt_replace_ortho2 = 0;
	my $cnt_homolog_confirmed = 0;
	my $cnt_replace_homolog1 = 0;
	my $cnt_replace_homolog2 = 0;
	my $cnt_false_ortho = 0;
	my $cnt_false_homolog = 0;
	my $cnt_NA_keep = 0;
	my $nrcheck = {};
	my $cnt_overlap = 0;
	foreach my $line (@L){
		my @A0 = split(/\t/, $line);
		my @A;
		for(my $i = 0; $i < 17; $i++){		# truncate column number to 17
			if(defined $A0[$i]){
				push(@A, $A0[$i]);
			}
			else{
				push(@A, "-");
			}
		}
		
#		$hBDH
#		$hash->{$mode}{db2q}{$A[1]}{$A[4]}{identity} = $A[10];
#		$hash->{$mode}{db2q}{$A[1]}{$A[4]}{pcov} = $A[12];
#		$hash->{$mode}{db2q}{$A[1]}{$A[4]}{score} = $A[14];
#		$hash->{$mode}{db2q}{$A[1]}{$A[4]}{Eval} = $A[16];
		
		my $each_r = "";
		my $pattern = "null";
		if($A[8]){
			my $tophit_of_db = "-";
			my $tophit_of_q = "-";
			if($A[8] eq '1-to-1 ortholog' && $A[2] ne '-'){
				if($hfitdb->{$A[0]}){
					$tophit_of_db = $hfitdb->{$A[0]};
				}
				if($hfitq->{$A[2]}){
					$tophit_of_q = $hfitq->{$A[2]};
				}
				
				if($hpartnerR2->{$A[0]}{$A[2]}){
					$tophit_of_db = $A[2];
					$pattern = "true";
					$cnt_ortho_confirmed++;
				}
				elsif($tophit_of_db ne '-'){
					if($hpartnerR2->{$A[0]}{$tophit_of_db}){
						if(! $hprevq->{ortho}{$tophit_of_db}){
							$A[2] = $tophit_of_db;
							$A[3] = $tophit_of_db;
							$A[7] = "blastn/p + same genomic position (Asm2sv)";
							$A[8] = "1-to-1 ortholog (revise)";
							$A[14] = $tophit_of_db.".refversion";
							$A[15] = $A[0];
							$pattern = "assign unlinked";
							$cnt_replace_ortho1++;
						}
						else{
							$A[2] = $tophit_of_db;
							$A[3] = $tophit_of_db;
							$A[7] = "blastn/p + same genomic position (Asm2sv)";
							$A[8] = "1-to-1 ortholog (revise)";
							$A[14] = $tophit_of_db.".refversion";
							$A[15] = $A[0];
							$pattern = "replace with other linked";
							$cnt_replace_ortho2++;
						}
					}
					else{
						$A[2] = "-";
						$A[3] = "-";
						$A[4] = "-";
						$A[5] = "-";
						$A[6] = $tophit_of_db;
						$A[8] = "partner has another pair (revise)";
						$A[14] = $A[0];
						$A[15] = "-";
						$pattern = "remove";
						$cnt_false_ortho++;
					}
				}
				else{
					$A[2] = "-";
					$A[3] = "-";
					$A[4] = "-";
					$A[5] = "-";
					$A[6] = "-";
					$A[7] = "-";
					$A[8] = "-";
					$A[14] = $A[0];
					$A[15] = "-";
					$pattern = "remove";
					$cnt_false_ortho++;
				}
				
				if($pattern ne 'null'){
					if($A[2] ne '-' && ! $nrcheck->{$A[2]}){
						$nrcheck->{$A[2]} = 1;
					}
					elsif($A[2] ne '-'){
						$nrcheck->{$A[2]} += 1;
						$cnt_overlap++;
					}
					
					$each_r = join("\t", @A);
					if($tophit_of_db ne '-' && defined $hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{identity}){
						$each_r .= "\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{identity}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{pcov}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{score}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{Eval}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{status};
					}
					else{
						$each_r .= "\t-\t-\t-\t-\t-";
					}
					if($tophit_of_db ne '-' && defined $hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{identity}){
						$each_r .= "\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{identity}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{pcov}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{score}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{Eval}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{status};
					}
					else{
						$each_r .= "\t-\t-\t-\t-\t-";
					}
					$each_r .= "\t$pattern";
				}
			}
			elsif($A[2] ne '-'){
				if($hfitdb->{$A[0]}){
					$tophit_of_db = $hfitdb->{$A[0]};
				}
				if($hfitq->{$A[2]}){
					$tophit_of_q = $hfitq->{$A[2]};
				}
				
				if($hpartnerR2->{$A[0]}{$A[2]}){
					$tophit_of_db = $A[2];
					$pattern = "true";
					$cnt_homolog_confirmed++;
				}
				elsif($tophit_of_db ne '-'){
					if($hpartnerR2->{$A[0]}{$tophit_of_db}){
						if(! $hprevq->{ortho}{$tophit_of_db}){
							$A[2] = $tophit_of_db;
							$A[3] = $tophit_of_db;
							$A[7] = "blastn/p + same genomic position (Asm2sv)";
							$A[8] = "1-to-1 ortholog (revise)";
							$A[14] = $tophit_of_db.".refversion";
							$A[15] = $A[0];
							$pattern = "assign unlinked";
							$cnt_replace_homolog1++;
						}
						else{
							$A[2] = $tophit_of_db;
							$A[3] = $tophit_of_db;
							$A[7] = "blastn/p + same genomic position (Asm2sv)";
							$A[8] = "1-to-1 ortholog (revise)";
							$A[14] = $tophit_of_db.".refversion";
							$A[15] = $A[0];
							$pattern = "replace with other linked";
							$cnt_replace_homolog2++;
						}
					}
					else{
						$A[2] = "-";
						$A[3] = "-";
						$A[4] = "-";
						$A[5] = "-";
						$A[6] = $tophit_of_db;
						$A[8] = "partner has another pair (revise)";
						$A[14] = $A[0];
						$A[15] = "-";
						$pattern = "remove";
						$cnt_false_homolog++;
					}
				}
				else{
					$A[2] = "-";
					$A[3] = "-";
					$A[4] = "-";
					$A[5] = "-";
					$A[6] = "-";
					$A[7] = "-";
					$A[8] = "-";
					$A[14] = $A[0];
					$A[15] = "-";
					$pattern = "remove";
					$cnt_false_homolog++;
				}
				
				if($pattern ne 'null'){
					if($A[2] ne '-' && ! $nrcheck->{$A[2]}){
						$nrcheck->{$A[2]} = 1;
					}
					elsif($A[2] ne '-'){
						$nrcheck->{$A[2]} += 1;
						$cnt_overlap++;
					}
					
					$each_r = join("\t", @A);
					if($tophit_of_db ne '-' && defined $hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{identity}){
						$each_r .= "\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{identity}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{pcov}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{score}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{Eval}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{status};
					}
					else{
						$each_r .= "\t-\t-\t-\t-\t-";
					}
					if($tophit_of_db ne '-' && defined $hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{identity}){
						$each_r .= "\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{identity}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{pcov}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{score}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{Eval}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{status};
					}
					else{
						$each_r .= "\t-\t-\t-\t-\t-";
					}
					$each_r .= "\t$pattern";
				}
			}
			elsif($A[7] ne '-' && $A[4] ne '-'){
				if($hfitdb->{$A[0]}){
					$tophit_of_db = $hfitdb->{$A[0]};
				}
				if($hfitq->{$A[4]}){
					$tophit_of_q = $hfitq->{$A[4]};
				}
				
				if($hpartnerR2->{$A[0]}{$A[4]}){
					if(! $hprevq->{ortho}{$A[4]}){
						$tophit_of_db = $A[4];
						$A[2] = $tophit_of_db;
						$A[3] = $tophit_of_db;
						$A[4] = "-";
						$A[5] = "-";
						$A[6] = "-";
						$A[7] = "blastn/p + same genomic position (Asm2sv)";
						$A[8] = "1-to-1 ortholog (revise)";
						$A[14] = $tophit_of_db.".refversion";
						$A[15] = $A[0];
						$pattern = "assign unlinked";
						$cnt_replace_ortho1++;
					}
					else{
						$tophit_of_db = $A[4];
						$A[2] = $tophit_of_db;
						$A[3] = $tophit_of_db;
						$A[4] = "-";
						$A[5] = "-";
						$A[6] = "-";
						$A[7] = "blastn/p + same genomic position (Asm2sv)";
						$A[8] = "1-to-1 ortholog (revise)";
						$A[14] = $tophit_of_db.".refversion";
						$A[15] = $A[0];
						$pattern = "replace with other linked";
						$cnt_replace_ortho2++;
					}
				}
				elsif($tophit_of_db ne '-'){
					if($hpartnerR2->{$A[0]}{$tophit_of_db}){
						if(! $hprevq->{ortho}{$tophit_of_db}){
							$A[2] = $tophit_of_db;
							$A[3] = $tophit_of_db;
							$A[4] = "-";
							$A[5] = "-";
							$A[6] = "-";
							$A[7] = "blastn/p + same genomic position (Asm2sv)";
							$A[8] = "1-to-1 ortholog (revise)";
							$A[14] = $tophit_of_db.".refversion";
							$A[15] = $A[0];
							$pattern = "assign unlinked";
							$cnt_replace_ortho1++;
						}
						else{
							$A[2] = $tophit_of_db;
							$A[3] = $tophit_of_db;
							$A[4] = "-";
							$A[5] = "-";
							$A[6] = "-";
							$A[7] = "blastn/p + same genomic position (Asm2sv)";
							$A[8] = "1-to-1 ortholog (revise)";
							$A[14] = $tophit_of_db.".refversion";
							$A[15] = $A[0];
							$pattern = "replace with other linked";
							$cnt_replace_ortho2++;
						}
					}
					else{
						$pattern = "keep";
						$cnt_NA_keep++;
					}
				}
				else{
					$pattern = "keep";
					$cnt_NA_keep++;
				}
				
				if($pattern ne 'keep'){
					if($A[2] ne '-' && ! $nrcheck->{$A[2]}){
						$nrcheck->{$A[2]} = 1;
					}
					elsif($A[2] ne '-'){
						$nrcheck->{$A[2]} += 1;
						$cnt_overlap++;
					}
					
					$each_r = join("\t", @A);
					if($tophit_of_db ne '-' && defined $hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{identity}){
						$each_r .= "\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{identity}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{pcov}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{score}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{Eval}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{status};
					}
					else{
						$each_r .= "\t-\t-\t-\t-\t-";
					}
					if($tophit_of_db ne '-' && defined $hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{identity}){
						$each_r .= "\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{identity}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{pcov}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{score}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{Eval}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{status};
					}
					else{
						$each_r .= "\t-\t-\t-\t-\t-";
					}
					$each_r .= "\t$pattern";
				}
				else{
					$each_r = join("\t", @A);
					$each_r .= "\t$A0[21]\t$A0[22]\t$A0[23]\t$A0[24]\t-";
					$each_r .= "\t$A0[25]\t$A0[26]\t$A0[27]\t$A0[28]\t-";
					$each_r .= "\t$pattern";
				}
			}
			elsif($A[7] ne '-' && $A[5] ne '-'){
				if($hfitdb->{$A[0]}){
					$tophit_of_db = $hfitdb->{$A[0]};
				}
				if($hfitq->{$A[5]}){
					$tophit_of_q = $hfitq->{$A[5]};
				}
				
				if($hpartnerR2->{$A[0]}{$A[5]}){
					if(! $hprevq->{ortho}{$A[5]}){
						$tophit_of_db = $A[5];
						$A[2] = $tophit_of_db;
						$A[3] = $tophit_of_db;
						$A[4] = "-";
						$A[5] = "-";
						$A[6] = "-";
						$A[7] = "blastn/p + same genomic position (Asm2sv)";
						$A[8] = "1-to-1 ortholog (revise)";
						$A[14] = $tophit_of_db.".refversion";
						$A[15] = $A[0];
						$pattern = "assign unlinked";
						$cnt_replace_ortho1++;
					}
					else{
						$tophit_of_db = $A[5];
						$A[2] = $tophit_of_db;
						$A[3] = $tophit_of_db;
						$A[4] = "-";
						$A[5] = "-";
						$A[6] = "-";
						$A[7] = "blastn/p + same genomic position (Asm2sv)";
						$A[8] = "1-to-1 ortholog (revise)";
						$A[14] = $tophit_of_db.".refversion";
						$A[15] = $A[0];
						$pattern = "replace with other linked";
						$cnt_replace_ortho2++;
					}
				}
				elsif($tophit_of_db ne '-'){
					if($hpartnerR2->{$A[0]}{$tophit_of_db}){
						if(! $hprevq->{ortho}{$tophit_of_db}){
							$A[2] = $tophit_of_db;
							$A[3] = $tophit_of_db;
							$A[4] = "-";
							$A[5] = "-";
							$A[6] = "-";
							$A[7] = "blastn/p + same genomic position (Asm2sv)";
							$A[8] = "1-to-1 ortholog (revise)";
							$A[14] = $tophit_of_db.".refversion";
							$A[15] = $A[0];
							$pattern = "assign unlinked";
							$cnt_replace_ortho1++;
						}
						else{
							$A[2] = $tophit_of_db;
							$A[3] = $tophit_of_db;
							$A[4] = "-";
							$A[5] = "-";
							$A[6] = "-";
							$A[7] = "blastn/p + same genomic position (Asm2sv)";
							$A[8] = "1-to-1 ortholog (revise)";
							$A[14] = $tophit_of_db.".refversion";
							$A[15] = $A[0];
							$pattern = "replace with other linked";
							$cnt_replace_ortho2++;
						}
					}
					else{
						$pattern = "keep";
						$cnt_NA_keep++;
					}
				}
				else{
					$pattern = "keep";
					$cnt_NA_keep++;
				}
				
				if($pattern ne 'keep'){
					if($A[2] ne '-' && ! $nrcheck->{$A[2]}){
						$nrcheck->{$A[2]} = 1;
					}
					elsif($A[2] ne '-'){
						$nrcheck->{$A[2]} += 1;
						$cnt_overlap++;
					}
					
					$each_r = join("\t", @A);
					if($tophit_of_db ne '-' && defined $hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{identity}){
						$each_r .= "\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{identity}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{pcov}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{score}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{Eval}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{status};
					}
					else{
						$each_r .= "\t-\t-\t-\t-\t-";
					}
					if($tophit_of_db ne '-' && defined $hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{identity}){
						$each_r .= "\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{identity}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{pcov}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{score}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{Eval}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{status};
					}
					else{
						$each_r .= "\t-\t-\t-\t-\t-";
					}
					$each_r .= "\t$pattern";
				}
				else{
					$each_r = join("\t", @A);
					$each_r .= "\t$A0[21]\t$A0[22]\t$A0[23]\t$A0[24]\t-";
					$each_r .= "\t$A0[25]\t$A0[26]\t$A0[27]\t$A0[28]\t-";
					$each_r .= "\t$pattern";
				}
			}
			elsif($A[7] ne '-' && $A[6] ne '-'){
				if($hfitdb->{$A[0]}){
					$tophit_of_db = $hfitdb->{$A[0]};
				}
				if($hfitq->{$A[6]}){
					$tophit_of_q = $hfitq->{$A[6]};
				}
				
				if($hpartnerR2->{$A[0]}{$A[6]}){
					if(! $hprevq->{ortho}{$A[6]}){
						$tophit_of_db = $A[6];
						$A[2] = $tophit_of_db;
						$A[3] = $tophit_of_db;
						$A[4] = "-";
						$A[5] = "-";
						$A[6] = "-";
						$A[7] = "blastn/p + same genomic position (Asm2sv)";
						$A[8] = "1-to-1 ortholog (revise)";
						$A[14] = $tophit_of_db.".refversion";
						$A[15] = $A[0];
						$pattern = "assign unlinked";
						$cnt_replace_ortho1++;
					}
					else{
						$tophit_of_db = $A[6];
						$A[2] = $tophit_of_db;
						$A[3] = $tophit_of_db;
						$A[4] = "-";
						$A[5] = "-";
						$A[6] = "-";
						$A[7] = "blastn/p + same genomic position (Asm2sv)";
						$A[8] = "1-to-1 ortholog (revise)";
						$A[14] = $tophit_of_db.".refversion";
						$A[15] = $A[0];
						$pattern = "replace with other linked";
						$cnt_replace_ortho2++;
					}
				}
				elsif($tophit_of_db ne '-'){
					if($hpartnerR2->{$A[0]}{$tophit_of_db}){
						if(! $hprevq->{ortho}{$tophit_of_db}){
							$A[2] = $tophit_of_db;
							$A[3] = $tophit_of_db;
							$A[4] = "-";
							$A[5] = "-";
							$A[6] = "-";
							$A[7] = "blastn/p + same genomic position (Asm2sv)";
							$A[8] = "1-to-1 ortholog (revise)";
							$A[14] = $tophit_of_db.".refversion";
							$A[15] = $A[0];
							$pattern = "assign unlinked";
							$cnt_replace_ortho1++;
						}
						else{
							$A[2] = $tophit_of_db;
							$A[3] = $tophit_of_db;
							$A[4] = "-";
							$A[5] = "-";
							$A[6] = "-";
							$A[7] = "blastn/p + same genomic position (Asm2sv)";
							$A[8] = "1-to-1 ortholog (revise)";
							$A[14] = $tophit_of_db.".refversion";
							$A[15] = $A[0];
							$pattern = "replace with other linked";
							$cnt_replace_ortho2++;
						}
					}
					else{
						$pattern = "keep";
						$cnt_NA_keep++;
					}
				}
				else{
					$pattern = "keep";
					$cnt_NA_keep++;
				}
				
				if($pattern ne 'keep'){
					if($A[2] ne '-' && ! $nrcheck->{$A[2]}){
						$nrcheck->{$A[2]} = 1;
					}
					elsif($A[2] ne '-'){
						$nrcheck->{$A[2]} += 1;
						$cnt_overlap++;
					}
					
					$each_r = join("\t", @A);
					if($tophit_of_db ne '-' && defined $hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{identity}){
						$each_r .= "\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{identity}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{pcov}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{score}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{Eval}."\t".$hBDH->{nuc}{db2q}{$A[0]}{$tophit_of_db}{status};
					}
					else{
						$each_r .= "\t-\t-\t-\t-\t-";
					}
					if($tophit_of_db ne '-' && defined $hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{identity}){
						$each_r .= "\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{identity}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{pcov}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{score}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{Eval}."\t".$hBDH->{prot}{db2q}{$A[0]}{$tophit_of_db}{status};
					}
					else{
						$each_r .= "\t-\t-\t-\t-\t-";
					}
					$each_r .= "\t$pattern";
				}
				else{
					$each_r = join("\t", @A);
					$each_r .= "\t$A0[21]\t$A0[22]\t$A0[23]\t$A0[24]\t-";
					$each_r .= "\t$A0[25]\t$A0[26]\t$A0[27]\t$A0[28]\t-";
					$each_r .= "\t$pattern";
				}
			}
			else{
				$each_r = join("\t", @A);
				$each_r .= "\t$A0[21]\t$A0[22]\t$A0[23]\t$A0[24]\t-";
				$each_r .= "\t$A0[25]\t$A0[26]\t$A0[27]\t$A0[28]\t-";
				$each_r .= "\tkeep";
			}
		}
		if($each_r){
			my @B = split(/\t/, $each_r);
			my $add_info = "";
			my $gene1_chr = "-";
			my $gene2_chr = "-";
			if($hgffdb->{$B[0]}){
				if($hgffdb->{$B[0]}{strand} eq '+'){
					$hgffdb->{$B[0]}{strand} = "plus";
				}
				elsif($hgffdb->{$B[0]}{strand} eq '-'){
					$hgffdb->{$B[0]}{strand} = "minus";
				}
				$gene1_chr = $hgffdb->{$B[0]}{Chr};
				if($hseqinfo->{$dbpref}{$gene1_chr}){
					$gene1_chr = $hseqinfo->{$dbpref}{$gene1_chr};
					$cnt_seqalias++;
				}
				else{
					$gene1_chr = "-";
				}
				$add_info .= "\t".$hgffdb->{$B[0]}{Chr}."\t".$gene1_chr."\t".$hgffdb->{$B[0]}{pos0}."\t".$hgffdb->{$B[0]}{pos1}."\t".$hgffdb->{$B[0]}{strand};
			}
			else{
				$add_info .= "\t-\t-\t-\t-";
			}
			if($B[2] ne '-' && $hgffq->{$B[2]}){
				if($hgffq->{$B[2]}{strand} eq '+'){
					$hgffq->{$B[2]}{strand} = "plus";
				}
				elsif($hgffq->{$B[2]}{strand} eq '-'){
					$hgffq->{$B[2]}{strand} = "minus";
				}
				$gene2_chr = $hgffq->{$B[2]}{Chr};
				if($hseqinfo->{$qpref}{$gene2_chr}){
					$gene2_chr = $hseqinfo->{$qpref}{$gene2_chr};
					$cnt_seqalias++;
				}
				else{
					$gene2_chr = "-";
				}
				$add_info .= "\t".$hgffq->{$B[2]}{Chr}."\t".$gene2_chr."\t".$hgffq->{$B[2]}{pos0}."\t".$hgffq->{$B[2]}{pos1}."\t".$hgffq->{$B[2]}{strand};
			}
			else{
				$add_info .= "\t-\t-\t-\t-";
			}
			
			if($gene1_chr ne '-' && $gene2_chr ne '-' && $gene1_chr ne $gene2_chr){
				$rev_contradict .= $each_r.$add_info."\n";
				$cnt_contrachr++;
			}
			
			my @Ba = @B;
			if($B[6] && $B[6] ne '-'){
				my @Hmtmp = split(/,/, $B[6]);
				my @rtmp;
				my @rtmpa;
				foreach my $hmgene (@Hmtmp){
					if($hmgene =~ /\;/){
						my @subhmtmp = split(/\;/, $hmgene);
						$hmgene = $subhmtmp[0];
					}
					if($hgffdb->{$hmgene}{Chr}){
						my $ginfo = $hmgene.";".$hgffdb->{$hmgene}{Chr}.";".$hgffdb->{$hmgene}{pos0}.";".$hgffdb->{$hmgene}{pos1}.";".$hgffdb->{$hmgene}{strand};
						push(@rtmp, $ginfo);			# Glyma.01G000174.v4.JE1;chr01;92087;92955;-,GlymaJER.01G000001.1;chr01;88551;93728;-
						
						if($halias->{$hmgene} && $halias->{$hmgene} ne '-'){
							$ginfo = $halias->{$hmgene}.";".$hgffdb->{$hmgene}{Chr}.";".$hgffdb->{$hmgene}{pos0}.";".$hgffdb->{$hmgene}{pos1}.";".$hgffdb->{$hmgene}{strand};
							push(@rtmpa, $ginfo);
						}
					}
					elsif($hgffq->{$hmgene}{Chr}){
						my $ginfo = $hmgene.";".$hgffq->{$hmgene}{Chr}.";".$hgffq->{$hmgene}{pos0}.";".$hgffq->{$hmgene}{pos1}.";".$hgffq->{$hmgene}{strand};
						push(@rtmp, $ginfo);			# Glyma.01G000174.v4.JE1;chr01;92087;92955;-,GlymaJER.01G000001.1;chr01;88551;93728;-
						push(@rtmpa, $ginfo);
					}
				}
				if(@rtmp){
					$B[6] = join(",", @rtmp);
				}
				if(@rtmpa){
					$Ba[6] = join(",", @rtmpa);
				}
			}
			if($B[9] && $B[9] ne '-'){
				my @Hmtmp = split(/,/, $B[9]);
				my @rtmp;
				my @rtmpa;
				foreach my $hmgene (@Hmtmp){
					if($hmgene =~ /\;/){
						my @subhmtmp = split(/\;/, $hmgene);
						$hmgene = $subhmtmp[0];
					}
					if($hgffdb->{$hmgene}{Chr}){
						my $ginfo = $hmgene.";".$hgffdb->{$hmgene}{Chr}.";".$hgffdb->{$hmgene}{pos0}.";".$hgffdb->{$hmgene}{pos1}.";".$hgffdb->{$hmgene}{strand};
						push(@rtmp, $ginfo);			# Glyma.01G000174.v4.JE1;chr01;92087;92955;-,GlymaJER.01G000001.1;chr01;88551;93728;-
						
						if($halias->{$hmgene} && $halias->{$hmgene} ne '-'){
							$ginfo = $halias->{$hmgene}.";".$hgffdb->{$hmgene}{Chr}.";".$hgffdb->{$hmgene}{pos0}.";".$hgffdb->{$hmgene}{pos1}.";".$hgffdb->{$hmgene}{strand};
							push(@rtmpa, $ginfo);
						}
					}
					elsif($hgffq->{$hmgene}{Chr}){
						my $ginfo = $hmgene.";".$hgffq->{$hmgene}{Chr}.";".$hgffq->{$hmgene}{pos0}.";".$hgffq->{$hmgene}{pos1}.";".$hgffq->{$hmgene}{strand};
						push(@rtmp, $ginfo);			# Glyma.01G000174.v4.JE1;chr01;92087;92955;-,GlymaJER.01G000001.1;chr01;88551;93728;-
						push(@rtmpa, $ginfo);
					}
				}
				if(@rtmp){
					$B[9] = join(",", @rtmp);
				}
				if(@rtmpa){
					$Ba[9] = join(",", @rtmpa);
				}
			}
			if($B[10] && $B[10] ne '-'){
				my @Hmtmp = split(/,/, $B[10]);
				my @rtmp;
				my @rtmpa;
				foreach my $hmgene (@Hmtmp){
					if($hmgene =~ /\;/){
						my @subhmtmp = split(/\;/, $hmgene);
						$hmgene = $subhmtmp[0];
					}
					if($hgffdb->{$hmgene}{Chr}){
						my $ginfo = $hmgene.";".$hgffdb->{$hmgene}{Chr}.";".$hgffdb->{$hmgene}{pos0}.";".$hgffdb->{$hmgene}{pos1}.";".$hgffdb->{$hmgene}{strand};
						push(@rtmp, $ginfo);			# Glyma.01G000174.v4.JE1;chr01;92087;92955;-,GlymaJER.01G000001.1;chr01;88551;93728;-
						
						if($halias->{$hmgene} && $halias->{$hmgene} ne '-'){
							$ginfo = $halias->{$hmgene}.";".$hgffdb->{$hmgene}{Chr}.";".$hgffdb->{$hmgene}{pos0}.";".$hgffdb->{$hmgene}{pos1}.";".$hgffdb->{$hmgene}{strand};
							push(@rtmpa, $ginfo);
						}
					}
					elsif($hgffq->{$hmgene}{Chr}){
						my $ginfo = $hmgene.";".$hgffq->{$hmgene}{Chr}.";".$hgffq->{$hmgene}{pos0}.";".$hgffq->{$hmgene}{pos1}.";".$hgffq->{$hmgene}{strand};
						push(@rtmp, $ginfo);			# Glyma.01G000174.v4.JE1;chr01;92087;92955;-,GlymaJER.01G000001.1;chr01;88551;93728;-
						push(@rtmpa, $ginfo);
					}
				}
				if(@rtmp){
					$B[10] = join(",", @rtmp);
				}
				if(@rtmpa){
					$Ba[10] = join(",", @rtmpa);
				}
			}
			
			$rev .= join("\t", @B).$add_info."\n";
			
			if($halias->{$Ba[0]} && $halias->{$Ba[0]} ne '-'){
				$Ba[0] = $halias->{$Ba[0]};
				
				if($halias->{$Ba[14]} && $halias->{$Ba[14]} ne '-'){
					$Ba[14] = $halias->{$Ba[14]};
				}
				if($halias->{$Ba[15]} && $halias->{$Ba[15]} ne '-'){
					$Ba[15] = $halias->{$Ba[15]};
				}
				
				$rev_alias .= join("\t", @Ba).$add_info."\n";
				$cnt_alias++;
			}
		}
	}
	
	my $ncnt_ortho = $cnt_ortho_confirmed + $cnt_replace_ortho1 + $cnt_replace_ortho2 + $cnt_homolog_confirmed + $cnt_replace_homolog1 + $cnt_replace_homolog2;
	print "! [$cnt_overlap] overlapping assigned partner\n";
	if($cnt_seqalias > 0){
		print "! [$cnt_contrachr] possibly contradicting seq ID (translocated)\n";
	}
	print "! [$ncnt_ortho] : partner candidates in total\n";
	print " - [$cnt_ortho_confirmed] : 1-to-1 ortholog (confirmed)\n";
	print " - [$cnt_replace_ortho1] : assigned unlinked 1-to-1 ortholog\n";
	print " - [$cnt_replace_ortho2] : replaced 1-to-1 ortholog\n";
	print " - [$cnt_homolog_confirmed] : homolog -> 1-to-1 ortholog confirmed\n";
	print " - [$cnt_replace_homolog1] : homolog -> assigned unlinked 1-to-1 ortholog\n";
	print " - [$cnt_replace_homolog2] : homolog -> replaced 1-to-1 ortholog\n";
	print "! [$cnt_false_ortho] : not partner ortholog in strict (revised)\n";
	print "! [$cnt_false_homolog] : not partner homolog in strict (revised)\n";
	print "! [$cnt_NA_keep] : other (kept)\n";
}

return ($rev, $rev_alias, $cnt_alias, $rev_contradict);
}


#-------------------------------------------------------------------------------
sub Fitinfo2hash{
my $file = shift;

open(my $fh, "<", $file) or die;
my $cnt = 0;
my $htmp1 = {};
my $htmp2 = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($cnt == 0){
		$cnt++;
		next;
	}
	
	my @A = split(/\t/, $line);
	if($A[1] ne '-'){
		$htmp1->{$A[0]}{$A[1]}{fit} = $A[2];
		$htmp1->{$A[0]}{$A[1]}{nscore} = $A[3];
		$htmp1->{$A[0]}{$A[1]}{pscore} = $A[4];
		$htmp1->{$A[0]}{$A[1]}{judge} = $A[5];		# Asm2sv judge
		
		$htmp2->{$A[1]}{$A[0]}{fit} = $A[2];
		$htmp2->{$A[1]}{$A[0]}{nscore} = $A[3];
		$htmp2->{$A[1]}{$A[0]}{pscore} = $A[4];
		$htmp2->{$A[1]}{$A[0]}{judge} = $A[5];		# Asm2sv judge
	}
	$cnt++;
}
close $fh;

$cnt = 0;
my @GID1 = keys(%{$htmp1});
my @GID2 = keys(%{$htmp2});
my $hfit1 = {};
my $hfit2 = {};
my $include_partly_hit = 0;
my $record = "";

foreach my $gid1 (@GID1){
	my $heach = $htmp1->{$gid1};
	my @HID = keys(%{$heach});
	my $tag_present = 0;
	my $tag_nscore = 0;
	my $tag_pscore = 0;
	foreach my $hid (@HID){
		if(! defined $heach->{$hid}{fit} || $heach->{$hid}{fit} eq '-'){
			$heach->{$hid}{fit} = "-";
		}
		if(! defined $heach->{$hid}{nscore} || $heach->{$hid}{nscore} eq '-'){
			$heach->{$hid}{nscore} = 0;
		}
		if(! defined $heach->{$hid}{pscore} || $heach->{$hid}{pscore} eq '-'){
			$heach->{$hid}{pscore} = 0;
		}
		if($heach->{$hid}{judge} && $heach->{$hid}{judge} =~ /present/){		# Asm2sv judge
			$tag_present++;
		}
		if($heach->{$hid}{nscore} && $heach->{$hid}{nscore} > 0){
			$tag_nscore++;
		}
		if($heach->{$hid}{pscore} && $heach->{$hid}{pscore} > 0){
			$tag_pscore++;
		}
	}
	
	my $AoH = [];
	foreach my $hid (@HID){
		if($heach->{$hid}{fit} eq '-'){
			next;
		}
		
		my @tmp;
		my $pusharray = 0;
		if($tag_present > 0 && $tag_pscore > 0 && $tag_nscore > 0){
			if($heach->{$hid}{judge} && $heach->{$hid}{judge} =~ /present/ && $heach->{$hid}{pscore} && $heach->{$hid}{pscore} > 0 && $heach->{$hid}{nscore} && $heach->{$hid}{nscore} > 0){
				my $calc_score = $heach->{$hid}{nscore} + $heach->{$hid}{pscore} - $heach->{$hid}{fit};
				push(@tmp, $hid);
				push(@tmp, $heach->{$hid}{fit});
				push(@tmp, $heach->{$hid}{nscore});
				push(@tmp, $heach->{$hid}{pscore});
				push(@tmp, $heach->{$hid}{judge});
				push(@tmp, $calc_score);
				$pusharray = 1;
			}
		}
		elsif($tag_present > 0 && $tag_pscore == 0 && $tag_nscore > 0 && $include_partly_hit eq '1'){
			if($heach->{$hid}{judge} && $heach->{$hid}{judge} =~ /present/ && $heach->{$hid}{nscore} && $heach->{$hid}{nscore} > 0){
				my $calc_score = $heach->{$hid}{nscore} - $heach->{$hid}{fit};
				push(@tmp, $hid);
				push(@tmp, $heach->{$hid}{fit});
				push(@tmp, $heach->{$hid}{nscore});
				push(@tmp, $heach->{$hid}{pscore});
				push(@tmp, $heach->{$hid}{judge});
				push(@tmp, $calc_score);
				$pusharray = 1;
			}
		}
		elsif($tag_present > 0 && $tag_pscore > 0 && $tag_nscore == 0 && $include_partly_hit eq '1'){
			if($heach->{$hid}{judge} && $heach->{$hid}{judge} =~ /present/ && $heach->{$hid}{pscore} && $heach->{$hid}{pscore} > 0){
				my $calc_score = $heach->{$hid}{pscore} - $heach->{$hid}{fit};
				push(@tmp, $hid);
				push(@tmp, $heach->{$hid}{fit});
				push(@tmp, $heach->{$hid}{nscore});
				push(@tmp, $heach->{$hid}{pscore});
				push(@tmp, $heach->{$hid}{judge});
				push(@tmp, $calc_score);
			}
		}
		elsif($tag_present > 0 && $include_partly_hit eq '1'){
			if($heach->{$hid}{judge} && $heach->{$hid}{judge} =~ /present/){
				my $calc_score = 0 - $heach->{$hid}{fit};
				push(@tmp, $hid);
				push(@tmp, $heach->{$hid}{fit});
				push(@tmp, $heach->{$hid}{nscore});
				push(@tmp, $heach->{$hid}{pscore});
				push(@tmp, $heach->{$hid}{judge});
				push(@tmp, $calc_score);
			}
		}
		if($pusharray eq '1'){
			push(@{$AoH}, \@tmp);
		}
	}
	if(@{$AoH}){
#		@{$AoH} = sort {$b->[2] <=> $a->[2]} @{$AoH};
#		@{$AoH} = sort {$b->[3] <=> $a->[3]} @{$AoH};
#		@{$AoH} = sort {$a->[1] <=> $b->[1]} @{$AoH};
		@{$AoH} = sort {$b->[5] <=> $a->[5]} @{$AoH};
		
		$hfit1->{$gid1} = $AoH->[0][0];
		
		foreach my $tmp (@{$AoH}){
			$record .= $gid1."\t".join("\t", @{$tmp})."\n";
		}
	}
}

foreach my $gid2 (@GID2){
	my $heach = $htmp2->{$gid2};
	my @HID = keys(%{$heach});
	my $tag_present = 0;
	my $tag_nscore = 0;
	my $tag_pscore = 0;
	foreach my $hid (@HID){
		if($heach->{$hid}{fit} eq '-'){
			next;
		}
		
		if(! defined $heach->{$hid}{fit} || $heach->{$hid}{fit} eq '-'){
			$heach->{$hid}{fit} = "-";
		}
		if(! defined $heach->{$hid}{nscore} || $heach->{$hid}{nscore} eq '-'){
			$heach->{$hid}{nscore} = 0;
		}
		if(! defined $heach->{$hid}{pscore} || $heach->{$hid}{pscore} eq '-'){
			$heach->{$hid}{pscore} = 0;
		}
		if($heach->{$hid}{judge} && $heach->{$hid}{judge} =~ /present/){
			$tag_present++;
		}
		if($heach->{$hid}{nscore} && $heach->{$hid}{nscore} > 0){
			$tag_nscore++;
		}
		if($heach->{$hid}{pscore} && $heach->{$hid}{pscore} > 0){
			$tag_pscore++;
		}
	}
	
	my $AoH = [];
	foreach my $hid (@HID){
		if($heach->{$hid}{fit} eq '-'){
			next;
		}
		
		my @tmp;
		my $pusharray = 0;
		if($tag_present > 0 && $tag_pscore > 0 && $tag_nscore > 0){
			if($heach->{$hid}{judge} && $heach->{$hid}{judge} =~ /present/ && $heach->{$hid}{pscore} && $heach->{$hid}{pscore} > 0 && $heach->{$hid}{nscore} && $heach->{$hid}{nscore} > 0){
				my $calc_score = $heach->{$hid}{nscore} + $heach->{$hid}{pscore} - $heach->{$hid}{fit};
				push(@tmp, $hid);
				push(@tmp, $heach->{$hid}{fit});
				push(@tmp, $heach->{$hid}{nscore});
				push(@tmp, $heach->{$hid}{pscore});
				push(@tmp, $heach->{$hid}{judge});
				push(@tmp, $calc_score);
				$pusharray = 1;
			}
		}
		elsif($tag_present > 0 && $tag_pscore == 0 && $tag_nscore > 0 && $include_partly_hit eq '1'){
			if($heach->{$hid}{judge} && $heach->{$hid}{judge} =~ /present/ && $heach->{$hid}{nscore} && $heach->{$hid}{nscore} > 0){
				my $calc_score = $heach->{$hid}{nscore} - $heach->{$hid}{fit};
				push(@tmp, $hid);
				push(@tmp, $heach->{$hid}{fit});
				push(@tmp, $heach->{$hid}{nscore});
				push(@tmp, $heach->{$hid}{pscore});
				push(@tmp, $heach->{$hid}{judge});
				push(@tmp, $calc_score);
				$pusharray = 1;
			}
		}
		elsif($tag_present > 0 && $tag_pscore > 0 && $tag_nscore == 0 && $include_partly_hit eq '1'){
			if($heach->{$hid}{judge} && $heach->{$hid}{judge} =~ /present/ && $heach->{$hid}{pscore} && $heach->{$hid}{pscore} > 0){
				my $calc_score = $heach->{$hid}{pscore} - $heach->{$hid}{fit};
				push(@tmp, $hid);
				push(@tmp, $heach->{$hid}{fit});
				push(@tmp, $heach->{$hid}{nscore});
				push(@tmp, $heach->{$hid}{pscore});
				push(@tmp, $heach->{$hid}{judge});
				push(@tmp, $calc_score);
			}
		}
		elsif($tag_present > 0 && $include_partly_hit eq '1'){
			if($heach->{$hid}{judge} && $heach->{$hid}{judge} =~ /present/){
				my $calc_score = 0 - $heach->{$hid}{fit};
				push(@tmp, $hid);
				push(@tmp, $heach->{$hid}{fit});
				push(@tmp, $heach->{$hid}{nscore});
				push(@tmp, $heach->{$hid}{pscore});
				push(@tmp, $heach->{$hid}{judge});
				push(@tmp, $calc_score);
			}
		}
		if($pusharray eq '1'){
			push(@{$AoH}, \@tmp);
		}
	}
	if(@{$AoH}){
#		@{$AoH} = sort {$b->[2] <=> $a->[2]} @{$AoH};
#		@{$AoH} = sort {$b->[3] <=> $a->[3]} @{$AoH};
#		@{$AoH} = sort {$a->[1] <=> $b->[1]} @{$AoH};
		@{$AoH} = sort {$b->[5] <=> $a->[5]} @{$AoH};
		
		$hfit2->{$gid2} = $AoH->[0][0];
		
		foreach my $tmp (@{$AoH}){
			$record .= $gid2."\t".join("\t", @{$tmp})."\n";
		}
	}
}

return ($hfit1, $hfit2, $record);
}


#-------------------------------------------------------------------------------
sub FitPos0{
my ($DBG, $p0, $p1, $rfile_bdhn, $rfile_bdhp, $asm2sv_tsv, $hgffdb, $hgffq, $hseqinfo, $dbpref, $qpref, $th_flanking, $increment, $strict_chr_identical, $i) = @_;

my $hASV = {};
if(-e $asm2sv_tsv){
	open(my $fh, "<", $asm2sv_tsv) or die;
	my $cnt = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		unless($line){
			next;
		}
		
		my @A = split(/\t/, $line);
		my $dbgene = $A[23];
		if(! $hASV->{$dbgene}){
			$hASV->{$dbgene}{dbseqid} = $A[2];
			$hASV->{$dbgene}{qseqid} = $A[8];
			$hASV->{$dbgene}{fpos0} = $A[9];
			$hASV->{$dbgene}{fpos1} = $A[10];
			$hASV->{$dbgene}{gpos0} = $A[12];
			$hASV->{$dbgene}{gpos1} = $A[13];
			$hASV->{$dbgene}{judge} = $A[22];
		}
	}
	close $fh;
}

#	0	ID 1=[Gmax_Enrei_pseudomol_v3.31]
#	1	gene_id [1]
#	2	Asm2sv_partner [1]
#	3	ID 2=[Gmax_Fukuyutaka_pseudomol_v3.12_nuclear_transcript]
#	4	gene_id [2]
#	5	Asm2sv_BDH_parter
#	6	[1] length
#	7	[2] length
#	8	[1] alignment length
#	9	[2] alignment length
#	10	[1] %identity
#	11	[2] %identity
#	12	[1] %coverage
#	13	[2] %coverage
#	14	[1] score
#	15	[2] score
#	16	[1] Eval
#	17	[2] Eval
#	18	k
#	19	l

my $hBDH = {};
if(-e $rfile_bdhn){
	open(my $fh, "<", $rfile_bdhn) or die;
	my $mode = "nuc";
	my $cnt = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		unless($line){
			next;
		}
		if($line =~ /\[1/ && $line =~ /\[2/){
			next;
		}
		
		my @A = split(/,/, $line);
		if($A[14] && $A[15] && $A[14] =~ /\d/ && $A[15] =~ /\d/){
			if(! $hBDH->{$mode}{$A[1]}{$A[4]}){
				$hBDH->{$mode}{$A[1]}{$A[4]} = $A[14] + $A[15];
				$cnt++;
			}
			else{
				my $sum_score = $A[14] + $A[15];
				if($hBDH->{$mode}{$A[1]}{$A[4]} < $sum_score){
					$hBDH->{$mode}{$A[1]}{$A[4]} = $sum_score;
				}
			}
		}
	}
	close $fh;
}
if(-e $rfile_bdhp){
	open(my $fh, "<", $rfile_bdhp) or die;
	my $mode = "prot";
	my $cnt = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		unless($line){
			next;
		}
		if($line =~ /\[1/ && $line =~ /\[2/){
			next;
		}
		
		my @A = split(/,/, $line);
		if($A[14] && $A[15] && $A[14] =~ /\d/ && $A[15] =~ /\d/){
			if(! $hBDH->{$mode}{$A[1]}{$A[4]}){
				$hBDH->{$mode}{$A[1]}{$A[4]} = $A[14] + $A[15];
				$cnt++;
			}
			else{
				my $sum_score = $A[14] + $A[15];
				if($hBDH->{$mode}{$A[1]}{$A[4]} < $sum_score){
					$hBDH->{$mode}{$A[1]}{$A[4]} = $sum_score;
				}
			}
		}
	}
	close $fh;
}

my $r = "";
my $cnt_pair = 0;
my $cnt_skip_idalias = 0;
for(my $j = $p0; $j < $p1; $j++){
	if($DBG->[$j]){
		my $dbgene = $DBG->[$j];
		my $dbsid = $hgffdb->{$dbgene}{Chr};
		my $hASV_dbgene = {};
		my $NrQ = [];
		if($hBDH->{nuc}{$dbgene} && $hBDH->{prot}{$dbgene} && $hASV->{$dbgene}){
			my $hash1 = $hBDH->{nuc}{$dbgene};
			my $hash2 = $hBDH->{prot}{$dbgene};
			$hASV_dbgene = $hASV->{$dbgene};
			
			my @A1 = keys(%{$hash1});
			@A1 = sort {$a cmp $b} @A1;
			my @A2 = keys(%{$hash2});
			@A2 = sort {$a cmp $b} @A2;
			
			my $hash = {};
			foreach my $str (@A1){
				$hash->{$str} = 1;
			}
			foreach my $str (@A2){
				$hash->{$str} = 1;
			}
			
			@{$NrQ} = keys(%{$hash});
			@{$NrQ} = sort {$a cmp $b} @{$NrQ};
		}
		
		if(@{$NrQ}){
			my $AoNrQ = [];
			foreach my $qcand (@{$NrQ}){
				if($hgffq->{$qcand} && defined $hASV->{$dbgene}{qseqid} && $hASV->{$dbgene}{qseqid} eq $hgffq->{$qcand}{Chr}){
					my $qsid = $hgffq->{$qcand}{Chr};
					if($hseqinfo->{$dbpref}{$dbsid} && $hseqinfo->{$qpref}{$qsid} && $hgffdb->{$dbgene}{scaffold} eq 'false' && $hgffq->{$qcand}{scaffold} eq 'false'){
						if($strict_chr_identical eq 'true' && $hseqinfo->{$dbpref}{$dbsid} ne $hseqinfo->{$qpref}{$qsid}){
							$cnt_skip_idalias++;
							next;
						}
					}
					
					my $qpos0_gff = $hgffq->{$qcand}{pos0};
					my $qpos1_gff = $hgffq->{$qcand}{pos1};
					
					for(my $t = 0; $t < $th_flanking; $t += $increment){
						if($hASV->{$dbgene}{gpos0} ne '-' && $hASV->{$dbgene}{gpos1} ne '-' && $hASV->{$dbgene}{gpos0} =~ /\d/ && $hASV->{$dbgene}{gpos1} =~ /\d/){
							my $qposg0_ASV = $hASV->{$dbgene}{gpos0} - $t;
							my $qposg1_ASV = $hASV->{$dbgene}{gpos1} + $t;
							if($qposg0_ASV <= $qpos0_gff && $qpos1_gff <= $qposg1_ASV){
								my @tmp;
								push(@tmp, $qcand);
								push(@tmp, $t);
								push(@{$AoNrQ}, \@tmp);
								last;
							}
						}
					}
				}
			}
			my $each_r;
			if(@{$AoNrQ}){
				@{$AoNrQ} = sort {$a->[1] <=> $b->[1]} @{$AoNrQ};
				foreach my $tmp (@{$AoNrQ}){
					unless($hBDH->{nuc}{$dbgene}{$tmp->[0]}){
						$hBDH->{nuc}{$dbgene}{$tmp->[0]} = "-";
					}
					unless($hBDH->{prot}{$dbgene}{$tmp->[0]}){
						$hBDH->{prot}{$dbgene}{$tmp->[0]} = "-";
					}
					
					$each_r .= "$dbgene\t$tmp->[0]\t$tmp->[1]\t$hBDH->{nuc}{$dbgene}{$tmp->[0]}\t$hBDH->{prot}{$dbgene}{$tmp->[0]}\t$hASV->{$dbgene}{judge}\n";
					$cnt_pair++;
				}
			}
			if(! $each_r){
				$each_r = "$dbgene\t-\t-\t-\t-\t$hASV->{$dbgene}{judge}\n";
			}
			$r .= $each_r;
		}
	}
}

if($strict_chr_identical eq 'true' && $cnt_skip_idalias > 0){
	print " - thread [$i] | [$cnt_pair] pairs, [$cnt_skip_idalias] unmatching seq ID alias\n";
}
else{
	print " - thread [$i] | [$cnt_pair] pairs\n";
}

return $r;
}


#-------------------------------------------------------------------------------
sub Read_BDconfident{
my ($hash, $file, $mode) = @_;

#print " - reading [$file]...\n";
open(my $fh, "<", $file) or die;
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	unless($line){
		next;
	}
	
	#	0	ID 1=[Gmax_Enrei_pseudomol_v3.31]
	#	1	gene_id [1]
	#	2	Asm2sv_partner [1]
	#	3	ID 2=[Gmax_Fukuyutaka_pseudomol_v3.12_nuclear_transcript]
	#	4	gene_id [2]
	#	5	Asm2sv_BDH_parter
	#	6	[1] length
	#	7	[2] length
	#	8	[1] alignment length
	#	9	[2] alignment length
	#	10	[1] %identity
	#	11	[2] %identity
	#	12	[1] %coverage
	#	13	[2] %coverage
	#	14	[1] score
	#	15	[2] score
	#	16	[1] Eval
	#	17	[2] Eval
	#	18	k
	#	19	l
	
	my @A = split(/,/, $line);
	unless($hash->{$mode}{$A[1]}{$A[4]}){
		$hash->{$mode}{db2q}{$A[1]}{$A[4]}{identity} = $A[10];
		$hash->{$mode}{db2q}{$A[1]}{$A[4]}{pcov} = $A[12];
		$hash->{$mode}{db2q}{$A[1]}{$A[4]}{score} = $A[14];
		$hash->{$mode}{db2q}{$A[1]}{$A[4]}{Eval} = $A[16];
		
		$hash->{$mode}{q2db}{$A[4]}{$A[1]}{identity} = $A[11];
		$hash->{$mode}{q2db}{$A[4]}{$A[1]}{pcov} = $A[13];
		$hash->{$mode}{q2db}{$A[4]}{$A[1]}{score} = $A[15];
		$hash->{$mode}{q2db}{$A[4]}{$A[1]}{Eval} = $A[17];
		
		if(defined $A[18] && defined $A[19] && $A[18] eq '0' && $A[19] eq '0'){
			$hash->{$mode}{db2q}{$A[1]}{$A[4]}{status} = "BDBH";
			$hash->{$mode}{q2db}{$A[4]}{$A[1]}{status} = "BDBH";
		}
		else{
			$hash->{$mode}{db2q}{$A[1]}{$A[4]}{status} = "BDH";
			$hash->{$mode}{q2db}{$A[4]}{$A[1]}{status} = "BDH";
		}
		$cnt++;
	}
}
close $fh;

print " - [$cnt] candidates in [$file]\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Gff_to_hash{
my $gff3 = shift;

#print "! reading [$gff3]...\n";
open(my $fh, "<", $gff3) or die;
my $hash = {};
my @GID;
my $cnt = 0;
my $gcnt = 0;
my $tcnt = 0;
my $pcnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	$line =~ s/\"//g;
	
	if($line =~ /\#/){
		next;
	}
	
	my @A = split(/\t/, $line);
	unless($A[8]){
		next;
	}
	my @tag = split(/\;/, $A[8]);
	my $gid;
	my $tid;
	my $name_evm;
	$cnt++;
	
	if($A[2] eq 'gene'){
		foreach my $val (@tag){
			if($val =~ /ID\=/ && $val !~ /Alt_ID\=/){
				$gid = $val;
				$gid =~ s/ID\=//;
			}
			elsif($val =~ /Name\=EVM/){
				$name_evm = $val;
				$name_evm =~ s/Name\=EVM//;
			}
		}
		
#		if($name_evm){
#			$A[8] = "ID\=".$gid.";Note\=EVM";
#		}
#		else{
			$A[8] = "ID\=".$gid;
#		}
		
		if($gid && ! $hash->{hgff}{$gid}){
			$hash->{hgff}{$gid}{Chr} = $A[0];		#1	seqid
			$hash->{hgff}{$gid}{pos0} = $A[3];		#2	pos0
			$hash->{hgff}{$gid}{pos1} = $A[4];		#3	pos1
			$hash->{hgff}{$gid}{strand} = $A[6];	#4	strand
			push(@GID, $gid);
			
			if($A[0] =~ /scaffold/ || $A[0] =~ /unanchored/ || $A[0] =~ /mitochon/ || $A[0] =~ /chloro/){
				$hash->{hgff}{$gid}{scaffold} = "true";
			}
			else{
				$hash->{hgff}{$gid}{scaffold} = "false";
			}
			$gcnt++;
		}
	}
	elsif($A[2] eq 'mRNA' || $A[2] eq 'transcript'){
		foreach my $val (@tag){
			if($val =~ /ID\=/ && $val !~ /Alt_ID\=/){
				$tid = $val;
				$tid =~ s/ID\=//;
			}
			elsif($val =~ /Parent\=/){
				$gid = $val;
				$gid =~ s/Parent\=//;
			}
		}
		
		if($gid && $tid){
			$hash->{t2g}{$tid} = $gid;
			$hash->{g2t}{$gid} .= $tid."\n";
			$tcnt++;
		}
	}
	elsif($A[2] eq 'cds' || $A[2] eq 'CDS'){
		foreach my $val (@tag){
			if($val =~ /Parent\=/){
				$tid = $val;
				$tid =~ s/Parent\=//;
			}
		}
		
		if($tid){
			$hash->{CDS}{$tid} = "true";
			$pcnt++;
		}
	}
	elsif($A[2] eq 'exon'){
		foreach my $val (@tag){
			if($val =~ /Parent\=/){
				$tid = $val;
				$tid =~ s/Parent\=//;
			}
		}
		
		if($tid){
			$hash->{CDS}{$tid} = "true";
			$pcnt++;
		}
	}
}
close $fh;

$hash->{GID} = \@GID;
print " - [$gcnt] genes, [$tcnt] transcripts in [$gff3]\n";

return $hash;
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

#-------------------------------------------------------------------------------
sub Gffread2seq{
my ($bin_gffread, $genome, $gff, $rfile_transcript, $rfile_cds, $rfile_protein) = @_;

my $cmd1 = "$bin_gffread $gff -g $genome -w $rfile_transcript > /dev/null 2>&1";
my $cmd2 = "$bin_gffread $gff -g $genome -x $rfile_cds > /dev/null 2>&1";
my $cmd3 = "$bin_gffread $gff -g $genome -y $rfile_protein > /dev/null 2>&1";

my $err = 0;
if(! -e $rfile_transcript){
	if(system($cmd1) != 0){
		print " - [$rfile_transcript] failed.\n";
		$err++;
	}
	else{
		print " - [$rfile_transcript]\n";
	}
}
else{
	print " - [$rfile_transcript] already exists.\n";
}
if(! -e $rfile_cds){
	if(system($cmd2) != 0){
		print " - [$rfile_cds] failed.\n";
		$err++;
	}
	else{
		print " - [$rfile_cds]\n";
	}
}
else{
	print " - [$rfile_cds] already exists.\n";
}
if(! -e $rfile_protein){
	if(system($cmd3) != 0){
		print " - [$rfile_protein] failed.\n";
		$err++;
	}
	else{
		print " - [$rfile_protein]\n";
	}
}
else{
	print " - [$rfile_protein] already exists.\n";
}

return $err;
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
sub TestBlat{
my $bin = shift;

my $tmp = "_test.txt";
my $err = 1;
system("$bin > $tmp 2>&1");
if(-e $tmp){
	open(my $fh, "<", $tmp);
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line){
			if($line =~ /Standalone BLAT/){
				$err = 0;
				last;
			}
		}
	}
	close $fh;
	system("rm $tmp");
}
if($err == 1){
	print "! [$bin] may not be exectable... (please check it)\n";
}
else{
	print "! [$bin] is exectable\n";
}

return $err;
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

