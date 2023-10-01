#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use Getopt::Long;
use FindBin;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "Psl_tigPlot3.pl version 1.42\n";
$version .= "last update: [2023\/1\/28]\n";
$version .= "copyright: ryoichi yano\n";

#print $version;

#-------------------------------------------------------------------------------

#---------------------------------------------------------//
my $help_command =<<"EOS";
Basic usage: 
  /path/to/Psl_tigPlot3.pl -d [db fasta] -q [query fasta] -s [bp scale] -o [Result_directory] -I [db seqID] -r [XX-YY] -t [AA-BB] --cpu [CPU thread]

--db or -d           : reference sequence (fasta, required)
--query or -q        : query sequence (fasta, required)
--wgpsldir or -w     : pre-analyzed whole genome alignment dir (optional)
--scale or -s        : scale (Mb, kb, bp)
--ID or -I           : sequence ID in reference fasta
--region or -r       : target sequence region in reference fasta
--target or -t       : specified region to be drawn in R-plot and PAV analysis
--match or -m        : percent similarity threshold for blast alignment (1-100)
--plot or -p         : output R-plot
--all or -a          : output R-plots for all candidate alignments
--keep or -k         : keep intermediate files
--cpu or -c          : CPU thread number (default 4)
--note or -n         : free description included in summary together with result
--output_dir or -o   : output directory
--help or -h         : display usage

EOS

#---------------------------------------------------------//
#gettin parameters from command line
my $dir;
my $dbfasta;
my $qfasta;
my $wgpsldir;
my $scale;
my $seqid;
my $tregion;
my $ablines;
my $tstrand;
my $neighbor_tlen;
my $mth;
my $allplot;
my $outplot;
my $skip_alignmethod;
my $keep;
my $cpu;
my $note;
my $help;

GetOptions('--db=s' => \$dbfasta, '--query=s' => \$qfasta, '--wgpsldir=s' => \$wgpsldir, '--scale=s' => \$scale, '--match=i' => \$mth, '--neighbor=i' => \$neighbor_tlen, '--ID=s' => \$seqid,  '--region=s' => \$tregion, '--target=s' => \$ablines, '--output_dir=s' => \$dir, '--cpu=i' => \$cpu, '--all' => \$allplot, '--skip_alignmethod' => \$skip_alignmethod, '--keep' => \$keep, '--plot' => \$outplot, '--note=s' => \$note, '--strand=s' => \$tstrand, '--help' => \$help);

my $lth = 1000;
my $tpos0;
my $tpos1;
my $nomask = "false";
#my $nomask = "true";

unless($dbfasta){
	print "! missing db fasta input\n";
	print "\n$help_command\n";
	goto END;
}
if($dbfasta && ! -e $dbfasta){
	print "! missing db fasta  input\n";
	print "\n$help_command\n";
	goto END;
}

unless($qfasta){
	print "! missing query fasta input\n";
	print "\n$help_command\n";
	goto END;
}
if($qfasta && ! -e $qfasta){
	print "! missing query fasta [$qfasta]\n";
	print "\n$help_command\n";
	goto END;
}

#if(! $wgpsldir || ! -d $wgpsldir){
	$wgpsldir = "null";
	$skip_alignmethod = "false";
#}
if(! $skip_alignmethod){
	$skip_alignmethod = "false";
}
else{
	$skip_alignmethod = "true";
}

unless($lth){
	print "! missing seq length threshold (bp)\n";
	print "\n$help_command\n";
	goto END;
}
if($lth =~ /[a-z]/i){
	print "! missing seq length threshold (bp)\n";
	print "\n$help_command\n";
	goto END;
}

unless($scale){
	$scale = "bp";
}
elsif($scale =~ /mb/i){
	$scale = "Mb";
}
elsif($scale =~ /kb/i){
	$scale = "kb";
}
elsif($scale =~ /b/i || $scale =~ /bp/i){
	$scale = "bp";
}
else{
	print "! missing scale info (Mb, kb, bp)\n";
	print "\n$help_command\n";
	goto END;
}

unless($seqid){
	print "! missing seq ID (chr ID)\n";
	print "\n$help_command\n";
	goto END;
}

unless($tregion){
	$nomask = "true";
	$tpos0 = "NA";
	$tpos1 = "NA";
#	print "! missing target genomic region [2] : specify region like \"10-20\"\n";
#	print "\n$help_command\n";
#	goto END;
}
else{
	my @P = split(/-/, $tregion);
	$tpos0 = $P[0];
	$tpos1 = $P[1];
	
	unless($tpos0){
		print "! missing tpos0\n";
		goto END;
	}
	elsif($tpos0 =~ /[a-z]/i){
		print "! abnormal input [$tpos0]\n";
		goto END;
	}
	
	unless($tpos1){
		print "! missing tpos1\n";
		goto END;
	}
	elsif($tpos1 =~ /[a-z]/i){
		print "! abnormal input [$tpos1]\n";
		goto END;
	}
}

my $ab0 = -1;
my $ab1 = -1;
if($ablines){
	my @P = split(/-/, $ablines);
	$ab0 = $P[0];
	$ab1 = $P[1];
	
	if($ab0 =~ /[a-z]/i){
		print "! abnormal input [$ab0]\n";
		goto END;
	}
	
	if($ab1 =~ /[a-z]/i){
		print "! abnormal input [$ab1]\n";
		goto END;
	}
}

unless($cpu){
	$cpu = 2;
}
$cpu = 2;

unless($mth){
	$mth = 95;
}
unless($allplot){
	$allplot = 0;
}
else{
	$allplot = 1;
}
unless($outplot){
	$outplot = 0;
	$allplot = 0;
}
elsif($outplot){
	$outplot = 1;
}
unless($keep){
	$keep = 0;
}
else{
	$keep = 1;
}
unless($note){
	$note = "-";
}

#$keep = 1;
#$outplot = 1;

if(! $tstrand){
	$tstrand = "null";
}
elsif($tstrand eq 'plus' || $tstrand eq '+'){
	$tstrand = "+";
}
elsif($tstrand eq 'minus' || $tstrand eq '-'){
	$tstrand = "-";
}
else{
	$tstrand = "null";
}

if($tstrand eq 'null'){
	if($neighbor_tlen){
		print "\n! disabled --neighbor as --strand is not specified\n";
	}
	$neighbor_tlen = 0;
}
if(! $neighbor_tlen || $neighbor_tlen =~ /\D/){
	if($tstrand ne 'null'){
		print "\n! disabled --strand as --neighbor is not specified\n";
	}
	$tstrand = "null";
}
elsif($neighbor_tlen > 10000){
	$neighbor_tlen = 10000;
}

my $flen = 0;

my $log = "\n---------------------------------------------------------------------------------------\n";
$log .= "! db fasta               : [$dbfasta]\n";
$log .= "! target fasta           : [$qfasta]\n";
#$log .= "! pre-analysis psl dir   : [$wgpsldir]\n";
$log .= "! skip align method      : [$skip_alignmethod]\n";
$log .= "! block length threshold : [$lth]\n";
$log .= "! db fasta seqID         : [$seqid]\n";
$log .= "! db fasta region        : [$tpos0] - [$tpos1]\n";
if($ab0 > 0 && $ab1 > 0){
	$log .= "! specified region       : [$ab0] - [$ab1]\n";
	$flen = abs($ab0 - $tpos0);
	$log .= "! flanking seq span      : [$flen] bp\n";
}
if($tstrand ne 'null'){
	$log .= "! strand                 : [$tstrand]\n";
	$log .= "! promoter/UTR seq span  : [$neighbor_tlen] bp\n";
}
$log .= "! entry name             : [$note]\n";

print $log;

my $q = "y";
unless($q){
#	print $analysis_log;
	
	print "\nOK? (Y/n): ";
	$q = <STDIN>;
	$q =~ s/\n//;
	$q =~ s/\r//;
	unless($q){
		$q = "y";
	}
}

if($q =~ /y/i || $q =~ /pipe/i){
	if($q =~ /y/i){
		print "\n";
	}
	
	if($scale eq 'Mb'){
		$tpos0 *= 1000000;
		$tpos1 *= 1000000;
		$ab0 *= 1000000;
		$ab1 *= 1000000;
	}
	elsif($scale eq 'kb'){
		$tpos0 *= 1000;
		$tpos1 *= 1000;
		$ab0 *= 1000000;
		$ab1 *= 1000000;
	}
	
	my $dbpref = $dbfasta;
	if($dbpref =~ /_tmp_/){
		my @DBPref = split(/_tmp_/, $dbpref);
		$dbpref = $DBPref[1];
	}
	if($dbpref =~ /\.fasta/){
		$dbpref =~ s/\.fasta//g;
	}
	elsif($dbpref =~ /\.fa/){
		$dbpref =~ s/\.fa//g;
	}
	
	my $qpref = $qfasta;
	if($qpref =~ /\.fasta/){
		$qpref =~ s/\.fasta//g;
	}
	elsif($qpref =~ /\.fa/){
		$qpref =~ s/\.fa//g;
	}
	if($qpref =~ /_tmp_/){
		my @tmp = split(/_tmp_/, $qpref);
		$qpref = $tmp[1];
	}
	
	unless($dir){
		$dir = $dbpref."_".$seqid."_".$tpos0."-".$tpos1."_vs_".$qpref;
	}
	unless(-e $dir){
		system("mkdir $dir");
	}
	chdir $dir;
	
	unless(-e $dbfasta){
		system("ln -s ../$dbfasta ./");
	}
	
	unless(-e $qfasta){
		system("ln -s ../$qfasta ./");
	}
	
	my $hqtarget = Open_fasta_as_hash($qfasta);
	my $qpsl = "../".$wgpsldir."/hint.psl";
	my $apref = "q-".$dbpref."_".$seqid."_".$tpos0."-".$tpos1."_target-".$qpref;
	my $afile = "specified_region_hits_".$dir.".tsv";
	my $rdfile = "raw_data_specified_region_hits_".$dir.".tsv";
	my $bfile = "specified_region_hits_db-".$dbpref."_q-".$qpref.".tsv";
	my $abhit_header = "directory\tdbfasta\tseqid\tpos0\tpos1\tsp0\tsp1\tqfasta\tqfasta seqid\thit pos0\thit pos1\tgroup\thit sp0\thit sp1\tbp (sp0-sp1)\t";
	$abhit_header .= "bp hit-align (sp0-sp1)\tbp hit-span (sp0-sp1)\tbp hit-span woN (sp0-sp1)\tbp insertion\talign ratio (1=not disrupted)\tinsert ratio (1=not disrupted)\t";
	$abhit_header .= "seq normality (1=not disrupted)\tjudge\tnote\t";
	$abhit_header .= "5'-promoter bp hit ($neighbor_tlen bp)\t5'-promoter align ratio ($neighbor_tlen bp)\t3'-UTR bp hit ($neighbor_tlen bp)\t3'-UTR align ratio ($neighbor_tlen bp)\t";
	$abhit_header .= "5' flanking bp hit ($flen bp)\t5' flanking align ratio ($flen bp)\t3' flanking bp hit ($flen bp)\t3' flanking align ratio ($flen bp)\n";
	
	my $method_eachalign = "true";
	if($method_eachalign eq 'true'){
		#-----------------------------------------------------//
		my $seqfasta = $seqid."_".$tpos0."-".$tpos1.".fasta";
		my $blast_result = $seqid."_".$tpos0."-".$tpos1."_".$dbpref."_".$qpref.".blast";
		
		print "! preparing fasta...\n";
		my $dbseqlen = Search_seq($dbfasta, $seqid, $tpos0, $tpos1, $seqfasta, "NA", "NA", "false");
		unless(-e $seqfasta){
			goto END;
		}
		if($tpos0 eq 'NA' && $tpos1 eq 'NA' && $dbseqlen){
			$tpos0 = 1;
			$tpos1 = $dbseqlen;
		}
		
		unless(-e $blast_result){
			my $makedb_tmp = "_tmp.".$qfasta."_makedb.txt";
			my $cmd1 = "makeblastdb -in $qfasta -dbtype nucl -hash_index 1> $makedb_tmp";
			print "! making blast db for [$qfasta]...\n";
			system($cmd1);
			
			my $cmd2 = "blastn -db $qfasta -query $seqfasta -outfmt 6 -num_threads $cpu -max_target_seqs 3 > $blast_result";
			print "! cmd=[$cmd2]\n";
			system($cmd2);
			
			Delete("$qfasta.nhd");
			Delete("$qfasta.nhi");
			Delete("$qfasta.nhr");
			Delete("$qfasta.nin");
			Delete("$qfasta.nog");
			Delete("$qfasta.nsd");
			Delete("$qfasta.nsi");
			Delete("$qfasta.nsq");
			Delete("$qfasta.ndb");
			Delete("$qfasta.not");
			Delete("$qfasta.ntf");
			Delete("$qfasta.nto");
			Delete("$makedb_tmp");
		}
		else{
			print "! [$blast_result] already exists. skipping blast...\n";
		}
		
		my $wth = abs($tpos1 - $tpos0);
	#	my $lth2 = abs($tpos1 - $tpos0) / 10;
		my $lth2 = 1000;
		my $AoCand = Read_blast2($seqid, $blast_result, $mth, $wth, $lth2, $tpos0, $tpos1);
		unless($AoCand){
			print "! missing blast hit for [$seqid] in [$qfasta]...\n";
			goto END;
		}
		
		my $numQ = @{$AoCand};
		my $nohit = "false";
		if($numQ == 0){
			print "! no blast hit...\n";
			$nohit = "true";
			goto JUMP;
		}
		
		my $qtarget = $qpref."_with_q-".$dbpref."_".$seqid."_".$tpos0."-".$tpos1."_selected.fasta";
		if($keep && -e $qtarget){
			goto SKIP;
		}
		elsif(-e $qtarget){
			system("rm $qtarget");
		}
		
		print "! preparing fasta...\n";
	#	for(my $i = 0; $i < $numQ; $i++){
		for(my $i = 0; $i < 5; $i++){
			unless($AoCand->[$i]){
				last;
			}
			
			my $cand = $AoCand->[$i];
			my $qseqid = $cand->[0];
			my $hpos0 = $cand->[1];
			my $hpos1 = $cand->[2];
			
			if($nomask eq 'true'){
				Search_seq($qfasta, $qseqid, "NA", "NA", $qtarget, "NA", "NA", "true");
			}
			else{
				Search_seq($qfasta, $qseqid, "NA", "NA", $qtarget, $hpos0, $hpos1, "true");		#masking non-blast-hit regions with "NNN"
			}
		}
		print "! output [$qtarget]...\n";
		
		SKIP:{
			my $skip = 1;
		}
		
	#	my $apref = "q-".$dbpref."_".$seqid."_target-".$qpref;
		my $maf = "last_".$apref.".maf";
		my $psl = "last_".$apref.".psl";
		my $lastal_e = 25;
		my $lastal_q = 3;
		my $lastal_j = 4;
		my $lastal_a = 1;
		my $lastal_b = 1;
		my $lastsp_s = 35;
		
		unless(-e $psl){
			print "\n";
			LASTINDEX($qtarget);
			my $cmd3 = "lastal -e $lastal_e -q $lastal_q -j $lastal_j -P $cpu -a $lastal_a -b $lastal_b $qtarget $seqfasta | last-split -s $lastsp_s > $maf";
			print "! cmd = [$cmd3]\n";
			system("$cmd3");
			
			my $cntfail_mafconv = 0;
			RedoMafconv:{
				my $RedoMafconv = 1;
			}
			
			my $cmd4 = "maf-convert psl $maf > $psl";
			print "! cmd = [$cmd4]\n";
			if(system("$cmd4") != 0){
				$cntfail_mafconv++;
				if($cntfail_mafconv <= 3){
					goto RedoMafconv;
				}
			}
			
			Delete("$qtarget.bck");
			Delete("$qtarget.des");
			Delete("$qtarget.prj");
			Delete("$qtarget.sds");
			Delete("$qtarget.ssp");
			Delete("$qtarget.suf");
			Delete("$qtarget.tis");
		}
		else{
			print "! [$psl] already exists, skipping...\n";
		}
		
		my $rh = Search_targetpos($psl, $lth, $seqid, $tpos0, $tpos1);
		if($rh->{data} eq 'false'){
			$nohit = "true";
			goto JUMP;
		}
		
		$rh->{stats_log} = $log."\n".$rh->{stats_log};
		
		if($keep eq '1'){
			my $log_dir = "log_stats";
			unless(-e $log_dir){
				system("mkdir $log_dir");
			}
			
			SAVE("./$log_dir/stats_$apref.log.txt", $rh->{stats_log});
			SAVE("./$log_dir/stats_$apref.count.tsv", $rh->{stats_count});
			SAVE("./$log_dir/stats_$apref.seqlength.tsv", $rh->{stats_pos});
			SAVE("./$log_dir/stats_$apref.hmatches.tsv", $rh->{stats_hmatches});
			SAVE("./$log_dir/stats_$apref.hmisMatches.tsv", $rh->{stats_hmisMatches});
			SAVE("./$log_dir/stats_$apref.hblockCount.tsv", $rh->{stats_hblockCount});
		}
		
		my $hposdata = $rh->{hposdata};
		my @TID = keys(%{$hposdata});
		@TID = sort {$a cmp $b} @TID;
		my $sum_matches_total = $rh->{sum_matches_total};
		
		my $AoAL = [];
		foreach my $tid (@TID){
			unless($tid =~ /All/){
				my @AL;
				push(@AL, $tid);
				push(@AL, $hposdata->{$tid}{hmatches});
				push(@AL, $hposdata->{$tid}{qpos_min});
				push(@AL, $hposdata->{$tid}{qpos_max});
				push(@AL, $hposdata->{$tid}{tpos_min});
				push(@AL, $hposdata->{$tid}{tpos_max});
				push(@AL, $hposdata->{$tid}{seqlen});
				push(@{$AoAL}, \@AL);
			}
		}
		@{$AoAL} = sort {$b->[1] <=> $a->[1]} @{$AoAL};
		
		my $abhit_summary;
		my $raw_abhit_data;
		my $pav_summary = "directory\tdbfasta\tID-region in dbfasta\tpav pos0\tpav pos1\tpav dist\tqfasta\tqfasta seqid\tpav pos0\tpav pos1\tpac dist\tnote\n";
		my $num_AoAL = @{$AoAL};
		if($num_AoAL == 0){
			$abhit_summary = "$dir\t$dbpref\t$seqid\t$tpos0\t$tpos1\t$ab0\t$ab1\t$qpref\tno hit\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t$note\t-\t-\t-\t-\t-\t-\t-\t-\n";
			APPEND($afile, $abhit_summary);
		#	SAVE($bfile, $pav_summary);
			
			if(! $keep){
				Delete($dbfasta);
				Delete($qfasta);
				Delete($qtarget);
				Delete($seqfasta);
				Delete($blast_result);
				Delete($maf);
				Delete($psl);
			}
			goto END;
		}
		
		print "\n! summary of candidate...\n";
		my $rank = 1;
		foreach my $AL (@{$AoAL}){
			print " [$rank] : [$AL->[0]] [$AL->[4]] - [$AL->[5]] ($AL->[1] bp)\n";
			$rank++;
		}
		print "\n";
		
		JUMP:{
			print "\n";
			my $jump = 1;
		}
		
		if(! -e $afile){
			if($abhit_summary){
				$abhit_summary = $abhit_header.$abhit_summary;
				SAVE($afile, $abhit_summary);
				$abhit_summary = "";
			}
			else{
				SAVE($afile, $abhit_header);
			}
		}
		
		if($nohit ne 'true'){
			print "! summarizing last alignment...\n";
			foreach my $AL (@{$AoAL}){
				my $tid = $AL->[0];
				
				my $tid_alt = $tid;
				if($tid =~ /chr00/){
					$tid_alt =~ s/chr00/unanchored_scaffolds/;
				}
				
				if($hqtarget->{$tid}{seq}){
					if($AL->[4] < 1){
						$AL->[4] = 1;
					}
					
					my $qtpref = $tid."_".$AL->[4]."_".$AL->[5];
					print "! [$seqid] vs [$tid] [$AL->[4]] - [$AL->[5]] ($AL->[1] bp) ...\n";
					
					my @MODES;
					$MODES[0] = "FR";
					$MODES[1] = "F";
					$MODES[2] = "R";
					
					my $Rrh_FR = {};
					my $Rrh_F = {};
					my $Rrh_R = {};
					
					foreach my $mode (@MODES){
						my $num_redo = 0;
						my $window_set_factor = 2;
						my $th_border = 500;
						
						Redo:{
							if($num_redo == 1){
								$window_set_factor = 1.5;
								$th_border = 1000;
							}
							elsif($num_redo == 2){
								$window_set_factor = 1;
								$th_border = 1500;
							}
							elsif($num_redo == 3){
								$window_set_factor = 0.5;
								$th_border = 2000;
							}
						}
						
						my $Rrh_tmp = R_plot($dir, $apref, $hposdata->{$tid}, $dbpref, $qpref, $seqid, $tid, $qtpref, $lth, $tpos0, $tpos1, $ab0, $ab1, 0, $outplot, $keep, $note, $AL->[6], $nomask, $AL->[4], $AL->[5], $hqtarget->{$tid}{seq}, $tstrand, $neighbor_tlen, $flen, $window_set_factor, $num_redo, $th_border, $mode);
						
						if($Rrh_tmp->{judge} && $Rrh_tmp->{judge} eq 'retry' && $num_redo < 3){
							print "! retry analysis ($num_redo)...\n";
							$num_redo++;
							goto Redo;
						}
						else{
							if($mode eq 'FR'){
								$Rrh_FR = $Rrh_tmp;
								
								if(defined $Rrh_tmp->{abhit_score} && $Rrh_tmp->{abhit_score} eq '100'){
									last;
								}
							}
							elsif($mode eq 'F'){
								$Rrh_F = $Rrh_tmp;
							}
							elsif($mode eq 'R'){
								$Rrh_R = $Rrh_tmp;
							}
						}
					}
					
					my $Rrh = SelectBest($Rrh_FR, $Rrh_F, $Rrh_R);
					if($Rrh->{abhit_summary}){
						if($Rrh->{raw}){
							$abhit_summary .= $Rrh->{abhit_summary}."\n";
							$raw_abhit_data .= $Rrh->{raw};
						}
						else{
							$abhit_summary .= $Rrh->{abhit_summary};
						}
					}
					if($Rrh->{pav_summary}){
						$pav_summary .= $Rrh->{pav_summary};
					}
					
					if($allplot eq '0'){
						last;
					}
				}
				elsif($hqtarget->{$tid_alt}{seq}){
					if($AL->[4] < 1){
						$AL->[4] = 1;
					}
					
					my $qtpref = $tid_alt."_".$AL->[4]."_".$AL->[5];
					print "! [$seqid] vs [$tid_alt] [$AL->[4]] - [$AL->[5]] ($AL->[1] bp) ...\n";
					
					my @MODES;
					$MODES[0] = "FR";
					$MODES[1] = "F";
					$MODES[2] = "R";
					
					my $Rrh_FR = {};
					my $Rrh_F = {};
					my $Rrh_R = {};
					
					foreach my $mode (@MODES){
						my $num_redo = 0;
						my $window_set_factor = 2;
						my $th_border = 500;
						
						Redo:{
							if($num_redo == 1){
								$window_set_factor = 1.5;
								$th_border = 1000;
							}
							elsif($num_redo == 2){
								$window_set_factor = 1;
								$th_border = 1500;
							}
							elsif($num_redo == 3){
								$window_set_factor = 0.5;
								$th_border = 2000;
							}
						}
						
						my $Rrh_tmp = R_plot($dir, $apref, $hposdata->{$tid_alt}, $dbpref, $qpref, $seqid, $tid_alt, $qtpref, $lth, $tpos0, $tpos1, $ab0, $ab1, 0, $outplot, $keep, $note, $AL->[6], $nomask, $AL->[4], $AL->[5], $hqtarget->{$tid_alt}{seq}, $tstrand, $neighbor_tlen, $flen, $window_set_factor, $num_redo, $th_border, $mode);
						
						if($Rrh_tmp->{judge} && $Rrh_tmp->{judge} eq 'retry' && $num_redo < 3){
							print "! retry analysis ($num_redo)...\n";
							$num_redo++;
							goto Redo;
						}
						else{
							if($mode eq 'FR'){
								$Rrh_FR = $Rrh_tmp;
								
								if(defined $Rrh_tmp->{abhit_score} && $Rrh_tmp->{abhit_score} eq '100'){
									last;
								}
							}
							elsif($mode eq 'F'){
								$Rrh_F = $Rrh_tmp;
							}
							elsif($mode eq 'R'){
								$Rrh_R = $Rrh_tmp;
							}
						}
					}
					
					my $Rrh = SelectBest($Rrh_FR, $Rrh_F, $Rrh_R);
					if($Rrh->{abhit_summary}){
						if($Rrh->{raw}){
							$abhit_summary .= $Rrh->{abhit_summary}."\n";
							$raw_abhit_data .= $Rrh->{raw};
						}
						else{
							$abhit_summary .= $Rrh->{abhit_summary};
						}
					}
					if($Rrh->{pav_summary}){
						$pav_summary .= $Rrh->{pav_summary};
					}
					
					if($allplot eq '0'){
						last;
					}
				}
				else{
					print "! missing seq for [$tid] or [$tid_alt] ...\n";
					$abhit_summary .= "$dir\t$dbpref\t$seqid\t$tpos0\t$tpos1\t$ab0\t$ab1\t$qpref\tmissing seq for [$tid] or [$tid_alt]\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t$note\n";
				}
			}
		}
		else{
			$abhit_summary .= "$dir\t$dbpref\t$seqid\t$tpos0\t$tpos1\t$ab0\t$ab1\t$qpref\tno hit\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t$note\t-\t-\t-\t-\t-\t-\t-\t-\n";
		}
		
		APPEND($afile, $abhit_summary);
		APPEND($rdfile, $raw_abhit_data);
	#	SAVE($bfile, $pav_summary);
		if($note && $note ne '-'){
			print "! [$note] done\n";
		}
		
		if(! $keep){
			Delete($dbfasta);
			Delete($qfasta);
			Delete($qtarget);
			Delete($seqfasta);
			Delete($blast_result);
			Delete($maf);
			Delete($psl);
		}
	}
	
	if(-e $dbfasta){
		system("rm $dbfasta");
	}
	if(-e $qfasta){
		system("rm $qfasta");
	}
}

END:{
	print "\n";
	my $end = 1;
}


################################################################################
#-------------------------------------------------------------------------------
sub SelectBest{
my $Rrh_FR = shift;
my $Rrh_F = shift;
my $Rrh_R = shift;

my $debug = 0;
if($debug eq '1'){
	print "\n--\n";
}

my $AoD = [];
if(defined $Rrh_FR->{abhit_summary} && defined $Rrh_FR->{abhit_score} && $Rrh_FR->{abhit_score} =~ /\d/){
	my $htmp = $Rrh_FR;
	my @S = split(/\n/, $htmp->{abhit_score});
	my @D = split(/\n/, $htmp->{abhit_summary});
	my $numS = @S;
	my $numD = @D;
	
	if($debug eq '1'){
		print "FR | $htmp->{abhit_summary}\n";
	}
	
	for(my $i = 0; $i < $numD; $i++){
		if(defined $S[$i] && defined $D[$i]){
			my @tmp;
			push(@tmp, $S[$i]);
			push(@tmp, $D[$i]);
			push(@tmp, $htmp->{judge});
			push(@tmp, 3);
			push(@tmp, "FR");
			push(@{$AoD}, \@tmp);
		}
	}
}
elsif($debug eq '1'){
	print "FR | missing entry\n";
}

if(defined $Rrh_F->{abhit_summary} && defined $Rrh_F->{abhit_score} && $Rrh_F->{abhit_score} =~ /\d/){
	my $htmp = $Rrh_F;
	my @S = split(/\n/, $htmp->{abhit_score});
	my @D = split(/\n/, $htmp->{abhit_summary});
	my $numS = @S;
	my $numD = @D;
	
	if($debug eq '1'){
		print "F | $htmp->{abhit_summary}\n";
	}
	
	for(my $i = 0; $i < $numD; $i++){
		if(defined $S[$i] && defined $D[$i]){
			my @tmp;
			push(@tmp, $S[$i]);
			push(@tmp, $D[$i]);
			push(@tmp, $htmp->{judge});
			push(@tmp, 2);
			push(@tmp, "F");
			push(@{$AoD}, \@tmp);
		}
	}
}
elsif($debug eq '1'){
	print "F | missing entry\n";
}

if(defined $Rrh_R->{abhit_summary} && defined $Rrh_R->{abhit_score} && $Rrh_R->{abhit_score} =~ /\d/){
	my $htmp = $Rrh_R;
	my @S = split(/\n/, $htmp->{abhit_score});
	my @D = split(/\n/, $htmp->{abhit_summary});
	my $numS = @S;
	my $numD = @D;
	
	if($debug eq '1'){
		print "R | $htmp->{abhit_summary}\n";
	}
	
	for(my $i = 0; $i < $numD; $i++){
		if(defined $S[$i] && defined $D[$i]){
			my @tmp;
			push(@tmp, $S[$i]);
			push(@tmp, $D[$i]);
			push(@tmp, $htmp->{judge});
			push(@tmp, 1);
			push(@tmp, "R");
			push(@{$AoD}, \@tmp);
		}
	}
}
elsif($debug eq '1'){
	print "R | missing entry\n";
}

my $num_AoD = @{$AoD};

if($num_AoD > 0){
	@{$AoD} = sort {$b->[3] <=> $a->[3]} @{$AoD};
	@{$AoD} = sort {$b->[0] <=> $a->[0]} @{$AoD};
	
	my $raw_data = "";
	for(my $i = 0; $i < $num_AoD; $i++){
		print "! [$i] = [$AoD->[$i][4]] | score = [$AoD->[$i][0]]\n";
		
		if($i == 0){
			$raw_data .= $AoD->[$i][1]."\ttop\n";
		}
		elsif($i == 0){
			$raw_data .= $AoD->[$i][1]."\t$i\n";
		}
	}
	
	my $tophit = {};
	$tophit->{abhit_score} = $AoD->[0][0];
	$tophit->{abhit_summary} = $AoD->[0][1];
	$tophit->{judge} = $AoD->[0][2];
	$tophit->{raw} = $raw_data;
	
	return $tophit;
}
else{
	if($debug eq '1'){
		print "! missing array of data...\n";
	}
	return $Rrh_FR;
}

}


#-------------------------------------------------------------------------------
sub Search_hintpsl{
my ($qpsl, $seqid, $tpos0, $tpos1, $ab0, $ab1, $lth) = @_;

open(my $fh, "<", $qpsl) or die;
my $hash = {};
my $sw = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	unless($line){
		next;
	}
	
	my @A = split(/\t/, $line);
	if($A[13] eq $seqid && $A[0] >= $lth){
		my $p0 = $A[15];
		my $p1 = $A[16];
		if($A[15] > $A[16]){
			$p1 = $A[15];
			$p0 = $A[16];
		}
		
		$A[15] = $p0;
		$A[16] = $p1;
		
		if($p0 <= $tpos0 && $tpos0 <= $p1){
			$hash->{$A[9]}{psl} .= join("\t", @A)."\n";
			$hash->{$A[9]}{matches} += $A[0];
			$sw = 1;
		}
		elsif($sw == 1 && $p0 <= $tpos0 && $p1 <= $tpos0){
			$hash->{$A[9]}{psl} .= join("\t", @A)."\n";
			$hash->{$A[9]}{matches} += $A[0];
		}
		elsif($sw == 1 && $p0 <= $tpos1 && $tpos1 <= $p1){
			$hash->{$A[9]}{psl} .= join("\t", @A)."\n";
			$hash->{$A[9]}{matches} += $A[0];
		}
		elsif($sw == 1 && $tpos1 < $p0 && $tpos1 < $p1){
			last;
		}
	}
}
close $fh;

my $rh = {};
unless($hash){
	return ("false", "NA", 0, 0, $rh);
}

my @Cands = keys(%{$hash});
@Cands = sort {$a cmp $b} @Cands;

my $AoC = [];
foreach my $cand (@Cands){
	my @tmp;
	push(@tmp, $cand);
	push(@tmp, $hash->{$cand}{matches});
	push(@{$AoC}, \@tmp);
}

@{$AoC} = sort {$b->[1] <=> $a->[1]} @{$AoC};
my $num_AoC = @{$AoC};
if($num_AoC == 0){
	return ("false", "NA", 0, 0, $rh);
}

print "\n! rank by seq ID\n";
my $n = 1;
foreach my $tmp (@{$AoC}){
	print " [$n] = [$tmp->[0]] [$tmp->[1]] bp\n";
	$n++;
	if($n == 6){
		last;
	}
}

my $tophit = $AoC->[0][0];
my @PSL = split(/\n/, $hash->{$tophit}{psl});
my $AoQ = [];
foreach my $line (@PSL){
	my @A = split(/\t/, $line);
	push(@{$AoQ}, \@A);
}
@{$AoQ} = sort {$a->[11] <=> $b->[11]} @{$AoQ};

my $prev_pos = 0;
$n = 0;
my $ngrp = 0;
my $htop = {};
foreach my $A (@{$AoQ}){
	print " $ngrp : $A->[11] - $A->[12] | $A->[15] - $A->[16]\n";
	if($n == 0){
		$htop->{$ngrp}{psl} .= join("\t", @{$A})."\n";
		$htop->{$ngrp}{matches} += $A->[0];
		$htop->{$ngrp}{Pos} .= $A->[11]."\n";
		$htop->{$ngrp}{Pos} .= $A->[12]."\n";
	}
	else{
		my $window = ($tpos1 - $tpos0) / 2;
		my $dist = abs($A->[11] - $prev_pos);
		if($dist > $window){
			$ngrp++;
		}
		
		$htop->{$ngrp}{psl} .= join("\t", @{$A})."\n";
		$htop->{$ngrp}{matches} += $A->[0];
		$htop->{$ngrp}{Pos} .= $A->[11]."\n";
		$htop->{$ngrp}{Pos} .= $A->[12]."\n";
	}
	$prev_pos = $A->[12];
	$n++;
}

my @GRP = keys(%{$htop});
@GRP = sort {$a <=> $b} @GRP;

my $AoC2 = [];
foreach my $ngrp (@GRP){
	my @tmp;
	push(@tmp, $ngrp);
	push(@tmp, $htop->{$ngrp}{matches});
	push(@tmp, $htop->{$ngrp}{Pos});
	push(@{$AoC2}, \@tmp);
}
@{$AoC2} = sort {$b->[1] <=> $a->[1]} @{$AoC2};

print "\n! rank by algnment block\n";
$n = 1;
foreach my $tmp (@{$AoC2}){
	my @POS = split(/\n/, $tmp->[2]);
	@POS = sort {$a <=> $b} @POS;
	my $npos = @POS;
	
	print " [$n] = [$tophit : $POS[0] - $POS[$npos - 1]] [$tmp->[1]] bp\n";
	$n++;
	if($n == 6){
		last;
	}
}

my $ntop = $AoC2->[0][0];
@PSL = split(/\n/, $htop->{$ntop}{psl});
my $AoAl = [];
my $AoFw = [];
my $AoRv = [];
my $AoAlsp = [];
my $AoFwsp = [];
my $AoRvsp = [];
my $hstrand = {};
my @dbposFw;
my @dbposRv;
my @qfposFw;
my @qfposRv;
foreach my $line (@PSL){
	my @A = split(/\t/, $line);
	my $strand = $A[8];
	my $blockSizes = $A[18];
	my $qStarts = $A[19];
	my $tStarts = $A[20];				# dbfasta
	
	my @BLK = split(/,/, $blockSizes);
	my @Qs = split(/,/, $qStarts);
	my @Ts = split(/,/, $tStarts);		# dbfasta
	my $num_BLK = @BLK;
	
	for(my $i = 0; $i < $num_BLK; $i++){
		my $blen = $BLK[$i];
		my $qp0 = $Qs[$i];
		my $tp0 = $Ts[$i];
		
		for(my $j = 0; $j < $blen; $j++){
			my $qpos = $qp0 + $j;
			my $tpos = $tp0 + $j;
			
			if($tpos0 <= $tpos && $tpos <= $tpos1){
				my @tmp;
				push(@tmp, $qpos);
				push(@tmp, $tpos);
				push(@{$AoAl}, \@tmp);
				$hstrand->{$strand} = 1;
				
				my $tpos_sub = $tpos - $tpos0;
				if($strand eq '+'){
					push(@{$AoFw}, \@tmp);
					push(@dbposFw, $tpos_sub);
					push(@qfposFw, $qpos);
				}
				elsif($strand eq '-'){
					push(@{$AoRv}, \@tmp);
					push(@dbposRv, $tpos_sub);
					push(@qfposRv, $qpos);
				}
			}
			if($ab0 <= $tpos && $tpos <= $ab1){
				my @tmp;
				push(@tmp, $qpos);
				push(@tmp, $tpos);
				push(@{$AoAlsp}, \@tmp);
				
				if($strand eq '+'){
					push(@{$AoFwsp}, \@tmp);
				}
				elsif($strand eq '-'){
					push(@{$AoRvsp}, \@tmp);
				}
			}
		}
	}
}

my @Strands = keys(%{$hstrand});
@Strands = sort {$a cmp $b} @Strands;
my $strands = join(",", @Strands);

print "\n! matching position (including flanking)\n";
my $num_Al = @{$AoAl};
@{$AoAl} = sort {$a->[1] <=> $b->[1]} @{$AoAl};
print " DB-based | $seqid: $AoAl->[0][1] - $AoAl->[$num_Al - 1][1] = $tophit : $AoAl->[0][0] - $AoAl->[$num_Al - 1][0] ($strands)\n";
@{$AoAl} = sort {$a->[0] <=> $b->[0]} @{$AoAl};
print " q-based  | $seqid: $AoAl->[0][1] - $AoAl->[$num_Al - 1][1] = $tophit : $AoAl->[0][0] - $AoAl->[$num_Al - 1][0] ($strands)\n";

my $num_Alsp = @{$AoAlsp};
if($AoAlsp->[0][0] && $AoAlsp->[$num_Alsp - 1][0]){
	print "\n! matching position (specified region)\n";
	@{$AoAlsp} = sort {$a->[1] <=> $b->[1]} @{$AoAlsp};
	print " DB-based | $seqid: $AoAlsp->[0][1] - $AoAlsp->[$num_Alsp - 1][1] = $tophit : $AoAlsp->[0][0] - $AoAlsp->[$num_Alsp - 1][0] ($strands)\n";
	@{$AoAlsp} = sort {$a->[0] <=> $b->[0]} @{$AoAlsp};
	print " q-based  | $seqid: $AoAlsp->[0][1] - $AoAlsp->[$num_Alsp - 1][1] = $tophit : $AoAlsp->[0][0] - $AoAlsp->[$num_Alsp - 1][0] ($strands)\n";
}

my $hpos0 = $AoAl->[0][0];
my $hpos1 = $AoAl->[$num_Al - 1][0];
if($hpos1 < $hpos0){
	$hpos1 = $AoAl->[0][0];
	$hpos0 = $AoAl->[$num_Al - 1][0];
}

# to use R_plot function, set @dbposFw to {qposFw} (looks opposite but it is not mistake)
$rh->{qposFw} = \@dbposFw;		#dbseq (= $dbpref)
$rh->{qposRv} = \@dbposRv;		#dbseq (= $dbpref)
$rh->{tposFw} = \@qfposFw;		#target seq (= $qpref)
$rh->{tposRv} = \@qfposRv;		#target seq (= $qpref)

return ("true", $tophit, $hpos0, $hpos1, $rh);
}


#-------------------------------------------------------------------------------
sub R_plot{
my $dir = shift;
my $prefix = shift;
my $rh = shift;
my $dbpref = shift;
my $qpref = shift;
my $seqid = shift;
my $qseqid = shift;
my $qtpref = shift;
my $lth = shift;
my $tpos0 = shift;		#dbseq (= $dbpref)
my $tpos1 = shift;		#dbseq (= $dbpref)
my $ab0 = shift;
my $ab1 = shift;
my $allplot = shift;
my $outplot = shift;
my $keep = shift;
my $note = shift;
my $hseqlen = shift;
my $nomask = shift;
my $qtpos0 = shift;		#target seq (= $qpref)
my $qtpos1 = shift;		#target seq (= $qpref)
my $qtseq = shift;		#target seq (= $qpref)
my $tstrand = shift;
my $neighbor_tlen = shift;
my $flen = shift;
my $window_set_factor = shift;
my $num_redo = shift;
my $th_border = shift;
my $mode = shift;

$allplot = 0;		# not use '1'

if(! defined $window_set_factor || ! $window_set_factor){
	$window_set_factor = 1;
}
elsif($window_set_factor < 0.1){
	$window_set_factor = 0.1;
}

my $qposFw = $rh->{qposFw};		#dbseq (= $dbpref)
my $qposRv = $rh->{qposRv};		#dbseq (= $dbpref)
my $tposFw = $rh->{tposFw};		#target seq (= $qpref)
my $tposRv = $rh->{tposRv};		#target seq (= $qpref)

my $num_qposFw = @{$qposFw};
my $num_qposRv = @{$qposRv};
my $num_tposFw = @{$tposFw};
my $num_tposRv = @{$tposRv};

unless($num_qposFw == $num_tposFw){
	print "\n! error | counts unmatched between query and target (Fw) : [$num_qposFw] [$num_tposFw]\n";
	return;
}
unless($num_qposRv == $num_tposRv){
	print "\n! error | counts unmatched between query and target (Rv) : [$num_qposRv] [$num_tposRv]\n";
	return;
}

my @Rqpos;
my @Rtpos;
if($mode eq 'FR' || $mode eq 'F'){
	if($qposFw && $tposFw){
		my $P0 = $qposFw;
		my $P1 = $tposFw;
		
		for(my $i = 0; $i < $num_qposFw; $i++){
			if(($P0->[$i] + $tpos0) >= $tpos0 && ($P0->[$i] + $tpos0) <= $tpos1){
				push(@Rqpos, $P0->[$i]);
				push(@Rtpos, $P1->[$i]);
			}
		}
	}
}

if($mode eq 'FR' || $mode eq 'R'){
	if($qposRv && $tposRv){
		my $P0 = $qposRv;
		my $P1 = $tposRv;
		
		for(my $i = 0; $i < $num_qposRv; $i++){
			if(($P0->[$i] + $tpos0) >= $tpos0 && ($P0->[$i] + $tpos0) <= $tpos1){
				push(@Rqpos, $P0->[$i]);
				push(@Rtpos, $P1->[$i]);
			}
		}
	}
}

my $window = int(($tpos1 - $tpos0) / $window_set_factor);
if($window > 100 * 1000){
	$window = 100 * 1000;
}

@Rtpos = sort {$a <=> $b} @Rtpos;
my $numRtpos = @Rtpos;
my $ngrp = 1;
my $hngrp = {};
for(my $i = 0; $i < $numRtpos - 1; $i++){
	my $j = $i + 1;
	my $p0 = $Rtpos[$i];
	my $p1 = $Rtpos[$j];
	
	if($p1 - $p0 <= $window){
		$hngrp->{dat}{$ngrp}{$p0} = 1;
		$hngrp->{dat}{$ngrp}{$p1} = 1;
		$hngrp->{grp}{$p0} = $ngrp;
		$hngrp->{grp}{$p1} = $ngrp;
		$hngrp->{cnt}{$ngrp} += 1;
	}
	else{
		$hngrp->{dat}{$ngrp}{$p0} = 1;
		$hngrp->{grp}{$p0} = $ngrp;
		$hngrp->{cnt}{$ngrp} += 1;
		
		$ngrp++;
		$hngrp->{dat}{$ngrp}{$p1} = 1;
		$hngrp->{grp}{$p1} = $ngrp;
		$hngrp->{cnt}{$ngrp} += 1;
	}
}

my $htmp1 = $hngrp->{dat};
my @nG = keys(%{$htmp1});
@nG = sort {$a <=> $b} @nG;

my $AoG = [];
foreach my $grp (@nG){
	my @G;
	push(@G, $grp);
	push(@G, $hngrp->{cnt}{$grp});
	push(@{$AoG}, \@G);
}
@{$AoG} = sort {$b->[1] <=> $a->[1]} @{$AoG};

my $abhit_summary;
my $pav_summary;
my $abhit_score;
my $abhit_judge;

my $Rrh = {};
$Rrh->{abhit_summary} = "";
$Rrh->{pav_summary} = "";
$Rrh->{judge} = "null";

foreach my $G (@{$AoG}){
	my $grp = $G->[0];
	
#	print "! group [$grp] [$G->[1]] bp ...\n";
	print "! $mode | group [$grp] ...\n";
	my $pfile_Fw = "cumpos_".$prefix."_".$qtpref."_".$grp."_Fw_lth".$lth.".tsv";
	my $pfile_Rv = "cumpos_".$prefix."_".$qtpref."_".$grp."_Rv_lth".$lth.".tsv";
	my $png = "plot_".$prefix."_".$qtpref."_".$grp.".png";
	my $png1 = "plot1_".$prefix."_".$qtpref."_".$grp.".png";
	my $png2 = "plot2_".$prefix."_".$qtpref."_".$grp.".png";
	my $png3 = "plot3_".$prefix."_".$qtpref."_".$grp.".png";
	my $png4 = "plot4_".$prefix."_".$qtpref."_".$grp.".png";

	my @each_Rqpos;
	my @each_Rtpos;
	my $AoP = [];
	my @P1ab;
	my $abhit_cnt = 0;
	my $n5hit_cnt = 0;			# hit count of 5' neighbor
	my $n3hit_cnt = 0;			# hit count of 3' neighbor
	
	my $ab0nb = $ab0 - $neighbor_tlen;
	my $ab1nb = $ab1 + $neighbor_tlen;
	my $prohit_cnt = 0;			# hit count of promoter
	my $utrhit_cnt = 0;			# hit count of 3'-UTR
	
	my $p0ab0_border = "null";
	my $p0ab1_border = "null";
	my $p1ab0_border = "null";
	my $p1ab1_border = "null";
	
	# determines alignment borders of genic region first ---------------------//
	my @P1regionPos;
	if($qposFw && $tposFw){
		my $P0 = $qposFw;		#dbseq (= $dbpref)
		my $P1 = $tposFw;		#target seq (= $qpref)
		
		for(my $i = 0; $i < $num_qposFw; $i++){
			if($hngrp->{grp}{$P1->[$i]} && $grp eq $hngrp->{grp}{$P1->[$i]}){
				my $modp0 = $P0->[$i] + $tpos0;
				if($tpos0 <= $modp0 && $modp0 <= $tpos1){
					if($ab0 - $th_border <= $modp0 && $modp0 <= $ab0 + $th_border){
						if($p1ab0_border eq 'null'){
							$p0ab0_border = $modp0;
							$p1ab0_border = $P1->[$i];
						}
						elsif( abs($modp0 - $ab0) < abs($p0ab0_border - $ab0) ){
							$p0ab0_border = $modp0;
							$p1ab0_border = $P1->[$i];
						}
					}
					if($ab1 - $th_border <= $modp0 && $modp0 <= $ab1 + $th_border){
						if($p1ab1_border eq 'null'){
							$p0ab1_border = $modp0;
							$p1ab1_border = $P1->[$i];
						}
						elsif( abs($modp0 - $ab1) < abs($p0ab1_border - $ab1) ){
							$p0ab1_border = $modp0;
							$p1ab1_border = $P1->[$i];
						}
					}
					if($ab0 <= $modp0 && $modp0 <= $ab1){
						push(@P1regionPos, $P1->[$i]);
					}
				}
			}
		}
	}

	if($qposRv && $tposRv){
		my $P0 = $qposRv;		#dbseq (= $dbpref)
		my $P1 = $tposRv;		#target seq (= $qpref)
		
		for(my $i = 0; $i < $num_qposRv; $i++){
			if($hngrp->{grp}{$P1->[$i]} && $grp eq $hngrp->{grp}{$P1->[$i]}){
				my $modp0 = $P0->[$i] + $tpos0;
				if($tpos0 <= $modp0 && $modp0 <= $tpos1){
					if($ab0 - $th_border <= $modp0 && $modp0 <= $ab0 + $th_border){
						if($p1ab0_border eq 'null'){
							$p0ab0_border = $modp0;
							$p1ab0_border = $P1->[$i];
						}
						elsif( abs($modp0 - $ab0) < abs($p0ab0_border - $ab0) ){
							$p0ab0_border = $modp0;
							$p1ab0_border = $P1->[$i];
						}
					}
					if($ab1 - $th_border <= $modp0 && $modp0 <= $ab1 + $th_border){
						if($p1ab1_border eq 'null'){
							$p0ab1_border = $modp0;
							$p1ab1_border = $P1->[$i];
						}
						elsif( abs($modp0 - $ab1) < abs($p0ab1_border - $ab1) ){
							$p0ab1_border = $modp0;
							$p1ab1_border = $P1->[$i];
						}
					}
					if($ab0 <= $modp0 && $modp0 <= $ab1){
						push(@P1regionPos, $P1->[$i]);
					}
				}
			}
		}
	}
	
	if($p1ab0_border eq 'null' || $p1ab1_border eq 'null'){
		if($num_redo < 3){
			$Rrh->{judge} = "retry";
			return $Rrh;
		}
	}
	
	my $num_P1regionPos = @P1regionPos;
	if($num_P1regionPos > 0){
		@P1regionPos = sort {$a <=> $b} @P1regionPos;
		if(defined $P1regionPos[0]){
			my $P1regionPos0 = $P1regionPos[0];
			if($p1ab0_border eq 'null'){
				$p1ab0_border = $P1regionPos0;
			}
			elsif($P1regionPos0 < $p1ab0_border){
				$p1ab0_border = $P1regionPos0;
			}
		}
		
		@P1regionPos = sort {$b <=> $a} @P1regionPos;
		if(defined $P1regionPos[0]){
			my $P1regionPos1 = $P1regionPos[0];
			if($p1ab1_border eq 'null'){
				$p1ab1_border = $P1regionPos1;
			}
			elsif($p1ab1_border < $P1regionPos1){
				$p1ab1_border = $P1regionPos1;
			}
		}
	}
	
	if(! defined $p1ab0_border){
		$p1ab0_border = "null";
	}
	if(! defined $p1ab1_border){
		$p1ab1_border = "null";
	}
	
	# then collect info ---------------------//
	my @P1promoterPos;
	my @P1utrPos;
	if($qposFw && $tposFw){
		my $P0 = $qposFw;		#dbseq (= $dbpref)
		my $P1 = $tposFw;		#target seq (= $qpref)
		
		my $r = "dbfasta\ttarget\n";
		for(my $i = 0; $i < $num_qposFw; $i++){
			if($hngrp->{grp}{$P1->[$i]} && $grp eq $hngrp->{grp}{$P1->[$i]}){
				my $modp0 = $P0->[$i] + $tpos0;
				if($tpos0 <= $modp0 && $modp0 <= $tpos1){
					$r .= $modp0."\t".$P1->[$i]."\n";
					push(@each_Rqpos, $modp0);
					push(@each_Rtpos, $P1->[$i]);
					
					my @ptmp;
					push(@ptmp, $modp0);
					push(@ptmp, $P1->[$i]);
					push(@{$AoP}, \@ptmp);
					
					if($p1ab0_border ne 'null' && $p1ab1_border ne 'null'){
						if($ab0 <= $modp0 && $modp0 <= $ab1 && $p1ab0_border <= $P1->[$i] && $P1->[$i] <= $p1ab1_border){
							push(@P1ab, $P1->[$i]);
							$abhit_cnt++;
						}
					}
					else{
						if($ab0 <= $modp0 && $modp0 <= $ab1){
							push(@P1ab, $P1->[$i]);
							$abhit_cnt++;
						}
					}
					
					if($flen > 0){
						if($modp0 < $ab0 && $modp0 < $ab1){
							$n5hit_cnt++;
						}
						elsif($modp0 > $ab0 && $modp0 > $ab1){
							$n3hit_cnt++;
						}
					}
					if($tstrand eq '+'){
						if($ab0nb <= $modp0 && $modp0 < $ab0){
							push(@P1promoterPos, $P1->[$i]);
							$prohit_cnt++;
						}
						elsif($ab1 < $modp0 && $modp0 <= $ab1nb){
							push(@P1utrPos, $P1->[$i]);
							$utrhit_cnt++;
						}
					}
					elsif($tstrand eq '-'){
						if($ab1 < $modp0 && $modp0 <= $ab1nb){
							push(@P1promoterPos, $P1->[$i]);
							$prohit_cnt++;
						}
						elsif($ab0nb <= $modp0 && $modp0 < $ab0){
							push(@P1utrPos, $P1->[$i]);
							$utrhit_cnt++;
						}
					}
				}
			}
		}

		open(my $rfh, ">", $pfile_Fw);
		print $rfh $r;
		close $rfh;
	#	print "! output [$rfile]\n";
	}

	if($qposRv && $tposRv){
		my $P0 = $qposRv;		#dbseq (= $dbpref)
		my $P1 = $tposRv;		#target seq (= $qpref)
		
		my $r = "dbfasta\ttarget\n";
		for(my $i = 0; $i < $num_qposRv; $i++){
			if($hngrp->{grp}{$P1->[$i]} && $grp eq $hngrp->{grp}{$P1->[$i]}){
				my $modp0 = $P0->[$i] + $tpos0;
				if($tpos0 <= $modp0 && $modp0 <= $tpos1){
					$r .= $modp0."\t".$P1->[$i]."\n";
					push(@each_Rqpos, $modp0);
					push(@each_Rtpos, $P1->[$i]);
					
					my @ptmp;
					push(@ptmp, $modp0);
					push(@ptmp, $P1->[$i]);
					push(@{$AoP}, \@ptmp);
					
					if($p1ab0_border ne 'null' && $p1ab1_border ne 'null'){
						if($ab0 <= $modp0 && $modp0 <= $ab1 && $p1ab0_border <= $P1->[$i] && $P1->[$i] <= $p1ab1_border){
							push(@P1ab, $P1->[$i]);
							$abhit_cnt++;
						}
					}
					else{
						if($ab0 <= $modp0 && $modp0 <= $ab1){
							push(@P1ab, $P1->[$i]);
							$abhit_cnt++;
						}
					}
					
					if($flen > 0){
						if($modp0 < $ab0 && $modp0 < $ab1){
							$n5hit_cnt++;
						}
						elsif($modp0 > $ab0 && $modp0 > $ab1){
							$n3hit_cnt++;
						}
					}
					if($tstrand eq '+'){
						if($ab0nb <= $modp0 && $modp0 < $ab0){
							push(@P1promoterPos, $P1->[$i]);
							$prohit_cnt++;
						}
						elsif($ab1 < $modp0 && $modp0 <= $ab1nb){
							push(@P1utrPos, $P1->[$i]);
							$utrhit_cnt++;
						}
					}
					elsif($tstrand eq '-'){
						if($ab1 < $modp0 && $modp0 <= $ab1nb){
							push(@P1promoterPos, $P1->[$i]);
							$prohit_cnt++;
						}
						elsif($ab0nb <= $modp0 && $modp0 < $ab0){
							push(@P1utrPos, $P1->[$i]);
							$utrhit_cnt++;
						}
					}
				}
			}
		}

		open(my $rfh, ">", $pfile_Rv);
		print $rfh $r;
		close $rfh;
	#	print "! output [$rfile]\n";
	}
	
	my $n5hit_ratio = sprintf("%.3f", $n5hit_cnt / $flen );
	my $n3hit_ratio = sprintf("%.3f", $n3hit_cnt / $flen );
	
	my $prohit_ratio = "-";
	my $utrhit_ratio = "-";
	my $prohit_rspan = "-";
	my $utrhit_rspan = "-";
	if($prohit_cnt > 0){
		$prohit_ratio = sprintf("%.3f", $prohit_cnt / $neighbor_tlen );
		
		@P1promoterPos = sort {$a <=> $b} @P1promoterPos;
		my $tmp_promutrnum = @P1promoterPos;
		
		if($p1ab0_border ne 'null' && $p1ab1_border ne 'null'){
			my @P2Gdist;
			push(@P2Gdist, abs($p1ab0_border - $P1promoterPos[0]));
			push(@P2Gdist, abs($p1ab0_border - $P1promoterPos[$tmp_promutrnum - 1]));
			push(@P2Gdist, abs($p1ab1_border - $P1promoterPos[0]));
			push(@P2Gdist, abs($p1ab1_border - $P1promoterPos[$tmp_promutrnum - 1]));
			@P2Gdist = sort {$a <=> $b} @P2Gdist;
			
			$prohit_rspan = sprintf("%.3f", abs($P1promoterPos[$tmp_promutrnum - 1] - $P1promoterPos[0]) / $neighbor_tlen );
			$prohit_rspan .= ";";
			$prohit_rspan .= sprintf("%.3f", ( abs($P1promoterPos[$tmp_promutrnum - 1] - $P1promoterPos[0]) + $P2Gdist[0] ) / $neighbor_tlen );
			$prohit_rspan .= ";".join(" ", @P2Gdist);
		}
		else{
			$prohit_rspan = sprintf("%.3f", abs($P1promoterPos[$tmp_promutrnum - 1] - $P1promoterPos[0]) / $neighbor_tlen );
			$prohit_rspan .= ";null;null";
		}
		$prohit_ratio .= ";".$prohit_rspan.";".$P1promoterPos[0].";".$P1promoterPos[$tmp_promutrnum - 1];
		# ratio_hit; ratio_span; rev_ratio_span; P2Gdist; span_pos0; span_pos1
	}
	if($utrhit_cnt > 0){
		$utrhit_ratio = sprintf("%.3f", $utrhit_cnt / $neighbor_tlen );
		
		@P1utrPos = sort {$a <=> $b} @P1utrPos;
		my $tmp_promutrnum = @P1utrPos;
		
		if($p1ab0_border ne 'null' && $p1ab1_border ne 'null'){
			my @P2Gdist;
			push(@P2Gdist, abs($p1ab0_border - $P1utrPos[0]));
			push(@P2Gdist, abs($p1ab0_border - $P1utrPos[$tmp_promutrnum - 1]));
			push(@P2Gdist, abs($p1ab1_border - $P1utrPos[0]));
			push(@P2Gdist, abs($p1ab1_border - $P1utrPos[$tmp_promutrnum - 1]));
			@P2Gdist = sort {$a <=> $b} @P2Gdist;
			
			$utrhit_rspan = sprintf("%.3f", abs($P1utrPos[$tmp_promutrnum - 1] - $P1utrPos[0]) / $neighbor_tlen );
			$utrhit_rspan .= ";";
			$utrhit_rspan .= sprintf("%.3f", ( abs($P1utrPos[$tmp_promutrnum - 1] - $P1utrPos[0]) + $P2Gdist[0] ) / $neighbor_tlen );
			$utrhit_rspan .= ";".join(" ", @P2Gdist);
		}
		else{
			$utrhit_rspan = sprintf("%.3f", abs($P1utrPos[$tmp_promutrnum - 1] - $P1utrPos[0]) / $neighbor_tlen );
			$utrhit_rspan .= ";null;null";
		}
		$utrhit_ratio .= ";".$utrhit_rspan.";".$P1utrPos[0].";".$P1utrPos[$tmp_promutrnum - 1];
		# ratio_hit; ratio_span; rev_ratio_span; P2Gdist; span_pos0; span_pos1
	}
	
	if($flen == 0){
		$n5hit_cnt = "-";
		$n3hit_cnt = "-";
		$n5hit_ratio = "-";
		$n3hit_ratio = "-";
	}
	if($tstrand eq 'null'){
		$prohit_cnt = "-";
		$utrhit_cnt = "-";
		$prohit_ratio = "-";
		$utrhit_ratio = "-";
	}
	
	my $abdist = abs($ab1 - $ab0);
	my $del_disruption_ratio = "-";
	if($abdist){
		$del_disruption_ratio = $abhit_cnt / $abdist;
	}
	
	my $P1ab_max = "-";
	my $P1ab_min = "-";
	my $qhit_len0 = "-";
	my $qhit_len1 = "-";
	my $bp_insert = "-";
	my $ins_disruption_ratio = "-";
	my $judge_tp0_to_tp1 = "-";
	if($abhit_cnt > 0){
		@P1ab = sort {$b <=> $a} @P1ab;
		$P1ab_max = $P1ab[0];
		@P1ab = sort {$a <=> $b} @P1ab;
		$P1ab_min = $P1ab[0];
		
		my $seq_p1max_to_min = substr($qtseq, $P1ab_min - 1, $P1ab_max - $P1ab_min + 1);
		$qhit_len0 = length($seq_p1max_to_min);
		$seq_p1max_to_min =~ s/N//gi;
		$qhit_len1 = length($seq_p1max_to_min);
		
		if($qhit_len1 - $abhit_cnt > 10){
			$bp_insert = $qhit_len1 - $abhit_cnt;
			$ins_disruption_ratio = $abhit_cnt / $qhit_len1;
		}
		
		if($ins_disruption_ratio ne '-'){
			if($qhit_len1 - $abhit_cnt > 5000){
				$judge_tp0_to_tp1 = "insertion >5kb";
			}
			elsif($qhit_len1 - $abhit_cnt > 4000){
				$judge_tp0_to_tp1 = "insertion >4kb";
			}
			elsif($qhit_len1 - $abhit_cnt > 3000){
				$judge_tp0_to_tp1 = "insertion >3kb";
			}
			elsif($qhit_len1 - $abhit_cnt > 2000){
				$judge_tp0_to_tp1 = "insertion >2kb";
			}
			elsif($qhit_len1 - $abhit_cnt > 1000){
				$judge_tp0_to_tp1 = "insertion >1kb";
			}
			elsif($qhit_len1 - $abhit_cnt > 500){
				$judge_tp0_to_tp1 = "insertion >500bp";
			}
			elsif($del_disruption_ratio ne '-'){
				if($del_disruption_ratio < 0.1){
					$judge_tp0_to_tp1 = "present (collapsed)";
				}
				elsif($del_disruption_ratio < 0.5){
					$judge_tp0_to_tp1 = "present (collapsed)";
				}
				elsif($del_disruption_ratio < 0.8){
					$judge_tp0_to_tp1 = "present (partly)";
				}
				else{
					$judge_tp0_to_tp1 = "present";
				}
			}
			else{
				$judge_tp0_to_tp1 = "present";
			}
		}
		elsif($del_disruption_ratio ne '-'){
			if($del_disruption_ratio < 0.1){
				$judge_tp0_to_tp1 = "present (collapsed)";
			}
			elsif($del_disruption_ratio < 0.5){
				$judge_tp0_to_tp1 = "present (collapsed)";
			}
			elsif($del_disruption_ratio < 0.8){
				$judge_tp0_to_tp1 = "present (partly)";
			}
			else{
				$judge_tp0_to_tp1 = "present";
			}
		}
	}
	else{
		@{$AoP} = sort {$a->[0] <=> $b->[0]} @{$AoP};
		my $numAoP = @{$AoP};
		my $P1ab_lower;
		my $P1ab_upper;
		for(my $i = 0; $i < $numAoP - 1; $i++){
			my $j = $i + 1;
			
			my $qp0 = $AoP->[$i][0];		#dbseq (= $dbpref)
			my $qp1 = $AoP->[$j][0];		#dbseq (= $dbpref)
			my $tp0 = $AoP->[$i][1];		#target seq (= $qpref)
			my $tp1 = $AoP->[$j][1];		#target seq (= $qpref)
			
			if($qp0 <= $ab0 && $ab1 <= $qp1){
				# ---------------[qp0]----[ab0]---------[ab1]----[qp1]----------------
				# ||||||||||||||||||                              |||||||||||||||||||||
				# ---------------[tp0]     "NNN" or absent       [tp1]----------------
				
				$P1ab_lower = $tp0;
				$P1ab_upper = $tp1;
				last;
			}
		}
		
		if($P1ab_lower && $P1ab_upper && $P1ab_lower =~ /[0-9]/ && $P1ab_upper =~ /[0-9]/){
			my $len_qtseq = length($qtseq);
			
			my $seq_tp0_to_tp1;
			if($P1ab_upper > $P1ab_lower){
				$seq_tp0_to_tp1 = substr($qtseq, $P1ab_lower - 1, $P1ab_upper - $P1ab_lower + 1);
			}
			else{
				$seq_tp0_to_tp1 = substr($qtseq, $P1ab_upper - 1, $P1ab_lower - $P1ab_upper + 1);
			}
			
			my $cptmp;
			if($seq_tp0_to_tp1){
				my $dist_P1ab_ul = abs($P1ab_lower - $P1ab_upper);
				my $ratio_ab01 = $dist_P1ab_ul / abs($ab1 - $ab0);
				
				if($dist_P1ab_ul < 10){
					$judge_tp0_to_tp1 = "absent [$P1ab_lower - $P1ab_upper; $dist_P1ab_ul bp]";
				}
				elsif($ratio_ab01 && $ratio_ab01 < 0.1){
					$judge_tp0_to_tp1 = "absent [$P1ab_lower - $P1ab_upper; $dist_P1ab_ul bp]";
				}
				else{
					$cptmp = $seq_tp0_to_tp1;
					$cptmp =~ s/n//gi;
					
					if(! $cptmp){
						$judge_tp0_to_tp1 = "NNN [$P1ab_lower - $P1ab_upper; $dist_P1ab_ul bp; determined=0 bp]";
					}
					else{
						$cptmp = length($cptmp);
						my $ratio_cptmp = $cptmp / abs($P1ab_lower - $P1ab_upper);
						
						if($ratio_cptmp >= 0.5 && $cptmp >= 2000){
							$judge_tp0_to_tp1 = "no hit but some [$P1ab_lower - $P1ab_upper; $dist_P1ab_ul bp; determined=$cptmp bp]";
						}
						else{
							$judge_tp0_to_tp1 = "NNN [$P1ab_lower - $P1ab_upper; $dist_P1ab_ul bp; determined=$cptmp bp]";
						}
					}
				}
			}
		}
	}
	
	if($ab0 > 0 && $ab1 > 0){
		if($p1ab0_border ne 'null' && $p1ab1_border ne 'null'){			# recalculate important values...
			my @tmp;
			push(@tmp, $p1ab0_border);
			push(@tmp, $p1ab1_border);
			@tmp = sort {$a <=> $b} @tmp;
			
			$p1ab0_border = $tmp[0];
			$p1ab1_border = $tmp[1];
			
			my $seq_p1span = substr($qtseq, $tmp[0] - 1, $tmp[1] - $tmp[0] + 1);
			$seq_p1span =~ s/N//gi;
			$qhit_len1 = length($seq_p1span);
			$qhit_len0 = abs($tmp[1] - $tmp[0] + 1);
			
			if($qhit_len1 - $abhit_cnt > 10){
				$bp_insert = $qhit_len1 - $abhit_cnt;
				$ins_disruption_ratio = $abhit_cnt / $qhit_len1;
			}
			else{
				$bp_insert = "-";
				$ins_disruption_ratio = "-";
			}
			
			if($ins_disruption_ratio ne '-'){
				if($qhit_len1 - $abhit_cnt > 5000){
					$judge_tp0_to_tp1 = "insertion >5kb";
				}
				elsif($qhit_len1 - $abhit_cnt > 4000){
					$judge_tp0_to_tp1 = "insertion >4kb";
				}
				elsif($qhit_len1 - $abhit_cnt > 3000){
					$judge_tp0_to_tp1 = "insertion >3kb";
				}
				elsif($qhit_len1 - $abhit_cnt > 2000){
					$judge_tp0_to_tp1 = "insertion >2kb";
				}
				elsif($qhit_len1 - $abhit_cnt > 1000){
					$judge_tp0_to_tp1 = "insertion >1kb";
				}
				elsif($qhit_len1 - $abhit_cnt > 500){
					$judge_tp0_to_tp1 = "insertion >500bp";
				}
				elsif($del_disruption_ratio ne '-'){
					if($del_disruption_ratio < 0.1){
						$judge_tp0_to_tp1 = "present (collapsed)";
					}
					elsif($del_disruption_ratio < 0.5){
						$judge_tp0_to_tp1 = "present (collapsed)";
					}
					elsif($del_disruption_ratio < 0.8){
						$judge_tp0_to_tp1 = "present (partly)";
					}
					else{
						$judge_tp0_to_tp1 = "present";
					}
				}
				else{
					$judge_tp0_to_tp1 = "present";
				}
			}
			elsif($del_disruption_ratio ne '-'){
				if($del_disruption_ratio < 0.1){
					$judge_tp0_to_tp1 = "present (collapsed)";
				}
				elsif($del_disruption_ratio < 0.5){
					$judge_tp0_to_tp1 = "present (collapsed)";
				}
				elsif($del_disruption_ratio < 0.8){
					$judge_tp0_to_tp1 = "present (partly)";
				}
				else{
					$judge_tp0_to_tp1 = "present";
				}
			}
			
			print " P1border0 -> $tmp[0] (alignment border)\n";
			print " P1border1 -> $tmp[1] (alignment border)\n";
		}
		else{
			$judge_tp0_to_tp1 .= " (missing border)";
			
			if($p1ab0_border eq 'null'){
				print " P1border0 -> $P1ab_min (alignment block)\n";
				$p1ab0_border = $P1ab_min;
			}
			if($p1ab1_border eq 'null'){
				print " P1border1 -> $P1ab_max (alignment block)\n";
				$p1ab1_border = $P1ab_max;
			}
		}
		
		my $freq_presence = 1;
		if($del_disruption_ratio ne '-' && $ins_disruption_ratio ne '-'){
			if($ins_disruption_ratio < $del_disruption_ratio){
				$freq_presence = $ins_disruption_ratio;
			}
			else{
				$freq_presence = $del_disruption_ratio;
			}
		}
		elsif($del_disruption_ratio ne '-' && $ins_disruption_ratio eq '-'){
			$freq_presence = $del_disruption_ratio;
		}
		elsif($del_disruption_ratio eq '-' && $ins_disruption_ratio ne '-'){
			$freq_presence = $ins_disruption_ratio;
		}
		
	#	my $abhit_header = "directory\tdbfasta\tseqid\tpos0\tpos1\tsp0\tsp1\tqfasta\tqfasta seqid\thit pos0\thit pos1\tgroup\thit sp0\thit sp1\tbp (sp0-sp1)\t";
	#	$abhit_header .= "bp hit-align (sp0-sp1)\tbp hit-span (sp0-sp1)\tbp hit-span woN (sp0-sp1)\tbp insertion\talign ratio (1=not disrupted)\tinsert ratio (1=not disrupted)\t";
	#	$abhit_header .= "seq normality (1=not disrupted)\tjudge\tnote\t";
	#	$abhit_header .= "5'-promoter bp hit ($neighbor_tlen bp)\t5'-promoter align ratio ($neighbor_tlen bp)\t3'-UTR bp hit ($neighbor_tlen bp)\t3'-UTR align ratio ($neighbor_tlen bp)\t";
	#	$abhit_header .= "5' flanking bp hit ($flen bp)\t5' flanking align ratio ($flen bp)\t3' flanking bp hit ($flen bp)\t3' flanking align ratio ($flen bp)\n";
		
		$abhit_summary .= $dir."\t".$dbpref."\t".$seqid."\t".$tpos0."\t".$tpos1."\t".$ab0."\t".$ab1."\t".$qpref."\t".$qseqid."\t".$qtpos0."\t".$qtpos1."\t".$grp."\t".$p1ab0_border."\t".$p1ab1_border."\t".$abdist."\t".$abhit_cnt."\t".$qhit_len0."\t".$qhit_len1."\t".$bp_insert."\t".$del_disruption_ratio."\t".$ins_disruption_ratio."\t".$freq_presence."\t".$judge_tp0_to_tp1."\t".$note."\t".$prohit_cnt."\t".$prohit_ratio."\t".$utrhit_cnt."\t".$utrhit_ratio."\t".$n5hit_cnt."\t".$n5hit_ratio."\t".$n3hit_cnt."\t".$n3hit_ratio."\n";
		
		$abhit_judge = $judge_tp0_to_tp1."\n";
		
		if($freq_presence ne '-'){
			$abhit_score .= $freq_presence."\n";
		}
		
		print " - [$note] [$mode] | $judge_tp0_to_tp1 | $freq_presence\n";
#		print "$abhit_summary\n";
#		print "$abhit_score";
	}
	else{
		$abhit_summary .= $dir."\t".$dbpref."\t".$seqid."\t".$tpos0."\t".$tpos1."\t-\t-\t".$qpref."\t".$qseqid."\t".$qtpos0."\t".$qtpos1."\t".$grp."\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t".$note."\t-\t-\t-\t-\t-\t-\t-\t-\n";
		$abhit_score .= "0\n";
		$abhit_judge = "not_found\n";
		print " - [$note] [$mode] | not_found\n";
	}
	
#	# not use---------------------------------------------//
#	my $nrlines = {};
#	@{$AoP} = sort {$a->[0] <=> $b->[0]} @{$AoP};
#	my $numAoP = @{$AoP};
#	for(my $i = 0; $i < $numAoP - 1; $i++){
#		my $j = $i + 1;
#		
#		my $qp0 = $AoP->[$i][0];
#		my $qp1 = $AoP->[$j][0];
#		my $tp0 = $AoP->[$i][1];
#		my $tp1 = $AoP->[$j][1];
#		
#		my $dq = $qp1 - $qp0;
#		my $dt = $tp1 - $tp0;
#		if($dt >= 100 && $qp0 != $qp1){
#			my $tmpl = $dir."\t".$dbpref."\t".$seqid."_".$tpos0."_".$tpos1."\t".$qp0."\t".$qp1."\t".$dq."\t".$qpref."\t".$qtpref."\t".$tp0."\t".$tp1."\t".$dt."\t".$note."\n";
#			unless($nrlines->{$tmpl}){
#				$pav_summary .= $tmpl;
#				$nrlines->{$tmpl} = 1;
#			}
#		}
#	}
#	
#	@{$AoP} = sort {$a->[1] <=> $b->[1]} @{$AoP};
#	for(my $i = 0; $i < $numAoP - 1; $i++){
#		my $j = $i + 1;
#		
#		my $qp0 = $AoP->[$i][0];
#		my $qp1 = $AoP->[$j][0];
#		my $tp0 = $AoP->[$i][1];
#		my $tp1 = $AoP->[$j][1];
#		
#		my $dq = $qp1 - $qp0;
#		my $dt = $tp1 - $tp0;
#		if($dq >= 100 && $tp0 != $tp1){
#			my $tmpl = $dir."\t".$dbpref."\t".$seqid."_".$tpos0."_".$tpos1."\t".$qp0."\t".$qp1."\t".$dq."\t".$qpref."\t".$qtpref."\t".$tp0."\t".$tp1."\t".$dt."\t".$note."\n";
#			unless($nrlines->{$tmpl}){
#				$pav_summary .= $tmpl;
#				$nrlines->{$tmpl} = 1;
#			}
#		}
#	}
#	# not use---------------------------------------------//
	
	@each_Rqpos = sort {$a <=> $b} @each_Rqpos;
	@each_Rtpos = sort {$a <=> $b} @each_Rtpos;
#	my $qpos_min = $each_Rqpos[0];
	my $qpos_min = $tpos0;
	my $tpos_min = $each_Rtpos[0];

	@each_Rqpos = sort {$b <=> $a} @each_Rqpos;
	@each_Rtpos = sort {$b <=> $a} @each_Rtpos;
#	my $qpos_max = $each_Rqpos[0];
	my $qpos_max = $tpos1;
	my $tpos_max = $each_Rtpos[0];
	
	if($nomask eq 'true'){
		$tpos_min = 1;
		$tpos_max = $hseqlen;
	}

	my $vline;
	$vline .= "abline(v\=$ab0, col\=\"black\")\n";
	$vline .= "abline(v\=$ab1, col\=\"black\")\n";
	my $hline;
	$hline .= "abline(h\=$ab0, col\=\"black\")\n";
	$hline .= "abline(h\=$ab1, col\=\"black\")\n";
	
	if($tstrand ne 'null'){
		$vline .= "abline(v\=$ab0nb, col\=\"black\", lty=3)\n";
		$vline .= "abline(v\=$ab1nb, col\=\"black\", lty=3)\n";
		$hline .= "abline(h\=$ab0nb, col\=\"black\", lty=3)\n";
		$hline .= "abline(h\=$ab1nb, col\=\"black\", lty=3)\n";
	}
	
	#------------------------------------------------------//
	my $script =<<"EOS";
fw <- read.table("$pfile_Fw", sep="\t", header=T)
rv <- read.table("$pfile_Rv", sep="\t", header=T)

png("$png1", width=2400, height=2400)
plot(fw[,1], fw[,2], pch=20, col="blue", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3, xlab="$seqid ($dbpref, $tpos0 - $tpos1)", ylab="$qtpref ($qpref)")
points(rv[,1], rv[,2], pch=20, col="red", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3)
abline(h=$tpos_min, col="black")
abline(h=$tpos_max, col="black")
abline(v=$qpos_min, col="black")
abline(v=$qpos_max, col="black")
$vline
dev.off()

png("$png2", width=2400, height=2400)
plot(fw[,1], fw[,2], pch=20, col="red", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3, xlab="$seqid ($dbpref, $tpos0 - $tpos1)", ylab="$qtpref ($qpref)")
points(rv[,1], rv[,2], pch=20, col="blue", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3)
abline(h=$tpos_min, col="black")
abline(h=$tpos_max, col="black")
abline(v=$qpos_min, col="black")
abline(v=$qpos_max, col="black")
$vline
dev.off()

png("$png3", width=2400, height=2400)
plot(fw[,2], fw[,1], pch=20, col="blue", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3, xlab="$qtpref ($qpref)", ylab="$seqid ($dbpref, $tpos0 - $tpos1)")
points(rv[,2], rv[,1], pch=20, col="red", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3)
abline(v=$tpos_min, col="black")
abline(v=$tpos_max, col="black")
abline(h=$qpos_min, col="black")
abline(h=$qpos_max, col="black")
$hline
dev.off()

png("$png4", width=2400, height=2400)
plot(fw[,2], fw[,1], pch=20, col="red", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3, xlab="$qtpref ($qpref)", ylab="$seqid ($dbpref, $tpos0 - $tpos1)")
points(rv[,2], rv[,1], pch=20, col="blue", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3)
abline(v=$tpos_min, col="black")
abline(v=$tpos_max, col="black")
abline(h=$qpos_min, col="black")
abline(h=$qpos_max, col="black")
$hline
dev.off()

EOS

	my $R_script = "R-plot_$prefix.$grp.txt";
	if($outplot eq '1'){
		open(my $sfh, ">", $R_script) or die;
		print $sfh $script;
		close $sfh;
		#print "! output [$R_script] ...\n";

		#my $sh_script = "R-WGCNA-clust.sh";
		#open(my $shfh, ">", $sh_script) or die;
		#print $shfh "Rscript --vanilla --slave $R_script";
		#close $shfh;

		my $Rcmd = "Rscript --vanilla --slave $R_script >/dev/null 2>&1";
		#print "! cmd=[$Rcmd]\n";
		system("$Rcmd");

		if(-e "Rplots.pdf"){
			system("rm Rplots.pdf");
		}
	}
	
	if($keep eq '0'){
		if(-e $pfile_Fw){
			system("rm $pfile_Fw");
		}
		if(-e $pfile_Rv){
			system("rm $pfile_Rv");
		}
		if(-e $R_script){
			system("rm $R_script");
		}
	}
	
	if($allplot eq '0'){
		last;
	}
}

$Rrh->{abhit_summary} = $abhit_summary;
$Rrh->{pav_summary} = $pav_summary;
$Rrh->{abhit_score} = $abhit_score;
$Rrh->{abhit_judge} = $abhit_judge;

return $Rrh;
}


#-------------------------------------------------------------------------------
sub Open_fasta_as_hash{
my $file = shift;

print "! reading [$file] as hash ...\n";
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
			
			my $len = length($seq);
			$hash->{$ID}{seq} = $seq;
			$hash->{$ID}{len} = $len;
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
		print "! duplicated seq ID : [$ID]\n";
	}
	
	my $len = length($seq);
	$hash->{$ID}{seq} = $seq;
	$hash->{$ID}{len} = $len;
	$total_len += $len;
	push(@SID, $ID);
}

my $numID = @SID;
print "! [$numID] sequences with [$total_len] bp\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub LASTINDEX{
my $fasta = shift;

my $reference = $fasta;
#$reference =~ s/\.fasta//;
#$reference =~ s/\.fa//;

print "! creating lastal index for [$reference]\n";
system("lastdb $reference $fasta");

}


#-------------------------------------------------------------------------------
sub Read_blast2{
my $qid = shift;
my $file = shift;
my $mth = shift;
my $wth = shift;
my $lth = shift;
my $tpos0 = shift;
my $tpos1 = shift;

open(my $fh, "<", $file) or die;
my $hash = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	#chr01_24000000-24050000 tig00000060     97.276  50030   464     756     1       49999   3145499 3096338 0.0     84007
	my @A = split(/\t/, $line);
	if($A[1] !~ /chr00/i && $A[2] >= $mth && $A[3] >= $lth){
		my $p0 = $A[6];
		my $p1 = $A[7];
		my $t0 = $A[8];
		my $t1 = $A[9];
		
		if($p0 > $p1){
			$A[6] = $p1;
			$A[7] = $p0;
		}
		if($t0 > $t1){
			$A[8] = $t1;
			$A[9] = $t0;
		}
		
		$hash->{$A[1]} .= join(",", @A)."\n";
	}
}
close $fh;

my @HID = keys(%{$hash});
my $AoB = [];
foreach my $hid (@HID){
	my $AoA = [];
	my @D = split(/\n/, $hash->{$hid});
	foreach my $line (@D){
		my @A = split(/,/, $line);
		push(@{$AoA}, \@A);
	}
	@{$AoA} = sort {$b->[11] <=> $a->[11]} @{$AoA};		#sort by score
	
	my $sl = {};
	my $AoA2 = [];
	foreach my $A (@{$AoA}){
		my $p0 = $A->[6];
		my $p1 = $A->[7];
		
		my $cnt_prevhit = 0;
		for(my $i = $p0; $i <= $p1; $i++){
			unless($sl->{$i}){
				$sl->{$i} = 1;
			}
			else{
				$cnt_prevhit++;
			}
		}
		
		if($cnt_prevhit < abs($A->[9] - $A->[8]) / 2){
			push(@{$AoA2}, $A);
		}
	}
	@{$AoA2} = sort {$a->[8] <=> $b->[8]} @{$AoA2};		#sort by position
	
	my $n = @{$AoA2};
	my $p0 = $AoA2->[0][6];
	my $p1 = $AoA2->[0][7];
	my $t0 = $AoA2->[0][8];
	my $t1 = $AoA2->[0][9];
	my $alen = $AoA2->[0][3];
	my $pcnt = $AoA2->[0][2];
	my $score = $AoA2->[0][11];
	my $subn = 1;
	if($n > 1){
		for(my $i = 0; $i < $n - 1; $i++){
			my $j = $i + 1;
			my $A1 = $AoA2->[$i];
			my $A2 = $AoA2->[$j];
			
			my $d = $A2->[8] - $A1->[9];
			if($d <= $wth){
#				print " [$i] : [$A1->[6] - $A1->[7]] | [$A1->[8] - $A1->[9]] [$A1->[2]]\n";
#				print " [$j] : [$A2->[6] - $A2->[7]] | [$A2->[8] - $A2->[9]] [$A2->[2]] connected\n\n";
				
				$alen += $A2->[3];
				$pcnt += $A2->[2];
				$score += $A2->[11];
				$subn += 1;
				
				if($p1 < $A2->[7]){
					$p1 = $A2->[7];
				}
				if($t1 < $A2->[9]){
					$t1 = $A2->[9];
				}
			}
			else{
#				print " [$i] : [$A1->[6] - $A1->[7]] | [$A1->[8] - $A1->[9]] [$A1->[2]]\n";
#				print " [$j] : [$A2->[6] - $A2->[7]] | [$A2->[8] - $A2->[9]] [$A2->[2]] split\n\n";
				
				if($p1 < $A1->[7]){
					$p1 = $A1->[7];
				}
				if($t1 < $A1->[9]){
					$t1 = $A1->[9];
				}
				
				$pcnt = $pcnt / $subn;
				my $point = $score / (101 - $pcnt);
				my @B;
				push(@B, $hid);
				push(@B, $t0);
				push(@B, $t1);
				push(@B, $alen);
				push(@B, $p0);
				push(@B, $p1);
				push(@B, $score);
				push(@B, $pcnt);
				push(@B, $subn);
				push(@B, $point);
				push(@{$AoB}, \@B);
				
				$alen = $A2->[3];
				$pcnt = $A2->[2];
				$score = $A2->[11];
				$subn = 1;
				$p0 = $A2->[6];
				$p1 = $A2->[7];
				$t0 = $A2->[8];
				$t1 = $A2->[9];
			}
		}
	}
	
	$pcnt = $pcnt / $subn;
	my $point = $score / (101 - $pcnt);
	my @B;
	push(@B, $hid);
	push(@B, $t0);
	push(@B, $t1);
	push(@B, $alen);
	push(@B, $p0);
	push(@B, $p1);
	push(@B, $score);
	push(@B, $pcnt);
	push(@B, $subn);
	push(@B, $point);
	push(@{$AoB}, \@B);
}

@{$AoB} = sort {$b->[9] <=> $a->[9]} @{$AoB};

my $tdist = abs($tpos1 - $tpos0);
my $cnt = 1;
my $AoB2 = [];
foreach my $B (@{$AoB}){
	if($cnt == 1){
		push(@{$AoB2}, $B);
		print " [$cnt] : [$B->[0]] [$B->[1]] - [$B->[2]] : [$B->[3]] bp, score = [$B->[6]], point = [$B->[9]]\n";
		$tdist -= $B->[3];
	}
	else{
		if($B->[3] <= ($tdist * 5)){
			push(@{$AoB2}, $B);
			print " [$cnt] : [$B->[0]] [$B->[1]] - [$B->[2]] : [$B->[3]] bp, score = [$B->[6]], point = [$B->[9]]\n";
		}
		else{
			print " [omit] : [$B->[0]] [$B->[1]] - [$B->[2]] : [$B->[3]] bp, score = [$B->[6]], point = [$B->[9]]\n";
		}
	}
	
	if($cnt == 10){
		last;
	}
	$cnt++;
}

return $AoB2;
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
sub Search_seq{
my $file = shift;
my $seqid = shift;
my $tpos0 = shift;
my $tpos1 = shift;
my $rfile = shift;
my $mask0 = shift;
my $mask1 = shift;
my $append = shift;

unless($append){
	$append = "false";
}

#print "! reading [$file]...\n";
open(my $fh, "<", $file) or die;
my $ID;
my $seq;
my $tseq;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ />/){
		if($seq && $ID eq $seqid){
			$tseq = $seq;
			last;
		}
		
		$ID = $line;
		$ID =~ s/>//;
		my @B = split(/\s/, $ID);
		$ID = $B[0];
		$seq = "";
	}
	else{
		$seq .= $line;
	}
}
close $fh;

if($seq && $ID eq $seqid){
	$tseq = $seq;
}

my $len = 0;
unless($tseq){
	print "! sequence not found for [$seqid]\n";
	return;
}
else{
	$len = length($tseq);
	
	if($tpos0 ne 'NA' && $tpos1 ne 'NA'){
		my $d = abs($tpos1 - $tpos0) + 1;
		$tseq = ">".$seqid."_".$tpos0."-".$tpos1."\n".substr($tseq, $tpos0 - 1, $d);
	}
	elsif($mask0 ne 'NA' && $mask1 ne 'NA'){
		my $s0 = substr($tseq, 0, $mask0 - 1);
		my $s1 = substr($tseq, $mask0 - 1, $mask1 - $mask0 + 1);
		my $s2 = substr($tseq, $mask1);
		$s0 =~ s/[a-z]/n/gi;
		$s2 =~ s/[a-z]/n/gi;
		$tseq = ">".$seqid."\n".$s0.$s1.$s2;
	}
	else{
		$tseq = ">".$seqid."\n".$tseq;
	}
}

if($tseq){
	if($append eq 'true'){
		open(my $rfh, ">>", $rfile);
		print $rfh $tseq;
		close $rfh;
	}
	else{
		open(my $rfh, ">", $rfile);
		print $rfh $tseq;
		close $rfh;
		print "! output [$rfile]\n";
	}
}

return $len;
}


#----------------------------------------------------------
sub Search_targetpos{
my $file = shift;
my $lth = shift;
my $seqid = shift;
my $tpos0 = shift;
my $tpos1 = shift;

print "! reading psl [$file]...\n";
open(my $fh, "<", $file) or die;
my @PSLdata;
my $hnrblock = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	$line =~ s/unanchored_scaffolds/chr00/;
	push(@PSLdata, $line);
	
	my @A = split(/\t/, $line);
	
	my $matches = $A[0];
	my $misMatches = $A[1];
	my $repMatches = $A[2];
	my $nCount = $A[3];
	my $qNumInsert = $A[4];
	my $qBaseInsert = $A[5];
	my $tNumInsert = $A[6];
	my $tBaseInsert = $A[7];
	my $strand = $A[8];
	my $qName = $A[9];
	my $qSize = $A[10];
	my $qStart = $A[11];
	my $qEnd = $A[12];
	my $tName = $A[13];
	my $tSize = $A[14];
	my $tStart = $A[15];
	my $tEnd = $A[16];
	my $blockCount = $A[17];
	my $blockSizes = $A[18];
	my $qStarts = $A[19];
	my $tStarts = $A[20];
	
	if($matches >= $lth){
		my @qpos = split(/,/, $qStarts);
		my @tpos = split(/,/, $tStarts);
		my @blks = split(/,/, $blockSizes);
		my $npart = @blks;
		
		my @nqpos;
		my @ntpos;
		for(my $i = 0; $i < $npart; $i++){
			my $q0 = $qpos[$i];
			my $t0 = $tpos[$i];
			my $b0 = $blks[$i];
			
			if($strand eq '-'){
				$q0 = $qSize - $q0 - $b0;
			}
			
			for(my $j = 0; $j < $b0; $j++){
				my $nq0 = $q0 + $j;
				my $nt0 = $t0 + $j;
				
				if(! defined $hnrblock->{$qName}{$tName}{$nt0}{bs}){
					$hnrblock->{$qName}{$tName}{$nt0}{bs} = $b0;
					$hnrblock->{$qName}{$tName}{$nt0}{qpos} = $nq0;
					$hnrblock->{$qName}{$tName}{$nt0}{qStart} = $q0;
					$hnrblock->{$qName}{$tName}{$nt0}{qjudge}{$q0} = 1;
#					print "define | $nt0 | $q0\n";
				}
				elsif($b0 > $hnrblock->{$qName}{$tName}{$nt0}{bs}){			# select alignment with max blocksize (remove overlapping)
					if($hnrblock->{$qName}{$tName}{$nt0}{qStart} && $hnrblock->{$qName}{$tName}{$nt0}{qStart} ne $q0){
						my $prev_q0 = $hnrblock->{$qName}{$tName}{$nt0}{qStart};
						$hnrblock->{$qName}{$tName}{$nt0}{qjudge}{$prev_q0} = 0;
#						print "overlap | $nt0 | $prev_q0 -> 0\n";
					}
					
					$hnrblock->{$qName}{$tName}{$nt0}{bs} = $b0;
					$hnrblock->{$qName}{$tName}{$nt0}{qpos} = $nq0;
					$hnrblock->{$qName}{$tName}{$nt0}{qStart} = $q0;
					$hnrblock->{$qName}{$tName}{$nt0}{qjudge}{$q0} = 1;
#					print "replace | $nt0 | $q0 -> 1\n";
				}
				else{
					$hnrblock->{$qName}{$tName}{$nt0}{qjudge}{$q0} = 0;
#					print "remove | $nt0 | $q0\n";
				}
			}
		}
	}
}
close $fh;

my $hstat = {};
my $hqsize = {};
my $htsize = {};
my $hpos = {};
my $hstrand = {};
my $hmatches = {};
my $hmisMatches = {};
my $hblockCount = {};
my $cnt = 0;
my $cnt_selected = 0;
my $bp_overlap = 0;
foreach my $line (@PSLdata){
	my @A = split(/\t/, $line);
	
	my $matches = $A[0];
	my $misMatches = $A[1];
	my $repMatches = $A[2];
	my $nCount = $A[3];
	my $qNumInsert = $A[4];
	my $qBaseInsert = $A[5];
	my $tNumInsert = $A[6];
	my $tBaseInsert = $A[7];
	my $strand = $A[8];
	my $qName = $A[9];
	my $qSize = $A[10];
	my $qStart = $A[11];
	my $qEnd = $A[12];
	my $tName = $A[13];
	my $tSize = $A[14];
	my $tStart = $A[15];
	my $tEnd = $A[16];
	my $blockCount = $A[17];
	my $blockSizes = $A[18];
	my $qStarts = $A[19];
	my $tStarts = $A[20];
	
	$hstat->{sum_matches}{$qName}{$tName} += $matches;
	$hstat->{sum_misMatches}{$qName}{$tName} += $misMatches;
	$hstat->{sum_repMatches}{$qName}{$tName} += $repMatches;
	$hstat->{sum_nCount}{$qName}{$tName} += $nCount;
	$hstat->{sum_blockCount}{$qName}{$tName} += $blockCount;
	$hstat->{sum_blockSizes}{$qName}{$tName} .= $blockSizes;
	
	$hstat->{sum_matches}{total} += $matches;
	$hstat->{sum_misMatches}{total} += $misMatches;
	$hstat->{sum_repMatches}{total} += $repMatches;
	$hstat->{sum_nCount}{total} += $nCount;
	$hstat->{sum_blockCount}{total} += $blockCount;
	$hstat->{sum_blockSizes}{total} .= $blockSizes;
	
	$hqsize->{$qName} = $qSize;
	$htsize->{$tName} = $tSize;
	$hmatches->{$qName}{$tName} += $matches;
	$hmisMatches->{$qName}{$tName} += $misMatches;
	$hblockCount->{$qName}{$tName} += $blockCount;
	
	if($matches >= $lth){
		my @qpos = split(/,/, $qStarts);
		my @tpos = split(/,/, $tStarts);
		my @blks = split(/,/, $blockSizes);
		my $npart = @blks;
		
		my @nqpos;
		my @ntpos;
		for(my $i = 0; $i < $npart; $i++){
			my $q0 = $qpos[$i];
			my $t0 = $tpos[$i];
			my $b0 = $blks[$i];
			
			if($strand eq '-'){
				$q0 = $qSize - $q0 - $b0;
			}
			
			for(my $j = 0; $j < $b0; $j++){
				my $nq0 = $q0 + $j;
				my $nt0 = $t0 + $j;
				
				if(defined $hnrblock->{$qName}{$tName}{$nt0}{qpos} && $hnrblock->{$qName}{$tName}{$nt0}{qpos} eq $nq0 && defined $hnrblock->{$qName}{$tName}{$nt0}{qjudge}{$q0} && $hnrblock->{$qName}{$tName}{$nt0}{qjudge}{$q0} == 1){			# best possible position
					push(@nqpos, $nq0);
					push(@ntpos, $nt0);
				}
				else{
					$bp_overlap++;
				}
			}
		}
		
		$qStarts = join(",", @nqpos);
		$tStarts = join(",", @ntpos);
		
		$hpos->{$qName}{$tName} .= $qStarts."\t".$tStarts."\n";
		$hstrand->{$qName}{$tName} .= $strand."\n";
		$cnt_selected++;
	}
	
	$cnt++;
}

if($cnt == 0){
	print "! no data found for [$seqid]\n";
	my $rhnull = {};
	$rhnull->{data} = "false";
	return $rhnull;
}

my @QID = keys(%{$hqsize});
my @TID = keys(%{$htsize});
@QID = sort {$a cmp $b} @QID;
@TID = sort {$a cmp $b} @TID;
my $numQ = @QID;
my $numT = @TID;

if($QID[0] =~ /00/){
	my $chr00 = shift(@QID);
	push(@QID, $chr00);
}
if($TID[0] =~ /00/){
	my $chr00 = shift(@TID);
	push(@TID, $chr00);
}

my $stats_log;
$stats_log .= "! total [$cnt] lines with \'chr\' ID in [$file]\n";
$stats_log .= "! [$hstat->{sum_matches}{total}] = total number of matching bases that aren't repeats\n";
$stats_log .= "! [$hstat->{sum_misMatches}{total}] = total number of bases that don't match\n";
$stats_log .= "! [$hstat->{sum_repMatches}{total}] = total number of matching bases that are part of repeats\n";
$stats_log .= "! [$hstat->{sum_nCount}{total}] = total number of 'N' bases\n";
$stats_log .= "! [$hstat->{sum_blockCount}{total}] = total number of blocks in the alignment\n";
$stats_log .= "! [$cnt_selected] lines with matches > [$lth] bp\n";
$stats_log .= "! [$numQ] = total number of seq ID in query\n";
$stats_log .= "! [$numT] = total number of seq ID in target\n";
$stats_log .= "! [$bp_overlap] = overlapping bases (removed)\n";

if($numQ < 30 && $numT < 30){
	$stats_log .= "! seq IDs in query  = [".join(",", @QID)."]\n";
	$stats_log .= "! seq IDs in target = [".join(",", @TID)."]\n";
}
#$stats_log .= "\n";

print "! aligned to [$numT] target sequences\n";
print "! [$hstat->{sum_matches}{total}] bp = total number of matching bases that aren't repeats\n";
print "! [$bp_overlap] = overlapping bases (removed)\n";

my $hqcumsize = {};
my $htcumsize = {};

my $stats_pos = "type\tseqID\tlength\tcumlative length\n";
my $qcumpos = 0;
foreach my $qid (@QID){
	$stats_pos .= "query\t".$qid."\t".$hqsize->{$qid}."\t".$qcumpos."\n";
	$hqcumsize->{$qid} = $qcumpos;
	$qcumpos += $hqsize->{$qid};
}
$stats_pos .= "query\t-\t-\t".$qcumpos."\n";

my $tcumpos = 0;
foreach my $tid (@TID){
	$stats_pos .= "target\t".$tid."\t".$htsize->{$tid}."\t".$tcumpos."\n";
	$htcumsize->{$tid} = $tcumpos;
	$tcumpos += $htsize->{$tid};
}
$stats_pos .= "target\t-\t-\t".$tcumpos."\n";

#print "! converting value for position...\n";
my @qposFw;
my @qposRv;
my @tposFw;
my @tposRv;
my $hposdata = {};
my $stats_count = "query,target,num block,num block (Fw),num block (Rv),num alignment\n";
my $stats_hmatches .= "query | target->\t".join("\t", @TID)."\n";
my $stats_hmisMatches .= "query | target->\t".join("\t", @TID)."\n";
my $stats_hblockCount .= "query | target->\t".join("\t", @TID)."\n";
for(my $q = 0; $q < $numQ; $q++){
	my $qid = $QID[$q];
	
	my @val0;
	my @val1;
	my @val2;
	for(my $t = 0; $t < $numT; $t++){
		my $tid = $TID[$t];
		
		if($hpos->{$qid}{$tid}){
			my @L = split(/\n/, $hpos->{$qid}{$tid});
			my @S = split(/\n/, $hstrand->{$qid}{$tid});
			my $numL = @L;
			my $numS = @S;
			
			unless($numL == $numS){
				print "! count unmatched between start_position and strand info for [$qid] [$tid]\n";
				die;
			}
			
			my $num_qpos = 0;
			my $num_tpos = 0;
			my $num_fwaln = 0;
			my $num_rvaln = 0;
			
			my @qpos_sub;
			my @qposFw_sub;
			my @qposRv_sub;
			my @tpos_sub;
			my @tposFw_sub;
			my @tposRv_sub;
			
			for(my $i = 0; $i < $numL; $i++){
				my @subL = split(/\t/, $L[$i]);
				my @qpos = split(/,/, $subL[0]);
				my @tpos = split(/,/, $subL[1]);
				my $strand = $S[$i];
				
				$num_qpos += @qpos;
				$num_tpos += @tpos;
				
				if($strand eq '+'){
					$num_fwaln += 1;
				}
				elsif($strand eq '-'){
					$num_rvaln += 1;
				}
				
				unless($num_qpos == $num_tpos){
					print "! total count of position data unmatched between [$qid $num_qpos] and [$tid $num_tpos]\n";
					die;
				}
				
				foreach my $p (@qpos){
					if($strand eq '+'){
						push(@qposFw_sub, $p);
					}
					elsif($strand eq '-'){
						push(@qposRv_sub, $p);
					}
					push(@qpos_sub, $p);
					
					$p += $hqcumsize->{$qid};
					
					if($strand eq '+'){
						push(@qposFw, $p);
					}
					elsif($strand eq '-'){
						push(@qposRv, $p);
					}
				}
				foreach my $p (@tpos){
					if($strand eq '+'){
						push(@tposFw_sub, $p);
					}
					elsif($strand eq '-'){
						push(@tposRv_sub, $p);
					}
					push(@tpos_sub, $p);
					
					$p += $htcumsize->{$tid};
					
					if($strand eq '+'){
						push(@tposFw, $p);
					}
					elsif($strand eq '-'){
						push(@tposRv, $p);
					}
				}
			}
			
			$hposdata->{$tid}{qposFw} = \@qposFw_sub;
			$hposdata->{$tid}{qposRv} = \@qposRv_sub;
			$hposdata->{$tid}{tposFw} = \@tposFw_sub;
			$hposdata->{$tid}{tposRv} = \@tposRv_sub;
			$hposdata->{$tid}{qcumpos} = $hqsize->{$qid};
			$hposdata->{$tid}{tcumpos} = $htsize->{$tid};
			$hposdata->{$tid}{hmatches} = $hmatches->{$qid}{$tid};
			$hposdata->{$tid}{stats_pos} .= "query\t-\t-\t0\n";
			$hposdata->{$tid}{stats_pos} .= "query\t-\t-\t".$hqsize->{$qid}."\n";
			$hposdata->{$tid}{stats_pos} .= "target\t-\t-\t0\n";
			$hposdata->{$tid}{stats_pos} .= "target\t-\t-\t".$htsize->{$tid};
			$hposdata->{$tid}{seqlen} = $htsize->{$tid};
			
			@qpos_sub = sort {$a <=> $b} @qpos_sub;
			@tpos_sub = sort {$a <=> $b} @tpos_sub;
			$hposdata->{$tid}{qpos_min} = $qpos_sub[0];
			$hposdata->{$tid}{tpos_min} = $tpos_sub[0];
			
			@qpos_sub = sort {$b <=> $a} @qpos_sub;
			@tpos_sub = sort {$b <=> $a} @tpos_sub;
			$hposdata->{$tid}{qpos_max} = $qpos_sub[0];
			$hposdata->{$tid}{tpos_max} = $tpos_sub[0];
			
#			print "! [$qid] [$tid] = [$numL] blocks (Fw=[$num_fwaln], Rv=[$num_rvaln]) with [$num_qpos] alignments\n";
			$stats_count .= $qid."\t".$tid."\t".$numL."\t".$num_fwaln."\t".$num_rvaln."\t".$num_qpos."\n";
			
			push(@val0, $hmatches->{$qid}{$tid});
			push(@val1, $hmisMatches->{$qid}{$tid});
			push(@val2, $hblockCount->{$qid}{$tid});
		}
	}
	
	$stats_hmatches .= $qid."\t".join("\t", @val0)."\n";
	$stats_hmisMatches .= $qid."\t".join("\t", @val1)."\n";
	$stats_hblockCount .= $qid."\t".join("\t", @val2)."\n";
}

$hposdata->{All}{qposFw} = \@qposFw;
$hposdata->{All}{qposRv} = \@qposRv;
$hposdata->{All}{tposFw} = \@tposFw;
$hposdata->{All}{tposRv} = \@tposRv;
$hposdata->{All}{qcumpos} = $qcumpos;
$hposdata->{All}{tcumpos} = $tcumpos;
$hposdata->{All}{stats_pos} = $stats_pos;

my $rh = {};
$rh->{QID} = \@QID;
$rh->{TID} = \@TID;
$rh->{hposdata} = $hposdata;
$rh->{stats_log} = $stats_log;
$rh->{stats_count} = $stats_count;
$rh->{stats_pos} = $stats_pos;
$rh->{stats_hmatches} = $stats_hmatches;
$rh->{stats_hmisMatches} = $stats_hmisMatches;
$rh->{stats_hblockCount} = $stats_hblockCount;
$rh->{sum_matches_total} = $hstat->{sum_matches}{total};
$rh->{data} = "true";

return $rh;
}


#----------------------------------------------------------
sub SAVE{
my $file = shift;
my $str = shift;

if($str){
	open(my $fh, ">", $file) or die;
	print $fh $str;
	close $fh;
}

#print "! output [$file]\n";

}


#----------------------------------------------------------
sub APPEND{
my $file = shift;
my $str = shift;

if($str){
	open(my $fh, ">>", $file) or die;
	print $fh $str;
	close $fh;
}

#print "! output [$file]\n";

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
		unless($file =~ /\.fai/ || $file =~ /\.fa\.n../ || $file =~ /\.fa\.p../ || $file =~ /\.fasta\.n../ || $file =~ /\.fasta\.p../ || $file =~ /\.ffa\.n../ || $file =~ /\.csv/ || $file =~ /_BLK/){
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

