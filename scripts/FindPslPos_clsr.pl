#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;

my $refgenome = shift;
my $dbgff = shift;
my $qpref = shift;
my $qfasta = shift;
my $qpsl = shift;
my $genelist = shift;
my $chrinfo = shift;
my $afile = shift;
my $ffile = shift;
my $rfile_hits0 = shift;		# "r1_hitalign_flk-".$kbflnk."kb.tsv";
my $rfile_hits1 = shift;		# "r1_hitsummary_flk-".$kbflnk."kb.tsv";
my $outplot_dir = shift;		# "plot";
my $dir = shift;
my $wpath = shift;
my $initspan = shift;
my $neighbor_tlen = shift;
my $p0 = shift;
my $p1 = shift;
my $np = shift;
my $outplot = shift;

if(! $refgenome || ! $wpath){
	goto END;
}

my $hID = Read_list($genelist);
my $hgffinfo = Read_gff($dbgff, $hID);
my $QL = $hgffinfo->{RL};
my $hlist = $hgffinfo->{hash};
my $dbg2t = $hgffinfo->{g2t};
my $dbt2g = $hgffinfo->{t2g};
my $numQL = @{$QL};
my $hqseq = Open_fasta_as_hash($qfasta);

my $ho = {};
$ho->{hsample}{0} = $refgenome;
$ho->{hsample}{1} = $qpref;
$ho->{name}{$refgenome} = 1;
$ho->{name}{$qpref} = 1;

my $hseqinfo = {};
if($chrinfo && $chrinfo ne 'null' && -e $chrinfo){
	$hseqinfo = Read_seqinfo($chrinfo);
}

my @ALNF;
$ALNF[0] = $qpsl;

my $htmp = FindPslPos_rev(\@ALNF, $QL, $p0, $p1, $np, $initspan, $neighbor_tlen, $hseqinfo, $ho, 0, $refgenome, $hqseq, $outplot, $outplot_dir, $dir, $wpath, $qpref);

my $header_summary = "directory\tdbfasta\tseqid\tpos0\tpos1\tsp0\tsp1\tqfasta\tqfasta seqid\thit pos0\thit pos1\tgroup\thit sp0\thit sp1\tbp (sp0-sp1)\tbp hit-align (sp0-sp1)\tbp hit-span (sp0-sp1)\tbp hit-span woN (sp0-sp1)\tbp insertion\talign ratio (1=not disrupted)\tinsert ratio (1=not disrupted)\tseq normality (1=not disrupted)\tjudge\tgene_id\t";
$header_summary .= "5'-promoter bp hit ($neighbor_tlen bp)\t5'-promoter align ratio ($neighbor_tlen bp)\t3'-UTR bp hit ($neighbor_tlen bp)\t3'-UTR align ratio ($neighbor_tlen bp)\t";
$header_summary .= "5' flanking bp hit ($initspan bp)\t5' flanking align ratio ($initspan bp)\t3' flanking bp hit ($initspan bp)\t3' flanking align ratio ($initspan bp)";

my $hits_header = "gene_id\tbase_ref\tchr\tpos0\tpos1\tflanking\tclass\tsearch_mode\tnmatches (bp)\taln_strand\thit_pos0\thit_pos1\tq_name\tq_seqID\tq_seqlen\tq_apos0\tq_apos1\tdb_name\tdb_seqID\tdb_seqlen\tdb_apos0\tdb_apos1\tblockCount\tblockSizes\tqStarts\tdbStarts\tlimit_qpos0\tlimit_qpos1\n";
my $summary_header = "gene_id\tbase_ref\tchr\tpos0\tpos1\tflanking\tclass\tqname\thitbp (ref-base)\thitbp (+flanking, ref-base)\tq-hitseq\thit_pos0\thit_pos1\tlimit_qpos0\tlimit_qpos1\tnum_candidate_chrs\n";
my $summary2_header = $header_summary."\n";

SAVE($rfile_hits0, $htmp->{hit}, $hits_header);
#SAVE($rfile_hits1, $htmp->{summary}, $summary_header);
SAVE($afile, $htmp->{summary2}, $summary2_header);
SAVE($ffile, $htmp->{failed_list}, "");

END:{
	my $end = 1;
}


######################################################################
#---------------------------------------------------------------------
sub SAVE{
my ($file, $str, $header) = @_;

if(! -e $file){
	if($str && $header){
		$str = $header.$str;
		open(my $rfh, ">", $file) or die;
		print $rfh $str;
		close $rfh;
	}
	elsif($str){
		open(my $rfh, ">", $file) or die;
		print $rfh $str;
		close $rfh;
	}
}
else{
	if($str){
		open(my $rfh, ">>", $file) or die;
		print $rfh $str;
		close $rfh;
	}
}

}


#---------------------------------------------------------------------
sub FindPslPos_rev{
my $ALNF = shift;
my $QL = shift;
my $p0 = shift;
my $p1 = shift;
my $np = shift;
my $initspan = shift;
my $neighbor_tlen = shift;
my $hseqinfo = shift;
my $ho = shift;
my $t = shift;
my $refgenome = shift;
my $hqseq = shift;
my $outplot = shift;
my $outplot_dir = shift;
my $dir = shift;
my $wpath = shift;
my $qpref = shift;

my $hsample = $ho->{hsample};

print " partition [$t] : search for [$np] entries...\n";
my $hash = {};
foreach my $file (@{$ALNF}){
	my @F0 = split(/q-/, $file);
	my @F1 = split(/_DB-/, $F0[1]);
	my $qname0 = $F1[0];
	my $dbname0 = $F1[1];
	$dbname0 =~ s/\.psl//;
	
	if($ho->{altname}{$qname0}){
		$qname0 = $ho->{altname}{$qname0};
	}
	if($ho->{altname}{$dbname0}){
		$dbname0 = $ho->{altname}{$dbname0};
	}
	
	if($ho->{name}{$qname0} && $ho->{name}{$dbname0}){
		open(my $fh, "<", $file);
		while(my $line = <$fh>){
			$line =~ s/\n//;
			$line =~ s/\r//;
			$hash->{$file} .= $line."\n";
		}
		close $fh;
	}
}

my $hit = "";
my $summary = "";
my $summary2 = "";
my $failed_list = "";
my $cntq = 0;
my $cnt_hit = 0;
my $cnt_selected = 0;
my $seg = 100;
my $add = $seg;
for(my $p = $p0; $p < $p1; $p++){
	if(! $QL->[$p]){
		next;
	}
	
	my @Q = split(/,/, $QL->[$p]);
	my @Ord = keys(%{$hsample});
	@Ord = sort {$a <=> $b} @Ord;
	my $numOrd = @Ord;
	
	my $hhitinfo = {};
	foreach my $norder (@Ord){
		if($norder == 0){
			next;
		}
		my $j = $norder;
		my $i = 0;
		
		my $psl = "";
		my $qname = "";
		my $dbname = "";
		my $mode = "false";
		
		if($hsample->{$i} && $hsample->{$j}){
			foreach my $file (@{$ALNF}){
				my @F0 = split(/q-/, $file);
				my @F1 = split(/_DB-/, $F0[1]);
				my $qname0 = $F1[0];
				my $dbname0 = $F1[1];
				$dbname0 =~ s/\.psl//;
				
				if($dbname0 eq $hsample->{$i} && $qname0 eq $hsample->{$j}){
					$psl = $file;
					$qname = $hsample->{$j};
					$dbname = $hsample->{$i};
					$mode = "fw";
					last;
				}
				else{
					if($ho->{altname}{$qname0}){
						$qname0 = $ho->{altname}{$qname0};
					}
					if($ho->{altname}{$dbname0}){
						$dbname0 = $ho->{altname}{$dbname0};
					}
					
					if($dbname0 eq $hsample->{$i} && $qname0 eq $hsample->{$j}){
						$psl = $file;
						$qname = $hsample->{$j};
						$dbname = $hsample->{$i};
						$mode = "fw";
						last;
					}
				}
			}
		}
		
		if($mode eq 'false'){
			print "! error : incompatible file name...\n";
			die;
		}
		
	#	print "! searching for [$psl]...\n";
		my @DL = split(/\n/, $hash->{$psl});
		my $sw = 0;
		foreach my $line (@DL){
			$line =~ s/\n//;
			$line =~ s/\r//;
			
	#		1. matches - Number of matching bases that aren't repeats.
	#		2. misMatches - Number of bases that don't match.
	#		3. repMatches - Number of matching bases that are part of repeats.
	#		4. nCount - Number of 'N' bases.
	#		5. qNumInsert - Number of inserts in query.
	#		6. qBaseInsert - Number of bases inserted into query.
	#		7. tNumInsert - Number of inserts in target.
	#		8. tBaseInsert - Number of bases inserted into target.
	#		9. strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
	#		10. qName - Query sequence name.
	#		11. qSize - Query sequence size.
	#		12. qStart - Alignment start position in query.
	#		13. qEnd - Alignment end position in query.
	#		14. tName - Target sequence name.
	#		15. tSize - Target sequence size.
	#		16. tStart - Alignment start position in target.
	#		17. tEnd - Alignment end position in target.
	#		18. blockCount - Number of blocks in the alignment.
	#		19. blockSizes - Comma-separated list of sizes of each block.
	#		20. qStarts - Comma-separated list of start position of each block in query.
	#		21. tStarts - Comma-separated list of start position of each block in target.
			
			my @L = split(/\t/, $line);
			if($L[9] =~ /scaffold/ || $L[13] =~ /scaffold/){
				next;
			}
			
			if($mode eq 'fw'){
				my $seqid9 = $L[9];
				my $seqid13 = $L[13];
				
				if($hseqinfo->{$qname}{$L[9]}){
					$seqid9 = $hseqinfo->{$qname}{$L[9]}." | ".$L[9];
#					$L[9] = $hseqinfo->{$qname}{$L[9]};
				}
				if($hseqinfo->{$dbname}{$L[13]}){
					$seqid13 = $hseqinfo->{$dbname}{$L[13]}." | ".$L[13];
#					$L[13] = $hseqinfo->{$dbname}{$L[13]};
				}
				
				my $qseqName = $seqid9;
				my $qSize = $L[10];
				my $qStart = $L[11];
				my $qEnd = $L[12];
				my $tseqName = $seqid13;
				my $tSize = $L[14];
				my $tStart = $L[15];
				my $tEnd = $L[16];
				my $blockCount = $L[17];
				my $blockSizes = $L[18];
				my $qStarts = $L[19];
				my $tStarts = $L[20];
				
#				if($L[9] eq $Q[1] && $L[13] eq $Q[1]){
				if($L[13] eq $Q[1]){		# only when chr ID match
					for(my $p = $Q[2] - $initspan; $p < $Q[3] + $initspan; $p += 1){
						if($L[15] <= $p && $p <= $L[16]){
							my $gene_region1 = "F";
							my $gene_region2 = "F";
							my @BLK = split(/,/, $blockSizes);
							my @QS = split(/,/, $qStarts);
							my @TS = split(/,/, $tStarts);
							my $numBLK = @BLK;
							
							if(! $hhitinfo->{$j}{$L[9]}{pos0}){
								for(my $tmp_q2 = $Q[2]; $tmp_q2 < $Q[3]; $tmp_q2++){
									if($hhitinfo->{$j}{$L[9]}{pos0}){
										last;
									}
									
									if($L[15] <= $tmp_q2 && $tmp_q2 <= $L[16]){
										for(my $k = 0; $k < $numBLK; $k++){
											if($L[8] eq '+'){
												my $blk_qstart = $QS[$k];
												my $blk_qend = $QS[$k] + $BLK[$k];
												my $blk_tstart = $TS[$k];
												my $blk_tend = $TS[$k] + $BLK[$k];
												
												if($TS[$k] <= $tmp_q2 && $tmp_q2 <= $blk_tend){
													my $coef = ($blk_qend - $blk_qstart) / ($blk_tend - $blk_tstart);
													my $intc = $blk_qstart - $coef * $blk_tstart;
													my $hitpos = $coef * $tmp_q2 + $intc;
													$gene_region1 = $hitpos." | $L[8] | ($blk_qend - $blk_qstart) / ($blk_tend - $blk_tstart) | $Q[2] | $tmp_q2";
													$hhitinfo->{$j}{$L[9]}{pos0} = $hitpos;
													last;
												}
											}
											elsif($L[8] eq '-'){
												my $blk_qstart = $qSize - $QS[$k] - $BLK[$k];
												my $blk_qend = $qSize - $QS[$k];
												my $blk_tstart = $TS[$k];
												my $blk_tend = $TS[$k] + $BLK[$k];
												
												if($TS[$k] <= $tmp_q2 && $tmp_q2 <= $blk_tend){
													my $coef = ($blk_qstart - $blk_qend) / ($blk_tend - $blk_tstart);
													my $intc = $blk_qstart - $coef * $blk_tstart;
													my $hitpos = $coef * $tmp_q2 + $intc;
													$gene_region1 = $hitpos." | $L[8] | ($blk_qstart - $blk_qend) / ($blk_tend - $blk_tstart) | $Q[2] | $tmp_q2";
													$hhitinfo->{$j}{$L[9]}{pos0} = $hitpos;
													last;
												}
											}
										}
									}
								}
							}
							if(! $hhitinfo->{$j}{$L[9]}{pos1}){
								for(my $tmp_q3 = $Q[3]; $tmp_q3 > $Q[2]; $tmp_q3--){
									if($hhitinfo->{$j}{$L[9]}{pos1}){
										last;
									}
									
									if($L[15] <= $tmp_q3 && $tmp_q3 <= $L[16]){
										for(my $k = 0; $k < $numBLK; $k++){
											if($L[8] eq '+'){
												my $blk_qstart = $QS[$k];
												my $blk_qend = $QS[$k] + $BLK[$k];
												my $blk_tstart = $TS[$k];
												my $blk_tend = $TS[$k] + $BLK[$k];
												
												if($TS[$k] <= $tmp_q3 && $tmp_q3 <= $blk_tend){
													my $coef = ($blk_qend - $blk_qstart) / ($blk_tend - $blk_tstart);
													my $intc = $blk_qstart - $coef * $blk_tstart;
													my $hitpos = $coef * $tmp_q3 + $intc;
													$gene_region2 = $hitpos." | $L[8] | ($blk_qend - $blk_qstart) / ($blk_tend - $blk_tstart) | $Q[3] | $tmp_q3";
													$hhitinfo->{$j}{$L[9]}{pos1} = $hitpos;
													last;
												}
											}
											elsif($L[8] eq '-'){
												my $blk_qstart = $qSize - $QS[$k] - $BLK[$k];
												my $blk_qend = $qSize - $QS[$k];
												my $blk_tstart = $TS[$k];
												my $blk_tend = $TS[$k] + $BLK[$k];
												
												if($TS[$k] <= $tmp_q3 && $tmp_q3 <= $blk_tend){
													my $coef = ($blk_qstart - $blk_qend) / ($blk_tend - $blk_tstart);
													my $intc = $blk_qstart - $coef * $blk_tstart;
													my $hitpos = $coef * $tmp_q3 + $intc;
													$gene_region2 = $hitpos." | $L[8] | ($blk_qstart - $blk_qend) / ($blk_tend - $blk_tstart) | $Q[3] | $tmp_q3";
													$hhitinfo->{$j}{$L[9]}{pos1} = $hitpos;
													last;
												}
											}
										}
									}
								}
							}
							
							$hhitinfo->{$j}{$L[9]}{eachseq_nmatch} += $L[0];
							$hhitinfo->{$j}{$L[9]}{allpos} .= $qStart.",".$qEnd.",";
							$hhitinfo->{$j}{$L[9]}{data} .= $Q[0]."\t".$refgenome."\t".$Q[1]."\t".$Q[2]."\t".$Q[3]."\t".$initspan."\tref\t".$mode."\t".$L[0]."\t".$L[8]."\t".$gene_region1."\t".$gene_region2."\t";
							$hhitinfo->{$j}{$L[9]}{data} .= $qname."\t".$qseqName."\t".$qSize."\t".$qStart."\t".$qEnd."\t";
							$hhitinfo->{$j}{$L[9]}{data} .= $dbname."\t".$tseqName."\t".$tSize."\t".$tStart."\t".$tEnd."\t";
							$hhitinfo->{$j}{$L[9]}{data} .= $blockCount."\t".$blockSizes."\t".$qStarts."\t".$tStarts."\n";
							$cnt_hit++;
							last;
						}
					}
				}
			}
		}
	}
	foreach my $norder (@Ord){
		if($norder == 0){
			next;
		}
		my $j = $norder;
		
		if($hhitinfo->{$j}){
			my $htmp = $hhitinfo->{$j};
			my @HCRS = keys(%{$htmp});
			@HCRS = sort {$a cmp $b} @HCRS;
			my $num_HCRS = @HCRS;
			
			my $AoHCRS = [];
			foreach my $htseq (@HCRS){
				unless($hhitinfo->{$j}{$htseq}{eachseq_nmatch}){
					$hhitinfo->{$j}{$htseq}{eachseq_nmatch} = 0;
				}
				
				my @tmp;
				push(@tmp, $htseq);
				push(@tmp, $hhitinfo->{$j}{$htseq}{eachseq_nmatch});
				push(@{$AoHCRS}, \@tmp);
			}
			
			@{$AoHCRS} = sort {$b->[1] <=> $a->[1]} @{$AoHCRS};
			
			my $tophseq = $AoHCRS->[0][0];
			
			unless($hhitinfo->{$j}{$tophseq}{pos0}){
				$hhitinfo->{$j}{$tophseq}{pos0} = "null";
			}
			unless($hhitinfo->{$j}{$tophseq}{pos1}){
				$hhitinfo->{$j}{$tophseq}{pos1} = "null";
			}
			
			if($hhitinfo->{$j}{$tophseq}{allpos} && $hhitinfo->{$j}{$tophseq}{allpos} ne 'null'){
				my @AllPos0 = split(/,/, $hhitinfo->{$j}{$tophseq}{allpos});
				my @AllPos;
				my $cnt_pos = 0;
				foreach my $posx (@AllPos0){
					push(@AllPos, $posx);
					$cnt_pos++;
				}
				@AllPos = sort {$a <=> $b} @AllPos;
				my $min_pos = $AllPos[0];
				my $max_pos = $AllPos[$cnt_pos - 1];
				my $n50_pos = $AllPos[int($cnt_pos / 2)];
				my $limit_pos0 = $n50_pos - ($initspan * 1.2);
				my $limit_pos1 = $n50_pos + ($initspan * 1.2);
				
				if($hhitinfo->{$j}{$tophseq}{pos0} && $hhitinfo->{$j}{$tophseq}{pos0} ne 'null'){
					$limit_pos0 = $hhitinfo->{$j}{$tophseq}{pos0} - ($initspan * 1);
					if($limit_pos0 < 0){
						$limit_pos0 = 0;
					}
				}
				if($hhitinfo->{$j}{$tophseq}{pos1} && $hhitinfo->{$j}{$tophseq}{pos1} ne 'null'){
					$limit_pos1 = $hhitinfo->{$j}{$tophseq}{pos1} + ($initspan * 1);
				}
				
				my @HL = split(/\n/, $hhitinfo->{$j}{$tophseq}{data});
				my $AoHL0 = [];
				foreach my $line (@HL){
					my @A = split(/\t/, $line);
					push(@{$AoHL0}, \@A);
				}
				@{$AoHL0} = sort {$a->[15] <=> $b->[15]} @{$AoHL0};		# sort by qstart
				
				my $blk_qspan_th = 100 * 1000;		# if distance is > 100 kb, seprate alignment blocks within the tophit chromosome
				my $prev_qpos1 = 0;
				my $hqblk = {};
				my $num_qblk = 0;
				foreach my $A (@{$AoHL0}){
					if( ($A->[15] - $prev_qpos1) <= $blk_qspan_th){
						$hqblk->{$num_qblk}{data} .= join("\t", @{$A})."\n";
						$hqblk->{$num_qblk}{nmatch} += $A->[8];
					}
					else{
						$num_qblk++;
						$hqblk->{$num_qblk}{data} .= join("\t", @{$A})."\n";
						$hqblk->{$num_qblk}{nmatch} += $A->[8];
					}
					$prev_qpos1 = $A->[16];
				}
				
				my @NQBLK = keys(%{$hqblk});
				@NQBLK = sort {$a <=> $b} @NQBLK;
				
				my $Aoqtmp = [];
				foreach my $nqblk (@NQBLK){
					my @tmp;
					push(@tmp, $hqblk->{$nqblk}{data});
					push(@tmp, $hqblk->{$nqblk}{nmatch});
					push(@{$Aoqtmp}, \@tmp);
				}
				@{$Aoqtmp} = sort {$b->[1] <=> $a->[1]} @{$Aoqtmp};
				
				@HL = split(/\n/, $Aoqtmp->[0][0]);			# select alignment block with top score within the  tophit chromosome
				
				my $AoHL = [];
				my $sum_nmatch_re = 0;
				my @AllPosS;
				foreach my $line (@HL){
					my @A = split(/\t/, $line);
	#				if($limit_pos0 <= $A[15] && $A[15] <= $limit_pos1 && $limit_pos0 <= $A[16] && $A[16] <= $limit_pos1){
						push(@{$AoHL}, \@A);
						$sum_nmatch_re += $A[8];
						push(@AllPosS, $A[15]);
						push(@AllPosS, $A[16]);
						$cnt_selected++;
	#				}
				}
				
				my $cnt_posre = @AllPosS;
				@AllPosS = sort {$a <=> $b} @AllPosS;
				my $min_pos_re = $AllPosS[0];
				my $max_pos_re = $AllPosS[$cnt_posre - 1];
				@{$AoHL} = sort {$a->[20] <=> $b->[20]} @{$AoHL};		# sort by tstart
				
				my $gpos0 = $Q[2];
				my $gpos1 = $Q[3];
				my $gpos0f = $Q[2] - $initspan;
				my $gpos1f = $Q[3] + $initspan;
				my $gpos0n1 = $Q[2] - $neighbor_tlen;
				my $gpos1n1 = $Q[2];
				my $gpos0n2 = $Q[3];
				my $gpos1n2 = $Q[3] + $neighbor_tlen;
				my $gpos0n3 = $Q[2] - $initspan;
				my $gpos1n3 = $Q[2];
				my $gpos0n4 = $Q[3];
				my $gpos1n4 = $Q[3] + $initspan;
				my $gphitbp0 = 0;
				my $gphitbp1 = 0;
				my $gphitbp2 = 0;
				my $gphitbp3 = 0;
				my $gphitbp4 = 0;
				my $gphitbp5 = 0;
				my @Gbp0;
				my @Gbp1;
				my @Gbp2;
				my @Gbp3;
				my @Gbp4;
				my @Gbp5;
				my @Qbp0;
				my @Qbp1;
				my @Qbp2;
				my @Qbp3;
				my @Qbp4;
				my @Qbp5;
				my $dotdata_flanking_fw = "qpos,refpos\n";
				my $dotdata_flanking_rv = "qpos,refpos\n";
				
				foreach my $A (@{$AoHL}){
					$hit .= join("\t", @{$A})."\t".$limit_pos0."\t".$limit_pos1."\n";
					
					# 4317,501,	49877752,49882073 (89229,84908),	77485,81802,
					
					my $strand = $A->[9];
					my @BLK = split(/,/, $A->[23]);
					my @QS = split(/,/, $A->[24]);
					my @GS = split(/,/, $A->[25]);		# TS
					my $numBLK = @BLK;
					
					for(my $k = 0; $k < $numBLK; $k++){
						my $gstart = $GS[$k];
						my $gend = $GS[$k] + $BLK[$k];
						my $tmp_gbp = $GS[$k];
						my $tmp_qbp = $QS[$k];
						
						if($A->[9] eq '-'){
							$tmp_qbp = $A->[14] - $tmp_qbp;		# 49966981 - 49877752 = 89229
						}
						
						for(my $gp = $gstart; $gp < $gend; $gp++){
							if($gpos0 <= $gp && $gp <= $gpos1){
								$gphitbp0++;
								push(@Gbp0, $tmp_gbp);
								push(@Qbp0, $tmp_qbp);
							}
							if($gpos0f <= $gp && $gp <= $gpos1f){
								$gphitbp1++;
								push(@Gbp1, $tmp_gbp);
								push(@Qbp1, $tmp_qbp);
								
								if($A->[9] eq '+'){
									$dotdata_flanking_fw .= "$tmp_qbp,$tmp_gbp\n";
								}
								elsif($A->[9] eq '-'){
									$dotdata_flanking_rv .= "$tmp_qbp,$tmp_gbp\n";
								}
							}
							if($gpos0n1 <= $gp && $gp < $gpos1n1){
								$gphitbp2++;
								push(@Gbp2, $tmp_gbp);
								push(@Qbp2, $tmp_qbp);
							}
							if($gpos0n2 < $gp && $gp <= $gpos1n2){
								$gphitbp3++;
								push(@Gbp3, $tmp_gbp);
								push(@Qbp3, $tmp_qbp);
							}
							if($gpos0n3 <= $gp && $gp < $gpos1n3){
								$gphitbp4++;
								push(@Gbp4, $tmp_gbp);
								push(@Qbp4, $tmp_qbp);
							}
							if($gpos0n4 < $gp && $gp <= $gpos1n4){
								$gphitbp5++;
								push(@Gbp5, $tmp_gbp);
								push(@Qbp5, $tmp_qbp);
							}
							
							$tmp_gbp++;
							if($A->[9] eq '+'){
								$tmp_qbp++;
							}
							elsif($A->[9] eq '-'){
								$tmp_qbp--;
							}
						}
					}
				}
				
				@Gbp0 = sort {$a <=> $b} @Gbp0;
				@Gbp1 = sort {$a <=> $b} @Gbp1;
				@Gbp2 = sort {$a <=> $b} @Gbp2;
				@Gbp3 = sort {$a <=> $b} @Gbp3;
				@Gbp4 = sort {$a <=> $b} @Gbp4;
				@Gbp5 = sort {$a <=> $b} @Gbp5;
				@Qbp0 = sort {$a <=> $b} @Qbp0;
				@Qbp1 = sort {$a <=> $b} @Qbp1;
				@Qbp2 = sort {$a <=> $b} @Qbp2;
				@Qbp3 = sort {$a <=> $b} @Qbp3;
				@Qbp4 = sort {$a <=> $b} @Qbp4;
				@Qbp5 = sort {$a <=> $b} @Qbp5;
				
				if(! $Qbp1[0] || ! $Qbp1[$gphitbp1 - 1]){
					if($outplot){
						print " partition [$t] : cannot draw plot for [$Q[0]] because of no alignment...\n";
						$failed_list .= $Q[0]."\tmissing Qbp1 (cannot draw plot)\n";
					}
					else{
						$failed_list .= $Q[0]."\tmissing Qbp1 (cannot draw plot)\n";
					}
					
					my $gpos0f = $Q[2] - $initspan;
					my $gpos1f = $Q[3] + $initspan;
					$summary2 .= $hsample->{$j}."\t".$refgenome."\t".$Q[1]."\t".$gpos0f."\t".$gpos1f."\t".$Q[2]."\t".$Q[3]."\t".$hsample->{$j}."\t-\t-\t-\t-\t";
					$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tabsent\t".$Q[0]."\t";
					$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\n";
					next;
				}
				
				#-----------------------------------------------------------------------//
				if($outplot){
					my $png1 = $Q[0]."_".$qpref."_vs_".$refgenome."_plot1.png";
					my $png2 = $Q[0]."_".$qpref."_vs_".$refgenome."_plot2.png";
					my $png3 = $Q[0]."_".$qpref."_vs_".$refgenome."_plot3.png";
					my $png4 = $Q[0]."_".$qpref."_vs_".$refgenome."_plot4.png";
					my $tmp_dotcsv_fw = $Q[0]."_".$qpref."_vs_".$refgenome."_fw.csv";
					my $tmp_dotcsv_rv = $Q[0]."_".$qpref."_vs_".$refgenome."_rv.csv";
					my $tmp_Rscript = "./$outplot_dir/".$Q[0]."_".$qpref."_vs_".$refgenome."_Rsc.txt";
					
					my $q2nb = $Q[2] - $neighbor_tlen;
					my $q3nb = $Q[3] + $neighbor_tlen;
					if($q2nb < 0){
						$q2nb = 0;
					}
					if($q3nb < 0){
						$q3nb = 0;
					}
					
					my $vline;
					$vline .= "abline(v\=$Q[2], col\=\"black\")";
					$vline .= "abline(v\=$Q[3], col\=\"black\")";
					$vline .= "abline(v\=$q2nb, col\=\"black\", lty=3)\n";
					$vline .= "abline(v\=$q3nb, col\=\"black\", lty=3)\n";
					my $hline;
					$hline .= "abline(h\=$Q[2], col\=\"black\")";
					$hline .= "abline(h\=$Q[3], col\=\"black\")";
					$hline .= "abline(h\=$q2nb, col\=\"black\", lty=3)\n";
					$hline .= "abline(h\=$q3nb, col\=\"black\", lty=3)\n";
					
					my $qpos_min = $Qbp1[0];
					my $qpos_max = $Qbp1[$gphitbp1 - 1];
					my $tpos_min = $Gbp1[0];
					my $tpos_max = $Gbp1[$gphitbp1 - 1];
					
					my $Rscript =<<"EOS";
workingDir = "$wpath/$dir/$outplot_dir"
setwd(workingDir)
getwd()

fw <- read.csv("$tmp_dotcsv_fw", header=T)
rv <- read.csv("$tmp_dotcsv_rv", header=T)

png("$png1", width=2400, height=2400)
plot(fw[,1], fw[,2], pch=20, col="blue", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3, xlab="$tophseq ($qpref)", ylab="$Q[0] ($refgenome, $Q[1], $Q[2] - $Q[3])")
points(rv[,1], rv[,2], pch=20, col="red", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3)
abline(h=$tpos_min, col="black")
abline(h=$tpos_max, col="black")
abline(v=$qpos_min, col="black")
abline(v=$qpos_max, col="black")
$hline
dev.off()

png("$png2", width=2400, height=2400)
plot(fw[,1], fw[,2], pch=20, col="red", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3, xlab="$tophseq ($qpref)", ylab="$Q[0] ($refgenome, $Q[1], $Q[2] - $Q[3])")
points(rv[,1], rv[,2], pch=20, col="blue", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3)
abline(h=$tpos_min, col="black")
abline(h=$tpos_max, col="black")
abline(v=$qpos_min, col="black")
abline(v=$qpos_max, col="black")
$hline
dev.off()

png("$png3", width=2400, height=2400)
plot(fw[,2], fw[,1], pch=20, col="blue", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3, xlab="$Q[0] ($refgenome, $Q[1], $Q[2] - $Q[3])", ylab="$tophseq ($qpref)")
points(rv[,2], rv[,1], pch=20, col="red", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3)
abline(v=$tpos_min, col="black")
abline(v=$tpos_max, col="black")
abline(h=$qpos_min, col="black")
abline(h=$qpos_max, col="black")
$vline
dev.off()

png("$png4", width=2400, height=2400)
plot(fw[,2], fw[,1], pch=20, col="red", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3, xlab="$Q[0] ($refgenome, $Q[1], $Q[2] - $Q[3])", ylab="$tophseq ($qpref)")
points(rv[,2], rv[,1], pch=20, col="blue", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3)
abline(v=$tpos_min, col="black")
abline(v=$tpos_max, col="black")
abline(h=$qpos_min, col="black")
abline(h=$qpos_max, col="black")
$vline
dev.off()

EOS
					
					if($dotdata_flanking_fw){
						open(my $rfh, ">", "./$outplot_dir/$tmp_dotcsv_fw");
						print $rfh $dotdata_flanking_fw;
						close $rfh;
					}
					if($dotdata_flanking_rv){
						open(my $rfh, ">", "./$outplot_dir/$tmp_dotcsv_rv");
						print $rfh $dotdata_flanking_rv;
						close $rfh;
					}
					if($Rscript){
						open(my $rfh, ">", $tmp_Rscript);
						print $rfh $Rscript;
						close $rfh;
					}
					
					my $Rcmd = "Rscript --vanilla --slave $tmp_Rscript >/dev/null 2>&1";
					#print "! cmd=[$Rcmd]\n";
					system("$Rcmd");
					
					if(-e "Rplots.pdf"){
						system("rm Rplots.pdf");
					}
				}
				#-----------------------------------------------------------------------//
				
				$summary .= $Q[0]."\t".$refgenome."\t".$Q[1]."\t".$Q[2]."\t".$Q[3]."\t".$initspan."\tref\t".$hsample->{$j}."\t".$gphitbp0."\t".$gphitbp1."\t".$tophseq."\t".$hhitinfo->{$j}{$tophseq}{pos0}."\t".$hhitinfo->{$j}{$tophseq}{pos1}."\t".$limit_pos0."\t".$limit_pos1."\t".$num_HCRS;
				$summary .= "\n";
				
				if(! $Qbp0[0] || ! $Qbp0[$gphitbp0 - 1]){
					$failed_list .= $Q[0]."\tmissing Qbp0 (found alignment within the flanking region but cannot find gene position)\n";
					
					my $gpos0f = $Q[2] - $initspan;
					my $gpos1f = $Q[3] + $initspan;
					$summary2 .= $hsample->{$j}."\t".$refgenome."\t".$Q[1]."\t".$gpos0f."\t".$gpos1f."\t".$Q[2]."\t".$Q[3]."\t".$hsample->{$j}."\t-\t-\t-\t-\t";
					$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tabsent\t".$Q[0]."\t";
					$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\n";
					next;
				}
				
				my $bp_sp01 = abs($Q[2] - $Q[3]) + 1;
				my $bp_qsp01 = abs($Qbp0[0] - $Qbp0[$gphitbp0 - 1]) + 1;
				my $bp_ins = $bp_qsp01 - $bp_sp01;
				if($bp_ins < 1){
					$bp_ins = "-";
				}
				my $ratio_galign = $gphitbp0 / $bp_sp01;
				my $ratio_seqnormality = $ratio_galign;
				
				my $ratio_ins = "-";
				if($bp_qsp01 > 0){
					$ratio_ins = $bp_sp01 / $bp_qsp01;
					
					if($ratio_ins >= 1){
						$ratio_ins = "-";
					}
					elsif($ratio_seqnormality > $ratio_ins){
						$ratio_seqnormality = $ratio_ins;
					}
				}
				
				my $bp_5neighbor = 0;
				my $bp_3neighbor = 0;
				my $bp_5flanking = 0;
				my $bp_3flanking = 0;
				my $ratio_5neighbor = 0;
				my $ratio_3neighbor = 0;
				my $ratio_5flanking = 0;
				my $ratio_3flanking = 0;
				
				if($gphitbp2 > 1 && abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]) > 0){
					$bp_5neighbor = $gphitbp2;
					$ratio_5neighbor = $gphitbp2 / abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]);
					
					if($Q[4] eq '-'){
						$bp_3neighbor = $gphitbp2;
						$ratio_3neighbor = $gphitbp2 / abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]);
					}
				}
				if($gphitbp3 > 1 && abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]) > 0){
					$bp_3neighbor = $gphitbp3;
					$ratio_3neighbor = $gphitbp3 / abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]);
					
					if($Q[4] eq '-'){
						$bp_5neighbor = $gphitbp3;
						$ratio_5neighbor = $gphitbp3 / abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]);
					}
				}
				if($gphitbp4 > 1 && abs($Gbp4[0] - $Gbp4[$gphitbp4 - 1]) > 0){
					$bp_5flanking = $gphitbp4;
					$ratio_5flanking = $gphitbp4 / abs($Gbp4[0] - $Gbp4[$gphitbp4 - 1]);
					
					if($Q[4] eq '-'){
						$bp_3flanking = $gphitbp4;
						$ratio_3flanking = $gphitbp4 / abs($Gbp4[0] - $Gbp4[$gphitbp4 - 1]);
					}
				}
				if($gphitbp5 > 1 && abs($Gbp5[0] - $Gbp5[$gphitbp5 - 1]) > 0){
					$bp_3flanking = $gphitbp5;
					$ratio_3flanking = $gphitbp5 / abs($Gbp5[0] - $Gbp5[$gphitbp5 - 1]);
					
					if($Q[4] eq '-'){
						$bp_5flanking = $gphitbp5;
						$ratio_5flanking = $gphitbp5 / abs($Gbp5[0] - $Gbp5[$gphitbp5 - 1]);
					}
				}
				
				my $protein_coding = "FALSE";
				if($Q[6] && $Q[6] > 0){
					$protein_coding = "TRUE";
				}
				
				my $qseq_woN = "-";
				my $bp_qseq_woN = "-";
				if(! defined $hqseq->{$tophseq}{len}){
					print "! error : missing seq info for [$tophseq]\n";
					die;
				}
				if($hqseq->{$tophseq}{len} > $Qbp0[$gphitbp0 - 1]){
					$qseq_woN = substr($hqseq->{$tophseq}{seq}, $Qbp0[0] - 1, $bp_qsp01);
					$qseq_woN =~ s/N//gi;
					$bp_qseq_woN = length($qseq_woN);
				}
				else{
					$qseq_woN = substr($hqseq->{$tophseq}{seq}, $Qbp0[0] - 1);
					$qseq_woN =~ s/N//gi;
					$bp_qseq_woN = length($qseq_woN);
				}
				
				my $judge = "present";
				if($bp_ins ne '-'){
					if($bp_ins > 5000){
						$judge = "insertion > 5kb";
					}
					elsif($bp_ins > 4000){
						$judge = "insertion > 4kb";
					}
					elsif($bp_ins > 3000){
						$judge = "insertion > 3kb";
					}
					elsif($bp_ins > 2000){
						$judge = "insertion > 2kb";
					}
					elsif($bp_ins > 1000){
						$judge = "insertion > 1kb";
					}
					elsif($bp_ins > 500){
						$judge = "insertion > 500bp";
					}
				}
				elsif($ratio_seqnormality < 0.1){
					$judge = "present (collapsed)";
				}
				elsif($ratio_seqnormality < 0.5){
					$judge = "present (collapsed)";
				}
				elsif($ratio_seqnormality < 0.8){
					$judge = "present (partly)";
				}
				
				$summary2 .= $hsample->{$j}."\t".$refgenome."\t".$Q[1]."\t".$gpos0f."\t".$gpos1f."\t".$Q[2]."\t".$Q[3]."\t".$hsample->{$j}."\t".$tophseq."\t".$Qbp1[0]."\t".$Qbp1[$gphitbp1 - 1]."\t1\t";
				$summary2 .= $Qbp0[0]."\t".$Qbp0[$gphitbp0 - 1]."\t".$bp_sp01."\t".$gphitbp0."\t".$bp_qsp01."\t".$bp_qseq_woN."\t".$bp_ins."\t".$ratio_galign."\t".$ratio_ins."\t".$ratio_seqnormality."\t".$judge."\t".$Q[0]."\t";
				$summary2 .= $bp_5neighbor."\t".$ratio_5neighbor."\t".$bp_3neighbor."\t".$ratio_3neighbor."\t".$bp_5flanking."\t".$ratio_5flanking."\t".$bp_3flanking."\t".$ratio_3flanking."\n";
			}
			else{
				$summary .= $Q[0]."\t".$refgenome."\t".$Q[1]."\t".$Q[2]."\t".$Q[3]."\t".$initspan."\tref\t".$hsample->{$j}."\t-\t-\t-\t-\t-\t-\t-\t-\n";
				
				my $gpos0f = $Q[2] - $initspan;
				my $gpos1f = $Q[3] + $initspan;
				$summary2 .= $hsample->{$j}."\t".$refgenome."\t".$Q[1]."\t".$gpos0f."\t".$gpos1f."\t".$Q[2]."\t".$Q[3]."\t".$hsample->{$j}."\t-\t-\t-\t-\t";
				$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tabsent\t".$Q[0]."\t";
				$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\n";
			}
		}
		else{
			$summary .= $Q[0]."\t".$refgenome."\t".$Q[1]."\t".$Q[2]."\t".$Q[3]."\t".$initspan."\tref\t".$hsample->{$j}."\t-\t-\t-\t-\t-\t-\t-\t-\n";
			
			my $gpos0f = $Q[2] - $initspan;
			my $gpos1f = $Q[3] + $initspan;
			$summary2 .= $hsample->{$j}."\t".$refgenome."\t".$Q[1]."\t".$gpos0f."\t".$gpos1f."\t".$Q[2]."\t".$Q[3]."\t".$hsample->{$j}."\t-\t-\t-\t-\t";
			$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tabsent\t".$Q[0]."\t";
			$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\n";
		}
	}
	
	$cntq++;
	if($cntq == $seg){
		print " partition [$t] : [$cntq] [$cnt_selected]\n";
		$seg += $add;
	}
}

print "! partition [$t] : [$cnt_hit] candidates, [$cnt_selected] hits\n";

my $rh = {};
$rh->{hit} = $hit;
$rh->{summary} = $summary;
$rh->{summary2} = $summary2;
$rh->{failed_list} = $failed_list;

return $rh;
}


#-------------------------------------------------------------------------------
sub Read_list{
my $file = shift;

print "! reading ID list [$file]...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	unless($line){
#		print "! missing line\n";
		next;
	}
	elsif($line =~ /\#/){
		next;
	}
	elsif($line =~ /gid/ && $line =~ /chr/i){
		next;
	}
	
	my @A;
	if($line =~ /,/){
		@A = split(/,/, $line);
	}
	elsif($line =~ /\t/){
		@A = split(/\t/, $line);
	}
	else{
		$A[0] = $line;
	}
	
	if($A[0]){
		$hash->{$A[0]} = 1;
		$cnt++;
	}
}
close $fh;

print " [$cnt] entries\n";

return $hash;
}


#---------------------------------------------------------------------
sub Read_seqinfo{
my $file = shift;

print "! reading [$file]...\n";
open(my $fh, "<", $file);
my $hash = {};
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	$line =~ s/\"//g;
	
	my @A = split(/\t/, $line);
	my @B = split(/,/, $A[2]);
	$hash->{$A[0]}{$A[2]} = $A[1];
	$hash->{$A[0]}{$B[0]} = $A[1];
#	print "$A[0] $B[0] $A[1]\n";
	$cnt++;
}
close $fh;

print "! [$cnt] lines\n";

return $hash;
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

print "! total [$numID] sequence ID, [$total_len] bp\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Read_gff{
my $gff3 = shift;
my $hID = shift;

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
		
		if($hID->{$gid}){
			push(@GID, $gid);
			$hash->{$gid}{data} = $gid.",".$A[0].",".$A[3].",".$A[4].",".$A[6];
		}
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
		
		if($hID->{$gid}){
			$g2t->{$gid} .= $tid."\n";
			$t2g->{$tid} = $gid;
			$hash->{$gid}{num_variant} += 1;
		}
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
		if($hID->{$gid}){
			$hash->{$gid}{num_CDS} += 1;
		}
	}
}

print " [$gcnt] genes\n";

#my $r = "gid,chr,pos0,pos1,strand,num_variant,total_num_CDS\n";
my $r = "";
my $numpc = 0;
my $numnc = 0;
my $numpc_tr = 0;
my $hcnt = {};
foreach my $gid (@GID){
	unless($hash->{$gid}{num_variant}){
		$hash->{$gid}{num_variant} = 0;
	}
	unless($hash->{$gid}{num_CDS}){
		$numnc++;
		$hash->{$gid}{num_CDS} = 0;
		$hash->{$gid}{protein_coding} = "false";
	}
	else{
		$numpc++;
		$numpc_tr += $hash->{$gid}{num_variant};
		$hash->{$gid}{protein_coding} = "true";
	}
	
	my $tmp = $hash->{$gid}{data}.",".$hash->{$gid}{num_variant}.",".$hash->{$gid}{num_CDS};
	my @B = split(/,/, $tmp);
	$hcnt->{$B[1]} += 1;
	$r .= $tmp."\n";
}

print " [$numpc] protein-coding, [$numpc_tr] transcripts\n";
print " [$numnc] non-coding\n";

#my @Chrs = keys(%{$hcnt});
#@Chrs = sort {$a cmp $b} @Chrs;
#
#$r .= "\nChr,num gene\n";
#foreach my $sid (@Chrs){
#	$r .= $sid.",".$hcnt->{$sid}."\n";
#}

my @RL = split(/\n/, $r);
my $rh = {};
$rh->{RL} = \@RL;
$rh->{hash} = $hash;
$rh->{g2t} = $g2t;
$rh->{t2g} = $t2g;

return $rh;
}



