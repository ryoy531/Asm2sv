#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use Getopt::Long;
use FindBin;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "unitesummary.pl version 1.42\n";
$version .= "last update: [2023\/2\/7]\n";
$version .= "copyright: ryoichi yano\n";

#print $version;

#-------------------------------------------------------------------------------

my $data_path = getcwd;
unless($data_path){
	print "! cannot find data PATH...\n";
	goto END;
}

my $file_str = shift;
my $refgenome = shift;
my $chrinfo = shift;
my $qorf_str = shift;
my $zip = shift;

unless($file_str){
	print "! missing sv file string...\n";
	goto END;
}
unless($qorf_str){
	print "! missing qorf file string...\n";
	goto END;
}
unless($refgenome){
	$refgenome = "null";
}
unless($chrinfo){
	$chrinfo = "null";
}

my @F0 = split(/,/, $file_str);
my @Q0 = split(/,/, $qorf_str);
my $numF0 = @F0;
my $numQ0 = @Q0;

my $AoF = [];
my $err = 0;
for(my $i = 0; $i < $numF0; $i++){
	if(! $Q0[$i] || ! -e $Q0[$i]){
		$Q0[$i] = "null";
	}
	if($F0[$i] && -e $F0[$i]){
		my @tmp;
		push(@tmp, $F0[$i]);
		push(@tmp, $Q0[$i]);
		push(@{$AoF}, \@tmp);
	}
}

print "\n";
my $hseqinfo = {};
if($chrinfo ne 'null' && -e $chrinfo){
	$hseqinfo = Read_seqinfo($chrinfo);
}

print "! collecting data...\n";
my $hash = {};
foreach my $tmp (@{$AoF}){
	my $file = $tmp->[0];
	
	$hash = Collectdata($file, $hash, $refgenome, $hseqinfo);
	if($hash->{err}){
		print "! abort script due to error in [$chrinfo]...\n";
		goto END;
	}
}

print "! summarizing data...\n";
my $hginfo = $hash->{ginfo};
my @GID = keys(%{$hginfo});
my $AoG1 = [];
my $AoG2 = [];
foreach my $gid (@GID){
	my @tmp;
	push(@tmp, $gid);
	push(@tmp, $hginfo->{$gid}{Chr});
	push(@tmp, $hginfo->{$gid}{sp0});		# gene position in GFF
	
	if($hginfo->{$gid}{Chr} !~ /scaffold/){
		push(@{$AoG1}, \@tmp);
	}
	else{
		push(@{$AoG2}, \@tmp);
	}
}

my $hsortedGID = {};
my $norder = 0;
my $tmp_seqid = "-";
my @GID2;
if(@{$AoG1}){
	@{$AoG1} = sort {$a->[2] <=> $b->[2]} @{$AoG1};		# sort by gene position
	@{$AoG1} = sort {$a->[1] cmp $b->[1]} @{$AoG1};		# sort by chr ID
	
	foreach my $tmp (@{$AoG1}){
		if($tmp_seqid ne $tmp->[1]){
			$tmp_seqid = $tmp->[1];
			$norder++;
		}
		$hsortedGID->{$norder}{gid} .= $tmp->[0]."\n";
		$hsortedGID->{$norder}{seqid} = $tmp->[1];
		$hsortedGID->{$norder}{num} += 1;
		push(@GID2, $tmp->[0]);
	}
}
if(@{$AoG2}){
	@{$AoG2} = sort {$a->[2] <=> $b->[2]} @{$AoG2};		# sort by gene position
	@{$AoG2} = sort {$a->[1] cmp $b->[1]} @{$AoG2};		# sort by chr ID
	
	foreach my $tmp (@{$AoG2}){
		if($tmp_seqid ne $tmp->[1]){
			$tmp_seqid = $tmp->[1];
			$norder++;
		}
		$hsortedGID->{$norder}{gid} .= $tmp->[0]."\n";
		$hsortedGID->{$norder}{seqid} = $tmp->[1];
		$hsortedGID->{$norder}{num} += 1;
		push(@GID2, $tmp->[0]);
	}
}

print "! [$norder] sequence ID\n";

my $hpav = $hash->{sv};
my @Samples = keys(%{$hpav});
@Samples = sort {$a cmp $b} @Samples;

my @SHeader;
my @Sorted_samples;
foreach my $id (@Samples){
	if($refgenome ne 'null'){
		if($id eq $refgenome){
			push(@SHeader, "$id [reference]");
			push(@Sorted_samples, $id);
		}
	}
}
foreach my $id (@Samples){
	if($refgenome ne 'null'){
		if($id ne $refgenome){
			push(@SHeader, $id);
			push(@Sorted_samples, $id);
		}
	}
	else{
		push(@SHeader, $id);
		push(@Sorted_samples, $id);
	}
}
@Samples = @Sorted_samples;

#my $r1 = "gene_id,Chr,sp0,sp1,cnt 0-0.3,cnt 0.3-0.6,cnt 0.6-0.9,cnt >=0.9,reference,".join(",", @Samples)."\n";
my $r1 = "gene_id,Chr,pos0,pos1,cnt 0-0.3,cnt 0.3-0.6,cnt 0.6-0.9,cnt >=0.9,".join(",", @SHeader)."\n";
my $r1B = "gene_id,Chr,pos0,pos1,cnt 0-0.6,cnt 0.6-1.2,cnt 1.2-0.1.8,cnt >=1.8,".join(",", @SHeader)."\n";
my $r1_nmatrix = "gene_id,".join(",", @SHeader)."\n";
my $r1B_nmatrix = "gene_id,".join(",", @SHeader)."\n";
#my $r1 = "gene_id,reference,".join(",", @Samples)."\n";
$r1 =~ s/\.genome//g;
$r1 =~ s/_GCA_022114995_renamed//g;
my $r2 = $r1;
my $r3 = $r1;
my $r4 = $r1;
my $r5 = $r1;
my $r6 = $r1;
my $r7 = $r1;
my $r8 = $r1;
my $r9 = $r1;
my $r10 = $r1;
my $r11 = $r1;
my $r11B = $r1B;
my $r12 = $r1;
my $r15 = $r1;
my $r17 = $r1;
my $r19 = $r1;
my $r101 = $r1;
my $r101M = $r1;
my $r101B = $r1;
my $r102 = $r1;
my $r103 = $r1;
my $r104 = $r1;
my $r105 = $r1;
my $r201 = $r1;
my $r201M = $r1;
my $r201B = $r1;
my $r201i = $r1;
my $r202 = $r1;
my $r203 = $r1;
my $r204 = $r1;
my $r205 = $r1;
my $r301 = $r1;
my $r301M = $r1;
my $r301B = $r1;
my $r302 = $r1;
my $r303 = $r1;
my $r304 = $r1;
my $r305 = $r1;
my $r11_nmatrix = $r1_nmatrix;
my $r15_nmatrix = $r1_nmatrix;
my $r17_nmatrix = $r1_nmatrix;
my $r19_nmatrix = $r1_nmatrix;
my $r11B_nmatrix = $r1B_nmatrix;
my $r11R_nmatrix = $r1B_nmatrix;
my $r11N_nmatrix = $r1B_nmatrix;
my $r11Z_nmatrix = $r1B_nmatrix;
my $r11L_nmatrix = $r1B_nmatrix;
my $r11RZ_nmatrix = $r1B_nmatrix;
my $r11RL_nmatrix = $r1B_nmatrix;

#my $r_hist = "gene_id,Chr,sp0,sp1,gene 0-0.25,gene 0.25-0.50,gene 0.50-0.75,gene 0.75-1.0,gene 1.0-,.,protein 0-0.25,protein 0.25-0.50,protein 0.50-0.75,protein 0.75-1.0,protein 1.0-,.,promoter 0-0.25,promoter 0.25-0.50,promoter 0.50-0.75,promoter 0.75-1.0,promoter 1.0-,.,3'-UTR 0-0.25,3'-UTR 0.25-0.50,3'-UTR 0.50-0.75,3'-UTR 0.75-1.0,3'-UTR 1.0-\n";
my $r_hist = "gene_id,Chr,pos0,pos1,refgenome_judge,";
$r_hist .= "gene 0-0.3,gene 0.3-0.6,gene 0.6-0.9,gene >=0.9,.,";
$r_hist .= "protein 0-0.3,protein 0.3-0.6,protein 0.6-0.9,protein >=0.9,.,";
$r_hist .= "promoter 0-0.3,promoter 0.3-0.6,promoter 0.6-0.9,promoter >=0.9,.,";
$r_hist .= "3'-UTR 0-0.3,3'-UTR 0.3-0.6,3'-UTR 0.6-0.9,3'-UTR >=0.9\n";

my $cnt0 = 0;
my $cnt1 = 0;
my $cnt_geno1 = 0;
my $cnt_geno1_rm = 0;
my $cnt_geno2 = 0;
my $cnt_geno3 = 0;
my $cnt_geno4 = 0;
my $cnt_geno5 = 0;
my $cnt_geno6 = 0;
my $cnt_geno7 = 0;
my $cnt_geno8 = 0;
my $cnt_geno9 = 0;
my $cnt_geno10 = 0;
my $cnt_geno11 = 0;
my $cnt_geno11_rm = 0;
my $cnt_geno12 = 0;
my $cnt_geno13 = 0;
my $cnt_geno14 = 0;
my $cnt_geno15 = 0;
my $cnt_geno16 = 0;
my $cnt_geno17 = 0;
my $cnt_geno18 = 0;
my $cnt_geno19 = 0;
my $cnt_geno101 = 0;
my $cnt_geno102 = 0;
my $cnt_geno103 = 0;
my $cnt_geno104 = 0;
my $cnt_geno105 = 0;
my $cnt_geno201 = 0;
my $cnt_geno202 = 0;
my $cnt_geno203 = 0;
my $cnt_geno204 = 0;
my $cnt_geno205 = 0;
my $cnt_geno301 = 0;
my $cnt_geno302 = 0;
my $cnt_geno303 = 0;
my $cnt_geno304 = 0;
my $cnt_geno305 = 0;
my $cnt_reffail = 0;
my $cnt_diffchr = 0;
my $cnt_trans = 0;
my $hgeneclass = {};
my $hgeneinfo = {};
my $hhpos = {};
my $nrange = 20;
my $interblk_th = 1000 * 1000;
my $blknumgid_th = int($nrange / 4);

for(my $nc = 1; $nc <= $norder; $nc++){
	print " - $hsortedGID->{$nc}{seqid}\n";
	my @eachGID = split(/\n/, $hsortedGID->{$nc}{gid});
	
	for(my $ig = 0; $ig < $hsortedGID->{$nc}{num}; $ig++){
		my $gid = $eachGID[$ig];
		my $hjudge = {};
		my $num_conserved = 0;
		my $num_conserved_prom = 0;
		my $num_conserved_utr3 = 0;
		my $num_nonzero1 = 0;
		my $num_nonzero2 = 0;
		my $num_nonzero3 = 0;
		my $num_nonzero4 = 0;
		my $num_nonzero5 = 0;
		my $num_nodata1 = 0;
		my $num_nodata2 = 0;
		my $num_nodata3 = 0;
		my $num_nodata4 = 0;
		my $num_nodata5 = 0;
		my $num_diffchr = 0;
		my $num_trans = 0;
		my $num_genotyped = 0;
		my $num_misgenotype = 0;
		my @Geno1;
		my @Geno1B;
		my @Geno2;
		my @Geno3;
		my @Geno4;
		my @Geno5;
		my @Geno101;
		my @Geno101B;
		my @Geno102;
		my @Geno103;
		my @Geno104;
		my @Geno105;
		my @Geno201;
		my @Geno201B;
		my @Geno201i;
		my @Geno202;
		my @Geno203;
		my @Geno204;
		my @Geno205;
		my @Geno301;
		my @Geno301B;
		my @Geno302;
		my @Geno303;
		my @Geno304;
		my @Geno305;
	#	push(@Geno1, 1);	#reference
	#	push(@Geno2, 1);	#reference
	#	push(@Geno3, 1);	#reference
	#	push(@Geno4, 1);	#reference
	#	push(@Geno5, 1);	#reference
		
		#-----------------------------------------------------//
		my $failprot_refgenome = 0;
		foreach my $id (@Samples){
			if($refgenome ne 'null'){
				if($id eq $refgenome){
					if(! $hash->{svprot}{$id}{$gid} || $hash->{svprot}{$id}{$gid} eq '-'){
						$failprot_refgenome = 1;
					}
					elsif($hash->{svprot}{$id}{$gid} < 0.97){
						$failprot_refgenome = 1;
					}
					last;
				}
			}
		}
		
		my $judge_refgenome = "-";
		foreach my $id (@Samples){
			if(! $hash->{sv}{$id}{$gid} || $hash->{sv}{$id}{$gid} eq '-'){
				$hash->{sv}{$id}{$gid} = 0;
				$hash->{raw}{$id}{$gid} = 0;
				$num_nodata1++;
			}
			elsif($hash->{sv}{$id}{$gid} > 1){
				$hash->{sv}{$id}{$gid} = 1;
			}
			
			if(! $hash->{sv2}{$id}{$gid} || $hash->{sv2}{$id}{$gid} eq '-'){
				$hash->{sv2}{$id}{$gid} = 0;
			}
			
			if(! $hash->{svprot}{$id}{$gid} || $hash->{svprot}{$id}{$gid} eq '-'){
				$hash->{svprot}{$id}{$gid} = 0;
				$num_nodata3++;
			}
			elsif($hash->{svprot}{$id}{$gid} > 1){
				$hash->{svprot}{$id}{$gid} = 1;
			}
			
			if(! $hash->{svprom}{$id}{$gid} || $hash->{svprom}{$id}{$gid} eq '-'){
				$hash->{svprom}{$id}{$gid} = 0;
				$num_nodata4++;
			}
#			elsif($hash->{svprom}{$id}{$gid} > 1){
#				$hash->{svprom}{$id}{$gid} = 1;
#			}
			
			if(! $hash->{svutr3}{$id}{$gid} || $hash->{svutr3}{$id}{$gid} eq '-'){
				$hash->{svutr3}{$id}{$gid} = 0;
				$num_nodata5++;
			}
#			elsif($hash->{svutr3}{$id}{$gid} > 1){
#				$hash->{svutr3}{$id}{$gid} = 1;
#			}
			
			if($hash->{svjudge2}{$id}{$gid} && $hash->{svjudge2}{$id}{$gid} eq 'genotyped'){
				$num_genotyped++;
			}
			else{
				$num_misgenotype++;
			}
			
			my $fgeno = $hash->{geno_including_flanking}{$id}{$gid};		# true if flanking seq len >=5000 and sv ratio > 0
			$hjudge->{$fgeno} += 1;
			
			if($hash->{sv}{$id}{$gid} >= 0.95){
				if($refgenome ne 'null'){
					if($id ne $refgenome){
						$num_conserved++;
					}
				}
				else{
					$num_conserved++;
				}
			}
			if(0.95 <= $hash->{svprom}{$id}{$gid} && $hash->{svprom}{$id}{$gid} >= 1.05){
				if($refgenome ne 'null'){
					if($id ne $refgenome){
						$num_conserved_prom++;
					}
				}
				else{
					$num_conserved_prom++;
				}
			}
			if(0.95 <= $hash->{svutr3}{$id}{$gid} && $hash->{svutr3}{$id}{$gid} <= 1.05){
				if($refgenome ne 'null'){
					if($id ne $refgenome){
						$num_conserved_utr3++;
					}
				}
				else{
					$num_conserved_utr3++;
				}
			}
			if($hash->{sv}{$id}{$gid} > 0){
				$num_nonzero1++;
			}
			if($hash->{raw}{$id}{$gid} > 0){
				$num_nonzero2++;
			}
			if($hash->{svprot}{$id}{$gid} > 0){
				$num_nonzero3++;
			}
			if($hash->{svprom}{$id}{$gid} > 0){
				$num_nonzero4++;
			}
			if($hash->{svutr3}{$id}{$gid} > 0){
				$num_nonzero5++;
			}
			
			$hash->{sv}{$id}{$gid} = sprintf("%.3f", $hash->{sv}{$id}{$gid});
			$hash->{sv2}{$id}{$gid} = sprintf("%.3f", $hash->{sv2}{$id}{$gid});
			$hash->{raw}{$id}{$gid} = sprintf("%.3f", $hash->{raw}{$id}{$gid});
			$hash->{svprot}{$id}{$gid} = sprintf("%.3f", $hash->{svprot}{$id}{$gid});
			$hash->{svprom}{$id}{$gid} = sprintf("%.3f", $hash->{svprom}{$id}{$gid});
			$hash->{svutr3}{$id}{$gid} = sprintf("%.3f", $hash->{svutr3}{$id}{$gid});
			
			push(@Geno1, $hash->{sv}{$id}{$gid});
			push(@Geno1B, $hash->{sv2}{$id}{$gid});
			push(@Geno2, $hash->{raw}{$id}{$gid});
			push(@Geno3, $hash->{svprot}{$id}{$gid});
			push(@Geno4, $hash->{svprom}{$id}{$gid});
			push(@Geno5, $hash->{svutr3}{$id}{$gid});
			
			if($refgenome ne 'null'){
				if($id eq $refgenome){
					$judge_refgenome = $fgeno;
				}
				
				my $refchr = $hash->{refchr}{$refgenome}{$gid};
				my $refpos0 = $hash->{refpos0}{$refgenome}{$gid};
				my $refpos1 = $hash->{refpos1}{$refgenome}{$gid};
				
				my $judge_diffchr = "null";
				my $judge_trans = "null";
				if($hash->{hconvid}{$id}{$gid} && $hash->{hconvid}{$id}{$gid} ne '-' && $hash->{hconvid}{$id}{$gid} !~ /scaffold/){		# absent if chrinfo tsv is absent
					my $hitconvid = $hash->{hconvid}{$id}{$gid};
					my $hitseq = $hash->{hseqid}{$id}{$gid};
					my $hitpos0 = $hash->{hpos0}{$id}{$gid};
					my $hitpos1 = $hash->{hpos1}{$id}{$gid};
					
					if($hitconvid ne $refchr){
						$judge_diffchr = "true";
						$hash->{transloc_diffchr_info}{$id}{$gid} = $hitconvid." [".$hitseq.";".$hitpos0."-".$hitpos1."]";
					}
					else{
						$judge_diffchr = "false";
						
						my $ignb0 = $ig - $nrange;
						my $ignb1 = $ig + $nrange;
						if($ignb0 < 0){
							$ignb0 = 0;
						}
						if($ignb1 > $hsortedGID->{$nc}{num}){
							$ignb1 = $hsortedGID->{$nc}{num};
						}
						my @NBhits;
						for(my $ignb = $ignb0; $ignb < $ignb1; $ignb++){
							unless($eachGID[$ignb]){
								next;
							}
							
							my $nb_gid = $eachGID[$ignb];		#neighboring gene
							if($hash->{hconvid}{$id}{$nb_gid}){
								my $nb_hitconvid = $hash->{hconvid}{$id}{$nb_gid};
								my $nb_hitseq = $hash->{hseqid}{$id}{$nb_gid};
								my $nb_hitpos0 = $hash->{hpos0}{$id}{$nb_gid};
								my $nb_hitpos1 = $hash->{hpos1}{$id}{$nb_gid};
								
								if($nb_hitconvid ne '-' && $nb_hitconvid eq $refchr && $nb_hitpos0 ne '-' && $nb_hitpos1 ne '-'){
									push(@NBhits, $nb_hitpos0);
									push(@NBhits, $nb_hitpos1);
								}
							}
						}
						if(@NBhits){
							@NBhits = sort {$a <=> $b} @NBhits;
							
							my $prev_nb_hitpos = -1000000;
							my $hblk = {};
							my $nblk = 0;
							foreach my $nb_hitpos (@NBhits){
								if( abs($nb_hitpos - $prev_nb_hitpos) > $interblk_th){
									$nblk++;
									$hblk->{$nblk}{min_pos} = $nb_hitpos;
									$hblk->{$nblk}{max_pos} = $nb_hitpos;
									$hblk->{$nblk}{num} += 1;
								}
								else{
									$hblk->{$nblk}{max_pos} = $nb_hitpos;
									$hblk->{$nblk}{num} += 1;
								}
								$prev_nb_hitpos = $nb_hitpos;
							}
							
							my @NBLKs = keys(%{$hblk});
							@NBLKs = sort {$a <=> $b} @NBLKs;
							foreach my $nblk (@NBLKs){
								if($hblk->{$nblk}{min_pos} <= $hitpos0 && $hitpos1 <= $hblk->{$nblk}{max_pos}){
									if($hblk->{$nblk}{num} >= $blknumgid_th){		# within gene block
										$judge_trans = "false";
									}
									else{
										$judge_trans = "true";
									}
								}
							}
						}
					}
				}
				
				if($judge_diffchr eq 'true' || $judge_trans eq 'true'){
					if($judge_diffchr eq 'true'){
						$num_diffchr++;
						
						push(@Geno201, $hash->{sv}{$id}{$gid});
						push(@Geno201B, $hash->{sv2}{$id}{$gid});
						push(@Geno202, $hash->{raw}{$id}{$gid});
						push(@Geno203, $hash->{svprot}{$id}{$gid});
						push(@Geno204, $hash->{svprom}{$id}{$gid});
						push(@Geno205, $hash->{svutr3}{$id}{$gid});
						push(@Geno201i, $hash->{transloc_diffchr_info}{$id}{$gid});
					}
					else{
						push(@Geno201, 0);		# pseudo value
						push(@Geno201B, 0);		# pseudo value
						push(@Geno202, 0);		# pseudo value
						push(@Geno203, 0);		# pseudo value
						push(@Geno204, 0);		# pseudo value
						push(@Geno205, 0);		# pseudo value
						push(@Geno201i, 0);		# pseudo value
					}
					if($judge_trans eq 'true'){
						$num_trans++;
						
						push(@Geno301, $hash->{sv}{$id}{$gid});
						push(@Geno301B, $hash->{sv2}{$id}{$gid});
						push(@Geno302, $hash->{raw}{$id}{$gid});
						push(@Geno303, $hash->{svprot}{$id}{$gid});
						push(@Geno304, $hash->{svprom}{$id}{$gid});
						push(@Geno305, $hash->{svutr3}{$id}{$gid});
					}
					else{
						push(@Geno301, 0);		# pseudo value
						push(@Geno301B, 0);		# pseudo value
						push(@Geno302, 0);		# pseudo value
						push(@Geno303, 0);		# pseudo value
						push(@Geno304, 0);		# pseudo value
						push(@Geno305, 0);		# pseudo value
					}
					
					push(@Geno101, $hash->{sv}{$id}{$gid});
					push(@Geno101B, $hash->{sv2}{$id}{$gid});
					push(@Geno102, $hash->{raw}{$id}{$gid});
					push(@Geno103, $hash->{svprot}{$id}{$gid});
					push(@Geno104, $hash->{svprom}{$id}{$gid});
					push(@Geno105, $hash->{svutr3}{$id}{$gid});
				}
				else{
					push(@Geno101, 0);		# pseudo value
					push(@Geno101B, 0);		# pseudo value
					push(@Geno102, 0);		# pseudo value
					push(@Geno103, 0);		# pseudo value
					push(@Geno104, 0);		# pseudo value
					push(@Geno105, 0);		# pseudo value
					
					push(@Geno201, 0);		# pseudo value
					push(@Geno201B, 0);		# pseudo value
					push(@Geno202, 0);		# pseudo value
					push(@Geno203, 0);		# pseudo value
					push(@Geno204, 0);		# pseudo value
					push(@Geno205, 0);		# pseudo value
					push(@Geno201i, 0);		# pseudo value
					
					push(@Geno301, 0);		# pseudo value
					push(@Geno301B, 0);		# pseudo value
					push(@Geno302, 0);		# pseudo value
					push(@Geno303, 0);		# pseudo value
					push(@Geno304, 0);		# pseudo value
					push(@Geno305, 0);		# pseudo value
				}
			}
		}
		#-----------------------------------------------------//
		
		if($hjudge->{true}){
			$hgeneinfo->{num_true}{$gid} = $hjudge->{true};
		}
		else{
			$hgeneinfo->{num_true}{$gid} = 0;
		}
		($hgeneinfo->{hist_value1}{$gid}, $hgeneinfo->{hist_header1}) = Array2histcsv(\@Geno1);
		($hgeneinfo->{hist_value3}{$gid}, $hgeneinfo->{hist_header3}) = Array2histcsv(\@Geno3);
		($hgeneinfo->{hist_value4}{$gid}, $hgeneinfo->{hist_header4}) = Array2histcsv(\@Geno4);
		($hgeneinfo->{hist_value5}{$gid}, $hgeneinfo->{hist_header5}) = Array2histcsv(\@Geno5);
		($hgeneinfo->{hist_value1B}{$gid}, $hgeneinfo->{hist_header1B}) = Array2histcsv2(\@Geno1B);
		
		$r_hist .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$judge_refgenome.",".$hgeneinfo->{hist_value1}{$gid}.",.,".$hgeneinfo->{hist_value3}{$gid}.",.,".$hgeneinfo->{hist_value4}{$gid}.",.,".$hgeneinfo->{hist_value5}{$gid}."\n";
		
		if($refgenome ne 'null' && $judge_refgenome ne 'true'){		# skip record if refgenome is not genotyped (assume 1 in all genes but some may fail)
			$cnt_reffail++;
			next;
		}
		
		if($refgenome ne 'null'){
			if($num_diffchr > 0 || $num_trans > 0){
				($hgeneinfo->{hist_value101}{$gid}, $hgeneinfo->{hist_header101}) = Array2histcsv(\@Geno101);
				($hgeneinfo->{hist_value101B}{$gid}, $hgeneinfo->{hist_header101B}) = Array2histcsv(\@Geno101B);
				($hgeneinfo->{hist_value103}{$gid}, $hgeneinfo->{hist_header103}) = Array2histcsv(\@Geno103);
				($hgeneinfo->{hist_value104}{$gid}, $hgeneinfo->{hist_header104}) = Array2histcsv(\@Geno104);
				($hgeneinfo->{hist_value105}{$gid}, $hgeneinfo->{hist_header105}) = Array2histcsv(\@Geno105);
				
#				if($num_nonzero1 > 0 && $num_nodata1 == 0){
				if($num_nonzero1 > 0){
					$r101 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value101}{$gid}.",".join(",", @Geno101)."\n";			# disrupt
					$r101B .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value101B}{$gid}.",".join(",", @Geno101B)."\n";		# indel
					$hgeneclass->{transloc_val}{$gid} = 1;
					$cnt_geno101++;
					
					if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
						$r101M .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value101}{$gid}.",".join(",", @Geno101)."\n";
					}
				}
				if($num_nonzero2 > 0){
					$r102 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno102)."\n";
					$hgeneclass->{transloc_raw}{$gid} = 1;
					$cnt_geno102++;
				}
#				if($num_nonzero3 > 0 && $num_nodata3 == 0){
				if($num_nonzero3 > 0 && $failprot_refgenome == 0){
					$r103 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value103}{$gid}.",".join(",", @Geno103)."\n";
					$hgeneclass->{transloc_prot}{$gid} = 1;
					$cnt_geno103++;
				}
#				if($num_nonzero4 > 0 && $num_nodata4 == 0){
				if($num_nonzero4 > 0){
					$r104 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value104}{$gid}.",".join(",", @Geno104)."\n";
					$hgeneclass->{transloc_prom}{$gid} = 1;
					$cnt_geno104++;
				}
#				if($num_nonzero5 > 0 && $num_nodata5 == 0){
				if($num_nonzero5 > 0){
					$r105 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value105}{$gid}.",".join(",", @Geno105)."\n";
					$hgeneclass->{transloc_utr3}{$gid} = 1;
					$cnt_geno105++;
				}
				
				if($num_diffchr > 0){
					$cnt_diffchr++;
					
					($hgeneinfo->{hist_value201}{$gid}, $hgeneinfo->{hist_header201}) = Array2histcsv(\@Geno201);		# disrupt
					($hgeneinfo->{hist_value201B}{$gid}, $hgeneinfo->{hist_header201B}) = Array2histcsv(\@Geno201B);	# indel
					($hgeneinfo->{hist_value203}{$gid}, $hgeneinfo->{hist_header203}) = Array2histcsv(\@Geno203);
					($hgeneinfo->{hist_value204}{$gid}, $hgeneinfo->{hist_header204}) = Array2histcsv(\@Geno204);
					($hgeneinfo->{hist_value205}{$gid}, $hgeneinfo->{hist_header205}) = Array2histcsv(\@Geno205);
					
	#				if($num_nonzero1 > 0 && $num_nodata1 == 0){
					if($num_nonzero1 > 0){
						$r201 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value201}{$gid}.",".join(",", @Geno201)."\n";
						$r201B .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value201B}{$gid}.",".join(",", @Geno201B)."\n";
						$r201i .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value201}{$gid}.",".join(",", @Geno201i)."\n";
						$hgeneclass->{transloc_val}{$gid} = 1;
						$cnt_geno201++;
						
						if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
							$r201M .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value201}{$gid}.",".join(",", @Geno201)."\n";
						}
					}
					if($num_nonzero2 > 0){
						$r202 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno202)."\n";
						$hgeneclass->{transloc_raw}{$gid} = 1;
						$cnt_geno202++;
					}
	#				if($num_nonzero3 > 0 && $num_nodata3 == 0){
					if($num_nonzero3 > 0 && $failprot_refgenome == 0){
						$r203 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value203}{$gid}.",".join(",", @Geno203)."\n";
						$hgeneclass->{transloc_prot}{$gid} = 1;
						$cnt_geno203++;
					}
	#				if($num_nonzero4 > 0 && $num_nodata4 == 0){
					if($num_nonzero4 > 0){
						$r204 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value204}{$gid}.",".join(",", @Geno204)."\n";
						$hgeneclass->{transloc_prom}{$gid} = 1;
						$cnt_geno204++;
					}
	#				if($num_nonzero5 > 0 && $num_nodata5 == 0){
					if($num_nonzero5 > 0){
						$r205 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value205}{$gid}.",".join(",", @Geno205)."\n";
						$hgeneclass->{transloc_utr3}{$gid} = 1;
						$cnt_geno205++;
					}
				}
				else{
					if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
						$r201M .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,\n";
					}
				}
				
				if($num_trans > 0){
					$cnt_trans++;
					
					($hgeneinfo->{hist_value301}{$gid}, $hgeneinfo->{hist_header301}) = Array2histcsv(\@Geno301);		# disrupt
					($hgeneinfo->{hist_value301B}{$gid}, $hgeneinfo->{hist_header301B}) = Array2histcsv(\@Geno301B);	# indel
					($hgeneinfo->{hist_value303}{$gid}, $hgeneinfo->{hist_header303}) = Array2histcsv(\@Geno303);
					($hgeneinfo->{hist_value304}{$gid}, $hgeneinfo->{hist_header304}) = Array2histcsv(\@Geno304);
					($hgeneinfo->{hist_value305}{$gid}, $hgeneinfo->{hist_header305}) = Array2histcsv(\@Geno305);
					
	#				if($num_nonzero1 > 0 && $num_nodata1 == 0){
					if($num_nonzero1 > 0){
						$r301 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value301}{$gid}.",".join(",", @Geno301)."\n";
						$r301B .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value301B}{$gid}.",".join(",", @Geno301B)."\n";
						$hgeneclass->{transloc_val}{$gid} = 1;
						$cnt_geno301++;
						
						if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
							$r301M .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value301}{$gid}.",".join(",", @Geno301)."\n";
						}
					}
					if($num_nonzero2 > 0){
						$r302 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno302)."\n";
						$hgeneclass->{transloc_raw}{$gid} = 1;
						$cnt_geno302++;
					}
	#				if($num_nonzero3 > 0 && $num_nodata3 == 0){
					if($num_nonzero3 > 0 && $failprot_refgenome == 0){
						$r303 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value303}{$gid}.",".join(",", @Geno303)."\n";
						$hgeneclass->{transloc_prot}{$gid} = 1;
						$cnt_geno303++;
					}
	#				if($num_nonzero4 > 0 && $num_nodata4 == 0){
					if($num_nonzero4 > 0){
						$r304 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value304}{$gid}.",".join(",", @Geno304)."\n";
						$hgeneclass->{transloc_prom}{$gid} = 1;
						$cnt_geno304++;
					}
	#				if($num_nonzero5 > 0 && $num_nodata5 == 0){
					if($num_nonzero5 > 0){
						$r305 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value305}{$gid}.",".join(",", @Geno305)."\n";
						$hgeneclass->{transloc_utr3}{$gid} = 1;
						$cnt_geno305++;
					}
				}
				else{
					if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
						$r301M .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,\n";
					}
				}
			}
			else{
				if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
					$r101M .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,\n";
					$r201M .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,\n";
					$r301M .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,\n";
				}
			}
		}
		
#		if($num_nonzero1 > 0 && $num_nodata1 == 0){
		if($num_nonzero1 > 0){
			$r1 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1}{$gid}.",".join(",", @Geno1)."\n";
			$r1B .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1B}{$gid}.",".join(",", @Geno1B)."\n";
			$hgeneclass->{val}{$gid} = 1;
			$cnt_geno1++;
		}
		else{
			$cnt_geno1_rm++;
		}
		if($num_nonzero2 > 0){
			$r2 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno2)."\n";
			$hgeneclass->{raw}{$gid} = 1;
			$cnt_geno2++;
		}
#		if($num_nonzero3 > 0 && $num_nodata3 == 0){
		if($num_nonzero3 > 0 && $failprot_refgenome == 0){
			$r5 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value3}{$gid}.",".join(",", @Geno3)."\n";
			$hgeneclass->{prot}{$gid} = 1;
			$cnt_geno5++;
		}
#		if($num_nonzero4 > 0 && $num_nodata4 == 0){
		if($num_nonzero4 > 0){
			$r7 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value4}{$gid}.",".join(",", @Geno4)."\n";
			$hgeneclass->{prom}{$gid} = 1;
			$cnt_geno7++;
		}
#		if($num_nonzero5 > 0 && $num_nodata5 == 0){
		if($num_nonzero5 > 0){
			$r9 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value5}{$gid}.",".join(",", @Geno5)."\n";
			$hgeneclass->{utr3}{$gid} = 1;
			$cnt_geno9++;
		}
		
		if($num_conserved > 0 && $num_misgenotype == 0){		# double check of no data with $num_misgenotype == 0
			if($num_nonzero1 > 0 && $num_nodata1 == 0){
				$r11 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1}{$gid}.",".join(",", @Geno1)."\n";
				$r11_nmatrix .= $gid.",".join(",", @Geno1)."\n";
				$r11B .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1B}{$gid}.",".join(",", @Geno1B)."\n";
				$r11B_nmatrix .= $gid.",".join(",", @Geno1B)."\n";
				my $hrzval = Calc_normalizedvals(\@Geno1B);
				my $Geno1R = $hrzval->{R};
				my $Geno1N = $hrzval->{N};
				my $Geno1Z = $hrzval->{Z};
				my $Geno1L = $hrzval->{L};
				my $Geno1RZ = $hrzval->{RZ};
				my $Geno1RL = $hrzval->{RL};
				$r11R_nmatrix .= $gid.",".join(",", @{$Geno1R})."\n";
				$r11N_nmatrix .= $gid.",".join(",", @{$Geno1N})."\n";
				$r11Z_nmatrix .= $gid.",".join(",", @{$Geno1Z})."\n";
				$r11L_nmatrix .= $gid.",".join(",", @{$Geno1L})."\n";
				$r11RZ_nmatrix .= $gid.",".join(",", @{$Geno1RZ})."\n";
				$r11RL_nmatrix .= $gid.",".join(",", @{$Geno1RL})."\n";
				$hgeneclass->{val_1cnsv}{$gid} = 1;
				$cnt_geno11++;
			}
			else{
				$cnt_geno11_rm++;
			}
			if($num_nonzero2 > 0){
				$r12 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno2)."\n";
				$hgeneclass->{raw_1cnsv}{$gid} = 1;
				$cnt_geno12++;
			}
			if($num_nonzero3 > 0 && $num_nodata3 == 0 && $failprot_refgenome == 0){
				$r15 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value3}{$gid}.",".join(",", @Geno3)."\n";
				$r15_nmatrix .= $gid.",".join(",", @Geno3)."\n";
				$hgeneclass->{prot_1cnsv}{$gid} = 1;
				$cnt_geno15++;
			}
		}
		if($num_conserved_prom > 0 && $num_nonzero4 > 0 && $num_nodata4 == 0){
			$r17 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value4}{$gid}.",".join(",", @Geno4)."\n";
			$r17_nmatrix .= $gid.",".join(",", @Geno4)."\n";
			$hgeneclass->{prom_1cnsv}{$gid} = 1;
			$cnt_geno17++;
		}
		if($num_conserved_utr3 > 0 && $num_nonzero5 > 0 && $num_nodata5 == 0){
			$r19 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value5}{$gid}.",".join(",", @Geno5)."\n";
			$r19_nmatrix .= $gid.",".join(",", @Geno5)."\n";
			$hgeneclass->{utr3_1cnsv}{$gid} = 1;
			$cnt_geno19++;
		}
		
		if(! $hjudge->{false}){			# may be removed
			if($num_nonzero1 > 0 && $num_nodata1 == 0){
				$r3 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1}{$gid}.",".join(",", @Geno1)."\n";
				$hgeneclass->{val_allgeno}{$gid} = 1;
				$cnt_geno3++;
			}
			if($num_nonzero2 > 0){
				$r4 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno2)."\n";
				$hgeneclass->{raw_allgeno}{$gid} = 1;
				$cnt_geno4++;
			}
			if($num_nonzero3 > 0 && $num_nodata3 == 0 && $failprot_refgenome == 0){
				$r6 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value3}{$gid}.",".join(",", @Geno3)."\n";
				$hgeneclass->{prot_allgeno}{$gid} = 1;
				$cnt_geno6++;
			}
			if($num_nonzero4 > 0 && $num_nodata4 == 0){
				$r8 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value4}{$gid}.",".join(",", @Geno4)."\n";
				$hgeneclass->{prom_allgeno}{$gid} = 1;
				$cnt_geno8++;
			}
			if($num_nonzero5 > 0 && $num_nodata5 == 0){
				$r10 .= $gid.",".$hginfo->{$gid}{Chr}.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value5}{$gid}.",".join(",", @Geno5)."\n";
				$hgeneclass->{utr3_allgeno}{$gid} = 1;
				$cnt_geno10++;
			}
			$cnt1++;
		}
		$cnt0++;
	}
}

print "! [$cnt0] gene ID\n";
print "! SV data statistics\n";
print "  - all count\n";
print "    [$cnt_geno1] selected, [$cnt_geno1_rm] removed due to missing value in >=1 sample\n";
print "  - conserved in at least one sample\n";
print "    [$cnt_geno11] selected, [$cnt_geno11_rm] removed due to missing value in >=1 sample\n";

if($refgenome ne 'null'){
	print "! [$cnt_reffail] skipped as failed to detect in the analysis of reference genome (ref vs ref)\n";
	
	if($chrinfo ne 'null'){
		print "! [$cnt_diffchr] translocating to different chromosome\n";
		print "! [$cnt_trans] translocating to different region of the same chromosome\n";
	}
}

print "! saving data...\n";
my $r1_less500 = "";
my $r2_less500 = "";
my $r5_less500 = "";
my $r7_less500 = "";
my $r9_less500 = "";
my $r3_less500 = "";
my $r4_less500 = "";
my $r6_less500 = "";
my $r8_less500 = "";
my $r10_less500 = "";
my $r11_less500 = "";
my $r12_less500 = "";
my $r15_less500 = "";
my $r17_less500 = "";
my $r19_less500 = "";
($r1_less500, $hgeneclass) = Select_less500($r1, $hgeneclass, "val_less500");
($r2_less500, $hgeneclass) = Select_less500($r2, $hgeneclass, "raw_less500");
($r5_less500, $hgeneclass) = Select_less500($r5, $hgeneclass, "prot_less500");
($r7_less500, $hgeneclass) = Select_less500($r7, $hgeneclass, "promoter_less500");
($r9_less500, $hgeneclass) = Select_less500($r9, $hgeneclass, "utr3_less500");
#($r3_less500, $hgeneclass) = Select_less500($r3, $hgeneclass, "val_allgeno_less500");
#($r4_less500, $hgeneclass) = Select_less500($r4, $hgeneclass, "raw_allgeno_less500");
#($r6_less500, $hgeneclass) = Select_less500($r6, $hgeneclass, "prot_allgeno_less500");
#($r8_less500, $hgeneclass) = Select_less500($r8, $hgeneclass, "promoter_allgeno_less500");
#($r10_less500, $hgeneclass) = Select_less500($r10, $hgeneclass, "utr3_allgeno_less500");
($r11_less500, $hgeneclass) = Select_less500($r11, $hgeneclass, "val_1cnsv_less500");
($r12_less500, $hgeneclass) = Select_less500($r12, $hgeneclass, "raw_1cnsv_less500");
($r15_less500, $hgeneclass) = Select_less500($r15, $hgeneclass, "prot_1cnsv_less500");
($r17_less500, $hgeneclass) = Select_less500($r17, $hgeneclass, "promoter_1cnsv_less500");
($r19_less500, $hgeneclass) = Select_less500($r19, $hgeneclass, "utr3_1cnsv_less500");

my @Class = keys(%{$hgeneclass});
@Class = sort {$a cmp $b} @Class;

my $nsample = @Samples;
my $geneinfo = "gene_id,num genotyped with >10kb flanking,".join(",", @Class).",".$hgeneinfo->{hist_header1}."\n";
foreach my $gid (@GID2){
	my @Info;
	foreach my $class (@Class){
		if($hgeneclass->{$class}{$gid}){
			push(@Info, 1);
		}
		else{
			push(@Info, 0);
		}
	}
	$geneinfo .= $gid.",".$hgeneinfo->{num_true}{$gid}.",".join(",", @Info).",".$hgeneinfo->{hist_value1}{$gid}."\n";
}

my $dir = "combined_asm2sv";
unless(-e $dir){
	system("mkdir $dir");
}

my $rfile1 = "$dir/val_disrupt_q-".$nsample.".csv";
my $rfile1B = "$dir/val_indel_q-".$nsample.".csv";
my $rfile2 = "$dir/raw_q-".$nsample.".csv";
my $rfile3 = "$dir/val_disrupt_allgeno_q-".$nsample.".csv";
my $rfile4 = "$dir/raw_allgeno_q-".$nsample.".csv";
my $rfile5 = "$dir/prot_q-".$nsample.".csv";
my $rfile7 = "$dir/promoter_indel_q-".$nsample.".csv";
my $rfile9 = "$dir/3UTR_indel_q-".$nsample.".csv";
my $rfile6 = "$dir/prot_allgeno_q-".$nsample.".csv";
my $rfile8 = "$dir/promoter_indel_allgeno_q-".$nsample.".csv";
my $rfile10 = "$dir/3UTR_indel_allgeno_q-".$nsample.".csv";
my $rfile11 = "$dir/val_disrupt_1cnsv_q-".$nsample.".csv";
my $rfile11B = "$dir/val_indel_1cnsv_q-".$nsample.".csv";
my $rfile12 = "$dir/raw_1cnsv_q-".$nsample.".csv";
my $rfile15 = "$dir/prot_1cnsv_q-".$nsample.".csv";
my $rfile17 = "$dir/promoter_indel_1cnsv_q-".$nsample.".csv";
my $rfile19 = "$dir/3UTR_indel_1cnsv_q-".$nsample.".csv";
my $rfile1l = "$dir/val_disrupt_q-".$nsample."_less500.csv";
my $rfile2l = "$dir/raw_q-".$nsample."_less500.csv";
my $rfile5l = "$dir/prot_q-".$nsample."_less500.csv";
my $rfile7l = "$dir/promoter_indel_q-".$nsample."_less500.csv";
my $rfile9l = "$dir/3UTR_indel_q-".$nsample."_less500.csv";
my $rfile3l = "$dir/val_disrupt_allgeno_q-".$nsample."_less500.csv";
my $rfile4l = "$dir/raw_allgeno_q-".$nsample."_less500.csv";
my $rfile6l = "$dir/prot_allgeno_q-".$nsample."_less500.csv";
my $rfile8l = "$dir/promoter_indel_allgeno_q-".$nsample."_less500.csv";
my $rfile10l = "$dir/3UTR_indel_allgeno_q-".$nsample."_less500.csv";
my $rfile11l = "$dir/val_disrupt_1cnsv_q-".$nsample."_less500.csv";
my $rfile12l = "$dir/raw_1cnsv_q-".$nsample."_less500.csv";
my $rfile15l = "$dir/prot_1cnsv_q-".$nsample."_less500.csv";
my $rfile17l = "$dir/promoter_indel_1cnsv_q-".$nsample."_less500.csv";
my $rfile19l = "$dir/3UTR_indel_1cnsv_q-".$nsample."_less500.csv";
my $rginfo = "$dir/geneinfo_q-".$nsample.".csv";
my $rfile11_nmatrix = "$dir/val_disrupt_1cnsv_q-".$nsample."_num-matrix.csv";
my $rfile11B_nmatrix = "$dir/val_indel_1cnsv_q-".$nsample."_num-matrix.csv";
my $rfile11V_nmatrix = "$dir/forClust_indel_1cnsv_q-".$nsample."_native.csv";
my $rfile11R_nmatrix = "$dir/forClust_indel_1cnsv_q-".$nsample."_relative.csv";
my $rfile11N_nmatrix = "$dir/forClust_indel_1cnsv_q-".$nsample."_averaged.csv";
my $rfile11L_nmatrix = "$dir/forClust_indel_1cnsv_q-".$nsample."_log.csv";
my $rfile11Z_nmatrix = "$dir/forClust_indel_1cnsv_q-".$nsample."_zscore.csv";
my $rfile11RL_nmatrix = "$dir/forClust_indel_1cnsv_q-".$nsample."_log-relative.csv";
my $rfile11RZ_nmatrix = "$dir/forClust_indel_1cnsv_q-".$nsample."_zscore-relative.csv";
my $rfile15_nmatrix = "$dir/prot_1cnsv_q-".$nsample."_num-matrix.csv";
my $rfile17_nmatrix = "$dir/promoter_indel_1cnsv_q-".$nsample."_num-matrix.csv";
my $rfile19_nmatrix = "$dir/3UTR_indel_1cnsv_q-".$nsample."_num-matrix.csv";

SAVE($rfile1, $r1);
SAVE($rfile1B, $r1B);
SAVE($rfile2, $r2);
#SAVE($rfile3, $r3);
#SAVE($rfile4, $r4);
SAVE($rfile5, $r5);
#SAVE($rfile6, $r6);
SAVE($rfile7, $r7);
#SAVE($rfile8, $r8);
SAVE($rfile9, $r9);
#SAVE($rfile10, $r10);
SAVE($rfile11, $r11);
SAVE($rfile11B, $r11B);
SAVE($rfile12, $r12);
SAVE($rfile15, $r15);
SAVE($rfile17, $r17);
SAVE($rfile19, $r19);
SAVE($rfile1l, $r1_less500);
SAVE($rfile2l, $r2_less500);
SAVE($rfile5l, $r5_less500);
#SAVE($rfile3l, $r3_less500);
#SAVE($rfile4l, $r4_less500);
#SAVE($rfile6l, $r6_less500);
SAVE($rfile11l, $r11_less500);
SAVE($rfile12l, $r12_less500);
SAVE($rfile15l, $r15_less500);
SAVE($rfile17l, $r17_less500);
SAVE($rfile19l, $r19_less500);
SAVE($rfile11_nmatrix, $r11_nmatrix);
SAVE($rfile11V_nmatrix, $r11B_nmatrix);
SAVE($rfile11B_nmatrix, $r11B_nmatrix);
SAVE($rfile11R_nmatrix, $r11R_nmatrix);
#SAVE($rfile11N_nmatrix, $r11N_nmatrix);
SAVE($rfile11L_nmatrix, $r11L_nmatrix);
#SAVE($rfile11Z_nmatrix, $r11Z_nmatrix);
SAVE($rfile11RL_nmatrix, $r11RL_nmatrix);
#SAVE($rfile11RZ_nmatrix, $r11RZ_nmatrix);
SAVE($rfile15_nmatrix, $r15_nmatrix);
SAVE($rfile17_nmatrix, $r17_nmatrix);
SAVE($rfile19_nmatrix, $r19_nmatrix);
#SAVE($rginfo, $geneinfo);
SAVE($rginfo, $r_hist);

my $combined_predictorf = "$dir/combined_predict_orf_q-".$nsample.".fasta";
if(-e $combined_predictorf){
	system("rm $combined_predictorf");
}
foreach my $tmp (@{$AoF}){
	my $qorf = $tmp->[1];
	
	if(-e $qorf){
		system("cat $qorf >> $combined_predictorf");
	}
}

if($refgenome ne 'null'){
	if($cnt_diffchr > 0 || $cnt_trans > 0){
		my $rfile101 = "$dir/val_disrupt_transloc_q-".$nsample.".csv";
		my $rfile101M = "$dir/val_disrupt_1cnsv_transloc_q-".$nsample.".csv";
		my $rfile101B = "$dir/val_indel_transloc_q-".$nsample.".csv";
		my $rfile102 = "$dir/raw_transloc_q-".$nsample.".csv";
		my $rfile103 = "$dir/prot_transloc_q-".$nsample.".csv";
		my $rfile104 = "$dir/promoter_indel_transloc_q-".$nsample.".csv";
		my $rfile105 = "$dir/3UTR_indel_transloc_q-".$nsample.".csv";
		
		SAVE($rfile101, $r101);
		SAVE($rfile101M, $r101M);
		SAVE($rfile101B, $r101B);
		SAVE($rfile102, $r102);
		SAVE($rfile103, $r103);
		SAVE($rfile104, $r104);
		SAVE($rfile105, $r105);
		
		if($cnt_diffchr > 0){
			my $rfile201 = "$dir/val_disrupt_transloc-diffchr_q-".$nsample.".csv";
			my $rfile201M = "$dir/val_disrupt_1cnsv_transloc-diffchr_q-".$nsample.".csv";
			my $rfile201B = "$dir/val_indel_transloc-diffchr_q-".$nsample.".csv";
			my $rfile202 = "$dir/raw_transloc-diffchr_q-".$nsample.".csv";
			my $rfile203 = "$dir/prot_transloc-diffchr_q-".$nsample.".csv";
			my $rfile204 = "$dir/promoter_indel_transloc-diffchr_q-".$nsample.".csv";
			my $rfile205 = "$dir/3UTR_indel_transloc-diffchr_q-".$nsample.".csv";
			my $rfile201i = "$dir/info_indel_transloc-diffchr_q-".$nsample.".csv";
			
			SAVE($rfile201, $r201);
			SAVE($rfile201M, $r201M);
			SAVE($rfile201B, $r201B);
			SAVE($rfile202, $r202);
			SAVE($rfile203, $r203);
			SAVE($rfile204, $r204);
			SAVE($rfile205, $r205);
			SAVE($rfile201i, $r201i);
		}
		if($cnt_trans > 0){
			my $rfile301 = "$dir/val_disrupt_transloc-samechr_q-".$nsample.".csv";
			my $rfile301M = "$dir/val_disrupt_1cnsv_transloc-samechr_q-".$nsample.".csv";
			my $rfile301B = "$dir/val_indel_transloc-samechr_q-".$nsample.".csv";
			my $rfile302 = "$dir/raw_transloc-samechr_q-".$nsample.".csv";
			my $rfile303 = "$dir/prot_transloc-samechr_q-".$nsample.".csv";
			my $rfile304 = "$dir/promoter_indel_transloc-samechr_q-".$nsample.".csv";
			my $rfile305 = "$dir/3UTR_indel_transloc-samechr_q-".$nsample.".csv";
			
			SAVE($rfile301, $r301);
			SAVE($rfile301M, $r301M);
			SAVE($rfile301B, $r301B);
			SAVE($rfile302, $r302);
			SAVE($rfile303, $r303);
			SAVE($rfile304, $r304);
			SAVE($rfile305, $r305);
		}
	}
}

if($cnt_geno11 > 0){
	R_PCAplot($data_path, $rfile11V_nmatrix, "native");
	R_PCAplot($data_path, $rfile11R_nmatrix, "relative");
#	R_PCAplot($data_path, $rfile11N_nmatrix, "averaged");
	R_PCAplot($data_path, $rfile11L_nmatrix, "log2-transformed[native]");
#	R_PCAplot($data_path, $rfile11Z_nmatrix, "zscore");
	R_PCAplot($data_path, $rfile11RL_nmatrix, "log2-transformed[relative]");
#	R_PCAplot($data_path, $rfile11RZ_nmatrix, "zscore-relative");
}

if($zip){
	print "! --zip is invoked, zipping [$dir] ...\n";
	my $rzip = $dir.".zip";
	if(-e $rzip){
		system("rm $rzip");
	}
	system("zip -r $rzip $dir/ > /dev/null 2>&1");
}

END:{
#	print "\n! End of script.\n\n";
	my $end = 1;
}


################################################################################
#-------------------------------------------------------------------------------
sub R_PCAplot{
my $data_path = shift;
my $file0 = shift;
my $description = shift;

my $R_file0 = $file0;
$R_file0 =~ s/\.csv/_0\.txt/;

my $img_file0 = $file0;
my $img_file1 = $file0;
my $img_file2 = $file0;
my $img_file3 = $file0;
my $img_file4 = $file0;
$img_file0 =~ s/\.csv/_PCAplot\.png/;
$img_file1 =~ s/\.csv/_flushclust\.png/;
$img_file2 =~ s/\.csv/_flushclust_hang\.png/;
$img_file3 =~ s/\.csv/_hclust\.png/;
$img_file4 =~ s/\.csv/_hclust_hang\.png/;

my $script0 =<<EOS;
library(stats)
workingDir = "$data_path"
setwd(workingDir)
EOS

$script0 .=<<EOS;
#----------------------flashclust----------------------#
library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE)
Data <- read.csv("$file0", header=T, row.names=1)
datExpr <- as.data.frame(t(Data))
sampleTree <- flashClust(dist(datExpr), method = "complete")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
png("$img_file1", width=3000, height=1200)
plot(sampleTree, main = "Sample clustering, $description (flashClust)", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
png("$img_file2", width=3000, height=1200)
plot(sampleTree, hang=-1, main = "Sample clustering, $description (flashClust)", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
EOS

open(my $rfh0, ">", $R_file0) or die;
print $rfh0 $script0;
close $rfh0;

print "\n! generating clustering based on [$file0] ...\n";
if(system("Rscript --vanilla --slave $R_file0 > /dev/null 2>&1") == 0){
	if(-e $img_file1 && -e $img_file2){
		print "! output [$img_file1]\n";
		print "! output [$img_file2]\n";
	}
	else{
		print "! executed but missing image files...\n";
	}
}
else{
	print "! failed...\n";
}

#system("rm $R_file0");
if(-e "Rplots.pdf"){
	system("rm Rplots.pdf");
}
if(-e "Rplots1.pdf"){
	system("rm Rplots1.pdf");
}

}


#-------------------------------------------------------------------------------
sub Calc_normalizedvals{
my $A = shift;

my $max = 0;
my $avr = 0;
my $n = @{$A};
foreach my $val (@{$A}){
	if($max < $val){
		$max = $val;
	}
	$avr += $val;
}
$avr = $avr / $n;

my $sd = 0;
foreach my $val (@{$A}){
	$sd += ($val - $avr) * ($val - $avr);
}
$sd = sqrt($sd / $n);

my $R = [];
my $N = [];
my $Z = [];
my $L = [];
foreach my $val (@{$A}){
	my $rval = $val / $max;
	push(@{$R}, $rval);
	
	my $nval = $val / $avr;
	push(@{$N}, $nval);
	
	my $zval = 0;
	if($sd > 0){
		$zval = ($val - $avr) / $sd;
	}
	push(@{$Z}, $zval);
	
	my $logval = 0;
	if($val < 1 / 1024){
		$logval = -10;
	}
	elsif($val > 1024){
		$logval = 10;
	}
	else{
		$logval = log($val) / log(2);
	}
	push(@{$L}, $logval);
}

my $rmax = 0;
my $ravr = 0;
foreach my $val (@{$R}){
	if($rmax < $val){
		$rmax = $val;
	}
	$ravr += $val;
}
$ravr = $ravr / $n;

my $rsd = 0;
foreach my $val (@{$R}){
	$rsd += ($val - $ravr) * ($val - $ravr);
}
$rsd = sqrt($rsd / $n);

my $RZ = [];
my $RL = [];
foreach my $val (@{$R}){
	my $zval = 0;
	if($rsd > 0){
		$zval = ($val - $ravr) / $rsd;
	}
	push(@{$RZ}, $zval);
	
	my $logval = 0;
	if($val < 1 / 1024){
		$logval = -10;
	}
	elsif($val > 1024){
		$logval = 10;
	}
	else{
		$logval = log($val) / log(2);
	}
	push(@{$RL}, $logval);
}

my $rh = {};
$rh->{R} = $R;
$rh->{N} = $N;
$rh->{Z} = $Z;
$rh->{L} = $L;
$rh->{RZ} = $RZ;
$rh->{RL} = $RL;

return $rh;
}


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
	
	# example of chrinfo (tsv)
	#	prefix		convert_ID	original_ID
	#	Gmax_189	chr01		Gm01
	#	Gmax_189	chr02		Gm02
	#	Gmax_189	chr03		Gm03
	
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
sub Array2histcsv{
my $A = shift;

my $hash = {};
foreach my $val (@{$A}){
	for(my $i = 0; $i <= 1; $i += 0.3){
		my $p0 = $i;
		my $p1 = $i + 0.3;
		
		if($i < 1){
			if($p0 <= $val && $val < $p1){
				$hash->{$p0} += 1;
				last;
			}
		}
		else{
			if($p0 <= $val){
				$hash->{$p0} += 1;
				last;
			}
		}
	}
}

my @Hist;
my @Header;
for(my $i = 0; $i <= 1; $i += 0.3){
	my $p0 = $i;
	my $p1 = $i + 0.3;
	if($i < 1){
		push(@Header, "$p0 - $p1");
	}
	else{
		push(@Header, "$p0 -");
	}
	
	if($hash->{$p0}){
		push(@Hist, $hash->{$p0});
	}
	else{
		push(@Hist, 0);
	}
}

my $hist = join(",", @Hist);
my $header = join(",", @Header);

return ($hist, $header);
}


#-------------------------------------------------------------------------------
sub Array2histcsv2{
my $A = shift;

my $hash = {};
foreach my $val (@{$A}){
	for(my $i = 0; $i <= 2; $i += 0.6){
		my $p0 = $i;
		my $p1 = $i + 0.6;
		
		if($i < 1){
			if($p0 <= $val && $val < $p1){
				$hash->{$p0} += 1;
				last;
			}
		}
		else{
			if($p0 <= $val){
				$hash->{$p0} += 1;
				last;
			}
		}
	}
}

my @Hist;
my @Header;
for(my $i = 0; $i <= 2; $i += 0.6){
	my $p0 = $i;
	my $p1 = $i + 0.6;
	if($i < 1){
		push(@Header, "$p0 - $p1");
	}
	else{
		push(@Header, "$p0 -");
	}
	
	if($hash->{$p0}){
		push(@Hist, $hash->{$p0});
	}
	else{
		push(@Hist, 0);
	}
}

my $hist = join(",", @Hist);
my $header = join(",", @Header);

return ($hist, $header);
}


#-------------------------------------------------------------------------------
sub Select_less500{
my $lines = shift;
my $hgeneclass = shift;
my $str = shift;

my @L = split(/\n/, $lines);
my $cnt = 0;
my $AoA = [];
my $sorted = "";
foreach my $l (@L){
	if($cnt == 0){
		$sorted .= $l."\n";
	}
	else{
		my @A = split(/,/, $l);
		my $numA = @A;
		my $sum = 0;
		for(my $i = 9; $i < $numA; $i++){
			$sum += $A[$i];
		}
		my @tmp;
		push(@tmp, $sum);
		push(@tmp, \@A);
		push(@tmp, $A[0]);
		push(@{$AoA}, \@tmp);
	}
	$cnt++;
}

@{$AoA} = sort {$a->[0] <=> $b->[0]} @{$AoA};

for(my $j = 0; $j < 500; $j++){
	if($AoA->[$j][1]){
		my $A = $AoA->[$j][1];
		my $gid = $AoA->[$j][2];
		$sorted .= join(",", @{$A})."\n";
		$hgeneclass->{$str}{$gid} = 1;
	}
}

return ($sorted, $hgeneclass);
}


#-------------------------------------------------------------------------------
sub Collectdata{
my $file = shift;
my $hash = shift;
my $refgenome = shift;
my $hseqinfo = shift;

print " reading [$file]...\n";
if(-e $file){
	open(my $fh, "<", $file) or die;
	my $cnt = 0;
	my $neighbor_tlen = "null";
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		unless($line){
			next;
		}
		if($cnt == 0){
			$hash->{header} = $line;
			
			my @A = split(/\t/, $line);
			foreach my $val (@A){
				if($val =~ /\(/ && $val =~ /bp\)/){
					my @tmp0 = split(/\(/, $val);
					if($tmp0[1]){
						my @tmp1 = split(/\)/, $tmp0[1]);
						if($tmp1[0]){
							$tmp1[0] =~ s/bp//;
							$tmp1[0] =~ s/\s//g;
							if($tmp1[0] =~ /\d/ && $tmp1[0] !~ /\D/){
								$neighbor_tlen = $tmp1[0];
							}
						}
					}
				}
			}
		}
		else{
			my @A = split(/\t/, $line);
			
			my $gid = $A[23];
	#		my $sv = $A[21];				# value representing the degree of disruption: [0] < deletion or insertion < [1.0]
			my $aln_ratio = $A[19];
			my $ins_ratio = $A[20];
			my $promoter = $A[25];
			my $utr3 = $A[27];
			my $id = $A[7];
			
			if($A[32] !~ /true/i){			# skip non-coding transcript
				$hash->{geno_including_flanking}{$id}{$gid} = "false";
				next;
			}
			
			if($promoter =~ /\;/){
				# ratio_hit; ratio_span; rev_ratio_span; P2Gdist; span_pos0; span_pos1
				my @tmpA = split(/\;/, $promoter);
				if($tmpA[2] ne 'null'){
					$promoter = $tmpA[2];
				}
				elsif($tmpA[1] ne 'null'){
					$promoter = $tmpA[1];
				}
				else{
					$promoter = $tmpA[0];
				}
			}
			if($utr3 =~ /\;/){
				# ratio_hit; ratio_span; rev_ratio_span; P2Gdist; span_pos0; span_pos1
				my @tmpA = split(/\;/, $utr3);
				if($tmpA[2] ne 'null'){
					$utr3 = $tmpA[2];
				}
				elsif($tmpA[1] ne 'null'){
					$utr3 = $tmpA[1];
				}
				else{
					$utr3 = $tmpA[0];
				}
			}
			
	#		my $sv2 = $sv;					# value representing deletion and insertion separately: [0] < deletion < [1.0] < insertion < [2.0]
			my $dbspan_woN = $A[14];
			my $hitaln = $A[15];
			my $hitspan_woN = $A[17];
			my $bpinsert = $A[18];
			
			my $sv = "-";
			my $sv2 = "-";
			my $sv3 = "-";
			if(! defined $aln_ratio || $aln_ratio eq '-' || $hitspan_woN eq '-'){
				$sv = 0;
				$sv2 = 0;
				$sv3 = 0;
			}
			else{
				my $aln_ratio = $hitaln / $dbspan_woN;
				my $alnspan_ratio = $hitspan_woN / $dbspan_woN;
				
				if(defined $bpinsert && $bpinsert ne '-'){
					my $ins_rvratio = $dbspan_woN / $hitspan_woN;
					if($alnspan_ratio > $ins_rvratio){
						$sv2 = $hitspan_woN / $dbspan_woN;
						$sv3 = $sv2;
					}
					else{
						$sv2 = $alnspan_ratio;
						$sv3 = $aln_ratio;
					}
					$sv = $alnspan_ratio - $bpinsert / $dbspan_woN;
				}
				else{
					$sv = $alnspan_ratio;
					$sv2 = $alnspan_ratio;
					$sv3 = $aln_ratio;
				}
			}
			
			if($refgenome ne 'null' && $refgenome eq $id){
				$hash->{refchr}{$id}{$gid} = $A[2];
				$hash->{refpos0}{$id}{$gid} = $A[5];
				$hash->{refpos1}{$id}{$gid} = $A[6];
			}
			
			$hash->{sv}{$id}{$gid} = $sv;
			$hash->{sv2}{$id}{$gid} = $sv2;
			$hash->{svprom}{$id}{$gid} = $promoter;
			$hash->{svutr3}{$id}{$gid} = $utr3;
#			$hash->{svprot}{$id}{$gid} = $sv;
			
			$hash->{svjudge1}{$id}{$gid} = $A[22];					# native judge
			$hash->{svjudge2}{$id}{$gid} = "null";
			if($A[22] !~ /missing border/ && $A[22] !~ /absent/ && $A[22] !~ /NNN/ && $A[22] !~ /no hit but some/){
				if($A[22] =~ /present/ || $A[22] =~ /insertion/){
					$hash->{svjudge2}{$id}{$gid} = "genotyped";		# judge classified
				}
			}
			
			if($A[20] ne '-' && $A[20] ne '0' && $A[22] =~ /insertion/){
				$hash->{raw}{$id}{$gid} = 1 / $A[20];
			}
			else{
				$hash->{raw}{$id}{$gid} = $sv;
			}
			
			my $rprot = 0;
			if($A[34] && $A[34] ne '-' && $A[34] && $A[35] ne '-' && $A[35] ne '.'){
#				$rprot = $A[34] / $A[35];
				my @DBplen = split(/,/, $A[33]);
				my @DBplenr;
				foreach my $val (@DBplen){
					if(defined $val && $val =~ /\d/ && $val !~ /\D/ && $val > 0){
						push(@DBplenr, $val);
					}
				}
				@DBplenr = sort {$b <=> $a} @DBplenr;
				
				if($DBplenr[0]){
					my $protlen_db = $DBplenr[0];
					my $protlen_q = $A[34];
					my $alnlen = $A[35];
					my $pcnt_similarity = $A[36];
					my $pcnt_gap = $A[37];
					
					my $dbprot_factor = log($protlen_db) / log(10);
					my $difflen_ratio = abs($protlen_db - $protlen_q) / $protlen_db;
					
					$rprot = ($pcnt_similarity - $pcnt_gap) / 100 - ($difflen_ratio ** $dbprot_factor);
					if($rprot < 0){
						$rprot = 0;
					}
				}
			}
			$hash->{rprot}{$id}{$gid} = $rprot;
			$hash->{svprot}{$id}{$gid} = $rprot;
			
#			if(defined $sv && $sv ne '-' && $rprot > $sv){
#				$hash->{svprot}{$id}{$gid} = $rprot;
#			}
			
			if($refgenome ne 'null' && $hseqinfo->{$id}){
				if(defined $A[8] && defined $A[9] && $A[9] ne '-' && $A[9] ne 'null' && $A[9] !~ /\D/ && $A[9] =~ /\d/ && defined $A[10] && $A[10] ne '-' && $A[10] ne 'null' && $A[10] !~ /\D/ && $A[10] =~ /\d/){
					if($A[8] !~ /scaffold/i){
						if($hseqinfo->{$id}{$A[8]}){
							$hash->{hconvid}{$id}{$gid} = $hseqinfo->{$id}{$A[8]};
						}
						else{
							print "! missing seq ID alias for [$A[8]] of [$id]\n";
							$hash->{err} = 1;
							return $hash;
						}
					}
					else{
						if($hseqinfo->{$id}{$A[8]}){
							$hash->{hconvid}{$id}{$gid} = $hseqinfo->{$id}{$A[8]};
						}
						else{
							$hash->{hconvid}{$id}{$gid} = "-";
						}
					}
					$hash->{hseqid}{$id}{$gid} = $A[8];
					$hash->{hpos0}{$id}{$gid} = $A[9];
					$hash->{hpos1}{$id}{$gid} = $A[10];
				}
				else{
					$hash->{hconvid}{$id}{$gid} = "-";
					$hash->{hseqid}{$id}{$gid} = "-";
					$hash->{hpos0}{$id}{$gid} = "-";
					$hash->{hpos1}{$id}{$gid} = "-";
				}
			}
			
			my $geno_including_flanking = "false";
			if($A[39] && $A[39] ne '-' && $A[39] >= 5000 && $sv && $sv ne '-' && $sv > 0){
				$geno_including_flanking = "true";
			}
			$hash->{geno_including_flanking}{$id}{$gid} = $geno_including_flanking;
			
			if(! $hash->{ginfo}{$gid}){
				$hash->{ginfo}{$gid}{Chr} = $A[2];
				$hash->{ginfo}{$gid}{pos0} = $A[3];
				$hash->{ginfo}{$gid}{pos1} = $A[4];
				$hash->{ginfo}{$gid}{sp0} = $A[5];
				$hash->{ginfo}{$gid}{sp1} = $A[6];
				$hash->{ginfo}{$gid}{flen} = $A[39];
				$hash->{ginfo}{$gid}{line} = $cnt;
			}
			else{
				my $judge = 0;
				if($hash->{ginfo}{$gid}{Chr} ne $A[2]){
					print "! error at line [$cnt] : distinct dbfasta for [$gid] : [$hash->{ginfo}{$gid}{Chr} (line $hash->{ginfo}{$gid}{line}) ne $A[2]]\n";
					$judge++;
				}
#				if($hash->{ginfo}{$gid}{pos0} ne $A[3]){
#					print "! error at line [$cnt] : distinct pos0 for [$gid] : [$hash->{ginfo}{$gid}{pos0} (line $hash->{ginfo}{$gid}{line}) ne $A[3]]\n";
#					$judge++;
#				}
#				if($hash->{ginfo}{$gid}{pos1} ne $A[4]){
#					print "! error at line [$cnt] : distinct pos1 for [$gid] : [$hash->{ginfo}{$gid}{pos1} (line $hash->{ginfo}{$gid}{line}) ne $A[4]]\n";
#					$judge++;
#				}
				if($hash->{ginfo}{$gid}{sp0} ne $A[5]){
					print "! error at line [$cnt] : distinct sp0 for [$gid] : [$hash->{ginfo}{$gid}{sp0} (line $hash->{ginfo}{$gid}{line}) ne $A[5]]\n";
					$judge++;
				}
				if($hash->{ginfo}{$gid}{sp1} ne $A[6]){
					print "! error at line [$cnt] : distinct sp1 for [$gid] : [$hash->{ginfo}{$gid}{sp1} (line $hash->{ginfo}{$gid}{line}) ne $A[6]]\n";
					$judge++;
				}
				if($judge > 0){
					die;
				}
			}
			
			if(! $hash->{dbfasta}){
				$hash->{dbfasta} = $A[1];
			}
#			elsif($hash->{dbfasta} ne $A[1]){
#				print "! error : distinct dbfasta [$hash->{dbfasta} ne $A[1]]\n";
#				die;
#			}
		}
		$cnt++;
	}
	close $fh;
}

return $hash;
}


#-------------------------------------------------------------------------------
sub SearchSummary{
my ($kw1, $kw2) = @_;

my $tmp = "_tmp_searchsummary.txt";
system("find . | grep rev_summary_genome2pav_ | grep -v old | grep -v plot | grep -v original | grep .tsv > $tmp");

print "\n! searching for asm2pav summary...\n";
open(my $fh, "<", $tmp) or die;
my @F;
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($kw1 && $kw2){
		if($line =~ /$kw1/i && $line !~ /$kw2/i){
			push(@F, $line);
			print " $line\n";
			$cnt++;
		}
	}
	elsif($kw1 && ! $kw2){
		if($line =~ /$kw1/i){
			push(@F, $line);
			print " $line\n";
			$cnt++;
		}
	}
	elsif(! $kw1 && $kw2){
		if($line !~ /$kw2/i){
			push(@F, $line);
			print " $line\n";
			$cnt++;
		}
	}
	else{
		push(@F, $line);
		print " $line\n";
		$cnt++;
	}
}
close $fh;

system("rm $tmp");

print "! [$cnt] files found\n";

return \@F;
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


#----------------------------------------------------------
sub APPEND{
my $file = shift;
my $str = shift;

open(my $fh, ">>", $file) or die;
print $fh $str;
close $fh;

#print "! output [$file]\n";

}


#-------------------------------------------------------------------------------
sub Delete{
my $file = shift;

if(-e $file){
	system("rm $file");
}

}





