#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use FindBin;

my $bin_gth = shift;
my $bin_miniprot = shift;
my $bin_samtools = shift;
my $bin_gffread = shift;
my $bin_blat = shift;
my $bin_blat2hint = shift;
my $bin_matcher = shift;
my $gid = shift;
my $tmp_dbprot = shift;
my $tmp_dbtranscript = shift;
my $tmp_dbcds = shift;
my $tmp_qfasta = shift;
my $tmp_qgff = shift;
my $tmp_qprot = shift;
my $tmp_align = shift;
my $tmpdir = shift;
my $qpref = shift;
my $qpos0 = shift;
my $combined_qprot = shift;
my $combined_qgff = shift;

my $script_path = $FindBin::Bin;
my $bin_FindCDS_mRNA = $script_path."/FindCDS_mRNA.pl";
my $bin_gfftoprot = $script_path."/gfftoprot.pl";

if(! $bin_FindCDS_mRNA || ! -e $bin_FindCDS_mRNA){
	print "! missing FindCDS_mRNA.pl ...\n";
	goto END;
}
if(! $bin_gfftoprot || ! -e $bin_gfftoprot){
	print "! missing gfftoprot.pl ...\n";
	goto END;
}

if(! $bin_gth || ! -e $bin_gth){
	print "! missing gth...\n";
	goto END;
}
if(! $bin_miniprot || ! -e $bin_miniprot){
	print "! missing miniprot...\n";
	goto END;
}
if(! $bin_samtools || ! -e $bin_samtools){
	print "! missing samtools...\n";
	goto END;
}
if(! $bin_gffread || ! -e $bin_gffread){
	print "! missing gffread...\n";
	goto END;
}
if(! $bin_blat || ! -e $bin_blat){
	print "! missing blat...\n";
	goto END;
}
if(! $bin_matcher || ! -e $bin_matcher){
	print "! missing matcher...\n";
	goto END;
}
if(! $tmp_dbprot || ! -e $tmp_dbprot){
	print "! missing dbprot\n";
	goto END;
}
if(! $tmp_dbtranscript || ! -e $tmp_dbtranscript){
	print "! missing dbtranscript\n";
	goto END;
}
if(! $tmp_dbcds || ! -e $tmp_dbcds){
	print "! missing dbcds\n";
	goto END;
}
if(! $tmp_qgff){
	print "! missing qgff file name\n";
	goto END;
}
if(! $tmp_qprot){
	print "! missing qprot file name\n";
	goto END;
}
if(! $tmp_align){
	print "! missing matcher align file name\n";
	goto END;
}
if(! $tmpdir || ! -e $tmpdir){
	print "! missing temp directory\n";
	goto END;
}

if($gid){
	print "\n! findorf.pl | gid=[$gid]\n";
}
if(-e $combined_qprot){
	my $judge_done = Search_previous($combined_qprot, $gid, $qpref);
	
	if($judge_done && $judge_done eq 'true'){
		print "! already done for [$gid], skip...\n";
		goto END;
	}
}

my $cmd = "";
my $minIdentity = 95;
if(-e $tmp_dbprot && -e $tmp_qfasta){
	my $judge = 0;
	my $AoHit = [];
	
	if(! -e $tmp_align){
		print "\n! miniprot...\n";
		my $method2 = "miniprot";
		my @CMD2;
		$CMD2[0] = "$bin_miniprot --gff $tmp_qfasta $tmp_dbprot > $tmp_qgff 2> /dev/null";
#		$CMD2[1] = "$bin_gffread $tmp_qgff -g $tmp_qfasta -y $tmp_qprot";
		$CMD2[1] = "$bin_gfftoprot $tmp_qfasta $tmp_qgff $tmp_qprot";
		
		foreach my $cmd (@CMD2){
			print "! cmd=[$cmd]\n";
			system($cmd);
		}
		Rm_ifempty($tmp_qprot);
		
		if(-e $tmp_qprot){
			my $AoR = Revseq2array($tmp_qprot);
			
			my $num_seq = 0;
			foreach my $rtmp (@{$AoR}){
				$num_seq++;
				print "! [$method2] [$num_seq] | start...\n";
				
				my $cand_fasta = $rtmp->[0];
				my $split_status = $rtmp->[1];
				SAVE($tmp_qprot, $cand_fasta);
				
				if(! -e $tmp_qprot){
					next;
				}
				
				Rm($tmp_qgff);
				my @CMD2R;
				$CMD2R[0] = "$bin_miniprot --gff $tmp_qfasta $tmp_qprot > $tmp_qgff 2> /dev/null";
#				$CMD2R[1] = "$bin_gffread $tmp_qgff -g $tmp_qfasta -y $tmp_qprot";
				$CMD2R[1] = "$bin_gfftoprot $tmp_qfasta $tmp_qgff $tmp_qprot";
				
				print "! [$method2] [$num_seq] (each redo)\n";
				foreach my $cmd (@CMD2R){
					print "! cmd=[$cmd]\n";
					system($cmd);
				}
				
				my $num_eachseq = 0;
				if(-e $tmp_qprot){
					$num_eachseq = CheckFasta($tmp_qprot);
				}
				if($num_eachseq == 0){
					print "! [$method2] [$num_seq] | no sequence in the file, skip...\n";
					Rm($tmp_qprot);
					next;
				}
				
				my $cmd1 = "$bin_matcher -asequence $tmp_dbprot -bsequence $tmp_qprot -outfile $tmp_align";
				print "! [$method2] [$num_seq] | cmd=[$cmd1]\n";
				system($cmd1);
				
				if(-e $tmp_align){
					my $hmatcher = Read_matcher($tmp_align);
					if($hmatcher->{score} && $hmatcher->{score} ne '-'){
						my $hitseq = PickupFasta($tmp_qprot, $gid, $qpref, $method2);
						
						my ($gff_lines, $gff_native, $gff_transcript) = ReadGff($tmp_qgff, $gid, $qpref, $hitseq->[0], $hitseq->[5], $method2, $qpos0);
						$hitseq->[0] = $qpref."_".$gid;
#						ADD("_gff_native.gff", $gff_native);
						
						my @tmp;
						push(@tmp, $hitseq->[0]);		# revID
						push(@tmp, $hitseq->[1]);		# seq
						push(@tmp, $method2);
						push(@tmp, $hmatcher->{score} * $hitseq->[2]);
						push(@tmp, $hmatcher->{native});
						push(@tmp, $gff_lines);
						push(@tmp, $gff_transcript);
						push(@tmp, $hmatcher->{identity});
						push(@{$AoHit}, \@tmp);
						
						print "! [$method2] [$num_seq] | [$hitseq->[0]] score=[$hmatcher->{score}]\n";
					}
				}
				
#				Rm($tmp_qprot);
				Rm($tmp_align);
			}
		}
	}
	
	if(! -e $tmp_align){
		print "\n! miniprot (extended)...\n";
		my $method2 = "miniprot";
		my @CMD2;
		$CMD2[0] = "$bin_miniprot -O 0 -E 0 -J 0 -L 6 -B 50 --gff $tmp_qfasta $tmp_dbprot > $tmp_qgff 2> /dev/null";
#		$CMD2[1] =  "$bin_gffread $tmp_qgff -g $tmp_qfasta -y $tmp_qprot";
		$CMD2[1] =  "$bin_gfftoprot $tmp_qfasta $tmp_qgff $tmp_qprot";
		
		foreach my $cmd (@CMD2){
			print "! cmd=[$cmd]\n";
			system($cmd);
		}
		Rm_ifempty($tmp_qprot);
		
		if(-e $tmp_qprot){
			my $AoR = Revseq2array($tmp_qprot);
			
			my $num_seq = 0;
			foreach my $rtmp (@{$AoR}){
				$num_seq++;
				print "! [$method2] [$num_seq] | start...\n";
				
				my $cand_fasta = $rtmp->[0];
				my $split_status = $rtmp->[1];
				SAVE($tmp_qprot, $cand_fasta);
				
				if(! -e $tmp_qprot){
					next;
				}
				
				Rm($tmp_qgff);
				my @CMD2R;
				$CMD2R[0] = "$bin_miniprot -O 0 -E 0 -J 0 -L 6 -B 50 --gff $tmp_qfasta $tmp_qprot > $tmp_qgff 2> /dev/null";
#				$CMD2R[1] = "$bin_gffread $tmp_qgff -g $tmp_qfasta -y $tmp_qprot";
				$CMD2R[1] = "$bin_gfftoprot $tmp_qfasta $tmp_qgff $tmp_qprot";
				
				print "! [$method2] [$num_seq] (each redo)\n";
				foreach my $cmd (@CMD2R){
					print "! cmd=[$cmd]\n";
					system($cmd);
				}
				
				my $num_eachseq = 0;
				if(-e $tmp_qprot){
					$num_eachseq = CheckFasta($tmp_qprot);
				}
				if($num_eachseq == 0){
					print "! [$method2] [$num_seq] | no sequence in the file, skip...\n";
					Rm($tmp_qprot);
					next;
				}
				
				my $cmd1 = "$bin_matcher -asequence $tmp_dbprot -bsequence $tmp_qprot -outfile $tmp_align";
				print "! [$method2] [$num_seq] | cmd=[$cmd1]\n";
				system($cmd1);
				
				if(-e $tmp_align){
					my $hmatcher = Read_matcher($tmp_align);
					if($hmatcher->{score} && $hmatcher->{score} ne '-'){
						my $hitseq = PickupFasta($tmp_qprot, $gid, $qpref, $method2);
						
						my ($gff_lines, $gff_native, $gff_transcript) = ReadGff($tmp_qgff, $gid, $qpref, $hitseq->[0], $hitseq->[5], $method2, $qpos0);
						$hitseq->[0] = $qpref."_".$gid;
#						ADD("_gff_native.gff", $gff_native);
						
						my @tmp;
						push(@tmp, $hitseq->[0]);		# revID
						push(@tmp, $hitseq->[1]);		# seq
						push(@tmp, $method2);
						push(@tmp, $hmatcher->{score} * $hitseq->[2]);
						push(@tmp, $hmatcher->{native});
						push(@tmp, $gff_lines);
						push(@tmp, $gff_transcript);
						push(@tmp, $hmatcher->{identity});
						push(@{$AoHit}, \@tmp);
						
						print "! [$method2] [$num_seq] | [$hitseq->[0]] score=[$hmatcher->{score}]\n";
					}
				}
				
#				Rm($tmp_qprot);
				Rm($tmp_align);
			}
		}
	}
	
	if(@{$AoHit}){
		@{$AoHit} = sort {$b->[3] <=> $a->[3]} @{$AoHit};
#		@{$AoHit} = sort {$b->[7] <=> $a->[7]} @{$AoHit};
		
		my $num_hit = @{$AoHit};
		my $nrseq = {};
		my $iskip = 0;
		for(my $i = 0; $i < 3; $i++){
			unless($AoHit->[$i]){
				next;
			}
			
			my $candhit = $AoHit->[$i];
			
			if($nrseq->{$candhit->[1]}){
				$iskip++;
				next;
			}
			$nrseq->{$candhit->[1]} = 1;
			
			my $j = $i + 1 - $iskip;
			$candhit->[0] .= ".t".$j;
			
			print "! [$i] hit = [$candhit->[0]] [$candhit->[2]]\n";
			if($i == 0){
				my $candhit_fasta = ">".$candhit->[0]."\n".$candhit->[1]."\n";
				SAVE($tmp_qprot, $candhit_fasta);
				SAVE($tmp_align, $candhit->[4]);
				
				BindFasta($combined_qprot, $qpref, $gid, $candhit->[0], $candhit->[1], $candhit->[2]);
				BindGff($combined_qgff, $candhit->[5]);
			}
			else{
				my $nisoform = ".t".$j;
				$candhit->[6] =~ s/\.t1/$nisoform/g;
				
				BindFasta($combined_qprot, $qpref, $gid, $candhit->[0], $candhit->[1], $candhit->[2]);
				BindGff($combined_qgff, $candhit->[6]);
			}
		}
		$judge = 1;
	}
	
	if($judge eq '0'){
		Rm($tmp_align);
		
		print "\n! trying genome threader...\n";
		my @CMD1;
		$CMD1[0] = "$bin_gth -force -genomic $tmp_qfasta -protein $tmp_dbprot -gff3out -skipalignmentout -o $tmp_qgff > /dev/null 2>&1";
#		$CMD1[1] = "$bin_gffread $tmp_qgff -g $tmp_qfasta -y $tmp_qprot";
		$CMD1[1] = "$bin_gfftoprot $tmp_qfasta $tmp_qgff $tmp_qprot";
		my $method1 = "gth";
		
		foreach my $cmd (@CMD1){
			print "! cmd=[$cmd]\n";
			system($cmd);
		}
		Rm_ifempty($tmp_qprot);
		
		if(-e $tmp_qprot){
			my $AoR = Revseq2array($tmp_qprot);
			
			my $num_seq = 0;
			foreach my $rtmp (@{$AoR}){
				$num_seq++;
				print "! [$method1] [$num_seq] | start...\n";
				
				my $cand_fasta = $rtmp->[0];
				my $split_status = $rtmp->[1];
				SAVE($tmp_qprot, $cand_fasta);
				
				if(! -e $tmp_qprot){
					next;
				}
				
				Rm($tmp_qgff);
				my @CMD2R;
				$CMD2R[0] = "$bin_gth -force -genomic $tmp_qfasta -protein $tmp_qprot -gff3out -skipalignmentout -o $tmp_qgff > /dev/null 2>&1";
#				$CMD2R[1] = "$bin_gffread $tmp_qgff -g $tmp_qfasta -y $tmp_qprot";
				$CMD2R[1] = "$bin_gfftoprot $tmp_qfasta $tmp_qgff $tmp_qprot";
				
				print "! [$method1] [$num_seq] (each redo)\n";
				foreach my $cmd (@CMD2R){
					print "! cmd=[$cmd]\n";
					system($cmd);
				}
				
				my $num_eachseq = 0;
				if(-e $tmp_qprot){
					$num_eachseq = CheckFasta($tmp_qprot);
				}
				if($num_eachseq == 0){
					print "! [$method1] [$num_seq] | no sequence in the file, skip...\n";
					Rm($tmp_qprot);
					next;
				}
				
				my $cmd1 = "$bin_matcher -asequence $tmp_dbprot -bsequence $tmp_qprot -outfile $tmp_align";
				print "! [$method1] [$num_seq] | cmd=[$cmd1]\n";
				system($cmd1);
				
				if(-e $tmp_align){
					my $hmatcher = Read_matcher($tmp_align);
					if($hmatcher->{score} && $hmatcher->{score} ne '-'){
						my $hitseq = PickupFasta($tmp_qprot, $gid, $qpref, $method1);
						
						my ($gff_lines, $gff_native) = ReadGff($tmp_qgff, $gid, $qpref, $hitseq->[0], $hitseq->[5], $method1, $qpos0);
						$hitseq->[0] = $qpref."_".$gid.".t1";
						
						my @tmp;
						push(@tmp, $hitseq->[0]);		# revID
						push(@tmp, $hitseq->[1]);		# seq
						push(@tmp, $method1);
						push(@tmp, $hmatcher->{score} * $hitseq->[2]);
						push(@tmp, $hmatcher->{native});
						push(@tmp, $gff_lines);
						push(@tmp, "null");
						push(@tmp, $hmatcher->{identity});
						push(@{$AoHit}, \@tmp);
						
						print "! [$method1] [$num_seq] | [$hitseq->[0]] score=[$hmatcher->{score}]\n";
					}
				}
				
#				Rm($tmp_qprot);
				Rm($tmp_align);
			}
		}
		
		if(@{$AoHit}){
			@{$AoHit} = sort {$b->[3] <=> $a->[3]} @{$AoHit};
#			@{$AoHit} = sort {$b->[7] <=> $a->[7]} @{$AoHit};
			my $tophit = $AoHit->[0];
			
			print "! tophit = [$tophit->[0]] [$tophit->[2]]\n";
			my $tophit_fasta = ">".$tophit->[0]."\n".$tophit->[1]."\n";
			SAVE($tmp_qprot, $tophit_fasta);
			SAVE($tmp_align, $tophit->[4]);
			
			BindFasta($combined_qprot, $qpref, $gid, $tophit->[0], $tophit->[1], $tophit->[2]);
			BindGff($combined_qgff, $tophit->[5]);
			$judge = 1;
		}
	}
	
	if($judge eq '0'){
		Rm($tmp_align);
		
		print "\n! trying blat...\n";
		my $psl_transcript = "$tmpdir/_qtranscript_".$gid.".psl";
		my $egtf_transcript = "$tmpdir/_qtranscript_".$gid.".E.gtf";
		my $egff_transcript = "$tmpdir/_qtranscript_".$gid.".E.gff";
		my $eprot_transcript = "$tmpdir/_qtranscript_".$gid.".E.fasta";
		my $prefix = $gid.".E.transcript";
		my $method3 = "manual";
		
		Rm($egtf_transcript);
		Rm($egff_transcript);
		
		my @CMD1;
		$CMD1[0] = "$bin_blat -t=dna -q=rna -noHead -minIdentity=$minIdentity -maxIntron=10000 -minScore=30 $tmp_qfasta $tmp_dbtranscript $psl_transcript";
		$CMD1[1] = "$bin_blat2hint --in=$psl_transcript --out=$egtf_transcript --ep_cutoff=0";
		
		foreach my $cmd (@CMD1){
			print "! cmd=[$cmd]\n";
			system($cmd);
		}
		
		my $strand = Select_tophit($psl_transcript, $psl_transcript);
		if(defined $strand){
			my $blat_status = Add_feature_info($egtf_transcript, $strand, $egff_transcript);
			if($blat_status ne 'null'){
				my @CMD2;
				$CMD2[0] = "$bin_gffread $egff_transcript -g $tmp_qfasta -w $eprot_transcript";
				foreach my $cmd (@CMD2){
					print "! cmd=[$cmd]\n";
					system($cmd);
				}
				
				my $rh = Getseqinfo($eprot_transcript);
				my $hseq = $rh->{data};
				my @SeqID = keys(%{$hseq});
				@SeqID = sort {$a cmp $b} @SeqID;
				
				foreach my $sid (@SeqID){
					my $cmd = "perl $bin_FindCDS_mRNA $hseq->{$sid}{seq} $prefix > $tmp_qprot";
		#			print "! cmd=[$cmd]\n";
					system($cmd);
					
					if(-e $tmp_qprot){
						my $AoR = Revseq2array($tmp_qprot);
						
						my $num_seq = 0;
						foreach my $rtmp (@{$AoR}){
							my $cand_fasta = $rtmp->[0];
							my $split_status = $rtmp->[1];
							SAVE($tmp_qprot, $cand_fasta);
							
							if(! -e $tmp_qprot){
								next;
							}
							$num_seq++;
							
							my $cmd1 = "$bin_matcher -asequence $tmp_dbprot -bsequence $tmp_qprot -outfile $tmp_align";
							print "! [$method3] [$num_seq] | cmd=[$cmd1]\n";
							system($cmd1);
							
							if(-e $tmp_align){
								my $hmatcher = Read_matcher($tmp_align);
								if($hmatcher->{score} && $hmatcher->{score} ne '-'){
									my $hitseq = PickupFasta($tmp_qprot, $gid, $qpref, $method3);
									
									my @tmp;
									push(@tmp, $hitseq->[0]);		# revID
									push(@tmp, $hitseq->[1]);		# seq
									push(@tmp, $method3);
									push(@tmp, $hmatcher->{score} * $hitseq->[2]);
									push(@tmp, $hmatcher->{native});
									push(@{$AoHit}, \@tmp);
									
									print "! [$method3] [$num_seq] | [$hitseq->[0]] score=[$hmatcher->{score}]\n";
								}
							}
							
							Rm($tmp_qprot);
							Rm($tmp_align);
						}
					}
				}
				
				if(@{$AoHit}){
					@{$AoHit} = sort {$b->[3] <=> $a->[3]} @{$AoHit};
					my $tophit = $AoHit->[0];
					
					print "! tophit = [$tophit->[0]] [$tophit->[2]]\n";
					my $tophit_fasta = ">".$tophit->[0]."\n".$tophit->[1]."\n";
					SAVE($tmp_qprot, $tophit_fasta);
					SAVE($tmp_align, $tophit->[4]);
					
					BindFasta($combined_qprot, $qpref, $gid, $tophit->[0], $tophit->[1], $tophit->[2]);
				}
			}
		}
	}
	
	my $keep = 0;
	if($keep eq '0'){
		Rm($tmp_dbprot);
		Rm("$tmp_dbprot.md5");
		Rm("$tmp_dbprot.suf");
		Rm($tmp_qfasta);
		Rm("$tmp_qfasta.fai");
		Rm($tmp_qgff);
	}
	if(-e "$tmp_dbprot.protein.bck"){
		system("rm $tmp_dbprot.protein.*");
	}
	if(-e "$tmp_qfasta.dna.tis"){
		system("rm $tmp_qfasta.dna.*");
	}
}

END:{
	print "\n";
}


################################################################################
#-------------------------------------------------------------------------------
sub Search_previous{
my ($combined_qprot, $gid, $qpref) = @_;

open(my $fh, "<", $combined_qprot);
my $judge = "false";
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /\>/){
		my $ID = $line;
		$ID =~ s/\>//;
		my @tmp = split(/_/, $ID);
		if($tmp[0] && $tmp[1] && $tmp[0] eq $qpref && $tmp[1] eq $gid){
			$judge = "true";
			last;
		}
	}
}
close $fh;

return $judge;
}


#-------------------------------------------------------------------------------
sub Revseq2array{
my $tmp_qprot = shift;

open(my $fh, "<", $tmp_qprot);
my $ID;
my $seq;
my $AoR = [];
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /\>/){
		if($seq){
			my $native_seq = $seq;
			$seq =~ s/\./*/g;
			my @tmp = split(/\*/, $seq);
			$seq = $tmp[0];
			
			my $status = 0;
			if($tmp[1]){
				$status = 1;
			}
			
			if($seq){
				my @CHAR = split(//, $seq);
				if($CHAR[0] ne 'M'){
					$status = 1;
					my @SPM = split(/M/, $seq);
					
					if($SPM[0]){
						shift(@SPM);
						my $num_SPM = @SPM;
						for(my $i = 0; $i < $num_SPM; $i++){
							my @eachSPM;
							for(my $j = $i; $j < $num_SPM; $j++){
								push(@eachSPM, $SPM[$j]);
							}
							my $each_tmpseq = "M".join("M", @eachSPM);
							
							my $r = ">".$ID."\n".$each_tmpseq."\n";
							my @rtmp;
							push(@rtmp, $r);
							push(@rtmp, $status);
							push(@{$AoR}, \@rtmp);
						}
					}
				}
				else{
					my $r = ">".$ID."\n".$seq."\n";
					my @rtmp;
					push(@rtmp, $r);
					push(@rtmp, $status);
					push(@{$AoR}, \@rtmp);
				}
			}
		}
		
		$ID = $line;
		$ID =~ s/\>//;
		my @tmp = split(/\s/, $ID);
		$ID = $tmp[0];
		$seq = "";
	}
	else{
		$seq .= $line;
	}
}
close $fh;

if($seq){
	my $native_seq = $seq;
	$seq =~ s/\./*/g;
	my @tmp = split(/\*/, $seq);
	$seq = $tmp[0];
	
	my $status = 0;
	if($tmp[1]){
		$status = 1;
	}
	
	if($seq){
		my @CHAR = split(//, $seq);
		if($CHAR[0] ne 'M'){
			$status = 1;
			my @SPM = split(/M/, $seq);
			
			if($SPM[0]){
				shift(@SPM);
				my $num_SPM = @SPM;
				for(my $i = 0; $i < $num_SPM; $i++){
					my @eachSPM;
					for(my $j = $i; $j < $num_SPM; $j++){
						push(@eachSPM, $SPM[$j]);
					}
					my $each_tmpseq = "M".join("M", @eachSPM);
					
					my $r = ">".$ID."\n".$each_tmpseq."\n";
					my @rtmp;
					push(@rtmp, $r);
					push(@rtmp, $status);
					push(@{$AoR}, \@rtmp);
				}
			}
		}
		else{
			my $r = ">".$ID."\n".$seq."\n";
			my @rtmp;
			push(@rtmp, $r);
			push(@rtmp, $status);
			push(@{$AoR}, \@rtmp);
		}
	}
}

return $AoR;
}


#-------------------------------------------------------------------------------
sub CheckFasta{
my $tmp_qprot = shift;

open(my $fh, "<", $tmp_qprot);
my $ID;
my $seq;
my $num_seq = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /\>/){
		if($seq){
			$num_seq++;
		}
		
		$ID = $line;
		$ID =~ s/\>//;
		my @tmp = split(/\s/, $ID);
		$ID = $tmp[0];
		$seq = "";
	}
	else{
		$seq .= $line;
	}
}
close $fh;

if($seq){
	$num_seq++;
}

return $num_seq;
}


#-------------------------------------------------------------------------------
sub PickupFasta{
my ($tmp_qprot, $gid, $qpref, $method) = @_;

open(my $fh, "<", $tmp_qprot);
my $ID;
my $seq;
my $AoR = [];
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /\>/){
		if($seq){
			my @CHAR = split(//, $seq);
			my $aa_1st = $CHAR[0];
			my $len = length($seq);
			my $revID = $qpref."_".$gid."_".$ID;
#			my $revID = $ID;
			
			my $score = $len;
			if($aa_1st && $aa_1st =~ /M/i){
				$score *= 1000000000;
			}
			
			my @tmp;
			push(@tmp, $revID);
			push(@tmp, $seq);
			push(@tmp, $score);
			push(@tmp, $len);
			push(@tmp, $aa_1st);
			push(@tmp, $ID);
			push(@{$AoR}, \@tmp);
		}
		
		$ID = $line;
		$ID =~ s/\>//;
		my @tmp = split(/\s/, $ID);
		$ID = $tmp[0];
		$seq = "";
	}
	else{
		$seq .= $line;
	}
}
close $fh;

if($seq){
	my @CHAR = split(//, $seq);
	my $aa_1st = $CHAR[0];
	my $len = length($seq);
	my $revID = $qpref."_".$gid."_".$ID;
#	my $revID = $ID;
	
	my $score = $len;
	if($aa_1st && $aa_1st =~ /M/i){
		$score *= 1000000000;
	}
	
	my @tmp;
	push(@tmp, $revID);
	push(@tmp, $seq);
	push(@tmp, $score);
	push(@tmp, $len);
	push(@tmp, $aa_1st);
	push(@tmp, $ID);
	push(@{$AoR}, \@tmp);
}

@{$AoR} = sort {$b->[2] <=> $a->[2]} @{$AoR};

return $AoR->[0];		# return only tophit
}


#-------------------------------------------------------------------------------
sub BindFasta{
my ($combined_qprot, $qpref, $gid, $ID, $seq, $method) = @_;

if($seq =~ /\./){
	my @tmp = split(/\./, $seq);
	if($tmp[0]){
		$seq = $tmp[0];
	}
	else{
		$seq = "";
	}
}
if($seq){
	my $r = ">".$ID."\n".$seq."\n";
	
	open(my $rfh, ">>", $combined_qprot);
	print $rfh $r;
	close $rfh;
}

}


#-------------------------------------------------------------------------------
sub BindGff{
my ($combined_qgff, $lines) = @_;

open(my $rfh, ">>", $combined_qgff);
print $rfh $lines;
close $rfh;

}


#-------------------------------------------------------------------------------
sub ReadGff{
my $file = shift;
my $gid = shift;		# e.g. GlymaJER.01G057564.1
my $qpref = shift;		# e.g. Gmax_508_v4.0
my $revID = shift;		# e.g. Gmax_508_v4.0_GlymaJER.01G057564.1_mRNA1_gth
my $ID = shift;
my $method = shift;		# e.g. gth
my $qpos0 = shift;

my $tmp_prefix1 = $qpref."_".$gid."_";
my $tmp_prefix2 = "_".$method;
my $tidx = $revID;
$tidx =~ s/$tmp_prefix1//;			# e.g. mRNA1_gth
$tidx =~ s/$tmp_prefix2//;			# e.g. mRNA1

my $rev_gid = $qpref."_".$gid;
my $rev_tid = $rev_gid.".t1";

open(my $fh, "<", $file);
my $hash = {};
my $CDS = [];
my $Exon = [];
my $Other = [];
my $strand = "";
my $native = "# ".$revID." | ".$ID."\n";
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
#	Gm01	gth	gene	3	3433	.	+	.	ID=gene1
#	Gm01	gth	mRNA	3	3433	.	+	.	ID=mRNA1;Parent=gene1;Target=GlymaJER.01G057564.1.t1 2 1068 +
#	Gm01	gth	exon	3	121	1	+	.	Parent=mRNA1
#	Gm01	gth	five_prime_cis_splice_site	122	123	0.05	+	.	Parent=mRNA1
#	Gm01	gth	three_prime_cis_splice_site	250	251	0.05	+	.	Parent=mRNA1
#	Gm01	gth	exon	252	2967	0.998	+	.	Parent=mRNA1
#	Gm01	gth	CDS	1565	2967	.	+	0	ID=CDS1;Parent=mRNA1
#	Gm01	gth	five_prime_cis_splice_site	2968	2969	0.05	+	.	Parent=mRNA1
#	Gm01	gth	three_prime_cis_splice_site	3068	3069	0.05	+	.	Parent=mRNA1
#	Gm01	gth	exon	3070	3433	1	+	.	Parent=mRNA1
#	Gm01	gth	CDS	3070	3433	.	+	1	ID=CDS1;Parent=mRNA1
	
	if($line){
		$native .= $line."\n";
		
		my @CHAR = split(//, $line);
		my @A = split(/\t/, $line);
		if($A[8] && $CHAR[0] !~ /\#/){
			my @tag = split(/\;/, $A[8]);
			
			$A[3] += $qpos0 - 1;
			$A[4] += $qpos0 - 1;
			
			if($A[2] eq 'gene'){
				$A[8] = "ID=".$rev_gid;
				$hash->{gene} = join("\t", @A)."\n";
				$strand = $A[6];
			}
			elsif($A[2] eq 'mRNA' || $A[2] eq 'transcript'){
				my $tid;
				foreach my $val (@tag){
					if($val =~ /ID\=/){
						$tid = $val;
						$tid =~ s/ID\=//;
						last;
					}
				}
				if(! $tid || $tid ne $tidx){
					next;
				}
				
				$A[2] = "transcript";
				$A[8] = "ID=".$rev_tid."\;Parent=".$rev_gid;
				if($tag[2] =~ /Target\=/){
					$A[8] .= $tag[2];
				}
				$hash->{transcript} = join("\t", @A)."\n";
				
				$A[2] = "gene";
				$A[8] = "ID=".$rev_gid;
				$hash->{gene_sub} = join("\t", @A)."\n";
			}
			else{
				my $tid;
				foreach my $val (@tag){
					if($val =~ /Parent\=/){
						$tid = $val;
						$tid =~ s/Parent\=//;
						last;
					}
				}
				if(! $tid || $tid ne $tidx){
					next;
				}
				
				if($A[2] eq 'CDS'){
					push(@{$CDS}, \@A);
				}
				elsif($A[2] eq 'exon'){
					push(@{$Exon}, \@A);
				}
				else{
					push(@{$Other}, \@A);
				}
			}
		}
	}
}
close $fh;

if($strand eq '+'){
	@{$CDS} = sort {$a->[3] <=> $b->[3]} @{$CDS};
	@{$Exon} = sort {$a->[3] <=> $b->[3]} @{$Exon};
	@{$Other} = sort {$a->[3] <=> $b->[3]} @{$Other};
}
else{
	@{$CDS} = sort {$b->[3] <=> $a->[3]} @{$CDS};
	@{$Exon} = sort {$b->[3] <=> $a->[3]} @{$Exon};
	@{$Other} = sort {$b->[3] <=> $a->[3]} @{$Other};
}

my $num_CDS = @{$CDS};
my $num_exon = @{$Exon};
my $num_other = @{$Other};

my $AoA = [];
for(my $i = 0; $i < $num_CDS; $i++){
	my $j = $i + 1;
	my $A = $CDS->[$i];
	$A->[8] = "ID=".$rev_tid.".CDS".$j."\;Parent=".$rev_tid;
	push(@{$AoA}, $A);
}
for(my $i = 0; $i < $num_exon; $i++){
	my $j = $i + 1;
	my $A = $Exon->[$i];
	$A->[8] = "ID=".$rev_tid.".exon".$j."\;Parent=".$rev_tid;
	push(@{$AoA}, $A);
}
for(my $i = 0; $i < $num_other; $i++){
	my $A = $Other->[$i];
	$A->[8] = "Parent=".$rev_tid;
	push(@{$AoA}, $A);
}

@{$AoA} = sort {$a->[3] <=> $b->[3]} @{$AoA};

if(! $hash->{gene}){
	$hash->{gene} = $hash->{gene_sub};
}

my $rgff = $hash->{gene}.$hash->{transcript};
my $rgff_transcript = $hash->{transcript};
foreach my $A (@{$AoA}){
	$rgff .= join("\t", @{$A})."\n";
	$rgff_transcript .= join("\t", @{$A})."\n";
}

$native .= "-------selected-------\n";
$native .= $rgff;
$native .= "-------selected-------\n\n";

return ($rgff, $native, $rgff_transcript);
}


#-------------------------------------------------------------------------------
sub Read_matcher{
my $file = shift;

open(my $fh, "<", $file);
my $hash = {};
my $native;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line){
		$native .= $line."\n";
	}
	else{
		$native .= "\n";
	}
	
	# Aligned_sequences: 2
	# 1: MELO.jh000002.1.t1
	# 2: mRNA1
	# Matrix: EBLOSUM62
	# Gap_penalty: 14
	# Extend_penalty: 4
	#
	# Length: 282
	# Identity:     282/282 (100.0%)
	# Similarity:   282/282 (100.0%)
	# Gaps:           0/282 ( 0.0%)
	# Score: 1491
	
	my @A;
	my $sw = 0;
	if($line =~ /Length:/){
		@A = split(/Length:/, $line);
		$A[1] =~ s/\s//g;
		
		if($A[1]){
			$hash->{len} = $A[1];
		}
		else{
			$hash->{len} = "-";
		}
	}
	elsif($line =~ /Identity:/){
		@A = split(/Identity:/, $line);
		$A[1] =~ s/\s//g;
		$A[1] =~ s/\)//g;
		$A[1] =~ s/\%//g;
		
		my @B = split(/\(/, $A[1]);
		if($B[0]){
			$hash->{identity_str} = $B[0];
		}
		else{
			$hash->{identity_str} = "-";
		}
		if($B[1]){
			$hash->{identity} = $B[1];
		}
		else{
			$hash->{identity} = 0;
		}
	}
	elsif($line =~ /Similarity:/){
		@A = split(/Similarity:/, $line);
		$A[1] =~ s/\s//g;
		$A[1] =~ s/\)//g;
		$A[1] =~ s/\%//g;
		
		my @B = split(/\(/, $A[1]);
		if($B[0]){
			$hash->{similarity_str} = $B[0];
		}
		else{
			$hash->{similarity_str} = "-";
		}
		if($B[1]){
			$hash->{similarity} = $B[1];
		}
		else{
			$hash->{similarity} = "-";
		}
	}
	elsif($line =~ /Gaps:/){
		@A = split(/Gaps:/, $line);
		$A[1] =~ s/\s//g;
		$A[1] =~ s/\)//g;
		$A[1] =~ s/\%//g;
		
		my @B = split(/\(/, $A[1]);
		if($B[0]){
			$hash->{gap_str} = $B[0];
		}
		else{
			$hash->{gap_str} = "-";
		}
		if($B[1]){
			$hash->{gap} = $B[1];
		}
		else{
			$hash->{gap} = "-";
		}
	}
	elsif($line =~ /Score:/){
		@A = split(/Score:/, $line);
		$A[1] =~ s/\s//g;
		
		if($A[1]){
			$hash->{score} = $A[1];
		}
		else{
			$hash->{score} = "-";
		}
	}
}
close $fh;

$hash->{native} = $native;

return $hash;
}


#-------------------------------------------------------------------------------
sub Getseqinfo{
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
my $Mstart = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ />/){
		if($seq){
			$hash->{$ID}{seq} = $seq;
			$hash->{$ID}{len} = length($seq);
			
			my $aa1 = substr($seq, 0, 1);
			if($aa1 =~ /m/i){
				$Mstart->{$ID} = 1;
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
	$hash->{$ID}{seq} = $seq;
	$hash->{$ID}{len} = length($seq);
	
	my $aa1 = substr($seq, 0, 1);
	if($aa1 =~ /m/i){
		$Mstart->{$ID} = 1;
	}
}

my $rhash = {};
$rhash->{cnt} = $cnt;
$rhash->{data} = $hash;
$rhash->{Mstart} = $Mstart;

return $rhash;
}


#-------------------------------------------------------------------------------
sub Add_feature_info{
my $hint = shift;
my $strand = shift;
my $rhint = shift;

open(my $fh, "<", $hint) or return "null";
my $hash = {};
my $ep = {};
my $exon = {};
my $intron = {};
my $nr = {};
my $ch = {};
my $stc = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	my @A = split(/\t/, $line);
	my @tag = split(/\;/, $A[8]);
	
	if($tag[0] =~ /mult/){
		next;
	}
	$tag[0] =~ s/grp\=//;
	
	if($strand->{$tag[0]}){
		$A[6] = $strand->{$tag[0]};
	}
	else{
		print "! error : missing strand information for [$tag[0]] | [$line]\n";
		die;
	}
	
	$hash->{$tag[0]} .= join("\t", @A)."\n";
}
close $fh;

my @ID = keys(%{$hash});
@ID = sort {$a cmp $b} @ID;

foreach my $id (@ID){
	my @L = split(/\n/, $hash->{$id});
	my $AoA = [];
	my $cntex = 0;
	foreach my $line (@L){
		my @A = split(/\t/, $line);
		
		if($A[2] eq 'ep' || $A[2] eq 'exon'){
			$cntex++;
		}
		
		push(@{$AoA}, \@A);
	}
	
	@{$AoA} = sort {$a->[3] <=> $b->[3]} @{$AoA};
	my $numA = @{$AoA};
	
	my $numex = 0;
	for(my $i = 0; $i < $numA; $i++){
		my $A = $AoA->[$i];
		
		if($A->[2] eq 'ep'){
			$A->[2] = "exon";
			$ep->{$id} .= join("\t", @{$A})."\n";
			$nr->{$id} .= $A->[3]."\n";
			$nr->{$id} .= $A->[4]."\n";
			$ch->{$id} = $A->[0];
		}
		elsif($A->[2] eq 'exon'){
			$exon->{$id} .= join("\t", @{$A})."\n";
			$nr->{$id} .= $A->[3]."\n";
			$nr->{$id} .= $A->[4]."\n";
			$ch->{$id} = $A->[0];
		}
		elsif($A->[2] eq 'intron'){
			$intron->{$id} .= join("\t", @{$A})."\n";
			$nr->{$id} .= $A->[3]."\n";
			$nr->{$id} .= $A->[4]."\n";
			$ch->{$id} = $A->[0];
		}
	}
}

my $hmin = {};
my $hmax = {};
foreach my $id (@ID){
	my @Pos = split(/\n/, $nr->{$id});
	my $min;
	my $max;
	foreach my $p (@Pos){
		unless($min && $max){
			$min = $p;
			$max = $p;
		}
		else{
			if($min > $p){
				$min = $p;
			}
			if($max < $p){
				$max = $p;
			}
		}
	}
	
	if($min){
		$hmin->{$id} = $min;
	}
	if($max){
		$hmax->{$id} = $max;
	}
}

my $r;
#my $tss = "transcription_start_site";		# transcription start site		-> cause error
#my $tts = "transcription_end_site";		# transcription terminal site	-> cause error
my $tss = "tss";		# transcription start site 
my $tts = "tts";		# transcription terminal site 
foreach my $id (@ID){
	my $eachr;
	if($strand->{$id} eq '+'){
		if($hmin->{$id}){
			$eachr .= $ch->{$id}."\tb2h\t".$tss."\t".$hmin->{$id}."\t".$hmin->{$id}."\t0\t".$strand->{$id}."\t\.\tgrp\=".$id."\;pri\=4\;src\=E\n";
		}
	}
	elsif($strand->{$id} eq '-'){
		if($hmin->{$id}){
			$eachr .= $ch->{$id}."\tb2h\t".$tts."\t".$hmin->{$id}."\t".$hmin->{$id}."\t0\t".$strand->{$id}."\t\.\tgrp\=".$id."\;pri\=4\;src\=E\n";
		}
	}
	
	if($strand->{$id} eq '+' && $stc->{$id}){
		$eachr .= $stc->{$id};
	}
	
	if($ep->{$id}){
		$eachr .= $ep->{$id};
	}
	if($exon->{$id}){
		$eachr .= $exon->{$id};
	}
	if($intron->{$id}){
		$eachr .= $intron->{$id};
	}
	
	if($strand->{$id} eq '-' && $stc->{$id}){
		$eachr .= $stc->{$id};
	}
	
	if($strand->{$id} eq '+'){
		if($hmax->{$id}){
			$eachr .= $ch->{$id}."\tb2h\t".$tts."\t".$hmax->{$id}."\t".$hmax->{$id}."\t0\t".$strand->{$id}."\t\.\tgrp\=".$id."\;pri\=4\;src\=E\n";
		}
	}
	elsif($strand->{$id} eq '-'){
		if($hmax->{$id}){
			$eachr .= $ch->{$id}."\tb2h\t".$tss."\t".$hmax->{$id}."\t".$hmax->{$id}."\t0\t".$strand->{$id}."\t\.\tgrp\=".$id."\;pri\=4\;src\=E\n";
		}
	}
	
	my $revr = $eachr;
	$revr =~ s/src\=E/src\=M/g;
	
#	$r .= $eachr.$revr;
	$r .= $eachr;
}

if($r){
	open(my $rfh, ">", $rhint) or die;
	print $rfh $r;
	close $rfh;
}

return "OK";
}


#-------------------------------------------------------------------------------
sub Select_tophit{
my $psl = shift;
my $rpsl = shift;

print "\n! reading [$psl]...\n";
open(my $fh, "<", $psl) or return 'null';
my $hash = {};
my $data = {};
my $strand = {};
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	my @A = split(/\t/, $line);
	my $match = $A[0];
	my $mismatch = $A[1];
	my $score = $match - $mismatch;
	
	unless($hash->{$A[9]}){
		$data->{$A[9]} = $line;
		$strand->{$A[9]} = $A[8];
		$hash->{$A[9]} = $score;
	}
	elsif($score > $hash->{$A[9]}){
		$data->{$A[9]} = $line;
		$strand->{$A[9]} = $A[8];
		$hash->{$A[9]} = $score;
	}
	$cnt++;
}
close $fh;

my @ID = keys(%{$data});
@ID = sort {$a cmp $b} @ID;
my $numID = @ID;

my $r;
foreach my $id (@ID){
	$r .= $data->{$id}."\n";
}

if($r){
	open(my $rfh, ">", $rpsl) or die;
	print $rfh $r;
	close $rfh;
	
	return $strand;
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
else{
	print "! missing string to save [$file]...\n";
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
}

}


#-------------------------------------------------------------------------------
sub Rm_ifempty{
my $file = shift;

if(-e $file){
	open(my $fh, "<", $file) or die;
	my $dat = "";
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line){
			$dat .= $line."\n";
		}
	}
	close $fh;
	
	if(! $dat){
		system("rm $file");		# empty file
	}
}

}


#-------------------------------------------------------------------------------
sub Rm{
my $file = shift;

if(-e $file){
	system("rm $file");
}

}

