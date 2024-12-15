#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use FindBin;

my $each_rvcf = shift;
my $candidate = shift;
my $vcf = shift;
my $header = shift;
my $dbfasta = shift;
my $sid = shift;
my $sample_uniqueID = shift;
my $asm2sv_tsv = shift;
my $resume_from_native = shift;
my $neighbor_bp = shift;
my $chrinfo = shift;
my $dbprefix = shift;
my $qprefix = shift;
my $t = shift;

if(! defined $each_rvcf || ! defined $t){
	goto END;
}

my ($hdbseq, $SID) = Open_fasta_as_hash($dbfasta);
my $num_SID = @{$SID};
my $htmp_chrinfo = Read_seqinfo($chrinfo);
my $hseqinfo = $htmp_chrinfo->{seqID};

print " thread [$t] : [$sid]...\n";
my $hAoPL = {};
my $AoGID = [];
my $hoverlap = {};
if(-e $asm2sv_tsv){
	print " thread [$t] : [$sid] | reading [$asm2sv_tsv] ...\n";
	open(my $fh, "<", $asm2sv_tsv) or die;
	my $cnt_seqnormal = 0;
	my $cnt_unmatchseq = 0;
	my $htmp = {};
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
			my $rgpos0 = $A[5];
			my $rgpos1 = $A[6];
			my $qsid = $A[8];
			my $qpos0 = $A[9];
			my $qpos1 = $A[10];
			my $normality_score = $A[21];
			my $asm2sv_judge = $A[22];
			
			if($rpos0 =~ /\D/ || $rpos1 =~ /\D/ || $qpos0 =~ /\D/ || $qpos1 =~ /\D/ || $rsid ne $sid){
				next;
			}
			
			if(! $hseqinfo->{$dbprefix}{$rsid}){
				print "! error: missing seqname alias for DB [$dbprefix] [$rsid]...\n";
				print "$asm2sv_tsv\n$line\n";
				die;
			}
			if(! $hseqinfo->{$qprefix}{$qsid}){
				print "! error: missing seqname alias for q [$qprefix] [$qsid]...\n";
				print "$asm2sv_tsv\n$line\n";
				die;
			}
			my $rsid_alias = $hseqinfo->{$dbprefix}{$rsid};
			my $qsid_alias = $hseqinfo->{$qprefix}{$qsid};
			
			if($rsid_alias ne $qsid_alias){
				print " - skip [$gid] as unmatched Chr alias | DB [$rsid_alias] != q [$qsid_alias]\n";
				$cnt_unmatchseq++;
				next;
			}
			if($normality_score && $normality_score =~ /\d/ && $normality_score > 0.99 && $normality_score < 1.01 && $asm2sv_judge && $asm2sv_judge eq 'present'){
				$htmp->{$A[2]}{normal} .= $line."\n";
				$cnt_seqnormal++;
			}
			else{
				$htmp->{$A[2]}{cand} .= $line."\n";
			}
			
			for(my $p = $rgpos0; $p <= $rgpos1; $p++){
				$hoverlap->{$p} .= $gid."\n";
			}
			
			my @tmp;
			push(@tmp, $gid);
			push(@tmp, $rpos0);
			push(@tmp, $rpos1);
			push(@{$AoGID}, \@tmp);
		}
	}
	close $fh;
	
	if($htmp->{$sid}{normal}){
		my @L = split(/\n/, $htmp->{$sid}{normal});
		my $AoPL = [];
		my $cnt_tmp = 0;
		foreach my $line (@L){
			my @A = split(/\t/, $line);
			push(@{$AoPL}, \@A);
			$cnt_tmp++;
		}
		@{$AoPL} = sort {$a->[3] <=> $b->[3]} @{$AoPL};
		
		$hAoPL->{$sid}{normal} = $AoPL;
	}
	
	if($htmp->{$sid}{cand}){
		my @L = split(/\n/, $htmp->{$sid}{cand});
		my $AoPL = [];
		my $cnt_tmp = 0;
		foreach my $line (@L){
			my @A = split(/\t/, $line);
			push(@{$AoPL}, \@A);
			$cnt_tmp++;
		}
		@{$AoPL} = sort {$a->[3] <=> $b->[3]} @{$AoPL};
		
		$hAoPL->{$sid}{cand} = $AoPL;
	}
}

@{$AoGID} = sort {$a->[1] <=> $b->[1]} @{$AoGID};

my @Gpos = keys(%{$hoverlap});
@Gpos = sort {$a <=> $b} @Gpos;
my $prev_gpos;
my $n_gregion = 0;
my $hngrp2gene = {};
my $hgene2ngrp = {};
foreach my $rpos (@Gpos){
	if(! $prev_gpos){
		$n_gregion++;
		if($hoverlap->{$rpos}){
			my @tmpGID = split(/\n/, $hoverlap->{$rpos});
			foreach my $gid (@tmpGID){
				unless($hngrp2gene->{$n_gregion}{pos0}){
					$hngrp2gene->{$n_gregion}{pos0} = $rpos;
				}
				$hngrp2gene->{$n_gregion}{pos1} = $rpos;
				$hngrp2gene->{$n_gregion}{gene}{$gid} = 1;
				$hgene2ngrp->{$gid} = $n_gregion;
#				print " - $rpos = [$gid] ($n_gregion)\n";
			}
		}
		$prev_gpos = $rpos;
	}
	else{
		if($rpos != $prev_gpos + 1){
			$n_gregion++;
		}
		if($hoverlap->{$rpos}){
			my @tmpGID = split(/\n/, $hoverlap->{$rpos});
			foreach my $gid (@tmpGID){
				unless($hngrp2gene->{$n_gregion}{pos0}){
					$hngrp2gene->{$n_gregion}{pos0} = $rpos;
				}
				$hngrp2gene->{$n_gregion}{pos1} = $rpos;
				$hngrp2gene->{$n_gregion}{gene}{$gid} = 1;
				$hgene2ngrp->{$gid} = $n_gregion;
#				print " - $rpos = [$gid] ($n_gregion)\n";
			}
		}
		$prev_gpos = $rpos;
	}
}

print " thread [$t] : [$sid] [$n_gregion] nr gene region considering overlap\n";

my $revheader1 =<<"EOS";
##fileformat=VCFv4.2
##source=Sniffles2_2.2
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Breakend; Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=ID,Number=1,Type=String,Description="Individual sample SV ID for multi-sample output">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=GT,Description="Genotype filter">
##FILTER=<ID=SUPPORT_MIN,Description="Minimum read support filter">
##FILTER=<ID=STDEV_POS,Description="SV Breakpoint standard deviation filter">
##FILTER=<ID=STDEV_LEN,Description="SV length standard deviation filter">
##FILTER=<ID=COV_MIN,Description="Minimum coverage filter">
##FILTER=<ID=COV_CHANGE,Description="Coverage change filter">
##FILTER=<ID=COV_CHANGE_FRAC,Description="Coverage fractional change filter">
##FILTER=<ID=MOSAIC_AF,Description="Mosaic maximum allele frequency filter">
##FILTER=<ID=ALN_NM,Description="Length adjusted mismatch filter">
##FILTER=<ID=STRAND,Description="Strand support filter">
##FILTER=<ID=SVLEN_MIN,Description="SV length filter">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">
##INFO=<ID=MOSAIC,Number=0,Type=Flag,Description="Structural variation classified as putative mosaic">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">
##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">
##INFO=<ID=SUPPORT_INLINE,Number=1,Type=Integer,Description="Number of reads supporting an INS/DEL SV (non-split events only)">
##INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description="Number of soft-clipped reads putatively supporting the long insertion SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">
##INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Standard deviation of structural variation start position">
##INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Standard deviation of structural variation length">
##INFO=<ID=COVERAGE,Number=.,Type=Float,Description="Coverages near upstream, start, center, end, downstream of structural variation">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Strands of supporting reads for structural variant">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count, summed up over all samples">
##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="List of read support for all samples">
##INFO=<ID=CONSENSUS_SUPPORT,Number=1,Type=Integer,Description="Number of reads that support the generated insertion (INS) consensus sequence">
##INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads (if enabled with --output-rnames)">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=NM,Number=.,Type=Float,Description="Mean number of query alignment length adjusted mismatches of supporting reads">
##INFO=<ID=PHASE,Number=.,Type=String,Description="Phasing information derived from supporting reads, represented as list of: HAPLOTYPE,PHASESET,HAPLOTYPE_SUPPORT,PHASESET_SUPPORT,HAPLOTYPE_FILTER,PHASESET_FILTER">
EOS

my $revheader2 = "";
my $revheader3 = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".$sample_uniqueID."\n";
my $rvcf = "";
if(! $resume_from_native && -e $header){
	open(my $fh, "<", $header) or die;
	my $nr = {};
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line){
			if(! $nr->{$line} && $line !~ /command\=/ && $line !~ /fileDate/ && $line !~ /\#CHROM/){
				if($line !~ /contig\=/i){
					$revheader1 .= $line."\n";
				}
				else{
					$revheader2 .= $line."\n";
				}
				$nr->{$line} = 1;
			}
		}
	}
	close $fh;
}

foreach my $sid (@{$SID}){
	$revheader2 .="##contig=<ID=".$sid.",length=".$hdbseq->{$sid}{len}.">\n";
}

my $list_remove = "";
my $salt = "ACGTacgt";
if(-e $candidate){
	print " thread [$t] : [$sid] | reading [$candidate] ...\n";
	open(my $fh, "<", $candidate) or die;
	my $hvardata = {};
	my $nrl = {};
	my $cnt_seqnormal = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line && $line !~ /\#/){
			my @A = split(/\t/, $line);
			my $numA = @A;
			if($numA < 9 || $A[1] =~ /\D/ || ! $A[9]){
				next;
			}
			elsif($A[0] && $A[0] ne $sid){
				next;
			}
			
			my @tag = split(/\:/, $A[9]);		# 1/1:60:0:226
			if(! $tag[3] || $tag[3] =~ /\D/){
				next;
			}
			if($nrl->{$line}){
				next;
			}
			$nrl->{$line} = 1;
			
			my $lenR = length($A[3]);
			my $lenA = length($A[4]);
			my $lendf = abs($lenR - $lenA);
			
#			if($A[1] eq '31483465'){
#				print "$line\n";
#			}
			
			my $judge_within_normal_region = "false";
			my $judge_within_indel_region = "false";
			my $info_norm = "-\t-\t-";
			my $info_indel = "-\t-\t-\t-\t-";
			if($hAoPL->{$A[0]}{normal}){
				my $AoPL = $hAoPL->{$A[0]}{normal};
				foreach my $PL (@{$AoPL}){
					for(my $p = $A[1]; $p <= $A[1] + $lenR; $p++){
						if($PL->[5] <= $p && $p <= $PL->[6]){
							$judge_within_normal_region = "true";
							$info_norm = $PL->[5]."\t".$PL->[6]."\t".$PL->[23];
							last;
						}
					}
					if($judge_within_normal_region eq 'true'){
						last;
					}
				}
			}
			
			my $ngrp;
			if($hAoPL->{$A[0]}{cand}){
				my $AoPL = $hAoPL->{$A[0]}{cand};
				foreach my $PL (@{$AoPL}){
					for(my $p = $A[1]; $p <= $A[1] + $lenR; $p++){
						if($PL->[5] - $neighbor_bp <= $p && $p <= $PL->[6] + $neighbor_bp){
							$judge_within_indel_region = "true";
							
							if(! defined $hgene2ngrp->{$PL->[23]}){
								print "! error: unable to find ngrp for [$PL->[23]]...\n";
								die;
							}
							$ngrp = $hgene2ngrp->{$PL->[23]};
							
							my $htmpg = $hngrp2gene->{$ngrp}{gene};
							my @memgene = keys(%{$htmpg});
							@memgene = sort {$a cmp $b} @memgene;
							
							$info_indel = $PL->[5]."\t".$PL->[6]."\t".join(",", @memgene)."\t".$ngrp."\ttrue";
							last;
						}
					}
					if($judge_within_indel_region eq 'true'){
						last;
					}
				}
			}
			
			if($judge_within_normal_region eq 'true' && $judge_within_indel_region eq 'false'){
				$list_remove .= $A[0]."\t".$A[1]."\t".$A[2]."\t".$info_norm."\t".$judge_within_normal_region."\t".$info_indel."\trm\n";
				$cnt_seqnormal++;
				next;
			}
			elsif($judge_within_indel_region eq 'false'){
				$list_remove .= $A[0]."\t".$A[1]."\t".$A[2]."\t".$info_norm."\t".$judge_within_normal_region."\t".$info_indel."\trm\n";
				$cnt_seqnormal++;
				next;
			}
			else{
				$list_remove .= $A[0]."\t".$A[1]."\t".$A[2]."\t".$info_norm."\t".$judge_within_normal_region."\t".$info_indel."\tpass\n";
				$hngrp2gene->{$ngrp}{vcf} .= $line."\n";
			}
			
			my $tmpREF = $A[3];
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
			$hvardata->{$A[0]}{$varid}{vcf} = $line;
			$hvardata->{$A[0]}{$varid}{pos0} = $A[1];
			$hvardata->{$A[0]}{$varid}{lenR} = $lenR;
			$hvardata->{$A[0]}{$varid}{lenA} = $lenA;
			$hvardata->{$A[0]}{$varid}{lendf} = $lendf;
		}
	}
	close $fh;
	
	my $cnt_selected = 0;
	my $cnt_overlap = 0;
	for(my $ngrp = 1; $ngrp <= $n_gregion; $ngrp++){
		if($hngrp2gene->{$ngrp}{vcf}){
			my @L = split(/\n/, $hngrp2gene->{$ngrp}{vcf});
			my $AoAL = [];
			my $hsubov = {};
			my $hsvsize = {};
			foreach my $line (@L){
				my @A = split(/\t/, $line);
				my $lenR = length($A[3]);
				my $lenA = length($A[4]);
				my $eachpend = $A[1] + $lenR;
				my $lendfx = $lenR;
				if($lenA > $lenR){
					$lendfx = $lenA;
				}
#				print " thread [$t] : [$sid] | $A[1] - $eachpend ($lendfx)\n";
				
				$hsvsize->{$A[1]}{data} = $line;
				$hsvsize->{$A[1]}{lendf} = $lendfx;
				
				for(my $p = $A[1]; $p < $eachpend; $p++){
					$hsubov->{$p}{cnt} += 1;
				}
				push(@{$AoAL}, \@A);
			}
			@{$AoAL} = sort {$a->[1] <=> $b->[1]} @{$AoAL};
			
			# seq1	42184785	Sniffles2.DEL.0S0	ATTAAGATGTTTTTTATCAACAAGTTGTTGATTTTAAGTGATGTTGAAGTTGATTCTT	A	60	PASS	PRECISE;SVTYPE=DEL;SVLEN=-56;END=42184842;SUPPORT=100;COVERAGE=100,100,100,100,100;STRAND=+;AF=1.000;STDEV_LEN=0.000;STDEV_POS=0.000	GT:GQ:DR:DV	1/1:60:0:100
			
			my @SubP = keys(%{$hsubov});
			@SubP = sort {$a <=> $b} @SubP;
			
			my $hsubgrp = {};
			my $nsubg = 0;
			my $prev_subp = 0;
			my $nl = 0;
			foreach my $p (@SubP){
				if($nl == 0){
					$hsubov->{$p}{grp} = $nsubg;
					$hsubgrp->{$nsubg}{$p} = 1;
				}
				else{
					if($p - $prev_subp != 1){
						$nsubg++;
					}
					$hsubov->{$p}{grp} = $nsubg;
					$hsubgrp->{$nsubg}{$p} = 1;
				}
				$prev_subp = $p;
				$nl++;
			}
			
			my $hsvjudge = {};
			my @NSUBG = keys(%{$hsubgrp});
			@NSUBG = sort {$a <=> $b} @NSUBG;
			my $AoAL2 = [];
			foreach my $isubg (@NSUBG){
				my $htmp = $hsubgrp->{$isubg};
				my @SubP2 = keys(%{$htmp});
				@SubP2 = sort {$a <=> $b} @SubP2;
				my $num_SubP2 = @SubP2;
#				print " thread [$t] : [$sid] | $SubP2[0] - $SubP2[$num_SubP2 - 1] ($isubg)\n";
				
				my $AoSBP = [];
				foreach my $p (@SubP2){
					if($hsvsize->{$p}{lendf}){
						my @tmp;
						push(@tmp, $p);
						push(@tmp, $hsvsize->{$p}{lendf});
						push(@tmp, $hsvsize->{$p}{data});
						push(@{$AoSBP}, \@tmp);
					}
				}
				@{$AoSBP} = sort {$b->[1] <=> $a->[1]} @{$AoSBP};
				
				my $itop = 0;
				foreach my $A (@{$AoSBP}){
					if($itop == 0){
						$hsvjudge->{$A->[0]} = "pass";		# keep largest SV among overlaping
						my @Atop = split(/\t/, $A->[2]);
						push(@{$AoAL2}, \@Atop);
					}
					else{
						$hsvjudge->{$A->[0]} = "rm";
					}
					$itop++;
				}
			}
			
			my $disable_merge_by_geneunit = 1;
			if(! $disable_merge_by_geneunit){
				my $prev_rend = 0;
				my $start_rpos = 0;
				my $seqR = "";
				my $seqA = "";
				$nl = 0;
				foreach my $A (@{$AoAL2}){
					if(! $hsvjudge->{$A->[1]} || $hsvjudge->{$A->[1]} eq 'rm'){
						$cnt_overlap++;
						next;
					}
					
					my $lenR = length($A->[3]);
					my $lenA = length($A->[4]);
					my $lendf = $lenA - $lenR;
					
					if($nl == 0){
						$seqR .= $A->[3];
						$seqA .= $A->[4];
						$start_rpos = $A->[1];
					}
					else{
						my $lendfx = $A->[1] - $prev_rend;
						if($lendfx < 0){
							print " thread [$t] : [$sid] | error...\n";
							die;
						}
						$seqR .= substr($hdbseq->{$sid}{seq}, $prev_rend - 1, $A->[1] - $prev_rend);
						$seqA .= substr($hdbseq->{$sid}{seq}, $prev_rend - 1, $A->[1] - $prev_rend);
						$seqR .= $A->[3];
						$seqA .= $A->[4];
					}
					
					$prev_rend = $A->[1] + $lenR;
					$cnt_selected++;
					$nl++;
				}
				
				$rvcf .= $sid."\t".$start_rpos."\tSniffles2.combined\t".$seqR."\t".$seqA."\t60\tPASS\t.\tGT:GQ:DR:DV\t1/1:60:0:100\n";
			}
			else{
				foreach my $A (@{$AoAL2}){
					if(! $hsvjudge->{$A->[1]} || $hsvjudge->{$A->[1]} eq 'rm'){
						$cnt_overlap++;
						next;
					}
					
					$rvcf .= join("\t", @{$A})."\n";
				}
			}
		}
	}
	
	print "! thread [$t] : [$sid] [$n_gregion] integrated SV from [$cnt_selected] nr-variants, [$cnt_overlap] removed\n";
}

if($list_remove){
	open(my $rfh, ">", "./log/list_removed.$sid.tsv");
	print $rfh $list_remove;
	close $rfh;
}
if($t == 0){
	$rvcf = $revheader1.$revheader2.$revheader3.$rvcf;
}

SAVE($each_rvcf, $rvcf);

END:{
	my $end = 1;
}


#################################################################################
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

