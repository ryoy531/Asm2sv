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

#---------------------------------------------------------//
my $data_path = getcwd;
unless($data_path){
	print "! cannot find data PATH...\n";
	goto END;
}
my $script_path = $FindBin::Bin;
my $bin_interpro2list = $script_path."/interpro2list.pl";
my $goobo = $script_path."/go.obo.gz";

#---------------------------------------------------------//
my $file_str = shift;
my $reffasta = shift;
my $refgff = shift;
my $datalist = shift;
my $chrinfo = shift;
my $qorf_str = shift;
my $zip = shift;
my $interpro_result = shift;
my $IDalias = shift;
my $svscore_pcnt_th = shift;

unless($file_str){
	print "! missing sv file string...\n";
	goto END;
}
unless($qorf_str){
	print "! missing qorf file string...\n";
	goto END;
}
if(! $reffasta || ! -e $reffasta){
	print "! missing reference fasta...\n";
	goto END;
}
if(! $refgff || ! -e $refgff){
	print "! missing reference Gff...\n";
	goto END;
}
unless($datalist){
	print "! missing datalist...\n";
	goto END;
}
unless($chrinfo){
	print "! missing chrinfo...\n";
	goto END;
}
if(! $interpro_result || ! -e $interpro_result){
	$interpro_result = "null";
}
if(! $IDalias || ! -e $IDalias){
	$IDalias = "null";
}
if(! defined $svscore_pcnt_th){
	$svscore_pcnt_th = 70;
}
elsif($svscore_pcnt_th =~ /\D/){
	$svscore_pcnt_th = 70;
}
elsif($svscore_pcnt_th =~ /\d/){
	if($svscore_pcnt_th < 5 || $svscore_pcnt_th > 95){
		$svscore_pcnt_th = 70;
	}
}

my @PREF = split(/\//, $datalist);
my $num_PREF = @PREF;
$PREF[$num_PREF - 1] =~ s/\.csv//;
$PREF[$num_PREF - 1] =~ s/\.tsv//;
my $prefix = $PREF[$num_PREF - 1];

my $dir = "combined_asm2sv";
if($prefix){
	$dir .= "_".$prefix;
}
unless(-e $dir){
	system("mkdir $dir");
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
my $halias_name = {};
my $hnorder_name = {};
my $norder_status = 0;
if($chrinfo ne 'null' && -e $chrinfo){
	my $htmp_chrinfo = Read_seqinfo($chrinfo);
	$hseqinfo = $htmp_chrinfo->{seqID};
	
	if($htmp_chrinfo->{alias_name}){
		$halias_name = $htmp_chrinfo->{alias_name};
	}
	if($htmp_chrinfo->{norder_status}){
		$hnorder_name = $htmp_chrinfo->{norder};
		$norder_status = 1;
	}
}

my $refgenome = "null";
my $hfasta = {};
if($reffasta ne 'null' && -e $reffasta){
	my @REFS = split(/\//, $reffasta);
	my $nREFS = @REFS;
	$refgenome = $REFS[$nREFS - 1];
	$refgenome =~ s/\.fasta//;
	$refgenome =~ s/\.fa//;
	
	$hfasta = Stats_fasta($reffasta);
}

my $hgff = {};
if($refgff ne 'null' && -e $refgff){
	my $htmp_gff = Gff_to_hash($refgff);
	$hgff = $htmp_gff->{hgff};
}

print "! collecting data...\n";
my $hash = {};
my $neighbor_tlen = 0;
foreach my $tmp (@{$AoF}){
	my $file = $tmp->[0];
	
	$hash = Collectdata($file, $hash, $refgenome, $hseqinfo);
	if($hash->{err}){
		print "! abort script due to error in [$chrinfo]...\n";
		goto END;
	}
}

my $hnrtlen = $hash->{neighbor_tlen};
my @NRtlen = keys(%{$hnrtlen});
my $num_NRtlen = @NRtlen;
if($num_NRtlen != 1){
	print "! abort script due to multiple neighbor_tlen (results may be mixture of different parameter settings)...\n";
	goto END;
}
else{
	$neighbor_tlen = $NRtlen[0];
}
print "! --neighbor = [$neighbor_tlen] bp\n";

my $hgoipr = {};
my $gfile0 = "./$dir/GOterm_list.csv";
my $gfile1 = "./$dir/GOnum_list.csv";
if(-e $interpro_result && -e $IDalias && -e $goobo){
	print "! preparing InterPro info based on [$interpro_result]...\n";
	my $cmd = "$bin_interpro2list $refgff $interpro_result $goobo $IDalias $gfile0 $gfile1 > ./$dir/log_interpro2list.txt 2>&1";
	if(system($cmd) != 0){
		print "! failed.\n";
		$gfile0 = "null";
		$gfile1 = "null";
	}
	elsif(-e $gfile1){
		print "! done.\n";
		$hgoipr = Read_GOdata($gfile1);
	}
}
else{
	print "! skip $bin_interpro2list due to the absence of below files...\n";
	if(! -e $interpro_result){
		print " - [$interpro_result] <not_found>\n";
	}
	if(! -e $IDalias){
		print " - [$IDalias] <not_found>\n";
	}
	if(! -e $goobo){
		print " - [$goobo] <not_found>\n";
	}
}

print "\n! now summarizing data...\n";
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
if($norder_status == 1){
	my $AoN = [];
	my $nrmiss = {};
	foreach my $id (@Samples){
		if(defined $hnorder_name->{$id} && $hnorder_name->{$id} !~ /\D/ && $hnorder_name->{$id} =~ /\d/ && ! $nrmiss->{0}{$id}){
			my @tmp;
			push(@tmp, $id);
			push(@tmp, $hnorder_name->{$id});
			push(@{$AoN}, \@tmp);
		}
		else{
			$nrmiss->{$id} = 1;
		}
	}
	@{$AoN} = sort {$a->[1] <=> $b->[1]} @{$AoN};
	
	my @Sorted_samples0;
	foreach my $tmp (@{$AoN}){
		push(@Sorted_samples0, $tmp->[0]);
	}
	foreach my $id (@Samples){
		if($nrmiss->{$id}){
			push(@Sorted_samples0, $id);
		}
	}
	@Samples = @Sorted_samples0;
}

foreach my $id (@Samples){
	if($refgenome ne 'null'){
		if($id eq $refgenome){
			if($halias_name && $halias_name->{$id}){
				push(@SHeader, "$halias_name->{$id} [reference]");
			}
			else{
				push(@SHeader, "$id [reference]");
			}
			push(@Sorted_samples, $id);
		}
	}
}
foreach my $id (@Samples){
	if($refgenome ne 'null'){
		if($id ne $refgenome){
			if($halias_name && $halias_name->{$id}){
				push(@SHeader, $halias_name->{$id});
			}
			else{
				push(@SHeader, $id);
			}
			push(@Sorted_samples, $id);
		}
	}
	else{
		push(@SHeader, $id);
		push(@Sorted_samples, $id);
	}
}
@Samples = @Sorted_samples;
my $num_samples = @Samples;

#my $r1 = "gene_id,Chr,sp0,sp1,cnt 0-0.3,cnt 0.3-0.6,cnt 0.6-0.9,cnt >=0.9,reference,".join(",", @Samples)."\n";
my $r1 = "gene_id,Chr,pos0,pos1,cnt 0-0.3,cnt 0.3-0.6,cnt 0.6-0.9,cnt >=0.9,".join(",", @SHeader)."\n";
my $r1B = "gene_id,Chr,pos0,pos1,cnt 0-0.6,cnt 0.6-1.2,cnt 1.2-1.8,cnt >=1.8,".join(",", @SHeader)."\n";
my $r1C = "gene_id,class,Chr,pos0,pos1,cnt 0-0.6,cnt 0.6-1.2,cnt 1.2-1.8,cnt >=1.8,".join(",", @SHeader)."\n";
#my $r1_nmatrix = "gene_id,".join(",", @SHeader)."\n";
my $r1_nmatrix = "gene_id,Chr,pos0,pos1,cumpos,".join(",", @SHeader)."\n";
my $r1B_nmatrix = "gene_id,Chr,pos0,pos1,cumpos,".join(",", @SHeader)."\n";
#my $r1 = "gene_id,reference,".join(",", @Samples)."\n";
$r1 =~ s/\.genome//g;
$r1 =~ s/_GCA_022114995_renamed//g;
my $r2 = $r1;
my $r3 = $r1;
my $r4 = $r1;
my $r5 = $r1;
my $r6 = $r1;
my $r7 = $r1;
my $r7D = $r1;
my $r8 = $r1;
my $r9 = $r1;
my $r9D = $r1;
my $r10 = $r1;
my $r1A = $r1B;
my $r1AD = $r1B;
my $r11 = $r1;
my $r11B = $r1B;
my $r11C = $r1C;
my $r11G = $r1C;
my $r12 = $r1;
my $r15 = $r1;
my $r17 = $r1B;
my $r19 = $r1B;
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
#---------------------------------------
my $r1R_nmatrix = $r1B_nmatrix;
my $r1N_nmatrix = $r1B_nmatrix;
my $r1Z_nmatrix = $r1B_nmatrix;
my $r1L_nmatrix = $r1B_nmatrix;
my $r1RZ_nmatrix = $r1B_nmatrix;
my $r1RL_nmatrix = $r1B_nmatrix;
#---------------------------------------
my $r1A_nmatrix = $r1B_nmatrix;
my $r1AD_nmatrix = $r1B_nmatrix;
my $r1AR_nmatrix = $r1B_nmatrix;
my $r1AN_nmatrix = $r1B_nmatrix;
my $r1AZ_nmatrix = $r1B_nmatrix;
my $r1AL_nmatrix = $r1B_nmatrix;
my $r1ARZ_nmatrix = $r1B_nmatrix;
my $r1ARL_nmatrix = $r1B_nmatrix;
#---------------------------------------
my $r7_nmatrix = $r1B_nmatrix;
my $r7R_nmatrix = $r1B_nmatrix;
my $r7N_nmatrix = $r1B_nmatrix;
my $r7Z_nmatrix = $r1B_nmatrix;
my $r7L_nmatrix = $r1B_nmatrix;
my $r7RZ_nmatrix = $r1B_nmatrix;
my $r7RL_nmatrix = $r1B_nmatrix;
#---------------------------------------
my $r9_nmatrix = $r1B_nmatrix;
my $r9R_nmatrix = $r1B_nmatrix;
my $r9N_nmatrix = $r1B_nmatrix;
my $r9Z_nmatrix = $r1B_nmatrix;
my $r9L_nmatrix = $r1B_nmatrix;
my $r9RZ_nmatrix = $r1B_nmatrix;
my $r9RL_nmatrix = $r1B_nmatrix;
#---------------------------------------
my $r11B_nmatrix = $r1B_nmatrix;
my $r11R_nmatrix = $r1B_nmatrix;
my $r11N_nmatrix = $r1B_nmatrix;
my $r11Z_nmatrix = $r1B_nmatrix;
my $r11L_nmatrix = $r1B_nmatrix;
my $r11RZ_nmatrix = $r1B_nmatrix;
my $r11RL_nmatrix = $r1B_nmatrix;
#---------------------------------------
my $r11C_nmatrix = $r1B_nmatrix;
my $r11CR_nmatrix = $r1B_nmatrix;
my $r11CN_nmatrix = $r1B_nmatrix;
my $r11CZ_nmatrix = $r1B_nmatrix;
my $r11CL_nmatrix = $r1B_nmatrix;
my $r11CRZ_nmatrix = $r1B_nmatrix;
my $r11CRL_nmatrix = $r1B_nmatrix;
#---------------------------------------
my $r11G_nmatrix = $r1B_nmatrix;
my $r11GR_nmatrix = $r1B_nmatrix;
my $r11GN_nmatrix = $r1B_nmatrix;
my $r11GZ_nmatrix = $r1B_nmatrix;
my $r11GL_nmatrix = $r1B_nmatrix;
my $r11GRZ_nmatrix = $r1B_nmatrix;
my $r11GRL_nmatrix = $r1B_nmatrix;
#---------------------------------------

my $r_gntable = "gene_id,Chr,pos0,pos1,cumpos,class,num_genotype,ratio_minor,num_each,".join(",", @SHeader)."\n";


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
		my $gid_seqid = $hginfo->{$gid}{Chr};
		
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
		my @Geno1A;
		my @Geno1AD;
		my @Geno2;
		my @Geno3;
		my @Geno4;
		my @Geno5;
		my @Geno4D;
		my @Geno5D;
		my @GenoTable;
		my @GenoTablep;
		my @GenoTableu;
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
			
			if(! $hash->{gntype}{$id}{$gid} || $hash->{gntype}{$id}{$gid} eq '-'){
				$hash->{gntype}{$id}{$gid} = "NA";
			}
			if(! $hash->{gntype_promoter}{$id}{$gid} || $hash->{gntype_promoter}{$id}{$gid} eq '-'){
				$hash->{gntype_promoter}{$id}{$gid} = "NA";
			}
			if(! $hash->{gntype_utr}{$id}{$gid} || $hash->{gntype_utr}{$id}{$gid} eq '-'){
				$hash->{gntype_utr}{$id}{$gid} = "NA";
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
				$hash->{svprom_disrupt}{$id}{$gid} = 0;
				$num_nodata4++;
			}
#			elsif($hash->{svprom}{$id}{$gid} > 1){
#				$hash->{svprom}{$id}{$gid} = 1;
#			}
			
			if(! $hash->{svutr3}{$id}{$gid} || $hash->{svutr3}{$id}{$gid} eq '-'){
				$hash->{svutr3}{$id}{$gid} = 0;
				$hash->{svutr3_disrupt}{$id}{$gid} = 0;
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
			
			if(! defined $hash->{geno_including_flanking}{$id}{$gid}){
				print "! error: missing {geno_including_flanking} infor for $id, $gid\n";
				die;
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
			
			$hash->{sv2_all}{$id}{$gid} = ($hash->{sv2}{$id}{$gid} + $hash->{svprom}{$id}{$gid} + $hash->{svutr3}{$id}{$gid}) / 3;
			if($hash->{sv2_all}{$id}{$gid} > 1){
				$hash->{sv_all}{$id}{$gid} = 1 / $hash->{sv2_all}{$id}{$gid};
			}
			else{
				$hash->{sv_all}{$id}{$gid} = $hash->{sv2_all}{$id}{$gid};
			}
			
			$hash->{sv}{$id}{$gid} = sprintf("%.4f", $hash->{sv}{$id}{$gid});
			$hash->{sv2}{$id}{$gid} = sprintf("%.4f", $hash->{sv2}{$id}{$gid});
			$hash->{sv2_all}{$id}{$gid} = sprintf("%.4f", $hash->{sv2_all}{$id}{$gid});
			$hash->{raw}{$id}{$gid} = sprintf("%.4f", $hash->{raw}{$id}{$gid});
			$hash->{svprot}{$id}{$gid} = sprintf("%.4f", $hash->{svprot}{$id}{$gid});
			$hash->{svprom}{$id}{$gid} = sprintf("%.4f", $hash->{svprom}{$id}{$gid});
			$hash->{svutr3}{$id}{$gid} = sprintf("%.4f", $hash->{svutr3}{$id}{$gid});
			$hash->{svprom_disrupt}{$id}{$gid} = sprintf("%.4f", $hash->{svprom_disrupt}{$id}{$gid});
			$hash->{svutr3_disrupt}{$id}{$gid} = sprintf("%.4f", $hash->{svutr3_disrupt}{$id}{$gid});
			
			push(@Geno1, $hash->{sv}{$id}{$gid});				# sv = disrupt score
			push(@Geno1B, $hash->{sv2}{$id}{$gid});				# sv2 = indel score
			push(@Geno1A, $hash->{sv2_all}{$id}{$gid});			# indel score
			push(@Geno1AD, $hash->{sv_all}{$id}{$gid});			# disrupt score
			push(@Geno2, $hash->{raw}{$id}{$gid});
			push(@Geno3, $hash->{svprot}{$id}{$gid});
			push(@Geno4, $hash->{svprom}{$id}{$gid});
			push(@Geno5, $hash->{svutr3}{$id}{$gid});
			push(@Geno4D, $hash->{svprom_disrupt}{$id}{$gid});
			push(@Geno5D, $hash->{svutr3_disrupt}{$id}{$gid});
			push(@GenoTable, $hash->{gntype}{$id}{$gid});
			push(@GenoTablep, $hash->{gntype_promoter}{$id}{$gid});
			push(@GenoTableu, $hash->{gntype_utr}{$id}{$gid});
			
			if($refgenome ne 'null'){
				if($id eq $refgenome){
					$judge_refgenome = $fgeno;
				}
				
				my $refchr = $hash->{refchr}{$refgenome}{$gid};
				my $refpos0 = $hash->{refpos0}{$refgenome}{$gid};
				my $refpos1 = $hash->{refpos1}{$refgenome}{$gid};
				
				if($hash->{hconvid}{$refgenome}{$gid}){
					$refchr = $hash->{hconvid}{$refgenome}{$gid};
				}
				
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
		($hgeneinfo->{hist_value4}{$gid}, $hgeneinfo->{hist_header4}) = Array2histcsv2(\@Geno4);
		($hgeneinfo->{hist_value5}{$gid}, $hgeneinfo->{hist_header5}) = Array2histcsv2(\@Geno5);
		($hgeneinfo->{hist_value1B}{$gid}, $hgeneinfo->{hist_header1B}) = Array2histcsv2(\@Geno1B);
		($hgeneinfo->{hist_value1A}{$gid}, $hgeneinfo->{hist_header1A}) = Array2histcsv2(\@Geno1A);
		($hgeneinfo->{hist_value1AD}{$gid}, $hgeneinfo->{hist_header1AD}) = Array2histcsv(\@Geno1AD);
		
		$hgeneinfo->{gntable_info}{$gid} = Array2histcsv3(\@GenoTable);
		$hgeneinfo->{gntable_infop}{$gid} = Array2histcsv3(\@GenoTablep);
		$hgeneinfo->{gntable_infou}{$gid} = Array2histcsv3(\@GenoTableu);
		
		$hginfo->{$gid}{pos0_p} = $hginfo->{$gid}{pos0};
		$hginfo->{$gid}{pos1_p} = $hginfo->{$gid}{pos0};
		$hginfo->{$gid}{pos0_u} = $hginfo->{$gid}{pos1};
		$hginfo->{$gid}{pos1_u} = $hginfo->{$gid}{pos1};
		my $cumpos_p = $hginfo->{$gid}{pos0_p} + $hfasta->{$gid_seqid}{cumlen0};
		my $cumpos_g = $hginfo->{$gid}{pos0} + $hfasta->{$gid_seqid}{cumlen0};
		my $cumpos_u = $hginfo->{$gid}{pos0_u} + $hfasta->{$gid_seqid}{cumlen0};
		
		if($hgff->{$gid}{strand}){
			if($hgff->{$gid}{strand} eq '+'){
				$hginfo->{$gid}{pos0_p} = $hginfo->{$gid}{pos0} - $neighbor_tlen - 1;
				$hginfo->{$gid}{pos1_p} = $hginfo->{$gid}{pos0} - 1;
				$hginfo->{$gid}{pos0_u} = $hginfo->{$gid}{pos1} + 1;
				$hginfo->{$gid}{pos1_u} = $hginfo->{$gid}{pos1} + $neighbor_tlen + 1;
				
				if($hginfo->{$gid}{pos0_p} < 1){
					$hginfo->{$gid}{pos0_p} = 1;
				}
				if($hginfo->{$gid}{pos1_p} < 1){
					$hginfo->{$gid}{pos1_p} = 1;
				}
				
				$cumpos_p = $hginfo->{$gid}{pos0_p} + $hfasta->{$gid_seqid}{cumlen0};
				$cumpos_g = $hginfo->{$gid}{pos0} + $hfasta->{$gid_seqid}{cumlen0};
				$cumpos_u = $hginfo->{$gid}{pos0_u} + $hfasta->{$gid_seqid}{cumlen0};
				
				$r_gntable .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$cumpos_p.",p,".$hgeneinfo->{gntable_infop}{$gid}.",".join(",", @GenoTablep)."\n";
				$r_gntable .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",g,".$hgeneinfo->{gntable_info}{$gid}.",".join(",", @GenoTable)."\n";
				$r_gntable .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$cumpos_u.",u,".$hgeneinfo->{gntable_infou}{$gid}.",".join(",", @GenoTableu)."\n";
			}
			elsif($hgff->{$gid}{strand} eq '-'){
				$hginfo->{$gid}{pos0_u} = $hginfo->{$gid}{pos0} - $neighbor_tlen - 1;
				$hginfo->{$gid}{pos1_u} = $hginfo->{$gid}{pos0} - 1;
				$hginfo->{$gid}{pos0_p} = $hginfo->{$gid}{pos1} + 1;
				$hginfo->{$gid}{pos1_p} = $hginfo->{$gid}{pos1} + $neighbor_tlen + 1;
				
				if($hginfo->{$gid}{pos0_u} < 1){
					$hginfo->{$gid}{pos0_u} = 1;
				}
				if($hginfo->{$gid}{pos1_u} < 1){
					$hginfo->{$gid}{pos1_u} = 1;
				}
				
				$cumpos_p = $hginfo->{$gid}{pos0_p} + $hfasta->{$gid_seqid}{cumlen0};
				$cumpos_g = $hginfo->{$gid}{pos0} + $hfasta->{$gid_seqid}{cumlen0};
				$cumpos_u = $hginfo->{$gid}{pos0_u} + $hfasta->{$gid_seqid}{cumlen0};
				
				$r_gntable .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$cumpos_u.",p,".$hgeneinfo->{gntable_infou}{$gid}.",".join(",", @GenoTableu)."\n";
				$r_gntable .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",g,".$hgeneinfo->{gntable_info}{$gid}.",".join(",", @GenoTable)."\n";
				$r_gntable .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$cumpos_p.",u,".$hgeneinfo->{gntable_infop}{$gid}.",".join(",", @GenoTablep)."\n";
			}
		}
		
		$r_hist .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$judge_refgenome.",".$hgeneinfo->{hist_value1}{$gid}.",.,".$hgeneinfo->{hist_value3}{$gid}.",.,".$hgeneinfo->{hist_value4}{$gid}.",.,".$hgeneinfo->{hist_value5}{$gid}."\n";
		
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
					$r101 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value101}{$gid}.",".join(",", @Geno101)."\n";			# disrupt
					$r101B .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value101B}{$gid}.",".join(",", @Geno101B)."\n";		# indel
					$hgeneclass->{transloc_val}{$gid} = 1;
					$cnt_geno101++;
					
					if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
						$r101M .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value101}{$gid}.",".join(",", @Geno101)."\n";
					}
				}
				if($num_nonzero2 > 0){
					$r102 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno102)."\n";
					$hgeneclass->{transloc_raw}{$gid} = 1;
					$cnt_geno102++;
				}
#				if($num_nonzero3 > 0 && $num_nodata3 == 0){
				if($num_nonzero3 > 0 && $failprot_refgenome == 0){
					$r103 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value103}{$gid}.",".join(",", @Geno103)."\n";
					$hgeneclass->{transloc_prot}{$gid} = 1;
					$cnt_geno103++;
				}
#				if($num_nonzero4 > 0 && $num_nodata4 == 0){
				if($num_nonzero4 > 0){
					$r104 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value104}{$gid}.",".join(",", @Geno104)."\n";
					$hgeneclass->{transloc_prom}{$gid} = 1;
					$cnt_geno104++;
				}
#				if($num_nonzero5 > 0 && $num_nodata5 == 0){
				if($num_nonzero5 > 0){
					$r105 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value105}{$gid}.",".join(",", @Geno105)."\n";
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
						$r201 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value201}{$gid}.",".join(",", @Geno201)."\n";
						$r201B .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value201B}{$gid}.",".join(",", @Geno201B)."\n";
						$r201i .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value201}{$gid}.",".join(",", @Geno201i)."\n";
						$hgeneclass->{transloc_val}{$gid} = 1;
						$cnt_geno201++;
						
						if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
							$r201M .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value201}{$gid}.",".join(",", @Geno201)."\n";
						}
					}
					if($num_nonzero2 > 0){
						$r202 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno202)."\n";
						$hgeneclass->{transloc_raw}{$gid} = 1;
						$cnt_geno202++;
					}
	#				if($num_nonzero3 > 0 && $num_nodata3 == 0){
					if($num_nonzero3 > 0 && $failprot_refgenome == 0){
						$r203 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value203}{$gid}.",".join(",", @Geno203)."\n";
						$hgeneclass->{transloc_prot}{$gid} = 1;
						$cnt_geno203++;
					}
	#				if($num_nonzero4 > 0 && $num_nodata4 == 0){
					if($num_nonzero4 > 0){
						$r204 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value204}{$gid}.",".join(",", @Geno204)."\n";
						$hgeneclass->{transloc_prom}{$gid} = 1;
						$cnt_geno204++;
					}
	#				if($num_nonzero5 > 0 && $num_nodata5 == 0){
					if($num_nonzero5 > 0){
						$r205 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value205}{$gid}.",".join(",", @Geno205)."\n";
						$hgeneclass->{transloc_utr3}{$gid} = 1;
						$cnt_geno205++;
					}
				}
				else{
					if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
						$r201M .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,\n";
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
						$r301 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value301}{$gid}.",".join(",", @Geno301)."\n";
						$r301B .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value301B}{$gid}.",".join(",", @Geno301B)."\n";
						$hgeneclass->{transloc_val}{$gid} = 1;
						$cnt_geno301++;
						
						if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
							$r301M .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value301}{$gid}.",".join(",", @Geno301)."\n";
						}
					}
					if($num_nonzero2 > 0){
						$r302 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno302)."\n";
						$hgeneclass->{transloc_raw}{$gid} = 1;
						$cnt_geno302++;
					}
	#				if($num_nonzero3 > 0 && $num_nodata3 == 0){
					if($num_nonzero3 > 0 && $failprot_refgenome == 0){
						$r303 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value303}{$gid}.",".join(",", @Geno303)."\n";
						$hgeneclass->{transloc_prot}{$gid} = 1;
						$cnt_geno303++;
					}
	#				if($num_nonzero4 > 0 && $num_nodata4 == 0){
					if($num_nonzero4 > 0){
						$r304 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value304}{$gid}.",".join(",", @Geno304)."\n";
						$hgeneclass->{transloc_prom}{$gid} = 1;
						$cnt_geno304++;
					}
	#				if($num_nonzero5 > 0 && $num_nodata5 == 0){
					if($num_nonzero5 > 0){
						$r305 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value305}{$gid}.",".join(",", @Geno305)."\n";
						$hgeneclass->{transloc_utr3}{$gid} = 1;
						$cnt_geno305++;
					}
				}
				else{
					if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
						$r301M .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,\n";
					}
				}
			}
			else{
				if($num_conserved > 0 && $num_misgenotype == 0 && $num_nodata1 == 0){
					$r101M .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,\n";
					$r201M .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,\n";
					$r301M .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,\n";
				}
			}
		}
		
		my $hrzval = Calc_normalizedvals(\@Geno1B);
		my $Geno1R = $hrzval->{R};
		my $Geno1N = $hrzval->{N};
		my $Geno1Z = $hrzval->{Z};
		my $Geno1L = $hrzval->{L};
		my $Geno1RZ = $hrzval->{RZ};
		my $Geno1RL = $hrzval->{RL};
		
		my $hrzvalA = Calc_normalizedvals(\@Geno1A);
		my $Geno1AR = $hrzvalA->{R};
		my $Geno1AN = $hrzvalA->{N};
		my $Geno1AZ = $hrzvalA->{Z};
		my $Geno1AL = $hrzvalA->{L};
		my $Geno1ARZ = $hrzvalA->{RZ};
		my $Geno1ARL = $hrzvalA->{RL};
		
#		if($num_nonzero1 > 0 && $num_nodata1 == 0){
		if($num_nonzero1 > 0){
			$r1 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1}{$gid}.",".join(",", @Geno1)."\n";
			$r1_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno1)."\n";
			$r1B .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1B}{$gid}.",".join(",", @Geno1B)."\n";
			$r1B_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno1B)."\n";
			$r1A .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1A}{$gid}.",".join(",", @Geno1A)."\n";
			$r1A_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno1A)."\n";
			$r1AD .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1A}{$gid}.",".join(",", @Geno1AD)."\n";
			$r1AD_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno1AD)."\n";
			
			$r1R_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1R})."\n";
			$r1N_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1N})."\n";
			$r1Z_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1Z})."\n";
			$r1L_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1L})."\n";
			$r1RZ_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1RZ})."\n";
			$r1RL_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1RL})."\n";
			$r1AR_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1AR})."\n";
			$r1AN_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1AN})."\n";
			$r1AZ_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1AZ})."\n";
			$r1AL_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1AL})."\n";
			$r1ARZ_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1ARZ})."\n";
			$r1ARL_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1ARL})."\n";
			$hgeneclass->{val}{$gid} = 1;
			$cnt_geno1++;
		}
		else{
			$cnt_geno1_rm++;
		}
		if($num_nonzero2 > 0){
			$r2 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno2)."\n";
			$hgeneclass->{raw}{$gid} = 1;
			$cnt_geno2++;
		}
#		if($num_nonzero3 > 0 && $num_nodata3 == 0){
		if($num_nonzero3 > 0 && $failprot_refgenome == 0){
			$r5 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value3}{$gid}.",".join(",", @Geno3)."\n";
			$hgeneclass->{prot}{$gid} = 1;
			$cnt_geno5++;
		}
#		if($num_nonzero4 > 0 && $num_nodata4 == 0){
		if($num_nonzero4 > 0){
			$r7 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$hgeneinfo->{hist_value4}{$gid}.",".join(",", @Geno4)."\n";
			$r7D .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$hgeneinfo->{hist_value4}{$gid}.",".join(",", @Geno4D)."\n";
			$r7_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno4)."\n";
			$hgeneclass->{prom}{$gid} = 1;
			$cnt_geno7++;
		}
#		if($num_nonzero5 > 0 && $num_nodata5 == 0){
		if($num_nonzero5 > 0){
			$r9 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$hgeneinfo->{hist_value5}{$gid}.",".join(",", @Geno5)."\n";
			$r9D .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$hgeneinfo->{hist_value5}{$gid}.",".join(",", @Geno5D)."\n";
			$r9_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno5)."\n";
			$hgeneclass->{utr3}{$gid} = 1;
			$cnt_geno9++;
		}
		
		if($num_conserved > 0 && $num_misgenotype == 0){		# double check of no data with $num_misgenotype == 0
			if($num_nonzero1 > 0 && $num_nodata1 == 0){
				$r11 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1}{$gid}.",".join(",", @Geno1)."\n";
				$r11_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno1)."\n";
				$r11B .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1B}{$gid}.",".join(",", @Geno1B)."\n";
				$r11B_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno1B)."\n";
				$r11C .= $gid.",gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1B}{$gid}.",".join(",", @Geno1B)."\n";
				$r11C_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno1B)."\n";
				$r11G .= $gid.",gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1B}{$gid}.",".join(",", @Geno1B)."\n";
				$r11G_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno1B)."\n";
				
				$r11R_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1R})."\n";
				$r11N_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1N})."\n";
				$r11Z_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1Z})."\n";
				$r11L_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1L})."\n";
				$r11RZ_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1RZ})."\n";
				$r11RL_nmatrix .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1RL})."\n";
				$r11CR_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1R})."\n";
				$r11CN_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1N})."\n";
				$r11CZ_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1Z})."\n";
				$r11CL_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1L})."\n";
				$r11CRZ_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1RZ})."\n";
				$r11CRL_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1RL})."\n";
				$r11GR_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1R})."\n";
				$r11GN_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1N})."\n";
				$r11GZ_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1Z})."\n";
				$r11GL_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1L})."\n";
				$r11GRZ_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1RZ})."\n";
				$r11GRL_nmatrix .= $gid."_gene,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno1RL})."\n";
				$hgeneclass->{val_allgenotyped}{$gid} = 1;
				$cnt_geno11++;
			}
			else{
				$cnt_geno11_rm++;
			}
			if($num_nonzero2 > 0){
				$r12 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno2)."\n";
				$hgeneclass->{raw_allgenotyped}{$gid} = 1;
				$cnt_geno12++;
			}
#			if($num_nonzero3 > 0 && $num_nodata3 == 0 && $failprot_refgenome == 0){
			if($num_nonzero3 > 0 && $failprot_refgenome == 0){
				$r15 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value3}{$gid}.",".join(",", @Geno3)."\n";
				$r15_nmatrix .= $gid.",".join(",", @Geno3)."\n";
				$r11C .= $gid.",protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value3}{$gid}.",".join(",", @Geno3)."\n";
				$r11C_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno3)."\n";
				$r11G .= $gid.",protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value3}{$gid}.",".join(",", @Geno3)."\n";
				$r11G_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @Geno3)."\n";
				
				my $hrzval = Calc_normalizedvals(\@Geno3);
				my $Geno3R = $hrzval->{R};
				my $Geno3N = $hrzval->{N};
				my $Geno3Z = $hrzval->{Z};
				my $Geno3L = $hrzval->{L};
				my $Geno3RZ = $hrzval->{RZ};
				my $Geno3RL = $hrzval->{RL};
#				$r11CR_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3R})."\n";
				$r11CN_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3N})."\n";
				$r11CZ_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3Z})."\n";
				$r11CL_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3L})."\n";
				$r11CRZ_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3RZ})."\n";
				$r11CRL_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3RL})."\n";
				$r11GR_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3R})."\n";
				$r11GN_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3N})."\n";
				$r11GZ_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3Z})."\n";
				$r11GL_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3L})."\n";
				$r11GRZ_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3RZ})."\n";
				$r11GRL_nmatrix .= $gid."_protein,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno3RL})."\n";
				$hgeneclass->{prot_allgenotyped}{$gid} = 1;
				$cnt_geno15++;
			}
		}
		if($num_conserved_prom > 0 && $num_nonzero4 > 0 && $num_nodata4 == 0){
			$r17 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$hgeneinfo->{hist_value4}{$gid}.",".join(",", @Geno4)."\n";
			$r17_nmatrix .= $gid.",".join(",", @Geno4)."\n";
			$r11C .= $gid.",promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$hgeneinfo->{hist_value4}{$gid}.",".join(",", @Geno4)."\n";
			$r11C_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$cumpos_p.",".join(",", @Geno4)."\n";
			
			my $hrzval = Calc_normalizedvals(\@Geno4);
			my $Geno4R = $hrzval->{R};
			my $Geno4N = $hrzval->{N};
			my $Geno4Z = $hrzval->{Z};
			my $Geno4L = $hrzval->{L};
			my $Geno4RZ = $hrzval->{RZ};
			my $Geno4RL = $hrzval->{RL};
			$r7R_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno4R})."\n";
			$r7N_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno4N})."\n";
			$r7Z_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno4Z})."\n";
			$r7L_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno4L})."\n";
			$r7RZ_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno4RZ})."\n";
			$r7RL_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno4RL})."\n";
			
			$r11CR_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$cumpos_p.",".join(",", @{$Geno4R})."\n";
			$r11CN_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$cumpos_p.",".join(",", @{$Geno4N})."\n";
			$r11CZ_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$cumpos_p.",".join(",", @{$Geno4Z})."\n";
			$r11CL_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$cumpos_p.",".join(",", @{$Geno4L})."\n";
			$r11CRZ_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$cumpos_p.",".join(",", @{$Geno4RZ})."\n";
			$r11CRL_nmatrix .= $gid."_promoter,".$gid_seqid.",".$hginfo->{$gid}{pos0_p}.",".$hginfo->{$gid}{pos1_p}.",".$cumpos_p.",".join(",", @{$Geno4RL})."\n";
			$hgeneclass->{prom_allgenotyped}{$gid} = 1;
			$cnt_geno17++;
		}
		if($num_conserved_utr3 > 0 && $num_nonzero5 > 0 && $num_nodata5 == 0){
			$r19 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$hgeneinfo->{hist_value5}{$gid}.",".join(",", @Geno5)."\n";
			$r19_nmatrix .= $gid.",".join(",", @Geno5)."\n";
			$r11C .= $gid.",utr,".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$hgeneinfo->{hist_value5}{$gid}.",".join(",", @Geno5)."\n";
			$r11C_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$cumpos_u.",".join(",", @Geno5)."\n";
			
			my $hrzval = Calc_normalizedvals(\@Geno5);
			my $Geno5R = $hrzval->{R};
			my $Geno5N = $hrzval->{N};
			my $Geno5Z = $hrzval->{Z};
			my $Geno5L = $hrzval->{L};
			my $Geno5RZ = $hrzval->{RZ};
			my $Geno5RL = $hrzval->{RL};
			$r9R_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno5R})."\n";
			$r9N_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno5N})."\n";
			$r9Z_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno5Z})."\n";
			$r9L_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno5L})."\n";
			$r9RZ_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno5RZ})."\n";
			$r9RL_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$cumpos_g.",".join(",", @{$Geno5RL})."\n";
			
			$r11CR_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$cumpos_u.",".join(",", @{$Geno5R})."\n";
			$r11CN_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$cumpos_u.",".join(",", @{$Geno5N})."\n";
			$r11CZ_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$cumpos_u.",".join(",", @{$Geno5Z})."\n";
			$r11CL_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$cumpos_u.",".join(",", @{$Geno5L})."\n";
			$r11CRZ_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$cumpos_u.",".join(",", @{$Geno5RZ})."\n";
			$r11CRL_nmatrix .= $gid."_utr,".$gid_seqid.",".$hginfo->{$gid}{pos0_u}.",".$hginfo->{$gid}{pos1_u}.",".$cumpos_u.",".join(",", @{$Geno5RL})."\n";
			$hgeneclass->{utr3_allgenotyped}{$gid} = 1;
			$cnt_geno19++;
		}
		
		if(! $hjudge->{false}){			# may be removed
			if($num_nonzero1 > 0 && $num_nodata1 == 0){
				$r3 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value1}{$gid}.",".join(",", @Geno1)."\n";
				$hgeneclass->{val_allgeno}{$gid} = 1;
				$cnt_geno3++;
			}
			if($num_nonzero2 > 0){
				$r4 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",-,-,-,-,".join(",", @Geno2)."\n";
				$hgeneclass->{raw_allgeno}{$gid} = 1;
				$cnt_geno4++;
			}
			if($num_nonzero3 > 0 && $num_nodata3 == 0 && $failprot_refgenome == 0){
				$r6 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value3}{$gid}.",".join(",", @Geno3)."\n";
				$hgeneclass->{prot_allgeno}{$gid} = 1;
				$cnt_geno6++;
			}
			if($num_nonzero4 > 0 && $num_nodata4 == 0){
				$r8 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value4}{$gid}.",".join(",", @Geno4)."\n";
				$hgeneclass->{prom_allgeno}{$gid} = 1;
				$cnt_geno8++;
			}
			if($num_nonzero5 > 0 && $num_nodata5 == 0){
				$r10 .= $gid.",".$gid_seqid.",".$hginfo->{$gid}{pos0}.",".$hginfo->{$gid}{pos1}.",".$hgeneinfo->{hist_value5}{$gid}.",".join(",", @Geno5)."\n";
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

my $svscore_nsample_th = int($num_samples * $svscore_pcnt_th / 100);

print "! saving data...\n";
my $pseudogoipr = {};
my $r1_lessConserved = "";
my $r1AD_lessConserved = "";
my $r2_lessConserved = "";
my $r5_lessConserved = "";
my $r7_lessConserved = "";
my $r9_lessConserved = "";
#my $r3_lessConserved = "";
#my $r4_lessConserved = "";
#my $r6_lessConserved = "";
#my $r8_lessConserved = "";
#my $r10_lessConserved = "";
my $r11_lessConserved = "";
my $r12_lessConserved = "";
my $r15_lessConserved = "";
my $r17_lessConserved = "";
my $r19_lessConserved = "";
($r1_lessConserved, $hgeneclass) = Select_lessConverved($r1, $hgeneclass, "val_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
($r1AD_lessConserved, $hgeneclass) = Select_lessConverved($r1AD, $hgeneclass, "val_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
($r2_lessConserved, $hgeneclass) = Select_lessConverved($r2, $hgeneclass, "raw_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
($r5_lessConserved, $hgeneclass) = Select_lessConverved($r5, $hgeneclass, "prot_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
($r7_lessConserved, $hgeneclass) = Select_lessConverved($r7D, $hgeneclass, "promoter_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
($r9_lessConserved, $hgeneclass) = Select_lessConverved($r9D, $hgeneclass, "utr3_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
#($r3_lessConserved, $hgeneclass) = Select_lessConverved($r3, $hgeneclass, "val_allgeno_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
#($r4_lessConserved, $hgeneclass) = Select_lessConverved($r4, $hgeneclass, "raw_allgeno_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
#($r6_lessConserved, $hgeneclass) = Select_lessConverved($r6, $hgeneclass, "prot_allgeno_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
#($r8_lessConserved, $hgeneclass) = Select_lessConverved($r8, $hgeneclass, "promoter_allgeno_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
#($r10_lessConserved, $hgeneclass) = Select_lessConverved($r10, $hgeneclass, "utr3_allgeno_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
($r11_lessConserved, $hgeneclass) = Select_lessConverved($r11, $hgeneclass, "val_allgenotyped_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
($r12_lessConserved, $hgeneclass) = Select_lessConverved($r12, $hgeneclass, "raw_allgenotyped_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
($r15_lessConserved, $hgeneclass) = Select_lessConverved($r15, $hgeneclass, "prot_allgenotyped_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
($r17_lessConserved, $hgeneclass) = Select_lessConverved($r17, $hgeneclass, "promoter_allgenotyped_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);
($r19_lessConserved, $hgeneclass) = Select_lessConverved($r19, $hgeneclass, "utr3_allgenotyped_less_$svscore_nsample_th", $pseudogoipr, $svscore_nsample_th);

my $r1_lessConservedwoTRs = "";
my $r1AD_lessConservedwoTRs = "";
my $r2_lessConservedwoTRs = "";
my $r5_lessConservedwoTRs = "";
my $r7_lessConservedwoTRs = "";
my $r9_lessConservedwoTRs = "";
#my $r3_lessConservedwoTRs = "";
#my $r4_lessConservedwoTRs = "";
#my $r6_lessConservedwoTRs = "";
#my $r8_lessConservedwoTRs = "";
#my $r10_lessConservedwoTRs = "";
my $r11_lessConservedwoTRs = "";
my $r12_lessConservedwoTRs = "";
my $r15_lessConservedwoTRs = "";
my $r17_lessConservedwoTRs = "";
my $r19_lessConservedwoTRs = "";
($r1_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r1, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
($r1AD_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r1AD, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
($r2_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r2, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
($r5_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r5, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
($r7_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r7D, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
($r9_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r9D, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
#($r3_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r3, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
#($r4_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r4, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
#($r6_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r6, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
#($r8_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r8, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
#($r10_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r10, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
($r11_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r11, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
($r12_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r12, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
($r15_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r15, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
($r17_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r17, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);
($r19_lessConservedwoTRs, $hgeneclass) = Select_lessConverved($r19, $hgeneclass, "null", $hgoipr, $svscore_nsample_th);

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

my $rfile1 = "$dir/gene_disrupt_q-".$nsample.".csv";
my $rfile1B = "$dir/gene_indel_q-".$nsample.".csv";
my $rfile1A = "$dir/sum_indel_q-".$nsample.".csv";
my $rfile1AD = "$dir/sum_disrupt_q-".$nsample.".csv";
my $rfile2 = "$dir/raw_q-".$nsample.".csv";
my $rfile3 = "$dir/gene_disrupt_allgeno_q-".$nsample.".csv";
my $rfile4 = "$dir/raw_allgeno_q-".$nsample.".csv";
my $rfile5 = "$dir/protein_q-".$nsample.".csv";
my $rfile7 = "$dir/promoter_indel_q-".$nsample.".csv";
my $rfile7D = "$dir/promoter_disrupt_q-".$nsample.".csv";
my $rfile9 = "$dir/3UTR_indel_q-".$nsample.".csv";
my $rfile9D = "$dir/3UTR_disrupt_q-".$nsample.".csv";
my $rfile6 = "$dir/protein_allgeno_q-".$nsample.".csv";
my $rfile8 = "$dir/promoter_indel_allgeno_q-".$nsample.".csv";
my $rfile10 = "$dir/3UTR_indel_allgeno_q-".$nsample.".csv";
my $rfile11 = "$dir/gene_disrupt_allgenotyped_q-".$nsample.".csv";
my $rfile11B = "$dir/gene_indel_allgenotyped_q-".$nsample.".csv";
my $rfile11C = "$dir/all_indel_allgenotyped_q-".$nsample.".csv";
my $rfile12 = "$dir/raw_allgenotyped_q-".$nsample.".csv";
my $rfile15 = "$dir/protein_allgenotyped_q-".$nsample.".csv";
my $rfile17 = "$dir/promoter_indel_allgenotyped_q-".$nsample.".csv";
my $rfile19 = "$dir/3UTR_indel_allgenotyped_q-".$nsample.".csv";
my $rfile1l = "$dir/gene_disrupt_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile1ADl = "$dir/sum_disrupt_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile2l = "$dir/raw_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile5l = "$dir/protein_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile7l = "$dir/promoter_disrupt_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile9l = "$dir/3UTR_disrupt_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile3l = "$dir/gene_disrupt_allgeno_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile4l = "$dir/raw_allgeno_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile6l = "$dir/protein_allgeno_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile8l = "$dir/promoter_indel_allgeno_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile10l = "$dir/3UTR_indel_allgeno_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile11l = "$dir/gene_disrupt_allgenotyped_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile12l = "$dir/raw_allgenotyped_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile15l = "$dir/protein_allgenotyped_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile17l = "$dir/promoter_indel_allgenotyped_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile19l = "$dir/3UTR_indel_allgenotyped_q-".$nsample."_lessConserved_th-".$svscore_nsample_th.".csv";
my $rfile1lTRs = "$dir/gene_disrupt_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile1ADlTRs = "$dir/sum_disrupt_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile2lTRs = "$dir/raw_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile5lTRs = "$dir/protein_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile7lTRs = "$dir/promoter_disrupt_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile9lTRs = "$dir/3UTR_disrupt_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile3lTRs = "$dir/gene_disrupt_allgeno_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile4lTRs = "$dir/raw_allgeno_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile6lTRs = "$dir/protein_allgeno_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile8lTRs = "$dir/promoter_indel_allgeno_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile10lTRs = "$dir/3UTR_indel_allgeno_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile11lTRs = "$dir/gene_disrupt_allgenotyped_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile12lTRs = "$dir/raw_allgenotyped_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile15lTRs = "$dir/protein_allgenotyped_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile17lTRs = "$dir/promoter_indel_allgenotyped_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rfile19lTRs = "$dir/3UTR_indel_allgenotyped_q-".$nsample."_lessConserved_th-".$svscore_nsample_th."_woTransposon.csv";
my $rginfo = "$dir/geneinfo_q-".$nsample.".csv";
#-------------------------------
my $rfile1_nmatrix = "$dir/gene_disrupt_q-".$nsample."_num-matrix.csv";
my $rfile1AD_nmatrix = "$dir/sum_disrupt_q-".$nsample."_num-matrix.csv";
my $rfile1B_nmatrix = "$dir/gene_indel_q-".$nsample."_num-matrix.csv";
my $rfile1V_nmatrix = "$dir/forClust_gene_indel_q-".$nsample."_native.csv";
my $rfile1R_nmatrix = "$dir/forClust_gene_indel_q-".$nsample."_relative.csv";
my $rfile1N_nmatrix = "$dir/forClust_gene_indel_q-".$nsample."_averaged.csv";
my $rfile1L_nmatrix = "$dir/forClust_gene_indel_q-".$nsample."_log.csv";
my $rfile1Z_nmatrix = "$dir/forClust_gene_indel_q-".$nsample."_zscore.csv";
my $rfile1RL_nmatrix = "$dir/forClust_gene_indel_q-".$nsample."_log-relative.csv";
my $rfile1RZ_nmatrix = "$dir/forClust_gene_indel_q-".$nsample."_zscore-relative.csv";
my $rfile1AV_nmatrix = "$dir/forClust_sum_indel_q-".$nsample."_native.csv";
my $rfile1AR_nmatrix = "$dir/forClust_sum_indel_q-".$nsample."_relative.csv";
my $rfile1AL_nmatrix = "$dir/forClust_sum_indel_q-".$nsample."_log.csv";
#-------------------------------
my $rfile11_nmatrix = "$dir/gene_disrupt_allgenotyped_q-".$nsample."_num-matrix.csv";
my $rfile11B_nmatrix = "$dir/gene_indel_allgenotyped_q-".$nsample."_num-matrix.csv";
my $rfile11V_nmatrix = "$dir/forClust_gene_indel_allgenotyped_q-".$nsample."_native.csv";
my $rfile11R_nmatrix = "$dir/forClust_gene_indel_allgenotyped_q-".$nsample."_relative.csv";
my $rfile11N_nmatrix = "$dir/forClust_gene_indel_allgenotyped_q-".$nsample."_averaged.csv";
my $rfile11L_nmatrix = "$dir/forClust_gene_indel_allgenotyped_q-".$nsample."_log.csv";
my $rfile11Z_nmatrix = "$dir/forClust_gene_indel_allgenotyped_q-".$nsample."_zscore.csv";
my $rfile11RL_nmatrix = "$dir/forClust_gene_indel_allgenotyped_q-".$nsample."_log-relative.csv";
my $rfile11RZ_nmatrix = "$dir/forClust_gene_indel_allgenotyped_q-".$nsample."_zscore-relative.csv";
#-------------------------------
my $rfile7_nmatrix = "$dir/promoter_disrupt_q-".$nsample."_num-matrix.csv";
my $rfile7B_nmatrix = "$dir/promoter_indel_q-".$nsample."_num-matrix.csv";
my $rfile7V_nmatrix = "$dir/forClust_promoter_indel_q-".$nsample."_native.csv";
my $rfile7R_nmatrix = "$dir/forClust_promoter_indel_q-".$nsample."_relative.csv";
my $rfile7N_nmatrix = "$dir/forClust_promoter_indel_q-".$nsample."_averaged.csv";
my $rfile7L_nmatrix = "$dir/forClust_promoter_indel_q-".$nsample."_log.csv";
my $rfile7Z_nmatrix = "$dir/forClust_promoter_indel_q-".$nsample."_zscore.csv";
my $rfile7RL_nmatrix = "$dir/forClust_promoter_indel_q-".$nsample."_log-relative.csv";
my $rfile7RZ_nmatrix = "$dir/forClust_promoter_indel_q-".$nsample."_zscore-relative.csv";
#-------------------------------
my $rfile9_nmatrix = "$dir/3UTR_disrupt_q-".$nsample."_num-matrix.csv";
my $rfile9B_nmatrix = "$dir/3UTR_indel_q-".$nsample."_num-matrix.csv";
my $rfile9V_nmatrix = "$dir/forClust_3UTR_indel_q-".$nsample."_native.csv";
my $rfile9R_nmatrix = "$dir/forClust_3UTR_indel_q-".$nsample."_relative.csv";
my $rfile9N_nmatrix = "$dir/forClust_3UTR_indel_q-".$nsample."_averaged.csv";
my $rfile9L_nmatrix = "$dir/forClust_3UTR_indel_q-".$nsample."_log.csv";
my $rfile9Z_nmatrix = "$dir/forClust_3UTR_indel_q-".$nsample."_zscore.csv";
my $rfile9RL_nmatrix = "$dir/forClust_3UTR_indel_q-".$nsample."_log-relative.csv";
my $rfile9RZ_nmatrix = "$dir/forClust_3UTR_indel_q-".$nsample."_zscore-relative.csv";
#-------------------------------
my $rfile11C_nmatrix = "$dir/all_indel_allgenotyped_q-".$nsample."_num-matrix.csv";
my $rfile11CV_nmatrix = "$dir/forClust_all_indel_allgenotyped_q-".$nsample."_native.csv";
my $rfile11CR_nmatrix = "$dir/forClust_all_indel_allgenotyped_q-".$nsample."_relative.csv";
my $rfile11CN_nmatrix = "$dir/forClust_all_indel_allgenotyped_q-".$nsample."_averaged.csv";
my $rfile11CL_nmatrix = "$dir/forClust_all_indel_allgenotyped_q-".$nsample."_log.csv";
my $rfile11CZ_nmatrix = "$dir/forClust_all_indel_allgenotyped_q-".$nsample."_zscore.csv";
my $rfile11CRL_nmatrix = "$dir/forClust_all_indel_allgenotyped_q-".$nsample."_log-relative.csv";
my $rfile11CRZ_nmatrix = "$dir/forClust_all_indel_allgenotyped_q-".$nsample."_zscore-relative.csv";
#-------------------------------
my $rfile11G_nmatrix = "$dir/protein-gene_indel_allgenotyped_q-".$nsample."_num-matrix.csv";
my $rfile11GV_nmatrix = "$dir/forClust_protein-gene_indel_allgenotyped_q-".$nsample."_native.csv";
my $rfile11GR_nmatrix = "$dir/forClust_protein-gene_indel_allgenotyped_q-".$nsample."_relative.csv";
my $rfile11GN_nmatrix = "$dir/forClust_protein-gene_indel_allgenotyped_q-".$nsample."_averaged.csv";
my $rfile11GL_nmatrix = "$dir/forClust_protein-gene_indel_allgenotyped_q-".$nsample."_log.csv";
my $rfile11GZ_nmatrix = "$dir/forClust_protein-gene_indel_allgenotyped_q-".$nsample."_zscore.csv";
my $rfile11GRL_nmatrix = "$dir/forClust_protein-gene_indel_allgenotyped_q-".$nsample."_log-relative.csv";
my $rfile11GRZ_nmatrix = "$dir/forClust_protein-gene_indel_allgenotyped_q-".$nsample."_zscore-relative.csv";
#-------------------------------
my $rfile15_nmatrix = "$dir/protein_allgenotyped_q-".$nsample."_num-matrix.csv";
my $rfile17_nmatrix = "$dir/promoter_indel_allgenotyped_q-".$nsample."_num-matrix.csv";
my $rfile19_nmatrix = "$dir/3UTR_indel_allgenotyped_q-".$nsample."_num-matrix.csv";
my $rgntable = "$dir/genotype_table_q-".$nsample.".csv";

SAVE($rfile1, $r1);
SAVE($rfile1B, $r1B);
SAVE($rfile1A, $r1A);
SAVE($rfile1AD, $r1AD);
SAVE($rfile2, $r2);
#SAVE($rfile3, $r3);
#SAVE($rfile4, $r4);
SAVE($rfile5, $r5);
#SAVE($rfile6, $r6);
SAVE($rfile7, $r7);
SAVE($rfile7D, $r7D);
#SAVE($rfile8, $r8);
SAVE($rfile9, $r9);
SAVE($rfile9D, $r9D);
#SAVE($rfile10, $r10);
#SAVE($rfile11, $r11);
#SAVE($rfile11B, $r11B);
#SAVE($rfile11C, $r11C);
#SAVE($rfile12, $r12);
#SAVE($rfile15, $r15);
#SAVE($rfile17, $r17);
#SAVE($rfile19, $r19);

if(-e $gfile1){
	SAVE($rfile1l, $r1_lessConserved);
	SAVE($rfile1ADl, $r1AD_lessConserved);
	SAVE($rfile2l, $r2_lessConserved);
	SAVE($rfile5l, $r5_lessConserved);
	SAVE($rfile7l, $r7_lessConserved);
	SAVE($rfile9l, $r9_lessConserved);
	#SAVE($rfile3l, $r3_lessConserved);
	#SAVE($rfile4l, $r4_lessConserved);
	#SAVE($rfile6l, $r6_lessConserved);
#	SAVE($rfile11l, $r11_lessConserved);
#	SAVE($rfile12l, $r12_lessConserved);
#	SAVE($rfile15l, $r15_lessConserved);
#	SAVE($rfile17l, $r17_lessConserved);
#	SAVE($rfile19l, $r19_lessConserved);
	SAVE($rfile1lTRs, $r1_lessConservedwoTRs);
	SAVE($rfile1ADlTRs, $r1AD_lessConservedwoTRs);
	SAVE($rfile2lTRs, $r2_lessConservedwoTRs);
	SAVE($rfile5lTRs, $r5_lessConservedwoTRs);
	SAVE($rfile7lTRs, $r7_lessConservedwoTRs);
	SAVE($rfile9lTRs, $r9_lessConservedwoTRs);
	#SAVE($rfile4lTRs, $r4_lessConservedwoTRs);
	#SAVE($rfile6lTRs, $r6_lessConservedwoTRs);
#	SAVE($rfile11lTRs, $r11_lessConservedwoTRs);
#	SAVE($rfile12lTRs, $r12_lessConservedwoTRs);
#	SAVE($rfile15lTRs, $r15_lessConservedwoTRs);
#	SAVE($rfile17lTRs, $r17_lessConservedwoTRs);
#	SAVE($rfile19lTRs, $r19_lessConservedwoTRs);
}
#-------------------------------
SAVE($rfile1_nmatrix, $r1_nmatrix);
SAVE($rfile1AD_nmatrix, $r1AD_nmatrix);
#SAVE($rfile1B_nmatrix, $r1B_nmatrix);
SAVE($rfile1V_nmatrix, $r1B_nmatrix);
SAVE($rfile1R_nmatrix, $r1R_nmatrix);
#SAVE($rfile1N_nmatrix, $r1N_nmatrix);
SAVE($rfile1L_nmatrix, $r1L_nmatrix);
#SAVE($rfile1Z_nmatrix, $r1Z_nmatrix);
#SAVE($rfile1RL_nmatrix, $r1RL_nmatrix);
#SAVE($rfile1RZ_nmatrix, $r1RZ_nmatrix);
SAVE($rfile1AV_nmatrix, $r1A_nmatrix);
SAVE($rfile1AR_nmatrix, $r1AR_nmatrix);
SAVE($rfile1AL_nmatrix, $r1AL_nmatrix);
#-------------------------------
SAVE($rfile7_nmatrix, $r7D);
#SAVE($rfile7B_nmatrix, $r7B_nmatrix);
SAVE($rfile7V_nmatrix, $r7_nmatrix);
SAVE($rfile7R_nmatrix, $r7R_nmatrix);
#SAVE($rfile7N_nmatrix, $r7N_nmatrix);
SAVE($rfile7L_nmatrix, $r7L_nmatrix);
#SAVE($rfile7Z_nmatrix, $r7Z_nmatrix);
#SAVE($rfile7RL_nmatrix, $r7RL_nmatrix);
#SAVE($rfile7RZ_nmatrix, $r7RZ_nmatrix);
#-------------------------------
SAVE($rfile9_nmatrix, $r9D);
#SAVE($rfile9B_nmatrix, $r9B_nmatrix);
SAVE($rfile9V_nmatrix, $r9_nmatrix);
SAVE($rfile9R_nmatrix, $r9R_nmatrix);
#SAVE($rfile9N_nmatrix, $r9N_nmatrix);
SAVE($rfile9L_nmatrix, $r9L_nmatrix);
#SAVE($rfile9Z_nmatrix, $r9Z_nmatrix);
#SAVE($rfile9RL_nmatrix, $r9RL_nmatrix);
#SAVE($rfile9RZ_nmatrix, $r9RZ_nmatrix);
#-------------------------------
#SAVE($rfile11_nmatrix, $r11_nmatrix);
#SAVE($rfile11B_nmatrix, $r11B_nmatrix);
#SAVE($rfile11V_nmatrix, $r11B_nmatrix);
#SAVE($rfile11R_nmatrix, $r11R_nmatrix);
#SAVE($rfile11N_nmatrix, $r11N_nmatrix);
#SAVE($rfile11L_nmatrix, $r11L_nmatrix);
#SAVE($rfile11Z_nmatrix, $r11Z_nmatrix);
#SAVE($rfile11RL_nmatrix, $r11RL_nmatrix);
#SAVE($rfile11RZ_nmatrix, $r11RZ_nmatrix);
#-------------------------------
#SAVE($rfile11C_nmatrix, $r11C_nmatrix);
#SAVE($rfile11CV_nmatrix, $r11C_nmatrix);
#SAVE($rfile11CR_nmatrix, $r11CR_nmatrix);
#SAVE($rfile11CN_nmatrix, $r11CN_nmatrix);
#SAVE($rfile11CL_nmatrix, $r11CL_nmatrix);
#SAVE($rfile11CZ_nmatrix, $r11CZ_nmatrix);
#SAVE($rfile11CRL_nmatrix, $r11CRL_nmatrix);
#SAVE($rfile11CRZ_nmatrix, $r11CRZ_nmatrix);
#-------------------------------
#SAVE($rfile11G_nmatrix, $r11G_nmatrix);
#SAVE($rfile11GV_nmatrix, $r11G_nmatrix);
#SAVE($rfile11GR_nmatrix, $r11GR_nmatrix);
#SAVE($rfile11GN_nmatrix, $r11GN_nmatrix);
#SAVE($rfile11GL_nmatrix, $r11GL_nmatrix);
#SAVE($rfile11GZ_nmatrix, $r11GZ_nmatrix);
#SAVE($rfile11GRL_nmatrix, $r11GRL_nmatrix);
#SAVE($rfile11GRZ_nmatrix, $r11GRZ_nmatrix);
#-------------------------------
#SAVE($rfile15_nmatrix, $r15_nmatrix);
#SAVE($rfile17_nmatrix, $r17_nmatrix);
#SAVE($rfile19_nmatrix, $r19_nmatrix);
#SAVE($rginfo, $geneinfo);
SAVE($rginfo, $r_hist);
SAVE($rgntable, $r_gntable);

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
		my $rfile101 = "$dir/gene_disrupt_transloc_q-".$nsample.".csv";
		my $rfile101M = "$dir/gene_disrupt_allgenotyped_transloc_q-".$nsample.".csv";
		my $rfile101B = "$dir/gene_indel_transloc_q-".$nsample.".csv";
		my $rfile102 = "$dir/raw_transloc_q-".$nsample.".csv";
		my $rfile103 = "$dir/protein_transloc_q-".$nsample.".csv";
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
			my $rfile201 = "$dir/gene_disrupt_transloc-diffchr_q-".$nsample.".csv";
			my $rfile201M = "$dir/gene_disrupt_allgenotyped_transloc-diffchr_q-".$nsample.".csv";
			my $rfile201B = "$dir/gene_indel_transloc-diffchr_q-".$nsample.".csv";
			my $rfile202 = "$dir/raw_transloc-diffchr_q-".$nsample.".csv";
			my $rfile203 = "$dir/protein_transloc-diffchr_q-".$nsample.".csv";
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
			my $rfile301 = "$dir/gene_disrupt_transloc-samechr_q-".$nsample.".csv";
			my $rfile301M = "$dir/gene_disrupt_allgenotyped_transloc-samechr_q-".$nsample.".csv";
			my $rfile301B = "$dir/gene_indel_transloc-samechr_q-".$nsample.".csv";
			my $rfile302 = "$dir/raw_transloc-samechr_q-".$nsample.".csv";
			my $rfile303 = "$dir/protein_transloc-samechr_q-".$nsample.".csv";
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
	R_PCAplot($data_path, $rfile1V_nmatrix, "native", $nsample);
	R_PCAplot($data_path, $rfile1R_nmatrix, "relative", $nsample);
#	R_PCAplot($data_path, $rfile1N_nmatrix, "averaged", $nsample);
	R_PCAplot($data_path, $rfile1L_nmatrix, "log2-transformed[native]", $nsample);
#	R_PCAplot($data_path, $rfile1Z_nmatrix, "zscore", $nsample);
#	R_PCAplot($data_path, $rfile1RL_nmatrix, "log2-transformed[relative]", $nsample);
#	R_PCAplot($data_path, $rfile1RZ_nmatrix, "zscore-relative", $nsample);
	R_PCAplot($data_path, $rfile1AV_nmatrix, "native", $nsample);
	R_PCAplot($data_path, $rfile1AR_nmatrix, "relative", $nsample);
	R_PCAplot($data_path, $rfile1AL_nmatrix, "log2-transformed[native]", $nsample);

	R_PCAplot($data_path, $rfile7V_nmatrix, "native", $nsample);
	R_PCAplot($data_path, $rfile7R_nmatrix, "relative", $nsample);
#	R_PCAplot($data_path, $rfile7N_nmatrix, "averaged", $nsample);
	R_PCAplot($data_path, $rfile7L_nmatrix, "log2-transformed[native]", $nsample);
#	R_PCAplot($data_path, $rfile7Z_nmatrix, "zscore", $nsample);
#	R_PCAplot($data_path, $rfile7RL_nmatrix, "log2-transformed[relative]", $nsample);
#	R_PCAplot($data_path, $rfile7RZ_nmatrix, "zscore-relative", $nsample);

	R_PCAplot($data_path, $rfile9V_nmatrix, "native", $nsample);
	R_PCAplot($data_path, $rfile9R_nmatrix, "relative", $nsample);
#	R_PCAplot($data_path, $rfile9N_nmatrix, "averaged", $nsample);
	R_PCAplot($data_path, $rfile9L_nmatrix, "log2-transformed[native]", $nsample);
#	R_PCAplot($data_path, $rfile9Z_nmatrix, "zscore", $nsample);
#	R_PCAplot($data_path, $rfile9RL_nmatrix, "log2-transformed[relative]", $nsample);
#	R_PCAplot($data_path, $rfile9RZ_nmatrix, "zscore-relative", $nsample);

#	R_PCAplot($data_path, $rfile11V_nmatrix, "native", $nsample);
#	R_PCAplot($data_path, $rfile11R_nmatrix, "relative", $nsample);
#	R_PCAplot($data_path, $rfile11N_nmatrix, "averaged", $nsample);
#	R_PCAplot($data_path, $rfile11L_nmatrix, "log2-transformed[native]", $nsample);
#	R_PCAplot($data_path, $rfile11Z_nmatrix, "zscore", $nsample);
#	R_PCAplot($data_path, $rfile11RL_nmatrix, "log2-transformed[relative]", $nsample);
#	R_PCAplot($data_path, $rfile11RZ_nmatrix, "zscore-relative", $nsample);

#	R_PCAplot($data_path, $rfile11CV_nmatrix, "native", $nsample);
#	R_PCAplot($data_path, $rfile11CR_nmatrix, "relative", $nsample);
#	R_PCAplot($data_path, $rfile11CN_nmatrix, "averaged", $nsample);
#	R_PCAplot($data_path, $rfile11CL_nmatrix, "log2-transformed[native]", $nsample);
#	R_PCAplot($data_path, $rfile11CZ_nmatrix, "zscore", $nsample);
#	R_PCAplot($data_path, $rfile11CRL_nmatrix, "log2-transformed[relative]", $nsample);
#	R_PCAplot($data_path, $rfile11CRZ_nmatrix, "zscore-relative", $nsample);

#	R_PCAplot($data_path, $rfile11GV_nmatrix, "native", $nsample);
#	R_PCAplot($data_path, $rfile11GR_nmatrix, "relative", $nsample);
#	R_PCAplot($data_path, $rfile11GN_nmatrix, "averaged", $nsample);
#	R_PCAplot($data_path, $rfile11GL_nmatrix, "log2-transformed[native]", $nsample);
#	R_PCAplot($data_path, $rfile11GZ_nmatrix, "zscore", $nsample);
#	R_PCAplot($data_path, $rfile11GRL_nmatrix, "log2-transformed[relative]", $nsample);
#	R_PCAplot($data_path, $rfile11GRZ_nmatrix, "zscore-relative", $nsample);
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
sub Read_GOdata{
my $file = shift;

my $transposons =<<"EOS";
GO:0015074
IPR043128
IPR000477
IPR001584
IPR012337
IPR036397
IPR041588
IPR041577
IPR005162
IPR041373
IPR019557
IPR013103
IPR002156
IPR004242
IPR029480
EOS

my @TRs = split(/\n/, $transposons);
my $hTRs = {};
foreach my $str (@TRs){
	$hTRs->{$str} = 1;
}

open(my $fh, "<", $file) or die;
my $hash = {};
my $cnt = 0;
my $cnt_tr = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	$line =~ s/\"//g;
	
	if($cnt == 0){
		$cnt++;
		next;
	}
	
#	ID,GO+IPR,IPR,Biological Process,Molecular Function,Cellular Component,IPR ID
#	Glyma.01G000100.v4.JE1,Glyma.01G000100.v4.JE1.t1,false,false,-,-,-,-
#	Glyma.01G000100.v4.JE1,Glyma.01G000100.v4.JE1.t2,false,false,-,-,-,-
#	Glyma.01G001000.v4.JE1,Glyma.01G001000.v4.JE1.t1,true,true,-,GO:0051879;GO:0051087;GO:0001671,-,IPR039981;IPR015310;IPR036338
	
	my @A = split(/,/, $line);
	my $numA = @A;
	for(my $i = 4; $i < $numA; $i++){
		if($A[$i] ne '-'){
			my @B = split(/\;/, $A[$i]);
			foreach my $str (@B){
				$str =~ s/\s//g;
				if($hTRs->{$str}){
					$hash->{transposon}{$A[0]} = "true";
					$cnt_tr++;
					last;
				}
			}
		}
	}
	$cnt++;
}
close $fh;

if($cnt_tr > 0){
	$hash->{status} = "true";
}

return $hash;
}


#-------------------------------------------------------------------------------
sub R_PCAplot{
my $data_path = shift;
my $file0 = shift;
my $description = shift;
my $num_sample = shift;

my $R_file0 = $file0;
$R_file0 =~ s/\.csv/_0\.txt/;

my $img_file0 = $file0;
my $img_file1 = $file0;
my $img_file2 = $file0;
my $img_file3 = $file0;
my $img_file4 = $file0;
$img_file0 =~ s/forClust_/plotClust_/;
$img_file1 =~ s/forClust_/plotClust_/;
$img_file2 =~ s/forClust_/plotClust_/;
$img_file3 =~ s/forClust_/plotClust_/;
$img_file4 =~ s/forClust_/plotClust_/;
$img_file0 =~ s/\.csv/_PCAplot\.png/;
$img_file1 =~ s/\.csv/_flushclust\.png/;
$img_file2 =~ s/\.csv/_flushclust_hang\.png/;
$img_file3 =~ s/\.csv/_hclust\.png/;
$img_file4 =~ s/\.csv/_hclust_hang\.png/;

my $spcols = "5";
if($num_sample > 1){
	$num_sample += 4;
	$spcols = "5:".$num_sample;
}

my $script0 =<<"EOS";
library(stats)
workingDir = "$data_path"
setwd(workingDir)
EOS

$script0 .=<<"EOS";
#----------------------flashclust----------------------#
library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE)
Data <- read.csv("$file0", header=T, row.names=1)
Data <- Data[, $spcols ]
datExpr <- as.data.frame(t(Data))
sampleTree <- flashClust(dist(datExpr), method = "complete")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
png("$img_file1", width=3000, height=1200)
plot(sampleTree, main = "Sample clustering, $description (flashClust)", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
EOS

my $disabled =<<"EOS";
png("$img_file2", width=3000, height=1200)
plot(sampleTree, hang=-1, main = "Sample clustering, $description (flashClust)", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
EOS

open(my $rfh0, ">", $R_file0) or die;
print $rfh0 $script0;
close $rfh0;

print "! generating clustering based on [$file0] ...\n";
if(system("Rscript --vanilla --slave $R_file0 > /dev/null 2>&1") == 0){
	if(-e $img_file1 && -e $img_file2){
		print "! output [$img_file1]\n";
		print "! output [$img_file2]\n";
	}
#	else{
#		print "! executed but missing image files...\n";
#	}
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
	if(0 <= $val && $val < 0.6){
		$hash->{C0} += 1;
	}
	elsif(0.6 <= $val && $val < 1.2){
		$hash->{C1} += 1;
	}
	elsif(1.2 <= $val && $val < 1.8){
		$hash->{C2} += 1;
	}
	elsif(1.8 <= $val){
		$hash->{C3} += 1;
	}
}

my @Hist;
my @Header;
push(@Header, "0 - 0.6");
push(@Header, "0.6 - 1.2");
push(@Header, "1.2 - 1.8");
push(@Header, "1.8 -");

if(defined $hash->{C0}){
	push(@Hist, $hash->{C0});
}
else{
	push(@Hist, 0);
}
if(defined $hash->{C1}){
	push(@Hist, $hash->{C1});
}
else{
	push(@Hist, 0);
}
if(defined $hash->{C2}){
	push(@Hist, $hash->{C2});
}
else{
	push(@Hist, 0);
}
if(defined $hash->{C3}){
	push(@Hist, $hash->{C3});
}
else{
	push(@Hist, 0);
}

my $hist = join(",", @Hist);
my $header = join(",", @Header);

return ($hist, $header);
}


#-------------------------------------------------------------------------------
sub Array2histcsv3{
my $A = shift;

my $hash = {};
foreach my $val (@{$A}){
	if(defined $val){
		$hash->{$val} += 1;
	}
}

my @NR = keys(%{$hash});
@NR = sort {$a cmp $b} @NR;

my @NR2;
foreach my $val (@NR){
	if($val eq 'Ref'){
		push(@NR2, $val);
	}
}
foreach my $val (@NR){
	if($val ne 'Ref'){
		push(@NR2, $val);
	}
}

my $AoA = [];
my $ngeno = 0;
my $ngeno_sample = 0;
my @items;
foreach my $val (@NR2){
	if($val ne '-'){
		my $str = $val."=".$hash->{$val};
		push(@items, $str);
		
		my @tmp;
		push(@tmp, $val);
		push(@tmp, $hash->{$val});
		push(@{$AoA}, \@tmp);
		$ngeno++;
		$ngeno_sample += $hash->{$val};
	}
	else{
		my $str = "NA=".$hash->{$val};
		push(@items, $str);
	}
}

@{$AoA} = sort {$a->[1] <=> $b->[1]} @{$AoA};
my $nminor_sample = $AoA->[0][1];
my $nminor_ratio = $nminor_sample / $ngeno_sample;
if($ngeno == 1){
	$nminor_sample = "-";
	$nminor_ratio = "-";
}

#my $header = "num_genotype,num_minor,ratio_minor,info";
my $gnstat = join("\;", @items);
#my $r = $ngeno.",".$nminor_sample.",".$nminor_ratio.",".$gnstat;
my $r = $ngeno.",".$nminor_ratio.",".$gnstat;

return $r;
}


#-------------------------------------------------------------------------------
sub Select_lessConverved{
my $lines = shift;
my $hgeneclass = shift;
my $str = shift;
my $hgoipr = shift;
my $svscore_nsample_th = shift;

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
		if(defined $hgoipr->{status}){
			if($hgoipr->{transposon}{$A[0]} && $hgoipr->{transposon}{$A[0]} eq 'true'){
				next;
			}
			else{
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
		}
		else{
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
	}
	$cnt++;
}

@{$AoA} = sort {$a->[0] <=> $b->[0]} @{$AoA};
my $num_gid = @{$AoA};

#for(my $j = 0; $j < 500; $j++){
for(my $j = 0; $j < $num_gid; $j++){
	if($AoA->[$j][1]){
		my $sum_score = $AoA->[$j][0];
		my $A = $AoA->[$j][1];
		my $gid = $AoA->[$j][2];
		
		if($sum_score <= $svscore_nsample_th){
			$sorted .= join(",", @{$A})."\n";
			if($str ne 'null'){
				$hgeneclass->{$str}{$gid} = 1;
			}
		}
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
								$hash->{neighbor_tlen}{$neighbor_tlen} += 1;
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
			my $gntype_promoter = "-";
			my $gntype_utr = "-";
			
			unless($A[32]){
				next;
			}
			if($A[32] !~ /true/i){			# skip non-coding transcript
				$hash->{geno_including_flanking}{$id}{$gid} = "false";
				next;
			}
			
			my $disrupt_promoter = 0;
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
				$disrupt_promoter = $promoter;
				
				if($promoter > 1){
					$disrupt_promoter = 1 / $promoter;
				}
				
				if($disrupt_promoter >= 0.98){
					$gntype_promoter = "Ref";
				}
				if($disrupt_promoter >= 0.98){
					$gntype_promoter = "Ref";
				}
				elsif($disrupt_promoter >= 0.90){
					if($promoter > 1){
						$gntype_promoter = "A1";
					}
					else{
						$gntype_promoter = "C1";
					}
				}
				elsif($disrupt_promoter >= 0.80){
					if($promoter > 1){
						$gntype_promoter = "A2";
					}
					else{
						$gntype_promoter = "C2";
					}
				}
				else{
					if($promoter > 1){
						$gntype_promoter = "B";
					}
					else{
						$gntype_promoter = "D";
					}
				}
			}
			
			my $disrupt_utr = 0;
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
				$disrupt_utr = $utr3;
				
				if($utr3 > 1){
					$disrupt_utr = 1 / $utr3;
				}
				
				if($disrupt_utr >= 0.98){
					$gntype_utr = "Ref";
				}
				elsif($disrupt_utr >= 0.90){
					if($utr3 > 1){
						$gntype_utr = "A1";
					}
					else{
						$gntype_utr = "C1";
					}
				}
				elsif($disrupt_utr >= 0.80){
					if($utr3 > 1){
						$gntype_utr = "A2";
					}
					else{
						$gntype_utr = "C2";
					}
				}
				else{
					if($utr3 > 1){
						$gntype_utr = "B";
					}
					else{
						$gntype_utr = "D";
					}
				}
			}
			
	#		my $sv2 = $sv;					# value representing deletion and insertion separately: [0] < deletion < [1.0] < insertion < [2.0]
			my $dbspan_woN = $A[14];
			my $hitaln = $A[15];
			my $hitspan_woN = $A[17];
			my $bpinsert = $A[18];
			
			my $gntype = "-";
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
						
			#			if($sv2 > 10){
			#				$sv2 = 10;
			#			}
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
				
				if($sv >= 0.98){
					$gntype = "Ref";
				}
				elsif($sv >= 0.90){
					if($sv2 > 1){
						$gntype = "A1";
					}
					else{
						$gntype = "C1";
					}
				}
				elsif($sv >= 0.80){
					if($sv2 > 1){
						$gntype = "A2";
					}
					else{
						$gntype = "C2";
					}
				}
				else{
					if($sv2 > 1){
						$gntype = "B";
					}
					else{
						$gntype = "D";
					}
				}
			}
			
			if($refgenome ne 'null' && $refgenome eq $id){
				$hash->{refchr}{$id}{$gid} = $A[2];
				$hash->{refpos0}{$id}{$gid} = $A[5];
				$hash->{refpos1}{$id}{$gid} = $A[6];
			}
			
			$hash->{gntype}{$id}{$gid} = $gntype;
			$hash->{gntype_promoter}{$id}{$gid} = $gntype_promoter;
			$hash->{gntype_utr}{$id}{$gid} = $gntype_utr;
			$hash->{sv}{$id}{$gid} = $sv;
			$hash->{sv2}{$id}{$gid} = $sv2;
			$hash->{svprom}{$id}{$gid} = $promoter;
			$hash->{svutr3}{$id}{$gid} = $utr3;
			$hash->{svprom_disrupt}{$id}{$gid} = $disrupt_promoter;
			$hash->{svutr3_disrupt}{$id}{$gid} = $disrupt_utr;
#			$hash->{svprot}{$id}{$gid} = $sv;
			
			$hash->{svjudge1}{$id}{$gid} = $A[22];					# native judge
			$hash->{svjudge2}{$id}{$gid} = "null";
#			if($A[22] !~ /missing border/ && $A[22] !~ /absent/ && $A[22] !~ /NNN/ && $A[22] !~ /no hit but some/ && $A[22] !~ /but not n-BDBH/){
#				if($A[22] =~ /present/ || $A[22] =~ /insertion/){
#					$hash->{svjudge2}{$id}{$gid} = "genotyped";		# judge classified
#				}
#			}
			if($A[22] !~ /NNN/ && $A[22] !~ /but not n-BDBH/){
				$hash->{svjudge2}{$id}{$gid} = "genotyped";		# judge classified
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
					my $diff_protgap = abs(100 - $A[37]);		# A[37] = rate of sequence range that covers DB protein
					
					my $dbprot_factor = log($protlen_db) / log(10);
					my $difflen_ratio = abs($protlen_db - $protlen_q) / $protlen_db;
					
					$rprot = ($pcnt_similarity - $diff_protgap) / 100 - ($difflen_ratio ** $dbprot_factor);
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
sub Stats_fasta{
my $file = shift;

print "! reading refseq [$file]...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
my $cnt = 0;
my $ID;
my $seq;
my $total_len = 0;
my @SID;
my @SID_ori;
my $norder = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	my @A = split(/\t/, $line);
	if($line =~ /\>/){
		unless($cnt == 0){
			if($hash->{$ID}){
				print "! duplicated seq ID : [$ID]\n";
			}
			
#			if($ID !~ /scaffold/ && $ID !~ /mitoch/ && $ID !~ /chloro/){
				my $len = length($seq);
	#			my $ID_ori = $ID;
	#			$ID =~ s/chr0//i;
	#			$ID =~ s/chr//i;
	#			$ID =~ s/Gm0//i;
	#			$ID =~ s/Gm//i;
	#			$hash->{$ID}{originalID} = $ID_ori;
				$hash->{$ID}{order} = $norder;
	#			$hash->{$ID}{seq} = $seq;
				$hash->{$ID}{len} = $len;
				$hash->{$ID}{cumlen0} = $total_len;
	#			$hash->{$ID_ori}{cumlen0} = $total_len;
				$total_len += $len;
				$hash->{$ID}{cumlen1} = $total_len;
	#			$hash->{$ID_ori}{cumlen1} = $total_len;
				push(@SID, $ID);
	#			push(@SID_ori, $ID_ori);
				$norder++;
				print " [$ID] = [$len] bp\n";
#			}
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
	
#	if($ID !~ /scaffold/ && $ID !~ /mitoch/ && $ID !~ /chloro/){
		my $len = length($seq);
#		my $ID_ori = $ID;
#		$ID =~ s/chr0//i;
#		$ID =~ s/chr//i;
#		$ID =~ s/Gm0//i;
#		$ID =~ s/Gm//i;
#		$hash->{$ID}{originalID} = $ID_ori;
		$hash->{$ID}{order} = $norder;
#		$hash->{$ID}{seq} = $seq;
		$hash->{$ID}{len} = $len;
		$hash->{$ID}{cumlen0} = $total_len;
#		$hash->{$ID_ori}{cumlen0} = $total_len;
		$total_len += $len;
		$hash->{$ID}{cumlen1} = $total_len;
#		$hash->{$ID_ori}{cumlen1} = $total_len;
		push(@SID, $ID);
#		push(@SID_ori, $ID_ori);
		$norder++;
		print " [$ID] = [$len] bp\n";
#	}
}

my $numID = @SID;
$hash->{SID} = \@SID;
#$hash->{SID_ori} = \@SID_ori;
$hash->{total_len} = $total_len;

print "! total [$numID] sequence ID, [$total_len] bp\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Gff_to_hash{
my $gff3 = shift;

print "! reading [$gff3]...\n";
open(my $fh, "<", $gff3) or die;
my $hash = {};
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
	
	if($A[0] =~ /scaffold/ || $A[0] =~ /unanchored/ || $A[0] =~ /mitochon/ || $A[0] =~ /chloro/){
#	if($A[0] =~ /mitochon/ || $A[0] =~ /chloro/){
		next;
	}
	
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
		
		if($gid){
			$hash->{hgff}{$gid}{Chr} = $A[0];		#1	seqid
			$hash->{hgff}{$gid}{pos0} = $A[3];		#2	pos0
			$hash->{hgff}{$gid}{pos1} = $A[4];		#3	pos1
			$hash->{hgff}{$gid}{strand} = $A[6];	#4	strand
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
}
close $fh;

#print "! [$gcnt] genes, [$tcnt] transcripts, [$pcnt] CDS lines\n";
print "! [$gcnt] genes\n";

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
else{
	print "! skip [$file] due to missing data lines\n";
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





