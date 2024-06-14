#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use FindBin;
use Getopt::Long;
#use File::HomeDir;
use Sys::Hostname;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "sv2assoc.pl version 1.01\n";
$version .= "last update: [2023\/3\/21]\n";
$version .= "copyright: ryoichi yano\n";

#print $version;

#-------------------------------------------------------------------------------

#---------------------------------------------------------//
my $example_csv = "example_phenodata.csv";
my $help_command =<<"EOS";
Basic usage: 
  /path/to/sv2assoc.pl -f [fasta] -d [sv table] -i [phenotype data] -o [Result_directory]

--fasta or -f         : reference fasta
--data or -d          : sv data
--i or -i             : phenotype data (csv)
--output_dir or -o    : output directory (automatically set if not specified)
--help or -h          : display usage

! "$example_csv" will be produced if not specified.

EOS

#---------------------------------------------------------//
my $host = hostname();
my $script_path = $FindBin::Bin;
my $wpath = getcwd();
#my $HOME = File::HomeDir->my_home;

#---------------------------------------------------------//
#gettin parameters from command line
my $dbfasta;
my $dbgff;
my $svtable;
my $phenodata;
my $dir;
my $gn_th = 50;
my $help;

GetOptions('--fasta=s' => \$dbfasta, '--gff=s' => \$dbgff, '--data=s' => \$svtable, '--i=s' => \$phenodata, '--o=s' => \$dir, '--help' => \$help);

my $status = "OK";
if(! $dbfasta || ! -e $dbfasta){
	print "! missing fasta\n";
	$status = "missing";
}
#if(! $dbgff || ! -e $dbgff){
#	print "! missing Gff\n";
#	$status = "missing";
#}
if(! $svtable || ! -e $svtable){
	print "! missing sv_genic\n";
	$status = "missing";
}
if(! $phenodata || ! -e $phenodata){
	print "! missing phenodata\n";
	$status = "missing";

	my $expheno =<<"EOS";
ID,pheno1,pheno2
data_type,case,numeric
group_low,-,-
group_high,-,-
sample_1,1,0.894393605
sample_2,1,0.172733425
sample_3,1,1.44403025
sample_4,1,4.2246643
sample_5,2,8.49038407
sample_6,2,7.43537367
EOS

	SAVE($example_csv, $expheno);
}
if($status eq 'missing'){
	print "\n$help_command";
	goto END;
}

my $prefix = $phenodata;
$prefix =~ s/\.csv//;
unless($dir){
	$dir = "sv2assoc_".$prefix;
}

my $log = "\n";
$log .= "! fasta        : [$dbfasta]\n";
#$log .= "! Gff          : [$dbgff]\n";
$log .= "! sv_genic     : [$svtable]\n";
$log .= "! phenodata    : [$phenodata]\n";
$log .= "! outdir       : [$dir]\n";

print $log,"\n";

my @APP;
push(@APP, "plink");

my $hbin = {};
foreach my $app (@APP){
	$hbin->{$app}[0] = $app;
}
my $hpath = {};
my $hbinpath = Search_app2($hbin, \@APP, $hpath, $script_path, "pipe", $phenodata);
if($hbinpath->{err}){
	print "! abort script due to missing program(s)..., please confirm program is in PATH.\n";
	goto END;
}

my $hpheno = Read_phenodata3($phenodata);
if(! $hpheno || $hpheno->{err}){
	goto END;
}

#---------------------------------------------------------------//
unless(-e $dir){
	system("mkdir $dir");
}
chdir $dir;
$wpath .= "/".$dir;

$dbfasta = LnFile($dbfasta);
$svtable = LnFile($svtable);
$phenodata = LnFile($phenodata);

for(my $i = 0; $i < $hpheno->{npheno}; $i++){
	if($hpheno->{pheno_datatype}{$i} eq 'null'){
		print "! missing data_type for [$hpheno->{pheno_name}{$i}] (numeric or case), skip...\n";
	}
	
	print " - phenotype = [$hpheno->{pheno_name}{$i}] ($hpheno->{pheno_datatype}{$i})\n";
	my $each_pheno_file = "phenovalue_".$hpheno->{pheno_name}{$i}.".phe";
	Prep_phenodata($hpheno, $i, $hpheno->{pheno_name}{$i}, $each_pheno_file);
}

my $nsample = $hpheno->{nsample};
my $num_pheno = $hpheno->{npheno};
my $tmpdir_Rplots = "_tmpdir_Rplots";
unless(-e $tmpdir_Rplots){
	system("mkdir $tmpdir_Rplots");
}

my ($PlotColors, $num_color) = Prep_plotColors(2);
#my $color_abline = "#708090";		# slategray
my $color_abline = "#808080";		# gray

my $rh = Read_gntable($svtable, $hpheno);
my $err = $rh->{err};
my $SID = $rh->{SID};
my $hash = $rh->{hash};
my $list_nosex = $rh->{list_nosex};
my $list_fam = $rh->{list_fam};
my $list_map = $rh->{list_map};

my @Pat;
push(@Pat, "pat1");
push(@Pat, "pat2");
push(@Pat, "pat3");

print "! preparing files...\n";
foreach my $pat (@Pat){
	my $upref = "gn_".$pat;
	my $rfile0 = $upref.".nosex";
	my $rfile1 = $upref.".fam";
	my $rfile2 = $upref.".map";
	my $rfile3 = $upref.".ped";
	my $log_plink = $upref.".log";
	
	my $combined_ped;
	foreach my $sid (@{$SID}){
		if($hash->{$sid}{$pat}){
			$combined_ped .= $hash->{$sid}{$pat}."\n";
		}
		else{
			print "! error : missing genotype data for [$sid]\n";
			goto END;
		}
	}
	SAVE($rfile3, $combined_ped);
}
print "! done.\n";

my $hfasta = Stats_fasta($dbfasta);

foreach my $pat (@Pat){
	print "\n============================================================\n";
	print "! [$pat]\n";
	my $upref = "gn_".$pat;
	my $rfile0 = $upref.".nosex";
	my $rfile1 = $upref.".fam";
	my $rfile2 = $upref.".map";
	my $rfile3 = $upref.".ped";
	my $log_plink = $upref.".log";
	
	SAVE($rfile0, $list_nosex);
	SAVE($rfile1, $list_fam);
	SAVE($rfile2, $list_map);
	
	my $cmd3 = "$hbinpath->{plink} --file $upref --make-bed --out $upref --allow-extra-chr";
	print "! cmd=[$cmd3]\n";
	if(system("$cmd3 > $log_plink 2>&1") != 0){
		print "! make-bed failed\n";
		next;
	}
	
	SAVE($rfile0, $list_nosex);
	SAVE($rfile1, $list_fam);
	SAVE($rfile2, $list_map);
	
	for(my $i = 0; $i < $hpheno->{npheno}; $i++){
		my $each_pheno_file = "phenovalue_".$hpheno->{pheno_name}{$i}.".phe";
		my $pref_assoc = "assoc_".$pat."_".$hpheno->{pheno_name}{$i};
		
		my $topx = 100;
		my $each_assoc = $pref_assoc.".qassoc";
		my $each_adjust = $pref_assoc.".qassoc.adjusted";
		my $each_summary = "summary_".$pref_assoc."_plink.csv";
		my $each_phenovals_top100 = "top100_phenovals_".$pref_assoc."_plink.csv";		# data for box plots
		my $each_topx_summary = "Top".$topx."_".$pref_assoc."_plink.csv";
		my $each_topx_gnmatrix1 = "Top".$topx."_".$pref_assoc."_genotype1.tsv";
		my $each_topx_gnmatrix2 = "Top".$topx."_".$pref_assoc."_genotype2.tsv";
		
		if($hpheno->{pheno_datatype}{$i} eq 'case'){
			$each_assoc = $pref_assoc.".assoc";
			$each_adjust = $pref_assoc.".assoc.adjusted";
		}
		
		my $hfplots = {};
		$hfplots->{each_plot1} = "plot_plink_".$pref_assoc."_-Log10P.png";
		$hfplots->{each_plot1_posi} = "plot_plink_".$pref_assoc."_-Log10P_posi.png";
		$hfplots->{each_plot1_nega} = "plot_plink_".$pref_assoc."_-Log10P_nega.png";
		$hfplots->{each_plot2} = "plot_plink_".$pref_assoc."_-Log10P_direction.png";
		
#		if($hpheno->{pheno_datatype}{$i} eq 'numeric'){
			if(-e $each_assoc && -e $each_adjust){
#				print "! [$each_assoc] [$each_adjust] already exist, skip plink...\n";
				system("rm $each_assoc");
				system("rm $each_adjust");
			}
#			else{
				my $sh = "cmd_assoc_".$hpheno->{pheno_name}{$i}.".sh";
				my $cmd = "$hbinpath->{plink} --bfile $upref --pheno $each_pheno_file --missing-code -9,NA,na --assoc --adjust --out $pref_assoc --allow-extra-chr";
				SAVE($sh, $cmd);
				print "! cmd=[$cmd]\n";
				if(system("bash $sh >> $log_plink 2>&1") != 0){
					print "! process failed, skip...\n";
					next;
				}
#			}
#		}
		
		my ($hdata, $htopx) = Summarize_qassoc($each_assoc, $each_adjust, $each_summary, $each_phenovals_top100, $hfasta, $hpheno->{nsample}, $gn_th, $topx, $hpheno->{pheno_datatype}{$i});
		if($hdata->{err}){
			goto JUMP;
		}
		
		R_GenPlot($hfasta, $hdata, $PlotColors, $num_color, $color_abline, $tmpdir_Rplots, $hfplots, 2400, 800, $pref_assoc, $wpath);
#		Extract_topX_candidate($htopx, $each_pheno_file, $tsv_gn1, $genebased_aa_csv, $topx, $each_topx_summary, $each_topx_gnmatrix1, $each_topx_gnmatrix2, $hpheno->{pheno_name}{$i});
	}
}

END:{
	print "\n! End of script.\n\n";
	my $end = 1;
}


################################################################################
#-------------------------------------------------------------------------------
sub Read_gntable{
my ($file, $hpheno) = @_;

my $hID2n = $hpheno->{ID2n};
my @SID = keys(%{$hID2n});
@SID = sort {$a cmp $b} @SID;

print "! reading [$file]...\n";
unless($file){
	return;
}

# gene_id	Chr	pos0	pos1	region	num_genotype	ratio_minor

open(my $fh, "<", $file) or return;
my $hash = {};
my $hID2i = {};
my $hi2ID = {};
my $cnt = 0;
my $err = 0;
my $list_nosex;
my $list_fam;
my $list_map;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	unless($line){
		next;
	}
	
	my @A = split(/,/, $line);
	if($cnt == 0){
		my $i = 0;
		foreach my $val (@A){
			if(defined $val){
				if($val =~ /\[/){
					my @tmp = split(/\[/, $val);
					$val = $tmp[0];
					$val =~ s/\s//;
				}
				
				if(defined $hID2n->{$val}){
					$hID2i->{$val} = $i;
					$hi2ID->{$i} = $val;
				}
			}
			$i++;
		}
		
		my $j = 0;
		foreach my $sid (@SID){
			if(defined $hID2i->{$sid}){
#				print " [$sid] = $hID2i->{$sid} <OK>\n";
				$list_nosex .= "0\t$sid\n";
				$list_fam .= "$sid 1 0 0 1 1\n";
				$hash->{$sid}{pat1} = "0 $sid 0 0 0 -9";
				$hash->{$sid}{pat2} = "0 $sid 0 0 0 -9";
				$hash->{$sid}{pat3} = "0 $sid 0 0 0 -9";
				$j++;
			}
			else{
				print " [$sid] = <missing>\n";
				$err = 1;
			}
		}
		
		if($err == 0){
			print " - [$j] sample name <all OK>\n";
		}
		else{
			last;
		}
	}
	else{
		if(defined $A[6] && $A[6] eq '-'){
			next;
		}
		
		my $seqnum = $A[1];
		$seqnum =~ s/chr0//i;
		$seqnum =~ s/chr//i;
		$list_map .= $seqnum."\t".$seqnum."_".$A[2]."__".$A[0]."__".$A[4]."\t0\t".$A[2]."\n";
		
		my @I = keys(%{$hi2ID});
		@I = sort {$a <=> $b} @I;
		
		foreach my $i (@I){
			my $sid = $hi2ID->{$i};
			if(! defined $A[$i] || $A[$i] eq 'NA' || $A[$i] eq 'na' || $A[$i] eq '0' || $A[$i] eq '-9'){
				$hash->{$sid}{pat1} .= " 0 0";
				$hash->{$sid}{pat2} .= " 0 0";
				$hash->{$sid}{pat3} .= " 0 0";
			}
			elsif($A[$i] eq 'Ref'){
				$hash->{$sid}{pat1} .= " AA AA";
				$hash->{$sid}{pat2} .= " AA AA";
				$hash->{$sid}{pat3} .= " AA AA";
			}
			elsif($A[$i] =~ /C/){
				if($A[$i] eq 'C1'){
					$hash->{$sid}{pat1} .= " C C";
				}
				elsif($A[$i] eq 'C2'){
					$hash->{$sid}{pat1} .= " G G";
				}
				else{
					$hash->{$sid}{pat1} .= " C C";
				}
				
				$hash->{$sid}{pat2} .= " C C";
				$hash->{$sid}{pat3} .= " AA AA";
			}
			elsif($A[$i] =~ /D/){
				$hash->{$sid}{pat1} .= " T T";
				$hash->{$sid}{pat2} .= " T T";
				$hash->{$sid}{pat3} .= " T T";
			}
			elsif($A[$i] =~ /A/){
				if($A[$i] eq 'A1'){
					$hash->{$sid}{pat1} .= " AAC AAC";
				}
				elsif($A[$i] eq 'A2'){
					$hash->{$sid}{pat1} .= " AGG AGG";
				}
				else{
					$hash->{$sid}{pat1} .= " AAC AAC";
				}
				
				$hash->{$sid}{pat2} .= " AAC AAC";
				$hash->{$sid}{pat3} .= " AA AA";
			}
			elsif($A[$i] =~ /B/){
				$hash->{$sid}{pat1} .= " AAT AAT";
				$hash->{$sid}{pat2} .= " AAT AAT";
				$hash->{$sid}{pat3} .= " AAT AAT";
			}
		}
	}
	$cnt++;
}
close $fh;

print "! total [$cnt] lines\n";

my $rh = {};
$rh->{err} = $err;
$rh->{SID} = \@SID;
$rh->{hash} = $hash;
$rh->{list_nosex} = $list_nosex;
$rh->{list_fam} = $list_fam;
$rh->{list_map} = $list_map;

return $rh;
}


#-------------------------------------------------------------------------------
sub Summarize_qassoc{
my ($each_assoc, $each_adjust, $each_summary, $each_phenovals_top100, $hfasta, $nsample, $gn_th, $topx, $data_type) = @_;

my $hash = {};
if(-e $each_adjust){
	print "! reading [$each_adjust]...\n";
	open(my $fh, "<", $each_adjust) or die;
	my $cnt = 0;
	my $cnt_var = 0;
	my $cnt_dup = 0;
	while(my $l = <$fh>){
		$l =~ s/\n//;
		$l =~ s/\r//;
		
		# (numeric)
#	 CHR                             SNP      UNADJ         GC       BONF       HOLM   SIDAK_SS   SIDAK_SD     FDR_BH     FDR_BY
#	   7                      7_25839789  5.843e-32  7.862e-32  2.028e-26  2.028e-26        INF        INF  2.028e-26  2.705e-25 
#	   7                        7_713475  6.436e-29  2.019e-08  2.234e-23  2.234e-23        INF        INF  7.701e-24  1.027e-22 
#	   6                       6_4255459  6.656e-29  1.207e-15   2.31e-23   2.31e-23        INF        INF  7.701e-24  1.027e-22 
#	   6                        6_762166  1.768e-28  2.631e-08  6.136e-23  6.136e-23        INF        INF  1.534e-23  2.045e-22 
		
		# (case)
#	 CHR                             SNP      UNADJ         GC       BONF       HOLM   SIDAK_SS   SIDAK_SD     FDR_BH     FDR_BY
#	   5                      5_42960102  1.829e-55    0.01553   9.19e-49   9.19e-49        INF        INF   9.19e-49  1.471e-47 
#	   5                      5_42960675  4.858e-55    0.01595  2.441e-48  2.441e-48        INF        INF  1.221e-48  1.954e-47 
#	   5                      5_42958943  1.021e-54    0.01627  5.131e-48  5.131e-48        INF        INF  1.603e-48  2.566e-47
		
		unless($l){
			next;
		}
		if($cnt > 0){
			my @A = split(/\s/, $l);
			my @B;
			foreach my $val (@A){
				if($val){
					push(@B, $val);
				}
			}
			my $numB = @B;
			
			if($numB == 10){
				unless($hash->{$B[1]}){
					$hash->{$B[1]}{UNADJ} = $B[2];
					$hash->{$B[1]}{GC} = $B[3];
					$hash->{$B[1]}{BONF} = $B[4];
					$hash->{$B[1]}{HOLM} = $B[5];
					$hash->{$B[1]}{SIDAK_SS} = $B[6];
					$hash->{$B[1]}{SIDAK_SD} = $B[7];
					$hash->{$B[1]}{FDR_BH} = $B[8];
					$hash->{$B[1]}{FDR_BY} = $B[9];
					$cnt_var++;
				}
				else{
					$cnt_dup++;
				}
			}
		}
		$cnt++;
	}
	close $fh;
	
	print "! [$cnt_var] variants, [$cnt_dup] duplicated entry\n";
}

my $hdata = {};
my $htopx = {};
$hdata->{err} = 0;
if(-e $each_assoc){
	print "! reading [$each_assoc]...\n";
	open(my $fh, "<", $each_assoc) or die;
	my $header = "CHR,CHR_No,SNP,bp,bp cumlative,NMISS,beta,SE,R2,T,P,minus log10P,unadj,GC,BONF,HOLM,SIDAK_SS,SIDAK_SD,FDR_BH,FDR_BY,gene_id,region\n";
	my $header2 = "CHR,bp,NMISS,beta,SE,R2,T,P,minus log10P,unadj,GC,BONF,HOLM,SIDAK_SS,SIDAK_SD,FDR_BH,FDR_BY\n";
	my $r = $header;
	my $cnt = 0;
	my $cnt_with_stat = 0;
	my $cnt_wo_stat = 0;
	my $AoA = [];
	while(my $l = <$fh>){
		$l =~ s/\n//;
		$l =~ s/\r//;
		
		unless($l){
			next;
		}
		if($cnt > 0){
			my @A = split(/\s/, $l);
			my @B;
			foreach my $val (@A){
				if($val){
					push(@B, $val);
				}
				elsif(defined $val && $val eq '0'){
					push(@B, $val);
				}
			}
			my $numB = @B;
			
			if($numB == 9 && $data_type eq 'numeric'){
				if(! $hfasta->{$B[0]}){
					print "! unable to find refseq information for [$B[0]] (1), stop script...\n";
					print "\n$data_type\n$l\n\n";
					$hdata->{err} = 1;
					close $fh;
					return $hdata;
				}
				
				my $ID_ori = $hfasta->{$B[0]}{originalID};
				my $cumlen = $hfasta->{$B[0]}{cumlen0} + $B[2];
				my $gnr = $B[3] / $nsample * 100;
				
				if($hash->{$B[1]} && $gnr >= $gn_th){
					my $mlog10p = - log($B[8]) / log(10);
					
					my $mlog10p_dr = $mlog10p;
					if($B[4] < 0){
						$mlog10p_dr = - $mlog10p;
					}
					
							# $numB == 9 (numeric)
					#	 CHR                             SNP         BP    NMISS       BETA         SE         R2        T            P 
					#	   1                         1_27761      27761       94     0.2153    0.09695    0.05087    2.221      0.02884 
					#	   1                         1_27781      27781       94    -0.3398     0.9658   0.001344  -0.3518       0.7258 
					#	   1                         1_27841      27841       93     0.1685     0.9712  0.0003306   0.1735       0.8627 
					
					my $gid_info = "-,-";
					if($B[1] =~ /__/){
						my @GID = split(/__/, $B[1]);
						if($GID[1] && $GID[2]){
							if($GID[2] eq 'p'){
								$gid_info = $GID[1].",promoter";
							}
							elsif($GID[2] eq 'g'){
								$gid_info = $GID[1].",gene";
							}
							elsif($GID[2] eq 'u'){
								$gid_info = $GID[1].",UTR";
							}
						}
					}
					
#					my $header = "CHR,CHR_No,SNP,bp,bp cumlative,NMISS,beta,SE,R2,T,P,minus log10P,unadj,GC,BONF,HOLM,SIDAK_SS,SIDAK_SD,FDR_BH,FDR_BY,gene_id,region\n";
					my $each_r = $ID_ori.",".$B[0].",".$B[1].",".$B[2].",".$cumlen.",".$B[3].",".$B[4].",".$B[5].",".$B[6].",".$B[7].",".$B[8].",".$mlog10p;
					$each_r .= ",".$hash->{$B[1]}{UNADJ}.",".$hash->{$B[1]}{GC}.",".$hash->{$B[1]}{BONF}.",".$hash->{$B[1]}{HOLM}.",".$hash->{$B[1]}{SIDAK_SS}.",".$hash->{$B[1]}{SIDAK_SD};
					$each_r .= ",".$hash->{$B[1]}{FDR_BH}.",".$hash->{$B[1]}{FDR_BY}.",".$gid_info;
					$r .= $each_r."\n";
					
#					my $header2 = "CHR,bp,NMISS,beta,SE,R2,T,P,minus log10P,unadj,GC,BONF,HOLM,SIDAK_SS,SIDAK_SD,FDR_BH,FDR_BY\n";
					my $each_info = $header2;
					$each_info .= $ID_ori.",".$B[2].",".$B[3].",".$B[4].",".$B[5].",".$B[6].",".$B[7].",".$B[8].",".$mlog10p;
					$each_info .= ",".$hash->{$B[1]}{UNADJ}.",".$hash->{$B[1]}{GC}.",".$hash->{$B[1]}{BONF}.",".$hash->{$B[1]}{HOLM}.",".$hash->{$B[1]}{SIDAK_SS}.",".$hash->{$B[1]}{SIDAK_SD};
					$each_info .= ",".$hash->{$B[1]}{FDR_BH}.",".$hash->{$B[1]}{FDR_BY}."\n";
					$cnt_with_stat++;
					
					my @tmp;
					push(@tmp, $ID_ori);
					push(@tmp, $B[2]);
					push(@tmp, $each_info);
					push(@tmp, $mlog10p);
					push(@tmp, $B[4]);
					push(@{$AoA}, \@tmp);
					
#					unless($hdata->{$ID_ori}){
#						$hdata->{$ID_ori} = "cumlen,mlog10p\n";
#					}
					
					if($B[4] > 0){
						$hdata->{$ID_ori}{mlog10p_posi} .= $cumlen.",".$mlog10p."\n";
					}
					elsif($B[4] < 0){
						$hdata->{$ID_ori}{mlog10p_nega} .= $cumlen.",".$mlog10p."\n";
					}
					
					$hdata->{$ID_ori}{mlog10p} .= $cumlen.",".$mlog10p."\n";
					$hdata->{$ID_ori}{mlog10p_dr} .= $cumlen.",".$mlog10p_dr."\n";
					$hdata->{num_var} += 1;
					
					if(! $hdata->{max}){
						$hdata->{max} = $mlog10p;
					}
					elsif($mlog10p > $hdata->{max}){
						$hdata->{max} = $mlog10p;
					}
					
					if(! $hdata->{max_beta}){
						$hdata->{max_beta} = $B[4];
					}
					elsif($B[4] > $hdata->{max_beta}){
						$hdata->{max_beta} = $B[4];
					}
					
					if(! $hdata->{min_beta}){
						$hdata->{min_beta} = $B[4];
					}
					elsif($B[4] < $hdata->{min_beta}){
						$hdata->{min_beta} = $B[4];
					}
				}
				else{
#					my $mlog10p = "NA";
#					my $tester = $B[8];
#					$tester =~ s/\.//g;
#					$tester =~ s/-//g;
#					if($tester !~ /\D/){
#						$mlog10p = - log($B[8]) / log(10);
#					}
#					$r .= $ID_ori.",".$B[0].",".$B[1].",".$B[2].",".$cumlen.",".$B[3].",".$B[4].",".$B[5].",".$B[6].",".$B[7].",".$B[8].",".$mlog10p;
#					$r .= ",-,-,-,-,-,-,-,-\n";
					$cnt_wo_stat++;
				}
			}
			elsif($data_type eq 'case'){
				if(! $hfasta->{$B[0]}){
					print "! unable to find refseq information for [$B[0]] (2), stop script...\n";
					print "\n$data_type\n$l\n\n";
					$hdata->{err} = 1;
					close $fh;
					return $hdata;
				}
				elsif($numB != 10){
					print "\n$data_type\n$l\n\n";
					die;
				}
				
				my $ID_ori = $hfasta->{$B[0]}{originalID};
				my $cumlen = $hfasta->{$B[0]}{cumlen0} + $B[2];
#				my $gnr = $B[3] / $nsample * 100;
				
#				if($hash->{$B[1]} && $gnr >= $gn_th){
				if($hash->{$B[1]} && $B[8] ne 'NA'){
					my $mlog10p = - log($B[8]) / log(10);
					
							# $numB == 10 (case)
					#	 CHR                             SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
					#	   1                          1_8865       8865    G  0.06452   0.1154    A        1.265       0.2607       0.5287 
					#	   1                         1_25940      25940   AG    0.375   0.4237    A       0.4702       0.4929        0.816 
					#	   1                         1_28655      28655    T    0.175  0.09091    C        3.534      0.06011        2.121 
					#	   1                         1_32929      32929    T   0.1953        0    C        35.44    2.637e-09           NA 
					
					my $gid_info = "-,-";
					if($B[1] =~ /__/){
						my @GID = split(/__/, $B[1]);
						if($GID[1] && $GID[2]){
							if($GID[2] eq 'p'){
								$gid_info = $GID[1].",promoter";
							}
							elsif($GID[2] eq 'g'){
								$gid_info = $GID[1].",gene";
							}
							elsif($GID[2] eq 'u'){
								$gid_info = $GID[1].",UTR";
							}
						}
					}
					
#					my $header = "CHR,CHR_No,SNP,bp,bp cumlative,NMISS,beta,SE,R2,T,P,minus log10P,unadj,GC,BONF,HOLM,SIDAK_SS,SIDAK_SD,FDR_BH,FDR_BY,gene_id,region\n";
					my $each_r = $ID_ori.",".$B[0].",".$B[1].",".$B[2].",".$cumlen.",-,-,-,-,-,".$B[8].",".$mlog10p;
					$each_r .= ",".$hash->{$B[1]}{UNADJ}.",".$hash->{$B[1]}{GC}.",".$hash->{$B[1]}{BONF}.",".$hash->{$B[1]}{HOLM}.",".$hash->{$B[1]}{SIDAK_SS}.",".$hash->{$B[1]}{SIDAK_SD};
					$each_r .= ",".$hash->{$B[1]}{FDR_BH}.",".$hash->{$B[1]}{FDR_BY}.",".$gid_info;
					$r .= $each_r."\n";
					
#					my $header2 = "CHR,bp,NMISS,beta,SE,R2,T,P,minus log10P,unadj,GC,BONF,HOLM,SIDAK_SS,SIDAK_SD,FDR_BH,FDR_BY\n";
					my $each_info = $header2;
					$each_info .= $ID_ori.",".$B[2].",-,-,-,-,-,".$B[8].",".$mlog10p;
					$each_info .= ",".$hash->{$B[1]}{UNADJ}.",".$hash->{$B[1]}{GC}.",".$hash->{$B[1]}{BONF}.",".$hash->{$B[1]}{HOLM}.",".$hash->{$B[1]}{SIDAK_SS}.",".$hash->{$B[1]}{SIDAK_SD};
					$each_info .= ",".$hash->{$B[1]}{FDR_BH}.",".$hash->{$B[1]}{FDR_BY}."\n";
					$cnt_with_stat++;
					
					my @tmp;
					push(@tmp, $ID_ori);
					push(@tmp, $B[2]);
					push(@tmp, $each_info);
					push(@tmp, $mlog10p);
					push(@tmp, "-");
					push(@{$AoA}, \@tmp);
					
#					unless($hdata->{$ID_ori}){
#						$hdata->{$ID_ori} = "cumlen,mlog10p\n";
#					}
					
					$hdata->{$ID_ori}{mlog10p} .= $cumlen.",".$mlog10p."\n";
					$hdata->{num_var} += 1;
					
					if(! $hdata->{max}){
						$hdata->{max} = $mlog10p;
					}
					elsif($mlog10p > $hdata->{max}){
						$hdata->{max} = $mlog10p;
					}
				}
				else{
#					my $mlog10p = "NA";
#					my $tester = $B[8];
#					$tester =~ s/\.//g;
#					$tester =~ s/-//g;
#					if($tester !~ /\D/){
#						$mlog10p = - log($B[8]) / log(10);
#					}
#					$r .= $ID_ori.",".$B[0].",".$B[1].",".$B[2].",".$cumlen.",".$B[3].",".$B[4].",".$B[5].",".$B[6].",".$B[7].",".$B[8].",".$mlog10p;
#					$r .= ",-,-,-,-,-,-,-,-\n";
					$cnt_wo_stat++;
				}
			}
		}
		$cnt++;
	}
	close $fh;
	
	if($data_type eq 'numeric'){
		if(abs($hdata->{max_beta}) > abs($hdata->{min_beta})){
			$hdata->{min_beta} = - $hdata->{max_beta};
		}
		elsif(abs($hdata->{max_beta}) < abs($hdata->{min_beta})){
			$hdata->{max_beta} = - $hdata->{min_beta};
		}
	}
	
	print "! [$cnt_with_stat] variants with additional stats (e.g. FDR_BH)\n";
	print "! [$cnt_wo_stat] skipped due to the low genotyping ratio or missing stats...\n";
	
	open(my $rfh, ">", $each_summary);
	print $rfh $r;
	close $rfh;
	print "! output [$each_summary]\n";
	
	@{$AoA} = sort {$b->[3] <=> $a->[3]} @{$AoA};
	my $num_AoA = @{$AoA};
	if($num_AoA > $topx){
		$num_AoA = $topx;
	}
	
	for(my $i = 0; $i < $num_AoA; $i++){
		my $sid = $AoA->[$i][0];
		my $bp = $AoA->[$i][1];
		$htopx->{$sid}{$bp}{info} = $AoA->[$i][2];
		$htopx->{$sid}{$bp}{mlog10p} = $AoA->[$i][3];
		$htopx->{$sid}{$bp}{beta} = $AoA->[$i][4];
		$htopx->{$sid}{$bp}{order} = $i;
		$htopx->{$sid}{$bp}{n} += 1;
	}
}

return ($hdata, $htopx);
}


#-------------------------------------------------------------------------------
sub Prep_phenodata{
my ($hash, $i, $phenoname, $rfile) = @_;

my $nsample = $hash->{nsample};
my $icol = $hash->{pheno_icol}{$i};

my $phenostr = "FID\tIID\t$phenoname\n";
for(my $j = 0; $j < $nsample; $j++){
	if(! defined $hash->{pheno_val}{$icol}{$j} || $hash->{pheno_val}{$icol}{$j} eq '-'){
		$hash->{pheno_val}{$i}{$j} = "NA";
	}
	my $name = $hash->{n2ID}{$j};
	my $val = $hash->{pheno_val}{$icol}{$j};
	$phenostr .= "$name\t1\t$val\n";
}

if($phenostr){
	open(my $rfh, ">", $rfile) or die;
	print $rfh $phenostr;
	close $rfh;
}

}


#-------------------------------------------------------------------------------
sub Read_phenodata3{
my $file = shift;

print "! reading [$file]...\n";
unless($file){
	return;
}

open(my $fh, "<", $file) or return;
my $hash = {};
my $cnt = 0;
my $num_sample = 0;
my $num_pheno = 0;
my $num_missing = 0;
my $numA = 0;
my $id_missing;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	unless($line){
		next;
	}
	else{
		my @A = split(/,/, $line);
		if($cnt == 0){
			$numA = @A;
			
			for(my $i = 1; $i < $numA; $i++){
#				print " [$A[$i]]\n";
				$hash->{pheno_name}{$num_pheno} = $A[$i];
				$hash->{pheno_icol}{$num_pheno} = $i;
				$hash->{pheno_i2n}{$i} = $num_pheno;
				$num_pheno++;
			}
		}
		elsif($cnt == 1){			# data type
			if($A[0] ne 'data_type'){
				print "! format error : expected 'data_type' but distinct string found [$A[0]]\n";
				die;
			}
			
			for(my $i = 1; $i < $numA; $i++){
				if(! $A[$i]){
					$A[$i] = "null";
				}
				elsif($A[$i] ne 'numeric' && $A[$i] ne 'case'){
					$A[$i] = "null";
				}
				
				my $icol = $hash->{pheno_i2n}{$i};
				$hash->{pheno_datatype}{$icol} = $A[$i];
			}
		}
		elsif($cnt == 2){			# group low threshold
			if($A[0] ne 'group_low'){
				print "! format error : expected 'group_low' but distinct string found [$A[0]]\n";
				die;
			}
			
			for(my $i = 1; $i < $numA; $i++){
				if(! defined $A[$i]){
					$A[$i] = "null";
				}
				
				my $icol = $hash->{pheno_i2n}{$i};
				my $dtype = $hash->{pheno_datatype}{$icol};
				if($dtype eq 'numeric'){
					my $tester = $A[$i];
					$tester =~ s/\.//;
					$tester =~ s/-//;
					if(! $tester || $tester =~ /\D/){
						$A[$i] = "null";
					}
				}
				$hash->{pheno_low}{$icol} = $A[$i];
			}
		}
		elsif($cnt == 3){			# group high threshold
			if($A[0] ne 'group_high'){
				print "! format error : expected 'group_high' but distinct string found [$A[0]]\n";
				die;
			}
			
			for(my $i = 1; $i < $numA; $i++){
				if(! defined $A[$i]){
					$A[$i] = "null";
				}
				
				my $icol = $hash->{pheno_i2n}{$i};
				my $dtype = $hash->{pheno_datatype}{$icol};
				if($dtype eq 'numeric'){
					my $tester = $A[$i];
					$tester =~ s/\.//;
					$tester =~ s/-//;
					if(! $tester || $tester =~ /\D/){
						$A[$i] = "null";
					}
				}
				$hash->{pheno_high}{$icol} = $A[$i];
			}
		}
		else{
			my $judge = "true";
			if(! $A[0] || $A[0] eq '-'){
				$judge = "false";
			}
			if($judge eq 'false'){
				next;
			}
			
			$hash->{n2ID}{$num_sample} = $A[0];
			$hash->{ID2n}{$A[0]} = $num_sample;
			
			for(my $i = 1; $i < $numA; $i++){
				if(! defined $A[$i]){
					$A[$i] = "-";
					 $num_missing++;
				}
				$hash->{pheno_val}{$i}{$num_sample} = $A[$i];
			}
#			print " [$A[0]]\n";
			$num_sample++;
		}
	}
	$cnt++;
}
close $fh;

$hash->{nsample} = $num_sample;
$hash->{npheno} = $num_pheno;

if($id_missing){
	print "! missing some sample names, abort script...\n";
	$hash->{err} = 1;
}
else{
	print " - [$num_sample] ID x [$num_pheno] phenotype found\n";
#	print " - [$num_missing] missing pheno value\n";
}

return $hash;
}


#----------------------------------------------------------------------
sub LnFile{
my $file = shift;

if($file && $file ne 'null'){
	my @tmp = split(/\//, $file);
	my $n = @tmp;
	
	if(-e $tmp[$n - 1]){
		system("rm $tmp[$n - 1]");
	}
	
	if(-e $file){
		system("ln -s $file ./");
	}
	elsif(-e "../$file"){
		system("ln -s ../$file ./");
	}
	
	return $tmp[$n - 1];
}

}


#-------------------------------------------------------------------------------
sub Stats_fasta{
my $file = shift;

print "! reading [$file] as base data...\n";
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
			
			if($ID !~ /scaffold/ && $ID !~ /mitoch/ && $ID !~ /chloro/){
				my $len = length($seq);
				my $ID_ori = $ID;
				$ID =~ s/chr0//i;
				$ID =~ s/chr//i;
				$ID =~ s/Gm0//i;
				$ID =~ s/Gm//i;
				$hash->{$ID}{originalID} = $ID_ori;
				$hash->{$ID}{order} = $norder;
	#			$hash->{$ID}{seq} = $seq;
				$hash->{$ID}{len} = $len;
				$hash->{$ID}{cumlen0} = $total_len;
				$hash->{$ID_ori}{cumlen0} = $total_len;
				$total_len += $len;
				$hash->{$ID}{cumlen1} = $total_len;
				$hash->{$ID_ori}{cumlen1} = $total_len;
				push(@SID, $ID);
				push(@SID_ori, $ID_ori);
				$norder++;
				print " [$ID_ori] = [$len] bp\n";
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
	if($hash->{$ID}){
		print "! duplicated seq ID : [$ID]\n";
	}
	
	if($ID !~ /scaffold/ && $ID !~ /mitoch/ && $ID !~ /chloro/){
		my $len = length($seq);
		my $ID_ori = $ID;
		$ID =~ s/chr0//i;
		$ID =~ s/chr//i;
		$ID =~ s/Gm0//i;
		$ID =~ s/Gm//i;
		$hash->{$ID}{originalID} = $ID_ori;
		$hash->{$ID}{order} = $norder;
	#	$hash->{$ID}{seq} = $seq;
		$hash->{$ID}{len} = $len;
		$hash->{$ID}{cumlen0} = $total_len;
		$hash->{$ID_ori}{cumlen0} = $total_len;
		$total_len += $len;
		$hash->{$ID}{cumlen1} = $total_len;
		$hash->{$ID_ori}{cumlen1} = $total_len;
		push(@SID, $ID);
		push(@SID_ori, $ID_ori);
		$norder++;
		print " [$ID_ori] = [$len] bp\n";
	}
}

my $numID = @SID;
$hash->{SID} = \@SID;
$hash->{SID_ori} = \@SID_ori;
$hash->{total_len} = $total_len;

print "! total [$numID] sequence ID, [$total_len] bp\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Mv_files{
my $A = shift;
my $dir = shift;

foreach my $file (@{$A}){
	if(-e $file){
		system("mv $file $dir");
	}
}

}


#-------------------------------------------------------------------------------
sub Rm_files{
my $A = shift;

foreach my $file (@{$A}){
	if(-e $file){
		system("rm $file");
	}
}

}


#-------------------------------------------------------------------------------
sub Open_fasta_as_hash{
my $file = shift;

print "! open [$file] as hash ...\n";
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
sub R_GenPlot{
my ($hfasta, $hdata, $PlotColors, $num_color, $color_abline, $tmpdir_Rplots, $hfplots, $png_width, $png_height, $pref_assoc, $wpath) = @_;

my $each_plot1 = $hfplots->{each_plot1};
my $each_plot1_posi = $hfplots->{each_plot1_posi};
my $each_plot1_nega = $hfplots->{each_plot1_nega};
my $each_plot2 = $hfplots->{each_plot2};

my $SID = $hfasta->{SID};
my $SID_ori = $hfasta->{SID_ori};
my $numSID = @{$SID};
my $hash = {};
my $ablines = "";
for(my $i = 0; $i < $numSID; $i++){
	my $sid = $SID->[$i];
	my $sid_ori = $SID_ori->[$i];
	if($sid_ori =~ /scaffold/){
		if($sid_ori =~ /chloro/ || $sid_ori =~ /mitoch/){
			next;
		}
	}
	
	my $cumlen0 = $hfasta->{$sid}{cumlen0};
	my $cumlen1 = $hfasta->{$sid}{cumlen1};
	
	if($i == 0){
		$ablines .= "abline(v=$cumlen0, col=\"$color_abline\", lwd=1)\n";
		$ablines .= "abline(v=$cumlen1, col=\"$color_abline\", lwd=1)\n";
	}
	elsif($i < $numSID - 1){
		$ablines .= "abline(v=$cumlen1, col=\"$color_abline\", lwd=1)\n";
	}
	
#	my $ncolor = $i % $num_color;
	my $ncolor = $i;
	if(! $PlotColors->[$ncolor]){
		print "! error : missing plot color for [$sid] i=[$i]...\n";
	}
	
	unless($hdata->{$sid_ori}){
		print "! missing data for [$sid]...\n";
	}
	else{
		$hash->{$ncolor}{mlog10p} .= $hdata->{$sid_ori}{mlog10p};
		
		if($hdata->{$sid_ori}{mlog10p_posi}){
			$hash->{$ncolor}{mlog10p_posi} .= $hdata->{$sid_ori}{mlog10p_posi};
		}
		if($hdata->{$sid_ori}{mlog10p_nega}){
			$hash->{$ncolor}{mlog10p_nega} .= $hdata->{$sid_ori}{mlog10p_nega};
		}
		if($hdata->{$sid_ori}{mlog10p_dr}){
			$hash->{$ncolor}{mlog10p_dr} .= $hdata->{$sid_ori}{mlog10p_dr};
		}
		$hash->{$ncolor}{color} = $PlotColors->[$ncolor];
	}
}

my @NC = keys(%{$hash});
@NC = sort {$a <=> $b} @NC;

my $Rscript =<<"EOS";
workingDir = "$wpath"
setwd(workingDir)
EOS

my $Rscript_plot1 = "";
my $Rscript_plot1_posi = "";
my $Rscript_plot1_nega = "";
my $Rscript_plot2 = "";
my $numNC = @NC;
for(my $i = 0; $i < $numNC; $i++){
	my $ncolor = $NC[$i];
	my $lines1 = "cumlen,value\n".$hash->{$ncolor}{mlog10p};
#	my $lines1_posi = "cumlen,value\n".$hash->{$ncolor}{mlog10p_posi};
#	my $lines1_nega = "cumlen,value\n".$hash->{$ncolor}{mlog10p_nega};
#	my $lines2 = "cumlen,value\n".$hash->{$ncolor}{mlog10p_dr};
	my $tmpcsv1 = $tmpdir_Rplots."/_plotdata_-log10p_".$pref_assoc."_".$ncolor.".csv";
	my $tmpcsv1_posi = $tmpdir_Rplots."/_plotdata_-log10p_posi_".$pref_assoc."_".$ncolor.".csv";
	my $tmpcsv1_nega = $tmpdir_Rplots."/_plotdata_-log10p_nega_".$pref_assoc."_".$ncolor.".csv";
	my $tmpcsv2 = $tmpdir_Rplots."/_plotdata_-log10p_dr_".$pref_assoc."_".$ncolor.".csv";
	
	open(my $rfh1, ">", $tmpcsv1) or die;
	print $rfh1 $lines1;
	close $rfh1;
	
#	open(my $rfh1_posi, ">", $tmpcsv1_posi) or die;
#	print $rfh1_posi $lines1_posi;
#	close $rfh1_posi;
#	
#	open(my $rfh1_nega, ">", $tmpcsv1_nega) or die;
#	print $rfh1_nega $lines1_nega;
#	close $rfh1_nega;
#	
#	open(my $rfh2, ">", $tmpcsv2) or die;
#	print $rfh2 $lines2;
#	close $rfh2;
	
	$Rscript .= "mlog10p_$ncolor <- as.matrix(read.csv(\"$tmpcsv1\", header=T))\n";
#	$Rscript .= "mlog10p_posi_$ncolor <- as.matrix(read.csv(\"$tmpcsv1_posi\", header=T))\n";
#	$Rscript .= "mlog10p_nega_$ncolor <- as.matrix(read.csv(\"$tmpcsv1_nega\", header=T))\n";
#	$Rscript .= "mlog10p_dr_$ncolor <- as.matrix(read.csv(\"$tmpcsv2\", header=T))\n";
	
	my $dotcolor_beta = "#898880";
	if($i%2 == 1){
		$dotcolor_beta = "#7e837f";
	}
	
	if(! defined $hdata->{min}){
		$hdata->{min} = - $hdata->{max};
	}
	
	if($i == 0){
		$Rscript_plot1 .= "plot(mlog10p_$ncolor".", pch=19, col=\"$hash->{$ncolor}{color}\", xlab =\"Genomic position\", ylab=\"-log10(p-value)\", cex.lab=1, xlim=c(0, $hfasta->{total_len}), ylim=c(0, $hdata->{max}))\n";
		$Rscript_plot1 .= "title(main=\"$pref_assoc\", cex.main=2)\n";
		
		$Rscript_plot1_posi .= "plot(mlog10p_posi_$ncolor".", pch=19, col=\"$hash->{$ncolor}{color}\", xlab =\"Genomic position\", ylab=\"-log10(p-value)\", cex.lab=1, xlim=c(0, $hfasta->{total_len}), ylim=c(0, $hdata->{max}))\n";
		$Rscript_plot1_posi .= "title(main=\"$pref_assoc (beta > 0)\", cex.main=2)\n";
		
		$Rscript_plot1_nega .= "plot(mlog10p_nega_$ncolor".", pch=19, col=\"$hash->{$ncolor}{color}\", xlab =\"Genomic negation\", ylab=\"-log10(p-value)\", cex.lab=1, xlim=c(0, $hfasta->{total_len}), ylim=c(0, $hdata->{max}))\n";
		$Rscript_plot1_nega .= "title(main=\"$pref_assoc (beta < 0)\", cex.main=2)\n";
		
		$Rscript_plot2 .= "plot(mlog10p_dr_$ncolor".", pch=19, col=\"$hash->{$ncolor}{color}\", xlab =\"Genomic position\", ylab=\"-log10(p-value) with direction\", cex.lab=1, xlim=c(0, $hfasta->{total_len}), ylim=c($hdata->{min}, $hdata->{max}))\n";
		$Rscript_plot2 .= "title(main=\"$pref_assoc\", cex.main=2)\n";
	}
	else{
		$Rscript_plot1 .= "points(mlog10p_$ncolor".", pch=19, col=\"$hash->{$ncolor}{color}\")\n";
		$Rscript_plot1_posi .= "points(mlog10p_posi_$ncolor".", pch=19, col=\"$hash->{$ncolor}{color}\")\n";
		$Rscript_plot1_nega .= "points(mlog10p_nega_$ncolor".", pch=19, col=\"$hash->{$ncolor}{color}\")\n";
		$Rscript_plot2 .= "points(mlog10p_dr_$ncolor".", pch=19, col=\"$hash->{$ncolor}{color}\")\n";
	}
}

my $bonf_th = - log(0.05 / $hdata->{num_var}) / log(10);

$Rscript .=<<"EOS";
png("$each_plot1", width=$png_width, height=$png_height)
$Rscript_plot1
$ablines
abline(h=0, col="$color_abline", lwd=1)
abline(h=$bonf_th, col="red", lwd=1, lty=5)
dev.off()

EOS

my $notuse =<<"EOS";
png("$each_plot1_posi", width=$png_width, height=$png_height)
$Rscript_plot1_posi
$ablines
abline(h=0, col="$color_abline", lwd=1)
abline(h=$bonf_th, col="red", lwd=1, lty=5)
dev.off()

png("$each_plot1_nega", width=$png_width, height=$png_height)
$Rscript_plot1_nega
$ablines
abline(h=0, col="$color_abline", lwd=1)
abline(h=$bonf_th, col="red", lwd=1, lty=5)
dev.off()

EOS

my $Rsc_file = $tmpdir_Rplots."/_Rscript_".$pref_assoc.".txt";
if($Rscript){
	open(my $rfh, ">", $Rsc_file) or die;
	print $rfh $Rscript;
	close $rfh;
}

system("Rscript --vanilla --slave $Rsc_file > /dev/null 2>&1");
#system("Rscript --vanilla --slave $Rsc_file");
print "! created manhattan plots\n";

#if(-e $Rsc_file){
#	system("rm $Rsc_file");
#}
if(-e "Rplots.pdf"){
	system("rm Rplots.pdf");
}

}


#-------------------------------------------------------------------------------
sub Prep_plotColors{
my $n = shift;

my @PlotColors;
my $num_color = 1;
if($n == 1){
	my @tmp;
	push(@tmp, "#808080");	# gray
	push(@tmp, "#a9a9a9");	# darkgray
	my $ntmp = @tmp;
	$num_color = $ntmp;
	
	for(my $i = 0; $i < 100; $i += $ntmp){
		for(my $j = 0; $j < $ntmp; $j++){
			push(@PlotColors, $tmp[$j]);
		}
	}
}
elsif($n == 2){
	my @tmp;
	push(@tmp, "#0068b7");
	push(@tmp, "#f39800");
	push(@tmp, "#898989");
	push(@tmp, "#c1ab05");
	push(@tmp, "#4496d3");
	push(@tmp, "#9cbb1c");
	push(@tmp, "#192f60");
	push(@tmp, "#6b3f31");
	push(@tmp, "#626063");
	push(@tmp, "#c4972f");
	push(@tmp, "#00608d");
	push(@tmp, "#387d39");
	push(@tmp, "#68a4d9");
	push(@tmp, "#f6ae54");
	push(@tmp, "#abb1b5");
	push(@tmp, "#ffdc00");
	push(@tmp, "#00afcc");
	push(@tmp, "#89c997");
	push(@tmp, "#434da2");
	push(@tmp, "#ac6b25");
	push(@tmp, "#504946");
	push(@tmp, "#e29676");
	my $ntmp = @tmp;
	$num_color = $ntmp;
	
	for(my $i = 0; $i < 100; $i += $ntmp){
		for(my $j = 0; $j < $ntmp; $j++){
			push(@PlotColors, $tmp[$j]);
		}
	}
}

return (\@PlotColors, $num_color);
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

print "! [$gcnt] genes, [$tcnt] transcripts, [$pcnt] CDS lines\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Check_previous{
my $file = shift;
my $id = shift;

my $sw = "false";
if(-e $file){
	open(my $fh, "<", $file) or die;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		if($id eq $line){
			$sw = "true";
		}
	}
	close $fh;
}

return $sw;
}


#-------------------------------------------------------------------------------
sub Batch_exe{
my $script = shift;
my $j = shift;
my $cnt = shift;

if($script && -e $script){
	print " thread [$j] : [$cnt] queries...\n";
	if(system("bash $script") == 0){
		print " thread [$j] : completed.\n";
		
		unless(-e "cmds"){
			system("mkdir cmds");
		}
		system("mv $script cmds");
	}
	else{
		print "! thread [$j] : failed\n";
		die;
	}
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


#---------------------------------------------------------------------
sub Search_app2{
my $hbin = shift;
my $APP = shift;
my $hpath = shift;
my $script_path = shift;
my $pipeline = shift;
my $genelist = shift;

$script_path =~ s/\/scripts//;
my @tmp = split(/\//, $genelist);
my $ntmp = @tmp;
$genelist = $tmp[$ntmp - 1];

my $tmp = "_Search_app_temp.$genelist.txt";
for(my $i = 0; $i < 100; $i++){
	if(-e $tmp){
		$tmp = "_Search_app_temp.$genelist.$i.txt";
	}
	if(! -e $tmp){
		last;
	}
}

unless($pipeline){
	print "! checking PATH of required programs...\n";
}
my $hash = {};
foreach my $app (@{$APP}){
	my $Bin = $hbin->{$app};
	if($hpath->{$app} && -e $hpath->{$app}){
		my @F = glob("$hpath->{$app}/*");
		foreach my $bin (@{$Bin}){
			foreach my $f (@F){
				if($f =~ /$bin/){
					$hash->{$bin} = $f;
				}
			}
		}
	}
	else{
		foreach my $bin (@{$Bin}){
			system("which $bin > $tmp 2>&1");
			unless($bin =~ /samtools/){
				system("which $bin.pl >> $tmp 2>&1");
			}
			system("find $script_path/bin | grep $bin >> $tmp 2>&1");
			
			open(my $fh, "<", $tmp);
			while(my $line = <$fh>){
				$line =~ s/\n//;
				$line =~ s/\r//;
				
				if($line && $line =~ /$bin/i && $line !~ /which: no/){
					if($line =~ /\//){
						my @D = split(/\//, $line);
						my $nD = @D;
						$D[$nD - 1] =~ s/\.pl//;
						if($D[$nD - 1] eq $bin){
							$hash->{$bin} = $line;
						}
					}
					else{
						$hash->{$bin} = $line;
					}
				}
			}
			close $fh;
			
			system("rm $tmp");
		}
	}
	
	foreach my $bin (@{$Bin}){
		if($hash->{$bin} && -e $hash->{$bin}){
			unless($pipeline){
				if($bin !~ /interpro/ && ! $pipeline){
					print " [$bin] = [$hash->{$bin}] <OK>\n";
				}
				elsif($bin =~ /interpro/){
					my $testcmd = "$hash->{$bin} --version > $tmp";
					system($testcmd);
					
					my $iprver;
					open(my $fh, "<", $tmp);
					while(my $line = <$fh>){
						$line =~ s/\n//;
						$line =~ s/\r//;
						
						if($line && $line =~ /InterProScan version/i){
							$line =~ s/InterProScan version//;
							$line =~ s/\s//g;
							my @tmp2 = split(/\-/, $line);
							if($tmp2[0] =~ /5\./){
								$iprver = $tmp2[0];
							}
						}
					}
					close $fh;
					
					system("rm $tmp");
					
					if(! $iprver || $iprver !~ /5\./){
						print " NOT found [$bin], please check PATH\n";
						$hash->{err} += 1;
					}
					elsif($iprver =~ /5\./){
						print " [$bin] = [$hash->{$bin}] <OK> (version $iprver)\n";
					}
				}
			}
		}
		else{
			print " NOT found [$bin], please check PATH\n";
			$hash->{err} += 1;
		}
	}
}

unless($pipeline){
	print "\n";
}

return $hash;
}



