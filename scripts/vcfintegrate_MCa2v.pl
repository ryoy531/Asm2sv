#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use bytes;
use Encode;
use FindBin;
use Getopt::Long;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "vcfintegrate_MCa2v.pl, version 1.01\n";
$version .= "last update: [2024\/4\/29]\n";

print "$version\n";

#-------------------------------------------------------------------------------

my $wpath = getcwd();
my $script_path = $FindBin::Bin;
my $bin_sub = $script_path."/sub_vcfMCa2v.pl";
my $bin_varcomb = $script_path."/sub_varcomb.pl";
if(! -e $bin_sub){
	print "! unable to find [$bin_sub], abort script...\n";
	goto END;
}
if(! -e $bin_varcomb){
	print "! unable to find [$bin_varcomb], abort script...\n";
	goto END;
}

my $hbin = {};
$hbin->{vg}[0] = "vg";
$hbin->{matcher}[0] = "matcher";

my $hbinpath = Search_app2($hbin, $script_path, 0);
print "\n";

my $datalist;
my $basevcf;
my $rfile;
my $policy;
my $cpu;
my $npart;
my $help;
my $rdir;
my $datapath;
my $debug;

GetOptions('--list=s' => \$datalist, '--base=s' => \$basevcf, '--path=s' => \$datapath, '--mode=s' => \$policy, '--thread=i' => \$cpu, '--output_dir=s' => \$rdir, '--help' => \$help);

my $err = 0;
if(! $basevcf){
	print "! --base [VCF file] not specified...\n";
}
elsif(! -e $basevcf){
	print "! missing --base [VCF file]...\n";
	$err++;
}
if(! $datalist || ! -e $datalist){
	print "! missing --list [datalist csv]\n";
	$err++;
}
if(! $policy){
	$policy = "mergedupSVs";
}
elsif($policy ne 'mergedupSVs' && $policy ne 'all' && $policy ne 'keep_larger2'){
	$policy = "mergedupSVs";
}
if(! $cpu || $cpu =~ /\D/){
	$cpu = 24;
}
if($err > 0){
	Help();
	goto END;
}
if(! $datapath || ! -e $datapath){
	$datapath = "null";
}

if(! $rdir){
	$rdir = "integratedVCF_for_pangenome";
}
unless(-e $rdir){
	system("mkdir $rdir");
}
#chdir $rdir;
#unless(-e "log"){
#	system("mkdir log");
#}

my ($VF, $VF_str, $err_list, $reffasta) = Read_datalist($datalist, $datapath, "check");
if($err_list > 0){
	print "! missing gene-SV VCF from Asm2sv, please check each result...\n";
	goto END;
}
if(! -e $reffasta){
	print "! missing reference fasta...\n";
	goto END;
}

my $num_VF = @{$VF};
my $list = "./".$rdir."/geneSV_list.tsv";
SAVE($list, $VF_str);

my ($hdbseq, $SID) = Open_fasta_as_hash($reffasta);
my $num_SID = @{$SID};
if(! $cpu || $cpu =~ /\D/){
	$cpu = $num_SID;
}
my $rprefix = "combined_geneSV_".$num_VF."genomes";
my $rfile0 = $rprefix.".vcf";
my $prefixb = Path2pref($basevcf);

print "\n! integrating data in each seq entry...\n";
my @eachRF1;
my @eachRF2;
my $d = int($num_SID / $cpu) + 1;
for(my $i = 0; $i < $d; $i++){
	my $p0 = $i * $cpu;
	my $p1 = ($i + 1) * $cpu;
	
	my $thrs = [];
	for(my $j = $p0; $j < $p1; $j++){
		if($SID->[$j]){
			if($SID->[$j] =~ /scaffold/i){
				if($SID->[$j] =~ /chloro/i || $SID->[$j] =~ /mitoch/i){
					next;
				}
			}
			
			my $ldir = "./".$rdir."/log";
			unless(-e $ldir){
				system("mkdir $ldir");
			}
			
			my $each_prefix = "each_multisample_".$SID->[$j];
			my $rfile1 = $ldir."/".$prefixb."_plus_".$each_prefix.".vcf";
			my $rfile2 = $ldir."/".$each_prefix.".vcf";
			my $cmd = "$bin_sub -r $reffasta -b $basevcf -l $list -s $SID->[$j] -o $each_prefix -w $rdir -m $hbinpath->{matcher} -p $policy";
			my $sh = $ldir."/_batch_".$SID->[$j].".sh";
			my $logfile = $ldir."/_log_".$SID->[$j].".txt";
			
			if(! -e $rfile1){
				my $thr = threads->new(\&Batch_exec, $cmd, $sh, $logfile, $SID->[$j], $j);
				push(@{$thrs}, $thr);
				push(@eachRF1, $rfile1);
				push(@eachRF2, $rfile2);
			}
			else{
				print " thread [$j] : [$SID->[$j]], already done, skip...\n";
				push(@eachRF1, $rfile1);
				push(@eachRF2, $rfile2);
			}
		}
	}
	foreach my $thr (@{$thrs}){
		$thr->join();
	}
}

my $rfile_combined = "./$rdir/".$prefixb."_plus_".$rprefix.".vcf";
my $exf0 = $prefixb."_plus_".$rprefix.".vcf";
my $exf1 = $prefixb."_plus_".$rprefix.".vcf.gz";
if(-e $rfile_combined){
	print " - once delete [$rfile_combined]\n";
	system("rm $rfile_combined");
}

my $rfile_combined2 = "./$rdir/".$rprefix.".vcf";
if(-e $rfile_combined2){
	print " - once delete [$rfile_combined2]\n";
	system("rm $rfile_combined2");
}

print "\n! combine results to [$rfile_combined]...\n";
my $nf = 0;
foreach my $rfile1 (@eachRF1){
	if($rfile1 =~ /scaffold/i){
		if($rfile1 =~ /chloro/i || $rfile1 =~ /mitoch/i){
			next;
		}
	}
	if(-e $rfile1){
		CombineVCF($rfile1, $rfile_combined, $nf);
	}
	else{
		print "! missing [$rfile1]...\n";
	}
	$nf++;
}
if(-e $rfile_combined){
	my $rfile_combined_gz = $rfile_combined.".gz";
	my $cmd = "bgzip -c $rfile_combined > $rfile_combined_gz";
	print "! bgzip [$rfile_combined]...\n";
	if(system($cmd) != 0){
		print "! failed.\n";
	}
	else{
		print "! done.\n";
	}
}

my $skip2 = 1;
if(! $skip2){
	print "\n! combine results to [$rfile_combined2]...\n";
	$nf = 0;
	foreach my $rfile2 (@eachRF2){
		if($rfile2 =~ /scaffold/i){
			if($rfile2 =~ /chloro/i || $rfile2 =~ /mitoch/i){
				next;
			}
		}
		if(-e $rfile2){
			CombineVCF($rfile2, $rfile_combined2, $nf);
		}
		else{
			print "! missing [$rfile2]...\n";
		}
		$nf++;
	}
	if(-e $rfile_combined2){
		my $rfile_combined2_gz = $rfile_combined2.".gz";
		my $cmd = "bgzip -c $rfile_combined2 > $rfile_combined2_gz";
		print "! bgzip [$rfile_combined2]...\n";
		if(system($cmd) != 0){
			print "! failed.\n";
		}
		else{
			print "! done.\n";
		}
	}
}

my $cmd_vgconst =<<"EOS";
$hbinpath->{vg} construct -r ../$reffasta -v $exf0 > examplePanRef.vg
$hbinpath->{vg} autoindex --workflow giraffe -r ../$reffasta -v $exf1 -p examplePanRef
$hbinpath->{vg} index -L -x examplePanRef.xg examplePanRef.vg
$hbinpath->{vg} snarls examplePanRef.xg > examplePanRef.snarls
$hbinpath->{vg} gbwt -o examplePanRef.gbwt -Z examplePanRef.giraffe.gbz
$hbinpath->{vg} stats -z examplePanRef.vg
EOS

SAVE("./$rdir/example_command_vgconst.sh", $cmd_vgconst);

END:{
	print "\n! End of script.\n\n";
}


################################################################################
#---------------------------------------------------------------------------------------
sub Help{
	my $help_command =<<"EOS";
Usage: 
  /path/to/Vcfintegrate_cactusMG_asm2sv.pl -b [base.vcf] -l [datalist csv] -m [SV reporting policy] -t [CPU thread] -o [output directory]

--list or -l          : data file list (PATH of fasta, GFF, psl, and other files)
                       *header should be "dbfasta" "dbgff" "genelist" "qfasta" "qpsl" "chralias"
                       
--base or -b          : base VCF (assume output of 'splitseq_run_cactus.pl')
--mode or -m          : 'mergedupSVs' (merge and simplify duplicating SVs at the same genomic location) or 'all' (native data)
--thread or -t        : CPU thread number (default 16)
--partition or -p     : num of parallel job (default 1)
EOS
	print "\n$help_command";
}


#-------------------------------------------------------------------------------
sub Read_datalist{
my $file = shift;
my $datapath = shift;
my $option = shift;
my $gid_alias = shift;

if(! $gid_alias || ! -e $gid_alias){
	$gid_alias = "null";
}

print "! reading data list [$file]...\n";
open(my $fh, "<", $file);
my @VF;
my $VF_str;
my @H;
my $cnt = 0;
my $err = 0;
my $nrfasta = {};
my $nrgff = {};
my $setref = "null";
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if(! $line || $line =~ /\#/){
		next;
	}
	else{
		my @A = split(/,/, $line);
		my $numA = @A;
		
		if($numA > 8){
			$numA = 8;
		}
		elsif($numA < 8){
			print "! expect column = 8 but only [$numA]...\n";
			$err++;
		}
		
		if($line =~ /dbfasta/ && $line =~ /qfasta/){
			# dbfasta	dbgff	genelist	qfasta	qgff	qpsl	outdir	custom_prefix
			
			for(my $i = 0; $i < $numA; $i++){
				push(@H, $A[$i]);
			}
			next;
		}
		else{
			$cnt++;
			my $missing = 0;
			for(my $i = 0; $i < $numA; $i++){
				if($H[$i] eq 'qpsl' || $H[$i] eq 'qgff' || $H[$i] eq 'chralias' || $H[$i] eq 'prev_result'){
					if(! $A[$i] || ! -e $A[$i]){
						$A[$i] = "null";
					}
				}
				elsif($H[$i] eq 'outdir'){
					if(! $A[$i] || $A[$i] eq '-' || $A[$i] eq 'NA' || $A[$i] eq 'na' || $A[$i] eq 'N/A' || $A[$i] eq 'n.a.'){
						$A[$i] = "null";
					}
				}
				elsif($H[$i] eq 'sample_prefix' || $H[$i] eq 'custom_prefix'){
					if(! $A[$i] || $A[$i] eq '-' || $A[$i] eq 'NA' || $A[$i] eq 'na' || $A[$i] eq 'N/A' || $A[$i] eq 'n.a.'){
						$A[$i] = "null";
					}
				}
				else{
					if($datapath eq 'null'){
						if(! -e $A[$i]){
							print " [$cnt] : missing [$A[$i]]\n";
							$missing++;
						}
					}
					else{
						if(! -e "$datapath/$A[$i]"){
							print " [$cnt] : missing [$datapath/$A[$i]]\n";
							$missing++;
						}
					}
				}
			}
			if($missing == 0){
				my $hash = {};
				for(my $i = 0; $i < $numA; $i++){
					$hash->{$H[$i]} = $A[$i];
				}
				
				my @tmp0 = split(/\//, $hash->{dbfasta});
				my $ntmp0 = @tmp0;
				my $dbpref = $tmp0[$ntmp0 - 1];
				$dbpref =~ s/\.fasta//;
				$dbpref =~ s/\.fa//;
				$hash->{dbpref} = $dbpref;
				
#				$nrfasta->{$tmp0[$ntmp0 - 1]} = 1;
				$nrfasta->{$dbpref} = 1;
				
				my @tmp1 = split(/\//, $hash->{qfasta});
				my $ntmp1 = @tmp1;
				my $qpref = $tmp1[$ntmp1 - 1];
				$qpref =~ s/\.fasta//;
				$qpref =~ s/\.fa//;
				$hash->{qpref} = $qpref;
				
				my @tmp2 = split(/\//, $hash->{dbgff});
				my $ntmp2 = @tmp2;
				my $dbgff = $tmp2[$ntmp2 - 1];
				
#				$nrgff->{$tmp2[$ntmp2 - 1]} = 1;
				$nrgff->{$dbgff} = 1;
				
				my @tmp3 = split(/\//, $hash->{qgff});
				my $ntmp3 = @tmp3;
				my $qgff = $tmp3[$ntmp3 - 1];
				
				if($hash->{outdir} eq 'null'){
					$hash->{outdir} = "result_asm2sv_".$dbpref."_vs_".$qpref;			# do not change this
				}
				
				if($option eq 'check'){
					if(! $hash->{custom_prefix}){
						$hash->{custom_prefix} = $qpref;
					}
					if($datapath eq 'null'){
						$hash->{geneSV} = $hash->{outdir}."/rev2sv/genebased_asm2sv_".$hash->{custom_prefix}.".vcf";
						$setref = $hash->{dbfasta};
						
						if(! -e $hash->{geneSV}){
							$hash->{geneSV} = $hash->{outdir}."/rev2sv/genebased_asm2sv_".$qpref.".vcf";
						}
					}
					else{
						$hash->{geneSV} = $datapath."/".$hash->{outdir}."/rev2sv/genebased_asm2sv_".$hash->{custom_prefix}.".vcf";
						$setref = $datapath."/".$hash->{dbfasta};
						
						if(! -e $hash->{geneSV}){
							$hash->{geneSV} = $datapath."/".$hash->{outdir}."/rev2sv/genebased_asm2sv_".$qpref.".vcf";
						}
					}
					
					if(-e $hash->{geneSV}){
						print " [$hash->{custom_prefix}] : v=[$hash->{geneSV}] <OK>\n";
						push(@VF, $hash->{geneSV});
						$VF_str .= "$hash->{custom_prefix}\t$hash->{geneSV}\n";
					}
					else{
						print " [$hash->{custom_prefix}] : missing [$hash->{geneSV}]...\n";
						$err = 1;
					}
				}
			}
			else{
				$err = 1;
			}
		}
	}
}
close $fh;

my @NRF = keys(%{$nrfasta});
my @NRG = keys(%{$nrgff});
my $num_NRF = @NRF;
my $num_NRG = @NRG;
my $unique_reffasta = "";
my $unique_refgff = "";

if($num_NRF == 1 && $num_NRG == 1){
	print "! dbfasta and dbgff are unique : [$NRF[0]] [$NRG[0]] <OK>\n";
	$unique_reffasta = $NRF[0];
	$unique_refgff = $NRG[0];
}
elsif($err == 0){
	print "\n! error : duplicated dbfasta and/or dbgff entries...\n";
	print "  - [$num_NRF] dbfasta : @NRF\n";
	print "  - [$num_NRG] dbgff : @NRG\n";
	$err = 1;
}
if($err == 0){
	print "! total [$cnt] samples\n";
}

return (\@VF, $VF_str, $err, $setref);
}


#--------------------------------------------------------------------------------
sub CombineVCF{
my ($file, $rfile_combined, $nf) = @_;

my $fh;
if($file =~ /\.gz/){
	open($fh, "gzip -dc $file |") or die;
}
else{
	open($fh, "<", $file) or die;
}

print " - collecting data [$file]...\n";
my $rvcf = "";
my $header = "";
my $cnt = 0;
my $seg = 100000;
my $add = $seg;
my $nsave = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if(! $line){
		next;
	}
	my @A = split(/\t/, $line);
	my $numA = @A;
	if($numA < 10){
		$header .= $line."\n";
	}
	elsif($line =~ /CHROM/){
		$header .= $line."\n";
	}
	else{
		$rvcf .= $line."\n";
		$cnt++;
		
		if($cnt == $seg){
			if($nsave == 0 && $nf == 0){
				$rvcf = $header.$rvcf;
				open(my $rfh, ">", $rfile_combined);
				print $rfh $rvcf;
				close $rfh;
			}
			else{
				open(my $rfh, ">>", $rfile_combined);
				print $rfh $rvcf;
				close $rfh;
			}
			$rvcf = "";
			$seg += $add;
			$nsave++;
		}
	}
}
close $fh;

if($nsave == 0 && $nf == 0){
	$rvcf = $header.$rvcf;
	open(my $rfh, ">", $rfile_combined);
	print $rfh $rvcf;
	close $rfh;
}
else{
	open(my $rfh, ">>", $rfile_combined);
	print $rfh $rvcf;
	close $rfh;
}

}


#--------------------------------------------------------------------------------
sub Batch_exec{
my ($cmd, $sh, $logfile, $sid, $t) = @_;

if($cmd && $sh && $logfile){
	print " thread [$t] : $sid ...\n";
	if(system("$cmd > $logfile 2>&1") != 0){
		print " thread [$t] : failed...\n";
	}
	else{
		print " thread [$t] : completed.\n";
	}
}

}


#-------------------------------------------------------------------------------
sub Open_fasta_as_hash{
my $file = shift;

print "! reading [$file]...\n";
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

print "! reference = [$file]\n";
print " - [$numID] sequence ID, [$total_len] bp\n";

return ($hash, \@SID);
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


#-------------------------------------------------------------------------------
sub Path2pref{
my $file = shift;

if($file && $file ne 'null'){
	my @tmp = split(/\//, $file);
	my $n = @tmp;
	$tmp[$n - 1] =~ s/\.gvcf\.gz//;
	$tmp[$n - 1] =~ s/\.gvcf//;
	$tmp[$n - 1] =~ s/\.vcf\.gz//;
	$tmp[$n - 1] =~ s/\.vcf//;
	return $tmp[$n - 1];
}

}


#---------------------------------------------------------------------------------------
sub Search_app2{
my ($hbin, $script_path, $xp_wrapper) = @_;

my @APP = keys(%{$hbin});
@APP = sort {$a cmp $b} @APP;
my $tmp = "_search_app_temp_";
my $tmp2 = "_vercheck_";

if($xp_wrapper == 0){
	print "! checking PATH of required programs...\n";
}

my $hash = {};
my $happ = {};
foreach my $app (@APP){
	my $Bin = $hbin->{$app};
	foreach my $bin (@{$Bin}){
		if(! $hash->{$bin} || ! -e $hash->{$bin}){
			system("which $bin > $tmp");
			system("find /usr/local/vg/bin/vg >> $tmp 2>&1");
			
			if(-e $tmp){
				open(my $fh, "<", $tmp);
				while(my $line = <$fh>){
					$line =~ s/\n//;
					$line =~ s/\r//;
					
					if($line && $line =~ /$bin/i){
						$hash->{$bin} = $line;
						last;
					}
				}
				close $fh;
				system("rm $tmp");
			}
		}
	}
	
	foreach my $bin (@{$Bin}){
		if($hash->{$bin} && -e $hash->{$bin}){
			if($xp_wrapper == 0){
				print " [$bin] = [$hash->{$bin}] <OK>\n";
			}
			$happ->{$app} = $hash->{$bin};
		}
		else{
			print " NOT found [$bin], please check PATH\n";
			$happ->{err} += 1;
		}
	}
}

return $happ;
}


