#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use threads;
use FindBin;
use Getopt::Long;
use Sys::Hostname;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "splitseq_run_cactus.pl, version 1.01\n";
$version .= "last update: [2024\/4\/29]\n";

print "$version\n";

#-------------------------------------------------------------------------------

my $host = hostname();
my $wpath = getcwd();
my $script_path = $FindBin::Bin;
my $status = "null";

my $hbin = {};
$hbin->{"cactus-pangenome"}[0] = "cactus-pangenome";

my $hbinpath = Search_app2($hbin, $script_path, 0);

if($hbinpath->{err}){
	print "! abort script due to missing program(s) or version inconsistency...\n";
	print "! may need [source /usr/local/cactus/cactus_env/bin/activate] to use specific python environment\n\n";
	goto END;
}

my $datalist;
my $chrinfo;
my $rfile;
my $cpu;
my $npart;
my $help;
my $rdir;

GetOptions('--list=s' => \$datalist, '--thread=i' => \$cpu, '--partition=i' => \$npart, '--chrinfo=s' => \$chrinfo, '--output=s' => \$rdir, '--resultVcf=s' => \$rfile, '--help' => \$help);

print "\n";
if(! $datalist || ! -e $datalist){
	print "! missing --list [datalist csv]\n";
	$status = "missing";
}
if(! $chrinfo || ! -e $chrinfo){
	print "! missing --chrinfo [csv file]\n";
	$status = "missing";
}
if($status eq 'missing'){
	Help();
	goto END;
}
if(! $cpu || $cpu =~ /\D/){
	$cpu = 16;
}
if(! $npart || $npart =~ /\D/){
	$npart = 1;
}
if(! $rdir){
	$rdir = "split_run_cactuspg";
}
unless(-e $rdir){
	system("mkdir $rdir");
}

my $htmp_chrinfo = Read_seqinfo($chrinfo);
my $hseqinfo = $htmp_chrinfo->{seqID};
my $hnrCID = $htmp_chrinfo->{nrCID};
my @nrCID = keys(%{$hnrCID});
@nrCID = sort {$a cmp $b} @nrCID;
my $num_CID = @nrCID;

my ($AoH, $err_list, $setref) = Read_datalist($datalist, "check");

print "\n! collecting seq entries...\n";
my $hdata = {};
foreach my $hash (@{$AoH}){
	$hdata = Prep_fasta_forPG($hash->{qfasta}, $hash->{custom_prefix}, $hseqinfo, \@nrCID, $rdir, $hdata);
}
chdir $rdir;

print "\n! run cactus-pangenome in each entry...\n";
my $mcpu = $cpu;
if($npart > 1){
	$mcpu = int($cpu / $npart);
	if($mcpu == 0){
		$mcpu = 1;
	}
	print "! npartition = [$npart] | will use [$mcpu] thread in each partition\n";
}

my $d = int($num_CID / $npart) + 1;
my @VF;
my $err = 0;
for(my $i = 0; $i < $d; $i++){
	my $p0 = $i * $npart;
	my $p1 = ($i + 1) * $npart;
	
	my $thrs = [];
	for(my $j = $p0; $j < $p1; $j++){
		if(! $nrCID[$j] || ! $hdata->{$nrCID[$j]}{entry}){
			next;
		}
		
		my $cid = $nrCID[$j];
		my $seqlist = "./data_".$cid."/list.txt";
		SAVE($seqlist, $hdata->{$cid}{entry});
		
		my $adir = "data_".$cid;
		my $js = "./data_".$cid."/js";
		my $vcf = "./data_".$cid."/data_".$cid.".vcf.gz";
		my $alog = "data_".$cid."/log.$cid.txt";
		my $sh = "data_".$cid."/cmd.$cid.sh";
		push(@VF, $vcf);
		
		if(-e $vcf){
			print " - [$cid] | [$vcf] already exists, skip...\n";
		}
		else{
			if(-e $js){
				system("rm -R $js");
			}
			
			my $cmd = "cd $wpath/$rdir\n";
			$cmd .= $hbinpath->{"cactus-pangenome"}." $js $seqlist --outDir $adir --outName $adir --reference $setref --vcf --giraffe --gfa --gbz --mgCores $mcpu --maxCores $mcpu > $alog 2>&1";
			SAVE($sh, $cmd);
			
			my $thr = threads->new(\&Batch, $sh, $cid, $alog);
			push(@{$thrs}, $thr);
		}
	}
	if(@{$thrs}){
		foreach my $thr (@{$thrs}){
			$err += $thr->join();
		}
	}
}
if($err > 0){
	print "! process failed, abort script...\n";
	goto END;
}
else{
	print "\n! combining VCF of all entries...\n";
	my $AoH = [];
	for(my $i = 0; $i < $d; $i++){
		my $p0 = $i * $npart;
		my $p1 = ($i + 1) * $npart;
		
		my $thrs = [];
		for(my $j = $p0; $j < $p1; $j++){
			if(! $VF[$j]){
				next;
			}
			
			my $thr = threads->new(\&ReadVcf, $VF[$j], 50);
			push(@{$thrs}, $thr);
		}
		if(@{$thrs}){
			foreach my $thr (@{$thrs}){
				my $rh = $thr->join();
				push(@{$AoH}, $rh);
			}
		}
	}
	
	my $rheader0;
	my $rheader1;
	my $rline;
	my $numF = @VF;
	my $cnt = 0;
	for(my $i = 0; $i < $numF; $i++){
		my $file = $VF[$i];
		my $hash = $AoH->[$i];
		if($i == 0){
			$rheader0 = $hash->{$file}{header};
			$rheader0 .= $hash->{$file}{header_contig};
			$rheader1 = $hash->{$file}{header_CHROM};
		}
		else{
			$rheader0 .= $hash->{$file}{header_contig};
		}
		$cnt += $hash->{$file}{cnt};
		$rline .= $hash->{$file}{data};
	}
	print "! total [$cnt] variants\n";
	
	if(! $rfile){
		$rfile = "cactus_pangenome.vcf";
	}
	
	my $rfinal = $rheader0.$rheader1.$rline;
	SAVE($rfile, $rfinal);
	print "! output result = [$rfile]\n";
}


END:{
	print "\n! End of script.\n\n";
}


########################################################################################
#---------------------------------------------------------------------------------------
sub Help{
	my $help_command =<<"EOS";
Usage: 
  /path/to/splitseq_run_cactus.pl -l [datalist csv] -c [chr alias info] -t [CPU thread]

--list or -l          : data file list (PATH of fasta, GFF, psl, and other files)
                       *header should be "dbfasta" "dbgff" "genelist" "qfasta" "qpsl" "chralias"
                       
--chralias or -c      : chromosome ID alias info (.tsv, recommended)
--resultVcf or -r     : name of result VCF file
--thread or -t        : CPU thread number (default 16)
--partition or -p     : num of parallel job (default 1)
EOS
	print "\n";
}


#---------------------------------------------------------------------------------------
sub ReadVcf{
my $file = shift;
my $thlen = shift;

my $hash = {};
my $cnt = 0;
if(-e $file){
	print " reading [$file]...\n";
	open(my $fh, "gzip -dc $file |");
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line =~ /\#/){
			if($line =~ /contig\=/){
				$hash->{$file}{header_contig} .= $line."\n";
			}
			elsif($line =~ /CHROM/){
				$hash->{$file}{header_CHROM} .= $line."\n";
			}
			else{
				$hash->{$file}{header} .= $line."\n";
			}
		}
		else{
			my @A = split(/\t/, $line);
			if($A[0] && $A[1] && $A[3] && $A[4]){
				my $lenR = length($A[3]);
				my $sw = 0;
				my @tmpA = split(/,/, $A[4]);
				foreach my $eachA (@tmpA){
					my $lenA = length($eachA);
					if( abs($lenR - $lenA) >= $thlen ){
						$sw = 1;
					}
				}
				
				if($sw eq '1'){
					my $varid = $A[0].$A[1].$A[3].$A[4];
					if(! $hash->{$file}{nr}{$varid}){
						$hash->{$file}{nr}{$varid} = 1;
						$hash->{$file}{data} .= $line."\n";
						$cnt++;
					}
				}
			}
			else{
				print "$line\n";
				last;
			}
		}
	}
	close $fh;
	
	print "! [$file] [$cnt] nrlines\n";
}

$hash->{$file}{cnt} = $cnt;

return $hash;
}


#---------------------------------------------------------------------------------------
sub Batch{
my ($sh, $cid, $alog) = @_;

print " - running [$cid]...\n";
my $err = 0;
if(system("bash $sh") != 0){
	print "! [$cid] failed. see log = [$alog]\n";
	$err++;
}
else{
	print "! [$cid] completed.\n";
}

return $err;
}


#-------------------------------------------------------------------------------
sub Prep_fasta_forPG{
my $file = shift;
my $prefix = shift;
my $hseqinfo = shift;
my $CID = shift;
my $rdir = shift;
my $hdata = shift;

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
				print "! duplicated seq ID : [$ID]\n";
			}
			if($hseqinfo->{$prefix}{$ID}){
				my $cid = $hseqinfo->{$prefix}{$ID};
				my $len = length($seq);
				$hash->{$cid} = $seq;
				$total_len += $len;
				push(@SID, $ID);
			}
			$seq = "";
		}
		
		$ID = $line;
		$ID =~ s/\>//;
		my @tmp = split(/\s/, $ID);
		$ID = $tmp[0];
		$seq = "";
		$cnt++;
	}
	else{
#		$line =~ s/\.//g;
#		if($line){
			$seq .= $line;
#		}
	}
}
close $fh;

if($seq){
	if($hash->{$ID}){
		print "! duplicated seq ID : [$ID]\n";
	}
	if($hseqinfo->{$prefix}{$ID}){
		my $cid = $hseqinfo->{$prefix}{$ID};
		my $len = length($seq);
		$hash->{$cid} = $seq;
		$total_len += $len;
		push(@SID, $ID);
	}
}

my $numID = @SID;

my $err = 0;
foreach my $cid (@{$CID}){
	if($hash->{$cid}){
		my $tdir = $rdir."/data_".$cid;
		unless(-e $tdir){
			system("mkdir $tdir");
		}
		my $tdir2 = $rdir."/data_".$cid."/fasta";
		unless(-e $tdir2){
			system("mkdir $tdir2");
		}
		my $path = $tdir2."/".$prefix."_".$cid.".fasta";
		my $fasta = ">".$prefix."_".$cid."\n".$hash->{$cid};
		
		open(my $rfh, ">", $path);
		print $rfh $fasta;
		close $rfh;
		
		$hdata->{$cid}{entry} .= $prefix."\t"."./data_".$cid."/fasta/".$prefix."_".$cid.".fasta\n";
	}
	else{
		print " - $prefix : unable to find seq for [$cid]...\n";
		$err++;
	}
}

if($err == 0){
	print " - $prefix = [$numID] seq\n";
}

$hdata->{err} += $err;

return $hdata;
}


#-------------------------------------------------------------------------------
sub Read_datalist{
my $file = shift;
my $option = shift;
my $gid_alias = shift;

if(! $gid_alias || ! -e $gid_alias){
	$gid_alias = "null";
}

print "! reading data list [$file]...\n";
open(my $fh, "<", $file);
my $AoA = [];
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
					if(! -e $A[$i]){
						print " [$cnt] : missing [$A[$i]]\n";
						$missing++;
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
				
				if($option eq 'check'){
					if(! $hash->{custom_prefix}){
						$hash->{custom_prefix} = $qpref;
					}
					if($dbpref eq $qpref){
						print " [$hash->{custom_prefix}] : q=[$hash->{qfasta}] <OK> *reference\n";
						$setref = $hash->{custom_prefix};
					}
					else{
						print " [$hash->{custom_prefix}] : q=[$hash->{qfasta}] <OK>\n";
					}
					push(@{$AoA}, $hash);
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

return ($AoA, $err, $setref);
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
	
	if(! $line || $line =~ /\#/){
		next;
	}
	elsif($line =~ /prefix/ || $line =~ /norder/){
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
	$hash->{nrCID}{$A[1]} = 1;
	
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
			if(-e $tmp){
				open(my $fh, "<", $tmp);
				while(my $line = <$fh>){
					$line =~ s/\n//;
					$line =~ s/\r//;
					
					if($bin eq 'cactus-pangenome'){
						if($line && $line =~ /cactus-pangenome/){
							my $cmdtest = "$line --help > $tmp2 2>&1";
							my $judge = "false";
							if(system($cmdtest) == 0){
								open(my $tfh, "<", $tmp2);
								while(my $line = <$tfh>){
									$line =~ s/\n//;
									$line =~ s/\r//;
									
									if($line){
										if($line =~ /usage/){
											$judge = "true";
										}
										if($line =~ /nomodule/i){
											print "$line\n";
											$judge = "false";
										}
									}
								}
								close $tfh;
								system("rm $tmp2");
							}
							if($judge eq 'true'){
								$hash->{$bin} = $line;
								last;
							}
						}
					}
					else{
						if($line && $line =~ /$bin/i){
							$hash->{$bin} = $line;
							last;
						}
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


#---------------------------------------------------------------------------------------
sub SAVE{
my ($file, $str) = @_;

if($str){
	open(my $fh, ">", $file);
	print $fh $str;
	close $fh;
}

}

