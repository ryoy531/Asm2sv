#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use FindBin;

my $random_IDlist = shift;
my $genelist = shift;
my $qfasta = shift;
my $dbgff = shift;
my $dbfasta = shift;
my $dbprotein = shift;
my $dbtranscript = shift;
my $dbcds = shift;
my $afile = shift;
my $rfile_final = shift;
my $n0 = shift;
my $n1 = shift;
my $tmpdir = shift;
my $log_orfsearch = shift;
my $script = shift;
my $combined_qprot = shift;
my $tmp_ralign = shift;
my $bin_findorf = shift;
my $i = shift;
my $combined_qgff = shift;
my $neighbor_tlen = shift;

if(! $random_IDlist || ! $genelist){
	goto END;
}
if(! $combined_qgff){
	$combined_qgff = $combined_qprot;
	$combined_qgff =~ s/\.fasta/\.gff/;
}
if(! defined $neighbor_tlen){
	$neighbor_tlen = 20000;
}

my $script_path = $FindBin::Bin;

my @APP;
push(@APP, "samtools");
push(@APP, "gffread");
push(@APP, "gth");
push(@APP, "blat");
push(@APP, "blat2hints");
push(@APP, "matcher");
push(@APP, "miniprot");

my $hbin = {};
foreach my $app (@APP){
	$hbin->{$app}[0] = $app;
}
my $hpath = {};
my $hbinpath = Search_app2($hbin, \@APP, $hpath, $script_path, "pipe", $genelist);
if($hbinpath->{err}){
	print "! abort script due to missing program(s)..., please confirm program is in PATH.\n";
	goto END;
}

my $qpref = $qfasta;
if($qpref =~ /\.fasta/){
	$qpref =~ s/\.fasta//g;
}
elsif($qpref =~ /\.fa/){
	$qpref =~ s/\.fa//g;
}

my $hID = Read_list($genelist);
my $hgffinfo = Read_gff($dbgff, $hID);
my $dbg2t = $hgffinfo->{g2t};
my $dbt2g = $hgffinfo->{t2g};
my $hlist = $hgffinfo->{hash};
my $randomID = Read_randomIDlist($random_IDlist);
my @G = split(/\n/, $randomID);

my $hqseq = Open_fasta_as_hash($qfasta);
my $hdbseq = Open_fasta_as_hash($dbfasta);
my $hdbprot = Open_fasta_as_hash($dbprotein);
my $hdbtranscript = Open_fasta_as_hash($dbtranscript);
my $hdbcds = Open_fasta_as_hash($dbcds);
my $hpav = Open_pav_results($afile);
my $hprevr = {};
if(-e $rfile_final){
	$hprevr = Check_prev_results($rfile_final);
}

my $ralign;
my $cmd;
my $cnt = 0;
my $cnt_skip = 0;
for(my $j = $n0; $j < $n1; $j++){
	if($G[$j]){
		my $gid = $G[$j];
		
		if($hpav->{$gid}{status} && $hpav->{$gid}{status} eq 'P'){
			my $dbsid = $hpav->{$gid}{dbsid};
			my $dbpos0 = $hpav->{$gid}{dbpos0};
			my $dbpos1 = $hpav->{$gid}{dbpos1};
			my $qsid = $hpav->{$gid}{qsid};
			my $qpos0 = $hpav->{$gid}{qpos0};
			my $qpos1 = $hpav->{$gid}{qpos1};
			
			my $tmp_dbprot = "$tmpdir/_dbprot_".$gid.".fasta";
			my $tmp_dbtranscript = "$tmpdir/_dbtranscript_".$gid.".fasta";
			my $tmp_dbcds = "$tmpdir/_dbcds_".$gid.".fasta";
			my $tmp_dbindex = "$tmpdir/_dbprot_".$gid.".fasta.protein";
			my $tmp_dbfasta = "$tmpdir/_dbseq_".$gid.".fasta";
			my $tmp_qfasta = "$tmpdir/_qseq_".$qpref."_".$gid.".fasta";
			my $tmp_qgff = "$tmpdir/_qprot_".$qpref."_".$gid.".gff";
			my $tmp_qprot = "$tmpdir/_qprot_".$qpref."_".$gid.".fasta";
			my $tmp_align = "$tmpdir/_align_".$gid.".txt";
			my $mk_process = "$tmpdir/_mk_".$gid.".txt";
			
			if($dbg2t->{$gid}){
				if($hprevr && $hprevr->{$gid}){
					$cnt_skip++;
				}
				else{
					my @dbTID = split(/\n/, $dbg2t->{$gid});
					unless(-e $tmp_dbprot){
						my $each_dbprot;
						foreach my $dbtid (@dbTID){
							if($hdbprot->{$dbtid}{seq}){
								$each_dbprot .= ">".$dbtid."\n".$hdbprot->{$dbtid}{seq}."\n";
							}
						}
						SAVE($tmp_dbprot, $each_dbprot);
					}
					unless(-e $tmp_dbtranscript){
						my $each_dbtranscript;
						foreach my $dbtid (@dbTID){
							if($hdbtranscript->{$dbtid}{seq}){
								$each_dbtranscript .= ">".$dbtid."\n".$hdbtranscript->{$dbtid}{seq}."\n";
							}
						}
						SAVE($tmp_dbtranscript, $each_dbtranscript);
					}
					unless(-e $tmp_dbcds){
						my $each_dbcds;
						foreach my $dbtid (@dbTID){
							if($hdbcds->{$dbtid}{seq}){
								$each_dbcds .= ">".$dbtid."\n".$hdbcds->{$dbtid}{seq}."\n";
							}
						}
						SAVE($tmp_dbcds, $each_dbcds);
					}
					
					my $qsid_seqlen = length($hqseq->{$qsid}{seq});
	#				my $qpos0_nb = $qpos0;
	#				my $qpos1_nb = $qpos1;
					my $qpos0_nb = $qpos0 - $neighbor_tlen;
					my $qpos1_nb = $qpos1 + $neighbor_tlen;
					if($qpos0_nb < 0){
						$qpos0_nb = 1;
					}
					if($qpos1_nb > $qsid_seqlen){
						$qpos1_nb = $qsid_seqlen;
					}
					
					my $sub_dbseq = Extract_and_save($tmp_dbfasta, $dbsid, $hdbseq->{$dbsid}{seq}, $dbpos0, $dbpos1);
					my $sub_qseq = Extract_and_save($tmp_qfasta, $qsid, $hqseq->{$qsid}{seq}, $qpos0_nb, $qpos1_nb);
					
					if($sub_dbseq && $sub_dbseq ne 'null' && $sub_qseq && $sub_qseq ne 'null'){
						my $len_sub_dbseq = length($sub_dbseq);
						my $len_sub_qseq = length($sub_qseq);
						
	#							$cmd .= "if test -e \"$mk_process\"; then\n";
	#							$cmd .= "\techo '$gid already processed' >> $log_orfsearch 2>&1\n";
	#							$cmd .= "fi\n";
						
	#							if($sub_qseq && -e $tmp_dbprot && -e $tmp_qfasta){
						if(-e $tmp_dbprot && -e $tmp_qfasta){
							if(-e $bin_findorf){
								$cmd .= "perl $bin_findorf $hbinpath->{gth} $hbinpath->{miniprot} $hbinpath->{samtools} $hbinpath->{gffread} $hbinpath->{blat} $hbinpath->{blat2hints} $hbinpath->{matcher} $gid $tmp_dbprot $tmp_dbtranscript $tmp_dbcds $tmp_qfasta $tmp_qgff $tmp_qprot $tmp_align $tmpdir $qpref $qpos0_nb $combined_qprot $combined_qgff >> $log_orfsearch 2>&1\n";
							}
							else{
								$cmd .= "if test ! -e \"$tmp_qgff\"; then\n";
								$cmd .= "\tif test -e \"$tmp_qfasta\"; then\n";
								$cmd .= "\t\t$hbinpath->{gth} -genomic $tmp_qfasta -protein $tmp_dbprot -gff3out -skipalignmentout -o $tmp_qgff >> $log_orfsearch 2>&1\n";
								$cmd .= "\t\t$hbinpath->{samtools} faidx $tmp_qfasta >> $log_orfsearch 2>&1\n";
								$cmd .= "\t\t$hbinpath->{gffread} $tmp_qgff -g $tmp_qfasta -y $tmp_qprot >> $log_orfsearch 2>&1\n";
								$cmd .= "\t\t$hbinpath->{matcher} -asequence $tmp_dbprot -bsequence $tmp_qprot -outfile $tmp_align >> $log_orfsearch 2>&1\n";
								
								$cmd .= "\t\tif test -e \"$tmp_dbprot\"; then\n";
								$cmd .= "\t\t\trm $tmp_dbprot\n";
								$cmd .= "\t\tfi\n";
								$cmd .= "\t\tif test -e \"$tmp_dbprot.md5\"; then\n";
								$cmd .= "\t\t\trm $tmp_dbprot.md5\n";
								$cmd .= "\t\tfi\n";
								$cmd .= "\t\tif test -e \"$tmp_dbprot.suf\"; then\n";
								$cmd .= "\t\t\trm $tmp_dbprot.\*\n";
								$cmd .= "\t\tfi\n";
								$cmd .= "\t\tif test -e \"$tmp_dbindex.bck\"; then\n";
								$cmd .= "\t\t\trm $tmp_dbindex.\*\n";
								$cmd .= "\t\tfi\n";
								$cmd .= "\t\tif test -e \"$tmp_qfasta\"; then\n";
								$cmd .= "\t\t\trm $tmp_qfasta\n";
								$cmd .= "\t\tfi\n";
								$cmd .= "\t\tif test -e \"$tmp_qfasta.fai\"; then\n";
								$cmd .= "\t\t\trm $tmp_qfasta.\*\n";
								$cmd .= "\t\tfi\n";
								$cmd .= "\t\tif test -e \"$tmp_qgff\"; then\n";
								$cmd .= "\t\t\trm $tmp_qgff\n";
								$cmd .= "\t\tfi\n";
								$cmd .= "\tfi\n";
								$cmd .= "fi\n";
							}
						}
						
						# time-consuming...
	#							my $len_NNN = Count_NNN($sub_qseq);
	#							$ralign .= $gid."\t".$tmp_align."\t".$len_NNN."\t".$tmp_qprot."\n";
						$ralign .= $gid."\t".$tmp_align."\t-\t".$tmp_qprot."\t-\t".$len_sub_dbseq."\t".$len_sub_qseq."\n";
						
	#							push(@gabbage, $tmp_dbprot);
	#							push(@gabbage, $tmp_qfasta);
	#							push(@gabbage, $tmp_qgff);
	#							push(@gabbage, $tmp_qprot);
	#							push(@gabbage, $tmp_align);
					}
				}
			}
		}
		$cnt++;
	}
}

print "! [$cnt] queries, [$cnt_skip] skipped\n";

if($ralign){
	ADD($tmp_ralign, $ralign);
}
if($cmd){
	SAVE($script, $cmd);
}

END:{
	print "\n";
}


######################################################################
#-------------------------------------------------------------------------------
sub Extract_and_save{
my $rfile = shift;
my $id = shift;
my $seq = shift;
my $p0 = shift;
my $p1 = shift;

if($seq){
	my $subseq = substr($seq, $p0 - 1, abs($p1 - $p0));
	my $subfasta = ">".$id."\n".$subseq;
	
	unless(-e $rfile){
		open(my $rfh, ">", $rfile) or return 'null';
		print $rfh $subfasta;
		close $rfh;
	}
	return $subseq;
}

}


#---------------------------------------------------------------------
sub SAVE{
my ($file, $str, $header) = @_;

if($str && $header){
	$str = $header.$str;
	open(my $rfh, ">", $file) or die;
	print $rfh $str;
	close $rfh;
	print "! output [$file]\n";
}
elsif($str){
	open(my $rfh, ">", $file) or die;
	print $rfh $str;
	close $rfh;
}

}


#---------------------------------------------------------------------
sub ADD{
my ($file, $str) = @_;

if($str){
	open(my $rfh, ">>", $file);
	print $rfh $str;
	close $rfh;
}

}


#-------------------------------------------------------------------------------
sub Read_randomIDlist{
my $file = shift;

open(my $fh, "<", $file) or die;
my $ID = "";
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	unless($line){
		next;
	}
	$ID .= $line."\n";
}
close $fh;

return $ID;
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
	elsif($line =~ /gene_id/ || $line =~ /\#/){
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


#-------------------------------------------------------------------------------
sub Check_prev_results{
my $file = shift;

open(my $fh, "<", $file) or die;
my $hash = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /directory/ && $line =~ /dbfasta/ && $line =~ /seqid/){
		next;
	}
	
	my @A = split(/\t/, $line);
	my $numA = @A;
	if($numA > 18){
		$hash->{$A[17]} = $line;
	}
}
close $fh;

return $hash;
}


#-------------------------------------------------------------------------------
sub Open_pav_results{
my $file = shift;

print "! reading PAV result [$file] ...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
my $cnt = 0;
my $cntP = 0;
my $cntA = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /directory/ && $line =~ /dbfasta/ && $line =~ /seqid/){
		next;
	}
	
	my @A = split(/\t/, $line);
	if(defined $A[18]){
		my $gid = $A[23];
		my $pav_status = $A[22];
		
		if($pav_status && $pav_status =~ /present/i){
			$hash->{$gid}{status} = "P";
			$hash->{$gid}{dbsid} = $A[2];
			$hash->{$gid}{dbpos0} = $A[5];
			$hash->{$gid}{dbpos1} = $A[6];
			$hash->{$gid}{qsid} = $A[8];
			$hash->{$gid}{qpos0} = $A[12];
			$hash->{$gid}{qpos1} = $A[13];
			$hash->{$gid}{data} = $line;
			$cntP++;
		}
		elsif($pav_status && $pav_status =~ /insert/i){
			$hash->{$gid}{status} = "P";
			$hash->{$gid}{dbsid} = $A[2];
			$hash->{$gid}{dbpos0} = $A[5];
			$hash->{$gid}{dbpos1} = $A[6];
			$hash->{$gid}{qsid} = $A[8];
			$hash->{$gid}{qpos0} = $A[12];
			$hash->{$gid}{qpos1} = $A[13];
			$hash->{$gid}{data} = $line;
			$cntP++;
		}
		else{
			$hash->{$gid}{status} = "A";
			$hash->{$gid}{data} = $line;
			$cntA++;
		}
		
		$cnt++;
	}
}
close $fh;

print "! total [$cnt] lines, P=[$cntP], A=[$cntA]\n";

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
			system("which $bin > $tmp");
			system("which $bin.pl >> $tmp");
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

