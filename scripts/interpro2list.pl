#!/usr/bin/perl -w
use strict;
use warnings;
use threads;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "Iprscan_purse.pl version 2.01\n";
$version .= "last update: [2021\/6\/9]\n";
$version .= "copyright: ryoichi yano\n";

#print "$version";
#-------------------------------------------------------------------------------

my $gff = shift;
my $dfile = shift;
my $goobo = shift;
my $IDfile = shift;
my $gfile0 = shift;
my $gfile1 = shift;

my $pth = 1;
my $q = "y";

my $err = 0;
if(! $gff || ! -e $gff){
	print "! missing Gff...\n";
	$err++;
}
if(! $dfile || ! -e $dfile){
	print "! missing InterProscan result...\n";
	$err++;
}
if(! $goobo || ! -e $goobo){
	print "! missing go.obo ...\n";
	$err++;
}
if(! $IDfile || ! -e $IDfile){
	$IDfile = "null";
}
if(! $gfile0){
	$gfile0 = "GOterm_list.csv";
}
if(! $gfile1){
	$gfile1 = "GOnum_list.csv";
}

if($err > 0){
	print "\n! usage : interpro2list.pl [gff] [InterProscan5 result] [/path/tp/go.obo.gz]\n";
	goto END;
}

if($q =~ /y/i){
	my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime;
	$year += 1900;
	$mon = Num2Month($mon);
	if($min < 10){
		$min = "0".$min;
	}
	if($sec < 10){
		$sec = "0".$sec;
	}
	my $timestamp = $year.$mon.$mday;
	
	my $hgff = GFF_to_hash($gff);
	my $t2g = $hgff->{t2g};
	my $g2t = $hgff->{g2t};
	
	my $IDs = [];
	@{$IDs} = keys(%{$t2g});
	@{$IDs} = sort {$a cmp $b} @{$IDs};
	
	my $hgoobo = GetGO($goobo);		#all GO category info + string
	my $dhash = Open_data_as_hash($dfile, $pth, $hgoobo);
	my $hash_all = $dhash->{all};
	my $hash_sep = $dhash->{sep};
	my $ps = $dhash->{program};
	my $GOn = $dhash->{GOnum};
	my $IPn = $dhash->{IPnum};
	my $hash_gocnt = {};
	
	my $reph = {};
	if($IDfile ne 'null' && -e $IDfile){
		$reph = Open_replaceID_info($IDfile);
	}
	
	my $GOnum_list = "ID,GO+IPR,IPR,Biological Process,Molecular Function,Cellular Component,IPR ID\n";
	my $GOinfo_list = $GOnum_list;
	my $GOsql_list = "ID\tGOIPR\tBiological Process\tMolecular Function\tCellular Component\tIPR ID\n";		#for mysql DB in Melonet-DB
	
	my $GOnum_list_repID = "ID,GO+IPR,IPR,Biological Process,Molecular Function,Cellular Component,IPR ID\n";
	my $GOinfo_list_repID = $GOnum_list;
	my $GOsql_list_repID = "ID\tGOIPR\tBiological Process\tMolecular Function\tCellular Component\tIPR ID\n";		#for mysql DB in Melonet-DB
	
	my $cnt_repID = 0;
	my $GOjudge_true = 0;
	my $IPRjudge_true = 0;
	my $gbdata = {};			#gene ID based data (not transcript_variant based)
	my $nrgid = {};				#non-redundant gene ID
	foreach my $id (@{$IDs}){
		#for gene_based id---------------------//
		unless($t2g->{$id}){
			print "! missing gid for [$id]\n";
			sleep(1);
			next;
		}
		
		my $gid = $t2g->{$id};
		
		my $rep_gid= $gid;
		my $rep_tid= $id;
		my $repID = "FALSE";
		if($reph->{$gid}){
			$rep_gid = $reph->{$gid};
			
			$rep_tid =~ s/$gid/$rep_gid/;
			$repID = "TRUE";
			
#			if($rep_tid eq $id || $rep_gid eq $gid){
#				print "$rep_tid $id | $rep_gid $gid\n";
#				sleep(1);
#				$repID = "FALSE";
#			}
			if($repID eq 'TRUE'){
				$cnt_repID++;
			}
		}
		
		my $GOjudge = "true";
		my $IPRjudge = "true";
		my $bpstr = "-";
		my $mfstr = "-";
		my $clstr = "-";
		my $IPstr = "-";
		if($dhash->{GO_bp}{$id}){
			$bpstr = $dhash->{GO_bp}{$id};
			$bpstr = NRGO($bpstr);
			
			my @each_GOn = split(/;/, $bpstr);
			foreach my $go (@each_GOn){
				if($go){
					$hash_gocnt->{$go} += 1;
				}
			}
			
			$gbdata->{$gid}{bp} .= $bpstr."\n";
		}
		elsif($dhash->{GO_bp}{$rep_tid}){
			$bpstr = $dhash->{GO_bp}{$rep_tid};
			$bpstr = NRGO($bpstr);
			
			my @each_GOn = split(/;/, $bpstr);
			foreach my $go (@each_GOn){
				if($go){
					$hash_gocnt->{$go} += 1;
				}
			}
			
			$gbdata->{$gid}{bp} .= $bpstr."\n";
		}
		if($dhash->{GO_mf}{$id}){
			$mfstr = $dhash->{GO_mf}{$id};
			$mfstr = NRGO($mfstr);
			
			my @each_GOn = split(/;/, $mfstr);
			foreach my $go (@each_GOn){
				if($go){
					$hash_gocnt->{$go} += 1;
				}
			}
			
			$gbdata->{$gid}{mf} .= $mfstr."\n";
		}
		elsif($dhash->{GO_mf}{$rep_tid}){
			$mfstr = $dhash->{GO_mf}{$rep_tid};
			$mfstr = NRGO($mfstr);
			
			my @each_GOn = split(/;/, $mfstr);
			foreach my $go (@each_GOn){
				if($go){
					$hash_gocnt->{$go} += 1;
				}
			}
			
			$gbdata->{$gid}{mf} .= $mfstr."\n";
		}
		if($dhash->{GO_cl}{$id}){
			$clstr = $dhash->{GO_cl}{$id};
			$clstr = NRGO($clstr);
			
			my @each_GOn = split(/;/, $clstr);
			foreach my $go (@each_GOn){
				if($go){
					$hash_gocnt->{$go} += 1;
				}
			}
			
			$gbdata->{$gid}{cl} .= $clstr."\n";
		}
		elsif($dhash->{GO_cl}{$rep_tid}){
			$clstr = $dhash->{GO_cl}{$rep_tid};
			$clstr = NRGO($clstr);
			
			my @each_GOn = split(/;/, $clstr);
			foreach my $go (@each_GOn){
				if($go){
					$hash_gocnt->{$go} += 1;
				}
			}
			
			$gbdata->{$gid}{cl} .= $clstr."\n";
		}
		if($dhash->{GO_ipr}{$id}){
			$IPstr = $dhash->{GO_ipr}{$id};
			$IPstr = NRGO($IPstr);
			
			my @each_GOn = split(/;/, $IPstr);
			foreach my $go (@each_GOn){
				if($go){
					$hash_gocnt->{$go} += 1;
				}
			}
			
			$gbdata->{$gid}{IP} .= $IPstr."\n";
		}
		elsif($dhash->{GO_ipr}{$rep_tid}){
			$IPstr = $dhash->{GO_ipr}{$rep_tid};
			$IPstr = NRGO($IPstr);
			
			my @each_GOn = split(/;/, $IPstr);
			foreach my $go (@each_GOn){
				if($go){
					$hash_gocnt->{$go} += 1;
				}
			}
			
			$gbdata->{$gid}{IP} .= $IPstr."\n";
		}
		
		if($bpstr eq '-' && $mfstr eq '-' && $clstr eq '-'){
			$GOjudge = "false";
		}
		if($IPstr eq '-'){
			$IPRjudge = "false";
		}
		
		my $GOIPRjudge = "true";
		if($bpstr eq '-' && $mfstr eq '-' && $clstr eq '-' && $IPstr eq '-'){
			$GOIPRjudge = "false";
		}
		
		if($GOjudge eq 'true' && $hash_all->{$id}){
			$gbdata->{$gid}{line} .= $hash_all->{$id};
		}
		
		$nrgid->{$gid} = $GOjudge;
		$gbdata->{$gid}{cnt} += 1;
		
		$GOnum_list .= $gid.",".$id.",".$GOIPRjudge.",".$IPRjudge.",".$bpstr.",".$mfstr.",".$clstr.",".$IPstr."\n";
		
		my $bpstr2 = "-";
		my $mfstr2 = "-";
		my $clstr2 = "-";
		my $IPstr2 = "-";
		if($dhash->{GO_bp2}{$id}){
			$bpstr2 = $dhash->{GO_bp2}{$id};
			$bpstr2 = NRGO2($bpstr2);
		}
		elsif($dhash->{GO_bp2}{$rep_tid}){
			$bpstr2 = $dhash->{GO_bp2}{$rep_tid};
			$bpstr2 = NRGO2($bpstr2);
		}
		if($dhash->{GO_mf2}{$id}){
			$mfstr2 = $dhash->{GO_mf2}{$id};
			$mfstr2 = NRGO2($mfstr2);
		}
		elsif($dhash->{GO_mf2}{$rep_tid}){
			$mfstr2 = $dhash->{GO_mf2}{$rep_tid};
			$mfstr2 = NRGO2($mfstr2);
		}
		if($dhash->{GO_cl2}{$id}){
			$clstr2 = $dhash->{GO_cl2}{$id};
			$clstr2 = NRGO2($clstr2);
		}
		elsif($dhash->{GO_cl2}{$rep_tid}){
			$clstr2 = $dhash->{GO_cl2}{$rep_tid};
			$clstr2 = NRGO2($clstr2);
		}
		if($dhash->{GO_ipr2}{$id}){
			$IPstr2 = $dhash->{GO_ipr2}{$id};
			$IPstr2 = NRGO2($IPstr2);
		}
		elsif($dhash->{GO_ipr2}{$rep_tid}){
			$IPstr2 = $dhash->{GO_ipr2}{$rep_tid};
			$IPstr2 = NRGO2($IPstr2);
		}
		
		my $bpstr_sql = $bpstr2;
		my $mfstr_sql = $mfstr2;
		my $clstr_sql = $clstr2;
		my $IPstr_sql = $IPstr;
		$bpstr_sql =~ s/\;/<sqlsep>/g;
		$mfstr_sql =~ s/\;/<sqlsep>/g;
		$clstr_sql =~ s/\;/<sqlsep>/g;
		$IPstr_sql =~ s/\;/<sqlsep>/g;
		
		$GOinfo_list .= $gid.",".$id.",".$GOIPRjudge.",".$IPRjudge.",".$bpstr2.",".$mfstr2.",".$clstr2.",".$IPstr2."\n";
		
		if($GOjudge eq 'true'){
			$GOjudge_true++;
		}
		if($IPRjudge eq 'true'){
			$IPRjudge_true++;
		}
	}
	
	my $num_total = @{$IDs};
	print "! [$num_total] entries. Of these,\n";
	print "! [$GOjudge_true] entries with GO annotation.\n";
	
	if($cnt_repID > 0){
		print "! [$cnt_repID] entries replaced ID according to [$IDfile].\n";
	}
	
	$GOnum_list = Sort_lines_by_columm($GOnum_list, 0, "csv");
	$GOinfo_list = Sort_lines_by_columm($GOinfo_list, 0, "csv");
	
	my @DD = split(/\//, $dfile);
	my $numDD = @DD;
	$dfile = $DD[$numDD - 1];
	
	SAVE($gfile0, $GOinfo_list);
	print "\n! output [$gfile0]\n";
	
	SAVE($gfile1, $GOnum_list);
	print "! output [$gfile1]\n";
}

END:{
	print "\n";
}


################################################################################
#-------------------------------------------------------------------------------
sub GFF_to_hash{
my $file = shift;

print "\n! reading [$file] to make G2T/T2G hash data...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
#	if($line =~ /Glyma/ && $line =~ /Wm82/){
#		$line =~ s/\.Wm82\.a2\.v1//g;
#		$line =~ s/\.Wm82\.a4\.v1//g;
#	}
	
	unless($line){
		next;
	}
#	elsif($line =~ /\#/){
#		next;
#	}
	
	my @A = split(/\t/, $line);
	
	if($A[0] =~ /\#/){
		next;
	}
	
	if($A[2] eq 'mRNA' || $A[2] eq 'transcript'){
		my @T = split(/\;/, $A[8]);
		my $gid;
		my $tid;
		foreach my $tag (@T){
			if($tag =~ /ID\=/){
				$tag =~ s/ID\=//;
				unless($tid){
					$tid = $tag;
				}
			}
			elsif($tag =~ /Parent\=/){
				$tag =~ s/Parent\=//;
				unless($gid){
					$gid = $tag;
				}
			}
			
			if($gid && $tid){
				last;
			}
		}
		
		if($gid && $tid){
			$hash->{g2t}{$gid} = $tid;
			$hash->{t2g}{$tid} = $gid;
			$cnt++;
		}
		elsif($tid){
			$gid = $tid;
			$hash->{g2t}{$gid} = $tid;
			$hash->{t2g}{$tid} = $gid;
			$cnt++;
		}
	}
}
close $fh;

print "! [$cnt] gene ID\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Sort_lines_by_columm{
my $lines = shift;
my $i = shift;
my $type = shift;

my @L = split(/\n/, $lines);

if(@L){
	my $header = shift(@L);

	my $AoA = [];
	foreach my $l (@L){
		my @A;
		if($type eq 'csv'){
			@A = split(/,/, $l);
		}
		elsif($type eq 'tsv'){
			@A = split(/\t/, $l);
		}
	#	print "$A[$i] | $l\n";
	#	sleep(1);
		
		push(@{$AoA}, \@A);
	}

	@{$AoA} = sort {$a->[$i] cmp $b->[$i]} @{$AoA};

	my $nl = $header."\n";
	foreach my $A (@{$AoA}){
		if($type eq 'csv'){
			$nl .= join(",", @{$A})."\n";
		}
		elsif($type eq 'tsv'){
			$nl .= join("\t", @{$A})."\n";
		}
	}

	return $nl;
}

}


#-------------------------------------------------------------------------------
sub Num2Month{
my $mon = shift;

if($mon == 0){$mon = "Jan";}
elsif($mon == 1){$mon = "Feb";}
elsif($mon == 2){$mon = "Mar";}
elsif($mon == 3){$mon = "Apr";}
elsif($mon == 4){$mon = "May";}
elsif($mon == 5){$mon = "Jun";}
elsif($mon == 6){$mon = "Jul";}
elsif($mon == 7){$mon = "Aug";}
elsif($mon == 8){$mon = "Sep";}
elsif($mon == 9){$mon = "Oct";}
elsif($mon == 10){$mon = "Nov";}
elsif($mon == 11){$mon = "Dec";}

return $mon;
}


#-------------------------------------------------------------------------------
sub NRGO{
my $str = shift;

my @tmp = split(/;/, $str);
my $hash = {};
foreach my $i (@tmp){
	$hash->{$i} = 1;
}

my @tmp2 = keys(%{$hash});
$str = join(";", @tmp2);

return $str;
}


#-------------------------------------------------------------------------------
sub NRGO2{
my $str = shift;

my @tmp = split(/;/, $str);
my $hash = {};
foreach my $i (@tmp){
	$hash->{$i} = 1;
}

my @tmp2 = keys(%{$hash});
$str = join(";", @tmp2);

return $str;
}


#-------------------------------------------------------------------------------
sub NRGO3{
my $str = shift;

my @tmp = split(/\n/, $str);
my $hash = {};
foreach my $subl (@tmp){
	my @tmp2 = split(/;/, $subl);
	foreach my $i (@tmp2){
		$hash->{$i} = 1;
	}
}

my @tmp2 = keys(%{$hash});
$str = join("; ", @tmp2);

return $str;
}


#-------------------------------------------------------------------------------
sub SAVE{
my $file = shift;
my $str = shift;

if($str){
	open(my $fh, ">", $file) or print "! unable to save [$file]\n";
	print $fh $str;
	close $fh;
}

}


#-------------------------------------------------------------------------------
sub Open_replaceID_info{
my $file = shift;

my $hash = {};
if($file && -e $file){
	print "! reading replace ID data [$file] ...\n";
	open(my $fh, "<", $file) or die;
	my $cnslist = 0;
	my $cnt = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		if($cnt == 0 && $line =~ /qID_consensus,qID_original/){
			$cnslist = 1;
		}
		
		if($line =~ /\t/){
			my @A = split(/\t/, $line);
			if($A[1]){
				$hash->{$A[0]} = $A[1];
			}
		}
		elsif($line =~ /,/){
			my @A = split(/,/, $line);
			if($cnslist == 0){
				if($A[1]){
					$hash->{$A[0]} = $A[1];
				}
			}
			elsif($cnslist == 1){
				if($A[3] && $A[4] && $A[3] ne '-' && $A[4] ne '-'){
					$hash->{$A[3]} = $A[4];
					$hash->{$A[4]} = $A[3];
				}
			}
		}
		$cnt++;
	}
	close $fh;
}

return $hash;
}


#-------------------------------------------------------------------------------
sub GetGO{
my $file = shift;

open(my $fh, "gzip -dc $file |") or die;
my $lines;
my $sw = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /[Term]/){
		$sw = 1;
	}
	if($sw == 1 && $line){
		if($line =~ /\[Term\]/ || $line =~ /id: / || $line =~ /name: / || $line =~ /namespace: /){
			$lines .= $line."\n";
			
#			[Term]
#			id: GO:0000001
#			name: mitochondrion inheritance
#			namespace: biological_process
		}
	}
}
close $fh;

my @GO = split(/\[Term\]/, $lines);

my $hash = {};
foreach my $lines (@GO){
	my @L = split(/\n/, $lines);
	my $id;
	my $str;
	my $cat;
	foreach my $line (@L){
		if($line =~ /\[Term\]/ || $line =~ /alt_id/){
			next;
		}
		elsif($line =~ /id: /){
			$id = $line;
			$id =~ s/id: //;
		}
		elsif($line =~ /name: /){
			$str = $line;
			$str =~ s/name: //;
		}
		elsif($line =~ /namespace: /){
			$cat = $line;
			$cat =~ s/namespace: //;
			$cat =~ s/biological_process/Biological Process/;
			$cat =~ s/molecular_function/Molecular Function/;
			$cat =~ s/cellular_component/Cellular Component/;
		}
	}
	
	if($id && $str && $cat){
		$hash->{$id}{str} = $str;
		$hash->{$id}{cat} = $cat;
	}
}

return $hash;
}


#-------------------------------------------------------------------------------
sub Open_data_as_hash{
my $file = shift;
my $pth = shift;
my $hgoobo = shift;

print "! reading [$file] ...\n";
open(my $fh, "<", $file) or die;
my $hash_all = {};
my $hash_sep = {};
my $hash_GO_tm = {};
my $hash_GO_iptm = {};
my $hash_GO_bp = {};
my $hash_GO_mf = {};
my $hash_GO_cl = {};
my $hash_GO_ipr = {};
my $hash_GO_bp2 = {};
my $hash_GO_mf2 = {};
my $hash_GO_cl2 = {};
my $hash_GO_ipr2 = {};
my $cnt = 0;
my $cnt_go = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	my @A = split(/\t/, $line);
	
#	if($A[0] =~ /Glyma/ && $A[0] =~ /Wm82/){
#		$A[0] =~ s/\.Wm82\.a2\.v1//g;
#		$A[0] =~ s/\.Wm82\.a4\.v1//g;
#	}
	$A[12] =~ s/,/ \|/g;
	
	my $t = $A[8]."t";
	if($t eq 't'){
		next;
	}
	
	unless($A[8] eq 'NA' || $A[11] eq 'NULL' || $A[12] eq 'NULL' || $A[8] eq '-' || $A[11] eq '-' || $A[12] eq '-'){
		my @name = split(/\|/, $A[0]);
		my $num_name = @name;
		if($num_name > 1 && $name[0] =~ /TR/ && $name[1] =~ /g/ && $name[2] =~ /m/){
			$A[0] = $name[2];
		}
		
		$hash_all->{$A[0]} .= $line."\n";		#hash of gene ID
		$hash_sep->{$A[3]} .= $line."\n";		#hash of IPRscan program
		
		my $pval;
		my @tmp;
		if($A[8] =~ /e/){
			@tmp = split(/e/, $A[8]);
			$tmp[1] =~ s/\+//;
			$tmp[1] =~ s/\-//;
			$pval = $tmp[0] * 10 ** (- $tmp[1]);
		}
		elsif($A[8] =~ /E/){
			@tmp = split(/E/, $A[8]);
			$tmp[1] =~ s/\+//;
			$tmp[1] =~ s/\-//;
			$pval = $tmp[0] * 10 ** (- $tmp[1]);
		}
		
		my $judge = "false";
		if($pval && $pth){
			if($pval < $pth){
				$judge = "true";
#				print "$pval $pth\n";
#				sleep(1);
			}
		}
		
		if($judge eq 'true'){
			if($A[13] && $A[13] =~ /\(GO:/ && $A[13] =~ /\)/){		#interproscan4
				my @GOs = split(/\), /, $A[13]);
				foreach my $go (@GOs){
					my $p0 = index($go, "(GO\:", 0);
					my $gonum = substr($go, $p0 + 1);
					$gonum =~ s/\)//;
					my $l = length($gonum);
#					$gonum = substr($gonum, 0, $l - 1);
					
					my $gostr = substr($go, 0, $p0 - 1);
					$gostr =~ s/ Biological Process:/Biological Process:/;
					$gostr =~ s/ Molecular Function:/Molecular Function:/;
					$gostr =~ s/ Cellular Component:/Cellular Component:/;
					
#					print "$gostr | $gonum\n";
#					sleep(1);
					
					$hash_GO_tm->{$gonum} = $gostr;
					$gostr =~ s/,//g;
					if($gostr =~ /Biological Process:/){
						$gostr =~ s/Biological Process: //;
						$hash_GO_bp->{$A[0]} .= $gonum.";";
						$hash_GO_bp2->{$A[0]} .= $gostr." (".$gonum.") [".$A[3]." ".$A[8]."];";
					}
					elsif($gostr =~ /Molecular Function:/){
						$gostr =~ s/Molecular Function: //;
						$hash_GO_mf->{$A[0]} .= $gonum.";";
						$hash_GO_mf2->{$A[0]} .= $gostr." (".$gonum.") [".$A[3]." ".$A[8]."];";
					}
					elsif($gostr =~ /Cellular Component:/){
						$gostr =~ s/Cellular Component: //;
						$hash_GO_cl->{$A[0]} .= $gonum.";";
						$hash_GO_cl2->{$A[0]} .= $gostr." (".$gonum.") [".$A[3]." ".$A[8]."];";
					}
				}
				$cnt_go++;
			}
			elsif($A[13] && $A[13] =~ /GO:/){		#interproscan5
				my @GOs = split(/\|/, $A[13]);		#GO:0016021|GO:0016192
				foreach my $go (@GOs){
					my $gonum = $go;
					
					if(! $hgoobo->{$go}{cat} || ! $hgoobo->{$gonum}{str}){
						print "! missing information for [$go]\n";
						sleep(1);
						die;
					}
					
					$hash_GO_tm->{$gonum} = $hgoobo->{$gonum}{cat}.": ".$hgoobo->{$gonum}{str};
					
					if($hgoobo->{$go}{cat} =~ /Biological Process/){
						$hash_GO_bp->{$A[0]} .= $gonum.";";
						$hash_GO_bp2->{$A[0]} .= $hgoobo->{$go}{str}." (".$gonum.") [".$A[3]." ".$A[8]."];";
					}
					elsif($hgoobo->{$go}{cat} =~ /Molecular Function/){
						$hash_GO_mf->{$A[0]} .= $gonum.";";
						$hash_GO_mf2->{$A[0]} .= $hgoobo->{$go}{str}." (".$gonum.") [".$A[3]." ".$A[8]."];";
					}
					elsif($hgoobo->{$go}{cat} =~ /Cellular Component/){
						$hash_GO_cl->{$A[0]} .= $gonum.";";
						$hash_GO_cl2->{$A[0]} .= $hgoobo->{$go}{str}." (".$gonum.") [".$A[3]." ".$A[8]."];";
					}
				}
				$cnt_go++;
			}
		}
		
		unless($A[11] eq 'NULL'){
			$hash_GO_iptm->{$A[11]} = $A[12];
			$hash_GO_ipr->{$A[0]} .= $A[11].";";
			$hash_GO_ipr2->{$A[0]} .= $A[12]." (".$A[11].");";
			
#			if($A[8]){
#				$hash_GO_ipr2->{$A[0]} .= $A[12]." (".$A[11].") [".$A[8]."];";
#			}
#			else{
#				$hash_GO_ipr2->{$A[0]} .= $A[12]." (".$A[11].");";
#			}
		}
	}
	$cnt++;
}
close $fh;

my @Ks0 = keys(%{$hash_sep});
my @GOnums = keys(%{$hash_GO_tm});
my @IPnums = keys(%{$hash_GO_iptm});
my $num_Ks0 = @Ks0;
my $num_GOnums = @GOnums;
my $num_IPnums = @IPnums;

print "! total [$cnt] line output by [$num_Ks0] programs\n";
print "! [$cnt_go] line with [$num_GOnums] GO annotation, [$num_IPnums] IPRID\n";

my $GO_bp = "GO number,class,description\n";
my $GO_mf = "";
my $GO_cc = "";
my $cnt_bp = 0;
my $cnt_mf = 0;
my $cnt_cc = 0;
foreach my $gonum (@GOnums){
	my $gostr = $hash_GO_tm->{$gonum};
	if($gostr =~ /Biological Process:/){
		$GO_bp .= $gonum.",Biological Process,".$gostr."\n";
		$cnt_bp++;
	}
	elsif($gostr =~ /Molecular Function:/){
		$GO_mf .= $gonum.",Molecular Function,".$gostr."\n";
		$cnt_mf++;
	}
	elsif($gostr =~ /Cellular Component:/){
		$GO_cc .= $gonum.",Cellular Component,".$gostr."\n";
		$cnt_cc++;
	}
}

print "! [$cnt_bp] Biological Process\n";
print "! [$cnt_mf] Molecular Function\n";
print "! [$cnt_cc] Cellular Component\n";

my $GO_ip = "IPRID,description\n";
my $cnt_ip = 0;
foreach my $ipnum (@IPnums){
	my $ipstr = $hash_GO_iptm->{$ipnum};
	$GO_ip .= $ipnum.",".$ipstr."\n";
	$cnt_ip++;
}

print "! [$cnt_ip] IPRID\n\n";

my $rhash = {};
$rhash->{all} = $hash_all;
$rhash->{sep} = $hash_sep;
$rhash->{GO_tm} = $hash_GO_tm;
$rhash->{GO_iptm} = $hash_GO_iptm;
$rhash->{GO_bp} = $hash_GO_bp;
$rhash->{GO_mf} = $hash_GO_mf;
$rhash->{GO_cl} = $hash_GO_cl;
$rhash->{GO_ipr} = $hash_GO_ipr;
$rhash->{GO_bp2} = $hash_GO_bp2;
$rhash->{GO_mf2} = $hash_GO_mf2;
$rhash->{GO_cl2} = $hash_GO_cl2;
$rhash->{GO_ipr2} = $hash_GO_ipr2;
$rhash->{program} = \@Ks0;
$rhash->{GOnum} = \@GOnums;
$rhash->{IPnum} = \@IPnums;
$rhash->{GO_description} = $GO_bp.$GO_mf.$GO_cc;
$rhash->{IP_description} = $GO_ip;

return $rhash;
}


#-------------------------------------------------------------------------------
sub GetID{
my $file = shift;

print "\n! reading [$file] ...\n";
open(my $fh, "<", $file) or die;
my $IDs = [];
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ />/){
		my $id = $line;
		$id =~ s/>//;
		
#		if($id =~ /Glyma/ && $id =~ /Wm82/){
#			$id =~ s/\.Wm82\.a2\.v1//g;
#			$id =~ s/\.Wm82\.a4\.v1//g;
#		}
		
		my @name0 = split(/\s/, $id);
		my @name = split(/\|/, $name0[0]);
		my $num_name = @name;
		if($num_name > 1 && $name[0] =~ /TR/ && $name[1] =~ /g/ && $name[2] =~ /m/){
			$id = $name[2];
		}
		
		if($id =~ /Glyma/ && $id =~ /pacid/){
			my @temp = split(/pacid/, $id);
			$id = $temp[0];
			$id =~ s/\s//g;
		}
		elsif($id =~ /MELO/ && $id =~ /gene\=/){
			$id = $name0[0];
		}
		elsif($id =~ /pacid\=/){
			$id = $name0[0];
			$id =~ s/\s//g;
		}
		else{
			$id = $name0[0];
			$id =~ s/\s//g;
		}
		
#		print "$line = $id\n";
#		my $c = <STDIN>;
		
		push(@{$IDs}, $id);
		$cnt++;
	}
}
close $fh;
print "! [$cnt] entries\n";

return $IDs;
}


