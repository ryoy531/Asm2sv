#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "Gff3_diffIDpos_2files.pl version 1.01\n";
$version .= "last update: [2023\/3\/26]\n";
$version .= "copyright: ryoichi yano [ryoichiy104\@gmail.com]\n";

#print $version;

#-------------------------------------------------------------------------------

# 同一のリファレンスfastaに紐づいた異なる2つのGFFにて、gene/transcript IDとposition(start, end)が異なるものを出力する
# アノテーションをアップデートした際に、情報が変わった遺伝子を探す

my $gff0 = shift;
my $gff1 = shift;
my $chrID_alias = shift;
my $rfile0 = shift;
my $rfile1 = shift;
my $rfilei = shift;

my $wpath = getcwd();

if(! $gff0 || ! -e $gff0){
	print "! missing gff0...\n";
	goto END;
}
if(! $gff1 || ! -e $gff1){
	print "! missing gff1...\n";
	goto END;
}
if(! $chrID_alias || ! -e $chrID_alias){
	print "! missing chrID_alias...\n";
	goto END;
}

my $q = "y";
if($q =~ /y/i || $q =~ /pipe/i){
	my $hrn = Read_aliaslist($chrID_alias);
	
	my @PRF0 = split(/\//, $gff0);
	my $num_PRF0 = @PRF0;
	my $prefix0 = $PRF0[$num_PRF0 - 1];
	$prefix0 =~ s/\.gff3//;
	$prefix0 =~ s/\.gff//;
	
	my @PRF1 = split(/\//, $gff1);
	my $num_PRF1 = @PRF1;
	my $prefix1 = $PRF0[$num_PRF1 - 1];
	$prefix1 =~ s/\.gff3//;
	$prefix1 =~ s/\.gff//;
	
	my $rh0 = Open_gff($gff0, "previous", $hrn->{f});
	my $rh1 = Open_gff($gff1, "new");
	my $t2g0 = $rh0->{t2g};
	my $t2g1 = $rh1->{t2g};
	my $g2t0 = $rh0->{g2t};
	my $g2t1 = $rh1->{g2t};
	
	my $hgene0 = $rh0->{hgene};
	my $hgene1 = $rh1->{hgene};
	my $GID0 = $rh0->{GID};
	my $GID1 = $rh1->{GID};
	my $htranscript0 = $rh0->{htranscript};
	my $htranscript1 = $rh1->{htranscript};
	
	my $cnt0 = 0;
	my $cnt_diff0 = 0;
	my $cnt_trdiff0 = 0;
	my $cnt_contrad0 = 0;
	my $miss_gff0 = "##gff-version 3\n";
	my $rdiff0 = "";
	my $rinfo = "gff0,gene_id,transcript_id,seqid,pos0,pos1,strand,_,gff1,gene_id,transcript_id,seqid,pos0,pos1,strand,judge\n";
	foreach my $gid0 (@{$GID0}){
		my $judge = "false";
		if(defined $hgene1->{$gid0}{seqid} && defined $hgene1->{$gid0}{pos0} && defined $hgene1->{$gid0}{pos1} && defined $hgene1->{$gid0}{strand}){
			if($hgene0->{$gid0}{seqid} eq $hgene1->{$gid0}{seqid} || $hgene1->{$gid0}{seqid} eq $hgene0->{$gid0}{seqid}){
				if($hgene0->{$gid0}{pos0} eq $hgene1->{$gid0}{pos0} && $hgene0->{$gid0}{pos1} eq $hgene1->{$gid0}{pos1} && $hgene0->{$gid0}{strand} eq $hgene1->{$gid0}{strand}){
					$judge = "match";
				}
			}
		}
		
		if($judge eq 'false'){
			$rdiff0 .= $gid0."\n";
			$cnt_diff0++;
			$rinfo .= "$gff0,$gid0,-,$hgene0->{$gid0}{seqid_ori},$hgene0->{$gid0}{pos0},$hgene0->{$gid0}{pos1},$hgene0->{$gid0}{strand},->,$gff1,-,-,-,-,-,-,NA\n";
			$miss_gff0 .= $hgene0->{$gid0}{data};
			
			if($g2t0->{$gid0}){
				my @TID0 = split(/\n/, $g2t0->{$gid0});
				foreach my $tid0 (@TID0){
					$rinfo .= "$gff0,$gid0,$tid0,$htranscript0->{$tid0}{seqid_ori},$htranscript0->{$tid0}{pos0},$htranscript0->{$tid0}{pos1},$htranscript0->{$tid0}{strand},->,$gff1,-,-,-,-,-,-,NA\n";
					$miss_gff0 .= $htranscript0->{$tid0}{data};
					$cnt_trdiff0++;
				}
			}
		}
		else{
			my $tmp_cnt_trdiff0 = 0;
			if($g2t0->{$gid0}){
				my @TID0 = split(/\n/, $g2t0->{$gid0});
				foreach my $tid0 (@TID0){
					my $trjudge = "false";
					if(defined $htranscript1->{$tid0}{seqid} && defined $htranscript1->{$tid0}{pos0} && defined $htranscript1->{$tid0}{pos1} && defined $htranscript1->{$tid0}{strand}){
						if($htranscript0->{$tid0}{pos0} eq $htranscript1->{$tid0}{pos0} && $htranscript0->{$tid0}{pos1} eq $htranscript1->{$tid0}{pos1} && $htranscript0->{$tid0}{strand} eq $htranscript1->{$tid0}{strand}){
							$trjudge = "match";
						}
					}
					if($trjudge eq 'false'){
						$cnt_trdiff0++;
						$tmp_cnt_trdiff0++;
					}
				}
			}
			
			if($tmp_cnt_trdiff0 > 0){
				$rinfo .= "$gff0,$gid0,-,$hgene0->{$gid0}{seqid_ori},$hgene0->{$gid0}{pos0},$hgene0->{$gid0}{pos1},$hgene0->{$gid0}{strand},->,$gff1,-,-,-,-,-,-,NA\n";
				$miss_gff0 .= $hgene0->{$gid0}{data};
				
				if($g2t0->{$gid0}){
					my @TID0 = split(/\n/, $g2t0->{$gid0});
					foreach my $tid0 (@TID0){
						$rinfo .= "$gff0,$gid0,$tid0,$htranscript0->{$tid0}{seqid_ori},$htranscript0->{$tid0}{pos0},$htranscript0->{$tid0}{pos1},$htranscript0->{$tid0}{strand},->,$gff1,-,-,-,-,-,-,NA\n";
						$miss_gff0 .= $htranscript0->{$tid0}{data};
					}
				}
				$cnt_contrad0++;
			}
			else{
				$rinfo .= "$gff0,$gid0,-,$hgene0->{$gid0}{seqid_ori},$hgene0->{$gid0}{pos0},$hgene0->{$gid0}{pos1},$hgene0->{$gid0}{strand},->,$gff1,$gid0,-,$hgene1->{$gid0}{seqid_ori},$hgene1->{$gid0}{pos0},$hgene1->{$gid0}{pos1},$hgene1->{$gid0}{strand},$judge\n";
				
				if($g2t0->{$gid0}){
					my @TID0 = split(/\n/, $g2t0->{$gid0});
					foreach my $tid0 (@TID0){
						$rinfo .= "$gff0,$gid0,$tid0,$htranscript0->{$tid0}{seqid_ori},$htranscript0->{$tid0}{pos0},$htranscript0->{$tid0}{pos1},$htranscript0->{$tid0}{strand},->,$gff1,$gid0,$tid0,$htranscript1->{$tid0}{seqid_ori},$htranscript1->{$tid0}{pos0},$htranscript1->{$tid0}{pos1},$htranscript1->{$tid0}{strand},$judge\n";
					}
				}
			}
		}
		$cnt0++;
	}
	
	my $cnt1 = 0;
	my $cnt_diff1 = 0;
	my $cnt_trdiff1 = 0;
	my $cnt_contrad1 = 0;
	my $miss_gff1 = "##gff-version 3\n";
	my $rdiff1 = "";
	foreach my $gid1 (@{$GID1}){
		my $judge = "false";
		if(defined $hgene0->{$gid1}{seqid} && defined $hgene0->{$gid1}{pos0} && defined $hgene0->{$gid1}{pos1} && defined $hgene0->{$gid1}{strand}){
			if($hgene0->{$gid1}{seqid} eq $hgene1->{$gid1}{seqid} || $hgene1->{$gid1}{seqid} eq $hgene0->{$gid1}{seqid}){
				if($hgene0->{$gid1}{pos0} eq $hgene1->{$gid1}{pos0} && $hgene0->{$gid1}{pos1} eq $hgene1->{$gid1}{pos1} && $hgene0->{$gid1}{strand} eq $hgene1->{$gid1}{strand}){
					$judge = "match";
				}
			}
		}
		
		if($judge eq 'false'){
			$rdiff1 .= $gid1."\n";
			$cnt_diff1++;
			$rinfo .= "$gff1,$gid1,-,$hgene1->{$gid1}{seqid_ori},$hgene1->{$gid1}{pos0},$hgene1->{$gid1}{pos1},$hgene1->{$gid1}{strand},->,$gff0,-,-,-,-,-,-,NA\n";
			$miss_gff1 .= $hgene1->{$gid1}{data};
			
			if($g2t1->{$gid1}){
				my @TID1 = split(/\n/, $g2t1->{$gid1});
				foreach my $tid1 (@TID1){
					$rinfo .= "$gff1,$gid1,$tid1,$htranscript1->{$tid1}{seqid_ori},$htranscript1->{$tid1}{pos0},$htranscript1->{$tid1}{pos1},$htranscript1->{$tid1}{strand},->,$gff0,-,-,-,-,-,-,NA\n";
					$miss_gff1 .= $htranscript1->{$tid1}{data};
					$cnt_trdiff1++;
				}
			}
		}
		else{
			my $tmp_cnt_trdiff1 = 0;
			if($g2t1->{$gid1}){
				my @TID1 = split(/\n/, $g2t1->{$gid1});
				foreach my $tid1 (@TID1){
					my $trjudge = "false";
					if(defined $htranscript0->{$tid1}{seqid} && defined $htranscript0->{$tid1}{pos0} && defined $htranscript0->{$tid1}{pos1} && defined $htranscript0->{$tid1}{strand}){
						if($htranscript0->{$tid1}{pos0} eq $htranscript1->{$tid1}{pos0} && $htranscript0->{$tid1}{pos1} eq $htranscript1->{$tid1}{pos1} && $htranscript0->{$tid1}{strand} eq $htranscript1->{$tid1}{strand}){
							$trjudge = "match";
						}
					}
					if($trjudge eq 'false'){
						$cnt_trdiff1++;
						$tmp_cnt_trdiff1++;
					}
				}
			}
			
			if($tmp_cnt_trdiff1 > 0){
				$rinfo .= "$gff1,$gid1,-,$hgene1->{$gid1}{seqid_ori},$hgene1->{$gid1}{pos0},$hgene1->{$gid1}{pos1},$hgene1->{$gid1}{strand},->,$gff0,-,-,-,-,-,-,NA\n";
				$miss_gff1 .= $hgene1->{$gid1}{data};
				
				if($g2t1->{$gid1}){
					my @TID1 = split(/\n/, $g2t1->{$gid1});
					foreach my $tid1 (@TID1){
						$rinfo .= "$gff1,$gid1,$tid1,$htranscript1->{$tid1}{seqid_ori},$htranscript1->{$tid1}{pos0},$htranscript1->{$tid1}{pos1},$htranscript1->{$tid1}{strand},->,$gff0,-,-,-,-,-,-,NA\n";
						$miss_gff1 .= $htranscript1->{$tid1}{data};
					}
				}
				$cnt_contrad1++;
			}
			else{
				$rinfo .= "$gff1,$gid1,-,$hgene1->{$gid1}{seqid_ori},$hgene1->{$gid1}{pos0},$hgene1->{$gid1}{pos1},$hgene1->{$gid1}{strand},->,$gff0,$gid1,-,$hgene0->{$gid1}{seqid_ori},$hgene0->{$gid1}{pos0},$hgene0->{$gid1}{pos1},$hgene0->{$gid1}{strand},$judge\n";
				
				if($g2t1->{$gid1}){
					my @TID1 = split(/\n/, $g2t1->{$gid1});
					foreach my $tid1 (@TID1){
						$rinfo .= "$gff1,$gid1,$tid1,$htranscript1->{$tid1}{seqid_ori},$htranscript1->{$tid1}{pos0},$htranscript1->{$tid1}{pos1},$htranscript1->{$tid1}{strand},->,$gff0,$gid1,$tid1,$htranscript0->{$tid1}{seqid_ori},$htranscript0->{$tid1}{pos0},$htranscript0->{$tid1}{pos1},$htranscript0->{$tid1}{strand},$judge\n";
					}
				}
			}
		}
		$cnt1++;
	}
	
	print "! [$cnt0] genes in [previous]\n";
	print "! [$cnt1] genes in [new]\n";
	print "! stats : present in [previous] but absent from [new]\n";
	print "  [$cnt_diff0] genes\n";
	print "  [$cnt_trdiff0] transcripts\n";
	print "  [$cnt_contrad0] unmatch between gene and transcript\n";
	print "! stats : present in [new]  but absent from [previous]\n";
	print "  [$cnt_diff1] genes\n";
	print "  [$cnt_trdiff1] transcripts\n";
	print "  [$cnt_contrad1] unmatch between gene and transcript\n";
	
	if(! $rfile0){
		$rfile0 = "diffIDpos_".$prefix0."_absentFrom-".$prefix1.".gff3";
	}
	if(! $rfile1){
		$rfile1 = "diffIDpos_".$prefix1."_absentFrom-".$prefix0.".gff3";
	}
	if(! $rfilei){
		$rfilei = "diffIDpos_summary_".$prefix1."_vs_".$prefix0.".csv";
	}
	SAVE($rfile0, $miss_gff0);
	SAVE($rfile1, $miss_gff1);
	SAVE($rfilei, $rinfo);
}

END:{
	my $end = 1;
}


################################################################################
#-------------------------------------------------------------------------------
sub Open_gff{
my $file = shift;
my $str = shift;
my $hrnf = shift;

print "! reading [$file] ($str)...\n";
open(my $fh, "<", $file) or die;
my $t2g = {};
my $g2t = {};
my $t2g_lines = "";
my $hgene = {};
my $htranscript = {};
my @GID;
my $gcnt = 0;
my @L;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if(! $line){
		next;
	}
	
	my @Char = split(//, $line);
	if($Char[0] =~ /\#/){
		next;
	}
	
	my @A = split(/\t/, $line);
	unless($line){
#		print "! missing line\n";
		next;
	}
	
	my $seqid_ori = $A[0];
	if($hrnf->{$A[0]}){
		$A[0] = $hrnf->{$A[0]};
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
		
		push(@GID, $gid);
		
		$hgene->{$gid}{data} = $line."\n";
		$hgene->{$gid}{seqid} = $A[0];
		$hgene->{$gid}{seqid_ori} = $seqid_ori;
		$hgene->{$gid}{pos0} = $A[3];
		$hgene->{$gid}{pos1} = $A[4];
		$hgene->{$gid}{strand} = $A[6];
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
			if($gid && $tid){
				last;
			}
		}
		
		$t2g->{$tid} = $gid;
		$t2g_lines .= $tid."\t".$gid."\n";
		$g2t->{$gid} .= $tid."\n";
		$hgene->{$gid}{num_variant} += 1;
		
		if($htranscript->{$tid}{data}){
			print "! caution : Gff format error at [$tid]...\n";
		}
		$htranscript->{$tid}{data} .= $line."\n";
		$htranscript->{$tid}{seqid} = $A[0];
		$htranscript->{$tid}{seqid_ori} = $seqid_ori;
		$htranscript->{$tid}{pos0} = $A[3];
		$htranscript->{$tid}{pos1} = $A[4];
		$htranscript->{$tid}{strand} = $A[6];
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
		
		if($tid){
			if(! $htranscript->{$tid}{data}){
				print "! caution : Gff format error at [$tid]...\n";
			}
			$htranscript->{$tid}{data} .= $line."\n";
		}
	}
	
	push(@L, $line);
}
close $fh;

foreach my $line (@L){
	my @A = split(/\t/, $line);
	my @tag = split(/\;/, $A[8]);
	
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
			print "! missing gid for [$tid]\n";
			sleep(1);
		}
		
		my $gid = $t2g->{$tid};
		$hgene->{$gid}{num_CDS} += 1;
		$htranscript->{$tid}{num_CDS} += 1;
		$htranscript->{$tid}{pos_CDS}{$A[3]} = $A[4];
	}
	elsif($A[2] eq 'exon'){
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
			print "! missing gid for [$tid]\n";
			sleep(1);
		}
		
		my $gid = $t2g->{$tid};
		$hgene->{$gid}{num_exon} += 1;
		$htranscript->{$tid}{num_exon} += 1;
		$htranscript->{$tid}{pos_exon}{$A[3]} = $A[4];
	}
}

print " [$gcnt] genes\n";

my $numpc = 0;
my $numnc = 0;
my $numpc_tr = 0;
my $hcnt = {};
foreach my $gid (@GID){
	unless($hgene->{$gid}{num_variant}){
		$hgene->{$gid}{num_variant} = 0;
	}
	unless($hgene->{$gid}{num_CDS}){
		$hgene->{$gid}{num_CDS} = 0;
		$numnc++;
	}
	else{
		$numpc++;
		$numpc_tr += $hgene->{$gid}{num_variant};
	}
}

print " [$numpc] protein-coding, [$numpc_tr] transcripts\n";
print " [$numnc] non-coding\n";

my $rh = {};
$rh->{hgene} = $hgene;
$rh->{htranscript} = $htranscript;
$rh->{t2g} = $t2g;
$rh->{g2t} = $g2t;
$rh->{GID} = \@GID;

return $rh;
}


#-------------------------------------------------------------------------------
sub Read_aliaslist{
my $file = shift;

#print "! reading [$file] ...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	unless($line){
		next;
	}
	my @A = split(/,/, $line);
	
	if($A[0] && $A[1] && $A[0] ne $A[1]){
		if($A[0] =~ /\.fa/ || $A[1] =~ /\.fa/){
			if(! $hash->{fasta0} && ! $hash->{fasta1}){
				$hash->{fasta0} = $A[0];
				$hash->{fasta1} = $A[1];
				
				$A[0] =~ s/\.fasta//;
				$A[1] =~ s/\.fasta//;
				$A[0] =~ s/\.fa//;
				$A[1] =~ s/\.fa//;
				$hash->{prefix0} = $A[0];
				$hash->{prefix1} = $A[1];
#				print " $A[0] -> $A[1] (fasta)\n";
			}
		}
		else{
			if(! $hash->{f}{$A[0]} && ! $hash->{r}{$A[1]}){
				$hash->{f}{$A[0]} = $A[1];
				$hash->{r}{$A[1]} = $A[0];
#				print " $A[0] -> $A[1]\n";
			}
			elsif($hash->{f}{$A[0]}){
#				print "! duplicated ID alias for [$A[0]]\n";
				$hash->{err} = 1;
				last;
			}
			elsif($hash->{r}{$A[1]}){
#				print "! duplicated ID alias for [$A[1]]\n";
				$hash->{err} = 1;
				last;
			}
		}
	}
}
close $fh;

return $hash;
}


#----------------------------------------------------------
sub SAVE{
my $file = shift;
my $str = shift;

open(my $fh, ">", $file) or die;
print $fh $str;
close $fh;

print "! output [$file]\n";

}


#----------------------------------------------------------
sub File_select{
my $keyword = shift;
my $str = shift;
my $akw = shift;

my @files2;
my $q1;

unless($keyword){
	print "\nEnter keyword for file search: ";
	$keyword = <STDIN>;
	$keyword =~ s/\n//;
	$keyword =~ s/\r//;
}

my @files = glob("*");

foreach my $file (@files){
	if($file =~ /$keyword/){
		unless($file =~ /\.fai/ || $file =~ /\.fa\.n../ || $file =~ /\.fa\.p../ || $file =~ /\.fasta\.n../ || $file =~ /\.fasta\.p../ || $file =~ /\.ffa\.n../ || $file =~ /diffID/){
			if($akw){
				if($file !~ /$akw/){
					push(@files2, $file);
				}
			}
			else{
				push(@files2, $file);
			}
		}
	}
}

print "\n---------------------------------------------------------------------\n";
my $cnt0 = 0;
foreach my $file (@files2){
	print "[$cnt0] $file\n";
	$cnt0++;
}
print "---------------------------------------------------------------------\n";

unless($q1){
	print "Select $str: ";
	$q1 = <STDIN>;
	$q1 =~ s/\n//;
	$q1 =~ s/\r//;
}
unless($q1){
	$q1 = 0;
}

if($files2[$q1] && -e $files2[$q1]){
	print "[$files2[$q1]]\n";
}

return $files2[$q1];
}



