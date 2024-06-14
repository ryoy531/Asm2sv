#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;

my $chrID_alias = shift;
my $aliasinfoi = shift;
my $rdir = shift;
my $prev_fasta = shift;
my $prev_gff = shift;
my $new_fasta = shift;
my $new_gff = shift;

if(! $chrID_alias || ! -e $chrID_alias){
	print "! error : missing chrID_alias...\n";
	goto END;
}
if(! $aliasinfoi || ! -e $aliasinfoi){
	print "! error : missing aliasinfoi...\n";
	goto END;
}
if(! $rdir || ! -e $rdir){
	print "! error : missing target directory...\n";
	goto END;
}
if(! $prev_fasta){
	print "! error : missing previous ref-fasta...\n";
	goto END;
}
if(! $prev_gff){
	print "! error : missing previous ref-Gff...\n";
	goto END;
}
if(! $new_fasta){
	print "! error : missing newious ref-fasta...\n";
	goto END;
}
if(! $new_gff){
	print "! error : missing newious ref-Gff...\n";
	goto END;
}

my $hrn = Read_aliaslist($chrID_alias);
if($hrn->{err}){
	goto END;
}

$prev_fasta = Path2name($prev_fasta);
$prev_gff = Path2name($prev_gff);
$new_fasta = Path2name($new_fasta);
$new_gff = Path2name($new_gff);

my $hrmID = Read_IDinfo($aliasinfoi, $prev_gff);

chdir $rdir;
my $gabdir = "log_previous";
unless(-e $gabdir){
	system("mkdir $gabdir");
}

my $hfile = SearchFiles($gabdir);
if($hfile->{err}){
	goto END;
}

my @Cat = keys(%{$hfile});
@Cat = sort {$a cmp $b} @Cat;

my $list_updated = "_list_updated.txt";
my $hupdated = {};
if(-e $list_updated){
	$hupdated = Read_updated($list_updated);
}

foreach my $tdir (@Cat){
#	print " [$tdir] = [$hfile->{$tdir}{cnt}]\n";
	my $htmp = $hfile->{$tdir}{file};
	my @F = keys(%{$htmp});
	@F = sort {$a cmp $b} @F;
	
	foreach my $file (@F){
		unless($hfile->{$file}){
			if($file =~ /\.log/ || $file =~ /final_log_stats_/){
				system("mv $htmp->{$file} $gabdir");
			}
			elsif($file =~ /list_absentJudge_/ && $file =~ /\.csv/){
				system("rm $htmp->{$file}");
			}
			elsif($file =~ /\.psl/){
				system("rm $htmp->{$file}");
			}
			elsif($tdir eq 'top'){
				if($file =~ /selected_predict_/){
					system("rm $htmp->{$file}");
				}
				elsif($file eq $prev_fasta || $file eq $prev_gff){
					system("rm $htmp->{$file}");
				}
				elsif($file =~ /\.psl/ || $file =~ /\.csv/){
					system("rm $htmp->{$file}");
				}
				elsif($file =~ /\.tsv/ && $file !~ /rev_summary_/){
					system("rm $htmp->{$file}");
				}
			}
			elsif($tdir eq 'round_0'){
				if($file =~ /\.fa/ || $file =~ /\.fasta/ || $file =~ /\.gff/ || $file =~ /\.csv/ || $file =~ /\.sh/){
					system("rm $htmp->{$file}");
				}
				elsif($file =~ /\.tsv/ && $file =~ /failed_ID_list/){
					system("rm $htmp->{$file}");
				}
				elsif($file =~ /\.tsv/ && $file =~ /randomID/){
					system("rm $htmp->{$file}");
				}
				elsif($file =~ /r1_hitalign_/ && $file =~ /\.tsv/){
					if(! $hupdated->{$tdir}{$file}){
						my $rh = Revise_record($htmp->{$file}, $hrmID, $new_fasta, $hrn->{f}, 0, "false");
						if($rh->{err}){
							print "! abort script...\n";
							last;
						}
						SAVE($htmp->{$file}, $rh->{lines});
						ADD($list_updated, "$tdir\t$file\t$rh->{cnt_match}\t$rh->{cnt_unmatch}\n");
					}
					else{
						print "! [$tdir] [$file] already updated\n";
					}
				}
				elsif($file =~ /specified_region_hits_/ && $file =~ /\.tsv/){
					if(! $hupdated->{$tdir}{$file}){
						my $rh = Revise_record($htmp->{$file}, $hrmID, $new_fasta, $hrn->{f}, 23, "false");
						if($rh->{err}){
							print "! abort script...\n";
							last;
						}
						SAVE($htmp->{$file}, $rh->{lines});
						ADD($list_updated, "$tdir\t$file\t$rh->{cnt_match}\t$rh->{cnt_unmatch}\n");
					}
					else{
						print "! [$tdir] [$file] already updated\n";
					}
				}
				elsif($file =~ /summary_genome2sv_/ && $file =~ /\.tsv/){
					if(! $hupdated->{$tdir}{$file}){
						my $rh = Revise_record($htmp->{$file}, $hrmID, $new_fasta, $hrn->{f}, 23, "34-38");
						if($rh->{err}){
							print "! abort script...\n";
							last;
						}
						SAVE($htmp->{$file}, $rh->{lines});
						ADD($list_updated, "$tdir\t$file\t$rh->{cnt_match}\t$rh->{cnt_unmatch}\n");
					}
					else{
						print "! [$tdir] [$file] already updated\n";
					}
				}
			}
			elsif($tdir =~ /round_/){
				if($file =~ /\.fa/ || $file =~ /\.fasta/ || $file =~ /\.gff/ || $file =~ /\.csv/ || $file =~ /\.sh/ || $file =~ /hist\.txt/ || $file =~ /log_genelist\.tsv/ || $file =~ /info_ralign\.tsv/){
					system("rm $htmp->{$file}");
				}
				elsif($file =~ /\.tsv/ && $file =~ /failed_ID_list/){
					system("rm $htmp->{$file}");
				}
				elsif($file =~ /\.tsv/ && $file =~ /randomID/){
					system("rm $htmp->{$file}");
				}
				elsif($file =~ /specified_region_hits_/ && $file =~ /\.tsv/){
					if(! $hupdated->{$tdir}{$file}){
						my $rh = Revise_record($htmp->{$file}, $hrmID, $new_fasta, $hrn->{f}, 23, "false");
						if($rh->{err}){
							print "! abort script...\n";
							last;
						}
						SAVE($htmp->{$file}, $rh->{lines});
						ADD($list_updated, "$tdir\t$file\t$rh->{cnt_match}\t$rh->{cnt_unmatch}\n");
					}
					else{
						print "! [$tdir] [$file] already updated\n";
					}
				}
				elsif($file =~ /summary_genome2sv_/ && $file =~ /\.tsv/){
					if(! $hupdated->{$tdir}{$file}){
						my $rh = Revise_record($htmp->{$file}, $hrmID, $new_fasta, $hrn->{f}, 23, "34-38");
						if($rh->{err}){
							print "! abort script...\n";
							last;
						}
						SAVE($htmp->{$file}, $rh->{lines});
						ADD($list_updated, "$tdir\t$file\t$rh->{cnt_match}\t$rh->{cnt_unmatch}\n");
					}
					else{
						print "! [$tdir] [$file] already updated\n";
					}
				}
			}
		}
	}
}

END:{
	my $end = 1;
}


################################################################################
#-------------------------------------------------------------------------------
sub Revise_record{
my $file = shift;
my $hrmID = shift;
my $new_fasta_pref = shift;
my $hrnf = shift;
my $i = shift;
my $prot_cols = shift;

$new_fasta_pref =~ s/\.fasta//;
$new_fasta_pref =~ s/\.fa//;

my @Pcols = split(/-/, $prot_cols);

#print "! reading [$file] ...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
my $cnt = 0;
my $cnt_match = 0;
my $cnt_unmatch = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	unless($line){
		next;
	}
	if($cnt == 0){
		$hash->{lines} .= $line."\n";
		$cnt++;
		next;
	}
	
	my @A = split(/\t/, $line);
	if(! $hrnf->{$A[2]}){
		print "! error : missing seqID alias for [$A[2]]...\n";
		$hash->{err} = 1;
		last;
	}
	if($A[$i] && $hrmID->{$A[$i]}){
		$cnt_unmatch++;
		next;
	}
	else{
		$A[1] = $new_fasta_pref;
		$A[2] = $hrnf->{$A[2]};
		
		if($prot_cols ne 'false'){
			for(my $j = $Pcols[0]; $j <= $Pcols[1]; $j++){
				$A[$j] = "-";
			}
		}
		
		$hash->{lines} .= join("\t", @A)."\n";
		$cnt_match++;
	}
	$cnt++;
}
close $fh;

$hash->{cnt_match} = $cnt_match;
$hash->{cnt_unmatch} = $cnt_unmatch;

print "! revised [$file] | [$cnt_match] kept, [$cnt_unmatch] removed\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub SearchFiles{
my $gabdir = shift;

my $tmp = "__tmp_searchfiles.txt";
system("find . | grep -v $gabdir > $tmp");

my $hash = {};
if(-e $tmp){
	open(my $fh, "<", $tmp);
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		unless($line){
			next;
		}
		my @D = split(/\//, $line);
		my $nD = @D;
		if($nD == 2){
			if($hash->{top}{file}{$D[$nD - 1]}){
				print "! duplicated file name : $D[$nD - 1]\n";
				$hash->{err} = 1;
				last;
			}
			$hash->{top}{file}{$D[$nD - 1]} = $line;
			$hash->{top}{cnt} += 1;
		}
		elsif($nD > 2){
			my $tdir = $D[1];
			if($hash->{$tdir}{file}{$D[$nD - 1]}){
				print "! duplicated file name : $D[$nD - 1]\n";
				$hash->{err} = 1;
				last;
			}
			$hash->{$tdir}{file}{$D[$nD - 1]} = $line;
			$hash->{$tdir}{cnt} += 1;
		}
#		else{
#			print "$line\n";
#		}
	}
	close $fh;
	
	system("rm $tmp");
}

return $hash;
}


#-------------------------------------------------------------------------------
sub Path2name{
my $file = shift;

my @PRF = split(/\//, $file);
my $num_PRF = @PRF;
my $prefix = $PRF[$num_PRF - 1];

return $prefix;
}


#-------------------------------------------------------------------------------
sub Prefix_Gff{
my $file = shift;

my @PRF = split(/\//, $file);
my $num_PRF = @PRF;
my $prefix = $PRF[$num_PRF - 1];
$prefix =~ s/\.gff3//;
$prefix =~ s/\.gff//;

return $prefix;
}


#-------------------------------------------------------------------------------
sub Read_IDinfo{
my $file = shift;
my $prev_gff = shift;

print "! reading [$file] ...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
my $cnt_match = 0;
my $cnt_unmatch = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	unless($line){
		next;
	}
	
	# gff0,gene_id,transcript_id,seqid,pos0,pos1,strand,_,gff1,gene_id,transcript_id,seqid,pos0,pos1,strand,judge
	
	my @A = split(/,/, $line);
	if(! $A[15]){
		next;
	}
	
	my @PRF0 = split(/\//, $A[0]);
	my $nprf = @PRF0;
	if($prev_gff eq $PRF0[$nprf - 1]){
		if($A[15] eq 'match'){
			$cnt_match++;
		}
		else{
			$hash->{$A[1]} = 1;
			$hash->{$A[2]} = 1;
			$cnt_unmatch++;
		}
	}
}
close $fh;

print "! [$cnt_match] entries matched (keep)\n";
print "! [$cnt_unmatch] unmatched (remove)\n";

return $hash;
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


#-------------------------------------------------------------------------------
sub Read_updated{
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
	my @A = split(/\t/, $line);
	if($A[0] && $A[1]){
		$hash->{$A[0]}{$A[1]} = 1;
	}
}
close $fh;

return $hash;
}


#----------------------------------------------------------
sub SAVE{
my $file = shift;
my $str = shift;

if($str){
	open(my $fh, ">", $file) or die;
	print $fh $str;
	close $fh;
	
#	print "! output [$file]\n";
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


