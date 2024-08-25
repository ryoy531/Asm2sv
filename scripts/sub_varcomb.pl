#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use FindBin;

my $bin_matcher = shift;
my $aseq = shift;
my $bseq = shift;
my $rfile = shift;
my $afile = shift;
my $limit_len = shift;

if(! $bin_matcher || ! -e $bin_matcher){
	die;
}
if(! $aseq || ! -e $aseq){
	die;
}
if(! $bseq || ! -e $bseq){
	die;
}
if(! $rfile){
	die;
}
if(! $afile){
	die;
}
if(! defined $limit_len){
	$limit_len = 100;
}

my $alen = Len_fasta($aseq);
my $blen = Len_fasta($bseq);
if($alen < $limit_len || $blen < $limit_len){
	open(my $rfh, ">", $afile);
	print $rfh "null";
	close $rfh;
}
else{
	print "! $alen bp | $blen bp\n";
	
	my $cmd = "$bin_matcher -asequence $aseq -bsequence $bseq -outfile $rfile > /dev/null 2>&1";
	if(system($cmd) != 0){
		print "! matcher failed | $aseq vs $bseq\n";
	#	Rm($aseq);
	#	Rm($bseq);
		die;
	}
	
	Judge_alignment($rfile, $afile, $alen, $blen);
}

END:{
	my $end = 1;
}


############################################################
#-----------------------------------------------------------
sub Judge_alignment{
my $file = shift;
my $afile = shift;
my $alen = shift;
my $blen = shift;

if($file && -e $file){
	open(my $fh, "<", $file);
	my $sw0 = 0;
	my $sw1 = 0;
	my $sw2 = 0;
	my @PosA;
	my @PosB;
	my $seqA;
	my $seqB;
	my $aln = "";
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line){
			if($sw0 == 0 && $line =~ /Score/){
				$sw0 = 1;
			}
			if($sw0 == 1 && $line =~ /=======================================/){
				$sw1 = 1;
			}
			if($sw1 == 1){
				if($line =~ /---------------------------------------/){
					last;
				}
				elsif($line =~ /:/){
					$line =~ s/\s//g;
					$aln .= $line;
					next;
				}
				
				my @Vals = split(/\s/, $line);
				if($sw2 == 0){
					if($line !~ /[a-z]/i){
						foreach my $val (@Vals){
							if(defined $val && $val =~ /\d/){
#								print "A | $val\n";
								push(@PosA, $val);
							}
						}
					}
					else{
						$line =~ s/\s//g;
						$line =~ s/-/N/g;
						$line =~ s/aseq//g;
						$seqA .= $line;
						$sw2 = 1;
					}
				}
				elsif($sw2 == 1){
					if($line !~ /[a-z]/i){
						foreach my $val (@Vals){
							if(defined $val && $val =~ /\d/){
#								print "B | $val\n";
								push(@PosB, $val);
							}
						}
						$sw2 = 0;
					}
					else{
						$line =~ s/\s//g;
						$line =~ s/-/N/g;
						$line =~ s/bseq//g;
						$seqB .= $line;
					}
				}
			}
		}
	}
	close $fh;
	
	my $ans = "null";
	my $len_match = length($aln);
	if($len_match > 0){
		my $nA = @PosA;
		my $nB = @PosB;
		my $lenA = length($seqA);
		my $lenB = length($seqB);
		my $ratioA = $len_match / $alen;
		my $ratioB = $len_match / $blen;
		
		my $woN_seqA = $seqA;
		my $woN_seqB = $seqB;
		$woN_seqA =~ s/N//g;
		$woN_seqB =~ s/N//g;
		my $woN_lenA = length($woN_seqA);
		my $woN_lenB = length($woN_seqB);
		
		if(@PosA){
			@PosA = sort{$a <=> $b} @PosA;
		}
		if(@PosB){
			@PosB = sort{$a <=> $b} @PosB;
		}
		
		print "$PosA[0] - $PosA[$nA - 1]\n";
		print "$PosB[0] - $PosB[$nB - 1]\n";
		print "A = $lenA ($ratioA) | woN = $woN_lenA\n";
		print "B = $lenB ($ratioB) | woN = $woN_lenB\n";
#		print "\n>A\n$seqA\n>B\n$seqB\n";
		
		if($ratioA > 0.9 && $ratioB > 0.9){
			$ans = $seqA;
			if($woN_lenB > $woN_lenA){
				$ans = $seqB;
			}
		}
	}
	
	open(my $rfh, ">", $afile);
	print $rfh $ans;
	close $rfh;
}

}


#-----------------------------------------------------------
sub Len_fasta{
my $file = shift;

if($file && -e $file){
	open(my $fh, "<", $file);
	my $seq = "";
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line){
			if($line =~ /\>/){
				next;
			}
			else{
				$line =~ s/\*//g;
				$seq .= $line;
			}
		}
	}
	close $fh;
	
	my $len = length($seq);
	return $len;
}
else{
	return 0;
}

}


#-----------------------------------------------------------
sub Rm{
my $file = shift;

if($file && -e $file){
	system("rm $file");
}

}

