#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;

my $fasta = shift;
my $gff3 = shift;
my $rfasta = shift;

if(! $fasta || ! -e $fasta){
	print "! gfftoprot.pl | missing input fasta...\n";
	goto END;
}
if(! $gff3 || ! -e $gff3){
	print "! gfftoprot.pl | missing input GFF3...\n";
	goto END;
}

print "! translating CDS from [$gff3]...\n";
my $hgff = Open_gff($gff3);
my $hfasta = Open_fasta_as_hash($fasta);
my $GID = $hgff->{GID};
my $g2t = $hgff->{g2t};
my $hash = $hgff->{data};

my $cnt_gid = 0;
my $cnt_tid = 0;
my $cnt_trans = 0;
my $cnt_err = 0;
my $combined_prot = "";
my $seg = 10000;
my $add = $seg;
foreach my $gid (@{$GID}){
	if($g2t->{$gid}){
		$cnt_gid++;
		my @TID = split(/\n/, $g2t->{$gid});
		
		foreach my $tid (@TID){
			if($hash->{num_CDS}{$tid} && $hash->{num_CDS}{$tid} > 0){
				$cnt_tid++;
				my $err = 0;
				my $strand = $hash->{strand}{$tid};
				my @L = split(/\n/, $hash->{CDS}{$tid});
				
				my $AoCDS = [];
				my $sid = "";
				foreach my $line (@L){
					my @A = split(/\t/, $line);
					unless($sid){
						$sid = $A[0];
					}
					elsif($sid ne $A[0]){
						$err++;
						last;
					}
					push(@{$AoCDS}, \@A);
				}
				if($err){
					$cnt_err++;
					next;
				}
				
				@{$AoCDS} = sort {$a->[3] <=> $b->[3]} @{$AoCDS};
				
				my $chrseq = $hfasta->{$sid}{seq};
				my $cds;
				foreach my $A (@{$AoCDS}){
					$cds .= substr($chrseq, $A->[3] - 1, $A->[4] - $A->[3] + 1);
				}
				$cds =~ tr/acgtn/ACGTN/;
				
				if($strand eq '-'){
					my $revcomp = $cds;
					$revcomp =~ tr/ACGT/TGCA/;
					$revcomp = reverse($revcomp);
					$cds = $revcomp;
				}
				my $len = length($cds);
				
				my $prot = "";
				for(my $k = 0; $k < $len; $k += 3){
					my $codon0 = "";
					if($k + 3 <= $len){
						$codon0 = substr($cds, $k, 3);
					}
					else{
						last;
					}
					
					my $aa0s;
					if($codon0 =~ /gct/i || $codon0 =~ /gcc/i || $codon0 =~ /gca/i || $codon0 =~ /gcg/i){$aa0s = "A";}
					elsif($codon0 =~ /tta/i || $codon0 =~ /ttg/i || $codon0 =~ /ctt/i || $codon0 =~ /ctc/i || $codon0 =~ /cta/i || $codon0 =~ /ctg/i){$aa0s = "L";}
					elsif($codon0 =~ /cgt/i || $codon0 =~ /cgc/i || $codon0 =~ /cga/i || $codon0 =~ /cgg/i || $codon0 =~ /aga/i || $codon0 =~ /agg/i){$aa0s = "R";}
					elsif($codon0 =~ /aaa/i || $codon0 =~ /aag/i){$aa0s = "K";}
					elsif($codon0 =~ /aat/i || $codon0 =~ /aac/i){$aa0s = "N";}
					elsif($codon0 =~ /atg/i){$aa0s = "M";}
					elsif($codon0 =~ /gat/i || $codon0 =~ /gac/i){$aa0s = "D";}
					elsif($codon0 =~ /ttt/i || $codon0 =~ /ttc/i){$aa0s = "F";}
					elsif($codon0 =~ /tgt/i || $codon0 =~ /tgc/i){$aa0s = "C";}
					elsif($codon0 =~ /cct/i || $codon0 =~ /ccc/i || $codon0 =~ /cca/i || $codon0 =~ /ccg/i){$aa0s = "P";}
					elsif($codon0 =~ /caa/i || $codon0 =~ /cag/i){$aa0s = "Q";}
					elsif($codon0 =~ /tct/i || $codon0 =~ /tcc/i || $codon0 =~ /tca/i || $codon0 =~ /tcg/i || $codon0 =~ /agt/i || $codon0 =~ /agc/i){$aa0s = "S";}
					elsif($codon0 =~ /gaa/i || $codon0 =~ /gag/i){$aa0s = "E";}
					elsif($codon0 =~ /act/i || $codon0 =~ /acc/i || $codon0 =~ /aca/i || $codon0 =~ /acg/i){$aa0s = "T";}
					elsif($codon0 =~ /ggt/i || $codon0 =~ /ggc/i || $codon0 =~ /gga/i || $codon0 =~ /ggg/i){$aa0s = "G";}
					elsif($codon0 =~ /tgg/i){$aa0s = "W";}
					elsif($codon0 =~ /cat/i || $codon0 =~ /cac/i){$aa0s = "H";}
					elsif($codon0 =~ /tat/i || $codon0 =~ /tac/i){$aa0s = "Y";}
					elsif($codon0 =~ /att/i || $codon0 =~ /atc/i || $codon0 =~ /ata/i){$aa0s = "I";}
					elsif($codon0 =~ /gtt/i || $codon0 =~ /gtc/i || $codon0 =~ /gta/i || $codon0 =~ /gtg/i){$aa0s = "V";}
					elsif($codon0 =~ /tag/i || $codon0 =~ /tga/i || $codon0 =~ /taa/i){$aa0s = "*";}
					
					if($aa0s){
						$prot .= $aa0s;
					}
					if($aa0s && $aa0s eq '*'){
						last;
					}
				}
				if($prot){
					$combined_prot .= ">".$tid."\n".$prot."\n";
					$cnt_trans++;
					
					if($cnt_trans == $seg){
#						print " [$cnt_trans]\n";
						$seg += $add;
					}
				}
			}
		}
	}
}

print " [$cnt_gid] gene ID\n";
print " [$cnt_tid] transcript with CDS\n";
print " [$cnt_err] abondaned due to contradicting seq ID\n";
print " [$cnt_trans] translated to protein seq\n";

unless($rfasta){
	$rfasta = $gff3;
	$rfasta =~ s/\.gff3//;
	$rfasta =~ s/\.gff//;
	$rfasta .= ".fasta";
}

SAVE($rfasta, $combined_prot);

END:{
	my $end = 1;
}



################################################################################
#-------------------------------------------------------------------------------
sub Open_fasta_as_hash{
my $file = shift;

#print "! reading [$file] ...\n";
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

#print "! total [$numID] sequence ID, [$total_len] bp\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Open_gff{
my $gff3 = shift;

#print "! reading [$gff3]...\n";
open(my $fh, "<", $gff3) or die;
my $t2g = {};
my $g2t = {};
my $hash = {};
my @GID;
my @GID2;
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
	if(! $A[2] || ! $A[8]){
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
		
		push(@GID, $gid);
		$hash->{gene}{$gid} = $line;
	}
	elsif($A[2] eq 'transcript' || $A[2] eq 'mRNA'){
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
		
		if(! $gid && $tid){
			$gid = $tid;
			push(@GID2, $gid);
		}
		
		$t2g->{$tid} = $gid;
		$g2t->{$gid} .= $tid."\n";
		$hash->{transcript}{$tid} .= $line."\n";
		$hash->{num_variant}{$gid} += 1;
		$hash->{strand}{$tid} = $A[6];
	}
	elsif($A[2] eq 'CDS'){
		my $tid;
		foreach my $val (@tag){
			if($val =~ /Parent\=/){
				$tid = $val;
				$tid =~ s/Parent\=//;
			}
			if($tid){
				last;
			}
		}
		
		unless($t2g->{$tid}){
			next;
		}
		
		my $gid = $t2g->{$tid};
		$hash->{transcript}{$tid} .= $line."\n";
		$hash->{CDS}{$tid} .= $line."\n";
		$hash->{num_CDS}{$gid} += 1;
		$hash->{num_CDS}{$tid} += 1;
	}
}
close $fh;

my $rh = {};
$rh->{data} = $hash;
$rh->{t2g} = $t2g;
$rh->{g2t} = $g2t;

if(@GID){
	$rh->{GID} = \@GID;
}
if(@GID2){
	$rh->{GID} = \@GID2;
}

return $rh;
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
my $str = shift;

print "! Enter keyword to select files: ";
my $keyword = <STDIN>;
$keyword =~ s/\n//;
$keyword =~ s/\r//;
unless($keyword){
	$keyword = "null";
}

print "! Enter anti-keyword to select files: ";
my $akeyword = <STDIN>;
$akeyword =~ s/\n//;
$akeyword =~ s/\r//;
unless($akeyword){
	$akeyword = "null";
}

my @files2;
my $q1;

my @files = glob("*");
foreach my $file (@files){
	if($file =~ /\.gff/){
		unless($file =~ /\.pl/ || $file =~ /\~/ || $file =~ /\.csv/){
			if($keyword ne 'null' && $akeyword ne 'null'){
				if($file =~ /$keyword/i && $file !~ /$akeyword/i){
					push(@files2, $file);
				}
			}
			elsif($keyword eq 'null' && $akeyword ne 'null'){
				if($file !~ /$akeyword/i){
					push(@files2, $file);
				}
			}
			elsif($keyword ne 'null' && $akeyword eq 'null'){
				if($file =~ /$keyword/i){
					push(@files2, $file);
				}
			}
			else{
				push(@files2, $file);
			}
		}
	}
}

print "\n";
print "--------------------------------------------------------------------------\n";
my $cnt0 = 0;
foreach my $file (@files2){
	print "[$cnt0] $file\n";
	$cnt0++;
}
print "--------------------------------------------------------------------------\n";

unless($q1){
	print "Select $str file (0,1,4 or all): ";
	$q1 = <STDIN>;
	$q1 =~ s/\n//;
	$q1 =~ s/\r//;
}
unless($q1){
	$q1 = 0;
}
elsif($q1 =~ /a/i){
	$q1 = "all";
}

my $file_str = "";
if($q1 eq 'all'){
	print "[all]\n";
	$file_str = join(",", @files2);
}
else{
	my @nums = split(/,/, $q1);
	my @tmp;
	foreach my $n (@nums){
		print "[$files2[$n]]\n";
		push(@tmp, $files2[$n]);
	}
	$file_str = join(",", @tmp);
}

return $file_str;
}


