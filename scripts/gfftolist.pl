#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;

my $gff3 = shift;
my $rcsv = shift;

if(! $gff3 || ! -e $gff3){
	print "! error : missing input GFF3...\n";
	goto END;
}
if(! $rcsv || $gff3 eq $rcsv){
	$rcsv = $gff3;
	$rcsv =~ s/\.gff3//;
	$rcsv =~ s/\.gff//;
	$rcsv = "summary_gene_".$rcsv.".csv";
}

my $rh = Open_gff($gff3);

SAVE($rcsv, $rh->{r});
#my $rt2g = "list_t2g_".$gff3.".tsv";
#SAVE($rt2g, $rh->{t2g_lines});

END:{
	my $end = 1;
}



################################################################################
#-------------------------------------------------------------------------------
sub Open_gff{
my $gff3 = shift;

print "! reading gff3...\n";
open(my $fh, "<", $gff3) or die;
my $t2g = {};
my $t2g_lines = "";
my $hash = {};
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
		
		$hash->{$gid}{data} = $gid.",".$A[0].",".$A[3].",".$A[4].",".$A[6];
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
		$hash->{$gid}{num_variant} += 1;
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
			print "! missing gid for [$tid]\n";
			sleep(1);
		}
		
		my $gid = $t2g->{$tid};
		$hash->{$gid}{num_CDS} += 1;
	}
}

print " [$gcnt] genes\n";

my $r = "gid,chr,pos0,pos1,strand,num_variant,total_num_CDS\n";
my $numpc = 0;
my $numnc = 0;
my $numpc_tr = 0;
my $hcnt = {};
foreach my $gid (@GID){
	unless($hash->{$gid}{num_variant}){
		$hash->{$gid}{num_variant} = 0;
	}
	unless($hash->{$gid}{num_CDS}){
		$hash->{$gid}{num_CDS} = 0;
		$numnc++;
	}
	else{
		$numpc++;
		$numpc_tr += $hash->{$gid}{num_variant};
	}
	
	my $tmp = $hash->{$gid}{data}.",".$hash->{$gid}{num_variant}.",".$hash->{$gid}{num_CDS};
	my @B = split(/,/, $tmp);
	$hcnt->{$B[1]} += 1;
	$r .= $tmp."\n";
}

print " [$numpc] protein-coding, [$numpc_tr] transcripts\n";
print " [$numnc] non-coding\n";

my @Chrs = keys(%{$hcnt});
@Chrs = sort {$a cmp $b} @Chrs;

$r .= "\nChr,num gene\n";
foreach my $sid (@Chrs){
	$r .= $sid.",".$hcnt->{$sid}."\n";
}

my $rh = {};
$rh->{r} = $r;
$rh->{pc} = $numpc;
$rh->{tr} = $numpc_tr;
$rh->{nc} = $numnc;
$rh->{t2g_lines} = $t2g_lines;

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


