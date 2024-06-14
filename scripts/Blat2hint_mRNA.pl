#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "Blat2hint_mRNA.pl version 1.01\n";
$version .= "last update: [2019\/2\/10]\n";
$version .= "copyright: ryoichi yano [ryoichiy104\@gmail.com]\n";

#-------------------------------------------------------------------------------

#mRNA seqをgenomeにアライメント -> pslフォーマットで出力 -> augustus trainingに必要なhint.E.gffファイルを生成

my $bin_blat = shift;
my $bin_blat2hint = shift;
my $dbfasta = shift;
my $qfasta = shift;
my $minIdentity = shift;
my $core = shift;
my $q = "y";

my $err = 0;
$err += FileCheck1($bin_blat, "blat");
$err += FileCheck1($bin_blat2hint, "blat2hints.pl");
$err += FileCheck1($dbfasta, "dbfasta");
$err += FileCheck1($qfasta, "qfasta");

if($err > 0){
	print "! abort script due to error...\n";
	goto END;
}

unless($dbfasta){
	$dbfasta = File_select(".fasta", "genome fasta");
}

unless($qfasta){
	$qfasta = File_select(".fasta", "mRNA fasta", $dbfasta);
}

unless($minIdentity){
	print "\nEnter minIdentity value in blat (default 99): ";
	$minIdentity = <STDIN>;
	$minIdentity =~ s/\n//;
	$minIdentity =~ s/\r//;
	unless($minIdentity){
		$minIdentity = 99;
	}
	unless($minIdentity =~ /[0-9]/){
		if($minIdentity =~ /[a-z]/i){
			print "! error: not numerous input [$minIdentity]\n";
			goto END;
		}
	}
	print "! [$minIdentity]\n";
}

unless($core){
	print "\nEnter CPU threads used for each analysis (default 8): ";
	$core = <STDIN>;
	$core =~ s/\n//;
	$core =~ s/\r//;
	unless($core){
		$core = 8;
	}
	unless($core =~ /[0-9]/){
		if($core =~ /[a-z]/i){
			print "! error: not numerous input [$core]\n";
			goto END;
		}
	}
	print "! [$core]\n";
}
if($core > 32){
	$core = 32;
}

my $analysis_log = "";
$analysis_log .= "\n------------------------------------------------------------------------------------\n";
$analysis_log .= "genome fasta             [$dbfasta]\n";
$analysis_log .= "mRNA/EST fasta           [$qfasta]\n";
$analysis_log .= "blat minIdentity         [$minIdentity]\n";
$analysis_log .= "CPU thread               [$core]\n";
$analysis_log .= "------------------------------------------------------------------------------------\n";

unless($q){
	print "\nOK? (Y/n): ";
	$q = <STDIN>;
	$q =~ s/\n//;
	$q =~ s/\r//;
	unless($q){
		$q = "y";
	}
}

if($q =~ /y/i){
	my $prefix_q = $qfasta;
	if($prefix_q =~ /\.fasta/){
		$prefix_q =~ s/\.fasta//;
	}
	elsif($prefix_q =~ /\.fa/){
		$prefix_q =~ s/\.fa//;
	}
	
	my $prefix_db = $dbfasta;
	if($prefix_db =~ /\.fasta/){
		$prefix_db =~ s/\.fasta//;
	}
	elsif($prefix_db =~ /\.fa/){
		$prefix_db =~ s/\.fa//;
	}
	
	my $rdir = "Blat2hint_mRNA_db-".$prefix_db."_q-".$prefix_q;
	unless(-e $rdir){
		system("mkdir $rdir");
	}
	
	#---------------------------------------------------------------------------//
	chdir $rdir;
	unless(-e $qfasta){
		system("ln -s ../$qfasta ./");
	}
	unless(-e $dbfasta){
		system("ln -s ../$dbfasta ./");
	}
	
	my $sdata = Getseqcount($qfasta);
	my $qcnt = $sdata->{cnt};
	my $hash_seq = $sdata->{seq};
	
	my $gdata = Getseqcount($dbfasta);
	my $hash_db = $gdata->{seq};
	
	my $result = "blat_db-".$prefix_db."_q-".$prefix_q.".psl";
#	my $rpsl = "blat_db-".$prefix_db."_q-".$prefix_q."_tophit.psl";
	my $rpsl = "alignment.psl";
	my $hint = "blat_db-".$prefix_db."_q-".$prefix_q."_tophit.E.gff";
	
	print "\n";
	unless(-e $result){
		print "! starting blat db=[$dbfasta] q=[$qfasta]...\n";
		my $splits = Split_query($qfasta, $core, $prefix_q);
		my $j = @{$splits};
		
		my $thrs = [];
		my @eachRs;
		my @SHs;
		for(my $k = 0; $k < $j; $k++){
			my $sh_file = "tmp_q-".$prefix_q.".$k.sh";
			my $each_result = "tmp_q-".$prefix_q.".$k.psl";
			my $sh = "$bin_blat -t=dna -q=rna -noHead -minIdentity=$minIdentity -maxIntron=150000 -minScore=30 $dbfasta $splits->[$k] $each_result >> blat_log.txt 2>&1\n";
			
			if(-e $each_result){
				print "! [$each_result] already exists. skipping...\n";
			}
			else{
				open(my $shfh, ">", $sh_file) or die;
				print $shfh $sh;
				close $shfh;
				
				print "! thread [$k] starting | blat ...\n";
				my $thr = threads->new(\&Bash_exe, $sh_file, $k);
				push(@{$thrs}, $thr);
				push(@eachRs, $each_result);
				push(@SHs, $sh_file);
			}
		}
		foreach my $thr (@{$thrs}){
			$thr->join();
		}
		
		my $cmd = "cat ".join(" ", @eachRs)." > $result";
		system("$cmd");
		print "! output [$result]\n";
		
		foreach my $sh_file (@SHs){
			system("rm $sh_file");
		}
		foreach my $each_result (@eachRs){
			system("rm $each_result");
		}
		foreach my $each_fasta (@{$splits}){
			system("rm $each_fasta");
		}
	}
	else{
		print "! [$result] already exists. skipping blat...\n";
	}
	
	my $strand = Select_tophit($result, $rpsl, $qcnt);
	
	unless(-e $hint){
		my $cmd = "$bin_blat2hint --in=$rpsl --out=$hint --maxgaplen=100";
	#	my $cmd = "$bin_blat2hint --in=$rpsl --out=$hint --remove_redundant=true";
	#	my $cmd = "$bin_blat2hint --in=$rpsl --out=$hint --minintronlen=20";		#did not improve result
		print "\n! cmd=[$cmd]\n";
		system("$cmd");
	}
	else{
		print "! [$hint] already exists. skip blat2hint...\n";
	}
	
	Add_feature_info($hint, $strand, "hint.E.gff", $hash_db);
	
	#---------------------------------------------------------------------------//
}

END:{
	my $end = "";
}


################################################################################
#-------------------------------------------------------------------------------
sub Getseqcount{
my $query = shift;

my $fh;
if($query =~ /\.gz/){
	open($fh, "gzip -dc $query |") or die;
}
else{
	open($fh, "<", $query) or die;
}
my $cnt = 0;
my $ID;
my $seq;
my $hash = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ />/){
		if($seq){
			$hash->{$ID} = $seq;
		}
		$ID = $line;
		$ID =~ s/\>//;
		
		my @tmp = split(/\s/, $ID);
		$ID = $tmp[0];
		
		$seq = "";
		$cnt++;
	}
	else{
		$seq .= $line;
	}
}
close $fh;

if($seq){
	$hash->{$ID} = $seq;
}

my $rhash = {};
$rhash->{cnt} = $cnt;
$rhash->{seq} = $hash;

return $rhash;
}


#-------------------------------------------------------------------------------
sub Add_feature_info{
my $hint = shift;
my $strand = shift;
my $rhint = shift;
my $hash_seq = shift;

print "\n! reading [$hint]...\n";
open(my $fh, "<", $hint) or die;
my $hash = {};
my $ep = {};
my $exon = {};
my $intron = {};
my $nr = {};
my $ch = {};
my $stc = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	my @A = split(/\t/, $line);
	my @tag = split(/\;/, $A[8]);
	
	if($tag[0] =~ /mult/){
		next;
	}
	$tag[0] =~ s/grp\=//;
	
	if($strand->{$tag[0]}){
		$A[6] = $strand->{$tag[0]};
	}
	else{
		print "! error : missing strand information for [$tag[0]] | [$line]\n";
		die;
	}
	
	$hash->{$tag[0]} .= join("\t", @A)."\n";
}
close $fh;

my @ID = keys(%{$hash});
@ID = sort {$a cmp $b} @ID;

foreach my $id (@ID){
#	my $Mpos0;
#	if($hash_seq->{$id}){
#		$Mpos0 = index($hash_seq->{$id}, "ATG");
#		print "[$Mpos0]\n";
#		sleep(1);
#	}
	
	my @L = split(/\n/, $hash->{$id});
	my $AoA = [];
	my $cntex = 0;
	foreach my $line (@L){
		my @A = split(/\t/, $line);
		
		if($A[2] eq 'ep' || $A[2] eq 'exon'){
			$cntex++;
		}
		
		push(@{$AoA}, \@A);
	}
	
	@{$AoA} = sort {$a->[3] <=> $b->[3]} @{$AoA};
	my $numA = @{$AoA};
	
	my $numex = 0;
	for(my $i = 0; $i < $numA; $i++){
		my $A = $AoA->[$i];
#		if($i == 0){
#			if(($A->[3] - 500) > 0){
#				$A->[3] -= 500;
#			}
#		}
#		elsif($i == ($numA - 1)){
#			$A->[4] += 500;
#		}
		
		#---------------------------------------------------//
		#The strategy of this section is not good way... It causes abnormal exons to fit the artificially introduced ATG codon
		my $disable1 = 1;
		if($disable1 eq '0'){
			if($A->[2] eq 'ep' || $A->[2] eq 'exon'){
				$numex++;
				my $chrseq = $hash_seq->{$A->[0]};
				my $d = $A->[4] - $A->[3] + 1;
				my $Mpos0;
				
				if($A->[6] eq '+' && $numex == 1){
					my $prev;
					for(my $i = 0; $i < 100; $i++){
						if($prev){
							$Mpos0 = index($chrseq, "ATG", $prev) + 1;
						}
						else{
							$Mpos0 = index($chrseq, "ATG", $A->[3]) + 1;
						}
						
						if($Mpos0 + 2 < $A->[4]){
							my $Mposg0 = $Mpos0;
							my $Mposg1 = $Mpos0 + 2;
							
		#					print "$id $Mposg0\n";
		#					print "$id $Mposg1\n";
		#					sleep(1);
							
							$stc->{$id} .= "$A->[0]\tb2h\tstart\t$Mposg0\t$Mposg1\t0\t$A->[6]\t\.\tgrp=$id\;pri=4\;src=E\n";
							$prev = $Mpos0;
						}
						else{
							last;
						}
					}
				}
				elsif($A->[6] eq '-' && $numex == $cntex){
					my $prev;
					for(my $i = 0; $i < 100; $i++){
						if($prev){
							$Mpos0 = rindex($chrseq, "CAT", $prev - 2) + 1;
						}
						else{
							$Mpos0 = rindex($chrseq, "CAT", $A->[4]) + 1;
						}
						
						if($Mpos0 + 2 < $A->[4]){
							my $Mposg0 = $Mpos0;
							my $Mposg1 = $Mpos0 + 2;
							
		#					print "$id $Mposg0\n";
		#					print "$id $Mposg1\n";
		#					sleep(1);
							
							$stc->{$id} .= "$A->[0]\tb2h\tstart\t$Mposg0\t$Mposg1\t0\t$A->[6]\t\.\tgrp=$id\;pri=4\;src=E\n";
							$prev = $Mpos0;
						}
						else{
							last;
						}
					}
				}
			}
		}
		#---------------------------------------------------//
		
		if($A->[2] eq 'ep'){
			$A->[2] = "exon";
			$ep->{$id} .= join("\t", @{$A})."\n";
			$nr->{$id} .= $A->[3]."\n";
			$nr->{$id} .= $A->[4]."\n";
			$ch->{$id} = $A->[0];
		}
		elsif($A->[2] eq 'exon'){
			$exon->{$id} .= join("\t", @{$A})."\n";
			$nr->{$id} .= $A->[3]."\n";
			$nr->{$id} .= $A->[4]."\n";
			$ch->{$id} = $A->[0];
		}
		elsif($A->[2] eq 'intron'){
			$intron->{$id} .= join("\t", @{$A})."\n";
			$nr->{$id} .= $A->[3]."\n";
			$nr->{$id} .= $A->[4]."\n";
			$ch->{$id} = $A->[0];
		}
	}
}

my $hmin = {};
my $hmax = {};
foreach my $id (@ID){
	my @Pos = split(/\n/, $nr->{$id});
	my $min;
	my $max;
	foreach my $p (@Pos){
		unless($min && $max){
			$min = $p;
			$max = $p;
		}
		else{
			if($min > $p){
				$min = $p;
			}
			if($max < $p){
				$max = $p;
			}
		}
	}
	
	if($min){
		$hmin->{$id} = $min;
	}
	if($max){
		$hmax->{$id} = $max;
	}
}

my $r;
#my $tss = "transcription_start_site";		# transcription start site		-> cause error
#my $tts = "transcription_end_site";		# transcription terminal site	-> cause error
my $tss = "tss";		# transcription start site 
my $tts = "tts";		# transcription terminal site 
foreach my $id (@ID){
	my $eachr;
	if($strand->{$id} eq '+'){
		if($hmin->{$id}){
			$eachr .= $ch->{$id}."\tb2h\t".$tss."\t".$hmin->{$id}."\t".$hmin->{$id}."\t0\t".$strand->{$id}."\t\.\tgrp\=".$id."\;pri\=4\;src\=E\n";
		}
	}
	elsif($strand->{$id} eq '-'){
		if($hmin->{$id}){
			$eachr .= $ch->{$id}."\tb2h\t".$tts."\t".$hmin->{$id}."\t".$hmin->{$id}."\t0\t".$strand->{$id}."\t\.\tgrp\=".$id."\;pri\=4\;src\=E\n";
		}
	}
	
	if($strand->{$id} eq '+' && $stc->{$id}){
		$eachr .= $stc->{$id};
	}
	
	if($ep->{$id}){
		$eachr .= $ep->{$id};
	}
	if($exon->{$id}){
		$eachr .= $exon->{$id};
	}
	if($intron->{$id}){
		$eachr .= $intron->{$id};
	}
	
	if($strand->{$id} eq '-' && $stc->{$id}){
		$eachr .= $stc->{$id};
	}
	
	if($strand->{$id} eq '+'){
		if($hmax->{$id}){
			$eachr .= $ch->{$id}."\tb2h\t".$tts."\t".$hmax->{$id}."\t".$hmax->{$id}."\t0\t".$strand->{$id}."\t\.\tgrp\=".$id."\;pri\=4\;src\=E\n";
		}
	}
	elsif($strand->{$id} eq '-'){
		if($hmax->{$id}){
			$eachr .= $ch->{$id}."\tb2h\t".$tss."\t".$hmax->{$id}."\t".$hmax->{$id}."\t0\t".$strand->{$id}."\t\.\tgrp\=".$id."\;pri\=4\;src\=E\n";
		}
	}
	
	my $revr = $eachr;
	$revr =~ s/src\=E/src\=M/g;
	
#	$r .= $eachr.$revr;
	$r .= $eachr;
}


open(my $rfh, ">", $rhint) or die;
print $rfh $r;
close $rfh;
print "! output [$rhint]\n";

}


#-------------------------------------------------------------------------------
sub Select_tophit{
my $psl = shift;
my $rpsl = shift;
my $qcnt = shift;

print "\n! reading [$psl]...\n";
open(my $fh, "<", $psl) or die;
my $hash = {};
my $data = {};
my $strand = {};
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	my @A = split(/\t/, $line);
	my $match = $A[0];
	my $mismatch = $A[1];
	my $score = $match - $mismatch;
	
	unless($hash->{$A[9]}){
		$data->{$A[9]} = $line;
		$strand->{$A[9]} = $A[8];
		$hash->{$A[9]} = $score;
	}
	elsif($score > $hash->{$A[9]}){
		$data->{$A[9]} = $line;
		$strand->{$A[9]} = $A[8];
		$hash->{$A[9]} = $score;
	}
	$cnt++;
}
close $fh;

my @ID = keys(%{$data});
@ID = sort {$a cmp $b} @ID;
my $numID = @ID;

print "! total [$qcnt] query sequence\n";
print "! [$numID] selected as tophit from [$cnt] lines\n";

my $map_ratio = sprintf("%.3f", $numID / $qcnt * 100);
print "! alignment ratio = [$map_ratio \%]\n";

my $r;
foreach my $id (@ID){
	$r .= $data->{$id}."\n";
}

open(my $rfh, ">", $rpsl) or die;
print $rfh $r;
close $rfh;
print "! output [$rpsl]\n";

return $strand;
}


#-------------------------------------------------------------------------------
sub Bash_exe{
my $sh_file = shift;
my $num = shift;

system("bash $sh_file\n");
print "! thread [$num] completed.\n";

return $num;
}


#-------------------------------------------------------------------------------
sub Split_query{
my $query = shift;
my $core = shift;
my $qstr = shift;

print "reading [$query]...\n";
my $fh;
if($query =~ /\.gz/){
	open($fh, "gzip -dc $query |") or die;
}
else{
	open($fh, "<", $query) or die;
}
my $cnt = 0;
my $ID;
my $seq;
my @IDs;
my $hash = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ />/){
		unless($cnt == 0){
			push(@IDs, $ID);
			$hash->{$ID} = $seq;
		}
		$ID = $line;
		$ID =~ s/>//;
		$seq = "";
		$cnt++;
	}
	else{
		$seq .= $line;
	}
}
close $fh;

if($seq){
	push(@IDs, $ID);
	$hash->{$ID} = $seq;
}

print "! [$cnt] sequences\n";

my $d = int($cnt / $core) + 1;
print "splitting query sequence into [$core] files...\n";
print "each file will contain [$d] sequences...\n";
my $seg = $d;
my $tmpfa;
my @rfiles;
my $n = 0;
for(my $i = 0; $i < $cnt; $i++){
	if($i == $seg){
		my $fasta = "$qstr.$n.fa";
		open(my $rfh, ">", $fasta) or die;
		print $rfh $tmpfa;
		close $rfh;
#		print "! [$fasta]\n";
		
		push(@rfiles, $fasta);
		
		$n++;
		$tmpfa = "";
		$seg += $d;
	}
	$tmpfa .= ">$IDs[$i]\n$hash->{$IDs[$i]}\n";
}

if($tmpfa){
	my $fasta = "$qstr.$n.fa";
	open(my $rfh, ">", $fasta) or die;
	print $rfh $tmpfa;
	close $rfh;
#	print "! [$fasta]\n";
	
	push(@rfiles, $fasta);
}
#print "! done\n\n";

return \@rfiles;
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


#-------------------------------------------------------------------------------
sub FileCheck1{
my $file = shift;
my $str = shift;

my $err = 0;
if(! $file){
	print "! missing $str ...\n";
	$err++;
}
elsif(! -e $file){
	print "! missing [$file] ...\n";
	$err++;
}

return $err;
}


