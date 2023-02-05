#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use FindBin;

my $bin_gth = shift;
my $bin_samtools = shift;
my $bin_gffread = shift;
my $bin_blat = shift;
my $bin_blat2hint = shift;
my $bin_matcher = shift;
my $gid = shift;
my $tmp_dbprot = shift;
my $tmp_dbtranscript = shift;
my $tmp_dbcds = shift;
my $tmp_qfasta = shift;
my $tmp_qgff = shift;
my $tmp_qprot = shift;
my $tmp_align = shift;
my $tmpdir = shift;
my $qpref = shift;
my $combined_qprot = shift;

my $script_path = $FindBin::Bin;
my $bin_FindCDS_mRNA = $script_path."/FindCDS_mRNA.pl";
if(! $bin_FindCDS_mRNA || ! -e $bin_FindCDS_mRNA){
	print "! missing FindCDS_mRNA.pl ...\n";
	goto END;
}

if(! $bin_gth || ! -e $bin_gth){
	print "! missing gth...\n";
	goto END;
}
if(! $bin_samtools || ! -e $bin_samtools){
	print "! missing samtools...\n";
	goto END;
}
if(! $bin_gffread || ! -e $bin_gffread){
	print "! missing gffread...\n";
	goto END;
}
if(! $bin_blat || ! -e $bin_blat){
	print "! missing blat...\n";
	goto END;
}
if(! $bin_matcher || ! -e $bin_matcher){
	print "! missing matcher...\n";
	goto END;
}
if(! $tmp_dbprot || ! -e $tmp_dbprot){
	print "! missing dbprot\n";
	goto END;
}
if(! $tmp_dbtranscript || ! -e $tmp_dbtranscript){
	print "! missing dbtranscript\n";
	goto END;
}
if(! $tmp_dbcds || ! -e $tmp_dbcds){
	print "! missing dbcds\n";
	goto END;
}
if(! $tmp_qgff){
	print "! missing qgff file name\n";
	goto END;
}
if(! $tmp_qprot){
	print "! missing qprot file name\n";
	goto END;
}
if(! $tmp_align){
	print "! missing matcher align file name\n";
	goto END;
}
if(! $tmpdir || ! -e $tmpdir){
	print "! missing temp directory\n";
	goto END;
}

if($gid){
	print "\n! findorf.pl | gid=[$gid]\n";
}
if(-e $combined_qprot){
	my $judge_done = Search_previous($combined_qprot, $gid, $qpref);
	
	if($judge_done && $judge_done eq 'true'){
		print "! already done for [$gid], skip...\n";
		goto END;
	}
}

my $cmd = "";
my $minIdentity = 95;
if(-e $tmp_dbprot && -e $tmp_qfasta){
	my $judge = 0;
	if(! -e $tmp_align){
		print "! trying genome threader...\n";
		my @CMD;
		$CMD[0] = "$bin_gth -force -genomic $tmp_qfasta -protein $tmp_dbprot -gff3out -skipalignmentout -o $tmp_qgff > /dev/null 2>&1";
		$CMD[1] = "$bin_gffread $tmp_qgff -g $tmp_qfasta -y $tmp_qprot";
		
		foreach my $cmd (@CMD){
			print "! cmd=[$cmd]\n";
			system($cmd);
		}
		Rm_ifempty($tmp_qprot);
		
		my $cmd3 = "$bin_matcher -asequence $tmp_dbprot -bsequence $tmp_qprot -outfile $tmp_align";
		if(-e $tmp_qprot){
			print "! cmd=[$cmd3]\n";
			system($cmd3);
		}
		if(-e $tmp_align){
			my $hmatcher = Read_matcher($tmp_align);
			if($hmatcher->{score} && $hmatcher->{score} ne '-'){
				print "\n! ORF found in [$tmp_qfasta] with score=[$hmatcher->{score}]\n";
				BindFasta($tmp_qprot, $combined_qprot, $gid, $qpref);
				$judge = 1;
			}
		}
	}
	
	if($judge eq '0'){
		Rm($tmp_align);
		
		print "\n! trying blat...\n";
		my $psl_transcript = "$tmpdir/_qtranscript_".$gid.".psl";
		my $egtf_transcript = "$tmpdir/_qtranscript_".$gid.".E.gtf";
		my $egff_transcript = "$tmpdir/_qtranscript_".$gid.".E.gff";
		my $eprot_transcript = "$tmpdir/_qtranscript_".$gid.".E.fasta";
		my $prefix = $gid.".E.transcript";
		
		Rm($egtf_transcript);
		Rm($egff_transcript);
		
		my @CMD1;
		$CMD1[0] = "$bin_blat -t=dna -q=rna -noHead -minIdentity=$minIdentity -maxIntron=10000 -minScore=30 $tmp_qfasta $tmp_dbtranscript $psl_transcript";
		$CMD1[1] = "$bin_blat2hint --in=$psl_transcript --out=$egtf_transcript --ep_cutoff=0";
		foreach my $cmd (@CMD1){
			print "! cmd=[$cmd]\n";
			system($cmd);
		}
		
		my $strand = Select_tophit($psl_transcript, $psl_transcript);
		if($strand && $strand ne 'null'){
			Add_feature_info($egtf_transcript, $strand, $egff_transcript);
			
			my @CMD2;
			$CMD2[0] = "$bin_samtools faidx $tmp_qfasta";
			$CMD2[1] = "$bin_gffread $egff_transcript -g $tmp_qfasta -w $eprot_transcript";
			foreach my $cmd (@CMD2){
				print "! cmd=[$cmd]\n";
				system($cmd);
			}
			
			my $rh = Getseqinfo($eprot_transcript);
			my $hseq = $rh->{data};
			my @SeqID = keys(%{$hseq});
			@SeqID = sort {$a cmp $b} @SeqID;
			foreach my $sid (@SeqID){
				my $cmd = "perl $bin_FindCDS_mRNA $hseq->{$sid}{seq} $prefix > $tmp_qprot";
	#			print "! cmd=[$cmd]\n";
				system($cmd);
				
				if(-e $tmp_qprot){
					my $cmd3 = "$bin_matcher -asequence $tmp_dbprot -bsequence $tmp_qprot -outfile $tmp_align";
					if(-e $tmp_qprot){
						print "! cmd=[$cmd3]\n";
						system($cmd3);
					}
					if(-e $tmp_align){
						my $hmatcher = Read_matcher($tmp_align);
						if($hmatcher->{score} && $hmatcher->{score} ne '-'){
							print "\n! ORF found in [$tmp_qfasta] with score=[$hmatcher->{score}]\n";
							BindFasta($tmp_qprot, $combined_qprot, $gid, $qpref);
							last;
						}
					}
				}
			}
		}
	}
	
	my $keep = 0;
	if($keep eq '0'){
		Rm($tmp_dbprot);
		Rm("$tmp_dbprot.md5");
		Rm("$tmp_dbprot.suf");
		Rm($tmp_qfasta);
		Rm("$tmp_qfasta.fai");
		Rm($tmp_qgff);
	}
	if(-e "$tmp_dbprot.protein.bck"){
		system("rm $tmp_dbprot.protein.*");
	}
	if(-e "$tmp_qfasta.dna.tis"){
		system("rm $tmp_qfasta.dna.*");
	}
}

END:{
	print "\n";
}


################################################################################
#-------------------------------------------------------------------------------
sub Search_previous{
my ($combined_qprot, $gid, $qpref) = @_;

open(my $fh, "<", $combined_qprot);
my $judge = "false";
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /\>/){
		my $ID = $line;
		$ID =~ s/\>//;
		my @tmp = split(/_/, $ID);
		if($tmp[0] && $tmp[1] && $tmp[0] eq $qpref && $tmp[1] eq $gid){
			$judge = "true";
			last;
		}
	}
}
close $fh;

return $judge;
}


#-------------------------------------------------------------------------------
sub BindFasta{
my ($tmp_qprot, $combined_qprot, $gid, $qpref) = @_;

open(my $fh, "<", $tmp_qprot);
my $ID;
my $seq;
my $r;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /\>/){
		if($seq){
			$r .= ">".$qpref."_".$gid."_".$ID."\n".$seq."\n";
		}
		
		$ID = $line;
		$ID =~ s/\>//;
		my @tmp = split(/\s/, $ID);
		$ID = $tmp[0];
		$seq = "";
	}
	else{
		$seq .= $line;
	}
}
close $fh;

if($seq){
	$r .= ">".$qpref."_".$gid."_".$ID."\n".$seq."\n";
}

open(my $rfh, ">>", $combined_qprot);
print $rfh $r;
close $rfh;

}


#-------------------------------------------------------------------------------
sub Read_matcher{
my $file = shift;

open(my $fh, "<", $file);
my $hash = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	# Aligned_sequences: 2
	# 1: MELO.jh000002.1.t1
	# 2: mRNA1
	# Matrix: EBLOSUM62
	# Gap_penalty: 14
	# Extend_penalty: 4
	#
	# Length: 282
	# Identity:     282/282 (100.0%)
	# Similarity:   282/282 (100.0%)
	# Gaps:           0/282 ( 0.0%)
	# Score: 1491
	
	my @A;
	my $sw = 0;
	if($line =~ /Length:/){
		@A = split(/Length:/, $line);
		$A[1] =~ s/\s//g;
		
		if($A[1]){
			$hash->{len} = $A[1];
		}
		else{
			$hash->{len} = "-";
		}
	}
	elsif($line =~ /Identity:/){
		@A = split(/Identity:/, $line);
		$A[1] =~ s/\s//g;
		$A[1] =~ s/\)//g;
		$A[1] =~ s/\%//g;
		
		my @B = split(/\(/, $A[1]);
		if($B[0]){
			$hash->{identity_str} = $B[0];
		}
		else{
			$hash->{identity_str} = "-";
		}
		if($B[1]){
			$hash->{identity} = $B[1];
		}
		else{
			$hash->{identity} = "-";
		}
	}
	elsif($line =~ /Similarity:/){
		@A = split(/Similarity:/, $line);
		$A[1] =~ s/\s//g;
		$A[1] =~ s/\)//g;
		$A[1] =~ s/\%//g;
		
		my @B = split(/\(/, $A[1]);
		if($B[0]){
			$hash->{similarity_str} = $B[0];
		}
		else{
			$hash->{similarity_str} = "-";
		}
		if($B[1]){
			$hash->{similarity} = $B[1];
		}
		else{
			$hash->{similarity} = "-";
		}
	}
	elsif($line =~ /Gaps:/){
		@A = split(/Gaps:/, $line);
		$A[1] =~ s/\s//g;
		$A[1] =~ s/\)//g;
		$A[1] =~ s/\%//g;
		
		my @B = split(/\(/, $A[1]);
		if($B[0]){
			$hash->{gap_str} = $B[0];
		}
		else{
			$hash->{gap_str} = "-";
		}
		if($B[1]){
			$hash->{gap} = $B[1];
		}
		else{
			$hash->{gap} = "-";
		}
	}
	elsif($line =~ /Score:/){
		@A = split(/Score:/, $line);
		$A[1] =~ s/\s//g;
		
		if($A[1]){
			$hash->{score} = $A[1];
		}
		else{
			$hash->{score} = "-";
		}
	}
}
close $fh;

return $hash;
}


#-------------------------------------------------------------------------------
sub Getseqinfo{
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
my $Mstart = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	if($line =~ />/){
		if($seq){
			$hash->{$ID}{seq} = $seq;
			$hash->{$ID}{len} = length($seq);
			
			my $aa1 = substr($seq, 0, 1);
			if($aa1 =~ /m/i){
				$Mstart->{$ID} = 1;
			}
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
	$hash->{$ID}{seq} = $seq;
	$hash->{$ID}{len} = length($seq);
	
	my $aa1 = substr($seq, 0, 1);
	if($aa1 =~ /m/i){
		$Mstart->{$ID} = 1;
	}
}

my $rhash = {};
$rhash->{cnt} = $cnt;
$rhash->{data} = $hash;
$rhash->{Mstart} = $Mstart;

return $rhash;
}


#-------------------------------------------------------------------------------
sub Add_feature_info{
my $hint = shift;
my $strand = shift;
my $rhint = shift;

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

}


#-------------------------------------------------------------------------------
sub Select_tophit{
my $psl = shift;
my $rpsl = shift;

print "\n! reading [$psl]...\n";
open(my $fh, "<", $psl) or return 'null';
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

my $r;
foreach my $id (@ID){
	$r .= $data->{$id}."\n";
}

open(my $rfh, ">", $rpsl) or die;
print $rfh $r;
close $rfh;

return $strand;
}


#-------------------------------------------------------------------------------
sub Rm_ifempty{
my $file = shift;

if(-e $file){
	open(my $fh, "<", $file) or die;
	my $dat = "";
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line){
			$dat .= $line."\n";
		}
	}
	close $fh;
	
	if(! $dat){
		system("rm $file");		# empty file
	}
}

}


#-------------------------------------------------------------------------------
sub Rm{
my $file = shift;

if(-e $file){
	system("rm $file");
}

}

