#!/usr/bin/perl -w
use strict;
use warnings;
use threads;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "Find_bestCDS.pl version 1.01\n";
$version .= "last update: [2019\/3\/7]\n";
$version .= "copyright: ryoichi yano\n";

#-------------------------------------------------------------------------------

#mRNA sequenceからprotein coding sequenceを探索してレポートする。
#Fix_ATGcodon_gff.plなどで使用。

my $mRNAseq = shift;		#input is not file but sequence string like 'TGAAATTCTCGAATAAGAAGACGTACCCACTTGCTGGAACGCTCAAATTTCTTCACCCCAGAGC'
my $prefix = shift;

unless($mRNAseq){
	goto END;
}

my $ATG = FindATG($mRNAseq);
my $AoT = [];
foreach my $p0 (@{$ATG}){
	my $proteinseq = Translate($mRNAseq, $p0);
	
	my $lnp = length($proteinseq);
#	if($lnp < 20){
	if($lnp < 10){
		next;
	}
	
	if($proteinseq){
		my $p = $p0 + 1;
		
		my @tmp;
		$tmp[0] = $lnp;
		$tmp[1] = $p;
		$tmp[2] = $proteinseq;
		push(@{$AoT}, \@tmp);
	}
}

@{$AoT} = sort {$b->[0] <=> $a->[0]} @{$AoT};

my @S;
my $fasta;
foreach my $tmp (@{$AoT}){
	my $sw = 0;
	foreach my $prev (@S){
		if($prev =~ /$tmp->[2]/i){
			$sw = 1;
			last;
		}
	}
	
#	if($tmp->[0] >= 20 && $sw eq '0'){
	if($tmp->[0] >= 10 && $sw eq '0'){
		my $id = "start".$tmp->[1]."_ln".$tmp->[0];
		
		if($prefix){
			$id = $prefix."_".$id;
		}
		
		$fasta .= ">".$id."\n".$tmp->[2]."\n";
		push(@S, $tmp->[2]);
	}
}

if($fasta){
	print $fasta;
}

END:{
	my $end = 1;
}

################################################################################
#-------------------------------------------------------------------------------
sub Translate{
my $seq = shift;
my $p = shift;

my $cds = substr($seq, $p);
my $ln = length($cds);

my $proteinseq;
for(my $i = 0; $i < $ln; $i += 3){
	my $codon = substr($cds, $i, 3);
	
	my $aa;
	if($codon){
		if($codon =~ /gct/i || $codon =~ /gcc/i || $codon =~ /gca/i || $codon =~ /gcg/i){$aa = "A";}
		elsif($codon =~ /tta/i || $codon =~ /ttg/i || $codon =~ /ctt/i || $codon =~ /ctc/i || $codon =~ /cta/i || $codon =~ /ctg/i){$aa = "L";}
		elsif($codon =~ /cgt/i || $codon =~ /cgc/i || $codon =~ /cga/i || $codon =~ /cgg/i || $codon =~ /aga/i || $codon =~ /agg/i){$aa = "R";}
		elsif($codon =~ /aaa/i || $codon =~ /aag/i){$aa = "K";}
		elsif($codon =~ /aat/i || $codon =~ /aac/i){$aa = "N";}
		elsif($codon =~ /atg/i){$aa = "M";}
		elsif($codon =~ /gat/i || $codon =~ /gac/i){$aa = "D";}
		elsif($codon =~ /ttt/i || $codon =~ /ttc/i){$aa = "F";}
		elsif($codon =~ /tgt/i || $codon =~ /tgc/i){$aa = "C";}
		elsif($codon =~ /cct/i || $codon =~ /ccc/i || $codon =~ /cca/i || $codon =~ /ccg/i){$aa = "P";}
		elsif($codon =~ /caa/i || $codon =~ /cag/i){$aa = "Q";}
		elsif($codon =~ /tct/i || $codon =~ /tcc/i || $codon =~ /tca/i || $codon =~ /tcg/i || $codon =~ /agt/i || $codon =~ /agc/i){$aa = "S";}
		elsif($codon =~ /gaa/i || $codon =~ /gag/i){$aa = "E";}
		elsif($codon =~ /act/i || $codon =~ /acc/i || $codon =~ /aca/i || $codon =~ /acg/i){$aa = "T";}
		elsif($codon =~ /ggt/i || $codon =~ /ggc/i || $codon =~ /gga/i || $codon =~ /ggg/i){$aa = "G";}
		elsif($codon =~ /tgg/i){$aa = "W";}
		elsif($codon =~ /cat/i || $codon =~ /cac/i){$aa = "H";}
		elsif($codon =~ /tat/i || $codon =~ /tac/i){$aa = "Y";}
		elsif($codon =~ /att/i || $codon =~ /atc/i || $codon =~ /ata/i){$aa = "I";}
		elsif($codon =~ /gtt/i || $codon =~ /gtc/i || $codon =~ /gta/i || $codon =~ /gtg/i){$aa = "V";}
		elsif($codon =~ /tag/i || $codon =~ /tga/i || $codon =~ /taa/i){$aa = "*";}
	}
	
	if($aa){
		$proteinseq .= $aa;
		if($aa eq '*'){
			last;
		}
	}
}

return $proteinseq;
}


#-------------------------------------------------------------------------------
sub FindATG{
my $seq = shift;

my $ln = length($seq);
my @Pos;
for(my $i = 0; $i < $ln; $i++){
	my $codon = substr($seq, $i, 3);
	
	if($codon =~ /atg/i){
		push(@Pos, $i);
#		print "$i $ln\n";
	}
}

return \@Pos;
}


