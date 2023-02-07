#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "psl2dot.pl version 1.02\n";
$version .= "last update: [2020\/1\/19]\n";
$version .= "copyright: ryoichi yano\n";

#print $version;

#-------------------------------------------------------------------------------

my $psl = shift;
my $lth = 100;
my $q = "pipe";
my $allq = 2;			#use all query sequences if set to 1. default = 0
my $lnq = 1000;			#threshold for query seq length
my $bpresolution = 100;
my $dir = shift;
my $cmmap = "null";

#--------------------------------------------------------------//
my $cwd = getcwd();
my $script_log = "";
unless($script_log){
	$script_log = "/media/hdd0/script_logs";
	unless(-e $script_log){
		$script_log = getcwd();
	}
}
#--------------------------------------------------------------//

#EST hintがオーバーラップする遺伝子や、元々は単一遺伝子であるはずがタンデムにスプリットされてアノテーションされた遺伝子を修正する。Search_confident_pep.pl の中で実行。

# linkage map input (output of AnchorScaff_xxx.pl)
# seqID	pos		marker_name		LG	cM
# chr01	45230	T000100097l		1	0.2754932
# chr01	45210	T000100097m		1	0.2754932
# chr01	45104	T000100097s		1	0.2754932
#
# interpretation of LG (linkage group) : LG = 1 -> chr01 or Chr01 or Gm01

if(! $q || $q ne 'pipe'){
	print $version;
}

unless($psl){
	$psl = File_select(".psl", "psl alignment file by LAST");
}

if(! $cmmap && ! $q){
	$cmmap = File_select_marker(".csv", "overlay linkage map data (output of \"AnchorScaff\" script, can be absent)");
}
else{
	$cmmap = "null";
}

unless($lth){
	print "\nEnter threshold for the length of each alignment block (default 2000 bp): ";
	$lth = <STDIN>;
	$lth =~ s/\n//;
	$lth =~ s/\r//;
	unless($lth){
		$lth = 2000;
	}
	if($lth =~ /\D/){
		print "! error: not numerous input [$lth]\n";
		goto END;
	}
	print "! [$lth] bp\n";
}

unless($q){
	unless($allq){
		print "\n----------------------------------------------------------------------\n";
		print "[0] query = chr and target = chr\n";
		print "[1] query = any and target = chr\n";
		print "[2] query = any and target = any\n";
		print "----------------------------------------------------------------------\n";
		print "select query strategy (default 2): ";
		$allq = <STDIN>;
		$allq =~ s/\n//;
		$allq =~ s/\r//;
		my $tmpq = $allq."t";
		if($tmpq eq 't'){
			$allq = 2;
		}
		if($allq ne '0' && $allq ne '1' && $allq ne '2'){
			$allq = 2;
		}
		print "! [$allq]\n";
	}
	
	if($allq ne '0' && ! $lnq){
		print "\nEnter threshold for query seq length (default 10 kb): ";
		$lnq = <STDIN>;
		$lnq =~ s/\n//;
		$lnq =~ s/\r//;
		unless($lnq){
			$lnq = 10;
		}
		if($lnq =~ /\D/){
			$lnq = 10;
		}
		$lnq *= 1000;
		print "! [$lnq] bp\n";
	}
	else{
		$lnq = 0;
	}
}
else{
	unless($allq){
		$allq = 0;
	}
	unless($lnq){
		$lnq = 0;
	}
}

unless($bpresolution){
	$bpresolution = 100;
}

my $analysis_log = "";
$analysis_log .= "\n------------------------------------------------------------------------------------\n";
$analysis_log .= "psl alignment file by LAST   [$psl]\n";
$analysis_log .= "overlay linkage map data     [$cmmap]\n";
$analysis_log .= "block length threshold (bp)  [$lth]\n";
$analysis_log .= "bp resolution in report      [$bpresolution]\n";
if($allq ne '0'){
	$analysis_log .= "query length threshold (bp)  [$lnq]\n";
}
$analysis_log .= "------------------------------------------------------------------------------------\n";

print $analysis_log,"\n";

my $pref = $psl;
$pref =~ s/\.psl//;
if($pref =~ /\//){
	my @PREF = split(/\//, $pref);
	my $nPREF = @PREF;
	$pref = $PREF[$nPREF - 1];
}

unless($dir){
	$dir = "_pslalign_lth".$lth."_".$pref;
}
#print "! cmd=[Psl_createPlot.pl $psl $lth y $allq $lnq $bpresolution $dir $cmmap]\n";

unless($q){
	print "\nOK? (Y/n): ";
	$q = <STDIN>;
	$q =~ s/\n//;
	$q =~ s/\r//;
	unless($q){
		$q = "y";
	}
}
#elsif($q =~ /pipe/i){
#	print $analysis_log,"\n";
#}

if($q =~ /y/i || $q =~ /pipe/i){
	if($q =~ /y/i){
		print "\n";
	}
	
	my $rh = Convert_to_cumpos($psl, $lth, $lnq, $allq);
#	$rh->{QID} = \@QID;
#	$rh->{TID} = \@TID;
#	$rh->{qposFw} = \@qposFw;
#	$rh->{qposRv} = \@qposRv;
#	$rh->{tposFw} = \@tposFw;
#	$rh->{tposRv} = \@tposRv;
#	$rh->{stats_log} = $stats_log;
#	$rh->{stats_count} = $stats_count;
#	$rh->{stats_pos} = $stats_pos;
#	$rh->{stats_hmatches} = $stats_hmatches;
#	$rh->{stats_hmisMatches} = $stats_hmisMatches;
#	$rh->{stats_hblockCount} = $stats_hblockCount;
	
	my $hcmmap = {};
	if($cmmap ne 'null' && -e $cmmap){
		$hcmmap = Read_AnchorScaff_output($cmmap);
	}
	
	unless(-e $dir){
		system("mkdir $dir");
	}
	chdir $dir;
	$cwd .= "/".$dir;
	
	SAVE("stats_log.txt", $rh->{stats_log});
	SAVE("stats_count.tsv", $rh->{stats_count});
	SAVE("stats_seqlength.tsv", $rh->{stats_pos});
	SAVE("stats_hmatches.tsv", $rh->{stats_hmatches});
	SAVE("stats_hmisMatches.tsv", $rh->{stats_hmisMatches});
	SAVE("stats_hblockCount.tsv", $rh->{stats_hblockCount});
	
	my $hposdata = $rh->{hposdata};
	my @ID = keys(%{$hposdata});
	@ID = sort {$a cmp $b} @ID;
	
	print "! making dot files...\n";
	my $cnt_plot = 0;
	my $R_script =<<"EOS";
workingDir = "$cwd"
setwd(workingDir)
getwd()

EOS
	
	my $summary = "query seq_id,length,qStart,qEnd,target (DB) seq_id,length,tStart,tEnd,matches,num block,ratio\n";
	foreach my $id (@ID){
		if($id eq 'All'){
			my $sh = R_plot($hposdata->{$id}, $id, $lth, $pref, $bpresolution, "All");
			$R_script .= $sh->{script};
			$cnt_plot++;
		}
		else{
			my $each_hposdata = $hposdata->{$id};
			my @eachID = keys(%{$each_hposdata});
			@eachID = sort {$a cmp $b} @eachID;
			
			my $th_alen = 100 * 1000;
			my $th_ratio = 0.6;
			foreach my $tid (@eachID){
				if($tid =~ /chloro/i || $tid =~ /mitoch/i){
					$th_alen = 50 * 1000;
					$th_ratio = 0.4;
				}
				
				my $ratio = $hposdata->{$id}{$tid}{matches} / $hposdata->{$id}{$tid}{qcumpos};
				if($hposdata->{$id}{$tid}{qcumpos} >= $th_alen && $ratio >= $th_ratio){
					my $sh = R_plot($hposdata->{$id}{$tid}, $id, $lth, $pref, $bpresolution, $tid, $hcmmap);
					$R_script .= $sh->{script};
					$summary .= "$id,$hposdata->{$id}{$tid}{qcumpos},$sh->{qpos0},$sh->{qpos1},$tid,$hposdata->{$id}{$tid}{tcumpos},$sh->{tpos0},$sh->{tpos1},$hposdata->{$id}{$tid}{matches},$hposdata->{$id}{$tid}{nblock},$ratio\n";
					$cnt_plot++;
				}
			}
		}
	}
	print "! done.\n";
	
	SAVE("stats_alignment_length.csv", $summary);
}

END:{
	if($q && $q !~ /pipe/i){
		print "\n! End of script.\n\n";
	}
}



################################################################################
#-------------------------------------------------------------------------------
sub Read_AnchorScaff_output{
my $file = shift;

print "! reading AnchorScaff output [$file]...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	my @A = split(/,/, $line);
	my $numA = @A;
	
	if($numA < 5){
		next;
	}
	
	# linkage map input (output of AnchorScaff_xxx.pl)
	# seqID	pos		marker_name		LG	cM
	# chr01	45230	T000100097l		1	0.2754932
	# chr01	45210	T000100097m		1	0.2754932
	# chr01	45104	T000100097s		1	0.2754932
	#
	# interpretation of LG (linkage group) : LG = 1 -> chr01 or Chr01 or Gm01
	
	my $seqNo = $A[0];
	$seqNo =~ s/chr0//i;
	$seqNo =~ s/chr//i;
	$seqNo =~ s/Gm0//;
	$seqNo =~ s/Gm//;
	
	if($seqNo eq $A[3]){
		my $marker = $A[2];
#		$marker =~ s/l//;
#		$marker =~ s/m//;
#		$marker =~ s/s//;
		
		$hash->{$A[0]}{$marker}{pos} += $A[1];
		$hash->{$A[0]}{$marker}{cnt} += 1;
		$hash->{$A[0]}{$marker}{cM} = $A[4];
		
		$cnt++;
	}
}
close $fh;

if($cnt > 0){
	my @SID = keys(%{$hash});
	@SID = sort {$a cmp $b} @SID;
	
	my $cnt_selected = 0;
	my $rh = {};
	foreach my $sid (@SID){
		my $htmp = $hash->{$sid};
		my @MK = keys(%{$htmp});
		@MK = sort {$a cmp $b} @MK;
		
		my @CM;
		foreach my $marker (@MK){
#			if($hash->{$sid}{$marker}{cnt} && $hash->{$sid}{$marker}{cnt} == 3){
			if($hash->{$sid}{$marker}{cnt} && $hash->{$sid}{$marker}{cnt} > 0){
				my $avrpos = int($hash->{$sid}{$marker}{pos} / $hash->{$sid}{$marker}{cnt});
				push(@CM, $hash->{$sid}{$marker}{cM});
				
				$rh->{$sid}{data} .= $marker."\t".$avrpos."\t".$hash->{$sid}{$marker}{cM}."\t".$sid."\n";
				$cnt_selected++;
			}
		}
		
		@CM = sort {$b <=> $a} @CM;
		$rh->{$sid}{maxcM} = $CM[0];
	}
	
	print "! [$cnt_selected] markers selected\n";
	return $rh;
}
else{
	print "! no marker data found\n";
}

}


#-------------------------------------------------------------------------------
sub R_plot{
my $rh = shift;
my $id = shift;
my $lth = shift;
my $pref = shift;
my $bpresolution = shift;
my $string = shift;
my $hcmmap = shift;

my $qcumpos = $rh->{qcumpos};
my $tcumpos = $rh->{tcumpos};
my $seqlen_str = $rh->{stats_pos};

my $qposFw = $rh->{qposFw};
my $qposRv = $rh->{qposRv};
my $tposFw = $rh->{tposFw};			# t = db fasta in lastal
my $tposRv = $rh->{tposRv};			# t = db fasta in lastal
my $bsizeFw = $rh->{bsizeFw};
my $bsizeRv = $rh->{bsizeRv};

my $num_qposFw = @{$qposFw};
my $num_qposRv = @{$qposRv};
my $num_tposFw = @{$tposFw};
my $num_tposRv = @{$tposRv};

unless($num_qposFw == $num_tposFw){
	print "\n! error | counts unmatched between query and target (Fw) : [$num_qposFw] [$num_tposFw]\n";
	return;
}
unless($num_qposRv == $num_tposRv){
	print "\n! error | counts unmatched between query and target (Rv) : [$num_qposRv] [$num_tposRv]\n";
	return;
}

my $dirmap = "linkage_map";
unless(-e $dirmap){
	system("mkdir $dirmap");
}

my $pfile_Fw;
my $pfile_Rv;
my $png;
my $png2;
my $png_map;
my $lgfile;
my $img_width = 3000;
my $img_height = 3000;
my $lwd = 1;
my $cex = 1;
if($string eq 'All'){
	$pfile_Fw = $string."_cumpos_Fw_".$pref.".tsv";
	$pfile_Rv = $string."_cumpos_Rv_".$pref.".tsv";
	$png = "plot_lth".$lth."_".$string."_".$pref.".png";
	$img_width = 3000;
	$img_height = 3000;
	$lwd = 1;
	$cex = 1;
}
else{
#	$pfile_Fw = $id."_vs_".$string."_cumpos_Fw_".$pref.".tsv";
#	$pfile_Rv = $id."_vs_".$string."_cumpos_Rv_".$pref.".tsv";
#	$png = "plot_".$id."_vs_".$string."_".$pref.".png";
	$pfile_Fw = $string."_vs_".$id."_cumpos_Fw_".$pref.".tsv";
	$pfile_Rv = $string."_vs_".$id."_cumpos_Rv_".$pref.".tsv";
	$png = "plot_lth".$lth."_".$string."_vs_".$id."_".$pref.".png";
	$png2 = "mplot_lth".$lth."_".$string."_vs_".$id."_".$pref.".png";
	$png_map = $dirmap."/map_".$string."_vs_".$id."_".$pref.".png";
	$img_width = 3000;
	$img_height = 3000;
	$lwd = 1;
	$cex = 2;
	
	if($hcmmap->{$id}{data}){
		$lgfile = $dirmap."/".$id."_linkage_map.tsv";
		$hcmmap->{$id}{data} = "marker\tposition\tcM\tseqID\n".$hcmmap->{$id}{data};
		
		open(my $lfh, ">", $lgfile);
		print $lfh $hcmmap->{$id}{data};
		close $lfh;
	}
}

my @QPOS;
my @TPOS;
if($qposFw && $tposFw){
	my $P0 = $qposFw;
	my $P1 = $tposFw;		# P1 = db fasta in lastal
	my $Bs = $bsizeFw;
	
	my $r = "query [x-axis]\ttarget (DB) [y-axis]\tblocksize\n";
	for(my $i = 0; $i < $num_qposFw; $i++){
#		$r .= $P0->[$i]."\t".$P1->[$i]."\t".$Bs->[$i]."\n";
		
		for(my $p = 0; $p <= $Bs->[$i]; $p += $bpresolution){
			my $dotpos0 = $P0->[$i] + $p;
			my $dotpos1 = $P1->[$i] + $p;
			$r .= $dotpos0."\t".$dotpos1."\t".$bpresolution."\n";
		}
		
		push(@QPOS, $P0->[$i]);
		push(@TPOS, $P1->[$i]);
	}

	open(my $rfh, ">", $pfile_Fw);
	print $rfh $r;
	close $rfh;
#	print "! output [$rfile]\n";
}

if($qposRv && $tposRv){
	my $P0 = $qposRv;
	my $P1 = $tposRv;
	my $Bs = $bsizeRv;
	
	my $r = "query [x-axis]\ttarget (DB) [y-axis]\tblocksize\n";
	for(my $i = 0; $i < $num_qposRv; $i++){
#		$r .= $P0->[$i]."\t".$P1->[$i]."\t".$Bs->[$i]."\n";
		
		for(my $p = 0; $p <= $Bs->[$i]; $p += $bpresolution){
			my $dotpos0 = $P0->[$i] - $p;
			my $dotpos1 = $P1->[$i] + $p;
			$r .= $dotpos0."\t".$dotpos1."\t".$bpresolution."\n";
		}
		
		push(@QPOS, $P0->[$i]);
		push(@TPOS, $P1->[$i]);
	}

	open(my $rfh, ">", $pfile_Rv);
	print $rfh $r;
	close $rfh;
#	print "! output [$rfile]\n";
}

my @SL = split(/\n/, $seqlen_str);
my $script_line;
foreach my $line (@SL){
	my @A = split(/\t/, $line);
	if($A[0] eq 'query'){
		$script_line .= "abline(v\=$A[3], col\=\"black\", lwd=".$lwd.")\n";
	}
	elsif($A[0] eq 'target'){
		$script_line .= "abline(h\=$A[3], col\=\"black\", lwd=".$lwd.")\n";
	}
}

#$script_line .= "abline(v\=$qcumpos, col\=\"black\")\n";
#$script_line .= "abline(h\=$tcumpos, col\=\"black\")\n";

my $subline_dist = 5000 * 1000;
my $subq = abs($qcumpos - 0) / $subline_dist;
my $subt = abs($tcumpos - 0) / $subline_dist;

unless($id =~ /all/i){
	for(my $i = 0; $i < $subq; $i++){
		my $subline_pos = $i * $subline_dist;
		$script_line .= "abline(v\=$subline_pos, col\=\"black\", lwd=".$lwd.")\n";
	}
	for(my $i = 0; $i < $subt; $i++){
		my $subline_pos = $i * $subline_dist;
		$script_line .= "abline(h\=$subline_pos, col\=\"black\", lwd=".$lwd.")\n";
	}
}

my $char_size = ", cex.axis=3, cex.lab=3";
if($id =~ /all/i){
	$char_size = "cex.axis=1, cex.lab=3";
}
else{
	$char_size = "cex.axis=3, cex.lab=3";
}

my $ylab = "";
my $xlab = "";
if($string eq 'All'){
	$ylab = "target (DB)";
	$xlab = "query";
}
else{
	$ylab = "target (DB) ".$string;
	$xlab = "query ".$id;
}

my $script;
if($lgfile && -e $lgfile){
	$script =<<"EOS";
fw <- read.table("$pfile_Fw", sep="\t", header=T)
rv <- read.table("$pfile_Rv", sep="\t", header=T)
lg <- read.table("$lgfile", sep="\t", header=T, row.names=1)

png("$png2", width=$img_width, height=$img_height)
plot(fw[,1], fw[,2], pch=20, col="blue", xlim=c(0, $qcumpos), ylim=c(0, $tcumpos), $char_size, ylab="$ylab", xlab="$xlab")
points(rv[,1], rv[,2], pch=20, col="red", xlim=c(0, $qcumpos), ylim=c(0, $tcumpos), $char_size)
$script_line
par(new = TRUE)
plot(lg[,1], lg[,2], pch=17, col="darkgreen", xlim=c(0, $qcumpos), ylim=c(0, $hcmmap->{$id}{maxcM}), $char_size, ylab="", xlab="", axes = FALSE, cex=4)
axis(4)
dev.off()

png("$png", width=$img_width, height=$img_height)
plot(fw[,1], fw[,2], pch=20, col="blue", xlim=c(0, $qcumpos), ylim=c(0, $tcumpos), $char_size, ylab="$ylab", xlab="$xlab")
points(rv[,1], rv[,2], pch=20, col="red", xlim=c(0, $qcumpos), ylim=c(0, $tcumpos), $char_size)
$script_line
dev.off()

png("$png_map", width=$img_width, height=$img_height)
plot(lg[,1], lg[,2], pch=17, col="darkgreen", xlim=c(0, $qcumpos), ylim=c(0, $hcmmap->{$id}{maxcM}), $char_size, ylab="linkage position (cM)", xlab="physical position (bp)", cex=4)
dev.off()

EOS
}
else{
	$script =<<"EOS";
fw <- read.table("$pfile_Fw", sep="\t", header=T)
rv <- read.table("$pfile_Rv", sep="\t", header=T)

png("$png", width=$img_width, height=$img_height)
plot(fw[,1], fw[,2], pch=20, col="blue", xlim=c(0, $qcumpos), ylim=c(0, $tcumpos), $char_size, ylab="$ylab", xlab="$xlab")
points(rv[,1], rv[,2], pch=20, col="red", xlim=c(0, $qcumpos), ylim=c(0, $tcumpos), $char_size)
$script_line
dev.off()

EOS
}

my $sh = {};
$sh->{script} = $script;

@QPOS = sort {$a <=> $b} @QPOS;
@TPOS = sort {$a <=> $b} @TPOS;
$sh->{qpos0} = $QPOS[0];
$sh->{tpos0} = $TPOS[0];

@QPOS = sort {$b <=> $a} @QPOS;
@TPOS = sort {$b <=> $a} @TPOS;
$sh->{qpos1} = $QPOS[0];
$sh->{tpos1} = $TPOS[0];

return $sh;
}


#----------------------------------------------------------
sub Convert_to_cumpos{
my $file = shift;
my $lth = shift;
my $lnq = shift;
my $allq = shift;

print "! reading psl [$file]...\n";
open(my $fh, "<", $file) or die;
my $hstat = {};
my $hqsize = {};
my $htsize = {};
my $hpos = {};
my $hstrand = {};
my $hmatches = {};
my $hmisMatches = {};
my $hblockCount = {};
my $hq_tid = {};
my $hq_tpos = {};
my $hallen = {};
my $cnt = 0;
my $cnt_selected = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
#	$line =~ s/unanchored_scaffolds/chr00/;
	
	my @A = split(/\t/, $line);
	
	my $matches = $A[0];
	my $misMatches = $A[1];
	my $repMatches = $A[2];
	my $nCount = $A[3];
	my $qNumInsert = $A[4];
	my $qBaseInsert = $A[5];
	my $tNumInsert = $A[6];
	my $tBaseInsert = $A[7];
	my $strand = $A[8];
	my $qName = $A[9];
	my $qSize = $A[10];
	my $qStart = $A[11];
	my $qEnd = $A[12];
	my $tName = $A[13];				# DB fasta in lastal
	my $tSize = $A[14];				# DB fasta in lastal
	my $tStart = $A[15];			# DB fasta in lastal
	my $tEnd = $A[16];				# DB fasta in lastal
	my $blockCount = $A[17];
	my $blockSizes = $A[18];
	my $qStarts = $A[19];
	my $tStarts = $A[20];
	
	my $sw = 0;
	if($qName =~ /mitochon/i || $tName =~ /mitochon/i){
		$sw = 0;
	}
	elsif($qName =~ /chloro/i || $tName =~ /chloro/i){
		$sw = 0;
	}
	else{
		if($allq eq '0' && $qName =~ /chr/i && $tName =~ /chr/i){
			$sw = 1;
		}
		elsif($allq eq '0' && $qName =~ /chr/ && $tName =~ /Gm/ && $tName !~ /scaffold/i){
			$sw = 1;
		}
		else{
			if($allq eq '1' && $qSize >= $lnq && $tName =~ /chr/i){
				$sw = 1;
			}
			elsif($allq eq '1' && $qSize >= $lnq && $tName =~ /Gm/ && $tName !~ /scaffold/i){
				$sw = 1;
			}
			elsif($allq eq '2' && $qSize >= $lnq){
				$sw = 1;
			}
		}
	}
	
	if($sw == 1){
		$hstat->{sum_matches}{$qName}{$tName} += $matches;
		$hstat->{sum_misMatches}{$qName}{$tName} += $misMatches;
		$hstat->{sum_repMatches}{$qName}{$tName} += $repMatches;
		$hstat->{sum_nCount}{$qName}{$tName} += $nCount;
		$hstat->{sum_blockCount}{$qName}{$tName} += $blockCount;
		$hstat->{sum_blockSizes}{$qName}{$tName} .= $blockSizes;
		
		$hstat->{sum_matches}{total} += $matches;
		$hstat->{sum_misMatches}{total} += $misMatches;
		$hstat->{sum_repMatches}{total} += $repMatches;
		$hstat->{sum_nCount}{total} += $nCount;
		$hstat->{sum_blockCount}{total} += $blockCount;
		$hstat->{sum_blockSizes}{total} .= $blockSizes;
		
		$hqsize->{$qName} = $qSize;
		$htsize->{$tName} = $tSize;
		$hmatches->{$qName}{$tName} += $matches;
		$hmisMatches->{$qName}{$tName} += $misMatches;
		$hblockCount->{$qName}{$tName} += $blockCount;
		$hq_tid->{$qName} .= $tName.",";
#		$hq_tpos->{$qName}{$tName}{pos0} = $tStart;
#		$hq_tpos->{$qName}{$tName}{pos1} = $tEnd;
		$hq_tpos->{$qName}{$tName}{POS} .= $tStarts.",";
		
		if($matches >= $lth){
			if($strand eq '-'){
				my $strlen = length($qStarts);
				$qStarts = substr($qStarts, 0, $strlen - 1);
				
				my @qpos = split(/,/, $qStarts);
				my @qposR;
				foreach my $qp (@qpos){
					$qp = $qSize - $qp;
					push(@qposR, $qp);
				}
				$qStarts = join(",", @qposR).",";
			}
			
			$hpos->{$qName}{$tName} .= $qStarts."\t".$tStarts."\t".$blockSizes."\n";
			$hstrand->{$qName}{$tName} .= $strand."\n";
			$hallen->{$qName}{$tName} += $matches;
			$cnt_selected++;
		}
		
		$cnt++;
	}
}
close $fh;

my @TID = keys(%{$htsize});
@TID = sort {$a cmp $b} @TID;

my @QID = keys(%{$hqsize});
if($allq eq '0'){
	@QID = sort {$a cmp $b} @QID;
}
else{
	my @QID2;
	foreach my $tid (@TID){
		my $AoQ = [];
		foreach my $qid (@QID){
			my @HID = split(/,/, $hq_tid->{$qid});
			my $AoH = [];
			foreach my $hid (@HID){
				if($hid){
					my @tmp;
					push(@tmp, $hid);
					push(@tmp, $hstat->{sum_matches}{$qid}{$hid});
					push(@{$AoH}, \@tmp);
				}
			}
			@{$AoH} = sort {$b->[1] <=> $a->[1]} @{$AoH};
			my $tophit_tid = $AoH->[0][0];
			
			if($tophit_tid eq $tid){
				my @TPOS = split(/,/, $hq_tpos->{$qid}{$tid}{POS});
				my $num_TPOS = @TPOS;
				my $sum = 0;
				foreach my $tpos (@TPOS){
					if($tpos){
						$sum += $tpos;
					}
				}
				my $tpos_avr = $sum / $num_TPOS;
				
				my @tmp;
				push(@tmp, $qid);
				push(@tmp, $hqsize->{$qid});
				push(@tmp, $tpos_avr);
				push(@{$AoQ}, \@tmp);
			}
		}
		@{$AoQ} = sort {$a->[2] <=> $b->[2]} @{$AoQ};
		
		foreach my $tmp (@{$AoQ}){
			if($tmp->[0] !~ /unanchored/i && $tmp->[0] !~ /chr00/i){
				push(@QID2, $tmp->[0]);
#				print "$tmp->[0] (1)\n";
			}
		}
	}
	
	foreach my $tid (@TID){
		my $AoQ = [];
		foreach my $qid (@QID){
			my @HID = split(/,/, $hq_tid->{$qid});
			my $AoH = [];
			foreach my $hid (@HID){
				if($hid){
					my @tmp;
					push(@tmp, $hid);
					push(@tmp, $hstat->{sum_matches}{$qid}{$hid});
					push(@{$AoH}, \@tmp);
				}
			}
			@{$AoH} = sort {$b->[1] <=> $a->[1]} @{$AoH};
			my $tophit_tid = $AoH->[0][0];
			
			if($tophit_tid eq $tid){
				my @TPOS = split(/,/, $hq_tpos->{$qid}{$tid}{POS});
				my $num_TPOS = @TPOS;
				my $sum = 0;
				foreach my $tpos (@TPOS){
					if($tpos){
						$sum += $tpos;
					}
				}
				my $tpos_avr = $sum / $num_TPOS;
				
				my @tmp;
				push(@tmp, $qid);
				push(@tmp, $hqsize->{$qid});
				push(@tmp, $tpos_avr);
				push(@{$AoQ}, \@tmp);
			}
		}
		@{$AoQ} = sort {$a->[2] <=> $b->[2]} @{$AoQ};
		
		foreach my $tmp (@{$AoQ}){
			if($tmp->[0] =~ /unanchored/i || $tmp->[0] =~ /chr00/i){
				push(@QID2, $tmp->[0]);
#				print "$tmp->[0] (2)\n";
			}
		}
	}
	
	@QID = @QID2;
}

my $numQ = @QID;
my $numT = @TID;

if($allq eq '0' && $QID[0] =~ /00/){
	my $chr00 = shift(@QID);
	push(@QID, $chr00);
}
if($TID[0] =~ /00/){
	my $chr00 = shift(@TID);
	push(@TID, $chr00);
}

my $stats_log;
$stats_log .= "! total [$cnt] lines with target ID in [$file]\n";
$stats_log .= "\n";
$stats_log .= "! [$hstat->{sum_matches}{total}] = total number of matching bases that aren't repeats\n";
$stats_log .= "! [$hstat->{sum_misMatches}{total}] = total number of bases that don't match\n";
$stats_log .= "! [$hstat->{sum_repMatches}{total}] = total number of matching bases that are part of repeats\n";
$stats_log .= "! [$hstat->{sum_nCount}{total}] = total number of 'N' bases\n";
$stats_log .= "! [$hstat->{sum_blockCount}{total}] = total number of blocks in the alignment\n";
$stats_log .= "\n";
$stats_log .= "! [$cnt_selected] lines with matches > [$lth] bp\n";
$stats_log .= "! [$numQ] = total number of seq ID in query\n";
$stats_log .= "! [$numT] = total number of seq ID in target\n";

if($numQ < 30 && $numT < 30){
	$stats_log .= "! seq IDs in query  = [".join(",", @QID)."]\n";
	$stats_log .= "! seq IDs in target = [".join(",", @TID)."]\n";
}
$stats_log .= "\n";

print $stats_log;

my $hqcumsize = {};
my $htcumsize = {};

my $stats_pos = "type\tseqID\tlength\tcumlative length\n";
my $qcumpos = 0;
foreach my $qid (@QID){
	$stats_pos .= "query\t".$qid."\t".$hqsize->{$qid}."\t".$qcumpos."\n";
	$hqcumsize->{$qid} = $qcumpos;
	$qcumpos += $hqsize->{$qid};
}
$stats_pos .= "query\t-\t-\t".$qcumpos."\n";

my $tcumpos = 0;
foreach my $tid (@TID){
	$stats_pos .= "target\t".$tid."\t".$htsize->{$tid}."\t".$tcumpos."\n";
	$htcumsize->{$tid} = $tcumpos;
	$tcumpos += $htsize->{$tid};
}
$stats_pos .= "target\t-\t-\t".$tcumpos."\n";

print "! converting value for position...\n";
my @qposFw;
my @qposRv;
my @tposFw;
my @tposRv;
my @bsizeFw;
my @bsizeRv;
my $hposdata = {};
my $stats_count = "query\ttarget\tnum block\tnum block (Fw)\tnum block (Rv)\tnum alignment\n";
my $stats_hmatches .= "query | target->\t".join("\t", @TID)."\n";
my $stats_hmisMatches .= "query | target->\t".join("\t", @TID)."\n";
my $stats_hblockCount .= "query | target->\t".join("\t", @TID)."\n";
for(my $q = 0; $q < $numQ; $q++){
	my $qid = $QID[$q];
	
	my @val0;
	my @val1;
	my @val2;
	for(my $t = 0; $t < $numT; $t++){
		my $tid = $TID[$t];
		
		if($hpos->{$qid}{$tid}){
			my @L = split(/\n/, $hpos->{$qid}{$tid});			# $qStarts."\t".$tStarts."\t".$blockSizes."\n";
			my @S = split(/\n/, $hstrand->{$qid}{$tid});
			my $numL = @L;
			my $numS = @S;
			
			unless($numL == $numS){
				print "! count unmatched between start_position and strand info for [$qid] [$tid]\n";
				die;
			}
			
			my $num_qpos = 0;
			my $num_tpos = 0;
			my $num_fwaln = 0;
			my $num_rvaln = 0;
			
			my @qposFw_sub;
			my @qposRv_sub;
			my @tposFw_sub;
			my @tposRv_sub;
			my @bsizeFw_sub;
			my @bsizeRv_sub;
			
			for(my $i = 0; $i < $numL; $i++){
				my @subL = split(/\t/, $L[$i]);
				my @qpos = split(/,/, $subL[0]);
				my @tpos = split(/,/, $subL[1]);
				my @bsize = split(/,/, $subL[2]);
				my $strand = $S[$i];
				
				$num_qpos += @qpos;
				$num_tpos += @tpos;
				
				if($strand eq '+'){
					$num_fwaln += 1;
				}
				elsif($strand eq '-'){
					$num_rvaln += 1;
				}
				
				unless($num_qpos == $num_tpos){
					print "! total count of position data unmatched between [$qid $num_qpos] and [$tid $num_tpos]\n";
					die;
				}
				
				foreach my $size (@bsize){
					if($strand eq '+'){
						push(@bsizeFw_sub, $size);
					}
					elsif($strand eq '-'){
						push(@bsizeRv_sub, $size);
					}
					
					if($strand eq '+'){
						push(@bsizeFw, $size);
					}
					elsif($strand eq '-'){
						push(@bsizeRv, $size);
					}
				}
				foreach my $p (@qpos){
					if($strand eq '+'){
						push(@qposFw_sub, $p);
					}
					elsif($strand eq '-'){
						push(@qposRv_sub, $p);
					}
					
					$p += $hqcumsize->{$qid};
					
					if($strand eq '+'){
						push(@qposFw, $p);
					}
					elsif($strand eq '-'){
						push(@qposRv, $p);
					}
				}
				foreach my $p (@tpos){
					if($strand eq '+'){
						push(@tposFw_sub, $p);
					}
					elsif($strand eq '-'){
						push(@tposRv_sub, $p);
					}
					
					$p += $htcumsize->{$tid};
					
					if($strand eq '+'){
						push(@tposFw, $p);
					}
					elsif($strand eq '-'){
						push(@tposRv, $p);
					}
				}
			}
			
			unless($hallen->{$qid}{$tid}){
				$hallen->{$qid}{$tid} = 0;
			}
			
			$hposdata->{$qid}{$tid}{qposFw} = \@qposFw_sub;
			$hposdata->{$qid}{$tid}{qposRv} = \@qposRv_sub;
			$hposdata->{$qid}{$tid}{tposFw} = \@tposFw_sub;
			$hposdata->{$qid}{$tid}{tposRv} = \@tposRv_sub;
			$hposdata->{$qid}{$tid}{bsizeFw} = \@bsizeFw_sub;
			$hposdata->{$qid}{$tid}{bsizeRv} = \@bsizeRv_sub;
			$hposdata->{$qid}{$tid}{qcumpos} = $hqsize->{$qid};
			$hposdata->{$qid}{$tid}{tcumpos} = $htsize->{$tid};
			$hposdata->{$qid}{$tid}{matches} = $hallen->{$qid}{$tid};
			$hposdata->{$qid}{$tid}{nblock} = $numL;
			$hposdata->{$qid}{$tid}{stats_pos} .= "query\t-\t-\t0\n";
			$hposdata->{$qid}{$tid}{stats_pos} .= "query\t-\t-\t".$hqsize->{$qid}."\n";
			$hposdata->{$qid}{$tid}{stats_pos} .= "target\t-\t-\t0\n";
			$hposdata->{$qid}{$tid}{stats_pos} .= "target\t-\t-\t".$htsize->{$tid};
			
#			print "! [$qid] [$tid] = [$numL] blocks (Fw=[$num_fwaln], Rv=[$num_rvaln]) with [$num_qpos] alignments\n";
			$stats_count .= $qid."\t".$tid."\t".$numL."\t".$num_fwaln."\t".$num_rvaln."\t".$num_qpos."\n";
			
			push(@val0, $hmatches->{$qid}{$tid});
			push(@val1, $hmisMatches->{$qid}{$tid});
			push(@val2, $hblockCount->{$qid}{$tid});
		}
		else{
			push(@val0, 0);
			push(@val1, 0);
			push(@val2, 0);
		}
	}
	
	$stats_hmatches .= $qid."\t".join("\t", @val0)."\n";
	$stats_hmisMatches .= $qid."\t".join("\t", @val1)."\n";
	$stats_hblockCount .= $qid."\t".join("\t", @val2)."\n";
}

$hposdata->{All}{qposFw} = \@qposFw;
$hposdata->{All}{qposRv} = \@qposRv;
$hposdata->{All}{tposFw} = \@tposFw;
$hposdata->{All}{tposRv} = \@tposRv;
$hposdata->{All}{bsizeFw} = \@bsizeFw;
$hposdata->{All}{bsizeRv} = \@bsizeRv;
$hposdata->{All}{qcumpos} = $qcumpos;
$hposdata->{All}{tcumpos} = $tcumpos;
$hposdata->{All}{stats_pos} = $stats_pos;

my $rh = {};
$rh->{QID} = \@QID;
$rh->{TID} = \@TID;
$rh->{hposdata} = $hposdata;
$rh->{stats_log} = $stats_log;
$rh->{stats_count} = $stats_count;
$rh->{stats_pos} = $stats_pos;
$rh->{stats_hmatches} = $stats_hmatches;
$rh->{stats_hmisMatches} = $stats_hmisMatches;
$rh->{stats_hblockCount} = $stats_hblockCount;

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
my $keyword = shift;
my $str = shift;

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
		unless($file =~ /\.fai/ || $file =~ /\.fa\.n../ || $file =~ /\.fa\.p../ || $file =~ /\.fasta\.n../ || $file =~ /\.fasta\.p../ || $file =~ /\.ffa\.n../ || $file =~ /\.csv/ || $file =~ /_BLK/){
			push(@files2, $file);
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


#----------------------------------------------------------
sub File_select_marker{
my $keyword = shift;
my $str = shift;

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
	if($file =~ /$keyword/ && $file =~ /aligned_/){
		unless($file =~ /\.fai/ || $file =~ /\.fa\.n../ || $file =~ /\.fa\.p../ || $file =~ /\.fasta\.n../ || $file =~ /\.fasta\.p../ || $file =~ /\.ffa\.n../ || $file =~ /_BLK/){
			push(@files2, $file);
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
my $tq1 = $q1."t";
if($tq1 eq 't'){
	print "[null]\n";
	return "null";
}

if($files2[$q1] && -e $files2[$q1]){
	print "[$files2[$q1]]\n";
	return $files2[$q1];
}
else{
	print "[null]\n";
	return "null";
}

}

