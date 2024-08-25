#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use Cwd;
use FindBin;
use Getopt::Long;
#use File::HomeDir;
use Sys::Hostname;

#-------------------------------------------------------------------------------
my $version = "\n";
$version .= "batch_asm2sv.pl version 1.42\n";
$version .= "last update: [2023\/1\/28]\n";
$version .= "copyright: ryoichi yano\n";

#print $version;

#-------------------------------------------------------------------------------

#---------------------------------------------------------//
my $help_command =<<"EOS";
Basic usage: 
  /path/to/batch_asm2sv.pl -d [reference fasta] -g [reference gff3] -l [gene list] -q [query fasta] -n [expand neighbor length] -o [Result_directory] --cpu [CPU thread]

--db or -d            : reference sequence (fasta, required)
--gff or -g           : reference gff (gff, required)
--list or -l          : list of gene to be analyzed (text, required)
--query or -q         : query fasta (fasta, required)
--neighbor or -n      : include [X] bp neighboring seq in addition to target (e.g. promoter)
--wo_orfsearch or -w  : without ORF prediction
--cpu or -c           : CPU thread number (default 8)
--output_dir or -o    : output directory
--plot                : output alignment plot then stop
--help or -h          : display usage

EOS

#---------------------------------------------------------//
my $host = hostname();
my $script_path = $FindBin::Bin;
my $wpath = getcwd();
#my $HOME = File::HomeDir->my_home;
my $bin_tigplot3 = $script_path."/Psl_tigPlot3.pl";
my $bin_findorf = $script_path."/findorf.pl";
my $bin_psl2dot = $script_path."/psl2dot.pl";
my $bin_BBDH = $script_path."/Blast_BDH.pl";
my $bin_NBDH = $script_path."/Blast_NBDH.pl";
my $bin_gff2prot = $script_path."/gfftoprot.pl";
my $ramsize = CheckRAM();		#Mb scale

#---------------------------------------------------------//
#gettin parameters from command line
my $dir;
my $dbfasta;
my $dbgff;
my $genelist;
my $qfasta;
my $qpsl;
my $outplot;
my $neighbor_tlen;
my $cpu;
my $cpu_gth;
my $initspan;
my $redo_examine;
my $redo_orfsearch;
my $skip_orfsearch;
my $cluster;
my $help;
my $test;

GetOptions('--db=s' => \$dbfasta, '--gff=s' => \$dbgff, '--list=s' => \$genelist, '--query=s' => \$qfasta, '--alignpsl=s' => \$qpsl, '--output_dir=s' => \$dir, '--neighbor=i' => \$neighbor_tlen, '--cpu=i' => \$cpu, '--xthread=i' => \$cpu_gth, '--span=i' => \$initspan, '--plot' => \$outplot, '--test' => \$test, '--1' => \$cluster, '--5' => \$redo_examine, '--7' => \$redo_orfsearch, '--wo_orfsearch' => \$skip_orfsearch, '--help' => \$help);

if($cluster){
	$host = "cluster1";
}

my $status = "OK";
if(! $dbfasta || ! -e $dbfasta){
	print "! missing db fasta\n";
	$status = "missing";
}
if(! $dbgff || ! -e $dbgff){
	print "! missing db gff\n";
	$status = "missing";
}
if(! $genelist || ! -e $genelist){
	print "! missing target gene ID list\n";
	$status = "missing";
}
if(! $qfasta){
	print "! missing query fasta\n";
	$status = "missing";
}
elsif(! -e $qfasta){
	print "! missing query fasta [$qfasta]\n";
	$status = "missing";
}
if($status eq 'missing'){
	print "\n$help_command";
	goto END;
}
if($redo_examine){
	$skip_orfsearch = 0;
	$redo_orfsearch = 0;
}

if(! $qpsl || ! -e $qpsl){
	$qpsl = "null";
}
$qpsl = "null";

if(! $initspan || $initspan =~ /\D/){
	$initspan = 20000;
}

unless($neighbor_tlen){
	$neighbor_tlen = 1000;
}
elsif($neighbor_tlen =~ /\D/){
	$neighbor_tlen = 1000;
}
elsif($neighbor_tlen < 0){
	$neighbor_tlen = 0;
}
elsif($neighbor_tlen > 10000){
	$neighbor_tlen = 10000;
}

unless($cpu){
	$cpu = 8;
}
elsif($cpu =~ /\D/){
	$cpu = 8;
}
elsif($cpu > 192){
	$cpu = 192;
}
unless($cpu_gth){
	$cpu_gth = $cpu;
}
elsif($cpu_gth > $cpu){
	$cpu_gth = $cpu;
}
elsif($cpu_gth < 0){
	$cpu_gth = 1;
}
unless($test){
	$test = 0;
}
my $cpu_blast = $cpu;

my $log = "\n";
$log .= "! reference fasta     : [$dbfasta]\n";
$log .= "! reference gff       : [$dbgff]\n";
$log .= "! target gene list    : [$genelist]\n";
$log .= "! query fasta         : [$qfasta]\n";
#$log .= "! query alignment psl : [$qpsl]\n";
$log .= "! neighbor seq target : [$neighbor_tlen] bp (e.g. promoter region)\n";

if($skip_orfsearch){
	$log .= "! --wo_orfsearch is invoked, skip ORF prediction...\n";
}

#print $log,"\n";

my @APP;
push(@APP, "samtools");
push(@APP, "gffread");
push(@APP, "gth");
push(@APP, "blat");
push(@APP, "blat2hints");
push(@APP, "matcher");
push(@APP, "miniprot");

my $hbin = {};
foreach my $app (@APP){
	$hbin->{$app}[0] = $app;
}
my $hpath = {};
my $hbinpath = Search_app2($hbin, \@APP, $hpath, $script_path, "pipe", $genelist, $qfasta);
if($hbinpath->{err}){
	print "! abort script due to missing program(s)..., please confirm program is in PATH.\n";
	goto END;
}

#---------------------------------------------------------------//
my $dbpref = $dbfasta;
if($dbpref =~ /\.fasta/){
	$dbpref =~ s/\.fasta//g;
}
elsif($dbpref =~ /\.fa/){
	$dbpref =~ s/\.fa//g;
}

my $qpref = $qfasta;
if($qpref =~ /\.fasta/){
	$qpref =~ s/\.fasta//g;
}
elsif($qpref =~ /\.fa/){
	$qpref =~ s/\.fa//g;
}
if($qpref =~ /\.contigs_split4/){
	my @QPF = split(/\.contigs_split/, $qpref);
	$qpref = $QPF[0];
}

unless($dir){
	$dir = "result_GenePAV_".$dbpref."_vs_".$qpref;
}
unless(-e $dir){
	system("mkdir $dir");
}
chdir $dir;

unless(-e $qpref){
	system("mkdir $qpref");
}

unless(-e $dbfasta){
	system("ln -s ../$dbfasta ./");
}
unless(-e $dbgff){
	system("ln -s ../$dbgff ./");
}
unless(-e $genelist){
	system("ln -s ../$genelist ./");
}
unless(-e $qfasta){
	system("cp ../$qfasta ./");
}
#if($qpsl ne 'null' && ! -e $qpsl){
#	system("cp ../$qpsl ./");
#}

unless(-e "cmds"){
	system("mkdir cmds");
}
unless(-e "log"){
	system("mkdir log");
}

my $combined_final = "summary_genome2sv_".$dir."_combined.tsv";
my $rfile_final = "summary_genome2sv_".$dir."_".$qpref.".tsv";
my $rfile_final_BDH = "summary_genome2sv_BDH_".$dir."_".$qpref.".tsv";
my $afile = "./$qpref/specified_region_hits_".$qpref.".tsv";
my $catafile = "specified_region_hits.tsv";
my $redo_log = "_redo_examine.log";

if($redo_examine && -e $redo_log){
	print "! --5 is invoked, but analysis appears to be already done. skip...\n";
	goto END;
}

my $hqseq = Open_fasta_as_hash($qfasta);
my $hdbseq = Open_fasta_as_hash($dbfasta);
my $HGFF = Gff_to_hash($dbgff);
my $hgff = $HGFF->{hgff};
my $dbt2g = $HGFF->{t2g};
my $dbg2t = $HGFF->{g2t};
my $dbt2CDS = $HGFF->{CDS};

print "! analyzing target gene list...\n";
my $log_genelist = "log_genelist.tsv";
my $hlist = Get_genelist($dbpref, $genelist, $log_genelist, $hgff, $dbg2t, $dbt2CDS);

my @G = keys(%{$hlist});
#@G = sort {$a cmp $b} @G;
my $num = @G;
my $d = int($num / $cpu) + 1;

#---------------------------------------------------------------//
my $disable_skipaln = 1;
my $add_examine = 0;

if(-e $afile && ! $outplot && $disable_skipaln eq '0'){
	print "! [$afile] already exists, skip step.1 ...\n";
}
elsif($redo_orfsearch && -e $afile){
	print "! --7 is invoked, skip step.1 ...\n";
	$redo_examine = 1;
}
else{
	print "! step.1 analyzing presence-absence in [$qfasta]...\n";
	#-----------------------------------------------------//
	my $wgpsldir = "null";
#	if(-e $bin_psl2dot && -e $qpsl){
#		$wgpsldir = "_psl2dotalign";
#		if(-d $wgpsldir){
#			system("rm -R $wgpsldir");
#		}
#		
#		print "! convert whole genome alignment [$qpsl] to hints...\n";
#		my $cmd1 = "$bin_psl2dot $qpsl $wgpsldir > /dev/null 2>&1";
#		print "! cmd=[$cmd1]\n";
#		system($cmd1);
#		SortPsl($qpsl, "$wgpsldir/hint.psl");
#	}
	
	my $hpreva = {};
	my @Gnotyet;
	my $numQL_Gnotyet = @G;
	if(-e $afile){
		print "! found previous result [$afile]\n";
		$hpreva = Check_prev_results($afile);
#		$hpreva = Update_prev_results($afile);		# delete entries without exact 'present' tag (including promoter/UTR)
		
		foreach my $gid (@G){
			if(! $hpreva->{$gid}){
				push(@Gnotyet, $gid);
			}
		}
		$numQL_Gnotyet = @Gnotyet;
		print "! [$numQL_Gnotyet] not yet analyzed...\n";
		
		if($numQL_Gnotyet > 0){
			$add_examine = 1;
		}
	}
	else{
		@Gnotyet = @G;
	}
	my $dorf = int($numQL_Gnotyet / $cpu) + 1;
	
	#-----------------------------------------------------//
	my @gabbage;
	my @I1;
	my @I3;
	my @I4;
	for(my $i = 0; $i < $cpu; $i++){
		my $n0 = $i * $dorf;
		my $n1 = ($i + 1) * $dorf;
		
		my $cmd;
		my $cnt = 0;
		my $cnt_skip = 0;
		my $log_pav = "thread_".$i."_".$dbpref."_".$qpref.".log";
		for(my $j = $n0; $j < $n1; $j++){
			if($Gnotyet[$j]){
				my $gid = $Gnotyet[$j];
				my $sid = $hlist->{$gid}{Chr};
				my $len_sid = $hdbseq->{$sid}{len};
				
				my $p0 = $hlist->{$gid}{pos0} - $initspan;
				my $p1 = $hlist->{$gid}{pos1} + $initspan;
				my $t0 = $hlist->{$gid}{pos0};
				my $t1 = $hlist->{$gid}{pos1};
				
				if($p0 < 1){
					$p0 = 1;
				}
				if($p1 > $len_sid){
					$p1 = $len_sid;
				}
				
				my $region = $p0."-".$p1;
				my $specified = $t0."-".$t1;
				
#				my $log_pav = $dbpref."_".$sid."_".$region."_".$qpref.".txt";
				my $dbfasta_tmp = $i."_tmp_".$dbfasta;
				unless(-e $dbfasta_tmp){
					system("ln -s $dbfasta $dbfasta_tmp");
					push(@gabbage, $dbfasta_tmp);
				}
				
				my $qfasta_tmp = $i."_tmp_".$qfasta;
				unless(-e $qfasta_tmp){
					system("ln -s $qfasta $qfasta_tmp");
					push(@gabbage, $qfasta_tmp);
				}
				
				if(! $hpreva->{$gid}){
#					print "$bin_tigplot3 -o $qpref -d $dbfasta_tmp -q $qfasta_tmp --scale bp -I $sid -r $region -t $specified -m 95 --note $gid >> $log_pav\n";
					if($outplot){
						$cmd .= "$bin_tigplot3 -o $qpref -d $dbfasta_tmp -q $qfasta_tmp --scale bp -I $sid -r $region -t $specified -m 95 --note $gid --neighbor $neighbor_tlen --plot --neighbor $neighbor_tlen --strand $hgff->{$gid}{strand} >> $log_pav 2>&1\n";
					}
					else{
						if($wgpsldir ne 'null'){
							$cmd .= "$bin_tigplot3 -o $qpref -d $dbfasta_tmp -q $qfasta_tmp -w $wgpsldir --scale bp -I $sid -r $region -t $specified -m 95 --note $gid --neighbor $neighbor_tlen --strand $hgff->{$gid}{strand}";
						}
						else{
							$cmd .= "$bin_tigplot3 -o $qpref -d $dbfasta_tmp -q $qfasta_tmp --scale bp -I $sid -r $region -t $specified -m 95 --note $gid --neighbor $neighbor_tlen --strand $hgff->{$gid}{strand}";
						}
						
#						if($initspan >= 50000){
#							$cmd .= " --skip_alignmethod";
#						}
						$cmd .= " >> $log_pav 2>&1\n";
#						$cmd .= "echo \'$gid\' >> hist.txt\n";
						$cmd .= "sleep 1\n";
					}
				}
				else{
#					$cmd .= "echo \'already done, skipping $identifer\'\n";
					$cnt_skip++;
				}
				$cnt++;
			}
		}
		
		if($cmd){
			$cmd .= "echo 'all_done' >> $log_pav\n";
			my $script;
			if($outplot){
				$script = "cmd_aln2sv_".$i.".sh";
				ADD($script, $cmd);
			}
			else{
				$script = "cmd_aln2sv_".$qpref."_".$i.".sh";
				SAVE($script, $cmd);
			}
			if(-e $log_pav){
				system("rm $log_pav");
			}
			push(@I1, $script);
			push(@I3, $cnt);
			push(@I4, $log_pav);
		}
		elsif($cnt_skip > 0){
			print "! thread [$i] : [$cnt_skip] queries skipped\n";
		}
	}
	
	if($host eq 'cluster1'){
		my $num_I1 = @I1;
		for(my $i = 0; $i < $num_I1; $i++){
			system("qsub $I1[$i]");
		}
	}
	else{
		my $num_I1 = @I1;
		my $thrs = [];
		for(my $i = 0; $i < $num_I1; $i++){
			my $thr = threads->new(\&Batch_exe, $I1[$i], $i, $I3[$i]);
			push(@{$thrs}, $thr);
			sleep(1);
		}
		foreach my $thr (@{$thrs}){
			$thr->join();
		}
	}
	
	if($host eq 'cluster1'){
		print "! waiting for job completion...\n";
		LOOPjob:{
			my $loopjob = 1;
		}
		my $num_I4 = @I4;
		my $ncompleted_pav = 0;
		foreach my $log_pav (@I4){
			$ncompleted_pav += JudgeJob($log_pav);
		}
		if($ncompleted_pav < $num_I4){
			sleep(30);
			goto LOOPjob;
		}
		else{
			for(my $i = 0; $i < $num_I4; $i++){
				if(-e $I1[$i]){
					my $mvcmd1 = "mv $I1[$i]"."* cmds";
					my $mvcmd2 = "mv $I4[$i]"."* log";
					system($mvcmd1);
					system($mvcmd2);
				}
			}
		}
	}
	
	print "! removing temp files...\n";
	Rm_files(\@I1);
	Mv_files(\@I4, "cmds");
	Rm_files(\@gabbage);
	
	if(-d $wgpsldir){
		system("rm -R $wgpsldir");
	}
}

if(-e $afile){
	my $pcnt_th = 90;
	my $pcov_th = 5;
	my $eval_th = "1e-60";
	my $rfile_nBDH1 = "./blast_nBDH/BDconfidenthit_pcnt-".$pcnt_th."_pcov-".$pcov_th."_Eval-".$eval_th.".csv";
	
	if(-e $rfile_nBDH1 && $add_examine){
		print " - once remove [blast_nBDH]\n";
		system("rm -R blast_nBDH");
	}
	if(! -e $rfile_nBDH1){
		my $dbgene_region = $dbpref."_gene_region.fasta";
		Convert2eachseq_Gff($dbgff, $hdbseq, $dbgene_region);
		
		my $qgene_region = $qpref."_gene_region.fasta";
		Convert2eachseq($afile, $hqseq, $qgene_region, $qpref);
		
		my $cpu_nblast = $cpu_blast;
		if($cpu_nblast > 16){
			$cpu_nblast = 16;
		}
		
		if(-e $qgene_region){
			print "! blast-n for BDH search with [$cpu_nblast] threads...\n";
			my $cmd_nBDH = "$bin_NBDH blast_nBDH $dbgff $dbgene_region $qgene_region $qpref 5 $pcnt_th $pcov_th $eval_th $cpu_nblast y > ./log/log_nBDH.txt 2>&1";
#			print "! cmd=[$cmd_nBDH]\n";
			if(system($cmd_nBDH) != 0){
				print "! blast_NBDH failed, abort script...\n";
				die;
			}
		}
		else{
			print "! missing [$qgene_region]...\n";
		}
	}
	else{
		print "! [$rfile_nBDH1] already exists, skip...\n";
	}
	
	if(-e $rfile_nBDH1){
		Revise_afile_nBDH($afile, $rfile_nBDH1);
	}
}

#if($outplot){
#	Combine_text($afile, $catafile);
#	
#	if(-e $qfasta){
#		system("rm $qfasta");
#	}
#	if(-d $qpref){
#		system("cp $qpref/*.png ./");
#		system("rm -R $qpref");
#	}
#	print "\n";
#	goto END;
#}

#---------------------------------------------------------------//
my $header_summary = "directory\tdbfasta\tseqid\tpos0\tpos1\tsp0\tsp1\tqfasta\tqfasta seqid\thit pos0\thit pos1\tgroup\thit sp0\thit sp1\tbp (sp0-sp1)\tbp hit-align (sp0-sp1)\tbp hit-span (sp0-sp1)\tbp hit-span woN (sp0-sp1)\tbp insertion\talign ratio (1=not disrupted)\tinsert ratio (1=not disrupted)\tseq normality (1=not disrupted)\tjudge\tgene_id\t";
$header_summary .= "5'-promoter bp hit ($neighbor_tlen bp)\t5'-promoter align ratio ($neighbor_tlen bp)\t3'-UTR bp hit ($neighbor_tlen bp)\t3'-UTR align ratio ($neighbor_tlen bp)\t";
$header_summary .= "5' flanking bp hit ($initspan bp)\t5' flanking align ratio ($initspan bp)\t3' flanking bp hit ($initspan bp)\t3' flanking align ratio ($initspan bp)\t";
$header_summary .= "protein coding\tdb protein length\tq protein length (predicted)\talign length\t\%similarity\t\%coverage\tscore\n";

if($skip_orfsearch || $outplot){
	if($skip_orfsearch){
		print "! --wo_orfsearch is invoked, skip ORF prediction...\n";
	}
	elsif($outplot){
		print "! --plot is invoked, skip ORF prediction...\n";
	}
	print "! preparing final result...\n";
	my $dbprotein = $dbpref."_protein.fasta";
	my $cmd3 = "$hbinpath->{gffread} $dbgff -g $dbfasta -y $dbprotein > /dev/null 2>&1";
	Cmdexe_unless_seq($cmd3, $dbprotein);
	my $hpav = Open_pav_results($afile);
	my $hdbprot = Open_fasta_as_hash($dbprotein);
	
	my $tmp_ralign = "info_ralign.tsv";
	if(-e $tmp_ralign){
		system("rm $tmp_ralign");
	}
	
	my $hprevr = {};
	if(-e $rfile_final){
		$hprevr = Check_prev_results($rfile_final);
	}
	
	my @G = keys(%{$hlist});
	@G = sort {$a cmp $b} @G;
	
	my $summary;
	foreach my $gid (@G){
		if($hprevr->{$gid}){
			$summary .= $hprevr->{$gid}."\n";
		}
		elsif($hpav->{$gid}{data}){
			my $len_dbprot = "";
			if($dbg2t->{$gid}){
				my @dbTID = split(/\n/, $dbg2t->{$gid});
				foreach my $dbtid (@dbTID){
					if($hdbprot->{$dbtid}{len}){
						$len_dbprot .= $hdbprot->{$dbtid}{len}.",";
					}
				}
			}
			unless($len_dbprot){
				$len_dbprot = "-";
			}
			
			my $len_qprot = "-";
			$summary .= $hpav->{$gid}{data}."\t".$hlist->{$gid}{protein_coding}."\t".$len_dbprot."\t".$len_qprot."\t-\t-\t-\t-\n";
		}
	}
	
	if($summary){
		if($outplot){
			if(-e $combined_final){
				ADD($combined_final, $summary);
				print "! add result to [$combined_final]\n";
			}
			else{
				$summary = $header_summary.$summary;
				SAVE($combined_final, $summary);
				print "! output [$combined_final]\n";
			}
		}
		else{
			$summary = $header_summary.$summary;
			SAVE($rfile_final, $summary);
			print "! output [$rfile_final]\n";
		}
	}
}
elsif(-e $afile){
	print "! step.2 verify gene presence based on protein alignment...\n";
	if($host ne 'cluster1'){
		my $max_cpu = int($ramsize / 5.30 / 1000);
		if(! $max_cpu){
			$max_cpu = 48;
		}
		if($cpu_gth ne $cpu){
			$max_cpu = $cpu_gth;
		}
		if($max_cpu > 48){
			$max_cpu = 48;
		}
		if($cpu > $max_cpu){
			print "! CPU threads for ORF search [$cpu] > [$max_cpu]\n";
			$cpu = $max_cpu;
		}
	}
	else{
		if($cpu_blast > 16){
			$cpu_blast = 16;
		}
	}
	
	my $hpav = Open_pav_results($afile);
	
	my $dbfai = $dbfasta.".fai";
	my $dbtranscript = $dbpref."_transcript.fasta";
	my $dbcds = $dbpref."_cds.fasta";
	my $dbprotein = $dbpref."_protein.fasta";
#	my $cmd0 = "$hbinpath->{samtools} faidx $dbfasta";
	my $cmd1 = "$hbinpath->{gffread} $dbgff -g $dbfasta -w $dbtranscript > /dev/null 2>&1";
	my $cmd2 = "$hbinpath->{gffread} $dbgff -g $dbfasta -x $dbcds > /dev/null 2>&1";
	my $cmd3 = "$hbinpath->{gffread} $dbgff -g $dbfasta -y $dbprotein > /dev/null 2>&1";
#	Cmdexe_unless_seq($cmd0, $dbfai);
	Cmdexe_unless_seq($cmd1, $dbtranscript);
	Cmdexe_unless_seq($cmd2, $dbcds);
	Cmdexe_unless_seq($cmd3, $dbprotein);
	my $hdbprot = Open_fasta_as_hash($dbprotein);
	my $hdbtranscript = Open_fasta_as_hash($dbtranscript);
	my $hdbcds = Open_fasta_as_hash($dbcds);
	
	print "! preparing dataset for ORF search... (may take a while)\n";
	my $hprevr = {};
	my @Gnotyet;
	my $numQL_Gnotyet = @G;
	my $num_prevorf = 0;
	my $combined_qprot = "predicted_".$qpref.".fasta";
	my $combined_qgff = "predicted_".$qpref.".gff";
	my $revised_qprot = "predicted_revised_".$qpref.".fasta";
	my $revised_qgff = "predicted_revised_".$qpref.".gff";
	
	if($redo_examine){
		print "! --5 is invoked\n";
		if(-e $rfile_final){
			print " - once remove [$rfile_final]\n";
			system("rm $rfile_final");
		}
		if(-e $rfile_final_BDH){
			print " - once remove [$rfile_final_BDH]\n";
			system("rm $rfile_final_BDH");
		}
		if(-e $combined_qprot){
			print " - once remove [$combined_qprot]\n";
			system("rm $combined_qprot");
		}
		if(-e $combined_qgff){
			print " - once remove [$combined_qgff]\n";
			system("rm $combined_qgff");
		}
		@Gnotyet = @G;
	}
	
	if(-e $rfile_final){
		print "! found previous result [$rfile_final]\n";
		my $hcheckorf = Check_prev_orfsearch($rfile_final);
		$hprevr = $hcheckorf->{hash};
		$num_prevorf = $hcheckorf->{norf};
		
		foreach my $gid (@G){
			if(! $hprevr->{$gid}){
				push(@Gnotyet, $gid);
			}
		}
		$numQL_Gnotyet = @Gnotyet;
		print " [$numQL_Gnotyet] not yet analyzed...\n";
		print " [$num_prevorf] with ORF information\n";
	}
	else{
		@Gnotyet = @G;
	}
	my $dorf = int($numQL_Gnotyet / $cpu) + 1;
	
	my $tmp_ralign = "info_ralign.tsv";
	if(-e $tmp_ralign){
		system("rm $tmp_ralign");
	}
	
	my $log_orfsearch_hist = "log_orfsearch_processed.txt";
	if(-e $log_orfsearch_hist){
		system("rm $log_orfsearch_hist");
	}
	
	my $tmpdir_sw = 0;
#	if(-e "/media/nvmer/tmp"){
#		$tmpdir = "/media/nvmer/tmp"."/_temp_asm2sv_".$qpref;
#		unless(-e $tmpdir){
#			system("mkdir $tmpdir");
#			$tmpdir_sw = 1;
#		}
#	}
	
	@G = keys(%{$hlist});
	my $ralign;
	my @empty;
	my @gabbage;
	my @I1;
	my @I4;
#	if($numQL_Gnotyet > 0 && $num_prevorf == 0){
	if($numQL_Gnotyet > 0 || $add_examine == 1){
		for(my $i = 0; $i < $cpu; $i++){
			my $n0 = $i * $dorf;
			my $n1 = ($i + 1) * $dorf;
			
			my $tmpdir = "./_temp_$i";
			unless(-e $tmpdir){
				system("mkdir $tmpdir");
			}
			
			my $cmd;
			my $cnt = 0;
			my $cnt_skip = 0;
			my $log_orfsearch = "thread_orfsearch_".$i.".log";
			my $script = "cmd_orfsearch_".$qpref."_".$i.".sh";
			for(my $j = $n0; $j < $n1; $j++){
				if($Gnotyet[$j]){
					my $gid = $Gnotyet[$j];
					
					if($hpav->{$gid}{status} && $hpav->{$gid}{status} eq 'P'){
						my $dbsid = $hpav->{$gid}{dbsid};
						my $dbpos0 = $hpav->{$gid}{dbpos0};
						my $dbpos1 = $hpav->{$gid}{dbpos1};
						my $qsid = $hpav->{$gid}{qsid};
						my $qpos0 = $hpav->{$gid}{qpos0};
						my $qpos1 = $hpav->{$gid}{qpos1};
						
						my $tmp_dbprot = "$tmpdir/_dbprot_".$gid.".fasta";
						my $tmp_dbtranscript = "$tmpdir/_dbtranscript_".$gid.".fasta";
						my $tmp_dbcds = "$tmpdir/_dbcds_".$gid.".fasta";
						my $tmp_dbindex = "$tmpdir/_dbprot_".$gid.".fasta.protein";
						my $tmp_dbfasta = "$tmpdir/_dbseq_".$gid.".fasta";
						my $tmp_qfasta = "$tmpdir/_qseq_".$qpref."_".$gid.".fasta";
						my $tmp_qgff = "$tmpdir/_qprot_".$qpref."_".$gid.".gff";
						my $tmp_qprot = "$tmpdir/_qprot_".$qpref."_".$gid.".fasta";
						my $tmp_align = "$tmpdir/_align_".$gid.".txt";
						my $mk_process = "$tmpdir/_mk_".$gid.".txt";
						
						if($gid =~ /\W/){
							my $modgid = $gid;
							$modgid =~ s/\W/_/g;
							$tmp_dbprot = "$tmpdir/_dbprot_".$modgid.".fasta";
							$tmp_dbtranscript = "$tmpdir/_dbtranscript_".$modgid.".fasta";
							$tmp_dbcds = "$tmpdir/_dbcds_".$modgid.".fasta";
							$tmp_dbindex = "$tmpdir/_dbprot_".$modgid.".fasta.protein";
							$tmp_dbfasta = "$tmpdir/_dbseq_".$modgid.".fasta";
							$tmp_qfasta = "$tmpdir/_qseq_".$qpref."_".$modgid.".fasta";
							$tmp_qgff = "$tmpdir/_qprot_".$qpref."_".$modgid.".gff";
							$tmp_qprot = "$tmpdir/_qprot_".$qpref."_".$modgid.".fasta";
							$tmp_align = "$tmpdir/_align_".$modgid.".txt";
							$mk_process = "$tmpdir/_mk_".$modgid.".txt";
						}
						
						if($dbg2t->{$gid}){
							if($hprevr && $hprevr->{$gid}){
								$cnt_skip++;
							}
							else{
								my @dbTID = split(/\n/, $dbg2t->{$gid});
								unless(-e $tmp_dbprot){
									my $each_dbprot;
									foreach my $dbtid (@dbTID){
										if($hdbprot->{$dbtid}{seq}){
											$each_dbprot .= ">".$dbtid."\n".$hdbprot->{$dbtid}{seq}."\n";
										}
									}
									SAVE($tmp_dbprot, $each_dbprot);
								}
								unless(-e $tmp_dbtranscript){
									my $each_dbtranscript;
									foreach my $dbtid (@dbTID){
										if($hdbtranscript->{$dbtid}{seq}){
											$each_dbtranscript .= ">".$dbtid."\n".$hdbtranscript->{$dbtid}{seq}."\n";
										}
									}
									SAVE($tmp_dbtranscript, $each_dbtranscript);
								}
								unless(-e $tmp_dbcds){
									my $each_dbcds;
									foreach my $dbtid (@dbTID){
										if($hdbcds->{$dbtid}{seq}){
											$each_dbcds .= ">".$dbtid."\n".$hdbcds->{$dbtid}{seq}."\n";
										}
									}
									SAVE($tmp_dbcds, $each_dbcds);
								}
								
								my $qsid_seqlen = length($hqseq->{$qsid}{seq});
	#							my $qpos0_nb = $qpos0;
	#							my $qpos1_nb = $qpos1;
								my $qpos0_nb = $qpos0 - $neighbor_tlen;
								my $qpos1_nb = $qpos1 + $neighbor_tlen;
								if($qpos0_nb < 0){
									$qpos0_nb = 1;
								}
								if($qpos1_nb > $qsid_seqlen){
									$qpos1_nb = $qsid_seqlen;
								}
								
								my $sub_dbseq = Extract_and_save($tmp_dbfasta, $dbsid, $hdbseq->{$dbsid}{seq}, $dbpos0, $dbpos1);
								my $sub_qseq = Extract_and_save($tmp_qfasta, $qsid, $hqseq->{$qsid}{seq}, $qpos0_nb, $qpos1_nb);
								my $len_sub_dbseq = length($sub_dbseq);
								my $len_sub_qseq = length($sub_qseq);
								
	#							$cmd .= "if test -e \"$mk_process\"; then\n";
	#							$cmd .= "\techo '$gid already processed' >> $log_orfsearch 2>&1\n";
	#							$cmd .= "fi\n";
								
	#							if($sub_qseq && -e $tmp_dbprot && -e $tmp_qfasta){
								if(-e $tmp_dbprot && -e $tmp_qfasta){
									$cmd .= "perl $bin_findorf $hbinpath->{gth} $hbinpath->{miniprot} $hbinpath->{samtools} $hbinpath->{gffread} $hbinpath->{blat} $hbinpath->{blat2hints} $hbinpath->{matcher} $gid $tmp_dbprot $tmp_dbtranscript $tmp_dbcds $tmp_qfasta $tmp_qgff $tmp_qprot $tmp_align $tmpdir $qpref $qpos0_nb $combined_qprot $combined_qgff >> $log_orfsearch 2>&1\n";
								}
								
								# time-consuming...
	#							my $len_NNN = Count_NNN($sub_qseq);
	#							$ralign .= $gid."\t".$tmp_align."\t".$len_NNN."\t".$tmp_qprot."\n";
								$ralign .= $gid."\t".$tmp_align."\t-\t".$tmp_qprot."\t-\t".$len_sub_dbseq."\t".$len_sub_qseq."\n";
								
	#							push(@gabbage, $tmp_dbprot);
	#							push(@gabbage, $tmp_qfasta);
	#							push(@gabbage, $tmp_qgff);
	#							push(@gabbage, $tmp_qprot);
	#							push(@gabbage, $tmp_align);
							}
						}
					}
					$cnt++;
				}
			}
			
			if($cnt_skip > 0){
				print "! [$cnt_skip] queries skipped\n";
			}
			
			if($cmd){
				if(-e $log_orfsearch){
					system("rm $log_orfsearch");
				}
				$cmd .= "echo 'all_done' >> $log_orfsearch\n";
				SAVE($script, $cmd);
				push(@I1, $script);
				push(@I4, $log_orfsearch);
			}
		}
		
		if($ralign){
			SAVE($tmp_ralign, $ralign);
		}
		else{
			print "! no candidate found...\n";
			if(! $outplot){
				system("touch $rfile_final");
			}
			else{
				if(-e $qfasta){
					system("rm $qfasta");
				}
#				if(-d $qpref){
#					system("rm -R $qpref");
#				}
			}
			
			print "! removing temp files...\n";
			for(my $i = 0; $i < $cpu; $i++){
				my $tmpdir = "./_temp_$i";
				if(-d $tmpdir){
					system("rm -R $tmpdir");
				}
			}
			
			Rm_files(\@I1);
			Mv_files(\@I4, "cmds");
			Rm_files(\@gabbage);
			
			goto END;
		}
		
		if($host eq 'cluster1'){
			print "! ORF search with [$cpu] nodes...\n";
			my $num_I1 = @I1;
			for(my $i = 0; $i < $num_I1; $i++){
				system("qsub $I1[$i]");
			}
		}
		else{
			print "! ORF search with [$cpu] threads...\n";
			my $num_I1 = @I1;
			my $thrs2 = [];
			for(my $i = 0; $i < $num_I1; $i++){
				my $thr = threads->new(\&Batch_exe2, $I1[$i], $i, $I4[$i]);
				push(@{$thrs2}, $thr);
			}
			foreach my $thr (@{$thrs2}){
				$thr->join();
			}
		}
	}
	
	if($host eq 'cluster1'){
		print "! waiting for job completion...\n";
		LOOPjob:{
			my $loopjob = 1;
		}
		my $num_I4 = @I4;
		my $ncompleted_gth = 0;
		foreach my $log_orfsearch (@I4){
			$ncompleted_gth += JudgeJob($log_orfsearch);
		}
		if($ncompleted_gth < $num_I4){
			sleep(30);
			goto LOOPjob;
		}
		else{
			for(my $i = 0; $i < $num_I4; $i++){
				if(-e $I1[$i]){
					my $mvcmd1 = "mv $I1[$i]"."* cmds";
					my $mvcmd2 = "mv $I4[$i]"."* log";
					system($mvcmd1);
					system($mvcmd2);
				}
			}
		}
	}
	
#	print "! preparing protein fasta data for BDH search...\n";
#	my $cmd13 = "$hbinpath->{gffread} $combined_qgff -g $qfasta -y $combined_qprot > /dev/null 2>&1";
	my $cmd13 = "$bin_gff2prot $qfasta $combined_qgff $combined_qprot";
	system($cmd13);
	my $num_converted_protseq = Count_seq_fasta($combined_qprot);
	
	my $pcnt_th = 50;
	my $pcov_th = 0.001;
	my $eval_th = "1e-30";
	my $rfile_BDH1 = "./blast_BDH/BDconfidenthit_pcnt-".$pcnt_th."_pcov-".$pcov_th."_Eval-".$eval_th.".csv";
	
	if($num_converted_protseq && $num_converted_protseq > 0){
		print "! blast-p for BDH search with [$cpu_blast] threads...\n";
		Purse_seq($dbprotein);			#remove "." attacched to the end of protein sequences
	#	Purse_seq($combined_qprot);		#remove "." attacched to the end of protein sequences
		
		if($redo_examine || $add_examine){
			if(-e $rfile_BDH1){
				print " - once remove [blast_BDH]\n";
				system("rm -R blast_BDH");
			}
		}
		if(! -e $rfile_BDH1){
			my $cmd_BDH = "$bin_BBDH blast_BDH $qfasta $dbfasta $combined_qgff $dbgff $combined_qprot $dbprotein 2 30 $pcnt_th $pcov_th $eval_th $cpu_blast y pipe n > ./log/log_BDH.txt 2>&1";
			if(system($cmd_BDH) != 0){
				print "! blast_BDH failed, abort script...\n";
				die;
			}
		}
		else{
			print "! [$rfile_BDH1] already exists, skip...\n";
		}
	}
	else{
		print "! no sequence found in [$combined_qprot]\n";
	}
	
	print "! revise Gff based on BDH result...\n";
	my $htopcand = ReviseGff_byBDH($rfile_BDH1, $combined_qgff, $revised_qgff);
#	$htopcand->{$gid}{len} = $Data[6];
#	$htopcand->{$gid}{similarity} = $Data[10];
#	$htopcand->{$gid}{gap} = 100 - $Data[12];
#	$htopcand->{$gid}{score} = $Data[14];
#	$htopcand->{$gid}{qprotlen} = $hBDH->{qprotlen}{$gid};
	
#	my $cmd14 = "$hbinpath->{gffread} $revised_qgff -g $qfasta -y $revised_qprot > /dev/null 2>&1";
	my $cmd14 = "$bin_gff2prot $qfasta $revised_qgff $revised_qprot";
	system($cmd14);
#	Purse_seq($revised_qprot);		#remove "." attacched to the end of protein sequences
	
	print "! summarizing results...\n";
#	my $hralign = Txt2hash($tmp_ralign);
#	my $summary = "directory\tdbfasta\tseqid\tpos0\tpos1\tsp pos0\tsp pos1\tqfasta\tqfasta seqid\tpos0\tpos1\tgroup\thit sp0\thit sp1\tbp (sp0-sp1)\tbp hit (sp0-sp1)\tratio\tjudge\tgene_id\tprotein coding\tq protein length (predicted)\talign length\t\%similarity\t\%gap\tscore\tnum_NNN (qseq)\n";
	my $summary;
	my $summary_BDH;
	my $cnt_topcand1 = 0;
	my $cnt_topcand2 = 0;
	my $cnt_topcand3 = 0;
	my $cnt_topcand4 = 0;
	
	@G = sort {$a cmp $b} @G;
	
	foreach my $gid (@G){
#		if($hprevr->{$gid}){
#			$summary .= $hprevr->{$gid}."\n";
#			$summary_BDH .= $hprevr->{$gid}."\n";
#		}
#		elsif($hpav->{$gid}{data}){
		if($hpav->{$gid}{data}){
			my $len_dbprot = "";
			if($dbg2t->{$gid}){
				my @dbTID = split(/\n/, $dbg2t->{$gid});
				foreach my $dbtid (@dbTID){
					if($hdbprot->{$dbtid}{len}){
						$len_dbprot .= $hdbprot->{$dbtid}{len}.",";
					}
				}
			}
			unless($len_dbprot){
				$len_dbprot = "-";
			}
			
			if($htopcand->{$gid}){
				$cnt_topcand1++;
				if($htopcand->{$gid}{coverage} && $htopcand->{$gid}{coverage} ne '-'){
					$summary .= $hpav->{$gid}{data}."\t".$hlist->{$gid}{protein_coding}."\t".$len_dbprot."\t".$htopcand->{$gid}{len}."\t".$htopcand->{$gid}{len}."\t".$htopcand->{$gid}{similarity}."\t".$htopcand->{$gid}{coverage}."\t".$htopcand->{$gid}{score}."\n";
					$summary_BDH .= $hpav->{$gid}{data}."\t".$hlist->{$gid}{protein_coding}."\t".$len_dbprot."\t".$htopcand->{$gid}{len}."\t".$htopcand->{$gid}{len}."\t".$htopcand->{$gid}{similarity}."\t".$htopcand->{$gid}{coverage}."\t".$htopcand->{$gid}{score}."\n";
					$cnt_topcand2++;
				}
				else{
					$summary .= $hpav->{$gid}{data}."\t".$hlist->{$gid}{protein_coding}."\t".$len_dbprot."\t".$htopcand->{$gid}{qprotlen}."\t.\t.\t.\t.\n";
					$summary_BDH .= $hpav->{$gid}{data}."\t".$hlist->{$gid}{protein_coding}."\t".$len_dbprot."\t".$htopcand->{$gid}{qprotlen}."\t.\t.\t.\t.\n";
					$cnt_topcand3++;
				}
			}
			else{
				$summary .= $hpav->{$gid}{data}."\t".$hlist->{$gid}{protein_coding}."\t".$len_dbprot."\t-\t.\t.\t.\t.\n";
				$summary_BDH .= $hpav->{$gid}{data}."\t".$hlist->{$gid}{protein_coding}."\t".$len_dbprot."\t-\t.\t.\t.\t.\n";
				$cnt_topcand4++;
			}
		}
	}
	
	if($summary){
		if($outplot){
			if(-e $combined_final){
				ADD($combined_final, $summary);
				print "! add result to [$combined_final]\n";
			}
			else{
				$summary = $header_summary.$summary;
				SAVE($combined_final, $summary);
				print "! output [$combined_final]\n";
			}
		}
		else{
			$summary = $header_summary.$summary;
			SAVE($rfile_final, $summary);
			print "! output [$rfile_final]\n";
			
			$summary_BDH = $header_summary.$summary_BDH;
			SAVE($rfile_final_BDH, $summary_BDH);
#			print "! output [$rfile_final_BDH]\n";
		}
	}
	
	print "! removing temp files...\n";
	for(my $i = 0; $i < $cpu; $i++){
		my $tmpdir = "./_temp_$i";
		if(-d $tmpdir){
			system("rm -R $tmpdir");
		}
	}
	
	Rm_files(\@I1);
	Mv_files(\@I4, "cmds");
	Rm_files(\@gabbage);
	RmBlasttmp("blast_BDH");
	RmBlasttmp("blast_nBDH");
	
	if($redo_examine){
		SAVE($redo_log, "--5 done.\n");
	}
}
else{
	print "! unable to find [$afile], no candidate found...\n";
	if(! $outplot){
		system("touch $rfile_final");
	}
	goto END;
}

if($outplot){
	Combine_text($afile, $catafile);
	
	if(-e $qfasta){
		system("rm $qfasta");
	}
	if(-d $qpref){
		system("cp $qpref/*.png ./");
		system("rm -R $qpref");
	}
}
else{
	if(-e $qfasta){
		system("rm $qfasta");
	}
}

SORT:{
	Sort_and_save($rfile_final);
}

END:{
	if($outplot && -e $genelist){
		system("rm $genelist");
	}
	
#	print "\n! End of script.\n\n";
	my $end = 1;
}


################################################################################
#-------------------------------------------------------------------------------
sub RmBlasttmp{
my $dir = shift;

if(-d $dir){
	my @GB0 = glob("$dir/*");
	foreach my $file (@GB0){
		if($file =~ /BDconfidenthit/ || $file =~ /bN_db-/ || $file =~ /bP_db-/){
			next;
		}
		system("rm $file");
	}
}

}


#-------------------------------------------------------------------------------
sub Revise_afile_nBDH{
my ($afile, $rfile_nBDH1) = @_;

my $hash = {};
if(-e $rfile_nBDH1){
	open(my $fh, "<", $rfile_nBDH1) or return;
	my $cnt = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		my @A = split(/,/, $line);
		if($line =~ /ID/ && $line =~ /Asm2sv_partner/){
			next;
		}
		
#		if($A[5] eq 'BDH' || $A[5] eq 'BDBH'){
		if($A[5] eq 'BDBH'){
			if(! $hash->{$A[1]}){
				$hash->{$A[1]} = 1;
				$cnt++;
			}
		}
	}
	close $fh;
#	print "! [$cnt] BDH found\n";
}

my $r = "";
if(-e $afile){
	open(my $fh, "<", $afile) or die;
	my $cnt_all = 0;
	my $cnt_present = 0;
	my $cnt_bdh = 0;
	my $cnt_cnv = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		my @A = split(/\t/, $line);
		if($A[0] eq 'directory'){
			$r .= $line."\n";
			next;
		}
		
		if($A[22] =~ /present/ || $A[22] =~ /insertion/){
			my $gid = $A[23];
			if($hash->{$gid}){
				if($A[22] =~ / but not n-BDBH/){
					$A[22] =~ s/ but not n-BDBH//g;
				}
				$r .= join("\t", @A)."\n";
				$cnt_bdh++;
			}
			else{
				if($A[22] =~ / but not n-BDBH/){
					$A[22] =~ s/ but not n-BDBH//g;
					$A[22] .= " but not n-BDBH";
				}
				else{
	#				$A[22] .= " but not n-BDH";
					$A[22] .= " but not n-BDBH";
				}
				$r .= join("\t", @A)."\n";
				$cnt_cnv++;
			}
			$cnt_present++;
		}
		else{
			$r .= $line."\n";
		}
		$cnt_all++;
	}
	close $fh;
	print " - [$cnt_all] gene entries in [$afile]...\n";
	print " - [$cnt_bdh] hit with BDBH in query genome\n";
	print " - [$cnt_cnv] failed to find expected blast-n hit\n";
}

if($r){
	open(my $rfh, ">", $afile);
	print $rfh $r;
	close $rfh;
}

}


#-------------------------------------------------------------------------------
sub Convert2eachseq{
my ($afile, $hqseq, $rfile, $qpref) = @_;

my $rfasta = "";
if(-e $afile){
	print "! preparing nucleotide sequences for each hit region...\n";
	open(my $fh, "<", $afile) or die;
	my $cnt = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		my @A = split(/\t/, $line);
		if($A[0] eq 'directory'){
			next;
		}
		
		if($A[22] ne 'absent' && $A[12] ne '-' && $A[13] ne '-' && $A[12] =~ /\d/ && $A[13] =~ /\d/){
			my $qseqid = $A[8];
			my $qpos0 = $A[12];
			my $qpos1 = $A[13];
			my $eachseq = substr($hqseq->{$qseqid}{seq}, $qpos0 - 1, $qpos1 - $qpos0 + 1);
			$rfasta .= ">".$qpref."__".$A[23]."\n".$eachseq."\n";
			$cnt++;
		}
	}
	close $fh;
	print " - [$cnt] found except 'absent' flag\n";
}

if($rfasta){
	open(my $rfh, ">", $rfile);
	print $rfh $rfasta;
	close $rfh;
}

}


#-------------------------------------------------------------------------------
sub Convert2eachseq_Gff{
my ($gff, $hdbseq, $rfile) = @_;

my $rfasta = "";
if(-e $gff){
	print "! preparing nucleotide sequences for each gene region...\n";
	open(my $fh, "<", $gff) or die;
	my $cnt = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		my @A = split(/\t/, $line);
		if($A[0] && $A[8] && $A[2] && $A[2] eq 'gene'){
			my @tag = split(/\;/, $A[8]);
			my $gid;
			foreach my $val (@tag){
				if($val =~ /ID\=/){
					$gid = $val;
					$gid =~ s/ID\=//;
					last;
				}
			}
			
			if($gid){
				my $qpos0 = $A[3];
				my $qpos1 = $A[4];
				my $eachseq = substr($hdbseq->{$A[0]}{seq}, $qpos0 - 1, $qpos1 - $qpos0 + 1);
				$rfasta .= ">".$gid."\n".$eachseq."\n";
				$cnt++;
			}
		}
	}
	close $fh;
	print " - [$cnt] gene entries\n";
}

if($rfasta){
	open(my $rfh, ">", $rfile);
	print $rfh $rfasta;
	close $rfh;
}

}


#-----------------------------------------------------------
sub ReviseGff_byBDH{
my ($bdh, $gff, $rgff) = @_;

open(my $fh1, "<", $bdh);
my $cnt = 0;
my $hBDH = {};
while(my $line = <$fh1>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($cnt == 0){
		$cnt++;
		next;
	}
	
	my @A = split(/,/, $line);
#	my $r1 = "ID 1=[".$gprefix1."],gene_id [1],Asm2sv_partner [1],ID 2=[".$prefix2."],gene_id [2],Asm2sv_BDH_parter,[1] length,[2] length,[1] alignment length,[2] alignment length,[1] \%identity,[2] \%identity,[1] \%coverage,[2] \%coverage,[1] score,[2] score,[1] Eval,[2] Eval,k,l\n";
	
	if($A[5]){
		if($A[5] eq 'BDH' && ! $hBDH->{tid_BDBH}{$A[0]}){
			if($A[14] ne '-' && $A[14] =~ /\d/){
				$hBDH->{tid}{$A[0]} = "BDH";
				$hBDH->{gid}{$A[1]} += 1;
				$hBDH->{query}{$A[1]} = $A[2];
				$hBDH->{data}{$A[1]}{$A[0]} = $line;
				$hBDH->{score_tmp}{$A[1]}{$A[0]} = $A[14];
			}
		}
		if($A[5] eq 'BDBH'){
			if($A[14] ne '-' && $A[14] =~ /\d/){
				$hBDH->{tid_BDBH}{$A[0]} = 1;
				$hBDH->{gid_BDBH}{$A[1]} += 1;
				
				$hBDH->{tid}{$A[0]} = "BDBH";
				$hBDH->{gid}{$A[1]} += 1;
				$hBDH->{query}{$A[1]} = $A[2];
				$hBDH->{data}{$A[1]}{$A[0]} = $line;
				$hBDH->{score_tmp}{$A[1]}{$A[0]} = $A[14] + 1000000;
			}
		}
	}
	if($A[6] ne '-' && $A[6] =~ /\d/){
		$hBDH->{qprotlen}{$A[1]} .= $A[6].",";
	}
	$cnt++;
}
close $fh1;

open(my $fh2, "<", $gff);
$cnt = 0;
my @GID;
my $hash = {};
my $g2t = {};
my $t2g = {};
my $AoGID = [];
while(my $line = <$fh2>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /\#/){
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
		
		my @tmp;
		push(@tmp, $gid);
		push(@tmp, $A[0]);
		push(@tmp, $A[3]);
		push(@{$AoGID}, \@tmp);
		$hash->{gene}{$gid} = $line;
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
		}
		
		$g2t->{$gid} .= $tid."\n";
		$t2g->{$tid} = $gid;
		$hash->{pos0}{$tid} = $A[3];
		$hash->{pos1}{$tid} = $A[4];
		$hash->{num_variant}{$gid} += 1;
		
		if($hBDH->{tid_BDBH}{$tid}){
			$A[8] .= "\;status=BDBH";
		}
		elsif($hBDH->{tid}{$tid}){
			$A[8] .= "\;status=BDH";
		}
		else{
			$A[8] .= "\;status=false";
		}
		$hash->{transcript}{$tid} .= join("\t", @A)."\n";
	}
	elsif($A[2] eq 'CDS'){
		my $tid;
		foreach my $val (@tag){
			if($val =~ /Parent\=/){
				$tid = $val;
				$tid =~ s/Parent\=//;
			}
		}
		
		if(! $tid || ! $t2g->{$tid}){
			next;
		}
		
		my $gid = $t2g->{$tid};
		$hash->{num_CDS}{$gid} += 1;
		$hash->{transcript}{$tid} .= $line."\n";
	}
	else{
		my $tid;
		foreach my $val (@tag){
			if($val =~ /Parent\=/){
				$tid = $val;
				$tid =~ s/Parent\=//;
			}
		}
		
		if(! $tid || ! $t2g->{$tid}){
			next;
		}
		
		my $gid = $t2g->{$tid};
		$hash->{transcript}{$tid} .= $line."\n";
	}
}
close $fh2;

@{$AoGID} = sort {$a->[2] <=> $b->[2]} @{$AoGID};
@{$AoGID} = sort {$a->[1] cmp $b->[1]} @{$AoGID};

my $r = "##gff-version 3\n";
my $cnt_bdh = 0;
my $cnt_homolog = 0;
my $cnt_discard = 0;
my $cnt_withtid = 0;
my $cnt_wotid = 0;
my $cnt_topcand = 0;
my $htopcand = {};
foreach my $tmp (@{$AoGID}){
	my $gid = $tmp->[0];
	
	unless($g2t->{$gid}){
		$cnt_wotid++;
		next;
	}
	$cnt_withtid++;
	
	my @TID = split(/\n/, $g2t->{$gid});
	my $AoS = [];
	foreach my $tid (@TID){
		if($hBDH->{data}{$gid}{$tid}){
			my @tmp;
			push(@tmp, $hBDH->{data}{$gid}{$tid});
			push(@tmp, $hBDH->{score_tmp}{$gid}{$tid});
			push(@tmp, $hBDH->{query}{$gid});
			push(@tmp, $tid);
			push(@tmp, $hBDH->{tid}{$tid});
			push(@{$AoS}, \@tmp);
		}
	}
	if(@{$AoS}){
		@{$AoS} = sort {$b->[1] <=> $a->[1]} @{$AoS};
		
	#	$summary .= $hpav->{$gid}{data}."\t".$hlist->{$gid}{protein_coding}."\t".$len_dbprot."\t".$len_qprot."\t".$htmp->{len}."\t".$htmp->{similarity}."\t".$htmp->{gap}."\t".$htmp->{score}."\n";
		
		my $asm2sv_query = $AoS->[0][2];
		my @Data = split(/,/, $AoS->[0][0]);
		$htopcand->{$asm2sv_query}{len} = $Data[6];
		$htopcand->{$asm2sv_query}{similarity} = $Data[10];
		$htopcand->{$asm2sv_query}{score} = $Data[14];
		$htopcand->{$asm2sv_query}{qprotlen} = $hBDH->{qprotlen}{$gid};
		if($Data[12] ne '-'){
			$htopcand->{$asm2sv_query}{coverage} = $Data[12];
		}
		else{
			$htopcand->{$asm2sv_query}{coverage} = "-";
		}
		$cnt_topcand++;
	}
	
	if($hBDH->{gid}{$gid}){
		$htopcand->{$gid}{BDH} = "true";
		my @Npos;
		my $each_gff = "";
		my $num_AoS = @{$AoS};
		my $num_BDBH = 0;
		my $numt = 1;
		foreach my $tmp (@{$AoS}){
			my $tid = $tmp->[3];
			
			my $numt0 = $tid;
			$numt0 =~ s/$gid//;
			my $numtr = ".t".$numt;
			my $tidr = $gid.$numtr;
			
			if($tmp->[4] eq 'BDBH'){
				push(@Npos, $hash->{pos0}{$tid});
				push(@Npos, $hash->{pos1}{$tid});
				$hash->{transcript}{$tid} =~ s/$tid/$tidr/g;
				$each_gff .= $hash->{transcript}{$tid};
				$num_BDBH++;
				$numt++;
			}
			elsif($tmp->[4] eq 'BDH' && $num_BDBH == 0){		# only keep BDBH
#			elsif($tmp->[4] eq 'BDH'){
				push(@Npos, $hash->{pos0}{$tid});
				push(@Npos, $hash->{pos1}{$tid});
				$hash->{transcript}{$tid} =~ s/$tid/$tidr/g;
				$each_gff .= $hash->{transcript}{$tid};
				$numt++;
			}
			else{
				$cnt_discard++;
			}
		}
		
		@Npos = sort {$a <=> $b} @Npos;
		my $npos0 = $Npos[0];
		@Npos = sort {$b <=> $a} @Npos;
		my $npos1 = $Npos[0];
		
		my @G = split(/\t/, $hash->{gene}{$gid});
		$G[3] = $npos0;
		$G[4] = $npos1;
		if($hBDH->{gid_BDBH}{$gid}){
			$G[8] .= "\;status=BDBH";
		}
		else{
			$G[8] .= "\;status=BDH";
		}
		$r .= join("\t", @G)."\n";
		$r .= $each_gff;
		$cnt_bdh++;
	}
	else{
		$htopcand->{$gid}{BDH} = "false";
		my @Npos;
		my $each_gff = "";
		foreach my $tid (@TID){
#			my $ntid = "homolog_".$tid;
#			$hash->{transcript}{$tid} =~ s/$tid/$ntid/g;
			
			push(@Npos, $hash->{pos0}{$tid});
			push(@Npos, $hash->{pos1}{$tid});
			$each_gff .= $hash->{transcript}{$tid};
		}
		
		@Npos = sort {$a <=> $b} @Npos;
		my $npos0 = $Npos[0];
		@Npos = sort {$b <=> $a} @Npos;
		my $npos1 = $Npos[0];
		
#		my $ngid = "homolog_".$gid;
		my @G = split(/\t/, $hash->{gene}{$gid});
		$G[3] = $npos0;
		$G[4] = $npos1;
#		$G[8] =~ s/$gid/$ngid/;
		$G[8] .= "\;status=false";
		$r .= join("\t", @G)."\n";
		$r .= $each_gff;
		$cnt_homolog++;
	}
}

print " - [$cnt_withtid / $cnt_wotid] with/without transcript\n";
print " - [$cnt_topcand] with candidate\n";
print " - [$cnt_bdh] bidirectional hit (BDH)\n";
print " - [$cnt_homolog] are not BDH\n";
#print " - [$cnt_discard] transcript entries discarded\n";

open(my $rfh, ">", $rgff);
print $rfh $r;
close $rfh;
print "! output [$rgff]\n";

return $htopcand;
}


#-----------------------------------------------------------
sub Purse_seq{
my $file = shift;

#>mRNA1 gene=species_seq1_G_594180

#print "! pursing [$file]...\n";
open(my $fh, "<", $file) or die;
my $id;
my $seq;
my $modseq;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /\>/){
		if($seq){
			$modseq .= ">".$id."\n".$seq."\n";
		}
		
		$id = $line;
		$id =~ s/\>//;
		my @tmp = split(/\s/, $id);
		$id = $tmp[0];
		$seq = "";
	}
	else{
		if($line =~ /\./){
			my @tmp = split(/\./, $line);
			if($tmp[0]){
				$line = $tmp[0];
			}
			else{
				$line = "";
			}
		}
		
		if($line){
			$seq .= $line;
		}
	}
}
close $fh;

if($seq){
	$modseq .= ">".$id."\n".$seq."\n";
}

if($modseq){
	open(my $rfh, ">", $file);
	print $rfh $modseq;
	close $rfh;
#	print "! output modified fasta : [$file]\n";
}

}


#-------------------------------------------------------------------------------
sub Read_and_makehist{
my ($file, $hist) = @_;

open(my $fh, "<", $file);
my $cnt = 0;
my $r;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	unless($line){
		next;
	}
	if($cnt == 0){
		$cnt++;
		next;
	}
	
	my @A = split(/\t/, $line);
	if($A[23]){
		$r .= $A[23]."\n";
	}
	$cnt++;
}
close $fh;

if($r){
	open(my $rfh, ">", $hist);
	print $rfh $r;
	close $rfh;
}

}


#-------------------------------------------------------------------------------
sub JudgeJob{
my $file = shift;

my $n = 0;
if(-e $file){
	open(my $fh, "<", $file);
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line && $line =~ /all_done/){
			$n = 1;
		}
	}
	close $fh;
}

return $n;
}


#-------------------------------------------------------------------------------
sub SortPsl{
my $psl = shift;
my $rfile = shift;

open(my $fh, "<", $psl);
my $AoA = [];
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	unless($line){
		next;
	}
	
	my @A = split(/\t/, $line);
	push(@{$AoA}, \@A);
	$cnt++;
}
close $fh;

@{$AoA} = sort {$a->[15] <=> $b->[15]} @{$AoA};
@{$AoA} = sort {$a->[13] cmp $b->[13]} @{$AoA};

my $r;
foreach my $A (@{$AoA}){
	$r .= join("\t", @{$A})."\n";
}

open(my $rfh, ">", $rfile);
print $rfh $r;
close $rfh;

}


#-------------------------------------------------------------------------------
sub Combine_text{
my $afile = shift;
my $rfile = shift;

open(my $fh, "<", $afile);
my $r;
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($cnt == 0){
		if(! -e $rfile){
			$r .= $line."\n";
		}
	}
	else{
		$r .= $line."\n";
	}
	$cnt++;
}
close $fh;

open(my $rfh, ">>", $rfile);
print $rfh $r;
close $rfh;

}


#-------------------------------------------------------------------------------
sub Sort_and_save{
my $file = shift;

#my $backup = $file.".backup";
#system("cp $file $backup");

if(-e $file){
	open(my $fh, "<", $file);
	my $AoA = [];
	my $AoAnc = [];
	my $cnt = 0;
	my $r;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		unless($line){
			next;
		}
		if($cnt == 0){
			$r .= $line."\n";
		}
		else{
			my @A = split(/\t/, $line);
			my @tmp;
			push(@tmp, $A[2]);
			push(@tmp, $A[3]);
			push(@tmp, $line);
			
			if($A[2] =~ /ch/i || $A[2] =~ /Chr/i){
				push(@{$AoA}, \@tmp);
			}
			else{
				push(@{$AoAnc}, \@tmp);
			}
		}
		$cnt++;
	}
	close $fh;

	if(@{$AoA}){
		@{$AoA} = sort {$a->[1] <=> $b->[1]} @{$AoA};
		@{$AoA} = sort {$a->[0] cmp $b->[0]} @{$AoA};
		
		foreach my $tmp (@{$AoA}){
			$r .= $tmp->[2]."\n";
		}
	}

	if(@{$AoAnc}){
		@{$AoAnc} = sort {$a->[1] <=> $b->[1]} @{$AoAnc};
		@{$AoAnc} = sort {$a->[0] cmp $b->[0]} @{$AoAnc};
		
		foreach my $tmp (@{$AoAnc}){
			$r .= $tmp->[2]."\n";
		}
	}

	if($r){
		if(-e $file){
			system("rm $file");
		}
		open(my $rfh, ">", $file);
		print $rfh $r;
		close $rfh;
	}
}

}


#-------------------------------------------------------------------------------
sub Count_seqlen{
my $file = shift;

open(my $fh, "<", $file);
my $seq;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line && $line !~ />/){
		$seq .= $line;
	}
}
close $fh;

my $len = 0;
if($seq){
	$len = length($seq);
}

return $len;
}


#-------------------------------------------------------------------------------
sub Count_NNN{
my $seq = shift;

$seq =~ s/\n//g;
$seq =~ s/A//gi;
$seq =~ s/C//gi;
$seq =~ s/G//gi;
$seq =~ s/T//gi;
$seq =~ s/U//gi;

my $len = 0;
if($seq){
	$len = length($seq);
}

return $len;
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
sub Txt2hash{
my $file = shift;

my $hash = {};
if(-e $file){
	open(my $fh, "<", $file);
	while(my $l = <$fh>){
		$l =~ s/\n//;
		$l =~ s/\r//;
		
		my @A = split(/\t/, $l);
		$hash->{$A[0]}{align} = $A[1];
		$hash->{$A[0]}{cnt_NNN} = $A[2];
		$hash->{$A[0]}{sub_qprot} = $A[3];
	}
	close $fh;
}

return $hash;
}


#-------------------------------------------------------------------------------
sub Mv_files{
my $A = shift;
my $dir = shift;

foreach my $file (@{$A}){
	if(-e $file){
		system("mv $file $dir");
	}
}

}


#-------------------------------------------------------------------------------
sub Rm_files{
my $A = shift;

foreach my $file (@{$A}){
	if(-e $file){
		system("rm $file");
	}
}

}


#-------------------------------------------------------------------------------
sub Extract_and_save{
my $rfile = shift;
my $id = shift;
my $seq = shift;
my $p0 = shift;
my $p1 = shift;

if($seq){
	my $subseq = substr($seq, $p0 - 1, abs($p1 - $p0));
	my $subfasta = ">".$id."\n".$subseq;
	
	unless(-e $rfile){
		open(my $rfh, ">", $rfile) or return 'null';
		print $rfh $subfasta;
		close $rfh;
	}
	return $subseq;
}

}


#-------------------------------------------------------------------------------
sub Cmdexe_unless_seq{
my $cmd = shift;
my $file = shift;

unless(-e $file){
	system($cmd);
}

}


#-------------------------------------------------------------------------------
sub Update_prev_results{
my $file = shift;

open(my $fh, "<", $file) or die;
my $hash = {};
my $selected = "";
my $cnt_all = 0;
my $cnt_selected = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /directory/ && $line =~ /dbfasta/ && $line =~ /seqid/){
		$selected .= $line."\n";
		next;
	}
	
	my @A = split(/\t/, $line);
	my $numA = @A;
	if($numA > 18){
#		if($A[22] eq 'present' && $A[19] ne '0'){
		if($A[22] =~ /present/){
			my $promoterutr_judge = 0;
			my $lowth = 0.95;
			my $upth = 1.05;
			
			if(defined $A[25] && $A[25] =~ /\;/){
				# ratio_hit; ratio_span; rev_ratio_span; P2Gdist; span_pos0; span_pos1
				my @tmpA = split(/\;/, $A[25]);
				if($tmpA[1] && $tmpA[2] && $tmpA[1] ne 'null' && $tmpA[2] ne 'null'){
					if($lowth <= $tmpA[1] && $lowth <= $tmpA[2] && $tmpA[1] <= $upth && $tmpA[2] <= $upth){
						$promoterutr_judge++;
					}
				}
			}
			if(defined $A[27] && $A[27] =~ /\;/){
				# ratio_hit; ratio_span; rev_ratio_span; P2Gdist; span_pos0; span_pos1
				my @tmpA = split(/\;/, $A[27]);
				if($tmpA[1] && $tmpA[2] && $tmpA[1] ne 'null' && $tmpA[2] ne 'null'){
					if($lowth <= $tmpA[1] && $lowth <= $tmpA[2] && $tmpA[1] <= $upth && $tmpA[2] <= $upth){
						$promoterutr_judge++;
					}
				}
			}
			
#			if($promoterutr_judge == 2){		# not_use
				$selected .= $line."\n";
				$hash->{$A[23]} = 1;
				$cnt_selected++;
#			}
		}
	}
	$cnt_all++;
}
close $fh;

my $cnt_rm = $cnt_all - $cnt_selected;
print "! [$cnt_selected] with exact 'P' tag (keep)\n";
print "! [$cnt_rm] once removed from [$file]\n";

for(my $i = 1; $i < 100; $i++){
	my $backup = $file.".backup".$i;
	if(! -e $backup){
		system("cp $file $backup");
		last;
	}
}

open(my $rfh, ">", $file);
print $rfh $selected;
close $rfh;

return $hash;
}


#-------------------------------------------------------------------------------
sub Check_prev_results{
my $file = shift;

open(my $fh, "<", $file) or die;
my $hash = {};
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /directory/ && $line =~ /dbfasta/ && $line =~ /seqid/){
		next;
	}
	
	my @A = split(/\t/, $line);
	my $numA = @A;
	if($numA > 18){
		$hash->{$A[23]} = $line;
	}
}
close $fh;

return $hash;
}


#-------------------------------------------------------------------------------
sub Check_prev_orfsearch{
my $file = shift;

open(my $fh, "<", $file) or die;
my $hash = {};
my $norf = 0;
my $nnot = 0;
my $revised = "";
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /directory/ && $line =~ /dbfasta/ && $line =~ /seqid/){
		$revised .= $line."\n";
		next;
	}
	
	my @A = split(/\t/, $line);
	my $numA = @A;
	if($A[38]){
#		if($A[38] eq '-' || $A[38] eq '.'){			# ORF information NOT attached
		if($A[38] eq '-'){							# ORF information NOT attached
			$nnot++;
		}
		else{										# ORF information already attached
			$norf++;
			$hash->{$A[23]} = $line;
		}
	}
	else{
		$nnot++;
	}
	
	if($A[22] =~ / but not n-BDBH/){
		$A[22] =~ s/ but not n-BDBH//g;
		$A[22] .= " but not n-BDBH";
	}
	$revised .= join("\t", @A)."\n";
}
close $fh;

if($revised){
	open(my $rfh, ">", $file);
	print $rfh $revised;
	close $rfh;
}

my $ratio = 0;
if( ($nnot + $norf) > 0){
	$ratio = $norf / ($nnot + $norf);
}

my $rh = {};
$rh->{norf} = $norf;
$rh->{ratio_norf} = $ratio;
$rh->{hash} = $hash;

return $rh;
}


#-------------------------------------------------------------------------------
sub Open_pav_results{
my $file = shift;

print "! open PAV result [$file] ...\n";
my $hash = {};
open(my $fh, "<", $file) or return $hash;
my $cnt = 0;
my $cntP = 0;
my $cntA = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /directory/ && $line =~ /dbfasta/ && $line =~ /seqid/){
		next;
	}
	
	my @A = split(/\t/, $line);
	if($A[18]){
		my $gid = $A[23];
		my $pav_status = $A[22];
		
		if($A[22] =~ / but not n-BDBH/){
			$A[22] =~ s/ but not n-BDBH//g;
			$A[22] .= " but not n-BDBH";
		}
		
		if($pav_status && $pav_status =~ /present/i){
			$hash->{$gid}{status} = "P";
			$hash->{$gid}{dbsid} = $A[2];
			$hash->{$gid}{dbpos0} = $A[5];
			$hash->{$gid}{dbpos1} = $A[6];
			$hash->{$gid}{qsid} = $A[8];
			$hash->{$gid}{qpos0} = $A[12];
			$hash->{$gid}{qpos1} = $A[13];
			$hash->{$gid}{data} = join("\t", @A);
			$cntP++;
		}
		elsif($pav_status && $pav_status =~ /insert/i){
			$hash->{$gid}{status} = "P";
			$hash->{$gid}{dbsid} = $A[2];
			$hash->{$gid}{dbpos0} = $A[5];
			$hash->{$gid}{dbpos1} = $A[6];
			$hash->{$gid}{qsid} = $A[8];
			$hash->{$gid}{qpos0} = $A[12];
			$hash->{$gid}{qpos1} = $A[13];
			$hash->{$gid}{data} = join("\t", @A);
			$cntP++;
		}
		else{
			$hash->{$gid}{status} = "A";
			$hash->{$gid}{data} = join("\t", @A);
			$cntA++;
		}
		
		$cnt++;
	}
}
close $fh;

print "! total [$cnt] lines, P=[$cntP], other=[$cntA]\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Count_seq_fasta{
my $file = shift;

#print "! reading [$file]...\n";
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
	$total_len += $len;
	push(@SID, $ID);
}

my $numID = @SID;

#print "! total [$numID] sequence ID, [$total_len] bp\n";

return $numID;
}


#-------------------------------------------------------------------------------
sub Open_fasta_as_hash{
my $file = shift;

print "! open [$file] as hash ...\n";
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

print "! total [$numID] sequence ID, [$total_len] bp\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Gff_to_hash{
my $gff3 = shift;

print "! reading [$gff3]...\n";
open(my $fh, "<", $gff3) or die;
my $hash = {};
my $cnt = 0;
my $gcnt = 0;
my $tcnt = 0;
my $pcnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	$line =~ s/\"//g;
	
	if($line =~ /\#/){
		next;
	}
	
	my @A = split(/\t/, $line);
	unless($A[8]){
		next;
	}
	my @tag = split(/\;/, $A[8]);
	my $gid;
	my $tid;
	my $name_evm;
	$cnt++;
	
#	if($A[0] =~ /scaffold/ || $A[0] =~ /unanchored/ || $A[0] =~ /mitochon/ || $A[0] =~ /chloro/){
#		next;
#	}
	
	if($A[2] eq 'gene'){
		foreach my $val (@tag){
			if($val =~ /ID\=/ && $val !~ /Alt_ID\=/){
				$gid = $val;
				$gid =~ s/ID\=//;
			}
			elsif($val =~ /Name\=EVM/){
				$name_evm = $val;
				$name_evm =~ s/Name\=EVM//;
			}
		}
		
#		if($name_evm){
#			$A[8] = "ID\=".$gid.";Note\=EVM";
#		}
#		else{
			$A[8] = "ID\=".$gid;
#		}
		
		if($gid){
			$hash->{hgff}{$gid}{Chr} = $A[0];		#1	seqid
			$hash->{hgff}{$gid}{pos0} = $A[3];		#2	pos0
			$hash->{hgff}{$gid}{pos1} = $A[4];		#3	pos1
			$hash->{hgff}{$gid}{strand} = $A[6];	#4	strand
			$gcnt++;
		}
	}
	elsif($A[2] eq 'mRNA' || $A[2] eq 'transcript'){
		foreach my $val (@tag){
			if($val =~ /ID\=/ && $val !~ /Alt_ID\=/){
				$tid = $val;
				$tid =~ s/ID\=//;
			}
			elsif($val =~ /Parent\=/){
				$gid = $val;
				$gid =~ s/Parent\=//;
			}
		}
		
		if($gid && $tid){
			$hash->{t2g}{$tid} = $gid;
			$hash->{g2t}{$gid} .= $tid."\n";
			$tcnt++;
		}
	}
	elsif($A[2] eq 'cds' || $A[2] eq 'CDS'){
		foreach my $val (@tag){
			if($val =~ /Parent\=/){
				$tid = $val;
				$tid =~ s/Parent\=//;
			}
		}
		
		if($tid){
			$hash->{CDS}{$tid} = "true";
			$pcnt++;
		}
	}
}
close $fh;

print "! [$gcnt] genes, [$tcnt] transcripts, [$pcnt] CDS lines\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Check_previous{
my $file = shift;
my $id = shift;

my $sw = "false";
if(-e $file){
	open(my $fh, "<", $file) or die;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		unless($line){
			next;
		}
		if($line =~ /\t/){
			my @A = split(/\t/, $line);
			$line = join("_", @A);
		}
		
		if($id eq $line){
			$sw = "true";
		}
	}
	close $fh;
}

return $sw;
}


#-------------------------------------------------------------------------------
sub Batch_exe{
my $script = shift;
my $j = shift;
my $cnt = shift;

if($script && -e $script){
	print " thread [$j] : [$cnt] queries...\n";
	if(system("bash $script") == 0){
		print " thread [$j] : completed.\n";
		
		unless(-e "cmds"){
			system("mkdir cmds");
		}
		system("mv $script cmds");
	}
	else{
		print "! thread [$j] : failed\n";
	}
}

}


#-------------------------------------------------------------------------------
sub Batch_exe2{
my $script = shift;
my $j = shift;
my $log_th = shift;

if($script && -e $script){
	print " thread [$j] : started...\n";
	if(system("bash $script > $log_th 2>&1") == 0){
		print " thread [$j] : completed.\n";
		
		unless(-e "cmds"){
			system("mkdir cmds");
		}
		system("mv $script cmds");
	}
	else{
		print "! thread [$j] : failed | log = [$log_th]\n";
	}
}

}


#-------------------------------------------------------------------------------
sub Get_genelist{
my $prefix = shift;
my $list = shift;
my $rfile = shift;
my $hgff = shift;
my $dbg2t = shift;
my $dbt2CDS = shift;

if(-e $list){
	my $hash = {};
	print "! reading [$list]...\n";
	open(my $fh, "<", $list) or die;
	my $cnt_p = 0;
	my $cnt_a = 0;
	my $log;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		if($line =~ /gid,chr,pos/){
			next;
		}
		
		my @A;
		if($line =~ /,/){
			@A = split(/,/, $line);
		}
		elsif($line =~ /\t/){
			@A = split(/\t/, $line);
		}
		else{
			$A[0] = $line;
		}
		my $gid = $A[0];
		
		if($gid && $hgff->{$gid}){
			$log .= $gid."\t".$hgff->{$gid}{Chr}."\t".$hgff->{$gid}{pos0}."\t".$hgff->{$gid}{pos1}."\n";
			
			if($dbg2t->{$gid}){
				my @TID = split(/\n/, $dbg2t->{$gid});
				foreach my $tid (@TID){
					if($dbt2CDS->{$tid}){
						$hgff->{$gid}{protein_coding} = "true";
						last;
					}
				}
			}
			
			unless($hgff->{$gid}{protein_coding}){
				$hgff->{$gid}{protein_coding} = "false";
			}
			
			$hash->{$gid} = $hgff->{$gid};
			$cnt_p++;
		}
		else{
			$cnt_a++;
		}
	}
	close $fh;
	
	print "! [$cnt_p] found in gff\n";
	print "! [$cnt_a] not found in gff\n";
	
	if($log){
		open(my $rfh, ">", $rfile);
		print $rfh $log;
		close $rfh;
	}
	
	return $hash;
}

}


#-------------------------------------------------------------------------------
sub Read_repeat_type{
my $file = shift;
my $rfile = shift;

print "! reading [$file]...\n";
my @TYP;
if(-e $file){
	open(my $fh, "<", $file) or die;
	my $cnt = 0;
	while(my $line = <$fh>){
		$line =~ s/\n//;
		$line =~ s/\r//;
		
		my @A = split(/\t/, $line);
		unless($line){
			next;
		}
		elsif($line =~ /PRIMARY_CNT_END/){
			last;
		}
		
		if($A[0]){
			push(@TYP, $A[0]);
		}
	}
	close $fh;
}

@TYP = sort {$a cmp $b} @TYP;

return \@TYP;
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


#----------------------------------------------------------
sub SAVE{
my $file = shift;
my $str = shift;

if($str){
	open(my $fh, ">", $file) or die;
	print $fh $str;
	close $fh;
	
	#print "! output [$file]\n";
}

}


#----------------------------------------------------------
sub APPEND{
my $file = shift;
my $str = shift;

open(my $fh, ">>", $file) or die;
print $fh $str;
close $fh;

#print "! output [$file]\n";

}


#-------------------------------------------------------------------------------
sub Delete{
my $file = shift;

if(-e $file){
	system("rm $file");
}

}


#---------------------------------------------------------------------
sub CheckRAM{

my $tmp = "_tmp_ramcheck.txt";
my $ramsize = 0;
if(system("free -h > $tmp") == 0){
	if(-e $tmp){
		open(my $fh, "<", $tmp);
		while(my $line = <$fh>){
			$line =~ s/\n//;
			$line =~ s/\r//;
			if($line && $line =~ /Mem:/i){
				my @A = split(/\s/, $line);
				foreach my $val (@A){
					if($val && $val =~ /G/i){
						$ramsize = $val;
						$ramsize =~ s/Gi//i;
						$ramsize =~ s/G//i;
						$ramsize *= 1024;				# MB scale
						last;
					}
					elsif($val && $val =~ /T/i){
						$ramsize = $val;
						$ramsize =~ s/Ti//i;
						$ramsize =~ s/T//i;
						$ramsize *= 1024 * 1024;		# MB scale
						last;
					}
				}
			}
		}
		close $fh;
		
		system("rm $tmp");
	}
}

unless($ramsize){
	$ramsize = 128000;
}

return $ramsize;
}


#---------------------------------------------------------------------
sub Search_app2{
my $hbin = shift;
my $APP = shift;
my $hpath = shift;
my $script_path = shift;
my $pipeline = shift;
my $genelist = shift;
my $qfasta = shift;

$script_path =~ s/\/scripts//;
my @tmp = split(/\//, $genelist);
my $ntmp = @tmp;
$genelist = $tmp[$ntmp - 1];

my $tmp = "_Search_app_temp.$genelist.$qfasta.txt";
for(my $i = 0; $i < 100; $i++){
	if(-e $tmp){
		$tmp = "_Search_app_temp.$genelist.$qfasta.$i.txt";
	}
	if(! -e $tmp){
		last;
	}
}

unless($pipeline){
	print "! checking PATH of required programs...\n";
}
my $hash = {};
foreach my $app (@{$APP}){
	my $Bin = $hbin->{$app};
	if($hpath->{$app} && -e $hpath->{$app}){
		my @F = glob("$hpath->{$app}/*");
		foreach my $bin (@{$Bin}){
			foreach my $f (@F){
				if($f =~ /$bin/){
					$hash->{$bin} = $f;
				}
			}
		}
	}
	else{
		foreach my $bin (@{$Bin}){
			system("which $bin > $tmp 2>&1");
			unless($bin =~ /samtools/){
				system("which $bin.pl >> $tmp 2>&1");
			}
			system("find $script_path/bin | grep $bin >> $tmp 2>&1");
			
			open(my $fh, "<", $tmp);
			while(my $line = <$fh>){
				$line =~ s/\n//;
				$line =~ s/\r//;
				
				if($line && $line =~ /$bin/i && $line !~ /which: no/){
					if($line =~ /\//){
						my @D = split(/\//, $line);
						my $nD = @D;
						$D[$nD - 1] =~ s/\.pl//;
						if($D[$nD - 1] eq $bin){
							$hash->{$bin} = $line;
						}
					}
					else{
						$hash->{$bin} = $line;
					}
				}
			}
			close $fh;
			
			system("rm $tmp");
		}
	}
	
	foreach my $bin (@{$Bin}){
		if($hash->{$bin} && -e $hash->{$bin}){
			unless($pipeline){
				if($bin !~ /interpro/ && ! $pipeline){
					print " [$bin] = [$hash->{$bin}] <OK>\n";
				}
				elsif($bin =~ /interpro/){
					my $testcmd = "$hash->{$bin} --version > $tmp";
					system($testcmd);
					
					my $iprver;
					open(my $fh, "<", $tmp);
					while(my $line = <$fh>){
						$line =~ s/\n//;
						$line =~ s/\r//;
						
						if($line && $line =~ /InterProScan version/i){
							$line =~ s/InterProScan version//;
							$line =~ s/\s//g;
							my @tmp2 = split(/\-/, $line);
							if($tmp2[0] =~ /5\./){
								$iprver = $tmp2[0];
							}
						}
					}
					close $fh;
					
					system("rm $tmp");
					
					if(! $iprver || $iprver !~ /5\./){
						print " NOT found [$bin], please check PATH\n";
						$hash->{err} += 1;
					}
					elsif($iprver =~ /5\./){
						print " [$bin] = [$hash->{$bin}] <OK> (version $iprver)\n";
					}
				}
			}
		}
		else{
			print " NOT found [$bin], please check PATH\n";
			$hash->{err} += 1;
		}
	}
}

unless($pipeline){
	print "\n";
}

return $hash;
}



