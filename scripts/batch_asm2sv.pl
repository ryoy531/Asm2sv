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
$version .= "batch_asm2sv.pl version 1.01\n";
$version .= "last update: [2020\/2\/12]\n";
$version .= "copyright: ryoichi yano [ryoichiy104\@gmail.com]\n";

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
my $skip_orfsearch;
my $cluster;
my $help;
my $test;

GetOptions('--db=s' => \$dbfasta, '--gff=s' => \$dbgff, '--list=s' => \$genelist, '--query=s' => \$qfasta, '--alignpsl=s' => \$qpsl, '--output_dir=s' => \$dir, '--neighbor=i' => \$neighbor_tlen, '--cpu=i' => \$cpu, '--xthread=i' => \$cpu_gth, '--span=i' => \$initspan, '--plot' => \$outplot, '--test' => \$test, '--1' => \$cluster, '--wo_orfsearch' => \$skip_orfsearch, '--help' => \$help);

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
elsif($cpu > 96){
	$cpu = 96;
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

my $hbin = {};
foreach my $app (@APP){
	$hbin->{$app}[0] = $app;
}
my $hpath = {};
my $hbinpath = Search_app2($hbin, \@APP, $hpath, $script_path, "pipe", $genelist);
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
my $afile = "./$qpref/specified_region_hits_".$qpref.".tsv";
my $catafile = "specified_region_hits.tsv";

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
@G = sort {$a cmp $b} @G;
my $num = @G;
my $d = int($num / $cpu) + 1;

#---------------------------------------------------------------//
my $disable_skipaln = 1;
if(-e $afile && ! $outplot && $disable_skipaln eq '0'){
	print "! [$afile] already exists...\n";
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
		
		foreach my $gid (@G){
			if(! $hpreva->{$gid}){
				push(@Gnotyet, $gid);
			}
		}
		$numQL_Gnotyet = @Gnotyet;
		print "! [$numQL_Gnotyet] not yet analyzed...\n";
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
				
				my $hist = "hist.txt";
				my $identifer = $qpref."_".$dbfasta."_".$sid."_".$region."_".$specified;
				my $acheck = Check_previous($hist, $identifer);
				
				if($acheck eq 'false'){
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
						
						$cmd .= "echo \'$identifer\' >> hist.txt\n";
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
$header_summary .= "protein coding\tdb protein length\tq protein length (predicted)\talign length\t\%similarity\t\%gap\tscore\n";

if($skip_orfsearch){
	print "! --wo_orfsearch is invoked, skip ORF prediction...\n";
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
			print "! CPU threads [$cpu] > [$max_cpu]\n";
			$cpu = $max_cpu;
		}
	}
	
	my $hpav = Open_pav_results($afile);
	my $hqseq = Open_fasta_as_hash($qfasta);
	
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
	
	print "! preparing dataset for genome threader... (may take a while)\n";
	my $hprevr = {};
	my @Gnotyet;
	my $numQL_Gnotyet = @G;
	my $num_prevorf = 0;
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
	
	my $log_gth_hist = "log_gth_processed.txt";
	if(-e $log_gth_hist){
		system("rm $log_gth_hist");
	}
	
	my $combined_qprot = "predicted_".$qpref.".fasta";
	
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
	if($numQL_Gnotyet > 0 && $num_prevorf == 0){
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
			my $log_gth = "thread_gth_".$i.".log";
			my $script = "cmd_gth_".$qpref."_".$i.".sh";
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
								
								my $sub_dbseq = Extract_and_save($tmp_dbfasta, $dbsid, $hdbseq->{$dbsid}{seq}, $dbpos0, $dbpos1);
								my $sub_qseq = Extract_and_save($tmp_qfasta, $qsid, $hqseq->{$qsid}{seq}, $qpos0, $qpos1);
								my $len_sub_dbseq = length($sub_dbseq);
								my $len_sub_qseq = length($sub_qseq);
								
	#							$cmd .= "if test -e \"$mk_process\"; then\n";
	#							$cmd .= "\techo '$gid already processed' >> $log_gth 2>&1\n";
	#							$cmd .= "fi\n";
								
	#							if($sub_qseq && -e $tmp_dbprot && -e $tmp_qfasta){
								if(-e $tmp_dbprot && -e $tmp_qfasta){
									if(-e $bin_findorf){
										$cmd .= "perl $bin_findorf $hbinpath->{gth} $hbinpath->{samtools} $hbinpath->{gffread} $hbinpath->{blat} $hbinpath->{blat2hints} $hbinpath->{matcher} $gid $tmp_dbprot $tmp_dbtranscript $tmp_dbcds $tmp_qfasta $tmp_qgff $tmp_qprot $tmp_align $tmpdir $qpref $combined_qprot >> $log_gth 2>&1\n";
									}
									else{
										$cmd .= "if test ! -e \"$tmp_qgff\"; then\n";
										$cmd .= "\tif test -e \"$tmp_qfasta\"; then\n";
										$cmd .= "\t\t$hbinpath->{gth} -genomic $tmp_qfasta -protein $tmp_dbprot -gff3out -skipalignmentout -o $tmp_qgff >> $log_gth 2>&1\n";
										$cmd .= "\t\t$hbinpath->{samtools} faidx $tmp_qfasta >> $log_gth 2>&1\n";
										$cmd .= "\t\t$hbinpath->{gffread} $tmp_qgff -g $tmp_qfasta -y $tmp_qprot >> $log_gth 2>&1\n";
										$cmd .= "\t\t$hbinpath->{matcher} -asequence $tmp_dbprot -bsequence $tmp_qprot -outfile $tmp_align >> $log_gth 2>&1\n";
										
										$cmd .= "\t\tif test -e \"$tmp_dbprot\"; then\n";
										$cmd .= "\t\t\trm $tmp_dbprot\n";
										$cmd .= "\t\tfi\n";
										$cmd .= "\t\tif test -e \"$tmp_dbprot.md5\"; then\n";
										$cmd .= "\t\t\trm $tmp_dbprot.md5\n";
										$cmd .= "\t\tfi\n";
										$cmd .= "\t\tif test -e \"$tmp_dbprot.suf\"; then\n";
										$cmd .= "\t\t\trm $tmp_dbprot.\*\n";
										$cmd .= "\t\tfi\n";
										$cmd .= "\t\tif test -e \"$tmp_dbindex.bck\"; then\n";
										$cmd .= "\t\t\trm $tmp_dbindex.\*\n";
										$cmd .= "\t\tfi\n";
										$cmd .= "\t\tif test -e \"$tmp_qfasta\"; then\n";
										$cmd .= "\t\t\trm $tmp_qfasta\n";
										$cmd .= "\t\tfi\n";
										$cmd .= "\t\tif test -e \"$tmp_qfasta.fai\"; then\n";
										$cmd .= "\t\t\trm $tmp_qfasta.\*\n";
										$cmd .= "\t\tfi\n";
										$cmd .= "\t\tif test -e \"$tmp_qgff\"; then\n";
										$cmd .= "\t\t\trm $tmp_qgff\n";
										$cmd .= "\t\tfi\n";
										$cmd .= "\tfi\n";
										$cmd .= "fi\n";
									}
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
				if(-e $log_gth){
					system("rm $log_gth");
				}
				$cmd .= "echo 'all_done' >> $log_gth\n";
				SAVE($script, $cmd);
				push(@I1, $script);
				push(@I4, $log_gth);
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
			print "! genome threader with [$cpu] nodes...\n";
			my $num_I1 = @I1;
			for(my $i = 0; $i < $num_I1; $i++){
				system("qsub $I1[$i]");
			}
		}
		else{
			print "! genome threader with [$cpu] threads...\n";
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
		foreach my $log_gth (@I4){
			$ncompleted_gth += JudgeJob($log_gth);
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
	
	print "! summarizing results...\n";
	my $hralign = Txt2hash($tmp_ralign);
#		my $summary = "directory\tdbfasta\tseqid\tpos0\tpos1\tsp pos0\tsp pos1\tqfasta\tqfasta seqid\tpos0\tpos1\tgroup\thit sp0\thit sp1\tbp (sp0-sp1)\tbp hit (sp0-sp1)\tratio\tjudge\tgene_id\tprotein coding\tq protein length (predicted)\talign length\t\%similarity\t\%gap\tscore\tnum_NNN (qseq)\n";
	my $summary;
	
	foreach my $gid (@G){
		my $len_NNN = 0;
		if($hralign->{$gid}{cnt_NNN}){
			$len_NNN = $hralign->{$gid}{cnt_NNN};
		}
		
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
			if($hralign->{$gid}{sub_qprot} && -e $hralign->{$gid}{sub_qprot}){
				$len_qprot = Count_seqlen($hralign->{$gid}{sub_qprot});
			}
			
			if($hralign->{$gid}{align} && -e $hralign->{$gid}{align}){
				my $htmp = Read_matcher($hralign->{$gid}{align});
				$summary .= $hpav->{$gid}{data}."\t".$hlist->{$gid}{protein_coding}."\t".$len_dbprot."\t".$len_qprot."\t".$htmp->{len}."\t".$htmp->{similarity}."\t".$htmp->{gap}."\t".$htmp->{score}."\n";
			}
			else{
				$summary .= $hpav->{$gid}{data}."\t".$hlist->{$gid}{protein_coding}."\t".$len_dbprot."\t".$len_qprot."\t.\t.\t.\t.\n";
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
	my $subseq = substr($seq, $p0, $p1);
	my $subfasta = ">".$id."\n".$subseq;
	
	unless(-e $rfile){
		open(my $rfh, ">", $rfile);
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
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	if($line =~ /directory/ && $line =~ /dbfasta/ && $line =~ /seqid/){
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
}
close $fh;

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
		
		if($pav_status && $pav_status =~ /present/i){
			$hash->{$gid}{status} = "P";
			$hash->{$gid}{dbsid} = $A[2];
			$hash->{$gid}{dbpos0} = $A[5];
			$hash->{$gid}{dbpos1} = $A[6];
			$hash->{$gid}{qsid} = $A[8];
			$hash->{$gid}{qpos0} = $A[12];
			$hash->{$gid}{qpos1} = $A[13];
			$hash->{$gid}{data} = $line;
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
			$hash->{$gid}{data} = $line;
			$cntP++;
		}
		else{
			$hash->{$gid}{status} = "A";
			$hash->{$gid}{data} = $line;
			$cntA++;
		}
		
		$cnt++;
	}
}
close $fh;

print "! total [$cnt] lines, P=[$cntP], A=[$cntA]\n";

return $hash;
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
	
	if($A[0] =~ /scaffold/ || $A[0] =~ /unanchored/ || $A[0] =~ /mitochon/ || $A[0] =~ /chloro/){
		next;
	}
	
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
		die;
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
		print "! thread [$j] : failed\n";
		die;
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

$script_path =~ s/\/scripts//;
my @tmp = split(/\//, $genelist);
my $ntmp = @tmp;
$genelist = $tmp[$ntmp - 1];

my $tmp = "_Search_app_temp.$genelist.txt";
for(my $i = 0; $i < 100; $i++){
	if(-e $tmp){
		$tmp = "_Search_app_temp.$genelist.$i.txt";
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



