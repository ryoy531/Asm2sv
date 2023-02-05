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
$version .= "batch_psl2sv_clsr.pl version 1.01\n";
$version .= "last update: [2020\/11\/22]\n";
$version .= "copyright: ryoichi yano [ryoichiy104\@gmail.com]\n";

#print $version;

#-------------------------------------------------------------------------------

#---------------------------------------------------------//
my $help_command =<<"EOS";
Basic usage: 
  /path/to/batch_psl2sv_clsr.pl -d [reference fasta] -g [reference gff3] -l [gene list] -q [query fasta] -a [psl] -c [chralias] -s [init span] -n [expand neighbor length] -o [Result_directory] --thread [CPU thread]

--db or -d            : reference sequence (fasta, required)
--gff or -g           : reference gff (gff, required)
--list or -l          : list of gene to be analyzed (text, required)
--query or -q         : query fasta (fasta, required)
--alignpsl or -a      : LAST alignment (psl, required)
--chralias or -c      : chromosome ID alias info (tsv, required)
--span or -s          : length of flanking region to be assessed
--neighbor or -n      : include [X] bp neighboring seq in addition to target (e.g. promoter)
--wo_orfsearch or -w  : without ORF prediction
--thread or -t        : CPU thread number (default 8)
--output_dir or -o    : output directory
--plot                : output alignment plot then stop
--help or -h          : display usage

EOS

#---------------------------------------------------------//
my $host = hostname();
my $script_path = $FindBin::Bin;
my $wpath = getcwd();
#my $HOME = File::HomeDir->my_home;
my $bin_findorf = $script_path."/findorf.pl";
my $bin_prepseq4annot = $script_path."/prepseq4annot.pl";
my $bin_FindPslPos_clsr = $script_path."/FindPslPos_clsr.pl";
my $ramsize = CheckRAM();		#Mb scale

#---------------------------------------------------------//
#gettin parameters from command line
my $dir;
my $dbfasta;
my $dbgff;
my $genelist;
my $qfasta;
my $qpsl;
my $chrinfo;
my $outplot;
my $neighbor_tlen;
my $cpu;
my $cpu_gth;
my $initspan;
my $skip_orfsearch;
my $help;
my $test;

GetOptions('--db=s' => \$dbfasta, '--gff=s' => \$dbgff, '--list=s' => \$genelist, '--query=s' => \$qfasta, '--alignpsl=s' => \$qpsl, '--chralias=s' => \$chrinfo, '--output_dir=s' => \$dir, '--neighbor=i' => \$neighbor_tlen, '--thread=i' => \$cpu, '--xthread=i' => \$cpu_gth, '--span=i' => \$initspan, '--plot' => \$outplot, '--wo_orfsearch' => \$skip_orfsearch, '--help' => \$help);

$host = "cluster1";

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
if(! $qpsl || ! -e $qpsl){
	print "! missing last alignment psl file\n";
	$status = "missing";
}
if(! $chrinfo || ! -e $chrinfo){
	$chrinfo = "null";
}
if($status eq 'missing'){
	print "\n$help_command";
	goto END;
}
if($outplot){
	print "\n! --plot is invoked\n";
}

if(! $initspan || $initspan =~ /\D/){
	$initspan = 50000;
}
my $kbflnk = int($initspan / 1000);

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
elsif($cpu > 128){
	$cpu = 128;
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
$log .= "! chr alias info      : [$chrinfo]\n";
$log .= "! query fasta         : [$qfasta]\n";
$log .= "! query alignment psl : [$qpsl]\n";
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
my $refgenome = $dbfasta;
if($refgenome =~ /\.fasta/){
	$refgenome =~ s/\.fasta//g;
}
elsif($refgenome =~ /\.fa/){
	$refgenome =~ s/\.fa//g;
}

my $qpref = $qfasta;
if($qpref =~ /\.fasta/){
	$qpref =~ s/\.fasta//g;
}
elsif($qpref =~ /\.fa/){
	$qpref =~ s/\.fa//g;
}

unless(-e $dir){
	system("mkdir $dir");
}
chdir $dir;

unless(-e "cmds"){
	system("mkdir cmds");
}
unless(-e "log"){
	system("mkdir log");
}

$dbfasta = LnFile($dbfasta);
$dbgff = LnFile($dbgff);
$genelist = LnFile($genelist);
$qfasta = LnFile($qfasta);
$qpsl = LnFile($qpsl);
$chrinfo = LnFile($chrinfo);

my $header_summary = "directory\tdbfasta\tseqid\tpos0\tpos1\tsp0\tsp1\tqfasta\tqfasta seqid\thit pos0\thit pos1\tgroup\thit sp0\thit sp1\tbp (sp0-sp1)\tbp hit-align (sp0-sp1)\tbp hit-span (sp0-sp1)\tbp hit-span woN (sp0-sp1)\tbp insertion\talign ratio (1=not disrupted)\tinsert ratio (1=not disrupted)\tseq normality (1=not disrupted)\tjudge\tgene_id\t";
$header_summary .= "5'-promoter bp hit ($neighbor_tlen bp)\t5'-promoter align ratio ($neighbor_tlen bp)\t3'-UTR bp hit ($neighbor_tlen bp)\t3'-UTR align ratio ($neighbor_tlen bp)\t";
$header_summary .= "5' flanking bp hit ($initspan bp)\t5' flanking align ratio ($initspan bp)\t3' flanking bp hit ($initspan bp)\t3' flanking align ratio ($initspan bp)";

my $combined_final = "summary_genome2sv_".$dir."_combined.tsv";
my $rfile_final = "summary_genome2sv_".$dir."_".$qpref.".tsv";
my $afile = "specified_region_hits_".$qpref.".tsv";
my $ffile = "failed_ID_list_".$qpref.".tsv";

print "! preparing info...\n";
my $hID = Read_list($genelist);
my $hgffinfo = Read_gff($dbgff, $hID);
my $QL = $hgffinfo->{RL};
my $hlist = $hgffinfo->{hash};
my $dbg2t = $hgffinfo->{g2t};
my $dbt2g = $hgffinfo->{t2g};
my $numQL = @{$QL};
my $hqseq = Open_fasta_as_hash($qfasta);

#---------------------------------------------------------------//
my $sw = 1;
if($outplot){
	if(-e $afile){
		$sw = 2;
	}
}
else{
	if(-e $afile){
		$sw = 0;
	}
}

if($sw > 0){
	print "! step.1 analyzing presence-absence in [$qfasta]...\n";
	my $ho = {};
	$ho->{hsample}{0} = $refgenome;
	$ho->{hsample}{1} = $qpref;
	$ho->{name}{$refgenome} = 1;
	$ho->{name}{$qpref} = 1;
	
	my @ALNF;
	$ALNF[0] = $qpsl;
	
	my $np = int($numQL / $cpu) + 1;
	
	my $hseqinfo = {};
	if($chrinfo && $chrinfo ne 'null' && -e $chrinfo){
		$hseqinfo = Read_seqinfo($chrinfo);
	}
	
	my $hits_header = "gene_id\tbase_ref\tchr\tpos0\tpos1\tflanking\tclass\tsearch_mode\tnmatches (bp)\taln_strand\thit_pos0\thit_pos1\tq_name\tq_seqID\tq_seqlen\tq_apos0\tq_apos1\tdb_name\tdb_seqID\tdb_seqlen\tdb_apos0\tdb_apos1\tblockCount\tblockSizes\tqStarts\tdbStarts\tlimit_qpos0\tlimit_qpos1\n";
	my $summary_header = "gene_id\tbase_ref\tchr\tpos0\tpos1\tflanking\tclass\tqname\thitbp (ref-base)\thitbp (+flanking, ref-base)\tq-hitseq\thit_pos0\thit_pos1\tlimit_qpos0\tlimit_qpos1\tnum_candidate_chrs\n";
	my $summary2_header = $header_summary."\n";

	my $rfile_hits0 = "r1_hitalign_flk-".$kbflnk."kb.tsv";
	my $rfile_hits1 = "r1_hitsummary_flk-".$kbflnk."kb.tsv";
	my $outplot_dir = "plot";
	if($outplot && ! -e $outplot_dir){
		system("mkdir $outplot_dir");
	}
	if(! $outplot){
		$outplot = 0;
	}
	if(! $chrinfo || ! -e $chrinfo){
		$chrinfo = "null";
	}
	
	print "! analyze genomic alignment with [$cpu] nodes...\n";
	my @SH;
	my @I1;
	my @I2;
	for(my $i = 0; $i < $cpu; $i++){
		my $p0 = $i * $np;
		my $p1 = ($i + 1) * $np;
		
		if($p0 > $numQL){
			next;
		}
		
		my $log_fpsl = "log_batch_findpslpos_$i.txt";
		my $cmd = "perl $bin_FindPslPos_clsr $refgenome $dbgff $qpref $qfasta $qpsl $genelist $chrinfo $afile $ffile $rfile_hits0 $rfile_hits1 $outplot_dir $dir $wpath $initspan $neighbor_tlen $p0 $p1 $np $outplot\n";
		$cmd .= "echo 'all_done' >> $log_fpsl\n";
		my $sh = "cmd_findpslpos_$qpref"."_$i.sh";
		SAVE($sh, $cmd, "");
		push(@I1, $sh);
		push(@I2, $log_fpsl);
		system("qsub $sh");
	}
	
	print "! waiting for job completion...\n";
	LOOPjob:{
		my $loopjob = 1;
	}
	my $numI1 = @I1;
	my $ncompleted_pav = 0;
	for(my $i = 0; $i < $numI1; $i++){
		$ncompleted_pav += JudgeJob($I2[$i]);
	}
	if($ncompleted_pav < $numI1){
		sleep(30);
		goto LOOPjob;
	}
	else{
		for(my $i = 0; $i < $numI1; $i++){
			if(-e $I1[$i]){
				my $mvcmd1 = "mv $I1[$i]"."* cmds";
				my $mvcmd2 = "mv $I2[$i]"."* log";
				system($mvcmd1);
				system($mvcmd2);
			}
		}
	}
}
else{
	print "! [$afile] already exits, skip 1st search...\n";
}
#---------------------------------------------------------------//
my $header_summary2 = "directory\tdbfasta\tseqid\tpos0\tpos1\tsp0\tsp1\tqfasta\tqfasta seqid\thit pos0\thit pos1\tgroup\thit sp0\thit sp1\tbp (sp0-sp1)\tbp hit-align (sp0-sp1)\tbp hit-span (sp0-sp1)\tbp hit-span woN (sp0-sp1)\tbp insertion\talign ratio (1=not disrupted)\tinsert ratio (1=not disrupted)\tseq normality (1=not disrupted)\tjudge\tgene_id\t";
$header_summary2 .= "5'-promoter bp hit ($neighbor_tlen bp)\t5'-promoter align ratio ($neighbor_tlen bp)\t3'-UTR bp hit ($neighbor_tlen bp)\t3'-UTR align ratio ($neighbor_tlen bp)\t";
$header_summary2 .= "5' flanking bp hit ($initspan bp)\t5' flanking align ratio ($initspan bp)\t3' flanking bp hit ($initspan bp)\t3' flanking align ratio ($initspan bp)\t";
$header_summary2 .= "protein coding\tdb protein length\tq protein length (predicted)\talign length\t\%similarity\t\%gap\tscore\n";

if($skip_orfsearch){
	print "! --wo_orfsearch is invoked, skip ORF prediction...\n";
	print "! preparing final result...\n";
	my $dbprotein = $refgenome."_protein.fasta";
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
				$summary = $header_summary2.$summary;
				SAVE($combined_final, $summary);
				print "! output [$combined_final]\n";
			}
		}
		else{
			$summary = $header_summary2.$summary;
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
	my $hdbseq = Open_fasta_as_hash($dbfasta);
	my $hqseq = Open_fasta_as_hash($qfasta);
	
	my @G = keys(%{$hlist});
	@G = sort {$a cmp $b} @G;
	
	my $dbfai = $dbfasta.".fai";
	my $dbtranscript = $refgenome."_transcript.fasta";
	my $dbcds = $refgenome."_cds.fasta";
	my $dbprotein = $refgenome."_protein.fasta";
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
			system("touch $rfile_final");
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
				$summary = $header_summary2.$summary;
				SAVE($combined_final, $summary);
				print "! output [$combined_final]\n";
			}
		}
		else{
			$summary = $header_summary2.$summary;
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
	system("touch $rfile_final");
	goto END;
}

SORT:{
	Sort_and_save($rfile_final);
}

END:{
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


#---------------------------------------------------------------------
sub Read_seqinfo{
my $file = shift;

print "! reading [$file]...\n";
open(my $fh, "<", $file);
my $hash = {};
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	$line =~ s/\"//g;
	
	my @A = split(/\t/, $line);
	my @B = split(/,/, $A[2]);
	$hash->{$A[0]}{$A[2]} = $A[1];
	$hash->{$A[0]}{$B[0]} = $A[1];
#	print "$A[0] $B[0] $A[1]\n";
	$cnt++;
}
close $fh;

print "! [$cnt] lines\n";

return $hash;
}


#---------------------------------------------------------------------
sub FindPslPos_rev{
my $ALNF = shift;
my $QL = shift;
my $p0 = shift;
my $p1 = shift;
my $np = shift;
my $initspan = shift;
my $neighbor_tlen = shift;
my $hseqinfo = shift;
my $ho = shift;
my $t = shift;
my $refgenome = shift;
my $hqseq = shift;
my $outplot = shift;
my $outplot_dir = shift;
my $dir = shift;
my $wpath = shift;
my $qpref = shift;

my $hsample = $ho->{hsample};

print " partition [$t] : search for [$np] entries...\n";
my $hash = {};
foreach my $file (@{$ALNF}){
	my @F0 = split(/q-/, $file);
	my @F1 = split(/_DB-/, $F0[1]);
	my $qname0 = $F1[0];
	my $dbname0 = $F1[1];
	$dbname0 =~ s/\.psl//;
	
	if($ho->{altname}{$qname0}){
		$qname0 = $ho->{altname}{$qname0};
	}
	if($ho->{altname}{$dbname0}){
		$dbname0 = $ho->{altname}{$dbname0};
	}
	
	if($ho->{name}{$qname0} && $ho->{name}{$dbname0}){
		open(my $fh, "<", $file);
		while(my $line = <$fh>){
			$line =~ s/\n//;
			$line =~ s/\r//;
			$hash->{$file} .= $line."\n";
		}
		close $fh;
	}
}

my $hit = "";
my $summary = "";
my $summary2 = "";
my $failed_list = "";
my $cntq = 0;
my $cnt_hit = 0;
my $cnt_selected = 0;
my $seg = 100;
my $add = $seg;
for(my $p = $p0; $p < $p1; $p++){
	if(! $QL->[$p]){
		next;
	}
	
	my @Q = split(/,/, $QL->[$p]);
	my @Ord = keys(%{$hsample});
	@Ord = sort {$a <=> $b} @Ord;
	my $numOrd = @Ord;
	
	my $hhitinfo = {};
	foreach my $norder (@Ord){
		if($norder == 0){
			next;
		}
		my $j = $norder;
		my $i = 0;
		
		my $psl = "";
		my $qname = "";
		my $dbname = "";
		my $mode = "false";
		
		if($hsample->{$i} && $hsample->{$j}){
			foreach my $file (@{$ALNF}){
				my @F0 = split(/q-/, $file);
				my @F1 = split(/_DB-/, $F0[1]);
				my $qname0 = $F1[0];
				my $dbname0 = $F1[1];
				$dbname0 =~ s/\.psl//;
				
				if($dbname0 eq $hsample->{$i} && $qname0 eq $hsample->{$j}){
					$psl = $file;
					$qname = $hsample->{$j};
					$dbname = $hsample->{$i};
					$mode = "fw";
					last;
				}
				else{
					if($ho->{altname}{$qname0}){
						$qname0 = $ho->{altname}{$qname0};
					}
					if($ho->{altname}{$dbname0}){
						$dbname0 = $ho->{altname}{$dbname0};
					}
					
					if($dbname0 eq $hsample->{$i} && $qname0 eq $hsample->{$j}){
						$psl = $file;
						$qname = $hsample->{$j};
						$dbname = $hsample->{$i};
						$mode = "fw";
						last;
					}
				}
			}
		}
		
		if($mode eq 'false'){
			print "! error : incompatible file name...\n";
			die;
		}
		
	#	print "! searching for [$psl]...\n";
		my @DL = split(/\n/, $hash->{$psl});
		my $sw = 0;
		foreach my $line (@DL){
			$line =~ s/\n//;
			$line =~ s/\r//;
			
	#		1. matches - Number of matching bases that aren't repeats.
	#		2. misMatches - Number of bases that don't match.
	#		3. repMatches - Number of matching bases that are part of repeats.
	#		4. nCount - Number of 'N' bases.
	#		5. qNumInsert - Number of inserts in query.
	#		6. qBaseInsert - Number of bases inserted into query.
	#		7. tNumInsert - Number of inserts in target.
	#		8. tBaseInsert - Number of bases inserted into target.
	#		9. strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
	#		10. qName - Query sequence name.
	#		11. qSize - Query sequence size.
	#		12. qStart - Alignment start position in query.
	#		13. qEnd - Alignment end position in query.
	#		14. tName - Target sequence name.
	#		15. tSize - Target sequence size.
	#		16. tStart - Alignment start position in target.
	#		17. tEnd - Alignment end position in target.
	#		18. blockCount - Number of blocks in the alignment.
	#		19. blockSizes - Comma-separated list of sizes of each block.
	#		20. qStarts - Comma-separated list of start position of each block in query.
	#		21. tStarts - Comma-separated list of start position of each block in target.
			
			my @L = split(/\t/, $line);
			if($L[9] =~ /scaffold/ || $L[13] =~ /scaffold/){
				next;
			}
			
			if($mode eq 'fw'){
				my $seqid9 = $L[9];
				my $seqid13 = $L[13];
				
				if($hseqinfo->{$qname}{$L[9]}){
					$seqid9 = $hseqinfo->{$qname}{$L[9]}." | ".$L[9];
#					$L[9] = $hseqinfo->{$qname}{$L[9]};
				}
				if($hseqinfo->{$dbname}{$L[13]}){
					$seqid13 = $hseqinfo->{$dbname}{$L[13]}." | ".$L[13];
#					$L[13] = $hseqinfo->{$dbname}{$L[13]};
				}
				
				my $qseqName = $seqid9;
				my $qSize = $L[10];
				my $qStart = $L[11];
				my $qEnd = $L[12];
				my $tseqName = $seqid13;
				my $tSize = $L[14];
				my $tStart = $L[15];
				my $tEnd = $L[16];
				my $blockCount = $L[17];
				my $blockSizes = $L[18];
				my $qStarts = $L[19];
				my $tStarts = $L[20];
				
#				if($L[9] eq $Q[1] && $L[13] eq $Q[1]){
				if($L[13] eq $Q[1]){		# only when chr ID match
					for(my $p = $Q[2] - $initspan; $p < $Q[3] + $initspan; $p += 1){
						if($L[15] <= $p && $p <= $L[16]){
							my $gene_region1 = "F";
							my $gene_region2 = "F";
							my @BLK = split(/,/, $blockSizes);
							my @QS = split(/,/, $qStarts);
							my @TS = split(/,/, $tStarts);
							my $numBLK = @BLK;
							
							if(! $hhitinfo->{$j}{$L[9]}{pos0}){
								for(my $tmp_q2 = $Q[2]; $tmp_q2 < $Q[3]; $tmp_q2++){
									if($hhitinfo->{$j}{$L[9]}{pos0}){
										last;
									}
									
									if($L[15] <= $tmp_q2 && $tmp_q2 <= $L[16]){
										for(my $k = 0; $k < $numBLK; $k++){
											if($L[8] eq '+'){
												my $blk_qstart = $QS[$k];
												my $blk_qend = $QS[$k] + $BLK[$k];
												my $blk_tstart = $TS[$k];
												my $blk_tend = $TS[$k] + $BLK[$k];
												
												if($TS[$k] <= $tmp_q2 && $tmp_q2 <= $blk_tend){
													my $coef = ($blk_qend - $blk_qstart) / ($blk_tend - $blk_tstart);
													my $intc = $blk_qstart - $coef * $blk_tstart;
													my $hitpos = $coef * $tmp_q2 + $intc;
													$gene_region1 = $hitpos." | $L[8] | ($blk_qend - $blk_qstart) / ($blk_tend - $blk_tstart) | $Q[2] | $tmp_q2";
													$hhitinfo->{$j}{$L[9]}{pos0} = $hitpos;
													last;
												}
											}
											elsif($L[8] eq '-'){
												my $blk_qstart = $qSize - $QS[$k] - $BLK[$k];
												my $blk_qend = $qSize - $QS[$k];
												my $blk_tstart = $TS[$k];
												my $blk_tend = $TS[$k] + $BLK[$k];
												
												if($TS[$k] <= $tmp_q2 && $tmp_q2 <= $blk_tend){
													my $coef = ($blk_qstart - $blk_qend) / ($blk_tend - $blk_tstart);
													my $intc = $blk_qstart - $coef * $blk_tstart;
													my $hitpos = $coef * $tmp_q2 + $intc;
													$gene_region1 = $hitpos." | $L[8] | ($blk_qstart - $blk_qend) / ($blk_tend - $blk_tstart) | $Q[2] | $tmp_q2";
													$hhitinfo->{$j}{$L[9]}{pos0} = $hitpos;
													last;
												}
											}
										}
									}
								}
							}
							if(! $hhitinfo->{$j}{$L[9]}{pos1}){
								for(my $tmp_q3 = $Q[3]; $tmp_q3 > $Q[2]; $tmp_q3--){
									if($hhitinfo->{$j}{$L[9]}{pos1}){
										last;
									}
									
									if($L[15] <= $tmp_q3 && $tmp_q3 <= $L[16]){
										for(my $k = 0; $k < $numBLK; $k++){
											if($L[8] eq '+'){
												my $blk_qstart = $QS[$k];
												my $blk_qend = $QS[$k] + $BLK[$k];
												my $blk_tstart = $TS[$k];
												my $blk_tend = $TS[$k] + $BLK[$k];
												
												if($TS[$k] <= $tmp_q3 && $tmp_q3 <= $blk_tend){
													my $coef = ($blk_qend - $blk_qstart) / ($blk_tend - $blk_tstart);
													my $intc = $blk_qstart - $coef * $blk_tstart;
													my $hitpos = $coef * $tmp_q3 + $intc;
													$gene_region2 = $hitpos." | $L[8] | ($blk_qend - $blk_qstart) / ($blk_tend - $blk_tstart) | $Q[3] | $tmp_q3";
													$hhitinfo->{$j}{$L[9]}{pos1} = $hitpos;
													last;
												}
											}
											elsif($L[8] eq '-'){
												my $blk_qstart = $qSize - $QS[$k] - $BLK[$k];
												my $blk_qend = $qSize - $QS[$k];
												my $blk_tstart = $TS[$k];
												my $blk_tend = $TS[$k] + $BLK[$k];
												
												if($TS[$k] <= $tmp_q3 && $tmp_q3 <= $blk_tend){
													my $coef = ($blk_qstart - $blk_qend) / ($blk_tend - $blk_tstart);
													my $intc = $blk_qstart - $coef * $blk_tstart;
													my $hitpos = $coef * $tmp_q3 + $intc;
													$gene_region2 = $hitpos." | $L[8] | ($blk_qstart - $blk_qend) / ($blk_tend - $blk_tstart) | $Q[3] | $tmp_q3";
													$hhitinfo->{$j}{$L[9]}{pos1} = $hitpos;
													last;
												}
											}
										}
									}
								}
							}
							
							$hhitinfo->{$j}{$L[9]}{eachseq_nmatch} += $L[0];
							$hhitinfo->{$j}{$L[9]}{allpos} .= $qStart.",".$qEnd.",";
							$hhitinfo->{$j}{$L[9]}{data} .= $Q[0]."\t".$refgenome."\t".$Q[1]."\t".$Q[2]."\t".$Q[3]."\t".$initspan."\tref\t".$mode."\t".$L[0]."\t".$L[8]."\t".$gene_region1."\t".$gene_region2."\t";
							$hhitinfo->{$j}{$L[9]}{data} .= $qname."\t".$qseqName."\t".$qSize."\t".$qStart."\t".$qEnd."\t";
							$hhitinfo->{$j}{$L[9]}{data} .= $dbname."\t".$tseqName."\t".$tSize."\t".$tStart."\t".$tEnd."\t";
							$hhitinfo->{$j}{$L[9]}{data} .= $blockCount."\t".$blockSizes."\t".$qStarts."\t".$tStarts."\n";
							$cnt_hit++;
							last;
						}
					}
				}
			}
		}
	}
	foreach my $norder (@Ord){
		if($norder == 0){
			next;
		}
		my $j = $norder;
		
		if($hhitinfo->{$j}){
			my $htmp = $hhitinfo->{$j};
			my @HCRS = keys(%{$htmp});
			@HCRS = sort {$a cmp $b} @HCRS;
			my $num_HCRS = @HCRS;
			
			my $AoHCRS = [];
			foreach my $htseq (@HCRS){
				unless($hhitinfo->{$j}{$htseq}{eachseq_nmatch}){
					$hhitinfo->{$j}{$htseq}{eachseq_nmatch} = 0;
				}
				
				my @tmp;
				push(@tmp, $htseq);
				push(@tmp, $hhitinfo->{$j}{$htseq}{eachseq_nmatch});
				push(@{$AoHCRS}, \@tmp);
			}
			
			@{$AoHCRS} = sort {$b->[1] <=> $a->[1]} @{$AoHCRS};
			
			my $tophseq = $AoHCRS->[0][0];
			
			unless($hhitinfo->{$j}{$tophseq}{pos0}){
				$hhitinfo->{$j}{$tophseq}{pos0} = "null";
			}
			unless($hhitinfo->{$j}{$tophseq}{pos1}){
				$hhitinfo->{$j}{$tophseq}{pos1} = "null";
			}
			
			if($hhitinfo->{$j}{$tophseq}{allpos} && $hhitinfo->{$j}{$tophseq}{allpos} ne 'null'){
				my @AllPos0 = split(/,/, $hhitinfo->{$j}{$tophseq}{allpos});
				my @AllPos;
				my $cnt_pos = 0;
				foreach my $posx (@AllPos0){
					push(@AllPos, $posx);
					$cnt_pos++;
				}
				@AllPos = sort {$a <=> $b} @AllPos;
				my $min_pos = $AllPos[0];
				my $max_pos = $AllPos[$cnt_pos - 1];
				my $n50_pos = $AllPos[int($cnt_pos / 2)];
				my $limit_pos0 = $n50_pos - ($initspan * 1.2);
				my $limit_pos1 = $n50_pos + ($initspan * 1.2);
				
				if($hhitinfo->{$j}{$tophseq}{pos0} && $hhitinfo->{$j}{$tophseq}{pos0} ne 'null'){
					$limit_pos0 = $hhitinfo->{$j}{$tophseq}{pos0} - ($initspan * 1);
					if($limit_pos0 < 0){
						$limit_pos0 = 0;
					}
				}
				if($hhitinfo->{$j}{$tophseq}{pos1} && $hhitinfo->{$j}{$tophseq}{pos1} ne 'null'){
					$limit_pos1 = $hhitinfo->{$j}{$tophseq}{pos1} + ($initspan * 1);
				}
				
				my @HL = split(/\n/, $hhitinfo->{$j}{$tophseq}{data});
				my $AoHL0 = [];
				foreach my $line (@HL){
					my @A = split(/\t/, $line);
					push(@{$AoHL0}, \@A);
				}
				@{$AoHL0} = sort {$a->[15] <=> $b->[15]} @{$AoHL0};		# sort by qstart
				
				my $blk_qspan_th = 100 * 1000;		# if distance is > 100 kb, seprate alignment blocks within the tophit chromosome
				my $prev_qpos1 = 0;
				my $hqblk = {};
				my $num_qblk = 0;
				foreach my $A (@{$AoHL0}){
					if( ($A->[15] - $prev_qpos1) <= $blk_qspan_th){
						$hqblk->{$num_qblk}{data} .= join("\t", @{$A})."\n";
						$hqblk->{$num_qblk}{nmatch} += $A->[8];
					}
					else{
						$num_qblk++;
						$hqblk->{$num_qblk}{data} .= join("\t", @{$A})."\n";
						$hqblk->{$num_qblk}{nmatch} += $A->[8];
					}
					$prev_qpos1 = $A->[16];
				}
				
				my @NQBLK = keys(%{$hqblk});
				@NQBLK = sort {$a <=> $b} @NQBLK;
				
				my $Aoqtmp = [];
				foreach my $nqblk (@NQBLK){
					my @tmp;
					push(@tmp, $hqblk->{$nqblk}{data});
					push(@tmp, $hqblk->{$nqblk}{nmatch});
					push(@{$Aoqtmp}, \@tmp);
				}
				@{$Aoqtmp} = sort {$b->[1] <=> $a->[1]} @{$Aoqtmp};
				
				@HL = split(/\n/, $Aoqtmp->[0][0]);			# select alignment block with top score within the  tophit chromosome
				
				my $AoHL = [];
				my $sum_nmatch_re = 0;
				my @AllPosS;
				foreach my $line (@HL){
					my @A = split(/\t/, $line);
	#				if($limit_pos0 <= $A[15] && $A[15] <= $limit_pos1 && $limit_pos0 <= $A[16] && $A[16] <= $limit_pos1){
						push(@{$AoHL}, \@A);
						$sum_nmatch_re += $A[8];
						push(@AllPosS, $A[15]);
						push(@AllPosS, $A[16]);
						$cnt_selected++;
	#				}
				}
				
				my $cnt_posre = @AllPosS;
				@AllPosS = sort {$a <=> $b} @AllPosS;
				my $min_pos_re = $AllPosS[0];
				my $max_pos_re = $AllPosS[$cnt_posre - 1];
				@{$AoHL} = sort {$a->[20] <=> $b->[20]} @{$AoHL};		# sort by tstart
				
				my $gpos0 = $Q[2];
				my $gpos1 = $Q[3];
				my $gpos0f = $Q[2] - $initspan;
				my $gpos1f = $Q[3] + $initspan;
				my $gpos0n1 = $Q[2] - $neighbor_tlen;
				my $gpos1n1 = $Q[2];
				my $gpos0n2 = $Q[3];
				my $gpos1n2 = $Q[3] + $neighbor_tlen;
				my $gpos0n3 = $Q[2] - $initspan;
				my $gpos1n3 = $Q[2];
				my $gpos0n4 = $Q[3];
				my $gpos1n4 = $Q[3] + $initspan;
				my $gphitbp0 = 0;
				my $gphitbp1 = 0;
				my $gphitbp2 = 0;
				my $gphitbp3 = 0;
				my $gphitbp4 = 0;
				my $gphitbp5 = 0;
				my @Gbp0;
				my @Gbp1;
				my @Gbp2;
				my @Gbp3;
				my @Gbp4;
				my @Gbp5;
				my @Qbp0;
				my @Qbp1;
				my @Qbp2;
				my @Qbp3;
				my @Qbp4;
				my @Qbp5;
				my $dotdata_flanking_fw = "qpos,refpos\n";
				my $dotdata_flanking_rv = "qpos,refpos\n";
				
				foreach my $A (@{$AoHL}){
					$hit .= join("\t", @{$A})."\t".$limit_pos0."\t".$limit_pos1."\n";
					
					# 4317,501,	49877752,49882073 (89229,84908),	77485,81802,
					
					my $strand = $A->[9];
					my @BLK = split(/,/, $A->[23]);
					my @QS = split(/,/, $A->[24]);
					my @GS = split(/,/, $A->[25]);		# TS
					my $numBLK = @BLK;
					
					for(my $k = 0; $k < $numBLK; $k++){
						my $gstart = $GS[$k];
						my $gend = $GS[$k] + $BLK[$k];
						my $tmp_gbp = $GS[$k];
						my $tmp_qbp = $QS[$k];
						
						if($A->[9] eq '-'){
							$tmp_qbp = $A->[14] - $tmp_qbp;		# 49966981 - 49877752 = 89229
						}
						
						for(my $gp = $gstart; $gp < $gend; $gp++){
							if($gpos0 <= $gp && $gp <= $gpos1){
								$gphitbp0++;
								push(@Gbp0, $tmp_gbp);
								push(@Qbp0, $tmp_qbp);
							}
							if($gpos0f <= $gp && $gp <= $gpos1f){
								$gphitbp1++;
								push(@Gbp1, $tmp_gbp);
								push(@Qbp1, $tmp_qbp);
								
								if($A->[9] eq '+'){
									$dotdata_flanking_fw .= "$tmp_qbp,$tmp_gbp\n";
								}
								elsif($A->[9] eq '-'){
									$dotdata_flanking_rv .= "$tmp_qbp,$tmp_gbp\n";
								}
							}
							if($gpos0n1 <= $gp && $gp < $gpos1n1){
								$gphitbp2++;
								push(@Gbp2, $tmp_gbp);
								push(@Qbp2, $tmp_qbp);
							}
							if($gpos0n2 < $gp && $gp <= $gpos1n2){
								$gphitbp3++;
								push(@Gbp3, $tmp_gbp);
								push(@Qbp3, $tmp_qbp);
							}
							if($gpos0n3 <= $gp && $gp < $gpos1n3){
								$gphitbp4++;
								push(@Gbp4, $tmp_gbp);
								push(@Qbp4, $tmp_qbp);
							}
							if($gpos0n4 < $gp && $gp <= $gpos1n4){
								$gphitbp5++;
								push(@Gbp5, $tmp_gbp);
								push(@Qbp5, $tmp_qbp);
							}
							
							$tmp_gbp++;
							if($A->[9] eq '+'){
								$tmp_qbp++;
							}
							elsif($A->[9] eq '-'){
								$tmp_qbp--;
							}
						}
					}
				}
				
				@Gbp0 = sort {$a <=> $b} @Gbp0;
				@Gbp1 = sort {$a <=> $b} @Gbp1;
				@Gbp2 = sort {$a <=> $b} @Gbp2;
				@Gbp3 = sort {$a <=> $b} @Gbp3;
				@Gbp4 = sort {$a <=> $b} @Gbp4;
				@Gbp5 = sort {$a <=> $b} @Gbp5;
				@Qbp0 = sort {$a <=> $b} @Qbp0;
				@Qbp1 = sort {$a <=> $b} @Qbp1;
				@Qbp2 = sort {$a <=> $b} @Qbp2;
				@Qbp3 = sort {$a <=> $b} @Qbp3;
				@Qbp4 = sort {$a <=> $b} @Qbp4;
				@Qbp5 = sort {$a <=> $b} @Qbp5;
				
				if(! $Qbp1[0] || ! $Qbp1[$gphitbp1 - 1]){
					if($outplot){
						print " partition [$t] : cannot draw plot for [$Q[0]] because of no alignment...\n";
						$failed_list .= $Q[0]."\tmissing Qbp1 (cannot draw plot)\n";
					}
					else{
						$failed_list .= $Q[0]."\tmissing Qbp1 (cannot draw plot)\n";
					}
					
					my $gpos0f = $Q[2] - $initspan;
					my $gpos1f = $Q[3] + $initspan;
					$summary2 .= $hsample->{$j}."\t".$refgenome."\t".$Q[1]."\t".$gpos0f."\t".$gpos1f."\t".$Q[2]."\t".$Q[3]."\t".$hsample->{$j}."\t-\t-\t-\t-\t";
					$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tabsent\t".$Q[0]."\t";
					$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\n";
					next;
				}
				
				#-----------------------------------------------------------------------//
				if($outplot){
					my $png1 = $Q[0]."_".$qpref."_vs_".$refgenome."_plot1.png";
					my $png2 = $Q[0]."_".$qpref."_vs_".$refgenome."_plot2.png";
					my $png3 = $Q[0]."_".$qpref."_vs_".$refgenome."_plot3.png";
					my $png4 = $Q[0]."_".$qpref."_vs_".$refgenome."_plot4.png";
					my $tmp_dotcsv_fw = $Q[0]."_".$qpref."_vs_".$refgenome."_fw.csv";
					my $tmp_dotcsv_rv = $Q[0]."_".$qpref."_vs_".$refgenome."_rv.csv";
					my $tmp_Rscript = "./$outplot_dir/".$Q[0]."_".$qpref."_vs_".$refgenome."_Rsc.txt";
					
					my $q2nb = $Q[2] - $neighbor_tlen;
					my $q3nb = $Q[3] + $neighbor_tlen;
					if($q2nb < 0){
						$q2nb = 0;
					}
					if($q3nb < 0){
						$q3nb = 0;
					}
					
					my $vline;
					$vline .= "abline(v\=$Q[2], col\=\"black\")";
					$vline .= "abline(v\=$Q[3], col\=\"black\")";
					$vline .= "abline(v\=$q2nb, col\=\"black\", lty=3)\n";
					$vline .= "abline(v\=$q3nb, col\=\"black\", lty=3)\n";
					my $hline;
					$hline .= "abline(h\=$Q[2], col\=\"black\")";
					$hline .= "abline(h\=$Q[3], col\=\"black\")";
					$hline .= "abline(h\=$q2nb, col\=\"black\", lty=3)\n";
					$hline .= "abline(h\=$q3nb, col\=\"black\", lty=3)\n";
					
					my $qpos_min = $Qbp1[0];
					my $qpos_max = $Qbp1[$gphitbp1 - 1];
					my $tpos_min = $Gbp1[0];
					my $tpos_max = $Gbp1[$gphitbp1 - 1];
					
					my $Rscript =<<"EOS";
workingDir = "$wpath/$dir/$outplot_dir"
setwd(workingDir)
getwd()

fw <- read.csv("$tmp_dotcsv_fw", header=T)
rv <- read.csv("$tmp_dotcsv_rv", header=T)

png("$png1", width=2400, height=2400)
plot(fw[,1], fw[,2], pch=20, col="blue", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3, xlab="$tophseq ($qpref)", ylab="$Q[0] ($refgenome, $Q[1], $Q[2] - $Q[3])")
points(rv[,1], rv[,2], pch=20, col="red", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3)
abline(h=$tpos_min, col="black")
abline(h=$tpos_max, col="black")
abline(v=$qpos_min, col="black")
abline(v=$qpos_max, col="black")
$hline
dev.off()

png("$png2", width=2400, height=2400)
plot(fw[,1], fw[,2], pch=20, col="red", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3, xlab="$tophseq ($qpref)", ylab="$Q[0] ($refgenome, $Q[1], $Q[2] - $Q[3])")
points(rv[,1], rv[,2], pch=20, col="blue", xlim=c($qpos_min, $qpos_max), ylim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3)
abline(h=$tpos_min, col="black")
abline(h=$tpos_max, col="black")
abline(v=$qpos_min, col="black")
abline(v=$qpos_max, col="black")
$hline
dev.off()

png("$png3", width=2400, height=2400)
plot(fw[,2], fw[,1], pch=20, col="blue", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3, xlab="$Q[0] ($refgenome, $Q[1], $Q[2] - $Q[3])", ylab="$tophseq ($qpref)")
points(rv[,2], rv[,1], pch=20, col="red", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_min, $tpos_max), cex.axis=2, cex.lab=1.5, cex=3)
abline(v=$tpos_min, col="black")
abline(v=$tpos_max, col="black")
abline(h=$qpos_min, col="black")
abline(h=$qpos_max, col="black")
$vline
dev.off()

png("$png4", width=2400, height=2400)
plot(fw[,2], fw[,1], pch=20, col="red", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3, xlab="$Q[0] ($refgenome, $Q[1], $Q[2] - $Q[3])", ylab="$tophseq ($qpref)")
points(rv[,2], rv[,1], pch=20, col="blue", ylim=c($qpos_min, $qpos_max), xlim=c($tpos_max, $tpos_min), cex.axis=2, cex.lab=1.5, cex=3)
abline(v=$tpos_min, col="black")
abline(v=$tpos_max, col="black")
abline(h=$qpos_min, col="black")
abline(h=$qpos_max, col="black")
$vline
dev.off()

EOS
					
					if($dotdata_flanking_fw){
						open(my $rfh, ">", "./$outplot_dir/$tmp_dotcsv_fw");
						print $rfh $dotdata_flanking_fw;
						close $rfh;
					}
					if($dotdata_flanking_rv){
						open(my $rfh, ">", "./$outplot_dir/$tmp_dotcsv_rv");
						print $rfh $dotdata_flanking_rv;
						close $rfh;
					}
					if($Rscript){
						open(my $rfh, ">", $tmp_Rscript);
						print $rfh $Rscript;
						close $rfh;
					}
					
					my $Rcmd = "Rscript --vanilla --slave $tmp_Rscript >/dev/null 2>&1";
					#print "! cmd=[$Rcmd]\n";
					system("$Rcmd");
					
					if(-e "Rplots.pdf"){
						system("rm Rplots.pdf");
					}
				}
				#-----------------------------------------------------------------------//
				
				$summary .= $Q[0]."\t".$refgenome."\t".$Q[1]."\t".$Q[2]."\t".$Q[3]."\t".$initspan."\tref\t".$hsample->{$j}."\t".$gphitbp0."\t".$gphitbp1."\t".$tophseq."\t".$hhitinfo->{$j}{$tophseq}{pos0}."\t".$hhitinfo->{$j}{$tophseq}{pos1}."\t".$limit_pos0."\t".$limit_pos1."\t".$num_HCRS;
				$summary .= "\n";
				
				if(! $Qbp0[0] || ! $Qbp0[$gphitbp0 - 1]){
					$failed_list .= $Q[0]."\tmissing Qbp0 (found alignment within the flanking region but cannot find gene position)\n";
					
					my $gpos0f = $Q[2] - $initspan;
					my $gpos1f = $Q[3] + $initspan;
					$summary2 .= $hsample->{$j}."\t".$refgenome."\t".$Q[1]."\t".$gpos0f."\t".$gpos1f."\t".$Q[2]."\t".$Q[3]."\t".$hsample->{$j}."\t-\t-\t-\t-\t";
					$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tabsent\t".$Q[0]."\t";
					$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\n";
					next;
				}
				
				my $bp_sp01 = abs($Q[2] - $Q[3]) + 1;
				my $bp_qsp01 = abs($Qbp0[0] - $Qbp0[$gphitbp0 - 1]) + 1;
				my $bp_ins = $bp_qsp01 - $bp_sp01;
				if($bp_ins < 1){
					$bp_ins = "-";
				}
				my $ratio_galign = $gphitbp0 / $bp_sp01;
				my $ratio_seqnormality = $ratio_galign;
				
				my $ratio_ins = "-";
				if($bp_qsp01 > 0){
					$ratio_ins = $bp_sp01 / $bp_qsp01;
					
					if($ratio_ins >= 1){
						$ratio_ins = "-";
					}
					elsif($ratio_seqnormality > $ratio_ins){
						$ratio_seqnormality = $ratio_ins;
					}
				}
				
				my $bp_5neighbor = 0;
				my $bp_3neighbor = 0;
				my $bp_5flanking = 0;
				my $bp_3flanking = 0;
				my $ratio_5neighbor = 0;
				my $ratio_3neighbor = 0;
				my $ratio_5flanking = 0;
				my $ratio_3flanking = 0;
				my $rspan_5neighbor = 0;
				my $rspan_3neighbor = 0;
				
				if($gphitbp2 > 1 && abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]) > 0){
					if($Q[4] eq '+'){
						$bp_5neighbor = $gphitbp2;
						$ratio_5neighbor = sprintf("%.3f", $gphitbp2 / abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]) );
						
						if(abs($Qbp2[0] - $Qbp2[$gphitbp2 - 1]) > 0){
							if(defined $Qbp0[0] && $Qbp0[0] > $Qbp2[$gphitbp2 - 1]){
								$rspan_5neighbor = sprintf("%.3f", abs($Qbp2[0] - $Qbp2[$gphitbp2 - 1]) / abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]) );
								$rspan_5neighbor .= ";";
								$rspan_5neighbor .= sprintf("%.3f", ( abs($Qbp2[0] - $Qbp2[$gphitbp2 - 1]) + abs($Qbp0[0] - $Qbp2[$gphitbp2 - 1]) ) / abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]) );
							}
							else{
								$rspan_5neighbor = sprintf("%.3f", abs($Qbp2[0] - $Qbp2[$gphitbp2 - 1]) / abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]) );
								$rspan_5neighbor .= ";null";
							}
							$ratio_5neighbor .= ";".$rspan_5neighbor.";null;".$Qbp2[0].";".$Qbp2[$gphitbp2 - 1];
						}
						else{
							$ratio_5neighbor .= ";null;null;null;null;null";
						}
						# ratio_hit; ratio_span; rev_ratio_span; P2Gdist; span_pos0; span_pos1
					}
					elsif($Q[4] eq '-'){
						$bp_3neighbor = $gphitbp2;
						$ratio_3neighbor = sprintf("%.3f", $gphitbp2 / abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]) );
						
						if(abs($Qbp2[0] - $Qbp2[$gphitbp2 - 1]) > 0){
							if(defined $Qbp0[0] && $Qbp0[0] > $Qbp2[$gphitbp2 - 1]){
								$rspan_3neighbor = sprintf("%.3f", abs($Qbp2[0] - $Qbp2[$gphitbp2 - 1]) / abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]) );
								$rspan_3neighbor .= ";";
								$rspan_3neighbor .= sprintf("%.3f", ( abs($Qbp2[0] - $Qbp2[$gphitbp2 - 1]) + abs($Qbp0[0] - $Qbp2[$gphitbp2 - 1]) ) / abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]) );
							}
							else{
								$rspan_3neighbor = sprintf("%.3f", abs($Qbp2[0] - $Qbp2[$gphitbp2 - 1]) / abs($Gbp2[0] - $Gbp2[$gphitbp2 - 1]) );
								$rspan_3neighbor .= ";null";
							}
							$ratio_3neighbor .= ";".$rspan_3neighbor.";null;".$Qbp2[0].";".$Qbp2[$gphitbp2 - 1];
						}
						else{
							$ratio_3neighbor .= ";null;null;null;null;null";
						}
						# ratio_hit; ratio_span; rev_ratio_span; P2Gdist; span_pos0; span_pos1
					}
				}
				if($gphitbp3 > 1 && abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]) > 0){
					if($Q[4] eq '+'){
						$bp_3neighbor = $gphitbp3;
						$ratio_3neighbor = sprintf("%.3f", $gphitbp3 / abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]) );
						
						if(abs($Qbp3[0] - $Qbp3[$gphitbp3 - 1]) > 0){
							if(defined $Qbp0[$gphitbp0 - 1] && $Qbp3[0] > $Qbp0[$gphitbp0 - 1]){
								$rspan_3neighbor = sprintf("%.3f", abs($Qbp3[0] - $Qbp3[$gphitbp3 - 1]) / abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]) );
								$rspan_3neighbor .= ";";
								$rspan_3neighbor .= sprintf("%.3f", ( abs($Qbp3[0] - $Qbp3[$gphitbp3 - 1]) + abs($Qbp3[0] - $Qbp0[$gphitbp0 - 1]) ) / abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]) );
							}
							else{
								$rspan_3neighbor = sprintf("%.3f", abs($Qbp3[0] - $Qbp3[$gphitbp3 - 1]) / abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]) );
								$rspan_3neighbor .= ";null";
							}
							$ratio_3neighbor .= ";".$rspan_3neighbor.";null;".$Qbp3[0].";".$Qbp3[$gphitbp3 - 1];
						}
						else{
							$ratio_3neighbor .= ";null;null;null;null;null";
						}
						# ratio_hit; ratio_span; rev_ratio_span; P2Gdist; span_pos0; span_pos1
					}
					elsif($Q[4] eq '-'){
						$bp_5neighbor = $gphitbp3;
						$ratio_5neighbor = sprintf("%.3f", $gphitbp3 / abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]) );
						
						if(abs($Qbp3[0] - $Qbp3[$gphitbp3 - 1]) > 0){
							if(defined $Qbp0[$gphitbp0 - 1] && $Qbp3[0] > $Qbp0[$gphitbp0 - 1]){
								$rspan_5neighbor = sprintf("%.3f", abs($Qbp3[0] - $Qbp3[$gphitbp3 - 1]) / abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]) );
								$rspan_5neighbor .= ";";
								$rspan_5neighbor .= sprintf("%.3f", ( abs($Qbp3[0] - $Qbp3[$gphitbp3 - 1]) + abs($Qbp3[0] - $Qbp0[$gphitbp0 - 1]) ) / abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]) );
							}
							else{
								$rspan_5neighbor = sprintf("%.3f", abs($Qbp3[0] - $Qbp3[$gphitbp3 - 1]) / abs($Gbp3[0] - $Gbp3[$gphitbp3 - 1]) );
								$rspan_5neighbor .= ";null";
							}
							$ratio_5neighbor .= ";".$rspan_5neighbor.";null;".$Qbp3[0].";".$Qbp3[$gphitbp3 - 1];
						}
						else{
							$ratio_5neighbor .= ";null;null;null;null;null";
						}
						# ratio_hit; ratio_span; rev_ratio_span; P2Gdist; span_pos0; span_pos1
					}
				}
				if($gphitbp4 > 1 && abs($Gbp4[0] - $Gbp4[$gphitbp4 - 1]) > 0){
					if($Q[4] eq '+'){
						$bp_5flanking = $gphitbp4;
						$ratio_5flanking = sprintf("%.3f", $gphitbp4 / abs($Gbp4[0] - $Gbp4[$gphitbp4 - 1]) );
					}
					elsif($Q[4] eq '-'){
						$bp_3flanking = $gphitbp4;
						$ratio_3flanking = sprintf("%.3f", $gphitbp4 / abs($Gbp4[0] - $Gbp4[$gphitbp4 - 1]) );
					}
				}
				if($gphitbp5 > 1 && abs($Gbp5[0] - $Gbp5[$gphitbp5 - 1]) > 0){
					if($Q[4] eq '+'){
						$bp_3flanking = $gphitbp5;
						$ratio_3flanking = sprintf("%.3f", $gphitbp5 / abs($Gbp5[0] - $Gbp5[$gphitbp5 - 1]) );
					}
					elsif($Q[4] eq '-'){
						$bp_5flanking = $gphitbp5;
						$ratio_5flanking = sprintf("%.3f", $gphitbp5 / abs($Gbp5[0] - $Gbp5[$gphitbp5 - 1]) );
					}
				}
				
				my $protein_coding = "FALSE";
				if($Q[6] && $Q[6] > 0){
					$protein_coding = "TRUE";
				}
				
				my $qseq_woN = "-";
				my $bp_qseq_woN = "-";
				if(! defined $hqseq->{$tophseq}{len}){
					print "! error : missing seq info for [$tophseq]\n";
					die;
				}
				if($hqseq->{$tophseq}{len} > $Qbp0[$gphitbp0 - 1]){
					$qseq_woN = substr($hqseq->{$tophseq}{seq}, $Qbp0[0] - 1, $bp_qsp01);
					$qseq_woN =~ s/N//gi;
					$bp_qseq_woN = length($qseq_woN);
				}
				else{
					$qseq_woN = substr($hqseq->{$tophseq}{seq}, $Qbp0[0] - 1);
					$qseq_woN =~ s/N//gi;
					$bp_qseq_woN = length($qseq_woN);
				}
				
				my $judge = "present";
				if($bp_ins ne '-'){
					if($bp_ins > 5000){
						$judge = "insertion > 5kb";
					}
					elsif($bp_ins > 4000){
						$judge = "insertion > 4kb";
					}
					elsif($bp_ins > 3000){
						$judge = "insertion > 3kb";
					}
					elsif($bp_ins > 2000){
						$judge = "insertion > 2kb";
					}
					elsif($bp_ins > 1000){
						$judge = "insertion > 1kb";
					}
					elsif($bp_ins > 500){
						$judge = "insertion > 500bp";
					}
				}
				elsif($ratio_seqnormality < 0.1){
					$judge = "present (collapsed)";
				}
				elsif($ratio_seqnormality < 0.5){
					$judge = "present (collapsed)";
				}
				elsif($ratio_seqnormality < 0.8){
					$judge = "present (partly)";
				}
				
				$summary2 .= $hsample->{$j}."\t".$refgenome."\t".$Q[1]."\t".$gpos0f."\t".$gpos1f."\t".$Q[2]."\t".$Q[3]."\t".$hsample->{$j}."\t".$tophseq."\t".$Qbp1[0]."\t".$Qbp1[$gphitbp1 - 1]."\t1\t";
				$summary2 .= $Qbp0[0]."\t".$Qbp0[$gphitbp0 - 1]."\t".$bp_sp01."\t".$gphitbp0."\t".$bp_qsp01."\t".$bp_qseq_woN."\t".$bp_ins."\t".$ratio_galign."\t".$ratio_ins."\t".$ratio_seqnormality."\t".$judge."\t".$Q[0]."\t";
				$summary2 .= $bp_5neighbor."\t".$ratio_5neighbor."\t".$bp_3neighbor."\t".$ratio_3neighbor."\t".$bp_5flanking."\t".$ratio_5flanking."\t".$bp_3flanking."\t".$ratio_3flanking."\n";
			}
			else{
				$summary .= $Q[0]."\t".$refgenome."\t".$Q[1]."\t".$Q[2]."\t".$Q[3]."\t".$initspan."\tref\t".$hsample->{$j}."\t-\t-\t-\t-\t-\t-\t-\t-\n";
				
				my $gpos0f = $Q[2] - $initspan;
				my $gpos1f = $Q[3] + $initspan;
				$summary2 .= $hsample->{$j}."\t".$refgenome."\t".$Q[1]."\t".$gpos0f."\t".$gpos1f."\t".$Q[2]."\t".$Q[3]."\t".$hsample->{$j}."\t-\t-\t-\t-\t";
				$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tabsent\t".$Q[0]."\t";
				$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\n";
			}
		}
		else{
			$summary .= $Q[0]."\t".$refgenome."\t".$Q[1]."\t".$Q[2]."\t".$Q[3]."\t".$initspan."\tref\t".$hsample->{$j}."\t-\t-\t-\t-\t-\t-\t-\t-\n";
			
			my $gpos0f = $Q[2] - $initspan;
			my $gpos1f = $Q[3] + $initspan;
			$summary2 .= $hsample->{$j}."\t".$refgenome."\t".$Q[1]."\t".$gpos0f."\t".$gpos1f."\t".$Q[2]."\t".$Q[3]."\t".$hsample->{$j}."\t-\t-\t-\t-\t";
			$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tabsent\t".$Q[0]."\t";
			$summary2 .= "-\t-\t-\t-\t-\t-\t-\t-\n";
		}
	}
	
	$cntq++;
	if($cntq == $seg){
		print " partition [$t] : [$cntq] [$cnt_selected]\n";
		$seg += $add;
	}
}

print "! partition [$t] : [$cnt_hit] candidates, [$cnt_selected] hits\n";

my $rh = {};
$rh->{hit} = $hit;
$rh->{summary} = $summary;
$rh->{summary2} = $summary2;
$rh->{failed_list} = $failed_list;

return $rh;
}


#-------------------------------------------------------------------------------
sub Read_list{
my $file = shift;

print "! reading ID list [$file]...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
my $cnt = 0;
while(my $line = <$fh>){
	$line =~ s/\n//;
	$line =~ s/\r//;
	
	unless($line){
#		print "! missing line\n";
		next;
	}
	elsif($line =~ /\#/){
		next;
	}
	elsif($line =~ /gid/ && $line =~ /chr/i){
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
	
	if($A[0]){
		$hash->{$A[0]} = 1;
		$cnt++;
	}
}
close $fh;

print " [$cnt] entries\n";

return $hash;
}


#-------------------------------------------------------------------------------
sub Read_gff{
my $gff3 = shift;
my $hID = shift;

print "! reading [$gff3]...\n";
open(my $fh, "<", $gff3) or die;
my $g2t = {};
my $t2g = {};
my $hash = {};
my @GID;
my $gcnt = 0;
my @L;
while(my $line = <$fh>){
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
		
		if($hID->{$gid}){
			push(@GID, $gid);
			$hash->{$gid}{data} = $gid.",".$A[0].",".$A[3].",".$A[4].",".$A[6];
		}
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
		
		if($hID->{$gid}){
			$g2t->{$gid} .= $tid."\n";
			$t2g->{$tid} = $gid;
			$hash->{$gid}{num_variant} += 1;
		}
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
			next;
		}
		
		my $gid = $t2g->{$tid};
		if($hID->{$gid}){
			$hash->{$gid}{num_CDS} += 1;
		}
	}
}

print " [$gcnt] genes\n";

#my $r = "gid,chr,pos0,pos1,strand,num_variant,total_num_CDS\n";
my $r = "";
my $numpc = 0;
my $numnc = 0;
my $numpc_tr = 0;
my $hcnt = {};
foreach my $gid (@GID){
	unless($hash->{$gid}{num_variant}){
		$hash->{$gid}{num_variant} = 0;
	}
	unless($hash->{$gid}{num_CDS}){
		$numnc++;
		$hash->{$gid}{num_CDS} = 0;
		$hash->{$gid}{protein_coding} = "false";
	}
	else{
		$numpc++;
		$numpc_tr += $hash->{$gid}{num_variant};
		$hash->{$gid}{protein_coding} = "true";
	}
	
	my $tmp = $hash->{$gid}{data}.",".$hash->{$gid}{num_variant}.",".$hash->{$gid}{num_CDS};
	my @B = split(/,/, $tmp);
	$hcnt->{$B[1]} += 1;
	$r .= $tmp."\n";
}

print " [$numpc] protein-coding, [$numpc_tr] transcripts\n";
print " [$numnc] non-coding\n";

#my @Chrs = keys(%{$hcnt});
#@Chrs = sort {$a cmp $b} @Chrs;
#
#$r .= "\nChr,num gene\n";
#foreach my $sid (@Chrs){
#	$r .= $sid.",".$hcnt->{$sid}."\n";
#}

my @RL = split(/\n/, $r);
my $rh = {};
$rh->{RL} = \@RL;
$rh->{hash} = $hash;
$rh->{g2t} = $g2t;
$rh->{t2g} = $t2g;

return $rh;
}


#----------------------------------------------------------------------
sub LnFile{
my $file = shift;

if($file && $file ne 'null'){
	my @tmp = split(/\//, $file);
	my $n = @tmp;
	
	if(-e $tmp[$n - 1]){
		system("rm $tmp[$n - 1]");
	}
	
	if(-e $file){
		system("ln -s $file ./");
	}
	elsif(-e "../$file"){
		system("ln -s ../$file ./");
	}
	
	return $tmp[$n - 1];
}

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
		if($A[38] eq '-'){						# ORF information NOT attached
			$nnot++;
		}
		else{									# ORF information already attached
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

print "! reading PAV result [$file] ...\n";
open(my $fh, "<", $file) or die;
my $hash = {};
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
	if(defined $A[18]){
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

print "! reading [$file] as hash ...\n";
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
sub Batch_prepseq{
my ($bin_prepseq4annot, $combined_preplog, $random_IDlist, $genelist, $qfasta, $dbgff, $dbfasta, $dbprotein, $dbtranscript, $dbcds, $afile, $rfile_final, $n0, $n1, $tmpdir, $log_gth, $script, $combined_qprot, $tmp_ralign, $bin_findorf, $i) = @_;

my $cmd = "$bin_prepseq4annot $random_IDlist $genelist $qfasta $dbgff $dbfasta $dbprotein $dbtranscript $dbcds $afile $rfile_final $n0 $n1 $tmpdir $log_gth $script $combined_qprot $tmp_ralign $bin_findorf $i >> $combined_preplog 2>&1";

my $cnt = $n1 - $n0;
print " thread [$i] : [$cnt] queries...\n";
if(system($cmd) == 0){
	print " thread [$i] : completed.\n";
}
else{
	print "! thread [$i] : failed\n";
	die;
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


#---------------------------------------------------------------------
sub SAVE{
my ($file, $str, $header) = @_;

if($str && $header){
	$str = $header.$str;
	open(my $rfh, ">", $file) or die;
	print $rfh $str;
	close $rfh;
	print "! output [$file]\n";
}
elsif($str){
	open(my $rfh, ">", $file) or die;
	print $rfh $str;
	close $rfh;
}

}


#---------------------------------------------------------------------
sub ADD{
my ($file, $str) = @_;

if($str){
	open(my $rfh, ">>", $file);
	print $rfh $str;
	close $rfh;
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
			system("which $bin.pl >> $tmp 2>&1");
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



