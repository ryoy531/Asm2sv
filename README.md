<a name="TOP"></a>
# Asm2sv
Asm2sv is an assembly-based comparative genomics pipeline designed to analyze gene-level structural variations (SV) between distinct genomes. In Asm2sv, one reference genome is set as a base then analyze gene-level SV in target genome(s). This pipeline is assumted for comparison of genomes in the same or closely related species. The basic idea of Asm2sv originates from the expectation that SV can occur flexibly between distinct genomes. For example, some may carry 5-kb insertion within a gene region while others do 1-kb deletion. Because it is difficult to estimate and compare the effect of such SVs on gene function between multiple genomes based solely on conventional variant call format (VCF), we developed the Asm2sv as an alternative method. By algorithmic genomic alignment analysis, the Asm2sv pipeline captures insertion, deletion, or translocation in each gene including its flanking region (promoter, 3'-UTR) then output numeric scores that represent the degree of conservation or disruption for the gene. The output SV scores can be united across multiple genomes to enable population-scale comparison. Asm2sv can also output VCF file that describes gene-level SV, which can be used for pangenome reference construction.

[System requirements](#System_requirements)  
[Software prerequisites](#Software_prerequisites)  
[Installation](#Installation)  
[Tutorial](#Tutorial)  
[Citation](#Citation)  

####
<a name="System_requirements"></a>
<h2>System requirements</h2>

Minimum (small genome)  
- CPU: 16 core/32 threads (e.g. AMD ryzen 5950x)
- RAM: 128 GB

Recommend (genome size >500 Mb)
- CPU: 64 core/128 threads (e.g. AMD ryzen threadripper)
- RAM: 512 GB

*SSD disk is recommended in both cases.  
####

<a name="Software_prerequisites"></a>
<h2>Software prerequisites</h2>

We have validated the Asm2sv pipeline with the following environment:  
- Ubuntu 20.04 LTS
- perl 5.30 or higher
- perl modules<sup>[*1]</sup>: 
  - strict
  - warnings
  - threads
  - Cwd
  - FindBin
  - Getopt::Long
  - Sys::Hostname
- Third-party programs<sup>[*2]</sup>
  - samtools (1.12 or higher)<sup>[*3]</sup>
  - gffread<sup>[*3]</sup>
  - blat<sup>[*3]</sup>
  - BLAST (2.10.1+ or higher)<sup>[*3]</sup>
  - LAST<sup>[*4]</sup>
  - matcher (from emboss)<sup>[*4]</sup>
  - minimap2 (2.26-r1175 or higher)<sup>[*3] [*a]</sup>
  - Sniffles2 (2.2 or higher)<sup>[*3] [*a]</sup>
  - genome threader<sup>[*5]</sup>
  - miniprot (0.12-r237 or higher)<sup>[*6]</sup>  
  - Minigraph-cactus (6.0.0)<sup>[*6] [*b]</sup>  
  - Vg (v1.54.0-1.55.0)<sup>[*6] [*b]</sup>  

<sup>[*1] It may be necessary to install some modules with `cpan install` command.</sup>  
<sup>[*2] They must be in your PATH.</sup>  
<sup>[*3] These can be installed with `miniconda3`. `bioconda` and `conda-forge` channels are required.</sup>  
<sup>[*4] These can be installed with `apt-get install` command.</sup>  
<sup>[*5] Included in the bin directory.</sup>  
<sup>[*6] Download from the github or web site of deverlopers.</sup>  
<sup>[*a] Required if using '--vcf' function.</sup>  
<sup>[*b] Required if preparing and using VCF for pangenome construction.</sup>  
  
<sup>[[Back to TOP]](#TOP)</sup>  
####

<a name="Installation"></a>
<h2 id="Installation">Installation</h2>

Download zip or type the following git command:
```
$ git clone https://github.com/ryoy531/Asm2sv.git
```
####
Change permission:
```
$ chmod 755 -R /path/to/Asm2sv
```
By the following command, you can check whether all required programs are in your PATH.
```
$ /path/to/Asm2sv check
```
If it returns the following message, installation is no problem.
```
! Found all (looks good).
```
  
<sup>[[Back to TOP]](#TOP)</sup>  
####

<a name="Tutorial"></a>
<h2>Tutorial</h2>

The fastest way to explain about the usage of Asm2sv is using tutorial data. Basic procedure to run Asm2sv is described in this section together with command lines.  
___
#### Step.1 Move to the directory `tutorial` then prepare a gene query list file from Gff.
```
$ cd tutorial  
$ ls      #check files  
  chrname_info.tsv  
  list_for_batch_exec.csv  
  reference.fasta, reference.gff  
  sample1.fasta, sample1.gff  
  sample2.fasta, sample2.gff  
  sample3.fasta, sample3.gff  
  sample4.fasta, sample4.gff  
  sample5.fasta, sample5.gff  
  sample6.fasta, sample6.gff  
  
$ /path/to/Asm2sv gfftolist -g reference.gff      #This will create query list csv
```
<sup>[Note] These genome fasta/GFF files originated from soybean genomes (see citaiton). These data was pruned to use in this tutorial.</sup>    

The above command will create a csv file named `summary_gene_reference.csv`. Columns in this data are described below:  
####
| column | string | meaning |
| ------ | ------ | ------- |
| 1 | gid | gene ID |
| 2 | chr | chromosome or sequence name |
| 3 | pos0 | start position |
| 4 | pos1 | end position |
| 5 | strand | direction |
| 6 | num_variant | number of transcript variant |
| 7 | total_num_CDS | number of CDS record |

<sup>[Note] You can pick up gene entries and delete others.</sup>    
####

___
#### Step.2A Run Asm2sv to compare reference and target genome (an example of one-by-one analysis).  
```
# command line example of comparing reference.fasta/GFF and sample1.fasta/GFF
$ /path/to/Asm2sv run -r reference.fasta -g reference.gff -l summary_gene_reference.csv -q sample1.fasta -s sample1.gff -o asm2sv_1 -t 16 -x 8 -n 5000  
```
<sup>*`-n 5000` means Asm2sv will analyze 5000 bp flanking region in each gene.<\sup>  
<sup>*To output VCF, add `--vcf -c chrname_info.tsv` to the command line. `chrname_info.tsv` is a file that describes "alias for sequence name" (see below).<\sup>  
<sup>*In case of using SGE grid, add `--1` to the command line.<\sup>  
####

The result file of the above command will be `./asm2sv_1/rev_summary_asm2sv_sample1.tsv`. Some important columns are described below:  
####
| column | string | meaning |
| ------ | ------ | ------- |
| 9 | qfasta seqid | sequence name in which candidate gene is found in target genome |
| 13 / 14 | hit sp0 / hit sp1 | best possible candidate gene is found in target genome between `hit sp0` to `hit sp1` |
| 15 | bp (sp0-sp1) | bp length of reference gene |
| 16 | bp hit-align (sp0-sp1) | bp alignment length |
| 17 | bp hit-span (sp0-sp1) | bp length between `hit sp0` to `hit sp1` |
| 18 | bp hit-span woN (sp0-sp1) | bp length between `hit sp0` to `hit sp1` except undetermined nucleotides |
| 19 | bp insertion | bp insertion length compared to reference gene |
| 20 | align ratio (1=not disrupted) | ratio of `column16 / column15` |
| 21 | insert ratio (1=not disrupted) | ratio of `column16 / column17` |
| 22 | seq normality (1=not disrupted) | a score representing the degree of gene conservation or disruption |
| 23 | judge | `present` if candidate gene is found in target genome |
| 24 | gene_id | gene ID of reference genome |
| 25 | 5'-promoter bp hit (5000 bp) | bp alignment length of promoter region |
| 26 | 5'-promoter align ratio (5000 bp) | additional SV information (delimited by semicolon) |
| 27 | 3'-UTR bp hit (5000 bp) | bp alignment length of UTR region |
| 28 | 3'-UTR align ratio (5000 bp) | additional SV information (delimited by semicolon) |
| 34 | db protein length | protein sequence length of reference gene |
| 35 | q protein length (predicted) | protein sequence length of best possible candidate gene in target genome |
| 36 | align length | length of protein alignment between reference and target genes |
| 37 | %similarity | %protein sequence similarity between reference and target genes |
| 38 | %gap | %gappped sequence between reference and target genes |
| 39 | score | matcher's score for protein sequence alignment |
####
- Among these columns, scores in the columns 20-23 are important as they represent the degree of conservation (or disruption) for the gene. If score is close to 1.0, it means gene sequence is conserved in target genome (judged as `present`). Meanings of other representative judge tags are as follows:  
 --  
 `present but not n-BDBH` = Candidate genomic region was found in the target genome for the query gene but it was not bidirectional best hit in BLASTn (possibly tandemly duplicated gene with similar function).  
 `present but missing ORF` = Bidirectional best hit (BDBH) genomic region was found in the target genome for the query gene but ORF prediction was faield possibly because it is basend on protein alignment.  
 `present (collapsed or partly)` = BDBH genomic region was found but it was partially collapsed.  
 `no hit but some...` = Candidate genomic region was only partially found but collapsed.  
 `absent or missing border` = Candidate genomic region was not found or collapsed.  
 --  
- In case of comparing multiple target genomes, these scores can be combined to generate a numeric genotype table (SV score table) through the `unite` command (see below).  
- Protein sequence is predicted in target genome regardless of its GFF (trial function).  
####

___
#### Step.2B Generate command lines for multiple target genomes (for parallel run of Asm2sv).  
By using the 'makecmd' option of Asm2sv, you can prepare command lines of multiple target queries at once. To use this, you need to prepare a list file (.csv) that describes information of path of data files for reference and target genomes.  
```
# command line example of makecmd  
$ /path/to/Asm2sv makecmd -l list_for_batch_exec.csv -c chrname_info.tsv --vcf -t 16 -x 8 -n 5000  
```
####
<sup>*`-c chrname_info.tsv --vcf` means Asm2sv will output VCF data in addition to numeric score data (optional).<\sup>  
  
This command will generate a bash script file named `cmd_asm2sv_list_for_batch_exec.sh` in which multiple command lines are described.  
```
$ cat cmd_asm2sv_list_for_batch_exec.sh  
  /path/to/Asm2sv run -r reference.fasta -g reference.gff -l summary_gene_reference.csv -q reference.fasta -t 16 -x 8 --neighbor 5000 -s reference.gff -o asm2sv_reference --dataID reference --vcf --chrinfo chrname_info.tsv  
  /path/to/Asm2sv run -r reference.fasta -g reference.gff -l summary_gene_reference.csv -q sample1.fasta -t 16 -x 8 --neighbor 5000 -s sample1.gff -o asm2sv_1 --dataID sample1 --vcf --chrinfo chrname_info.tsv  
  /path/to/Asm2sv run -r reference.fasta -g reference.gff -l summary_gene_reference.csv -q sample2.fasta -t 16 -x 8 --neighbor 5000 -s sample2.gff -o asm2sv_2 --dataID sample2 --vcf --chrinfo chrname_info.tsv  
  /path/to/Asm2sv run -r reference.fasta -g reference.gff -l summary_gene_reference.csv -q sample3.fasta -t 16 -x 8 --neighbor 5000 -s sample3.gff -o asm2sv_3 --dataID sample3 --vcf --chrinfo chrname_info.tsv  
  /path/to/Asm2sv run -r reference.fasta -g reference.gff -l summary_gene_reference.csv -q sample4.fasta -t 16 -x 8 --neighbor 5000 -s sample4.gff -o asm2sv_4 --dataID sample4 --vcf --chrinfo chrname_info.tsv  
  /path/to/Asm2sv run -r reference.fasta -g reference.gff -l summary_gene_reference.csv -q sample5.fasta -t 16 -x 8 --neighbor 5000 -s sample5.gff -o asm2sv_5 --dataID sample5 --vcf --chrinfo chrname_info.tsv  
  /path/to/Asm2sv run -r reference.fasta -g reference.gff -l summary_gene_reference.csv -q sample6.fasta -t 16 -x 8 --neighbor 5000 -s sample6.gff -o asm2sv_6 --dataID sample6 --vcf --chrinfo chrname_info.tsv  
```
####

Then just run a script to start analysis.  
```
$ bash cmd_asm2sv_list_for_batch_exec.sh
```
<sup>*In case of large genome like 1 Gb, one command will take ~48 hours. Please use high-performance PC or PC cluster.</sup>  
####

In `list_for_batch_exec.csv`, each column should include the information like below.
####
| column | string | meaning |
| ------ | ------ | ------- |
| 1 | dbfasta | reference fasta |
| 2 | dbgff | reference GFF |
| 3 | genelist | gene query list |
| 4 | qfasta | target fasta |
| 5 | qgff | target GFF |
| 6 | qpsl | `null` or specify .psl alignment file (optional) |
| 7 | outdir | output directory |
| 8 | custom_prefix | e.g. sample name |
####
  
In `chrname_info.tsv`, information should be described as follows.
| genome_prefix | alias_ID | original_ID |
| ------------- | -------- | ----------- |
| reference | chr01 | chr01 |
| reference | chr02 | chr02 |
| sample1 | chr01 | sample1_ch01 |
| sample1 | chr02 | sample1_ch02 |
| sample2 | chr01 | sample2_chr01 |
| sample2 | chr02 | sample2_chr02 |
| sample3 | chr01 | sample3ch1 |
| sample3 | chr02 | sample3ch2 |
| sample4 | chr01 | sample4chr1 |
| sample4 | chr02 | sample4chr2 |
| sample5 | chr01 | PREFIX1 |
| sample5 | chr02 | PREFIX2 |
| sample6 | chr01 | PREFIX202301 |
| sample6 | chr02 | PREFIX202302 |
####
<sup>*`genome_prefix` here is string except the file extention of fasta file (e.g. "sample1.fasta" -> "sample1").</sup>  
  
If all jobs are finished, the following files will be created.
```
# numeric score data for gene-SV
$ ls asm2sv_*/rev_summary_*  
  asm2sv_reference/rev_summary_asm2sv_reference.tsv  
  asm2sv_1/rev_summary_asm2sv_sample1.tsv
  asm2sv_2/rev_summary_asm2sv_sample2.tsv
  asm2sv_3/rev_summary_asm2sv_sample3.tsv
  asm2sv_4/rev_summary_asm2sv_sample4.tsv
  asm2sv_5/rev_summary_asm2sv_sample5.tsv  
  asm2sv_6/rev_summary_asm2sv_sample6.tsv

# VCF data for gene-SV
$ ls asm2sv_*/rev2sv/*.vcf | grep genebased  
  asm2sv_reference/rev2sv/genebased_asm2sv_reference.vcf
  asm2sv_1/rev2sv/genebased_asm2sv_sample1.vcf
  asm2sv_2/rev2sv/genebased_asm2sv_sample2.vcf
  asm2sv_3/rev2sv/genebased_asm2sv_sample3.vcf
  asm2sv_4/rev2sv/genebased_asm2sv_sample4.vcf
  asm2sv_5/rev2sv/genebased_asm2sv_sample5.vcf
  asm2sv_6/rev2sv/genebased_asm2sv_sample6.vcf
```
--  
[optional] To save computing  time, user can specify preliminary analyzed genomic alignment data as hint (.psl file). It can be obtained with the following commands (may require >1 TB RAM in case of >500 Mb genome). If absent, just describe `null` in the column 4.  
```
$ lastdb [reference prefix] [reference fasta]  
$ lastal [reference prefix] [target fasta] | last-split -s 35 > [alignment maf]  
$ maf-convert psl [alignment maf] > [alignment psl]
```
####

___
#### Step.3 Unite gene-SV numeric scores of multiple target genomes to generate an unified table.  

Following command will unify the result files of Step. 2B. New directory named `combined_asm2sv_list_for_batch_exec` will be created.  
```
$ /path/to/Asm2sv unite -i list_for_batch_exec.csv -c chrname_info.tsv
```
####
  
The `Asm2sv unite` command produces several kinds of output files. They are dependent on differenct summarizing policies.    
- `A1`: Score approaches ‘0’ dependent on the degree of disruption regardless of SV type (insertion or deletion).  
- `A2`: Lower score (0 to 1) for deletion, larger score (>1) for insertion.  
####

For example, the following data files are based on different combination of `policy A1 / A2` and `B1 / B2`.  
| file name | type | target |
| --------- | ----------- | ------ |
| sum_disrupt_q-7.csv | `A1` | gene, promoter, and UTR |
| gene_dirsupt_q-7.csv | `A1` | gene region |
| promoter_disrupt_q-7.csv | `A1` | promoter |
| 3UTR_disrupt_q-7.csv | `A1` | 3'-UTR |
| sum_indel_q-7.csv | `A2` | sum of gene, promoter, and UTR |
| gene_indel_q-7.csv | `A2` | gene region |
| promoter_indel_q-7.csv | `A2` | promoter |
| 3UTR_indel_q-7.csv | `A2` | 3'-UTR |
| protein_q-7.csv | `A2` | protein (*trial) |  
  
[Note] `A1` scores are assumed to be used to evaluate whether gene function is conserved between distinct genomes while `A2` scores are used to simply compare how sequences are different between genomes.  
  
<sup>[[Back to TOP]](#TOP)</sup>  
####

___
#### Step.4 (optional) Prepare VCF for pangenome reference construction.  

This step describes how to integrate Asm2sv's gene-SV VCF with those of cactus-pangenome then use it for pangenome reference construction. Because it relies on the functions of two third party tools `Minigraph-cactus` and `Vg`, user needs to cite these papers.  
  
First, run `splitseq_run_cactus.pl` from Asm2sv pipeline to prepare VCF with `cactus-pangenome`. This script splits reference fasta and obtain VCF in each sequence entry. 
```
$ /path/to/Asm2sv/scripts/splitseq_run_cactus.pl -l list_for_batch_exec.csv -c chrname_info.tsv -t 64 -p 2
```
####
<sup>*`-t` and `-p` specify the number of CPU threads and parallel partions, respectively.</sup>  
<sup>*This will create a VCF file such as `./split_run_cactuspg/cactus_pangenome.vcf`.</sup>  

Then, run `vcfintegrate_MCa2v.pl` to integrate Asm2sv's gene-SV VCF with those of cactus-pangenome. 
```
$ /path/to/Asm2sv/scripts/vcfintegrate_MCa2v.pl -b ./split_run_cactuspg/cactus_pangenome.vcf -l list_for_batch_exec.csv
```
####
<sup>*`-b` specifies the result VCF of first command.</sup>  
  
```
# check the result file
$ du -hs ./integratedVCF_for_pangenome/cactus_pangenome_plus_combined_geneSV_6genomes.vcf
  1.1M    ./integratedVCF_for_pangenome/cactus_pangenome_plus_combined_geneSV_6genomes.vcf
```
####

The resultant file `cactus_pangenome_plus_combined_geneSV_6genomes.vcf` includes not only SVs from cactus-pangenome but also those of Asm2sv. If they are overlapping, the script above keeps SVs of Asm2sv and discard another.
In the result directory, you can find a bash script file named `example_command_vgconst.sh`. By running it, user can construct pangenome reference with Vg.
```
$ cd integratedVCF_for_pangenome  
$ bash example_command_vgconst.sh
$ du -hs examplePanRef.*
  8.3M    examplePanRef.dist
  8.3M    examplePanRef.gbwt
  14M     examplePanRef.giraffe.gbz
  265M    examplePanRef.min
  8.0K    examplePanRef.snarls
  12M     examplePanRef.vg
  32M     examplePanRef.xg
```
####
Files named `examplePanRef.***` can be used for pangenome resequencing study with Vg.  
  
<sup>[[Back to TOP]](#TOP)</sup>  
####

<a name="Citation"></a>
<h2>Citation</h2>

This is pre-publication preview version that is aimed at evaluation before publication. If you access to applications and datasets, you agree to the following conditions.  
1. Please refrain from publication using the contens of this database and the genome reference dataset before our paper is published.
2. Please refrain from disclosure to third parties without our permission.
3. Pre-publication dataset may be updated or replaced without any notice.  



