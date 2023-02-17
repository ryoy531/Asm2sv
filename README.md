# Asm2sv
Asm2sv is an assembly-based comparative genomics pipeline that analyzes gene-level structural variations (SV). It is designed to analyze SV between reference and target genomes of the same or closely related species. The basic idea of Asm2sv originates from the expectation that different types of SV can be present within a gene between distinct genomes. For example, some may carry 5-kb insertion within a gene region while others do 1-kb deletion. Because it is difficult to compare such SVs across multiple genomes based on conventional variant call format (VCF), we developed the Asm2sv as an alternative method. By algorithmic genomic alignment analysis, the Asm2sv pipeline captures insertion, deletion, or translocation around each gene region then output numeric scores that represent the degree of conservation (or disruption) for the gene. The output SV scores can be united across multiple genomes to enable population-scale comparison. 
  
[System requirements](#System_requirements)  
[Software prerequisites](#Software_prerequisites)  
[Installation](#Installation)  
[Command options](#Command_options)  
[Tutorial](#Tutorial)  
####

<h2 id="System_requirements"># System requirements</h2>

Minimum (small genome)  
- CPU: 16 core/32 threads (e.g. AMD ryzen 5950x)
- RAM: >128 GB

Recommend (genome size >500 Mb)
- CPU: 32 core/64 threads (e.g. AMD ryzen threadripper)
- RAM: >256 GB

*SSD disk is highly recommended in both cases.  
####

<h2 id="Software_prerequisites"># Software prerequisites</h2>

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
  - matcher<sup>[*4]</sup>
  - genome threader<sup>[*5]</sup>

<sup>[*1] Some perl modules need to be installed with `cpan install` command.</sup>  
<sup>[*2] They must be in your PATH.</sup>  
<sup>[*3] These can be installed with `miniconda3` in which `bioconda` and `conda-forge` channels are added.</sup>  
<sup>[*4] These can be installed with `apt-get install` command.</sup>  
<sup>[*5] Included in the bin directory.</sup> 
####

<h2 id="Installation"># Installation</h2>

Download zip or type the following git command:
```
$ git clone https://github.com/ryoy531/Asm2sv.git
```
####

By the following command, you can check whether all required programs are in your PATH.
```
$ /path/to/Asm2sv check
```
If it returns the following message, installation is no problem.
```
! Found all (looks good).
```
####

<h2 id="Command_options"># Command options</h2>

Asm2sv has multiple command options as follows. Detailed usage is described in the Tutorial section below.  
####
- `gfftolist` produces a gene query list file based on reference Gff3 file. The output list file can be used in the `run` command.
```
$ /path/to/Asm2sv gfftolist -g [reference Gff3]
```
####  
- `run` analyzes gene-level SVs present in target genome using reference genome as a base.  
```
$ /path/to/Asm2sv run -d [reference fasta] -g [reference Gff3] -l [gene query list] -q [target fasta] -o [working directory] -t [CPU1] -x [CPU2] -n 5000
```
<sup>`-t` specifies CPU threads used for genomic alignment analysis.</sup>  
<sup>`-x` specifies CPU threads used for gene prediction with genome threader.</sup>  
<sup>`-n` specifies the length of promoter and 3'-UTR region for which SV scores will be calculated.</sup>  
####
- `makecmd` produces a suite of command lines for multiple target genomes.
```
$ /path/to/Asm2sv makecmd -i [list of target genome (.csv)] -t [CPU1] -x [CPU2] -n 5000
```
####  
- `unite` combines results of multiple target genomes to compare genotype in a population scale.
```
$ /path/to/Asm2sv unite -i [list of target genome (.csv)] -c [information of chromosome ID alias (.tsv)]
```
####
- `plot` produces .png image file(s) that show genomic alignment of specified gene region including its flanking sequences.
```
$ /path/to/Asm2sv run -d [reference fasta] -g [reference Gff3] -q [list of target fasta (.csv or .txt)] -i [gene query list] -o [working directory] -t [CPU1] -f [flanking seq length] -n 5000
```
<sup>`-f` specifies the length of flanking sequence that is shown in genomic alignment plot (.png image file) together with gene region.</sup>  
####
  
<h2 id="Tutorial"># Tutorial</h2>

Basic usage of Asm2sv is described below.  
___
#### Step.1 Move to the directory `tutorial` then prepare a gene query list file from Gff.
```
$ cd tutorial  
$ ls      #check files  
  chrname_info.tsv  
  list_for_batch_exec.csv  
  reference_genome.fasta  
  reference_genome.gff3  
  sample_genome_1.fasta  
  sample_genome_2.fasta  
  sample_genome_3.fasta 
  
$ /path/to/Asm2sv gfftolist -g reference_genome.gff3      #This will create query list csv
```
<sup>[Note] Some genome fasta files in the tutorial directory originate from publicly available data. We have modified them to reduce the file sizes.</sup>    

The above command will create a csv file named `summary_gene_reference_genome.csv`. Important columns are described below:  
####
| column | string | meaning |
| ------ | ------ | ------- |
| 1 | gid | gene ID |
| 2 | chr | chromosome or sequence name |
| 3 | pos0 | start position |
| 4 | pos1 | end position |
| 5 | strand | direction |

<sup>[Note] It is also able to use manually prepared list table.</sup>    
####

___
#### Step.2 Run Asm2sv to compare reference and target genomes (an example of one-to-one comparison).  
<sup>[Note] In case of comparing multiple target genomes, you may not need to perform this step. Please just read it.</sup>   
```
$ /path/to/Asm2sv run -d reference_genome.fasta -g reference_genome.gff3 -l summary_gene_reference_genome.csv -q sample_genome_1.fasta -o asm2sv_genome_1 -t 16 -x 16 -n 5000  
```
####

The result file will be `asm2sv_genome_1/rev_summary_genome2sv_sample_genome_1.tsv`. Some important columns are described below:  
####
| column | string | meaning |
| ------ | ------ | ------- |
| 9 | qfasta seqid | sequence name in which candidate gene is found in target genome |
| 13 / 14 | hit sp0 / hit sp1 | candidate gene of target genome locates between `hit sp0` to `hit sp1` |
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
| 26 | 5'-promoter align ratio (5000 bp) | additional Indel information (delimited by semicolon) |
| 27 | 3'-UTR bp hit (5000 bp) | bp alignment length of UTR region |
| 28 | 3'-UTR align ratio (5000 bp) | additional Indel information (delimited by semicolon) |
| 34 | db protein length | protein sequence length of reference gene |
| 35 | q protein length (predicted) | protein sequence length of candidate gene in target genome |
| 36 | align length | length of protein alignment between reference and target genes |
| 37 | %similarity | %protein sequence similarity between reference and target genes |
| 38 | %gap | %gappped sequence between reference and target genes |
| 39 | score | alignment score of matcher |

<sup>[Note] Protein sequence is predicted in target genome regardless of its Gff (not use it).</sup>    
####
- Among the above columns, scores in the columns 20-22 are important as they represent the degree of conservation or disruption for the gene. If score is close to 1.0, it means gene sequence is conserved in target genome.  
- In case of comparing multiple target genomes, these scores can be combined to generate a numeric genotype table (SV score table) through the `unite` command (see below).
####

___
#### Step.3 Prepare a bash script for batch execution of Asm2sv for multiple target genomes.  
<sup>[Note] In case of analyzing only one target genome, you don't need to perform this step. Please just read it.</sup>  
```
$ /path/to/Asm2sv makecmd -i list_for_batch_exec.csv -t 16 -x 16 -n 5000  
```
####

This command will generate a bash script file named `cmd_asm2sv_list_for_batch_exec.sh ` in which multiple command lines are described.  
```
$ cat cmd_asm2sv_list_for_batch_exec.sh  
  /usr/local/pipeline/r6c1_Asm2sv/Asm2sv run -d reference_genome.fasta -g reference_genome.gff3 -l summary_gene_reference_genome.csv -q reference_genome.fasta -t 16 --neighbor 5000 -o asm2sv_reference
  /usr/local/pipeline/r6c1_Asm2sv/Asm2sv run -d reference_genome.fasta -g reference_genome.gff3 -l summary_gene_reference_genome.csv -q sample_genome_1.fasta -t 16 --neighbor 5000 -o asm2sv_genome_1
  /usr/local/pipeline/r6c1_Asm2sv/Asm2sv run -d reference_genome.fasta -g reference_genome.gff3 -l summary_gene_reference_genome.csv -q sample_genome_2.fasta -t 16 --neighbor 5000 -o asm2sv_genome_2
  /usr/local/pipeline/r6c1_Asm2sv/Asm2sv run -d reference_genome.fasta -g reference_genome.gff3 -l summary_gene_reference_genome.csv -q sample_genome_3.fasta -t 16 --neighbor 5000 -o asm2sv_genome_3
```
####

Then just run a script.  
```
$ bash cmd_asm2sv_list_for_batch_exec.sh
```
<sup>[Note] In case of large genome like 1 Gb, one command will take ~48 hours. Please use many-core CPU or PC cluster.</sup>  
####

In the above command line of `Asm2sv makecmd`, a csv file named `list_for_batch_exec.csv` is specified as a list. Information in this file is described as below.
####
| column | string | meaning |
| ------ | ------ | ------- |
| 1 | dbfasta | reference fasta |
| 2 | dbgff | reference Gff3 |
| 3 | genelist | gene query list |
| 4 | qfasta | target fasta |
| 5 | qpsl | `null` or specify .psl alignment file (optional) |
| 6 | outdir | output directory |
####
Genomic alignment .psl file can be obtained with the following commands (may require >1 TB RAM in case of >500 Mb genome). If absent, just describe `null` in the column 4.  
```
$ lastdb [reference prefix] [reference fasta]  
$ lastal [reference prefix] [target fasta] | last-split -s 35 > [alignment maf]  
$ maf-convert psl [alignment maf] > [alignment psl]
```
####


___
#### Step.4 Unite the results of multiple target genomes to generate a genotype table.  
If the Step. 3 is successful, the following files will be created.
```
$ ls asm2sv_*/rev_summary_*  
  asm2sv_genome_1/rev_summary_genome2sv_sample_genome_1.tsv  
  asm2sv_genome_2/rev_summary_genome2sv_sample_genome_2.tsv  
  asm2sv_genome_3/rev_summary_genome2sv_sample_genome_3.tsv  
  asm2sv_reference/rev_summary_genome2sv_reference_genome.tsv
```
####

By the following command, they can be united to generate a genotype table. A directory named `combined_asm2sv` will be created.  
```
$ /path/to/Asm2sv unite -i list_for_batch_exec.csv -c chrname_info.tsv
```
####

Dependent on the data summarizing policies, there are several output files.  
- `A1`: Score approaches ‘0’ dependent on the degree of disruption regardless of indel.  
- `A2`: Lower score (0 to 1) for deletion, larger score (>1) for insertion.  
- `B1`: Output PAV scores for all genes (raw data including missing genotype).  
- `B2`: Output PAV scores when structural variations were precisely determined by genomic alignments in all target assemblies (without missing genotype).  
####

For example, the following data files are based on different combination of `policy A1 / A2` and `B1 / B2`.  
| file name | combination |
| --------- | ----------- |
| val_disrupt_q-4.csv | `A1` x `B1` |
| val_disrupt_1cnsv_q-4.csv | `A1` x `B2` |
| val_indel_q-4.csv | `A2` x `B1` |
| val_indel_1cnsv_q-4.csv | `A2` x `B2` |



