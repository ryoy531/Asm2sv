# Asm2sv
Asm2sv is an assembly-based comparative genomics pipeline that analyzes gene-level structural variations (SV). It is designed to analyze SV between reference and target genomes of the same or closely related species. The basic idea of Asm2sv originates from the expectation that different types of SV can be present within a gene between distinct genomes. For example, some may carry 5-kb insertion within a gene region while others do 1-kb deletion. Because it is difficult to compare such SVs across multiple genomes based on conventional variant call format (VCF), we developed the Asm2sv as an alternative method. By algorithmic genomic alignment analysis, the Asm2sv pipeline captures insertion, deletion, or translocation around each gene region then output numeric scores that represent the degree of conservation (or disruption) for the gene. That is, a score close to 1.0 means that sequence of the corresponding gene is conserved when comapred to reference genome. The output SV scores can be united across multiple genomes to enable population-scale comparison. 
<br></br>

## System requirements
Minimum (small genome)  
- CPU: 16 core/32 threads (e.g. AMD ryzen 5950x)
- RAM: >128 GB

Recommend (genome size >500 Mb)
- CPU: 32 core/64 threads (e.g. AMD ryzen threadripper)
- RAM: >256 GB

*SSD raid disk is desirable to speed up analysis.  
<br></br>

## Software prerequisites
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

<sup>[*1]</sup> Some perl modules need to be installed with `cpan install` command.  
<sup>[*2]</sup> They must be in your PATH.  
<sup>[*3]</sup> These can be installed with `miniconda3` in which `bioconda` and `conda-forge` channels are added.  
<sup>[*4]</sup> These can be installed with `apt-get install` command.  
<sup>[*5]</sup> Included in the bin directory.  
<br></br>

## Installation
Download zip or type the following git command:
```
git clone https://github.com/ryoy531/Asm2sv.git
```

By the following command, you can check whether all required programs are in your PATH.
```
/path/to/Asm2sv check
```
<br></br>

## Command options
Asm2sv has multiple command options as follows. Detailed usage is described in the Tutorial section below.  

***
- `gfftolist` produces a gene query list file based on reference Gff3 file. The output list file can be used in the `run` command.
```
/path/to/Asm2sv gfftolist -g [reference Gff3]
```
***
- `run` analyzes gene-level SVs present in target genome using reference genome as a base.  
```
/path/to/Asm2sv run -d [reference fasta] -g [reference Gff3] -l [gene query list] -q [target fasta] -o [working directory] -t [CPU1] -x [CPU2] -n 5000
```
<sup>`-t` specifies CPU threads used for genomic alignment analysis.</sup>  
<sup>`-x` specifies CPU threads used for gene prediction with genome threader.</sup>  
<sup>`-n` specifies the length of promoter and 3'-UTR region for which SV scores will be calculated.</sup>  
***
- `makecmd` produces a suite of command lines for multiple target genomes.
```
/path/to/Asm2sv makecmd -i [list of target genome (.csv)] -t [CPU1] -x [CPU2] -n 5000
```
***
- `unite` combines results of multiple target genomes to compare genotype in a population scale.
```
/path/to/Asm2sv unite -i [list of target genome (.csv)] -c [information of chromosome ID alias (.tsv)]
```
***
- `plot` produces .png image file(s) that show genomic alignment of specified gene region including its flanking sequences.
```
/path/to/Asm2sv run -d [reference fasta] -g [reference Gff3] -q [list of target fasta (.csv or .txt)] -i [gene query list] -o [working directory] -t [CPU1] -f [flanking seq length] -n 5000
```
<sup>`-f` specifies the length of flanking sequence that is shown in genomic alignment plot (.png image file) together with gene region.</sup>  
<br></br>

## Tutorial
Here, we would like to show the usages of command options based on tutorial dataset. Please obtain it via [Daizu-net](https://daizu-net.dna.affrc.go.jp/ap/top)


The genome sequence data included in this tutorial originate from previously published data of other research groups. To use them as tutorial data, we have modified data to reduce file sizes. Please note that these files are different from the original ones.
