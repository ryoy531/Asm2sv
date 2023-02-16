# Asm2sv
Asm2sv is an assembly-based comparative genomics pipeline that analyzes gene-level structural variations (SV) present between reference and target genomes (same species or closely related ones). The basic idea of Asm2sv originates from the expectation that different types of SV can be present within a gene between distinct genomes. For example, some may carry 5-kb insertion within a gene region while others do 1-kb deletion. Indeed, it is difficult to compare such SVs across multiple genomes based on conventional variant call format (VCF). By algorithmic genomic alignment analysis, the Asm2sv captures insertion, deletion, or translocation within gene region and its flanking region then output numeric scores that represent the degree of conservation (or disruption) in each gene. That is, a score close to 1.0 means that sequence of the corresponding gene is conserved when comapred to reference genome. The output SV scores can be united across multiple genomes to enable population-scale comparison. 

## System requirements
Minimum (small genome)  
- CPU: 16 core/32 threads (e.g. AMD ryzen 5950x)
- RAM: >128 GB

Recommend (genome size >500 Mb)
- CPU: 32 core/64 threads (e.g. AMD ryzen threadripper)
- RAM: >256 GB

*SSD raid disk is desirable to speed up analysis.  

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

<sup>[*1]</sup> Perl modules can be installed with `cpan install` command.  
<sup>[*2]</sup> They must be in your PATH.  
<sup>[*3]</sup> These can be installed with `miniconda3` in which `bioconda` and `conda-forge` channels are added.  
<sup>[*4]</sup> These can be installed with `apt-get install` command.  
<sup>[*5]</sup> Included in the bin directory.  

## Installation


## Commands and tutorial


```
 this will be highlighted in green
- this will be highlighted in red
```


The genome sequence data included in this tutorial originate from previously published data of other research groups. To use them as tutorial data, we have modified data to reduce file sizes. Please note that these files are different from the original ones.
