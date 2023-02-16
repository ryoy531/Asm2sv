# Asm2sv

Asm2sv is an assembly-based comparative genomics pipeline designed to investigate the structural genomic variations that are occurring in flexible manner between distinct genomes. By algorithmic genomic alignment analysis, it captures insertion, deletion, or translocation in each gene including the flanking region then output numeric scores that represent the degree of conservation (or disruption) for all genes. The resultant output SV scores can be united across multiple individuals to enable population-scale comparison. The basic idea of Asm2sv originates from the expectation that different types of SV can be present within a gene between distinct genomes. For example, some may carry 5-kb insertion within a gene region while others do 1-kb deletion. Because it is difficult to compare such SVs based on conventional variant call format (VCF), we developed the Asm2sv pipeline.


## System requirements
Minimum (small genome)  
- CPU: 16 core/32 threads
- RAM: 128 GB

Recommend (genome size >500 Mb)
- CPU: 32 core/64 threads
- RAM: 256 GB

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
- samtools (1.12 or higher)<sup>[*2]</sup>
- gffread<sup>[*2]</sup>
- blat<sup>[*2]</sup>
- BLAST (2.10.1+ or higher)<sup>[*2]</sup>
- LAST<sup>[*3]</sup>
- matcher<sup>[*3]</sup>
- genome threader<sup>[*4]</sup>

<sup>[*1]</sup> Perl modules can be installed with `cpan install` command.  



```
 this will be highlighted in green
- this will be highlighted in red
```


The genome sequence data included in this tutorial originate from previously published data of other research groups. To use them as tutorial data, we have modified data to reduce file sizes. Please note that these files are different from the original ones.
