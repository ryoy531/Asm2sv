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




```
 this will be highlighted in green
- this will be highlighted in red
```


The genome sequence data included in this tutorial originate from previously published data of other research groups. To use them as tutorial data, we have modified data to reduce file sizes. Please note that these files are different from the original ones.
