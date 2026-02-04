# colitis_memory
Code associated with the manuscript:
Clonal memory of colitis accumulates and promotes tumor growth. Surya Nagaraja, Lety Ojeda-Miron, Ruochi Zhang, Ena Oreskovic, Conrad Hock, Yan Hu, Daniel Zeve, Karina Sharma, Roni R. Hyman, Qiming Zhang, Andrew Castillo, David T. Breault, Ömer H. Yilmaz, Jason D. Buenrostro
in press at Nature, as of 2026-02-04

Preprint:
bioRxiv 2025.02.13.638099; doi: https://doi.org/10.1101/2025.02.13.638099

Single cell data is hosted by IGVF at:
https://tinyurl.com/nagaraja-buenrostro-2025

Bulk and spatial data available via GEO as part of series GSE316619.

For processing of raw fastq files and alignment, see SHARE-seq V2 pipeline:
https://github.com/masai1116/SHARE-seq-alignmentV2/

For footprinting and PRINT related software:
https://github.com/buenrostrolab/PRINT
https://github.com/buenrostrolab/scPrinter (version 1.2.0)

Please see tools directory for specific software and packages required. Input and output data used are included in the ref directory. Please see code comments for further instructions.

Versions for key tools:  
chromVAR 1.24.0  
BuenRTools 0.1.4  
cisTopic 0.3.0
motifMatchr 1.12.0
DeepLIFT custom implementation included in scprinter v1.2.0
tfmodisco-lite v2.2.1
memesuite v5.5.7
DESeq2 1.42.1
nf-core 4.0.0
MethylDackel 0.6.1
scikit 1.5.2
finemo 0.40
