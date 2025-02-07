# Ribosomal DNA arrays are the most H-DNA rich element in the human genome

# Introduction

This repository is part of the H-DNA paper `Ribosomal DNA arrays are the most H-DNA rich element in the human genome` by Georgakopoulos-Soares lab.

```
├── data
│   ├── centromeres
│   ├── chip_data
│   ├── DSBhotspots_hs1
│   ├── HDNA_data
│   ├── telomeres
│   └── transposable_elements
├── notebooks
└── scripts
```

- In the data directory, you will find relevant data, such as the H-DNA extractions in CHM13v2, transposable elements (Alu & SVA), DSB hotspots and chip seq data.
- In the scripts directory, you will find the pipeline we used to extract and systematically validate H-DNA motifs in CHM13v2 and ape genomes.
- In the notebooks directory, the basic data analysis code for the generated figures.

# Extraction

Mirror repeat is any sequence such that:

(x1...xN)s1s2...sM(xN...x1)

I used non-B gfa tool to extract mirror repeats. We required that any sequence to qualify as a mirror repeat, it's arm length must be at least 10 bp long
and it's corresponding spacer length at most 7 bp long.

Then, the H-DNA were extracted as a subset of mirror repeats abiding to the following constraints:

- AT Content less than 80%
- pyrimidine or pyrine rich (either at least 90%)

The non-B gfa tool is accessible on the following github link:

![non-B_gfa](https://github.com/abcsFrederick/non-B_gfa)

You will need to have this installed, in order to use the snakemake extraction pipeline.

# Analysis

The analysis was conducted in the following parts:

- Examine the overall density of mirror and H-DNA in human T2T satellite regions.
- Examine the rDNA regions and locate the H-DNA specific loci.
- Density in relationship to TSS/TES of protein coding and rRNA genes.
- Generalize the results to primate genomes.
- Examine the distribution of H-DNA motifs relative to PRMD9 alleles, for: predicted binding motifs (1), ChipSeq data (2), and PRMD9 Alleles.

# Updated Analysis

- EnrichmentRecombination contains the ChipSeq analysis on PRMD9B
- The Alu_SVA notebook contains the analysis with STR and G4 composition of H-DNA and the
distribution of H-DNA within transposable elements.
- Finally, the primatesHDNA notebook contains the updated figures on the diploid assemblies.

For any questions about the code, setup, extraction or the analysis, please send your email to

- Nikol Chantzi (nmc6088@psu.edu)
- Ilias Georgakopoulos-Soares (izg5139@psu.edu)

Thank you!
