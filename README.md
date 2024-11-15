# H-DNA & Mirror Repeats in Human Genome


# Extraction

Mirror repeat is any sequence such that:

(x1...xN)s1s2...sM(xN...x1)

I used non-B gfa tool to extract mirror repeats. We required that any sequence to qualify as a mirror repeat, it's arm length must be at least 10 bp long
and it's corresponding spacer length at most 7 bp long.

Then, the H-DNA were extracted as a subset of mirror repeats abiding to the following constraints:

- AT Content less than 80%
- pyrimidine or pyrine rich (either at least 90%)

# Analysis

The EDA was conducted in the following parts:

- Examine the overall density of mirror and H-DNA in human T2T satellite regions.
- Examine the rDNA regions and locate the H-DNA specific loci.
- Generalize the results to primate genomes.


# Updated Analysis

A new directory has been added with the new notebooks under the

```
notebooks/update
```

- EnrichmentRecombination contains the ChipSeq analysis on PRMD9B
- The Alu_SVA notebook contains the analysis with STR and G4 composition of H-DNA and the
distribution of H-DNA within transposable elements.
- Finally, the primatesHDNA notebook contains the updated figures on the diploid assemblies.


For any questions about the extraction or the analysis, please send your email to

- Nikol Chantzi (nmc6088@psu.edu)
- Ilias Georgakopoulos-Soares (izg5139@psu.edu)


