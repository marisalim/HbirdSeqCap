The code in this repository accompany the *in revision* paper:

*Pervasive genomic signatures of local adaptation to altitude across highland specialist Andean hummingbird populations*

Marisa C.W. Lim, Ke Bi, Christopher C. Witt, Catherine H. Graham, Liliana M. Davalos

## Variant calling pipeline
Note: this pipeline uses scripts written by [UC Berkeley MVZ/CGRL scientists](./CGRLScripts).

[1. Exon capture probe design](./1Probe_design): There were two probe design procedures - the first to capture previously identified candidate genes for high-altitude adaptation and the second to capture a 'random' set of exons across the entire hummingbird genome using the Anna's hummingbird (*Calypte anna*) as the reference.

[2. Clean raw_reads_and assemble](./2Clean_reads_assemble): Remove low quality reads, adapter sequences, PCR duplicates, and potential bacterial contaminant sequences; merge paired-end reads; and generate *de novo* assemblies as species-specific references

[3. Read alignment](./3Alignment): Map reads to species-specific reference assemblies and evaluate capture experiment (% reads retained after filtering but before alignment, length of mapped data, % of reads aligned to target region (specificity), % of targeted regions covered by at least one read (sensitivity), average coverage, variationin coverage, and % of sites retained at multiple coverage depths)

[4. SNP calling](./4SNPcalling): Call variants and filter loci with too much missing data, within 10bp of indels, that are not biallelic, and/or that have excessive heterozygosity

## Population genomic analyses

- [PCA](./5Pop_genomics_analyses/A_PCA)

- [Admixture](./5Pop_genomics_analyses/B_ngsAdmix)

- [Relatedness coefficient](./5Pop_genomics_analyses/C_NgsRelate)

- Gene flow estimates

    - Between population diversity ([Fst](./5Pop_genomics_analyses/D_Fstcalc), [Dxy](./5Pop_genomics_analyses/E_dxycalc))

    - [Isolation by distance - geodesic and least cost distance](./5Pop_genomics_analyses/F_IBD)


- [Within population diversity](./5Pop_genomics_analyses/G_WattersonsTheta)


- [Estimating Effective Migration Surfaces (EEMS)](./5Pop_genomics_analyses/H_EEMS)


- Test for natural selection with Latent Factor Mixed Models (LFMM)

    - [Call genotypes](./5Pop_genomics_analyses/I_ANGSD_genotype_calls)

    - [LFMM](./5Pop_genomics_analyses/J_LFMM)


- Shiny app

*work in progress*

Check out additional plots of the data:
- PCA - more axes
- SNP clines
- Candidate gene pathways
