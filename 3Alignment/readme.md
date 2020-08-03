# Alignment and exon capture evaluation script descriptions

## Alignment

*Map reads for each sample to respective species-specific reference fasta file using Ke Bi's [5-Alignment script](../CGRLScripts/5-Alignment).*

Required software:
- perl 5.24.0-threads
- samtools 1.3
- gatk 3.6
- picard 2.1.1
- novalign

I wrote batch wrappers `align_batch<#>_wrap.sh` to loop through Ke's script.

Alignment command:
```
perl 5-Alignment -f <path to species reference fasta> -r <path to cleaned reads to map> -o <output path> -G $GATK_HOME/GenomeAnalysisTK.jar -P $PICARD_HOME/picard.jar -c 0 -k 0 -m 1 -i 250 -v 100 -l 151 -t 180

```

Script flags:
- `-f` = reference fasta sequence file
- `-r` = path to cleaned read directory
- `-o` = path to results directory
- `-G` = path to GATK jar file
- `-P` = path to Picard jar file
- `-c` = only keep concordant mapping for PE reads? 1 = yes, 0 = no
- `-k` = keep original sam files? 1 = yes, 0 = no
- `-m` = method for alignment: 1) novalign, 2) bwa
- `-i` = if PE, average insert size for the libraries
- `-v` = if PE, STD insert size for the libraries (usually 0.1\*i)
- `-l` = read length
- `-t` = maximum alignment score acceptable for the best alignment

Wrapper scripts:
```
# batch scripts per species
sbatch align_batch2_wrap.sh #cviol
sbatch align_batch4_wrap.sh #ccoru batch1
sbatch align_batch5_wrap.sh #ccoru batch2
```

## Exon capture evaluation

*Calculate stats to evaluate how well the exon capture experiment worked with Ke Bi's [6-exonCaptureEvaluation script](../CGRLScripts/6-exonCaptureEvaluation).*

Required software:
- perl 5.18.4-threads
- samtools 1.3
- bedtools 2.25.0

I wrote wrapper scripts for each species, `exoncapeval_wrap_<species>.sh`, to loop through Ke's script.

Evaluation command:
```
perl 6-exonCaptureEvaluation Evaluation -genome <path to species reference fasta> -cleanDir <path to cleaned reads> -rawDir <path to raw reads> -bamDir <path to mapped reads> -resDir <output path> -bedFile <path to species reference bed file>
```

Script flags:
- `Evaluation` = calculate sensitivity and specificity
- `-genome` = reference genome
- `-cleanDir` = path to cleaned read directory
- `-rawDir` = path to raw read directory
- `-bamDir` = path to sorted bam files
- `-resDir` = path to results directory
- `-bedFile` = bed file that contains regions that probes are tiled at

Wrapper scripts:
```
sbatch exoncapeval_wrap_ccoru.sh
sbatch exoncapeval_wrap_cviol.sh
```
