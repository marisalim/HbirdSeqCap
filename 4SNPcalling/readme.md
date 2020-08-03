# SNP calling script descriptions

*Python scripts written for Python2.7.*

## Call SNPs

Required software
- perl 5.18.4-threads
- samtools 1.3.1
- bcftools 1.3.1

SNP call command:
```
samtools mpileup -s -Q 20 -A -t DP,SP -d 100000 -E -ugf <path to species reference fasta> <path to mapped reads bam files> | bcftools call -c -p 0.1 - > <output raw vcf file>
```

Wrapper scripts:
```
sbatch snprawfiltering_ccoru.sh
sbatch snprawfiltering_cviol.sh
```

## Filter raw vcf files

*Filter SNPs with Tyler Linderoth's [SNPcleaner225.pl script](../CGRLScripts/SNPcleaner225.pl).*

Required software:
- perl 5.18.4-threads

Filter command:
```
perl SNPcleaner225.pl -u 3 -k 27 -d 80 -a 0 -H 0.0001 -h 0 -B <path to bed file of sites to keep> -p <path to failed sites file> -v <path to raw vcf file>

```

Script flags:
- `-u` = minimum individual coverage threshold used in -k
- `-k` = min number of individuals with less than -u coverage
- `-d` = minimum site read depth
- `-a` = minimum number of alternate alleles for site
- `-H` = min p-value for exact test of excess of heterozygotes
- `-B` = name of BED file for kept SNP positions
- `-p` = name of file to dump sites that failed filters
- `-v` = process nonvariants

Wrapper scripts:
```
# run from Novalign_outs/ (where the vcf files are)
sbatch cleansnp_ccoru.sh
sbatch cleansnp_cviol.sh
```

After you get `sites_to_keep.bed` (in Novalign_outs directory), format file:

`cut -f1,2 t sites_to_keep.bed > sites_to_keep.keep`

## Remove sites

*The SNPcleaner step outputs a list of sites to remove due to poor quality. Go to the bed folder and find a file called 'combined_sites_to_remove.txt' (rename it with species abbreviation and copy to Novalign_outs). Remove the sites listed in the 'combined_sites_to_remove.txt' file from 'sites_to_keep.keep' file if they are there.*

1. delete header on the sites to remove file
2. concatenate the columns to compare between files (do this for sites_to_remove and sites_to_keep files) with `removesites.py`:

```
# per species
python ccoru_removesites.py
python cviol_removesites.py
```

Wrapper scripts for above python scripts:
```
# per species
sbatch removesites_wrap_ccoru.sh
sbatch removesites_wrap_cviol.sh
```

3. compare files
4. for output file, undo the concatenation and separate chr and site again for output file:
5. Create a file with a list of all the contigs in the keep list:

`cut -f1 speciesA.keep | sort -V | uniq | awk '{print $1":"}' > speciesA.rf`

The resulting 'speciesA.keep' and 'speciesA.rf' files will be used for ANGSD pipeline population genomic analyses!
