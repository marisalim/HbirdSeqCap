# Read cleaning and assembly script descriptions

*Note: I wrote all of the python scripts with append functions to save files. As a result, old outputs must be deleted if the script is rerun. Otherwise, the outputs will end up with old + new + new2 + ... in the results file. Append does not write over the file. These were all written for Python2.7.*

## Clean raw reads

*This step uses CGRL's [1-ScrubReads_d script](../CGRLScripts/1-ScrubReads_d).*

Required software:
- perl 5.24.0-threads
- cutadapt
- super-deduper
- bowtie2 2.2.7
- flash 1.2.11
- trimmomatic 0.36
- fastqc 0.11.3

I wrote a wrapper script, `scrubreads_wrap.sh`, to loop through samples with CGRL's script. I allocated 48 hours on HPC, but the fastqc step didn't finish, so I moved cleaned fastq files to new folder and ran `runfastqc.sh`.

Read cleaning command:
```
perl 1-ScrubReads_d cleanPE -f <path to raw data> -o /ScrubReads_Results/ -t $TRIMMOMATIC_HOME/trimmomatic-0.36.jar -c e_coli_K12.fasta -d all -M 0 -l 150 -z
```

Script flags:
- `cleanPE` = reads to clean are paired end (vs. single end)
- `-f` = path to raw reads
- `-o` = output folder
- `-t` = path to trimmomatic executable file
- `-c` = file path to check for *E. coli* contamination
- `-d` = the particular library that you want to process. If 'all' is specified, all libraries in `-f` folder will be processed sequentially.
- `-M` = keep unmerged overlapping PE reads? 1 = yes (for DGE) and 0 = no (for SNP calling)
- `-l` = read length
- `-z` = if z is supplied, use fastQC to evaluate cleaned sequence reads

Defaults used for other flags:
- percent similarity to contaminant file to consider a read to be a contaminant = 0.9 (90%)
- Trimmomatic trimming length cutoff = 36
- quality trimming threshold = 23
- maximum allowed ratio between the number of mismatched base pairs and the overlap length = 0.1 (1 in 10 bp can be mismatched in overlap)
- will get rid of reads with any runs of bases longer than 0.5\*read length (too long, possible chimeras)

Wrappers:
```
# these were run on an HPC that used the slurm queue system
sbatch scrubreads_wrap.sh

sbatch runfastqc.sh
```

## Assemble reads

*The purpose of this step was to create species-specific reference assemblies for 6 samples each per species (with most data but not too many unmerged reads based on fastqc reports). All samples will later be mapped to these references. The divergence between study samples and the available hummingbird reference genome (Anna's hummmingbird) is too large, such that the % of reads mapping to the Anna's hummingbird genome was too low for downstream population genomic analysis. Mapping to the species-specific references was much higher.*

*This step uses CGRL's [2-GenerateAssembliesPhylo script](../CGRLScripts/2-GenerateAssembliesPhylo).*

Required software:
- perl 5.24.0-threads
- spade 3.8.1

I wrote a wrapper script, `spades_wrap.sh`. 

Assembly command:
```
perl 2-GenerateAssembliesPhylo spades -reads readsforspades -out readsforspades/spades_outs -np 1
```

Script flags:
- `spades` = use spades for assembly
- `-reads` = input fastq files for spades (merged '_1.fq', '_2.fq', and '_u.fq' files from the scrub reads script output)
- `-out` = output folder
- `-np` = number of processors to use for assembly

## Generate species-specific reference

*A consequence of making separate species-specific references is that the target gene coordinates differ for each assembly. Consequently, to make the species-specific references we had to blast search for all the targeted genes in the assemblies for each species. At the end of this process, we will have 1 reference file per species. This will be achieved by concatenating across the 6 spades assemblies per species using CGRL's [3-FindingTargetsV9 script](../CGRLScripts/3-FindingTargetsV9).*

1. Extract candidate target gene sequences: `python Canna_candgene_ref.py`
2. Extract random marker target gene sequences: `python Canna_randmarks_ref.py`
3. Combine target sequences designed from Anna's hummingbird (outputs of steps 1 and 2) with target sequences designed from [transcriptomes](https://github.com/marisalim/Transcriptome_pipeline) into single fasta file. This will be used as `-t` input for step 4.
4. Concatenate assemblies and blast search: `sbatch findingtargs_wrap.sh`

Required software:
- perl 5.24.0-threads
- blat v35
- cap3
- cd-hit
- exonerate
- blast 2.2.31

Reference generation command:
```
perl 3-FindingTargetsV9 combineExon -t targeted_loci.fasta -a raw_assembly_dir -p 80 -b 1 -e 4 -f 500
```

Script flags:
- `combineExon` = reconstruct reference contigs by joining individual exons if they are from the same gene
- `-t` = target sequence files (combined results from `Canna_candgene_ref.py` and `Canna_randmarks_ref.py` into single file)
- `-a` = a folder with all final assemblies
- `-p` = how similar does the assembled contig have to be to the target (out of 100)?
- `-b` = merging individual assemblies? 1 = yes (for pop gen datasets) and 0 = no (for phylogenetic datasets)
- `-e` = target sequences could be one of the following: 1) individual exons, 2) cDNA, 3) transcripts, or 4) random (such as UCEs, no need for exon identification)
- `-f` = flanking bp to add

Defaults used for other flags:
- how much should you cluster the targets and assemblies at the get go = 0.98
- how much overlap is allowed between adjoining assembled contigs mapping to the same target = 0.3
- used in the initial BLAST step = 1e-10










