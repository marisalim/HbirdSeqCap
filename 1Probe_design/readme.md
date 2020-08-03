## Probe design script descriptions

*Note: I wrote all of the python scripts with append functions to save files. As a result, old outputs must be deleted if the script is rerun. Otherwise, the outputs will end up with old + new + new2 + ... in the results file. Append does not write over the file. These were all written for Python2.7.*

### Calculate length of CDSs from Anna's hummingbird genome.

- `exonlength.py`

Output CSV file has `seqid`, `CDS length`, `Attributes` columns and plots CDS length histograms.

### Design probes for candidate genes:

*Here are the collection of scripts I wrote to target candidate genes. They were used roughly in the order listed below, but some steps had to be re-done to refine the filtering parameters so this is not a linear pipeline.*

- `findcandidates.py`: look for list of candidate gene sequences from Anna's hummingbird genome by searching for gene name
- `cand_blast.sh`: look for missing candidate genes from Anna's hummingbird genome with blastp search. These genes could not be found by searching for gene name in genome annotation file (gff). Use blast to identify the sequence region in the Anna's genome (trying to get a hummingbird ref vs. any other bird sequence from NCBI database).
- `calclen_transcriptomecands.py`: calculates the length of candidate gene sequences from transcriptome data in cands_in_transcriptomes.txt
- `viewlendistr_cands.py`: after finding candidate genes, calculate the sequence length distribution. Use this information to estimate the total target length for candidate genes (needed when ordering probes from Nimblegen).
- `addbuffer_cands.py`: the exons are too short, ~100bp or less. Use this script to add a buffer (we used 100bp) to the start of the gene.

### Design probes for random markers (putative neutral sequence):

*Here are the collection of scripts I wrote to target random marker genes. They were used roughly in the order listed below, but some steps had to be re-done to refine the filtering parameters so this is not a linear pipeline.*

- `GCcontent.py`: for each exon in Anna's hummingbird genome, calculate the GC content (#G's + #C's)/total exon length.
- `viewGCcont.py`: Visualize GC% as histograms. Want to keep exons with GC% 35-70%.
- `fastaforRepeatMask.py`: use RepeatMasker to remove sequences with too many repeats. Check if there are enough genes left after this filter. This script extracts fasta file to give to [RepeatMasker](http://www.repeatmasker.org/cgi-bin/WEBRepeatMasker). The output fasta file (forRMasker3.fa) is too large (>10Mb) for RepeatMasker, so split it into 2 files (forRMasker.fa and forRMasker2.fa)
  - RepeatMasker settings as follows:
```
- Browse to load .fa file
- Search engine: abblast
- Speed/sensitivity: default
- DNA source: Chicken
- Return format: tar file
- Return method: email (enter email)
- Note that by default, in advanced settings, interspersed and simple repeats are masked (low complexity reads)
- Submit sequence
```
- `afterRMask.py`: remove repeat sequences from the GCcontent csv file (above150bp_35-70GC.csv). This script requires bash edits before running.
- `testrm.py`: this was another check script. Because there are 327 duplicate sequence IDs (same ID, multiple repeat descriptions in the RepeatMasker output), it looked like there was a discrepancy in how many CDSs I should have. Use this script to investigate discrepancy.
- `lenfilter.py`: choose random markers. Checks length not >600bp, remove sequences with N's, choose longest CDS per gene, calculate total target length for random markers, check if buffer needed (add 50 or 100bp buffer as needed to meet target length goal (5-6Mb) - this script doesn't do the buffer adding, see description below for `addbuffer_randmarks.py`). This script requires bash edits before running.
- `viewrandommarkers.py`: visualize random marker CDS length, GC%, and N's distribution
- `addbuffer_randmarks.py`: need to add 100bp (or as many as I can depending on start position) to the CDSs. Note that the start position must be 1, not 0. If it is 0, the `randmarks_re_GCandNsfilt.py` script will read 0 as having no start and trim the sequence length to 0. You need to save .csv to .bed to .bed1 at top of script, outside of script.
- `randmarks_re_GCandNsfilt.py`: Because I added more sequence to the random markers, the GC% distribution changes and some sequences may now have N’s. So redo the GC% and N count + filter. You need to save .csv to .bed to .bed1 at top of script, outside of script.

## Get sequences for qPCR positive and negative control primers
Uses a list of positive and negative control gene sequence start/stop positions, the reference genome (`Calypte_anna.fa`), and `getseqs_forcontrols.py`

```
cp ./poscontrols_captarg_adjusted_startstop.txt ./poscontrols_captarg_adjusted_startstop.bed #change to be .bed file
cp ./negcontrols.txt ./negcontrols.bed #change to be .bed file
# Delete comments and headers from the bed files
python getseqs_forcontrols.py
```

- output files are poscontrols.fa and negcontrols.fa
- check that positive controls are IN probe design set, check that negative controls are NOT IN probe design set >> these were used to check exon capture experiment prior to sequencing
- used Primer3 to design primers for positive and negative controls
