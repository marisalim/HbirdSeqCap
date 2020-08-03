#!/bin/bash 
#SBATCH -N 1 
#SBATCH -p LM 
#SBATCH --mem=128GB 
#SBATCH -t 24:00:00 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x 

source ~/.bashrc

# pop 1 and 2
realSFS cviol_forFstpop1_out.saf.idx cviol_forFstpop2_out.saf.idx > cviol_pop12.ml
realSFS fst index cviol_forFstpop1_out.saf.idx cviol_forFstpop2_out.saf.idx -sfs cviol_pop12.ml -fstout cviol_pop12_fstout
realSFS fst stats cviol_pop12_fstout.fst.idx > cviol_pop12_fstout.txt

# pop 1 and 3
realSFS cviol_forFstpop1_out.saf.idx cviol_forFstpop3_out.saf.idx > cviol_pop13.ml
realSFS fst index cviol_forFstpop1_out.saf.idx cviol_forFstpop3_out.saf.idx -sfs cviol_pop13.ml -fstout cviol_pop13_fstout
realSFS fst stats cviol_pop13_fstout.fst.idx > cviol_pop13_fstout.txt

# pop 1 and 4
realSFS cviol_forFstpop1_out.saf.idx cviol_forFstpop4_out.saf.idx > cviol_pop14.ml
realSFS fst index cviol_forFstpop1_out.saf.idx cviol_forFstpop4_out.saf.idx -sfs cviol_pop14.ml -fstout cviol_pop14_fstout
realSFS fst stats cviol_pop14_fstout.fst.idx > cviol_pop14_fstout.txt

# pop 1 and 5
realSFS cviol_forFstpop1_out.saf.idx cviol_forFstpop5_out.saf.idx > cviol_pop15.ml
realSFS fst index cviol_forFstpop1_out.saf.idx cviol_forFstpop5_out.saf.idx -sfs cviol_pop15.ml -fstout cviol_pop15_fstout
realSFS fst stats cviol_pop15_fstout.fst.idx > cviol_pop15_fstout.txt

# pop 1 and 6
realSFS cviol_forFstpop1_out.saf.idx cviol_forFstpop6_out.saf.idx > cviol_pop16.ml
realSFS fst index cviol_forFstpop1_out.saf.idx cviol_forFstpop6_out.saf.idx -sfs cviol_pop16.ml -fstout cviol_pop16_fstout
realSFS fst stats cviol_pop16_fstout.fst.idx > cviol_pop16_fstout.txt

# pop 2 and 3
realSFS cviol_forFstpop2_out.saf.idx cviol_forFstpop3_out.saf.idx > cviol_pop23.ml
realSFS fst index cviol_forFstpop2_out.saf.idx cviol_forFstpop3_out.saf.idx -sfs cviol_pop23.ml -fstout cviol_pop23_fstout
realSFS fst stats cviol_pop23_fstout.fst.idx > cviol_pop23_fstout.txt

# pop 2 and 4
realSFS cviol_forFstpop2_out.saf.idx cviol_forFstpop4_out.saf.idx > cviol_pop24.ml
realSFS fst index cviol_forFstpop2_out.saf.idx cviol_forFstpop4_out.saf.idx -sfs cviol_pop24.ml -fstout cviol_pop24_fstout
realSFS fst stats cviol_pop24_fstout.fst.idx > cviol_pop24_fstout.txt

# pop 2 and 5
realSFS cviol_forFstpop2_out.saf.idx cviol_forFstpop5_out.saf.idx > cviol_pop25.ml
realSFS fst index cviol_forFstpop2_out.saf.idx cviol_forFstpop5_out.saf.idx -sfs cviol_pop25.ml -fstout cviol_pop25_fstout
realSFS fst stats cviol_pop25_fstout.fst.idx > cviol_pop25_fstout.txt

# pop 2 and 6
realSFS cviol_forFstpop2_out.saf.idx cviol_forFstpop6_out.saf.idx > cviol_pop26.ml
realSFS fst index cviol_forFstpop2_out.saf.idx cviol_forFstpop6_out.saf.idx -sfs cviol_pop26.ml -fstout cviol_pop26_fstout
realSFS fst stats cviol_pop26_fstout.fst.idx > cviol_pop26_fstout.txt

# pop 3 and 4
realSFS cviol_forFstpop3_out.saf.idx cviol_forFstpop4_out.saf.idx > cviol_pop34.ml
realSFS fst index cviol_forFstpop3_out.saf.idx cviol_forFstpop4_out.saf.idx -sfs cviol_pop34.ml -fstout cviol_pop34_fstout
realSFS fst stats cviol_pop34_fstout.fst.idx > cviol_pop34_fstout.txt

# pop 3 and 5
realSFS cviol_forFstpop3_out.saf.idx cviol_forFstpop5_out.saf.idx > cviol_pop35.ml
realSFS fst index cviol_forFstpop3_out.saf.idx cviol_forFstpop5_out.saf.idx -sfs cviol_pop35.ml -fstout cviol_pop35_fstout
realSFS fst stats cviol_pop35_fstout.fst.idx > cviol_pop35_fstout.txt

# pop 3 and 6
realSFS cviol_forFstpop3_out.saf.idx cviol_forFstpop6_out.saf.idx > cviol_pop36.ml
realSFS fst index cviol_forFstpop3_out.saf.idx cviol_forFstpop6_out.saf.idx -sfs cviol_pop36.ml -fstout cviol_pop36_fstout
realSFS fst stats cviol_pop36_fstout.fst.idx > cviol_pop36_fstout.txt

# pop 4 and 5
realSFS cviol_forFstpop4_out.saf.idx cviol_forFstpop5_out.saf.idx > cviol_pop45.ml
realSFS fst index cviol_forFstpop4_out.saf.idx cviol_forFstpop5_out.saf.idx -sfs cviol_pop45.ml -fstout cviol_pop45_fstout
realSFS fst stats cviol_pop45_fstout.fst.idx > cviol_pop45_fstout.txt

# pop 4 and 6
realSFS cviol_forFstpop4_out.saf.idx cviol_forFstpop6_out.saf.idx > cviol_pop46.ml
realSFS fst index cviol_forFstpop4_out.saf.idx cviol_forFstpop6_out.saf.idx -sfs cviol_pop46.ml -fstout cviol_pop46_fstout
realSFS fst stats cviol_pop46_fstout.fst.idx > cviol_pop46_fstout.txt

# pop 5 and 6
realSFS cviol_forFstpop5_out.saf.idx cviol_forFstpop6_out.saf.idx > cviol_pop56.ml
realSFS fst index cviol_forFstpop5_out.saf.idx cviol_forFstpop6_out.saf.idx -sfs cviol_pop56.ml -fstout cviol_pop56_fstout
realSFS fst stats cviol_pop56_fstout.fst.idx > cviol_pop56_fstout.txt





