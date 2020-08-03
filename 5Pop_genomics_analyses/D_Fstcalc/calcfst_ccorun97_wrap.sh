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
realSFS ccoru_forFstpop1_out.saf.idx ccoru_forFstpop2_out.saf.idx > ccoru_pop12.ml
realSFS fst index ccoru_forFstpop1_out.saf.idx ccoru_forFstpop2_out.saf.idx -sfs ccoru_pop12.ml -fstout ccoru_pop12_fstout
realSFS fst stats ccoru_pop12_fstout.fst.idx > ccoru_pop12_fstout.txt

# pop 1 and 3
realSFS ccoru_forFstpop1_out.saf.idx ccoru_forFstpop3_out.saf.idx > ccoru_pop13.ml
realSFS fst index ccoru_forFstpop1_out.saf.idx ccoru_forFstpop3_out.saf.idx -sfs ccoru_pop13.ml -fstout ccoru_pop13_fstout
realSFS fst stats ccoru_pop13_fstout.fst.idx > ccoru_pop13_fstout.txt

# pop 1 and 4
realSFS ccoru_forFstpop1_out.saf.idx ccoru_forFstpop4_out.saf.idx > ccoru_pop14.ml
realSFS fst index ccoru_forFstpop1_out.saf.idx ccoru_forFstpop4_out.saf.idx -sfs ccoru_pop14.ml -fstout ccoru_pop14_fstout
realSFS fst stats ccoru_pop14_fstout.fst.idx > ccoru_pop14_fstout.txt

# pop 1 and 5
realSFS ccoru_forFstpop1_out.saf.idx ccoru_forFstpop5_out.saf.idx > ccoru_pop15.ml
realSFS fst index ccoru_forFstpop1_out.saf.idx ccoru_forFstpop5_out.saf.idx -sfs ccoru_pop15.ml -fstout ccoru_pop15_fstout
realSFS fst stats ccoru_pop15_fstout.fst.idx > ccoru_pop15_fstout.txt

# pop 1 and 6
realSFS ccoru_forFstpop1_out.saf.idx ccoru_forFstpop6_out.saf.idx > ccoru_pop16.ml
realSFS fst index ccoru_forFstpop1_out.saf.idx ccoru_forFstpop6_out.saf.idx -sfs ccoru_pop16.ml -fstout ccoru_pop16_fstout
realSFS fst stats ccoru_pop16_fstout.fst.idx > ccoru_pop16_fstout.txt

# pop 1 and 7
realSFS ccoru_forFstpop1_out.saf.idx ccoru_forFstpop7_out.saf.idx > ccoru_pop17.ml
realSFS fst index ccoru_forFstpop1_out.saf.idx ccoru_forFstpop7_out.saf.idx -sfs ccoru_pop17.ml -fstout ccoru_pop17_fstout
realSFS fst stats ccoru_pop17_fstout.fst.idx > ccoru_pop17_fstout.txt

# pop 1 and 8
realSFS ccoru_forFstpop1_out.saf.idx ccoru_forFstpop8_out.saf.idx > ccoru_pop18.ml
realSFS fst index ccoru_forFstpop1_out.saf.idx ccoru_forFstpop8_out.saf.idx -sfs ccoru_pop18.ml -fstout ccoru_pop18_fstout
realSFS fst stats ccoru_pop18_fstout.fst.idx > ccoru_pop18_fstout.txt

# pop 1 and 9
realSFS ccoru_forFstpop1_out.saf.idx ccoru_forFstpop9_out.saf.idx > ccoru_pop19.ml
realSFS fst index ccoru_forFstpop1_out.saf.idx ccoru_forFstpop9_out.saf.idx -sfs ccoru_pop19.ml -fstout ccoru_pop19_fstout
realSFS fst stats ccoru_pop19_fstout.fst.idx > ccoru_pop19_fstout.txt

# pop 2 and 3
realSFS ccoru_forFstpop2_out.saf.idx ccoru_forFstpop3_out.saf.idx > ccoru_pop23.ml
realSFS fst index ccoru_forFstpop2_out.saf.idx ccoru_forFstpop3_out.saf.idx -sfs ccoru_pop23.ml -fstout ccoru_pop23_fstout
realSFS fst stats ccoru_pop23_fstout.fst.idx > ccoru_pop23_fstout.txt

# pop 2 and 4
realSFS ccoru_forFstpop2_out.saf.idx ccoru_forFstpop4_out.saf.idx > ccoru_pop24.ml
realSFS fst index ccoru_forFstpop2_out.saf.idx ccoru_forFstpop4_out.saf.idx -sfs ccoru_pop24.ml -fstout ccoru_pop24_fstout
realSFS fst stats ccoru_pop24_fstout.fst.idx > ccoru_pop24_fstout.txt

# pop 2 and 5
realSFS ccoru_forFstpop2_out.saf.idx ccoru_forFstpop5_out.saf.idx > ccoru_pop25.ml
realSFS fst index ccoru_forFstpop2_out.saf.idx ccoru_forFstpop5_out.saf.idx -sfs ccoru_pop25.ml -fstout ccoru_pop25_fstout
realSFS fst stats ccoru_pop25_fstout.fst.idx > ccoru_pop25_fstout.txt

# pop 2 and 6
realSFS ccoru_forFstpop2_out.saf.idx ccoru_forFstpop6_out.saf.idx > ccoru_pop26.ml
realSFS fst index ccoru_forFstpop2_out.saf.idx ccoru_forFstpop6_out.saf.idx -sfs ccoru_pop26.ml -fstout ccoru_pop26_fstout
realSFS fst stats ccoru_pop26_fstout.fst.idx > ccoru_pop26_fstout.txt

# pop 2 and 7
realSFS ccoru_forFstpop2_out.saf.idx ccoru_forFstpop7_out.saf.idx > ccoru_pop27.ml
realSFS fst index ccoru_forFstpop2_out.saf.idx ccoru_forFstpop7_out.saf.idx -sfs ccoru_pop27.ml -fstout ccoru_pop27_fstout
realSFS fst stats ccoru_pop27_fstout.fst.idx > ccoru_pop27_fstout.txt

# pop 2 and 8
realSFS ccoru_forFstpop2_out.saf.idx ccoru_forFstpop8_out.saf.idx > ccoru_pop28.ml
realSFS fst index ccoru_forFstpop2_out.saf.idx ccoru_forFstpop8_out.saf.idx -sfs ccoru_pop28.ml -fstout ccoru_pop28_fstout
realSFS fst stats ccoru_pop28_fstout.fst.idx > ccoru_pop28_fstout.txt

# pop 2 and 9
realSFS ccoru_forFstpop2_out.saf.idx ccoru_forFstpop9_out.saf.idx > ccoru_pop29.ml
realSFS fst index ccoru_forFstpop2_out.saf.idx ccoru_forFstpop9_out.saf.idx -sfs ccoru_pop29.ml -fstout ccoru_pop29_fstout
realSFS fst stats ccoru_pop29_fstout.fst.idx > ccoru_pop29_fstout.txt

# pop 3 and 4
realSFS ccoru_forFstpop3_out.saf.idx ccoru_forFstpop4_out.saf.idx > ccoru_pop34.ml
realSFS fst index ccoru_forFstpop3_out.saf.idx ccoru_forFstpop4_out.saf.idx -sfs ccoru_pop34.ml -fstout ccoru_pop34_fstout
realSFS fst stats ccoru_pop34_fstout.fst.idx > ccoru_pop34_fstout.txt

# pop 3 and 5
realSFS ccoru_forFstpop3_out.saf.idx ccoru_forFstpop5_out.saf.idx > ccoru_pop35.ml
realSFS fst index ccoru_forFstpop3_out.saf.idx ccoru_forFstpop5_out.saf.idx -sfs ccoru_pop35.ml -fstout ccoru_pop35_fstout
realSFS fst stats ccoru_pop35_fstout.fst.idx > ccoru_pop35_fstout.txt

# pop 3 and 6
realSFS ccoru_forFstpop3_out.saf.idx ccoru_forFstpop6_out.saf.idx > ccoru_pop36.ml
realSFS fst index ccoru_forFstpop3_out.saf.idx ccoru_forFstpop6_out.saf.idx -sfs ccoru_pop36.ml -fstout ccoru_pop36_fstout
realSFS fst stats ccoru_pop36_fstout.fst.idx > ccoru_pop36_fstout.txt

# pop 3 and 7
realSFS ccoru_forFstpop3_out.saf.idx ccoru_forFstpop7_out.saf.idx > ccoru_pop37.ml
realSFS fst index ccoru_forFstpop3_out.saf.idx ccoru_forFstpop7_out.saf.idx -sfs ccoru_pop37.ml -fstout ccoru_pop37_fstout
realSFS fst stats ccoru_pop37_fstout.fst.idx > ccoru_pop37_fstout.txt

# pop 3 and 8
realSFS ccoru_forFstpop3_out.saf.idx ccoru_forFstpop8_out.saf.idx > ccoru_pop38.ml
realSFS fst index ccoru_forFstpop3_out.saf.idx ccoru_forFstpop8_out.saf.idx -sfs ccoru_pop38.ml -fstout ccoru_pop38_fstout
realSFS fst stats ccoru_pop38_fstout.fst.idx > ccoru_pop38_fstout.txt

# pop 3 and 9
realSFS ccoru_forFstpop3_out.saf.idx ccoru_forFstpop9_out.saf.idx > ccoru_pop39.ml
realSFS fst index ccoru_forFstpop3_out.saf.idx ccoru_forFstpop9_out.saf.idx -sfs ccoru_pop39.ml -fstout ccoru_pop39_fstout
realSFS fst stats ccoru_pop39_fstout.fst.idx > ccoru_pop39_fstout.txt

# pop 4 and 5
realSFS ccoru_forFstpop4_out.saf.idx ccoru_forFstpop5_out.saf.idx > ccoru_pop45.ml
realSFS fst index ccoru_forFstpop4_out.saf.idx ccoru_forFstpop5_out.saf.idx -sfs ccoru_pop45.ml -fstout ccoru_pop45_fstout
realSFS fst stats ccoru_pop45_fstout.fst.idx > ccoru_pop45_fstout.txt

# pop 4 and 6
realSFS ccoru_forFstpop4_out.saf.idx ccoru_forFstpop6_out.saf.idx > ccoru_pop46.ml
realSFS fst index ccoru_forFstpop4_out.saf.idx ccoru_forFstpop6_out.saf.idx -sfs ccoru_pop46.ml -fstout ccoru_pop46_fstout
realSFS fst stats ccoru_pop46_fstout.fst.idx > ccoru_pop46_fstout.txt

# pop 4 and 7
realSFS ccoru_forFstpop4_out.saf.idx ccoru_forFstpop7_out.saf.idx > ccoru_pop47.ml
realSFS fst index ccoru_forFstpop4_out.saf.idx ccoru_forFstpop7_out.saf.idx -sfs ccoru_pop47.ml -fstout ccoru_pop47_fstout
realSFS fst stats ccoru_pop47_fstout.fst.idx > ccoru_pop47_fstout.txt

# pop 4 and 8
realSFS ccoru_forFstpop4_out.saf.idx ccoru_forFstpop8_out.saf.idx > ccoru_pop48.ml
realSFS fst index ccoru_forFstpop4_out.saf.idx ccoru_forFstpop8_out.saf.idx -sfs ccoru_pop48.ml -fstout ccoru_pop48_fstout
realSFS fst stats ccoru_pop48_fstout.fst.idx > ccoru_pop48_fstout.txt

# pop 4 and 9
realSFS ccoru_forFstpop4_out.saf.idx ccoru_forFstpop9_out.saf.idx > ccoru_pop49.ml
realSFS fst index ccoru_forFstpop4_out.saf.idx ccoru_forFstpop9_out.saf.idx -sfs ccoru_pop49.ml -fstout ccoru_pop49_fstout
realSFS fst stats ccoru_pop49_fstout.fst.idx > ccoru_pop49_fstout.txt

# pop 5 and 6
realSFS ccoru_forFstpop5_out.saf.idx ccoru_forFstpop6_out.saf.idx > ccoru_pop56.ml
realSFS fst index ccoru_forFstpop5_out.saf.idx ccoru_forFstpop6_out.saf.idx -sfs ccoru_pop56.ml -fstout ccoru_pop56_fstout
realSFS fst stats ccoru_pop56_fstout.fst.idx > ccoru_pop56_fstout.txt

# pop 5 and 7
realSFS ccoru_forFstpop5_out.saf.idx ccoru_forFstpop7_out.saf.idx > ccoru_pop57.ml
realSFS fst index ccoru_forFstpop5_out.saf.idx ccoru_forFstpop7_out.saf.idx -sfs ccoru_pop57.ml -fstout ccoru_pop57_fstout
realSFS fst stats ccoru_pop57_fstout.fst.idx > ccoru_pop57_fstout.txt

# pop 5 and 8
realSFS ccoru_forFstpop5_out.saf.idx ccoru_forFstpop8_out.saf.idx > ccoru_pop58.ml
realSFS fst index ccoru_forFstpop5_out.saf.idx ccoru_forFstpop8_out.saf.idx -sfs ccoru_pop58.ml -fstout ccoru_pop58_fstout
realSFS fst stats ccoru_pop58_fstout.fst.idx > ccoru_pop58_fstout.txt

# pop 5 and 9
realSFS ccoru_forFstpop5_out.saf.idx ccoru_forFstpop9_out.saf.idx > ccoru_pop59.ml
realSFS fst index ccoru_forFstpop5_out.saf.idx ccoru_forFstpop9_out.saf.idx -sfs ccoru_pop59.ml -fstout ccoru_pop59_fstout
realSFS fst stats ccoru_pop59_fstout.fst.idx > ccoru_pop59_fstout.txt

# pop 6 and 7
realSFS ccoru_forFstpop6_out.saf.idx ccoru_forFstpop7_out.saf.idx > ccoru_pop67.ml
realSFS fst index ccoru_forFstpop6_out.saf.idx ccoru_forFstpop7_out.saf.idx -sfs ccoru_pop67.ml -fstout ccoru_pop67_fstout
realSFS fst stats ccoru_pop67_fstout.fst.idx > ccoru_pop67_fstout.txt

# pop 6 and 8
realSFS ccoru_forFstpop6_out.saf.idx ccoru_forFstpop8_out.saf.idx > ccoru_pop68.ml
realSFS fst index ccoru_forFstpop6_out.saf.idx ccoru_forFstpop8_out.saf.idx -sfs ccoru_pop68.ml -fstout ccoru_pop68_fstout
realSFS fst stats ccoru_pop68_fstout.fst.idx > ccoru_pop68_fstout.txt

# pop 6 and 9
realSFS ccoru_forFstpop6_out.saf.idx ccoru_forFstpop9_out.saf.idx > ccoru_pop69.ml
realSFS fst index ccoru_forFstpop6_out.saf.idx ccoru_forFstpop9_out.saf.idx -sfs ccoru_pop69.ml -fstout ccoru_pop69_fstout
realSFS fst stats ccoru_pop69_fstout.fst.idx > ccoru_pop69_fstout.txt

# pop 7 and 8
realSFS ccoru_forFstpop7_out.saf.idx ccoru_forFstpop8_out.saf.idx > ccoru_pop78.ml
realSFS fst index ccoru_forFstpop7_out.saf.idx ccoru_forFstpop8_out.saf.idx -sfs ccoru_pop78.ml -fstout ccoru_pop78_fstout
realSFS fst stats ccoru_pop78_fstout.fst.idx > ccoru_pop78_fstout.txt

# pop 7 and 9
realSFS ccoru_forFstpop7_out.saf.idx ccoru_forFstpop9_out.saf.idx > ccoru_pop79.ml
realSFS fst index ccoru_forFstpop7_out.saf.idx ccoru_forFstpop9_out.saf.idx -sfs ccoru_pop79.ml -fstout ccoru_pop79_fstout
realSFS fst stats ccoru_pop79_fstout.fst.idx > ccoru_pop79_fstout.txt

# pop 8 and 9
realSFS ccoru_forFstpop8_out.saf.idx ccoru_forFstpop9_out.saf.idx > ccoru_pop89.ml
realSFS fst index ccoru_forFstpop8_out.saf.idx ccoru_forFstpop9_out.saf.idx -sfs ccoru_pop89.ml -fstout ccoru_pop89_fstout
realSFS fst stats ccoru_pop89_fstout.fst.idx > ccoru_pop89_fstout.txt




