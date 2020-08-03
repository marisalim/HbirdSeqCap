#!/bin/bash 


#this took 5 hours to complete
#for i in {1..5}; do LFMM -x Ccoru_genotransposed.txt -v Ccoru_standardized_elevation.txt -K 2 -o run_K2_i5000b500_run$i -D 1 -i 5000 -b 500; done

#could not finish even 1 iteration with these -i and -b settings after 24 hours (timed out)
#try 2 K = 1 and 48 hours - timed out
#for i in {1..2}; do LFMM -x Ccoru_genotransposed.txt -v Ccoru_standardized_elevation.txt -K 1 -o run_K1_i500000b50000_run$i -D 1 -i 500000 -b 50000; done

#try again on Ke's computer
#LFMM -x Ccoru_genotransposed.txt -v Ccoru_standardized_elevation.txt -K 1 -o run_K1_i500000b50000_run1 -D 1 -i 500000 -b 50000

#try with fewer -i and -b
LFMM -x Ccoru_genotransposed.txt -v Ccoru_standardized_elevation.txt -K 1 -o run_K1_i250000b25000_run1 -D 1 -i 250000 -b 25000
