## This script is used to establish quantitative trait architectures and perform artificial selection experiments on the neutral populations generated using the Burnin.sh script.
## IMPORTANT: Change the output directory and in the first step and the variables in the second step before running.
## Run the script on server using: nohup bash Selection.sh > Selection.nohup &



##### Set variables here:

QTL_in=NQTLNQTL_REPLACE
SLN_NM=SLN_NM_REPLACE

# still need to change variables below to match.


## Step 1: Ceate directories to store the outputs. Make sure to change the directory names in the first three lines.

cd /users/r/b/rbrennan/evolve-resequence-simulation/simulations # Change this to a directory where you want to store all you simulation outputs.
mkdir "$QTL_in" # Change this to what you want to name this particular quantitative trait architectures and/or the experimental design that you are simulating.
cd "$QTL_in" # Same as above
for k in {1..100} # Number of simulation replicates that you want to create.
do
    mkdir 'SimRep'$k
    cd 'SimRep'$k

    # SET degree of selection! this must match below
    mkdir "$SLN_NM"

    cd ..
done
cd ..

## Step 2: Run burnin using SLiM 2. Variables inside the loop are all customizable and can be changed as desired.
## IMPORTANT: Number of generations has to be changed in the slim script.
## It is recommended if you copy this shell script for each of the trait architecture x experimental combinations that you want to test and make changes in the new script.
## For example, when simulating the scenario with 100 QTLs, copy this file and rename it as NQTL100.sh. Set the NQTL variable to be 100 and run it with "nohup bash NQTL100.sh > NQTL100.nohup &"
## Note: The number of generations in the selection experiment can only be changed in the Selection.slim script.
## Depending on the population size, more or fewer burnin generations might have been needed. If that's the case, edit the Selection.slim file as instructed in the file.

for k in {1..100} # Set the number of simulation replicates that you want to create.
do
    echo $k
            # Set the path to the SLiM program in the next line
            ~/bin/SLiM/slim \
            -d SimRepID=$k  \
            -d "BurninPath='/users/r/b/rbrennan/evolve-resequence-simulation/simulations/Burnin/'" \
            -d "BurninFilename='Burnin.txt'" \
            -d LCh=30000000 \
            -d RecRate=1e-8 \
            -d SampleSize=50 \
            -d NQTL=NQTL_REPLACE \
            -d ESMean=1.0 \
            -d "ESDist='f'" \
            -d LowFreq=F \
            -d FreqBound=0.05 \
            -d LowerPosBound=0 \
            -d UpperPosBound=29999999 \
            -d D=0.5 \
            -d OffspringSize=80000 \
            -d ParentSize=50 \
            -d Cutoff_Quantile=QUANTILE_REPLACE \
            -d "Quantile_label='SLN_NM_REPLACE'" \
            -d "OutPath='/users/r/b/rbrennan/evolve-resequence-simulation/simulations/NQTL20/'" \
            /users/r/b/rbrennan/evolve-resequence-simulation/SlimScripts/Selection.slim # Directory to the Selection.slim file included in the simulation tool.

        # reformat to save disc space
        cd ~/evolve-resequence-simulation/simulations/$QTL_in/SimRep$k/$SLN_NM/

        for i in $(ls | grep 'genomes' | cut -f 1-2 -d "_"); do

            # save num of genomes remaining in each:
            head -n 1 ${i}_genomes.txt | cut -f 4 -d " " | paste - >> number_genomes.txt

            sed -n '/^Mutations/,/^Genomes/p' ${i}_genomes.txt | grep -v "^Genomes:" | grep -v "^Mutations" | grep -v 'm3' > ${i}_mutations.txt

            rm ${i}_genomes.txt
        done

    done




# SimuRepID = Simulation Replicate ID
# ExpRepID = Experiment replicate ID
# Direction = direction of selection; set this at the line "for i in {T,F}" above
# BurninPath = path to the burnin files
# BurninFilename = name of the burnin file
# LCh = length of chromosome (CHANGE THIS ONLY WHEN THE BURNIN USES A DIFFERENT CHROMOSOME SIZE)
# RecRate = recombination rate (change the slim script if simulating multiple chromosomes)
# SampleSize = number of individuals to sample each generation
# NQTL = number of QTLs, (even number is recommended when the number is small)
# ESMean = absolute value of mean effect size
# ESDist = effect size distribution("f" for fixed or "e" for exponential),
# LowFreq = starting frequency preference (T if selecting for lower frequency, F if selecting for higher frequency or random)
# FreqBound = frequency bound (0.0~0.5, 0.0 if starting frequency is random)
# LowerPosBound = lower position bound (0 if random)
# UpperPosBound = upper position bound (LCh-1 if random)
# D = dominance coefficient (0.0~1.0, 1.0 for mutant being completely dominant and 0.0 for wildtype to being completely dominant.)
# Epistasis = F if there is no epistasis, T otherwise
# EpiSce = epistasis scenario; vector of size 9; first element must be 0; needs to be defined if Epistasis ==T; the nine value corresponds to phenotypes of genotypes in the following order: c(aabb, Aabb, AAbb, aaBb, AaBb, AABb, aaBB, AaBB, AABB); e.g. c(0,1,2,1,2,3,2,3,4) when there is no epistasis.
# PopSize = population size (can only downsample)
# SelectedSize = number of selected individuals in each generation
# OutPath = output path
