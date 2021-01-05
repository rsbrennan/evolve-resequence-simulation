## This script is used to generate neutral populations used for the selection experiment.
## IMPORTANT: Change the output directory in the first step and the variables in the second step before running.
## Run the script on server using: nohup bash Burnin.sh > Burnin.nohup &
## Note: the number of generations to be included in the burnin process depends on the population size (typically ten times the population size) and needs to be changed in the Burnin.slim file. Please refer to the notes in the .slim file.

## Step 1: Ceate directories to store the outputs. Make sure to change the directory name in the first line.

cd /users/r/b/rbrennan/evolve-resequence-simulation/simulations  # Change this to a directory where you want to store all you simulation outputs.
mkdir Burnin
cd Burnin
for k in {1..100} # Number of simulation replicates that you want to create.
do
    mkdir 'SimRep'$k
done
cd ..

## Step 2: Run burnin using SLiM 2. Variables inside the loop are all customizable and can be changed as desired.

for k in {1..100} # Set the number of simulation replicates that you want to create.
do
	SECONDS=0
    # Set the path to the SLiM program in the next line
    ~/bin/SLiM/slim \
    -d SimRepID=$k  \
    -d Mu=2e-8 \
    -d RecRate=1e-8 \
    -d LCh=30000000 \
    -d BurninSize=1000 \
    -d "BurninPath='/users/r/b/rbrennan/evolve-resequence-simulation/simulations/Burnin/'" \
    -d "BurninFilename='Burnin.txt'" \
    /users/r/b/rbrennan/evolve-resequence-simulation/SlimScripts/Burnin.slim # Directory to the Burnin.slim file included in the simulation tool.

    duration=$SECONDS
    echo "done with replicat ${k}"
    echo "Time to run replicate ${k}: $(($duration / 60)) minutes and $(($duration % 60)) seconds"
    echo " "

done

# SimuRepID = Simulation Replicate ID
# Mu = mutation rate
# RecRate = recombination rate (change the slim script if simulating multiple chromosomes)
# LCh = length of chromosome
# BurninSize = size of the burnin populations (IF THIS NEEDS TO BE CHANGED, MAKE SURE TO CHANGE MUTATION RATE, RECOMBINATION RATE, AND NUMBER OF GENERATIONS ACCORDINGLY)
# BurninPath = path to the burnin files
# BurninFilename = name of the burnin file

