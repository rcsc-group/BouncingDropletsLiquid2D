#!/bin/bash

# Author: Radu Cimpeanu
# Date: 31/07/2023

# Create test folder (for debugging only)
rm -r Impact_Test
cp -r MasterImpact/ Impact_Test
cd Impact_Test/		

# Compile code to create the executable (including visualisation)
qcc -O2 -w -fopenmp -Wall DropImpact.c -lm -o DropImpact -L$BASILISK/gl -lglutils -lfb_tiny

# Specify parallelisation features, here set to 1 to begin with
export OMP_NUM_THREADS=1

# Physical parameters:
# 1. rho liquid
# 2. rho gas
# 3. mu liquid
# 4. mu gas
# 5. sigma (surface tension coefficient)
# 6. g acceleration
# 7. drop radius
# 8. initial drop velocity
# 9. simulation end time
# 10. max level

# Run executable with above parameters
./DropImpact 1000.0 1.225 8.9e-4 1.825e-5 0.0749 9.81 0.3e-3 0.35 20.0 10

cd ..
