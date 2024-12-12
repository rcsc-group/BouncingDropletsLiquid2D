#!/bin/bash

# Author: Radu Cimpeanu
# Date: 31/07/2023

# Additional resolution levels, radii or velocities can be added below
for LEVEL in 9; do
	for RADIUS in 0.0002; do
		for VELOCITY in 0.2; do

			# Copy all files to renamed folder based on key parameters
			cp -r MasterImpact/ Bounce-R$RADIUS-V$VELOCITY-Level$LEVEL
			cd Bounce-R$RADIUS-V$VELOCITY-Level$LEVEL/

			# Compile code to create the executable (including visualisation)
			qcc -O2 -w -fopenmp -Wall DropImpact.c -lm -o DropImpact -L$BASILISK/gl -lglutils -lfb_tiny

			# Specify parallelisation features
			export OMP_NUM_THREADS=4

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
                        ./DropImpact 949.0 1.225 1.9e-2 1.825e-5 0.0206 0.001 $RADIUS $VELOCITY 10.001 $LEVEL

			cd ..
		done
	done
done
