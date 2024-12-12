# Droplet rebounds off fluid baths

Mathematical model and direct numerical simulation code infrastructure for drop impact onto deep liquid pools, supporting collaborative work with [Dr. Kat Phillips](https://warwick.ac.uk/fac/sci/camacs/people/phillips/) and [Prof. Paul Milewski](https://science.psu.edu/math/people/ppm5454).

The present code base complements the preprint available at [arXiv:2406.16750](https://arxiv.org/abs/2406.16750). It consists of two parts:
1. The Matlab-based reduced-dimensional model [code](https://github.com/rcsc-group/BouncingDropletsLiquid2D/tree/main/Model),
2. The C-based direct numerical simulation [code](https://github.com/rcsc-group/BouncingDropletsLiquid2D/tree/main/DNS),

which are compared as part of this body of work.

## DNS Installation
* The code relies on [Basilisk](<http://basilisk.fr/>) to model the Navier-Stokes equations. See the [installation page](<http://basilisk.fr/src/INSTALL>) for instructions. 
* Full visualisation capabilities have been used in order to generate animations. These may be switched off depending on the local architecture.
* The two-phase non-coalescing fluid volume implementation by V. Sanjay available [here](https://github.com/VatsalSy/Lifting-a-sessile-drop/blob/master/CaseI/two-phaseDOD.h) has been successfully employed in this study to limit numerical artifacts during contact time.

## Running the DNS code
Once the Basilisk structure is in place, the driver code here is built either as stand-alone test to be run via *run_test.sh* in which a single run can be customised, or as a loop over variables of interest (in this case drop radius $R$ and initial drop velocity $V_0$, with one of each values added to the *run_master_example.sh* for brevity. Other parameters can be varied through this shell script, with both physical and computational handles provided. 

in order to navigate parameter sweeps in velocity $V_0$ and resolution level, with one of each values added to the run_master_example.sh for brevity. Other parameters can be varied through this shell script, with both physical and computational handles provided. 
