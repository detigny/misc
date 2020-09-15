# misc
 Some miscellaneous python code I am using/writing for RAMSES simulations analysis
 
 Right now here we have:
 
 1) Routine for defining smbh pairs from ramses outputs as objects. From here some basic information can be readily extracted

 2) Routine for quick sink separation outputting without the use of libraries, and format corrected sink outputs for ramses run restarting. This script is    minimalistic as it was done to be as portable as possible.
 
 3) Rafikov.py is a code I'm writing to test some predictions from Rafikov's 2018 model of accretion disks evolution based on angular momentum transport. It includes functions that extract values such as surface density, Omega and dOmega/dR through .yt 
