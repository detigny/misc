# misc
 Some miscellaneous python code I am using/writing for RAMSES simulations analysis
 
 Right now here we have:
 
 1) SMBH.py is a code for defining smbh pairs from ramses outputs as objects. From here some basic information can be readily extracted

 2) scripy.py has routines for quick sink separation outputting without the use of libraries, and format corrected sink outputs for ramses run restarting (when downgrading from the expected .csv sink output format that is present in RAMSES right now, which is bugged). This script is minimalistic by design as it was done to be as portable as possible.
 
 3) rafikov.py is a code I'm writing to test some predictions from Rafikov's 2012 model of accretion disks evolution based on angular momentum transport. Specifically the evolution and behaviour of viscosity as a function of the fraction of pressure that depends on radiation/thermodynamics. It includes functions that extract values such as radial profiles for surface density, T, Omega and dOmega/dR through .yt, does analysis of Fourier modes and some other things.
