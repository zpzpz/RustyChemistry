# RustyChemistry
An open-source computational chemistry software package. This was developed as a side-project, and I have no current intentions to flesh it out. If you are interested in helping me with that, please let me know!!!

# Currently Implemented
Hartree-Fock/STO-3G energy calculations, for elements up to Helium
Basis set read from file (as long as its STO-3G xD)
Read from input files

# Planned
Finish periodic table for STO-3G
- Properly implement spherical harmonics
Implementation of other basis sets
- Loading from file 
- Proper integrations
Geometry optimization
Vibrational analysis (IR)

# How to use
`$ ./RustyChemistry input_file.txt [-v] [--verbose]`
There are some example input files included in the repo. 

# How to build
```
$ git clone https://github.com/zpzpz/RustyChemistry
$ cd RustyChemistry
$ cargo build
``` 

there may be some issues with OpenBLAS-- I don't remember how I solved this. Basically, the current rust implementation can cause some issues.

# Credits
I followed [Adam Baskerville's guides](https://adambaskerville.github.io/posts/HFIntegrals/) for most of this. The code is heavily inspired, but my own.

# LICENSE
GPL 3.0
