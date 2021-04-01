# Setting up a MD simulation with NAMD and VMD 

- Download the PDB structure from Protein Data Bank (https://www.rcsb.org/). 
  In the field "Enter search term" write  4AKE. Then, use the option "Download Files"
  and choose "PDB Format". The PDB file will be downloaded to your local machine.
  Then, you can transfer this file to Kebnekaise.

- Load the VMD modules on Kebnekaise:

ml icc/2018.3.222-GCC-7.3.0-2.30 ifort/2018.3.222-GCC-7.3.0-2.30 impi/2018.3.222; ml VMD/1.9.3-Python-2.7.15

- Start VMD on the terminal by "vmd" (assuming your path is the one where you put the 4ake.pdb file).
Open VMD:  File -> New Molecule (choose 4ake.pdb) and close the Molecule File Browser box dialog.

- On VMD Main, go to Extensions -> Tk Console write these commands:
  - (cd $PATH_TO_PDB_FILE if you are working in your local machine)
  - set chaina [atomselect top "protein and chain A"]
  - $chaina writepdb 4ake_chaina.pdb
  - Quit VMD (File -> Quit)
  - with your preferred text editor delete the atom 1656 OXT

- Download the CHARMM36 parameters file toppar_c36_jul20.tgz from MacKerell's lab (http://mackerell.umaryland.edu/charmm_ff.shtml)
  - extract the files with the command *tar zxvf toppar_c36_jul20.tgz* 
  - copy and paste the file *top_all36_prot.rtf*,  *toppar_water_ions.str* and *par_all36_prot.prm* to the same folder level (same $PATH) than the 4ake_chaina.pdb structure

- write these lines into a file called *4ake.pgn*:

  - > package require psfgen
  - > topology top_all36_prot.rtf
  - > pdbalias residue HIS HSE
  - > pdbalias atom ILE CD1 CD
  - > segment U {pdb 4ake_chaina.pdb}
  - > coordpdb 4ake_chaina.pdb U
  - > guesscoord
  - > writepdb 4ake_corr.pdb
  - > writepsf 4ake_corr.psf

- Correcting Structure. On a Kebnekaise Linux terminal  write:
  - > vmd -dispdev text -e 4ake.pgn  
  - type *exit* on VMD terminal to close it. You will obtain the .psf and .pdb with hidrogen atoms
  (4ake_corr.psf, 4ake_corr.psf)

- Solvation. Start VMD. On VMD Tk Console type:
  - > package require solvate
  - > solvate 4ake_corr.psf 4ake_corr.pdb -t 10 -o 4ake_wb 
  - Close VMD and check if the solvated protein molecule files (4ake_wb.psf, 4ake_corr.psf) were written.
  *-t 10* creates a water box whose sides are 10 Angstrom from the more distant protein atom.

- Ionization. On VMD, open the Tk console, use the following command to place Na and Cl ions at a 150mM 
  concentration:
  - autoionize -psf 4ake_wb.psf -pdb 4ake_wb.pdb -sc 0.15 -cation SOD
  - quit VMD
  - change the names of the resulting files *ionized.pdb* to *4ake_ion.pdb* and the same for *.psf file

- On VMD open *4ake_ion.psf* and *4ake_ion.pdb*. Then on Tk console type:
  - > set everyone [atomselect top all]
  - > measure minmax $everyone
  - > measure center $everyone
  - Keep a record of the center of mass position. Quit VMD. 


Cell vectors can be taken from the 4ake_ion.pdb file first line:
CRYST1   59.001   76.997   73.292  90.00  90.00  90.00 P 1           




 
