Please read the header to the file (line 1-94)

REMEMBER Rosetta is particially stochastic.  Make many structures and take the best!!!

Set the options in the file to those below.

--------------------------------------------------------------------------------------
## set CPUs to use
CPU={integer}

--------------------------------------------------------------------------------------
## To only score a file set . . . (only need 1 CPU)
RELAX=False
RIPDB={pdb.file}
MUTATE=False
MOSEQ=False
MOMSQ=False
MOSQE=False	
INTER=False

--------------------------------------------------------------------------------------
## To relax a file set . . . 
RELAX=True
RIPDB={pdb.file}
RINST={integer}   # number of total relaxations
MUTATE=False
MOSEQ=False
MOMSQ=False
MOSQE=False	
INTER=False

--------------------------------------------------------------------------------------
## To mutate a file using rosetta-fast-relax (each position selected will be mutated to best residue)
##   Additionally take the interface energy and compare and display the sequences and scores.
RELAX=True
RIPDB={pdb.file}
RINST={integer}   # number of total relaxations
RONSS={integer}   # number of top relaxations to take through to mutagenesis (should be less than 10% of RINST)
MUTATE=True
MIPDB=False
MIRES=False
MITHD=False
MOSEQ=True
MOMSQ=True
MOSQE=True	
MONST={integer}   # number of total mutation models.
INTER=True
IIFAC="A B"	      # define the interface (by group) to compare A and B to C, just define one group "A B" 
IOSQE=True

Then select the chain and positions to allow mutations (line 74)
The following mutates positions 7-11 on chain C
echo "
C 7 8 9 10 11 
" > mires.tx

--------------------------------------------------------------------------------------
## To thread a structure use the above options, but change 
MITHD={fasta.fa or file.pdb}
MITHC=C           # this is the chain to thread onto (line 79)

--------------------------------------------------------------------------------------
For additional options consult the file header.
