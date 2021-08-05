#!/bin/bash	

#author: James Ross : University of Leeds : March 2021 : james.f.ross@live.com
#		 This script prepares, runs and analyses a 'Rosetta fast relax with design' protocol.
#		 It accepts (as standard) a single input pdb to initite the process with RELAX (you will need to set up a loop over a list if you want multiple input pdbs)
#		 It can except a list of pdbs for the MUTATE function, supposing this is a pre-relaxed list (essential)
#
#		 Requirements:
#		 This is a bash script requiring Linux, mac or the 'windows subsystem for linux', in all cases the following must be on the PATH
#		 pymol   : version 2.4 or above
#		 python3 : with the seaborn package
#		 Rosetta : or
#				 : access to a HPC cluster managed by SGE with rosetta
#				 : a proxy jump setup is available
#				 : ssh public id must be shared to the cluster to allow the script to login without a password from your current machine.
#
#		 USAGE: This script is takes a single pdb and optionally relaxes, mutates and calculates an intereaction energy, providing sequence clustering and residue contributions
#	  	 USAGE: Fill out the options below.
#		 EXAMPLE: just run this script!
#		 
#		 ADVICE
# 		 It is advised that you provide a 'clean' pdb as input. Solvent removed, no alternative conformations.  Ligands and metal ions may require further support and scripting.
#		 Further information for ligands and metals should be sought from the rosetta help pages
# 		 Have you added "/usr/local/rosetta/main/source/bin/" to your .bashrc i.e.
# 		 export PATH=$PATH:/usr/local/rosetta/main/source/bin/
#		 For Local CPU management, please make sure that RONST, RONSS and MONST are equal or greater than the number of CPU's, which they should be anyway!!
#
#		 Warnings : This script removes the folders 'rosetta-0-initial', 'rosetta-1-relax', 'rosetta-2-mutate', 'rosetta-3-inter', 'rosetta-4-analysis' from the current directory
#					This script removes all '*.tx' files in the current directory
#					Makes use of 'nohup' and deletes nohup.out
#
# 	 	 IMPOVEMENTS TO MAKE:
#					Please send buggy copies of the script or suggestions to james.f.ross@live.com
#				
#		 WARNING !!!!
# 		 If your structure has come straight from MD and only contains chain 'P' but has multiple segment identifiers, the chains will be renamed to the final segment character.

########################################################################################################################################################################
# METHOD OPTIONS and OUTPUTS    #
########################################################################################################################################################################

# INPUT OPTIONS
# Option                    # Expected              # Comments

COMPUTE=LOCAL #arc4.leeds.ac.uk    # LOCAL Arc4 Arc3       # where to conduct the Rosetta calculations (remote options require '.leeds.ac.uk'
                          # if not local, it is assumed to be an arc cluster at the univerity of leeds, or another SGE cluster
                            # if not local, it is assumed you have already passed you ssh $PROXY id to the remote machine
HOURS=04:00:00              # HH:MM:SS for arc run  # this must exceed the maximum time, please wildly overestimate.
USER=chmjro                 # string                # your username on the remote machine, if required
EMAIL=chmjro@leeds.ac.uk     # email address         # your email address
remotedir=/nobackup/chmjro   # path                  # home directory on remote machine
PROXYJUMP=chmjro@remote-access.leeds.ac.uk           # FALSE or proxy jump user and address [user]@remote-access.leeds.ac.uk

CPU=10                      # Integer               # Number of CPUs to use  (make sure RINST and MONST are equal or greater than CPUs)
                            # only for local, HPC uses a defualt of 100, though 50 are usually given

########################
#       RELAX          #
########################
# It is highly recommended that structures are relaxed prior to any analysis with rosetta
# If you are starting from structures that have not been relaxed in rosetta, you MUST initiate this step
# you should generate a number of stuctures and pick the best to carry forward with the following options
RELAX=True                  # True False            # Relax structure : relax into rosetta forcefield, creates output folder 'rosetta-1-relax'
RIPDB=1coi.pdb              # Filename False        # Input pdb file. (if False you can provide a list of pre-relaxed structures for mutagenesis with MIPDB)
RIRMM=False              # move-map.file False    # provide a movemap for the relaxation if restraints are required.
RINST=10                    # Integer               # How many relaxations? : The number relaxations to make
RONSS=1                     # Integer               # How many results to carry through? : The number of best relaxations to carry through to the next stage 
                                                    # Typically do not exceed 10% of the total relaxations if making mutations
ROMET=True                # True False            # Produce output metrics and graphs for the relaxations. 		!!!!! INCOMPLETE !!!!!
ROPRS=False                 # True False            # Energy breakdown per residue output. - THIS PRODUCES A HUGE FILE!


########################
#       MUTATE         #
########################
# Mutagenesis of the protein structure using rosetta fast relax and design.
MUTATE=True                 # True False            # Mutate structure : Use rosetta fast relax and design, creates output folder 'rosetta-2-mutate'
MIPDB=False               # Filename False        # Input filelist : file with list of input pdbs, can contain a single pdb, 
                            #                         only use if no relaxation, otherwise "RONSS" determines input list.
                            #		                  i.e. 'False' is the default which takes inputs from the relaxation.
MIRES=False                 # Filename False        # Res file input : Provide a rosetta res file with mutation options
                            # 		                  if False provide mutation details below (mires.tx), or provide threading input with MITHD
MITHD=False                 # Filename False        # Threading a sequence to a structure : Provide a .pdb or .fa (fasta)
                            #                         Threading sequences MUST have the same number of residues as target structure, must have same numbering.
                            #                         must be the same chain, only one chain at present
                            #                         if Filename, provide threading details below (MITHC)
                            #                         if Filename, ensure 'mires.tx' (below) denotes the residue changes in order to generate output. 
MONST=10                # Integer               # many mutation runs? : The number of fast-relax-design runs to make PER relaxed structure (CPU must be a factor of this num
MOMET=True                 # True False            # Produce output metrics and graphs for the mutagenesis.  	!!!!! INCOMPLETE !!!!!
MOSEQ=True                  # True False            # Produce full output sequences (rosetta-full-sequence.fa)
MOMSQ=True                  # True False            # Produce position specific sequence outputs dependent on resfile, required for clustering (rosetta-spec-sequence.fa), single chain only
MOSQE=True                  # True False            # For each of the above sequence outputs, name - sequence - total energy
MOPRS=False                 # True False            # Energy breakdown per residue output. - THIS PRODUCES A HUGE FILE!

# STANDARD MUTATION (pick positions to mutate)
# if "MIRES=False" define mutational parameters here. Alternatively if set MITHD=filename.file providing a threading sequence (.fa) or structure (.pdb), use mires.tx to determine which residues to output.
# Provide chain and residue positions for mutagenesis, by default we will mutate to all amino acids except cysteine. Otherwise provide resfile.
# One chain per line, with residue positions, as below
rm mires.tx 2>/dev/null 
echo "
A 11 12 15 16 19 20
" > mires.tx
# PICK MUTATION TYPE (see here https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/resfiles)
PICKACID=ALLAAxc            # ALLAAxc (all amino acid, not cys) or POLAR (DEHKNQRST) or APOLAR (ACFGILMPVWY)
                            # or PIKAA-ACDEFGH (either native AA or ACDEFGH, choose yourself) or NOTAA-CW (any but not C or W)

# THREADING SEQUENCE TO STRUCTURE 
# set RONSS = 1
# set MONST = the number of structures for each threading, ie. if one mutation ~10, if 5 mutations about ~100, if 10+ mutations around ~1000, pick lowest
# if "MITHD=Filename" define threading parameters here.  
# Reqires the above mires.tx file to have the mutated positions, IF the 'combined-output.csv' is required
# residue numbers will be ignored, this replaces aminoacids sequenctially from source to target.
MITHC=B                     # Character False       # chain idenifier for replaced sequence. Not currently applicable for mutil chain replacements	

# DETERMINE STRINGENCY OF MOVEMAP ( how large is you 'flexible shell' around mutations)
flexsc=9                    # integer               # shell, in angstrom, around mutant residues which will have flexible sidechains (and backbone)
flexbb=9                    # integer               # shell, in angstrom, around flexible sidechains which will have flexible backbone
subjump=YES                 # YES or NO             # can subunits rigid-body move independent of each other !!!!! NO is NOT currently available !!!!

########################
#      INTERFACE       #
########################
# Analysis of the protein interface between chains.
INTER=True                  # True False            # Analyse interface : Use the InterfaceAnalyzer application to calculate ddG in Rosetta energy units (REU), creates output folder 'rosetta-3-inter'
IIFAC="A"                   # chain_names           # qoute chains of single group, if analaysing the interface between chain groups A and B vs C and D then use "A B"
IOSQE=True                  # True False            # For each of the analysed structures output, name - sequence - total energy - interface energy
IOMET=True                 # True False            # Produce output metrics and graphs for the Interface Energy.			!!!!! INCOMPLETE !!!!!
IOPRS=True                  # True False            # Energy breakdown per residue across interface output. - THIS PRODUCES A HUGE FILE!

########################
#       CLUSTER        #
########################
# Sequence clustering and energy analysis.
action=sequence              # either 'rmsd' or 'sequence' or 'both' or 'False'    !!!!!  RMSD MATRIX IS CURRENTLY DISABLED   !!!!!!!
                            # if sequence, will use the mutant sequence.
                            # if rmsd will use the whole protein or a selected range below			
# select a residue range for the rmsd
rrange=all                  # 1-25 or 'all'         # currently rrange is not working
chain=$IIFAC                # A or 'all'            # currently 'all' is not working
ndigits=5

# THIS IS THE END OF THE USER DEFINED INPUT, NOW PLEASE RUN THE SCRIPT.

########################################################################################################################################################################
# INTERNAL CHECKS AND FUNCTIONS    #
########################################################################################################################################################################

# allow temporary file preservation
DEBUG=True			# True False    # If True, then intermediate files are not cleaned up!

# find the current working directory
oridir=$(pwd)

# force passed unknown reidues
passunk=no			# yes no	# if yes, unknown residues will be deleted or cause error

########################
# PROGRAM DEPENDENCIES #
########################

########################
# rosetta			# 
if [ $COMPUTE = LOCAL ] ; then 
	if [ $(which score_jd2.*.linuxgccrelease | wc -l ) = 0 ] ; then 
		echo 'Rosetta not located on your path, seaching now . . .'
		rpath=$(locate get_fasta_from_pdb.py | head -1)
		if [ $(echo $rpath | wc -l ) = 1 ] ; then 
			rbpath=$(echo $rpath | sed 's#/tools/protein_tools/scripts/get_fasta_from_pdb.py#/main/source/bin/#g')
			rdpath=$(echo $rpath | sed 's#/tools/protein_tools/scripts/get_fasta_from_pdb.py#/main/database/#g')
			PATH=$PATH:$rbpath
		else	
			echo 'could not find rosetta'
			exit 0 
		fi
	fi

	# static or default versions?  Dependent on the Rosetta installation, executables can have either 'static' or 'default' in the executable name.
	if [ $( ls $rbpath/score_jd2.static* 2>/dev/null | wc -l ) = 1 ] ; then 
		version=static
		echo "using rosetta 'static' build"
	elif [ $( ls $rbpath/score_jd2.default* 2>/dev/null | wc -l ) = 1 ] ; then 
		version=default 
		echo "using rosetta 'default' build"
	else
		echo 'unable to find correct build of rosetta'
		exit 0
	fi
fi

########################
# pymol 			# check on PyMOl
rm nohup.out 2>/dev/null 
nohup pymol -qc 2>/dev/null 
test=$(grep "command not found" nohup.out | wc -l )
if [ $test = "0" ] ; then
	echo "PyMOL detected"
else	
	echo -en "\rPyMOL ***NOT*** detected"
	if [ $MIRES = False ] ; then
		echo "no resfile detected and no PyMOL program detected"
		echo "I cannot make the resfile without PYMOL_MAIN"
		exit 0
	fi
fi

########################
# python			# Check on python3 with seaborn


########################
#    OPTIONS CHECK     #
########################

# remove underscores from file names and fix stuff
if [ $RELAX = True ] ; then
	RIPDB=$(echo $RIPDB | sed 's/_/-/g')  2>/dev/null 
elif [ $MUTATE = True ] ; then
	if [ $(head -1 $MIPDB | rev | cut -c 1-3 | rev ) = pdb ] ; then 
		sed -i 's/_/-/g' $MIPDB
	else 
		MIPDB=$(echo $RIPDB | sed 's/_/-/g')  2>/dev/null 
	fi
fi 

# checking for contridictory or badly formatted inputs

if [ $RELAX = True ] ; then
	if [ -f $RIPDB ] ; then 
		echo "read initial pdb file as $RIPDB for relaxations"
		# how many chains?
		cnum=$(grep ' CA ' $RIPDB | cut -c 22 | sort | uniq  | sed '/^[[:space:]]*$/d' | wc -l )
		if [ $cnum = 0 ] ; then
			echo "Hey!! there aren't any chain names, please name your chains!"
			exit 0
		fi
		# how many segments?
		snum=$(grep ' CA ' $RIPDB | cut -c 73-76 | sort | uniq | sed '/^[[:space:]]*$/d' | wc -l )
		echo $cnum chains and $snum segments
		# remove any HET atoms!
		if [ $( grep HET $RIPDB | wc -l ) != 0 ] ; then 
			echo 'found HET atoms or segments, these have been removed'
		fi 
		grep -v HET $RIPDB > fix-$RIPDB
		if [ $snum -gt 0 ] ; then 
			if [ $cnum -lt $snum ] ; then 
				echo 'renaming chains due to an inconsistancy between number of chains ('$cnum') and number of segments ('$snum')
		to prevent this change occuring, please fix manually.'
				if [ $snum = $(grep ' CA ' $RIPDB | cut -c 76 | sort | uniq | wc -l ) ] ; then 
					rm  fixt-$RIPDB 2>/dev/null
					for segi in $(grep ' CA ' $RIPDB | cut -c 73-76 | sort | uniq) ; do 
						chain=$(echo $segi | cut -c 4)
						grep $segi fix-$RIPDB | sed "s/./$chain/22" >> fixt-$RIPDB
						echo TER >> fixt-$RIPDB 
					done
					mv fixt-$RIPDB fix-$RIPDB 
				else
					echo "please check input file chain and segment identifiers, they are currently incompatible with this script."
					exit 0
				fi
			else	
				echo "please check input file chain and segment identifiers, they are currently incompatible with this script."
				exit 0
			fi
		fi
		
		echo '
remove hydrogens 
save temp.pdb
quit' > rmhydro.pml
		sed -i "s/HSD/HIS/g" fix-$RIPDB
		pymol -qc fix-$RIPDB rmhydro.pml 1>/dev/null 
		
		RIPDBn=$(echo $RIPDB | sed 's/.pdb//g')
		mv temp.pdb $RIPDBn-relax.pdb
		RIPDB=$RIPDBn-relax.pdb
	else	
		echo "please check input file for relaxation, none found."
		exit 0
	fi
elif [ $MUTATE = True ] ; then 
	if [ -f $MIPDB ] ; then 
		echo "read initial as pdb list file $MIPDB with $(wc -l < $MIPDB) structures for mutations"
	else	
		echo "please check input file for mutations, none found."
		exit 0
	fi
elif [ $INTER = True ] ; then 
	if [ -f $IIPDB ] ; then 
		echo "read initial as pdb list file $IIPDB with $(wc -l < $IIPDB) structures for Interface analysis"
	else	
		echo "please check input file for interface analysis, none found."
		exit 0
	fi
fi
# check threading input.
if [ $MITHD != False ] ; then
	if [ $( echo $MITHD | awk -F'.' '{print $2}' ) = fa ] ; then 
		echo "detected threading with fasta file"
		echo "will thread sequences of $MITHD onto chain $MITHC of $RIPDB "
	elif [ $( echo $MITHD | awk -F'.' '{print $2}' ) = pdb ] ; then 
		echo "detected threading with pdb file,"
		echo "will thread chain $MITHC of $MITHD onto chain $MITHC of $RIPDB"
	else
		echo "detected threading but file type not recognised, please use .pdb or .fa"
		exit 0
	fi
	if [ ! -f "$MITHD" ]; then
    echo "$MITHD does not exist, exiting . . ."
	exit 0
	fi
fi

# check res file clasehes.
if [ $MIRES != False ] && [ $MITHD != False ] ; then 
	echo 'We cannnot have both a res file ('$MIRES') and a theading template file ('$MITHD')'
	exit = 0 
fi

# resifixer
resifix ()
{
# protein residues
sed -i 's/HIE/HIS/g;s/HSD/HIS/g;s/HSE/HIS/g;s/CYX/CYS/g' $1
# glycan residues
sed -i 's/0fA/FUC/g;s/0LB/GAL/g;s/0SA/SIA/g;s/2LB/GAL/g;s/3VB/NGA/g;
		s/4GB/BGC/g;s/WLB/GAL/g;s/WYB/NDG/g' $1
# glycan atom nomenclaure
sed -i 's/ C2N/ CN2/g;s/ CME SIA/CAN5 SIA/g;s/ CME/CAN2/g;s/ O2N/OCN2/g;s/ O5N/OCN5/g' $1
}

# checking for non protein residues
if [ $RELAX = True ] ; then
	tester=$RIPDB
elif [ $MUTATE = True ] ; then
	if [ $(head -1 $MIPDB | rev | cut -c 1-3 | rev ) = pdb ] ; then 
		tester=$(head -1 $MIPDB)
	else 
		tester=$MIPDB
	fi
fi 
cp $tester test.pdb

grep "ATOM\|HETATM" test.pdb | cut -c 18-20 | sort | uniq | sed '/ALA/d;/CYS/d;/ASP/d;/GLU/d;/PHE/d;
																/GLY/d;/HIS/d;/ILE/d;/LYS/d;/LEU/d;
																/MET/d;/ASN/d;/PRO/d;/GLN/d;/ARG/d;
																/SER/d;/THR/d;/VAL/d;/TRP/d;/TYR/d;
																/FUC/d;/GAL/d;/SIA/d;/NGA/d;/BGC/d;
																/NDG/d' > non-familiar.tx
if [ $( wc -l < non-familiar.tx) != 0 ] ; then 
echo "Non familiar residues detected! $(cat non-familiar.tx | tr '\n' ' ')"
echo "attempting fix"
resifix test.pdb

grep "ATOM\|HETATM" test.pdb | cut -c 18-20 | sort | uniq | sed '/ALA/d;/CYS/d;/ASP/d;/GLU/d;/PHE/d;
														/GLY/d;/HIS/d;/ILE/d;/LYS/d;/LEU/d;
														/MET/d;/ASN/d;/PRO/d;/GLN/d;/ARG/d;
														/SER/d;/THR/d;/VAL/d;/TRP/d;/TYR/d;
														/FUC/d;/GAL/d;/SIA/d;/NGA/d;/BGC/d;
														/NDG/d' > non-familiar.tx
fi
# are these glycans?
if [ $(grep "ATOM\|HETATM" test.pdb | cut -c 18-20 | sort | uniq | grep "FUC\|GAL\|SIA\|NGA\|BGC\|NDG\|Gal\|Glc\|Neu" | wc -l ) != 0 ] ; then 
	glycan=yes
	echo "Glycans detected!"
	echo "-include_sugars
-alternate_3_letter_codes pdb_sugar
-load_PDB_components false
-auto_detect_glycan_connections " > glycan.tx

fi 

if [ $( wc -l < non-familiar.tx) != 0 ] ; then 
echo Unknown residues detected 
if [ $passunk = yes ] ; then 
	echo "continuing, unknown residues will be deleted."
else	
	echo "unknown residues detected and not fixed! "
	echo "fix the pdb or activate the passunk option ingnore unknowns"
	echo "exiting . . . "
	exit 0
fi 
fi 


########################
#      transpose       #
########################

# pass first argument as $1 which should be any matrix
transpose ()
{
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $1 
}


########################
#       pdb2seq        #
########################

# pass first argument ($1) which should be a pdb file
# pass second argument ($2) which should be a chain character
seqgen ()
{
if [ -f "$2" ] ; then 
	grep 'CA' $1 | grep '^ATOM' | grep "[A-Z][A-Z][A-Z] $2 " | awk '{print $4}' > temp.tx
	sed 's/ALA/A/g; s/CYS/C/g; s/ASP/D/g; s/GLU/E/g; s/PHE/F/g; 
	     s/GLY/G/g; s/HIS/H/g; s/ILE/I/g; s/LYS/K/g; s/LEU/L/g; 
		 s/MET/M/g; s/ASN/N/g; s/PRO/P/g; s/GLN/Q/g; s/ARG/R/g; 
		 s/SER/S/g; s/THR/T/g; s/VAL/V/g; s/TRP/W/g; s/TYR/Y/g' < temp.tx  > transin.tx 
	transpose transin.tx | sed 's/ //g' > temp3.tx
	#sed 's/ //g'  transout.tx > temp3.tx
	echo ">"$1'_'$2 > head.tx
	cat head.tx temp3.tx > temp4.tx
	cat temp4.tx >> seqgen.fa
else
	for i in $(grep 'CA' $1 | grep '^ATOM' | cut -c 22 | sort | uniq) ; do 
		grep 'CA' $1 | grep '^ATOM' | grep "[A-Z][A-Z][A-Z] $i " | awk '{print $4}' > temp.tx
		sed 's/ALA/A/g; s/CYS/C/g; s/ASP/D/g; s/GLU/E/g; s/PHE/F/g; 
		     s/GLY/G/g; s/HIS/H/g; s/ILE/I/g; s/LYS/K/g; s/LEU/L/g; 
			 s/MET/M/g; s/ASN/N/g; s/PRO/P/g; s/GLN/Q/g; s/ARG/R/g; 
			 s/SER/S/g; s/THR/T/g; s/VAL/V/g; s/TRP/W/g; s/TYR/Y/g' < temp.tx  > transin.tx 
		transpose transin.tx | sed 's/ //g' > temp3.tx
		#sed 's/ //g'  transout.tx > temp3.tx
		echo ">"$1'_'$i > head.tx
		cat head.tx temp3.tx > temp4.tx
		cat temp4.tx >> seqgen.fa
	done
fi
}

########################
#   check CA RMSD      #   
########################
rmsmatrix ()
{

#isolate 
count=1
rm test.pdb 2>/dev/null 
for i in $(cat pdb.list) ; do 
echo MODEL    $count >> test.pdb
cat $i >> test.pdb
echo ENDMDL >> test.pdb
count=$(($count +1))
echo -en "\r$count"
done

echo "load test.pdb
intra_fit all " > rmscur.pml
rm rmsout*.tx 2>/dev/null
for i in $( seq 1 $(wc -l < pdb.list)) ; do
echo "intra_rms_cur polymer, $i" >> rmscur.pml
done
echo "save check.pdb, c. A & i. 1" >> rmscur.pml
rm check.pdb 2>/dev/null
pymol -qc test.pdb rmscur.pml > temp.tx

while [ ! -f check.pdb ] ; do
echo -en "\rgenerating rmsd's"
sleep 1
done

rm check.pdb
sleep 2

for i in $(seq 1 $( wc -l < pdb.list) ) ; do 
if [ $i = 1 ] ; then 
echo 0.000 > zero.tx
grep "state $i$" temp.tx | awk '{print $2}' > zerod.tx
cat zero.tx zerod.tx > zerot.tx ; mv zerot.tx rmsout$i.tx
else
j=$(( $i - 1 ))
grep "state $i$" temp.tx | awk '{print $2}' | sed ''$j' a 0.000' > rmsout$i.tx
fi
done

rm rms.tx  trms.tx 2>/dev/null
touch rms.tx  
for i in $( seq 1 $(wc -l < pdb.list)) ; do
echo -en "\rmatricising $i"
sed 's/-1.0/0.0/g;s/\[/ /g;s/\]//g' rmsout$i.tx | sed 's/\.//g' | sed -n '1,'$(wc -l < pdb.list)'p' > trms.tx
paste -d',' rms.tx trms.tx > ttrms.tx
mv ttrms.tx rms.tx
done
sed -i 's/^,//g' rms.tx
rm rmsout*.tx

transpose pdb.list | sed 's/ /,/g' > tpdb.list
cat tpdb.list rms.tx > rms2.tx
echo name > rmsh.tx
cat pdb.list >> rmsh.tx
paste -d',' rmsh.tx  rms2.tx | sed 's/.pdb//g' > matrix.csv
}


########################
#    Graph function    #
########################

########################
#   REMOTE COMPUTE     #
########################

if [ $PROXYJUMP = FALSE ] ; then
	PROXY=''
else
	PROXY=-oProxyJump=$PROXYJUMP
fi

date=$(echo $(date -R | cut -c 6-7,9-11,15-25 | sed 's/ /-/g;s/://g')-$RANDOM)

echo "
# Arc4 - Rosetta
# James Ross - j.f.ross@leeds.ac.uk - 2020

# specify environment and working directory
#$ -V -cwd
# Time request
#$ -l h_rt=$HOURS
# RAM request
#$ -l h_vmem=4G
# Request task arrays
#$ -t 1-XXX     # total jobs
#$ -tc YYY       # max simultaineous
# Notificiations
## -m be
## -M $EMAIL
# Specify Job
module load test rosetta/3.10

cd VVV" > remote.tx
echo '
state=$(echo 0000$SGE_TASK_ID | rev | cut -c 1-5 | rev )

echo "' >> remote.tx


########################################################################################################################################################################
#      INITIAL SCORE      #
###########################
if [ 1 = 0 ] ; then # just not sure if this is useful
rm -r $oridir/rosetta-0-initial score.sc 2>/dev/null 
mkdir $oridir/rosetta-0-initial
cd $oridir/rosetta-0-initial

echo ''
echo "
-database /usr/local/rosetta/main/database
-nstruct 1
-out:no_nstruct_label
-score:weights ref2015" > sflag.file


if [ $RELAX = True ] ; then
	echo '0.0 Initial Score'
	if [ $RIPDB != False ] ; then
		resifix $RIPDB
		echo '0.1 - Scoring initial pdb'
		
		# For local compute
		if [ $COMPUTE = LOCAL ] ; then 
			cp $oridir/$RIPDB .
			nohup score_jd2.$version.linuxgccrelease -in:file:s $RIPDB @sflag.file  2>/dev/null
			wait
		fi 

		# for remote compute
		if [ $COMPUTE != LOCAL ] ; then 
			date=$(echo $(date -R | cut -c 6-7,9-11,15-25 | sed 's/ /-/g;s/://g')$RANDOM)
			cp ../remote.tx .
			ssh $PROXY $USER@$COMPUTE mkdir $remotedir/$date-RosettaAutoScore
			
			scp $PROXY $oridir/$RIPDB $USER@$COMPUTE:$remotedir/$date-RosettaAutoScore/. 1>/dev/null
			echo "-database /apps/applications/rosetta/3.10/1/default/main/database
-nstruct 1
-out:no_nstruct_label
$(cat $oridir/glycan.tx)
-score:weights ref2015 \" > flag_"'$SGE_TASK_ID'".file
/apps/applications/rosetta/3.10/1/default/main/source/bin/score_jd2.default.linuxgccrelease -in:file:s $RIPDB @flag_"'$SGE_TASK_ID'".file
" >> remote.tx
			cat remote.tx | sed "s/XXX/1/g;s/YYY/1/g;s#VVV#$remotedir/$date-RosettaAutoScore#g" > run-rosetta-jd2-score.sh
			scp $PROXY run-rosetta-jd2-score.sh $USER@$COMPUTE:$remotedir/$date-RosettaAutoScore/. 1>/dev/null
			ssh $PROXY $USER@$COMPUTE " cd $remotedir/$date-RosettaAutoScore ; qsub run-rosetta-jd2-score.sh" > $date.tx 
			arcjobid=$( awk '{print $3}' $date.tx | awk -F'.' '{print $1}' )
			while [ $(ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid 2>/dev/null | wc -l ) != 0 ] ; do
				ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid | cut -c 41-43,104- > arcjobid.tx
				if [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == qw ] ; then 
					echo -en "\rsubmitted initial score to $COMPUTE, queueing "$(awk '{print $2}' arcjobid.tx)
				elif [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == r ] ; then 
					echo -en "\rsubmitted initial score to $COMPUTE, running last job                            "
				elif [ $( wc -l < arcjobid.tx) -gt 1 ] ; then 
					echo -en "\rsubmitted initial score to $COMPUTE, queueing, "$(grep qw arcjobid.tx | awk '{print $1}')", running "$(grep r arcjobid.tx | wc -l )
				fi
				sleep 10
			done
			scp $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoScore/score.sc scoret.sc 1>/dev/null
			sleep 1
			if [ ! -f score.sc ] ; then 
				mv scoret.sc score.sc
			else
				sed '1,2d' scoret.sc >> score.sc
			fi
			scp $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoScore/*.pdb . 1>/dev/null
			ssh $PROXY $USER@$COMPUTE rm -r $remotedir/$date-RosettaAutoScore
			if [ $( wc -l < score.sc ) -gt 2 ] ; then 
				echo -en "\rsubmitted initial score to $COMPUTE, job successful              "
				if [ $DEBUG != True ] ; then 
					rm *.tx 
				fi
			else
				echo -en "\rsubmitted initial score to $COMPUTE, job failed              "
				exit 0
			fi
			echo ''
		fi			
	fi

elif [ $MUTATE = True ] ; then 	 # <-- This has been disabled as probably not wanted
	echo '0.0 - Initial Score'
	cp ../$MIPDB .
	nuministruc=$(wc -l < $MIPDB)
	echo '0.1 - Scoring initial pdb list'
	
	# For local compute
	if [ $COMPUTE = LOCAL ] ; then 
		if [ $nuministruc -gt $CPU ] ; then 
			for i in $(seq 1 $CPU) ; do
				echo '
				for j in $(seq '$i' '$CPU' $(wc -l < '$MIPDB') ) ; do
				pdb=$(sed -n '"''"'$j'"'"'p'"'"' '$MIPDB')
				resifix $pdb
				cp $oridir/$pdb .
				nohup score_jd2.'$version'.linuxgccrelease -in:file:s $pdb @sflag.file  2>/dev/null 
				done' > run-cpu$i
				chmod 755 run-cpu$i
				./run-cpu$i &
			done
			wait
		else
			for i in $(cat $MIPDB) ; do
				cp $oridir/$i .
				nohup score_jd2.$version.linuxgccrelease -in:file:s $i @sflag.file  2>/dev/null
				wait
			done
		fi
	fi
	
	# for remote compute
	if [ $COMPUTE != LOCAL ] ; then 
		
		for i in $(cat $MIPDB) ; do 
			date=$(echo $(date -R | cut -c 6-7,9-11,15-25 | sed 's/ /-/g;s/://g')$RANDOM)
			cp ../remote.tx .
			ssh $PROXY $USER@$COMPUTE mkdir $remotedir/$date-RosettaAutoScore
			resifix $i
			scp $PROXY $oridir/$i $USER@$COMPUTE:$remotedir/$date-RosettaAutoScore/. 1>/dev/null
			echo "-database /apps/applications/rosetta/3.10/1/default/main/database
-nstruct 1
-out:no_nstruct_label
$(cat $oridir/glycan.tx)
-score:weights ref2015 \" > flag_"'$SGE_TASK_ID'".file
/apps/applications/rosetta/3.10/1/default/main/source/bin/score_jd2.default.linuxgccrelease -in:file:s $i @flag_"'$SGE_TASK_ID'".file
	" >> remote.tx
			cat remote.tx | sed "s/XXX/1/g;s/YYY/1/g;s#VVV#$remotedir/$date-RosettaAutoScore#g" > run-rosetta-jd2-score.sh
			scp $PROXY run-rosetta-jd2-score.sh $USER@$COMPUTE:$remotedir/$date-RosettaAutoScore/. 1>/dev/null
			ssh $PROXY $USER@$COMPUTE " cd $remotedir/$date-RosettaAutoScore ; qsub run-rosetta-jd2-score.sh" > $date.tx 
			arcjobid=$( awk '{print $3}' $date.tx | awk -F'.' '{print $1}' )
			while [ $(ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid 2>/dev/null | wc -l ) != 0 ] ; do
				ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid | cut -c 41-43,104- > arcjobid.tx
				if [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == qw ] ; then 
					echo -en "\rsubmitted initial score to $COMPUTE, queueing "$(awk '{print $2}' arcjobid.tx)
				elif [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == r ] ; then 
					echo -en "\rsubmitted initial score to $COMPUTE, running last job                 "
				elif [ $( wc -l < arcjobid.tx) -gt 1 ] ; then 
					echo -en "\rsubmitted initial score to $COMPUTE, queueing, "$(grep qw arcjobid.tx | awk '{print $1}')", running "$(grep r arcjobid.tx | wc -l )
				fi
				sleep 10
			done
			
			scp $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoScore/score.sc scoret.sc 1>/dev/null
			sleep 1
			if [ ! -f score.sc ] ; then 
				mv scoret.sc score.sc
			else
				sed '1,2d' scoret.sc >> score.sc
			fi
			scp $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoScore/*.pdb . 1>/dev/null
			ssh $PROXY $USER@$COMPUTE rm -r $remotedir/$date-RosettaAutoScore
			if [ $( wc -l < score.sc ) -gt 2 ] ; then 
				echo -en "\rsubmitted initial score to $COMPUTE, job successful          "
				if [ $DEBUG != True ] ; then
					rm *.tx 
				fi
			else
				echo -en "\rsubmitted initial score to $COMPUTE, job failed              "
				exit 0
			fi
			echo ''
		done
	fi

fi 
if [ $DEBUG != True ] ; then
	rm sflag.file
fi
cd $oridir
fi
########################################################################################################################################################################
#      BODY RELAX      #
########################################################################################################################################################################

if [ $RELAX = True ] ; then

	rm -r $oridir/rosetta-1-relax 2>/dev/null 
	mkdir $oridir/rosetta-1-relax
	cd $oridir/rosetta-1-relax
	#cp $oridir/$RIPDB rosetta-1-relax/$RIPDB

	echo '1.0 Relaxation'
	
	# For local compute	
	if [ $COMPUTE = LOCAL ] ; then 
	################################################	
	########################  CPU CONTROL
			
		# duplicate input files for each CPU
		rm ripdb.tx 2>/dev/null 
		for i in $(seq 1 $CPU) ; do 
			echo $oridir/$RIPDB >> ripdb.tx
		done
								# debug!
		cp $oridir/$RIPDB .
		RIPDBL=ripdb.tx


		# perform residue fix on pdbs
		echo "1.0a	fixing residue nomenclature"
		for i in $(cat $RIPDBL | sort | uniq) ; do 
			resifix $i 
		done

		echo '1.1 - Running local relaxations'
	################################################	
	########################  RUN RELAXATION
		
		# process the input list for the number of CPUs
		for i in $(seq 1 $CPU) ; do

			#make a temp directory per cpu
			mkdir $oridir/rosetta-1-relax/temp$i
			#make a flag file per cpu
			echo "
-out:path:all $oridir/rosetta-1-relax/temp$i
-relax:fast
$(cat $oridir/glycan.tx 2>/dev/null)
-database /usr/local/rosetta/main/database
-nstruct $(( $RINST / $CPU ))
-score:weights ref2015" > rflag$i.file
			# is there a move map restaint file?
			if [ $RIRMM != False ] ; then 
				cp $oridir/$RIRMM .
				echo "-in:file:movemap $RIRMM" >> rflag$i.file
			fi
			#make a rosetta executable per cpu and run them
			echo '
for j in $(seq '$i' '$CPU' $(wc -l < '$RIPDBL') ) ; do
pdb=$(sed -n '"''"'$j'"'"'p'"'"' '$RIPDBL')
nohup relax.'$version'.linuxgccrelease -in:file:s $pdb @rflag'$i'.file  2>/dev/null 
done' > run-cpu$i
			chmod 755 run-cpu$i
			./run-cpu$i &
		done
		wait
			#for i in $(cat ripdb.tx) ; do
		#	rm $i
		#done
		rm $RIPDB
		cd $oridir/rosetta-1-relax
		################################################	
		########################  MERGE OUTPUTS	
		echo '1.2 - merging outputs'
		# the output pdbs and score files are in seperate directories, with the same names, lets rename and merge
		for i in $(seq 1 $CPU) ; do
			#count files in destination directory
			#pdbcount=$(ls $oridir/rosetta-1-relax/*.pdb | wc -l )
			if [ $i = 1 ] ; then 
				#move everything from the first temp directory to the destination directory
				mv $oridir/rosetta-1-relax/temp$i/* $oridir/rosetta-1-relax
			else
				#move everything from the other temp directorys to the destination directory
				for j in $(ls temp$i/*.pdb) ; do 
					pdbcount=$(ls *.pdb | wc -l )
					#echo pdb count = $pdbcount
					#define old number (with leading 0's)
					oldnum=$(echo $j | rev | cut -c 1-8 | rev | cut -c 1-4)
					#define old number (without leading 0's)
					num=$(echo $j | rev | cut -c 1-8 | rev | cut -c 1-4 | sed 's/0//g')
					#define new number (without leading 0's)
					newnum=$(( 1 + $pdbcount))
					#define new number (with leading 0's)
					newnumname=$(echo 000$newnum | rev | cut -c 1-4 | rev)
					#define file name without number
					nonumname=$(echo $j | rev | cut -c 9- | awk -F'/' '{print $1}' | rev)

					# rename and move pdb and contents of score file
					mv $j $nonumname$newnumname.pdb
					grep $nonumname$oldnum temp$i/score.sc | sed "s/$nonumname$oldnum/$nonumname$newnumname/g" >> score.sc
					#echo check4
				done
				#echo check3
			fi
			#echo check2
			#count files in destination directory
		done
		if [ $DEBUG != True ] ; then
			rm -r temp* run-cpu* rflag*
		fi
		#echo check1
		cd $oridir
	fi
	
	# for remote compute
	if [ $COMPUTE != LOCAL ] ; then 
		echo '1.1 - Running remote relaxations'
		date=$(echo $(date -R | cut -c 6-7,9-11,15-25 | sed 's/ /-/g;s/://g')$RANDOM)
		cp ../remote.tx .
		ssh $PROXY $USER@$COMPUTE mkdir $remotedir/$date-RosettaAutoRelax
		resifix $RIPDB
		scp $PROXY $oridir/$RIPDB $USER@$COMPUTE:$remotedir/$date-RosettaAutoRelax/. 1>/dev/null
		echo "-database /apps/applications/rosetta/3.10/1/default/main/database
-out:suffix _"'$SGE_TASK_ID'"
-relax:fast
-out:suffix _"'$state'"
-out:file:scorefile score.sc
$(cat $oridir/glycan.tx 2>/dev/null)
-nstruct 1
-no_nstruct_label
-score:weights ref2015 \" > flag_"'$SGE_TASK_ID'".file
/apps/applications/rosetta/3.10/1/default/main/source/bin/relax.linuxgccrelease -in:file:s $RIPDB @flag_"'$SGE_TASK_ID'".file
" >> remote.tx
		cat remote.tx | sed "s/XXX/$RINST/g;s/YYY/100/g;s#VVV#$remotedir/$date-RosettaAutoRelax#g" > run-rosetta-relax.sh
		scp $PROXY run-rosetta-relax.sh $USER@$COMPUTE:$remotedir/$date-RosettaAutoRelax/. 1>/dev/null
		ssh $PROXY $USER@$COMPUTE " cd $remotedir/$date-RosettaAutoRelax ; qsub run-rosetta-relax.sh" > $date.tx 
		arcjobid=$( awk '{print $3}' $date.tx | awk -F'.' '{print $1}' )
		while [ $(ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid 2>/dev/null | wc -l ) != 0 ] ; do
			ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid | cut -c 41-43,104- > arcjobid.tx
			if [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == qw ] ; then 
				echo -en "\rsubmitted relax to $COMPUTE, queueing "$(awk '{print $2}' arcjobid.tx | awk -F':' '{print $1}')
			elif [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == r ] ; then 
				echo -en "\rsubmitted relax to $COMPUTE, running last job                   "
			elif [ $( wc -l < arcjobid.tx) -gt 1 ] ; then 
				echo -en "\rsubmitted relax to $COMPUTE, queueing "$(grep qw arcjobid.tx | awk '{print $2}')", running "$(grep r arcjobid.tx | wc -l )"        "
			fi
			sleep 30
		done
		scp $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoRelax/score.sc . 1>/dev/null
		scp $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoRelax/*.pdb . 1>/dev/null
		#ssh $PROXY $USER@$COMPUTE rm -r $remotedir/$date-RosettaAutoRelax
		if [ $( wc -l < score.sc ) -gt 2 ] ; then 
			echo -en "\rsubmitted relax to $COMPUTE, job successful              "
			if [ $DEBUG != True ] ; then
				rm *.tx 
			fi
		else
			echo -en "\rsubmitted relax to $COMPUTE, job failed              "
			exit 0
		fi
		echo ''

	fi
	
elif [ $MUTATE = True ] ; then
	cp ../*.pdb .
fi


########################################################################################################################################################################
#      BODY MUTATE     #
########################################################################################################################################################################
#echo checka
if [ $MUTATE = True ] ; then

	rm -r $oridir/rosetta-2-mutate 2>/dev/null 
	mkdir $oridir/rosetta-2-mutate
	cd $oridir/rosetta-2-mutate
	
	echo '2.0 Mutation'
################################################	
########################  INPUT PREPARATION
	echo "2.1 - Preparing mutation inputs"
	
####### 
####### 
####### CHECKING INPUT STURCUTURES
	# if an input pdb is not provided (default option is 'False'), generate an input list from the best relaxed structures
	if [ $MIPDB = False ] ; then 
		cat $oridir/rosetta-1-relax/score.sc | sed '1,2d' | awk '{print $22,$2}' | sort -n -k2 | tac | tail -$RONSS | awk '{print $1}' | sed 's/$/.pdb/g' > mipdb-auto.tx
		for i in $(cat mipdb-auto.tx) ; do 
			cp $oridir/rosetta-1-relax/$i .
		done
	else
		cp $oridir/$MIPDB mipdb-auto.tx
		for i in $(cat mipdb-auto.tx) ; do 
			cp $oridir/$i .
		done

	fi

####### 
####### 	
####### CHECKING RESFILE
	# if an res file is not provided (optional), generate an res file from the MIRES comments
	if [ $MIRES = False ] && [ $MITHD = False ] ; then 
		echo '	generating resfile'
		rm mires-auto.tx 2>/dev/null 
		# print the header of the res file
		echo 'NATRO
start' > mires-auto.tx
		# generate mutable selection
		for chain in $(cat $oridir/mires.tx | sed '/^[[:space:]]*$/d' | awk '{print $1}') ; do
			for res in $(grep $chain $oridir/mires.tx | cut -c 2-) ; do 
				echo "$res $chain $PICKACID" | sed s'/-/ /g' >> mires-auto.tx
			done
		done
		
####### generate movable selection - using pymol
		respdb=$(head -1 mipdb-auto.tx)
		if [ $RELAX = True ] ; then 
			cp $oridir/rosetta-1-relax/$respdb $oridir/$respdb
		fi
		echo "load $oridir/$respdb" > mires.pml
		sed '/^[[:space:]]*$/d' $oridir/mires.tx | sed "s/ /+/g;s/./ \& i. /2;s/^/+ c. /g" | tr '\n' ' ' | sed "s/^+ /select muts, /g
" >> mires.pml
		echo "
select shell, br. all near $flexsc of muts
select shelled, br. all near $flexbb of shell
select combined, muts + shell + shelled
save shell.pdb, shell
save combined.pdb, combined" >> mires.pml
		nohup pymol -qc mires.pml 2>/dev/null 
		wait
		cat shell.pdb | grep ' CA ' | awk '{print $6,$5}' | sed 's/$/ NATAA/g' >> mires-auto.tx
		rm shell.pdb
		cp mires-auto.tx mires-auto-mutate.tx
		MIRES=mires-auto.tx
####### generate movemap
		echo 'RESIDUE * NO' > move.map
		grep 'ATOM\|HETATM' $oridir/$respdb | cut -c 22-26 | uniq | nl > resnums.tx
		grep 'ATOM\|HETATM' combined.pdb | cut -c 22-26 | uniq > extnums.tx
		grep -f extnums.tx resnums.tx | awk '{print $1}' | sed 's/^/RESIDUE /g;s/$/ BBCHI/g' >> move.map
		echo "JUMP * $subjump" >> move.map
		
	fi

####### 
####### 
####### CHECKING SEQUENCE THREADING
	if [ $MITHD != False ] ; then
		respdb=$(sort mipdb-auto.tx | head -1)
		# print the header of the res file
		echo '	initiating threading'
		echo 'NATRO
start' > mires-auto.tx

####### check if the source file is a fasta or pdb
		if [ $( echo $MITHD | awk -F'.' '{print $2}' ) = fa ] ; then 
			cp $oridir/$MITHD seqgen.fa 
			for mltthdseq in $(grep '>' seqgen.fa | sed 's/>//g' ) ; do 
				# transpose horizontal source sequence to vertical sequence
				grep -A1 $mltthdseq seqgen.fa | sed -n '2p' | sed 's/\(.\{1\}\)/\1 /g' > transin.tx
				transpose transin.tx > srcseq-$mltthdseq.tx	
			done
		# if it is a pdb, generate a fasta file of the appropriate chain with seqgen function 
		elif [ $( echo $MITHD | awk -F'.' '{print $2}' ) = pdb ] ; then 
			mltthdseq=$( echo $MITHD | awk -F'.' '{print $1}' )
			rm seqgen.fa 2>/dev/null 
			seqgen $oridir/$MITHD $MITHC
			grep -A1 $mltthdseq seqgen.fa | sed -n '2p' seqgen.fa | sed 's/\(.\{1\}\)/\1 /g' > transin.tx
			transpose transin.tx > srcseq-$mltthdseq.tx	
		fi

		
####### generate a fasta file from the first of the relaxed structures as a template
		rm seqgen.fa 
		seqgen $oridir/$respdb $MITHC
		# transpose horizontal target sequence to vertical sequence
		sed -n '2p' seqgen.fa | sed 's/\(.\{1\}\)/\1 /g' > transin.tx
		transpose transin.tx > replaced-sequence.tx
		#mv transout.tx replaced-sequence.tx
		
####### generate residue numbering from the target structure
		cat $oridir/$respdb | grep ^ATOM | grep ' CA ' | grep "[A-Z][A-Z][A-Z] $MITHC" | awk '{print $6}' > number-sequence.tx
		paste number-sequence.tx replaced-sequence.tx > numrep.tx
		for mltthdseq in $(ls srcseq-*.tx | sed 's/srcseq-//g;s/.tx//g') ; do 
			paste number-sequence.tx srcseq-$mltthdseq.tx > nupsor.tx
			cp mires-auto.tx mires-auto-$mltthdseq.tx
			# use sdiff to compare vertical sequences and print mismatches to the res file
			sdiff numrep.tx nupsor.tx | grep '|' | awk '{print $4,$5}' | sed "s/ / $MITHC PIKAA /g" >> mires-auto-$mltthdseq.tx
			# check the number of positions and sequences for discontinuity
			ns=$(wc -l < number-sequence.tx)
			nr=$(wc -l < numrep.tx)
			np=$(wc -l < nupsor.tx)

			if [ "$ns" != "$nr" ] || [ "$ns" != "$np" ] ; then 
				echo " please check the source and template for threading. We have found different number of residues between selection and cannot thread.
				number of positions = $ns, target AA = $nr, source AA = $np"
				echo " exiting . . . " ; exit 0
			fi 
		
####### 	generate movable selection - using pymol
			echo "load $oridir/$respdb" > mires.pml
			sdiff numrep.tx nupsor.tx | grep '|' | awk '{print $4,"+"}' | tr '\n' ' ' | sed 's/ //g' | sed "s/^/select muts, c. $MITHC \& i. /g" | rev | cut -c2- | rev >> mires.pml
			echo "
select shell, br. all near $flexsc of muts
select shelled, br. all near $flexbb of shell
select combined, muts + shell + shelled
save shell.pdb, shell
save combined.pdb, combined" >> mires.pml
			nohup pymol -qc mires.pml 2>/dev/null 
			wait
			cat shell.pdb | grep ' CA ' | awk '{print $6,$5}' | sed 's/$/ NATAA/g' >> mires-auto-$mltthdseq.tx
			#cat mires-auto.tx
			rm shell.pdb mires.pml
			####### generate movemap
			echo 'RESIDUE * NO' > move.map
			grep 'ATOM\|HETATM' $oridir/$respdb | cut -c 22-26 | uniq | nl > resnums.tx
			grep 'ATOM\|HETATM' combined.pdb | cut -c 22-26 | uniq > extnums.tx
			grep -f extnums.tx resnums.tx | awk '{print $1}' | sed 's/^/RESIDUE /g;s/$/ BBCHI/g' >> move.map
			echo "JUMP * $subjump" >> move.map
		done
	fi

################################################
########################  CPU CONTROL
	echo "2.2 - Making mutations"
	
		# For local compute	
	if [ $COMPUTE = LOCAL ] ; then 
		
		# if there are more CPU's than input files, add additional copies of the input to the file list
		# this is a shite way of using multiple CPUs, please improve.
		if [ $CPU -gt $(wc -l < mipdb-auto.tx) ] ; then 
			# we currently only accept single models or more models than the number of CPU's
			if [ $(wc -l < mipdb-auto.tx) != 1 ] ; then
				echo "this script does not currently accept more CPU's than starting models (unless single model)"
				exit 0
			fi
			
			# for a single model input with mutltiple CPU, make a new line in the pdb list per CPU
			rm mipdb.txt 2>/dev/null 
			count=1
			for i in $(seq 2 $CPU) ; do 
				echo $count_$(head -1 mipdb-auto.tx) >> mipdb-auto.tx
				count=$(( $count + 1)) 
			done
			
		fi
		MIPDBL=mipdb-auto.tx
		sort mipdb-auto.tx | uniq > $oridir/mipdb-auto.tx
		
		for i in $(cat $MIPDBL ) ; do 
			resifix $i
		done
	################################################	
	########################  RUN MUTATION
		# process the input list for the number of CPUs
		
		for MIRES in $(ls mires-auto-*) ; do 															# <---- FIX THIS BETWEEN SCRIPTS
			MIRESname=$(echo $MIRES | sed 's/mires-auto-//g;s/mires-auto//g;s/.tx//g;s/fix-//g')		# <---- FIX THIS BETWEEN SCRIPTS
			for i in $(seq 1 $CPU) ; do
				#make a temp directory per cpu
				mkdir temp$i 2>/dev/null 
				#make a flag file per cpu
				#echo " cpu $i is making $MIRES"
				echo "
-out:path:all temp$i
-out:pdb
-out:suffix _$MIRESname
-out:file:scorefile score.sc
-relax:fast
$(cat $oridir/glycan.tx 2>/dev/null)
-database /usr/local/rosetta/main/database
-nstruct $(( $MONST / $CPU ))
-in:file:movemap move.map
-relax:respect_resfile 
-packing:resfile $MIRES                                                                                                                                                              
-score:weights ref2015" > mflag$i.file
				#make a rosetta executable per cpu and run them
				echo '
				for j in $(seq '$i' '$CPU' $(wc -l < '$MIPDBL') ) ; do
				pdb=$(sed -n '"''"'$j'"'"'p'"'"' '$MIPDBL')

				nohup relax.'$version'.linuxgccrelease -in:file:s $pdb @mflag'$i'.file 2>/dev/null 
				done' > run-cpu$i
				chmod 755 run-cpu$i
				./run-cpu$i &
			done
			wait
		done
	################################################	
	########################  MERGE OUTPUTS	
		echo '2.3 - consolidating outputs'
		# the output pdbs and score files are in seperate directories, with the same names, lets rename and merge
		
		# as we are not starting from the same initial pdb name, hopefully this will be simple
		rm sfcheck.tx 2>/dev/null 
		for i in $(seq 1 $CPU) ; do
			sed '1,2d' $oridir/rosetta-2-mutate/temp$i/score.sc | awk '{print $22}' | rev | cut -c 6- | rev | sort | uniq >> sfcheck.tx
		done
		sfchpre=$(wc -l < sfcheck.tx)
		sfchpost=$(cat sfcheck.tx | sort | uniq | wc -l )
		if [ $sfchpre = $sfchpost ] ; then
			rm *.pdb
			for i in $(seq 1 $CPU) ; do
				mv temp$i/*.pdb .
				if [ $i = 1 ] ; then 
					#move everything from the first temp directory to the destination directory
					mv temp$i/score.sc .
				else
				cat temp$i/score.sc | sed '1,2d' >> score.sc
				fi
			done
		#rm -r temp* run-cpu* mflag*
		else # this is a bit shit but the names should be unique (although naming problems will exist if over 9999 models total cat
			rm *.pdb
			for i in $(seq 1 $CPU) ; do
				#count files in destination directory
				pdbcount=$(ls *.pdb 2>/dev/null | wc -l ) 
				if [ $i = 1 ] ; then 
					#move the header of the score file to the destination directory
					head -2 $oridir/rosetta-2-mutate/temp1/score.sc > $oridir/rosetta-2-mutate/score.sc
				fi 
				#move everything from the temp directorys to the destination directory
				pdbcount2=0
				for j in $(ls $oridir/rosetta-2-mutate/temp$i/*.pdb) ; do 
					pdbcount2=$(( $pdbcount2 + 1 ))
					#define old number (with leading 0's)
					oldnum=$(echo $j | rev | cut -c 1-8 | rev | cut -c 1-4)
					#define old number (without leading 0's)
					num=$(echo $j | rev | cut -c 1-8 | rev | cut -c 1-4 | sed 's/0//g')
					#define new number (without leading 0's)
					newnum=$(( $pdbcount2 + $pdbcount))
					#define new number (with leading 0's)
					newnumname=$(echo 00000$newnum | rev | cut -c 1-6 | rev)
					#define file name without number
					nonumname=$(echo $j | rev | cut -c 9- | awk -F'/' '{print $1}' | rev)
					
					# rename and move pdb and contents of score file
					cp $oridir/rosetta-2-mutate/temp$i/$nonumname$oldnum.pdb $oridir/rosetta-2-mutate/$nonumname$newnumname.pdb
					grep $nonumname$oldnum $oridir/rosetta-2-mutate/temp$i/score.sc | sed "s/$nonumname$oldnum/$nonumname$newnumname/g" >> $oridir/rosetta-2-mutate/score.sc
					
				done
				
			done
		# clean up
		mv mires-auto.tx  mires-auto.txt
		if [ $DEBUG != True ] ; then
				rm -r temp* run-cpu* mflag* *.tx 2>/dev/null 
		fi
		fi
		for i in $(cat $MIPDBL 2>/dev/null ) ; do 
			rm $i 2>/dev/null 
		done
		ls *.pdb > mopdblist.tx  2>/dev/null 
		MOPDBL=mopdblist.tx
		mv mires-auto.txt  mires-auto.tx
		
		# mipdb location correction
		for i in $(cat $oridir/mipdb-auto.tx) ; do 
			cp $oridir/rosetta-1-relax/$i .
		done
		
	fi
		# for remote compute
	if [ $COMPUTE != LOCAL ] ; then 
		for i in $(cat mipdb-auto.tx) ; do 
			echo "2.3 - Running remote mutations for $i"
			date=$(echo $(date -R | cut -c 6-7,9-11,15-25 | sed 's/ /-/g;s/://g')$RANDOM)
			cp ../remote.tx .
			ssh $PROXY $USER@$COMPUTE mkdir $remotedir/$date-RosettaAutoMutate
			resifix $i
			scp $PROXY $oridir/rosetta-2-mutate/$i $USER@$COMPUTE:$remotedir/$date-RosettaAutoMutate/. 1>/dev/null
			scp $PROXY $oridir/rosetta-2-mutate/move.map $USER@$COMPUTE:$remotedir/$date-RosettaAutoMutate/. 1>/dev/null
			scp $PROXY $oridir/rosetta-2-mutate/mires-auto-mutate.tx $USER@$COMPUTE:$remotedir/$date-RosettaAutoMutate/res.file 1>/dev/null
			echo "-database /apps/applications/rosetta/3.10/1/default/main/database
-relax:fast
-out:suffix _"'$state'"
-out:file:scorefile score.sc
$(cat $oridir/glycan.tx 2>/dev/null )
-nstruct 1
-no_nstruct_label
-in:file:movemap move.map
-relax:respect_resfile
-packing:resfile res.file
-score:weights ref2015 \" > flag_"'$SGE_TASK_ID'".file
/apps/applications/rosetta/3.10/1/default/main/source/bin/relax.linuxgccrelease -in:file:s $i @flag_"'$SGE_TASK_ID'".file
	" >> remote.tx
			cat remote.tx | sed "s/XXX/$MONST/g;s/YYY/100/g;s#VVV#$remotedir/$date-RosettaAutoMutate#g" > run-rosetta-mutate.sh
			scp $PROXY run-rosetta-mutate.sh $USER@$COMPUTE:$remotedir/$date-RosettaAutoMutate/. 1>/dev/null
			ssh $PROXY $USER@$COMPUTE " cd $remotedir/$date-RosettaAutoMutate ; qsub run-rosetta-mutate.sh" > $date.tx 
			arcjobid=$( awk '{print $3}' $date.tx | awk -F'.' '{print $1}' )
			while [ $(ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid 2>/dev/null | wc -l ) != 0 ] ; do
				ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid | cut -c 41-43,104- > arcjobid.tx
				if [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == qw ] ; then 
					echo -en "\rsubmitted Mutate to $COMPUTE, queueing "$(awk '{print $2}' arcjobid.tx | awk -F':' '{print $1}')
				elif [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == r ] ; then 
					echo -en "\rsubmitted Mutate to $COMPUTE, running last job                   "
				elif [ $( wc -l < arcjobid.tx) -gt 1 ] ; then 
					echo -en "\rSubmitted mutate to $COMPUTE, queueing "$(grep qw arcjobid.tx | awk '{print $2}')", running "$(grep r arcjobid.tx | wc -l )"        "
				fi
				sleep 20
			done
                        scp $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoMutate/score.sc scoret.tx 1>/dev/null
                        if [ ! -f score.sc ] ; then
                                cp scoret.tx score.sc
                        else
                                sed '1,2d' scoret.tx >> score.sc
                        fi
			scp $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoMutate/*.pdb . 1>/dev/null
			#ssh $PROXY $USER@$COMPUTE rm -r $remotedir/$date-RosettaAutoMutate
			if [ $( wc -l < scoret.tx ) = $(( 2 + $MONST)) ] ; then 
				echo -en "\rsubmitted Mutate to $COMPUTE, job successful              "
			else
				echo -en "\rsubmitted Mutate to $COMPUTE, job failed              "
				exit 0
			fi
			rm scoret.tx
		echo ''
		done
		rm combined.pdb
		ls *.pdb > mopdblist.tx  2>/dev/null 
		MOPDBL=mopdblist.tx
		cp mires-auto.tx  mires-auto.txt
		cp mipdb-auto.tx $oridir/.
		

	fi
	if [ $DEBUG != True ] ; then
		rm -r temp* mflag* run-* nohup-out 2>/dev/null 
	fi
	

fi
cd $oridir

########################################################################################################################################################################
#    BODY INTERFACE    #
########################################################################################################################################################################

if [ $INTER = True ] ; then

	rm -r $oridir/rosetta-3-inter 2>/dev/null 
	mkdir $oridir/rosetta-3-inter
	cd  $oridir/rosetta-3-inter
	
	echo '3.0 Interface Analysis'
	
	# for local compute
	if [ $COMPUTE = LOCAL ] ; then 
		cp $oridir/rosetta-1-relax/*.pdb .
		cp $oridir/rosetta-2-mutate/*.pdb .
		ls -v *.pdb > pdb.list

		IIPDB=pdb.list

	################################################	
	########################  Run Interface analysis
		# process the input list for the number of CPUs
		for i in $(seq 1 $CPU) ; do
			#make a temp directory per cpu
			mkdir $oridir/rosetta-3-inter/temp$i
			#make a flag file per cpu
			echo "#specific options for InterfaceAnalyzer
-out:no_nstruct_label
-out:path:all $oridir/rosetta-3-inter/temp$i
-database /usr/local/rosetta/main/database
-fixedchains $IIFAC
-score:weights ref2015
$(cat $oridir/glycan.tx 2>/dev/null)
-compute_packstat true
-tracer_data_print false #make a score file with all the important info instead of just printing to the terminal
-out:file:score_only inter_score.sc #This will cause output of all of the info to a file called inter_score.sc
-pack_input false #will not relax the input interface residues (good if you have already minimized/packed your structure)
-pack_separated false #will also not pack the monomers to calculated dG bind.
-add_regular_scores_to_scorefile true #will run the rest of rosettas score function on your complex using score12
#these are some tweeks that we have found helpful
-atomic_burial_cutoff 0.01 #This is set to help rosetta identify buried polar atoms properly
-sasa_calculator_probe_radius 1.4 #This is the default water probe radius for SASA calculations, sometimes lowering the radius helps rosetta more accurately find buried polar atoms
-pose_metrics::interface_cutoff 8.0 # this defines how far away a CBeta atom can be from the other chain to be considered an interface residue" > iflag$i.file

			#make a rosetta executable per cpu and run them
			echo '
for j in $(seq '$i' '$CPU' $(wc -l < '$IIPDB') ) ; do
pdb=$(sed -n '"''"'$j'"'"'p'"'"' '$IIPDB')
nohup InterfaceAnalyzer.'$version'.linuxgccrelease -in:file:s $pdb @iflag'$i'.file 2>/dev/null
done' > run-cpu$i
			chmod 755 run-cpu$i
			./run-cpu$i &
		done
		wait

		for i in $(seq 1 $CPU) ; do
			if [ $i = 1 ] ; then 
				#move everything from the first temp directory to the destination directory
				cat $oridir/rosetta-3-inter/temp$i/inter_score.sc  > $oridir/rosetta-3-inter/inter_score.sc
			else
				cat $oridir/rosetta-3-inter/temp$i/inter_score.sc | sed '1,2d' >> $oridir/rosetta-3-inter/inter_score.sc
			fi
		done
		if [ $DEBUG != True ] ; then
			rm -r temp* iflag* run-* nohup-out pdb.list 2>/dev/null 
		fi
		
		# if a crash log is present.
		if [ -f ROSETTA_CRASH.log ] ; then 
			mv ROSETTA_CRASH.log roscrash.tx
			echo '	Crash report found during analysis, attempting repair'
			mkdir temp
			echo "#specific options for InterfaceAnalyzer
-out:no_nstruct_label
-out:path:all $oridir/rosetta-3-inter/temp
-database /usr/local/rosetta/main/database
-fixedchains $IIFAC
-score:weights ref2015
$(cat $oridir/glycan.tx 2>/dev/null)
-compute_packstat true
-tracer_data_print false #make a score file with all the important info instead of just printing to the terminal
-out:file:score_only inter_score.sc #This will cause output of all of the info to a file called inter_score.sc
-pack_input false #will not relax the input interface residues (good if you have already minimized/packed your structure)
-pack_separated false #will also not pack the monomers to calculated dG bind.
-add_regular_scores_to_scorefile true #will run the rest of rosettas score function on your complex using score12
#these are some tweeks that we have found helpful
-atomic_burial_cutoff 0.01 #This is set to help rosetta identify buried polar atoms properly
-sasa_calculator_probe_radius 1.4 #This is the default water probe radius for SASA calculations, sometimes lowering the radius helps rosetta more accurately find buried polar atoms
-pose_metrics::interface_cutoff 8.0 # this defines how far away a CBeta atom can be from the other chain to be considered an interface residue" > iflag.file

			for i in $(cat roscrash.tx | grep file | awk '{print $8}' | rev | awk -F '/' '{print $1}' | rev ) ; do 
				nohup InterfaceAnalyzer.$version.linuxgccrelease -in:file:s $oridir/rosetta-2-mutate/$i.pdb @iflag.file 2>/dev/null
			done
			cat $oridir/rosetta-3-inter/temp/inter_score.sc | sed '1,2d' >> $oridir/rosetta-3-inter/inter_score.sc
			
		fi
	fi
	
	# if remote compute
	if [ $COMPUTE != LOCAL ] ; then 
		echo "3.1 - Running remote interface analysis"
		date=$(echo $(date -R | cut -c 6-7,9-11,15-25 | sed 's/ /-/g;s/://g')$RANDOM)
		cp $oridir/rosetta-1-relax/*.pdb .
		cp $oridir/rosetta-2-mutate/*.pdb .
		ls -v *.pdb > pdb.list
		cp ../remote.tx .
		ssh $PROXY $USER@$COMPUTE mkdir $remotedir/$date-RosettaAutoInter
		
		
		scp $PROXY $oridir/rosetta-2-mutate/*.pdb $USER@$COMPUTE:$remotedir/$date-RosettaAutoInter/. 1>/dev/null
		scp $PROXY $oridir/rosetta-3-inter/pdb.list $USER@$COMPUTE:$remotedir/$date-RosettaAutoInter/. 1>/dev/null
		echo "-database /apps/applications/rosetta/3.10/1/default/main/database
-out:no_nstruct_label
-in:file:s "'$(sed -n ""$SGE_TASK_ID"p"'" pdb.list)
-fixedchains $IIFAC
-score:weights ref2015
$(cat $oridir/glycan.tx 2>/dev/null)
-compute_packstat true
-tracer_data_print false #make a score file with all the important info instead of just printing to the terminal
-out:file:score_only inter_score.sc #This will cause output of all of the info to a file called inter_score.sc
-pack_input false #will not relax the input interface residues (good if you have already minimized/packed your structure)
-pack_separated false #will also not pack the monomers to calculated dG bind.
-add_regular_scores_to_scorefile true #will run the rest of rosettas score function on your complex using score12
#these are some tweeks that we have found helpful
-atomic_burial_cutoff 0.01 #This is set to help rosetta identify buried polar atoms properly
-sasa_calculator_probe_radius 1.4 #This is the default water probe radius for SASA calculations, sometimes lowering the radius helps rosetta more accurately find buried polar atoms
-pose_metrics::interface_cutoff 8.0 # this defines how far away a CBeta atom can be from the other chain to be considered an interface residue\" > flag_"'$SGE_TASK_ID'".file

/apps/applications/rosetta/3.10/1/default/main/source/bin/InterfaceAnalyzer.linuxgccrelease @flag_"'$SGE_TASK_ID'".file
" >> remote.tx
		cat remote.tx | sed "s/XXX/$(wc -l < pdb.list)/g;s/YYY/100/g;s#VVV#$remotedir/$date-RosettaAutoInter#g" > run-rosetta-inter.sh
		scp $PROXY run-rosetta-inter.sh $USER@$COMPUTE:$remotedir/$date-RosettaAutoInter/. 1>/dev/null
		ssh $PROXY $USER@$COMPUTE " cd $remotedir/$date-RosettaAutoInter ; qsub run-rosetta-inter.sh" > $date.tx 
		arcjobid=$( awk '{print $3}' $date.tx | awk -F'.' '{print $1}' )
		while [ $(ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid 2>/dev/null | wc -l ) != 0 ] ; do
			ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid | cut -c 41-43,104- > arcjobid.tx
			if [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == qw ] ; then 
				echo -en "\rsubmitted InterfaceAnalyzer to $COMPUTE, queueing "$(awk '{print $2}' arcjobid.tx | awk -F':' '{print $1}')
			elif [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == r ] ; then 
				echo -en "\rsubmitted InterfaceAnalyzer to $COMPUTE, running last job                   "
			elif [ $( wc -l < arcjobid.tx) -gt 1 ] ; then 
				echo -en "\rsubmitted InterfaceAnalyzer to $COMPUTE, queueing "$(grep qw arcjobid.tx | awk '{print $2}')", running "$(grep r arcjobid.tx | wc -l )"        "
			fi
			sleep 30
		done
		scp $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoInter/inter_score.sc . 1>/dev/null
		ssh $PROXY $USER@$COMPUTE rm -r $remotedir/$date-RosettaAutoInter
		if [ $( wc -l < inter_score.sc ) -gt 2 ] ; then 
			echo -en "\rsubmitted InterfaceAnalyzer to $COMPUTE, job successful              "
		else
			echo -en "\rsubmitted InterfaceAnalyzer to $COMPUTE, job failed              "
			exit 0
		fi
		echo ''
	fi
	if [ $DEBUG != True ] ; then
		rm -r temp* iflag* run-* nohup-out 2>/dev/null 
	fi
	

fi
cd  $oridir

########################################################################################################################################################################
#    Body Analysis    #
########################################################################################################################################################################
echo "4.0 Analysis"

rm -r $oridir/rosetta-4-analysis 2>/dev/null 
mkdir $oridir/rosetta-4-analysis
cd  $oridir/rosetta-4-analysis

cenperes=False
if [ $ROPRS = True ] ; then 
enperes=1-relax
cenperes=True
fi 
if [ $MOPRS = True ] ; then 
enperes="$enperes 2-mutate"
cenperes=True
fi 
if [ $IOPRS = True ] ; then 
enperes="$enperes 3-inter"
cenperes=True
fi 

if [ $cenperes = True ] ; then 

	mkdir pdb
	for i in $enperes ; do 
		cp $oridir/rosetta-$i/*.pdb pdb/.
	done

	cd  pdb
	ls *.pdb > ../pdb.list
	cd  $oridir/rosetta-4-analysis

	# For local computation
	if [ $IOPRS = True ] && [ $COMPUTE = LOCAL ] ; then 
		echo "4.1 - Running energy breakdown	analysis"
		# process the input list for the number of CPUs
		for i in $(seq 1 $CPU) ; do
			#make a temp directory per cpu
			mkdir $oridir/rosetta-4-analysis/temp$i
			cd $oridir/rosetta-4-analysis/temp$i
			#make a flag file per cpu
			echo "#specific options for per-res
-database /usr/local/rosetta/main/database
-score:weights ref2015
$(cat $oridir/glycan.tx 2>/dev/null)
-out:file:score_only perres_score.sc" > rflag$i.file

			#make a rosetta executable per cpu and run them
			echo '
for j in $(seq '$i' '$CPU' $(wc -l < '../pdb.list') ) ; do
pdb=$(sed -n '"''"'$j'"'"'p'"'"' '../pdb.list')
nohup residue_energy_breakdown.'$version'.linuxgccrelease -in:file:s ../pdb/$pdb @rflag'$i'.file 2>/dev/null
done' > run-cpu$i
			chmod 755 run-cpu$i
			./run-cpu$i &
		done
		wait
		
		cd $oridir/rosetta-4-analysis
		for i in $(seq 1 $CPU) ; do
			if [ $CPU = 1 ] ; then
				mv temp$i/default.out perres.out
			else	
				cat temp$i/default.out >> perres.out	
			fi
			
		done
				
	fi

	# for remote compute
	if [ $IOPRS = True ] && [ $COMPUTE != LOCAL ] ; then 
		echo "4.1 - Running remote energy breakdown	analysis"
		date=$(echo $(date -R | cut -c 6-7,9-11,15-25 | sed 's/ /-/g;s/://g')$RANDOM)
		cp ../remote.tx .
		ssh $PROXY $USER@$COMPUTE mkdir $remotedir/$date-RosettaAutoPer

		scp $PROXY $oridir/rosetta-4-analysis/pdb/*.pdb $USER@$COMPUTE:$remotedir/$date-RosettaAutoPer/. 1>/dev/null
		scp $PROXY $oridir/rosetta-4-analysis/pdb.list $USER@$COMPUTE:$remotedir/$date-RosettaAutoPer/. 1>/dev/null		
		echo "-database /apps/applications/rosetta/3.10/1/default/main/database
-in:file:s "'$(sed -n ""$SGE_TASK_ID"p"'" pdb.list)
-score:weights ref2015 \" > flag_"'$SGE_TASK_ID'".file

/apps/applications/rosetta/3.10/1/default/main/source/bin/residue_energy_breakdown.linuxgccrelease @flag_"'$SGE_TASK_ID'".file
cat default.out >> perres.out ; rm default.out
	" >> remote.tx
		cat remote.tx | sed "s/XXX/$(wc -l < pdb.list)/g;s/YYY/100/g;s#VVV#$remotedir/$date-RosettaAutoPer#g" > run-rosetta-perres.sh
		scp $PROXY run-rosetta-perres.sh $USER@$COMPUTE:$remotedir/$date-RosettaAutoPer/. 1>/dev/null
		ssh $PROXY $USER@$COMPUTE " cd $remotedir/$date-RosettaAutoPer ; qsub run-rosetta-perres.sh" > $date.tx 
		arcjobid=$( awk '{print $3}' $date.tx | awk -F'.' '{print $1}' )
		while [ $(ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid 2>/dev/null | wc -l ) != 0 ] ; do
			ssh $PROXY $USER@$COMPUTE qstat | grep $arcjobid | cut -c 41-43,104- > arcjobid.tx
			if [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == qw ] ; then 
				echo -en "\rsubmitted energy breakdown analyze to $COMPUTE, queueing "$(awk '{print $2}' arcjobid.tx | awk -F':' '{print $1}')
			elif [ $( wc -l < arcjobid.tx) = 1 ] && [ $(awk '{print $1}' arcjobid.tx) == r ] ; then 
				echo -en "\rsubmitted energy breakdown analyze to $COMPUTE, running last job                   "
			elif [ $( wc -l < arcjobid.tx) -gt 1 ] ; then 
				echo -en "\rsubmitted energy breakdown analyze to $COMPUTE, queueing "$(grep qw arcjobid.tx | awk '{print $2}')", running "$(grep r arcjobid.tx | wc -l )"        "
			fi
			sleep 20
		done
		#echo "downloading energy breakdown, of size $(ssh $PROXY $USER@$COMPUTE ls -lh $remotedir/$date-RosettaAutoPer | grep perres.out | awk '{print $5}')"
		ssh $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoPer "sort perres.out | uniq > perres2.out ; mv perres2.out perres.out"
		scp $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoPer/perres.out .
		sort perres.out | uniq > perres2.out ; mv perres2.out perres.out
		ssh $PROXY $USER@$COMPUTE rm -r $remotedir/$date-RosettaAutoPer
		if [ $( wc -l < perres.out ) -gt 2 ] ; then 
			echo -en "\rsubmitted energy breakdown analyzer to $COMPUTE, job successful              "
		else
			echo -en "\rsubmitted energy breakdown analyzer to $COMPUTE, job failed              "
			exit 0
		fi
		echo ''
	fi	
	#cp $oridir/rosetta-4-analysis/pdb/perres.out $oridir/rosetta-4-analysis/perres.out
	cd $oridir/rosetta-4-analysis ; rm -r temp*
fi

echo '4.2 - collating sequence and energy data'
################################################	
########################  Generate total fasta sequence files. 
if [ $MOSEQ = True ] || [ $MOMSQ = True ] ; then 
	cd  $oridir/rosetta-4-analysis
fi

if [ $RELAX = False ] ; then 
	RIPDB=$(head -1 ../mipdb-auto.tx)
fi


if [ $MOSEQ = True ] ; then 
	#echo 'all chain initial'
	rm seqgen.fa  2>/dev/null 
	name=$(echo $RIPDB | sed 's/.pdb//g' )
	#echo 'file = '$RIPDB
	seqgen $oridir/$RIPDB
	#sed -i "s/$RIPDB/$name/g" seqgen.fa

	for i in $(cat $oridir/rosetta-2-mutate/$MOPDBL) ; do 
		#echo 'all chain mutant'
		name=$(echo $i | sed 's/.pdb//g' )
		seqgen $oridir/rosetta-2-mutate/$i 
		#sed -i "s/$i/$name/g" seqgen.fa
	
	done
	mv seqgen.fa rosetta-full-sequence.fa
	# generate mutant only fasta file for threading
for chain in $(cat rosetta-full-sequence.fa | grep '>' | rev | cut -c 1 | rev | sort | uniq) ; do
grep -A1 "_$chain$" rosetta-full-sequence.fa | grep -A1 mutate | grep -v '^--' > rosetta-thread-$chain-sequence.fa
kpvar=$(grep -A1 "_$chain$" rosetta-full-sequence.fa | grep -A1 mutate | grep -v '^--' | head -1 | awk -F'/' '{print $NF}')
rmvar=$(grep -A1 "_$chain$" rosetta-full-sequence.fa | grep -A1 mutate | grep -v '^--' | head -1 | sed "s/$kpvar//g" | cut -c 2-)
sed -i "s#$rmvar##g" rosetta-thread-$chain-sequence.fa
done
fi 

################################################	
########################  Generate specific fasta sequence files with total energy and interface energy (if appropriate)
if [ $MOMSQ = True ] && [ $MITHD = False ] ; then 
	#echo 'RIPDB = '$RIPDB
	#echo 'spec chain initial'
	rm seqgen.fa specseqpdb.tx 2>/dev/null 
	cat $oridir/rosetta-2-mutate/$MIRES | sed '1,2d' | grep -v 'NATAA$' | grep -v 'NATRO' > specseqgen.tx
	
	for i in $(seq 1 $(wc -l < specseqgen.tx) ) ; do 
		#sed -n ''$i'p' specseqgen.tx
		specseqpos=$(echo "---$(sed -n ''$i'p' specseqgen.tx | awk '{print $1}')" | rev | cut -c 1-4 | rev)
		specseqcha=$(sed -n ''$i'p' specseqgen.tx | awk '{print $2}')
		specseqpch=$(echo $specseqcha$specseqpos | sed 's/-/ /g')
		grep ' CA ' $oridir/$RIPDB | grep "$specseqpch" >> specseqpdb.tx
	done
	name=$(echo $RIPDB | sed 's/.pdb//g' )
	mv specseqpdb.tx $name
	#echo 'name ='$name 'chain ='$specseqcha
	seqgen $name $specseqcha
	rm $name

	for j in $(cat $oridir/rosetta-2-mutate/$MOPDBL) ; do 
		#echo 'spec chain mutant'
		for i in $(seq 1 $(wc -l < specseqgen.tx) ) ; do 
			#sed -n ''$i'p' specseqgen.tx
			specseqpos=$(echo "---$(sed -n ''$i'p' specseqgen.tx | awk '{print $1}')" | rev | cut -c 1-4 | rev)
			specseqcha=$(sed -n ''$i'p' specseqgen.tx | awk '{print $2}')
			specseqpch=$(echo $specseqcha$specseqpos | sed 's/-/ /g')
			grep ' CA ' $oridir/rosetta-2-mutate/$j | grep "$specseqpch" >> specseqpdb.tx
		done
		name=$(echo $j | sed 's/\// /g' | rev | cut -d' ' -f 1 | rev | sed 's/.pdb//g' )
		mv specseqpdb.tx $name
		seqgen $name $specseqcha
		rm $name
	done
	mv seqgen.fa rosetta-spec-sequence.fa
fi
if [ $MOMSQ = True ] && [ $MITHD != False ] ; then 
	rm seqgen.fa specseqpdb.tx 2>/dev/null 
	specseqcha=$(cut -c 1 $oridir/mires.tx | sed -n '2p')
	cut -c 3- $oridir/mires.tx | sed -n '2p'| tr ' ' '\n' | sed '$d' > mires-col.tx 
	for i in $(cat mires-col.tx ) ; do 
		#sed -n ''$i'p' specseqgen.tx
		specseqpos=$(echo "---$i" | rev | cut -c 1-4 | rev)
		specseqpch=$(echo $specseqcha$specseqpos | sed 's/-/ /g')
		grep ' CA ' $oridir/$RIPDB | grep "$specseqpch" >> specseqpdb.tx
	done
	name=$(echo $RIPDB | sed 's/.pdb//g' )
	mv specseqpdb.tx $name
	#echo 'name ='$name 'chain ='$specseqcha
	seqgen $name $specseqcha
	rm $name

	for j in $(cat $oridir/rosetta-2-mutate/$MOPDBL) ; do 
		#echo 'spec chain mutant'
		specseqcha=$(cut -c 1 $oridir/mires.tx | sed -n '2p')
		for i in $(cat mires-col.tx ) ; do 
			#sed -n ''$i'p' specseqgen.tx
			specseqpos=$(echo "---$i" | rev | cut -c 1-4 | rev)
			specseqpch=$(echo $specseqcha$specseqpos | sed 's/-/ /g')
			grep ' CA ' $oridir/rosetta-2-mutate/$j | grep "$specseqpch" >> specseqpdb.tx
		done
		name=$(echo $j | sed 's/\// /g' | rev | cut -d' ' -f 1 | rev | sed 's/.pdb//g' )
		mv specseqpdb.tx $name
		seqgen $name $specseqcha
		rm $name
	done
	mv seqgen.fa rosetta-spec-sequence.fa
fi

cd $oridir/rosetta-4-analysis

#prepare output for the relaxations
if [ $RELAX = True ] ; then
	echo 'name,total-REU,Sequence' > relaxed-output.csv
	if [ $INTER = True ] ; then 
		echo 'name,total-REU,Interface-REU,Sequence' > relaxed-output.csv
	fi
	for roname in $(sed '1,2d' $oridir/rosetta-1-relax/score.sc | awk '{print $22}') ; do 
		rotote=$(grep $roname$ $oridir/rosetta-1-relax/score.sc | awk '{print $2}')

		specseqcha=$(cut -c 1 $oridir/mires.tx | sed -n '2p')
		rosseq=$(sed -n '2p' rosetta-spec-sequence.fa)
		if [ $INTER = True ] ; then 
			rointe=$(grep $roname$ $oridir/rosetta-3-inter/inter_score.sc | awk '{print $6}')
			echo "$roname,$rotote,$rointe,$rosseq" >> relaxed-output.csv
		else
			echo "$roname,$rotote,$rosseq" >> relaxed-output.csv
		fi
	done
fi

#cat $oridir/mipdb-auto.tx

#prepare output for the mutations
if [ $MUTATE = True ] ; then
	rm mutate-output.csv 2>/dev/null
	echo 'name,total-REU,total-DDG,Sequence' > mutate-outputh.csv
	if [ $INTER = True ] ; then 
		echo 'name,total-REU,total-DDG,Interface-REU,inter-DDG,Sequence' > mutate-outputh.csv
	fi
	for i in $(cat $oridir/mipdb-auto.tx | sed 's/.pdb//g') ; do 
		if [ $RELAX = True ] ; then 
			vardir=rosetta-1-relax
		else
			vardir=rosetta-0-initial
		fi
		morote=$(grep $i$ $oridir/$vardir/score.sc | awk '{print $2}')
		motddg=0
		specseqcha=$(cut -c 1 $oridir/mires.tx | sed -n '2p')
		RIPDBn=$(echo $RIPDB | sed 's/.pdb//g')
		mosseq=$(grep -A1 "$RIPDBn\_$specseqcha$" rosetta-spec-sequence.fa | sed -n '2p' )
		if [ $INTER = True ] ; then 
			mointe=$(grep $i$ $oridir/rosetta-3-inter/inter_score.sc | awk '{print $6}')
			moiddg=0
			echo "$i,$morote,$motddg,$mointe,$moiddg,$mosseq" >> mutate-output.csv
		else
			echo "$i,$morote,$motddg,$mosseq" >> mutate-output.csv
		fi
		
		for moname in $(sed '1,2d' $oridir/rosetta-2-mutate/score.sc | grep $i | awk '{print $22}') ; do 
			motote=$(grep $moname$ $oridir/rosetta-2-mutate/score.sc | awk '{print $2}')
			motddg=$( echo "scale=3; $motote - $morote" | bc  )
			specseqcha=$(cut -c 1 $oridir/mires.tx | sed -n '2p')
			mosseq=$(grep -A1 "$moname\_$specseqcha$" rosetta-spec-sequence.fa | sed -n '2p' )
			if [ $INTER = True ] ; then 
				mointe=$(grep $moname$ $oridir/rosetta-3-inter/inter_score.sc | awk '{print $6}')
				moroin=$(grep $i$ $oridir/rosetta-3-inter/inter_score.sc | awk '{print $6}')
				moiddg=$( echo "scale=3; $mointe - $moroin" | bc  )
				echo "$moname,$motote,$motddg,$mointe,$moiddg,$mosseq" >> mutate-output.csv
			else
				echo "$moname,$motote,$motddg,$mosseq" >> mutate-output.csv
			fi
			
		done
	done
	sort -t',' -k1,1 mutate-output.csv > temp
	cat mutate-outputh.csv temp > mutate-output.csv
fi
if [ $DEBUG != True ] ; then
	rm *.tx
fi

########################################################################################################################################################################
#    Clustering    #
########################################################################################################################################################################



echo 5.0 Clustering
if [ $action = rmsd ] ; then
	rmsd=no						# yes	<--- needs fixing for rmsd matrix
	sequence=no							 		
	matrix="rmsd"
elif [ $action = sequence ] ; then 
	rmsd=no							 		
	sequence=yes
	matrix="sequence"
elif [ $action = both ] ; then 
	rmsd=no						# yes	<--- needs fixing for rmsd matrix
	sequence=yes
	matrix="rmsd sequence"
elif [ $action = False ] ; then 
	rmsd=no
	sequence=no
fi 
if [ $action != False ] ; then
	rm -r $oridir/rosetta-5-cluster 2>/dev/null 
	mkdir $oridir/rosetta-5-cluster
	cd  $oridir/rosetta-5-cluster

cp $oridir/rosetta-4-analysis/mutate-output.csv .
cp $oridir/rosetta-4-analysis/rosetta-spec-sequence.fa .

# define energy col if mutate 
if [ $MUTATE = True ] ; then 
energycol=total-REU
fi

# define energy col if inter
if [ $INTER = True ] ; then 
energycol=Interface-REU
fi 

###########################
## GENERATE RMSD SIMILARITY MATRIX
if [ $rmsd = yes ] ; then 
	ls -v *.pdb > pdb.list
	if [ $(wc -l < pdb.list) -gt 2000 ] ; then 
		echo "too many pdbs for rms matrix!"
		rmsd=no
	else
		rmsmatrix 1>/dev/null 
		mv matrix.csv rmsd-matrix.csv
	fi 
fi 

###########################
## GENERATE SEQUENCE SIMILARITY MATRIX AND ENERGY OVERLAY
if [ $sequence = yes ] ; then	
echo "
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import euclidean_distances
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.cluster import hierarchy


codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
         'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def one_hot_encode(seq):
    o = list(set(codes) - set(seq))
    s = pd.DataFrame(list(seq))    
    x = pd.DataFrame(np.zeros((len(seq),len(o)),dtype=int),columns=o)    
    a = s[0].str.get_dummies(sep=',')
    a = a.join(x)
    a = a.sort_index(axis=1)
    e = a.values.flatten()
    return e

def encode():
    mutseq = pd.read_csv('mutate-output.csv')
    
    sequences = mutseq['Sequence']
    
    encoded = np.array([one_hot_encode(seq) for seq in sequences])
    
    np.save('encoded_sequenced.npy', encoded)
    return mutseq 
'''
def cluster(distances, n_clusters):
    
    model = AgglomerativeClustering(linkage='average', n_clusters=n_clusters)
    
    model.fit(distances)
    
    return model.labels_    
'''    
df=encode()

X = np.load('encoded_sequenced.npy')
distances = euclidean_distances(X, X)

plt.figure(figsize = (12,10))
sns_heat = sns.heatmap(distances)
figure = sns_heat.get_figure()
figure.savefig('heat-native.png', dpi=800)
plt.clf()

row_linkage = hierarchy.linkage(distances, method='average')

cutree = [1, 3, 9, 27, 81, 243, 729]
my_array = []
for i in cutree :
    clusters = hierarchy.fcluster(row_linkage, i, 'maxclust')
    my_array.append(clusters)


my_array = np.array(my_array)

# Get sorted indices
sorted_indices = my_array[-1,:].argsort()

# Index distances with sorted indices
distances = distances[sorted_indices]

# Transpose and sort again 
distances = distances.T[sorted_indices]


np.savetxt('heat-clusters.csv', distances , delimiter=',',  fmt='%s')

plt.figure(figsize = (12,10))
sns_heat = sns.heatmap(distances)
figure = sns_heat.get_figure()
figure.savefig('heat-cluster.png', dpi=800) 

#print(my_array)

for i, n in enumerate(cutree):
    
    df[f'n{n}'] = my_array[i,:]
    df[f'e{n}'] = 0
    
    for j in range(1, n+1):
        
        mean = df['$energycol'].loc[df[f'n{n}'] == j].mean()
        df[f'e{n}'].loc[df[f'n{n}'] == j] = mean

df['mean_energy'] = 0

cluster_columns = [f'n{n}' for n in cutree]
energy_columns = [f'e{n}' for n in cutree]

cluster_columns = df[cluster_columns].to_numpy()
energy_columns = df[energy_columns].to_numpy()
 

matrix = []
for i, x in enumerate(cluster_columns):
    matrix_row = []
    for j, y in enumerate(cluster_columns):
        matrix_row.append(energy_columns[i, np.where(x == y)[0][-1]])
    matrix.append(matrix_row)

matrix = np.array(matrix)

# Index distances with sorted indices
matrix = matrix[sorted_indices]

# Transpose and sort again 
matrix = matrix.T[sorted_indices]


plt.figure(figsize = (12,10))
sns_heat = sns.heatmap(np.array(matrix))
figure = sns_heat.get_figure()
figure.savefig('heat-average.png', dpi=800)
" > mc.py
python mc.py 2>/dev/null
fi

fi



########################################################################################################################################################################
#    Clean up    #
########################################################################################################################################################################


echo '5.0 Complete!'

