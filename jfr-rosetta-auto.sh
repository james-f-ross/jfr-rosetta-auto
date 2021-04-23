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

COMPUTE=arc4.leeds.ac.uk    # LOCAL Arc4 Arc3       # where to conduct the Rosetta calculations (remote options require '.leeds.ac.uk'
                            # if not local, it is assumed to be an arc cluster at the univerity of leeds, or another SGE cluster
                            # if not local, it is assumed you have already passed you ssh $PROXY id to the remote machine
HOURS=04:00:00              # HH:MM:SS for arc run  # this must exceed the maximum time, please wildly overestimate.
USER=chmjro                 # string                # your username on the remote machine, if required
EMAIL=chmjro@leeds.ac.uk     # email address         # your email address
remotedir=/nobackup/chmjro   # path                  # home directory on remote machine
PROXYJUMP=FALSE           # FALSE or proxy jump user and address [user]@remote-access.leeds.ac.uk

CPU=1                      # Integer               # Number of CPUs to use  (make sure RINST and MONST are equal or greater than CPUs)
                            # only for local, HPC uses a defualt of 100, though 50 are usually given

########################
#       RELAX          #
########################
# It is highly recommended that structures are relaxed prior to any analysis with rosetta
# If you are starting from structures that have not been relaxed in rosetta, you MUST initiate this step
# you should generate a number of stuctures and pick the best to carry forward with the following options
RELAX=True                  # True False            # Relax structure : relax into rosetta forcefield, creates output folder 'rosetta-1-relax'
RIPDB=dockpdb.pdb               # Filename False        # Input pdb file. (if False you can provide a list of pre-relaxed structures for mutagenesis with MIPDB)
RINST=4                    # Integer               # How many relaxations? : The number relaxations to make
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
MONST=10                 # Integer               # many mutation runs? : The number of fast-relax-design runs to make PER relaxed structure (CPU must be a factor of this num
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
F 92 93 94 95
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
IIFAC="F"                   # chain_names           # qoute chains of single group, if analaysing the interface between chain groups A and B vs C and D then use "A B"
IOSQE=True                  # True False            # For each of the analysed structures output, name - sequence - total energy - interface energy
IOMET=True                 # True False            # Produce output metrics and graphs for the Interface Energy.			!!!!! INCOMPLETE !!!!!
IOPRS=True                  # True False            # Energy breakdown per residue across interface output. - THIS PRODUCES A HUGE FILE!

########################
#       CLUSTER        #
########################
# Sequence clustering and energy analysis.
action=both               # either 'rmsd' or 'sequence' or 'both' or 'False'
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
DEBUG=True                  # True False            # If True, then intermediate files are not cleaned up!

# find the current working directory
oridir=$(pwd)

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

# checking for contridictory or badly formatted inputs

if [ $RELAX = True ] ; then
	if [ -f $RIPDB ] ; then 
		echo "read initial pdb file as $RIPDB for relaxations"
		# how many chains?
		cnum=$(grep ' CA ' $RIPDB | cut -c 22 | sort | uniq  | sed '/^[[:space:]]*$/d' | wc -l )
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
						echo TER >> fixt-$RIPDB ; mv fixt-$RIPDB fix-$RIPDB 
					done
					
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

# remove underscores from file names
RIPDB=$(echo $RIPDB | sed 's/_/-/g')  2>/dev/null 
sed -i 's/_/-/g' $MIPDB  2>/dev/null 

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
#    Sequence Matrix   #   
########################

seqmatrix ()
{
# here we identify sequence similarities by a naive scoring function.
file=$1
echo "setting up . . ."
# check integrity of input
countnames=$(grep '>' $file | wc -l )
countsequences=$(grep -v '>' $file | wc -l )
if [ $countnames != $countsequences ] ; then
        echo "more than one chain per item, this code is for single chains"
        exit 0
fi

# check sequence lengths
clengthsequences=$(grep -v '>' $file | awk '{print length($0)}' | sort | uniq -c | wc -l )
if [ $clengthsequences != 1 ] ; then
        echo "sequences are of different lengths"
        exit 0
fi
lengthsequences=$( grep -v '>' $file | awk '{print length($0)}' | sort | uniq -c | awk '{print $2}')

# generate working files
echo name > names.tx
grep '>' $file >> names.tx
cp names.tx names2.tx
grep -v '>' $file | sed 's/\(.\{1\}\)/\1,/g' > sequence.tx

### find non-unique sequence positions
ffirst=1
# how many amino acids? 
AAs=$(head -1 sequence.tx | grep -o ',' | wc -l)
for i in $(seq 1 $AAs) ; do
	echo -en "\rreducing $i"
	if [ $(cut -d, -f $i sequence.tx | sort | uniq | wc -l) != 1 ] ; then 
		if [ $ffirst = 1 ] ; then 
			cut -d, -f $i sequence.tx > sequence3.tx
			ffirst=2
		else 
			cut -d, -f $i sequence.tx > sequence2.tx
			paste -d, sequence3.tx sequence2.tx > sequence4.tx
			mv sequence4.tx sequence3.tx
		fi
	fi 
done
mv sequence3.tx sequence.tx

for count in $(seq 1 $(grep -v '>' sequence.tx | wc -l) ) ; do
        echo -en "\rgenerating matrix $count"
        initial=$(sed -n ''$count'p' sequence.tx)
        ncount=$(( $count + 1 ))
        sed -n ''$ncount'p' names.tx > current.tx
        rm initial.tx 2>/dev/null
        for i in $(seq 1 $countsequences) ; do
                echo $initial >> initial.tx
        done

        awk '
BEGIN{
  FS=OFS=","
}
FNR==NR{
  for(i=1;i<=NF;i++){
    array[FNR,i]=$i
  }
  next
}
{
  for(i=1;i<=NF;i++){
    $i=array[FNR,i] $i
  }
}
1
' initial.tx  sequence.tx > acumu.tx

        # MERGE THE FILES
        sed -i 's/AA/0/g;s/CC/0/g;s/DD/0/g;s/EE/0/g;s/FF/0/g;
                                s/GG/0/g;s/HH/0/g;s/II/0/g;s/KK/0/g;s/LL/0/g;
                                s/MM/0/g;s/NN/0/g;s/PP/0/g;s/QQ/0/g;s/RR/0/g;
                                s/SS/0/g;s/TT/0/g;s/VV/0/g;s/WW/0/g;s/YY/0/g' acumu.tx
        sed -i 's/AP/1/g;s/AS/1/g;s/AT/1/g;s/AG/1/g' acumu.tx
        sed -i 's/CV/1/g' acumu.tx
        sed -i 's/DN/1/g;s/DE/1/g' acumu.tx
        sed -i 's/EN/1/g;s/EQ/1/g;s/EH/1/g;s/ED/1/g' acumu.tx
        sed -i 's/FW/1/g;s/FY/1/g;s/FL/1/g;s/FI/1/g;s/FM/1/g' acumu.tx
        sed -i 's/GA/1/g;s/GP/1/g;s/GG/1/g' acumu.tx
        sed -i 's/HK/1/g;s/HR/1/g;s/HE/1/g;s/HQ/1/g' acumu.tx
        sed -i 's/IF/1/g;s/IY/1/g;s/IL/1/g;s/IM/1/g;s/IV/1/g' acumu.tx
        sed -i 's/KR/1/g;s/KQ/1/g;s/KH/1/g' acumu.tx
        sed -i 's/LF/1/g;s/LY/1/g;s/LV/1/g;s/LI/1/g;s/LM/1/g' acumu.tx
        sed -i 's/MF/1/g;s/MY/1/g;s/MV/1/g;s/MI/1/g;s/ML/1/g' acumu.tx
        sed -i 's/NE/1/g;s/ND/1/g;s/NQ/1/g' acumu.tx
        sed -i 's/PA/1/g;s/PS/1/g;s/PT/1/g;s/PG/1/g' acumu.tx
        sed -i 's/QK/1/g;s/QE/1/g;s/QH/1/g;s/QN/1/g' acumu.tx
        sed -i 's/RK/1/g;s/RH/1/g' acumu.tx
        sed -i 's/SA/1/g;s/SP/1/g;s/ST/1/g;s/SG/1/g' acumu.tx
        sed -i 's/TA/1/g;s/TP/1/g;s/TS/1/g' acumu.tx
        sed -i 's/VL/1/g;s/VI/1/g;s/VM/1/g;s/VC/1/g' acumu.tx
        sed -i 's/WF/1/g;s/WY/1/g' acumu.tx
        sed -i 's/YW/1/g;s/YF/1/g;s/YL/1/g;s/YI/1/g;s/YM/1/g' acumu.tx
        sed -i 's/[A-Z][A-Z]/2/g' acumu.tx

        awk -F',' '{ for(i=1; i<=NF;i++) j+=$i; print j; j=0 }' acumu.tx > current2.tx
        cat current.tx current2.tx > current3.tx
        paste -d',' names2.tx current3.tx > temp.tx
        mv temp.tx names2.tx
done
mv names2.tx matrix.csv
echo ''
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
#   Cluster function   #  
########################.

printlink ()
{
echo "
# In[1]: Import
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pylab import savefig
from scipy.spatial import distance
from scipy.cluster import hierarchy

# In[2]: load sequence similarity matrix as df, set names
df = pd.read_csv('$met-matrix.csv')
df = df.set_index('name')
names=np.array(df.columns.values)

# In[3]: Plot original data as heatmap
plt.figure(figsize = (12,10))
sns_heat = sns.heatmap(df)
figure = sns_heat.get_figure()
figure.savefig('heat-native.png', dpi=800)

# In[7]:  Generate hierarchy linkage based on one of the following methods
# 'single', 'complete', 'average', 'ward'
row_linkage = hierarchy.linkage(distance.pdist(df), method='average')

# Generate Cluster map of hierarchy
# sns_clus = sns.clustermap(df, row_linkage=row_linkage, col_linkage=row_linkage, figsize=(10, 10))
# plt.savefig('heat-cluster.png', dpi=800)

# In[5]: Based on 'cutree' pull out linkage groupings
cutree = [3, 9, 27, 81, 243, 729, 2187]
my_array = []
for i in cutree :
    clusters = hierarchy.fcluster(row_linkage, i, 'maxclust')
    my_array.append(clusters)
np.savetxt('heat-clusters.csv', my_array, delimiter=',',  fmt='%s')
np.savetxt('heat-names.csv', names, fmt='%s') " > cluster-linkage.py
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
			
			scp $PROXY $oridir/$i $USER@$COMPUTE:$remotedir/$date-RosettaAutoScore/. 1>/dev/null
			echo "-database /apps/applications/rosetta/3.10/1/default/main/database
-nstruct 1
-out:no_nstruct_label
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

	echo '1.0 Relaxation'
	
	# For local compute	
	if [ $COMPUTE = LOCAL ] ; then 
	################################################	
	########################  CPU CONTROL
		echo '1.1 - Running local relaxations'
		
		# duplicate input files for each CPU
		rm ripdb.tx 2>/dev/null 
		for i in $(seq 1 $CPU) ; do 
			echo $oridir/$RIPDB >> ripdb.tx
		done
								# debug!
		cp $oridir/$RIPDB .
		RIPDBL=ripdb.tx

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
-database /usr/local/rosetta/main/database
-nstruct $(( $RINST / $CPU ))
-score:weights ref2015" > rflag$i.file
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
		
		scp $PROXY $oridir/$RIPDB $USER@$COMPUTE:$remotedir/$date-RosettaAutoRelax/. 1>/dev/null
		echo "-database /apps/applications/rosetta/3.10/1/default/main/database
-out:suffix _"'$SGE_TASK_ID'"
-relax:fast
-out:suffix _"'$state'"
-out:file:scorefile score.sc
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
		ssh $PROXY $USER@$COMPUTE rm -r $remotedir/$date-RosettaAutoRelax
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
				pdbcount=$(ls *.pdb | wc -l ) 2>/dev/null 
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
			
			scp $PROXY $oridir/rosetta-2-mutate/$i $USER@$COMPUTE:$remotedir/$date-RosettaAutoMutate/. 1>/dev/null
			scp $PROXY $oridir/rosetta-2-mutate/move.map $USER@$COMPUTE:$remotedir/$date-RosettaAutoMutate/. 1>/dev/null
			scp $PROXY $oridir/rosetta-2-mutate/mires-auto-mutate.tx $USER@$COMPUTE:$remotedir/$date-RosettaAutoMutate/res.file 1>/dev/null
			echo "-database /apps/applications/rosetta/3.10/1/default/main/database
-relax:fast
-out:suffix _"'$state'"
-out:file:scorefile score.sc
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
			ssh $PROXY $USER@$COMPUTE rm -r $remotedir/$date-RosettaAutoMutate
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
		#cp $oridir/rosetta-1-relax/*.pdb .
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
		echo "4.1 - Running remote energy breakdown	analysis"
		# process the input list for the number of CPUs
		for i in $(seq 1 $CPU) ; do
			#make a temp directory per cpu
			mkdir $oridir/rosetta-4-analysis/temp$i
			cd $oridir/rosetta-4-analysis/temp$i
			#make a flag file per cpu
			echo "#specific options for per-res
-database /usr/local/rosetta/main/database
-score:weights ref2015
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
-score:weights ref2015
-out:file:score_only perres.out \" > flag_"'$SGE_TASK_ID'".file

/apps/applications/rosetta/3.10/1/default/main/source/bin/residue_energy_breakdown.linuxgccrelease @flag_"'$SGE_TASK_ID'".file
cat default.out >> perres.out
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
		scp $PROXY $USER@$COMPUTE:$remotedir/$date-RosettaAutoPer/perres.out .
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

echo 'collating sequence and energy data'
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
	echo 'name,total-REU,total-DDG,Sequence' > mutate-output.csv
	if [ $INTER = True ] ; then 
		echo 'name,total-REU,total-DDG,Interface-REU,inter-DDG,Sequence' > mutate-output.csv
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
	sort -t',' -k1,1 mutate-output.csv > temp ; mv temp mutate-output.csv
fi
if [ $DEBUG != True ] ; then
	rm *.tx
fi

########################################################################################################################################################################
#    Clustering    #
########################################################################################################################################################################
echo 5.0 Clustering
if [ $action = rmsd ] ; then
	rmsd=yes
	sequence=no
	matrix="rmsd"
elif [ $action = sequence ] ; then 
	rmsd=no
	sequence=yes
	matrix="sequence"
elif [ $action = both ] ; then 
	rmsd=yes
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
fi 

cp $oridir/rosetta-3-inter/$pdb*.pdb .								 		#	<--- this probably need fixing
for i in $(cat $oridir/mipdb-auto.tx) ; do  
	cp $oridir/rosetta-1-relax/$i .
done
ls -v *.pdb > pdb.list
cp $oridir/rosetta-4-analysis/mutate-output.csv .
cp $oridir/rosetta-4-analysis/rosetta-spec-sequence.fa .

###########################
## GENERATE SEQUENCE FILE
if [ $sequence = yes ] ; then 
	rm mutseq.fa 2>/dev/null
	for i in $( seq 1 $( grep -v ",Sequence$" mutate-output.csv | wc -l ) ) ; do 
		name=$(grep -v ",Sequence$" mutate-output.csv | sed -n ''$i'p' | awk -F',' '{print $1}')
		sequ=$(grep -v ",Sequence$" mutate-output.csv | sed -n ''$i'p' | rev | awk -F',' '{print $1}' | rev )
		echo '>'$name >> mutseq.fa
		echo $sequ >> mutseq.fa
	done
	seqmatrix mutseq.fa
	sed 's/>//g' matrix.csv > sequence-matrix.csv
fi

###########################
## GENERATE RMSD SIMILARITY MATRIX
ls -v *.pdb > pdb.list
if [ $(wc -l < pdb.list) -gt 2000 ] ; then 
	echo "too many pdbs for rms matrix!"
	rmsd=no
fi

if [ $rmsd = yes ] ; then 
	rmsmatrix 1>/dev/null 
	mv matrix.csv rmsd-matrix.csv
fi

###########################
## CALCULATE CLUSTER LINKAGE

for met in $matrix ; do 

#for cluslink in single complete average ward ; do 
for cluslink in average ; do 
printlink
rm heat-clusters.csv heat-names.csv 2>/dev/null 
sed -i "s/average/$cluslink/g" cluster-linkage.py

echo '
Calculating linkage'
python3 cluster-linkage.py


########################
#    clean heat-names  #			<--- does this need fixing/removing?  do we still have #'.pdb_chain' identifiers?
########################
if [ $met = sequence ] ; then 
sed -i "s/.pdb_$chain//g" heat-names.csv
fi
########################
#   No Dendrogram      # 
########################
sed 's/,/ /g' heat-clusters.csv > heat-clusters-2.tx
transpose  heat-clusters-2.tx  | nl > heat-clusters-1.tx
colv=$(awk '{print NF}' heat-clusters-1.tx | sort -nu | tail -n 1)
sort -n -k$colv,$colv heat-clusters-1.tx | nl | sort -n -k2,2 | awk '{print $1}' > neworder1.tx
echo '0' > head.tx
cat head.tx neworder1.tx > neworder.tx
paste neworder.tx $met-matrix.csv | sed 's/,/ /g' | sort -n -k1,1 | cut -f 2- > newmat1.tx
transpose newmat1.tx | sed 's/ /,/g' > newmat1t.tx
paste -d',' neworder.tx newmat1t.tx | sort -n -k1,1 | cut -d',' -f 2- > clustered-matrix.csv

echo "
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pylab import savefig

# In[2]: load sequence similarity matrix as df, set names
df = pd.read_csv('clustered-matrix.csv')
df = df.set_index('name')
names=np.array(df.columns.values)

# In[3]: Plot original data as heatmap
plt.figure(figsize = (12,10))
sns_heat = sns.heatmap(df)
figure = sns_heat.get_figure()
figure.savefig('heat-cluster-matrix.png', dpi=800) "> cluster-matrix.py
echo 'Generate Cluster-map'
python3 cluster-matrix.py

########################
# average energy matrix #
########################
sed -i 's/,/ /g' heat-clusters.csv
transpose heat-clusters.csv | sed 's/ /,/g' > transout.tx

sed 's/-chain.//g;s/>//g' heat-names.csv > heat-names2.tx
#sed '1,2d'  no_pack_input_score.sc | awk '{print $42}' | sed 's/_0001$//g'  > namecheck.tx
if [ $INTER = True ] ; then 
	energyval=4
else
	energyval=2
fi
grep -v ",Sequence$" mutate-output.csv |  awk -F',' '{print $1}' > namecheck.tx
grep -v ",Sequence$" mutate-output.csv |  awk -F',' -v var=$energyval  '{print $var}'  > interenegy.tx

if [ $(sdiff heat-names2.tx namecheck.tx | grep "|\|>\|<" | wc -l ) != 0 ] ; then
        echo 'we have found that the total length or order of the data in the score files and matrix analysis differs'
        sdiff heat-names2.tx namecheck.tx | grep "|\|>\|<"
        sdiff heat-names2.tx namecheck.tx | grep "|\|>\|<" > issue.tx
        if [ $(awk '{print $1}' issue.tx) = '>' ] ; then
                issue=$(nl heat-names.csv | grep $(awk '{print $2}' issue.tx) | awk '{print $1}')
                sed -i ''$issue'd' heat-names.csv
                sed -i ''$issue'd' interenegy.tx
                if [ $(sdiff heat-names2.tx namecheck.tx | grep "|\|>\|<" | wc -l ) != 0 ] ; then
                        echo 'we could not fix the problem!'
                        exit 0
                fi
        fi
fi

if [ $(wc -l < heat-names2.tx) = $(wc -l < interenegy.tx) ] && [ $(wc -l < heat-names2.tx) = $(wc -l < transout.tx) ] ; then
        paste -d',' heat-names2.tx interenegy.tx transout.tx > naengr.tx
else
        echo ' heat-names2.tx interenegy.tx transout.tx are not equal in length!!!!'
        wc -l heat-names2.tx interenegy.tx transout.tx
        exit 0
fi



sort -nt, -k$(head naengr.tx | awk -F ',' '{print NF}' | sort -nu | head -n 1),$(head naengr.tx | awk -F ',' '{print NF}' | sort -nu | head -n 1) naengr.tx > naengrst.tx
#sed -i 's/,$//g' naengrst.tx
#namelen=$(awk -F',' '{print $1}' naengrst.tx | sort | uniq | wc -L)
#cut -c 1-$(($namelen-$ndigits)) naengrst.tx > front_name.tx
#cut -c $(($namelen-$ndigits+1))-$namelen naengrst.tx | awk -F',' '{print $1}' | sed 's/^/000000/g' | rev | cut -c 1-$ndigits | rev > back_name.tx
#awk -F',' '{print $2,$3,$4,$5,$6,$7,$8,$9,$10}' naengrst.tx > rest.tx

#paste  front_name.tx back_name.tx | awk '{print $1$2}' > dname.tx
#paste -d',' dname.tx rest.tx | sed 's/ /,/g;s/,*$//g' > naengrst.tx


####  DEBUG
#mv heat-cluster-matrix.png $cluslink-cluster-matrix.png
#done
#exit 0

# fix wild-type
for wtname in $(grep -v _.*_ naengrst.tx | awk -F',' '{print $1}') ; do 
	sed -i "s/$wtname,/$wtname-wt,/g" naengrst.tx
done
#echo wt-name is $wtname

# take total average
totave=$(awk -F',' '{print $2}' naengrst.tx | awk -F',' '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' | sed 's/ /,/g')

# how many columns
cols=$(head naengrst.tx | awk -F ',' '{print NF}' | sort -nu | head -n 1)

# make side-header
awk -F',' '{print $1}' naengrst.tx  > sider.tx
nl sider.tx | sort -V -k2,2 > nsider.tx  		# ADDED -V HERE FOR LOCAL
rm  assem2.csv assem2.tx  2>/dev/null
cp sider.tx assem2.tx
touch assem2.csv
total=$(wc -l < sider.tx)
count=1
scount=1
# make general average matrix
cp sider.tx assem.tx
date
for name in $(cat sider.tx) ; do

#name=$(head -1 sider.tx)
if [ $scount -le 1001 ] ; then

echo -en "\rstarting $name, $count of $total"
grep "$name," naengrst.tx > currentline.tx
rm assembly.tx assembly3.tx seen.tx 2>/dev/null
touch seen.tx
cp naengrst.tx naengrstemp.tx
for group in $(seq $cols -1 3 ) ; do

groupval=$(awk -F',' -v gpv=$group '{print $gpv}' currentline.tx)
groupaveval=$(awk -F',' -v col=$group -v gpv=$groupval '{if ($col == gpv) print $2;}' naengrstemp.tx | awk -F',' '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' | sed 's/ /,/g')

awk -F',' -v col=$group -v gpv=$groupval '{if ($col == gpv) print $1;}' naengrstemp.tx | sed "s/$/ $groupaveval/g" > assembly3.tx

grep -v -f seen.tx assembly3.tx >> assembly.tx

awk '{print $1}' assembly3.tx >> seen.tx
sort seen.tx | uniq > seen2.tx
mv seen2.tx seen.tx

#debug
#head assembly3.tx  assembly.tx seen.tx
#wc -l assembly3.tx  assembly.tx seen.tx
#echo group $group groupval $groupval groupaveval $groupaveval totave $totave

done

awk -F',' '{print $1}' naengrstemp.tx | sed "s/$/ $totave/g" > assembly3.tx
grep -v -f seen.tx assembly3.tx >> assembly.tx

nl assembly.tx | awk '{print $1,$2,$3}' | sort -V -k2,2 > nassembly.tx
paste nsider.tx  nassembly.tx | sort -n -k1,1 | awk '{print $5}' > assem.tx
paste -d',' assem2.tx assem.tx > assem3.tx
mv assem3.tx assem2.tx
count=$(($count + 1))

#debug
#cat  assembly3.tx  assembly.tx nassembly.tx  assem.tx assem3.tx assem2.tx
#wc -l assembly3.tx  assembly.tx nassembly.tx  assem.tx assem3.tx assem2.tx

fi
if [ $scount = 1000 ] ; then
sed -i 's/^,//g;s/,$//g' assem2.csv
sed -i 's/^,//g;s/,$//g' assem2.tx
paste -d',' assem2.csv assem2.tx > assem3.tx
mv assem3.tx assem2.csv
rm  assem2.tx ; touch assem2.tx
scount=0
fi
scount=$(($scount + 1))

done
sed -i 's/^,//g;s/,$//g' assem2.csv
sed -i 's/^,//g;s/,$//g' assem2.tx
paste -d',' assem2.csv assem2.tx > assem3.tx
mv assem3.tx assem2.csv

echo '
'
date
echo name > name.tx
cat name.tx sider.tx > header1.tx
transpose header1.tx | sed 's/ /,/g' > header2.tx
cat header2.tx assem2.csv > average-matrix.csv

sed -i "s/^,//g;s/,$//g" average-matrix.csv

echo "
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pylab import savefig

df = pd.read_csv('average-matrix.csv')
df = df.set_index('name')
names=np.array(df.columns.values)

plt.figure(figsize = (12,10))
sns_heat = sns.heatmap(df)
figure = sns_heat.get_figure()
figure.savefig('heat-average.png', dpi=800) " > average-matrix.py

python3 average-matrix.py

mv heat-native.png heat-native-$met.png
mv heat-cluster-matrix.png $cluslink-heat-cluster-matrix-$met.png
mv heat-average.png $cluslink-heat-average-$met.png
mv naengrstemp.tx $cluslink-naengrst-$met.csv
mv average-matrix.csv $cluslink-average-matrix-$met.csv

if [ $DEBUG != True ] ; then
	rm *.tx *.pdb
fi

done
done

########################################################################################################################################################################
#    Clean up    #
########################################################################################################################################################################


echo '5.0 Complete!'
#mv score.sc rosetta-0-initial/.
#rm run-cpu* *.tx nohup.out *flag[0-9]* pdb.list cluster.py 2>/dev/null 
