Step1: Create a fasta sequence for MMEJ deletions
For this create a fasta file for each MMEJ deletion event. The proxy of MMEJ deletion event is the sequence obtained between  'Start coordinate of the homolog1' to the 'Start coordinate  of the homolog2'. 

Step2:Remove microsatellites
This is peformed by using the tool Phobos where a deletion is filtered if the microsatllite percentage in it is > 80%

Step3: Remove MMEJ deletion events present in the reference genome
	For each type of MMEJ deletion event, create a fasta sequence that resembles a deletion. 
	Align these fasta sequences to the reference genome using bwa aln
	Create an index of the sam file using bwa samse
	Remove the MMEJ deletions with a complete match in the refence genome using the CIGAR field of a the sam file


