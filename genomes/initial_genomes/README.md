# Steps to running the pipeline to prep genomes for the Paipu pipeline
The following steps are described assuming the workdir = `/orange/kgraim/panmammalian/Panmammalian/genomes/initial_genomes`

1) Make input file with list of mammals you would like to prep genomes for called input.txt (script expects file to be in ). 
	Example:
	
	Ailuropoda melanoleuca
	Antilocapra americana
	Aotus nancymaae
	Bos mutus
	Bos taurus
	Bubalus bubalis
	Callithrix jacchus
	Canis lupus familiaris
	Capra hircus
	Carollia perspicillata

All of the following steps are executed within the master slurm script (`master_prep.sh`).
 NOTE: the pipeline expects input file (`input.txt`) to be created with the format described above and the `initial_genomes` directory. 

2. Input file `input.txt` is taken as input into `ncbi_queries/query.sh` to query reference genome information for each mammal. 
`ncbi_queries/query.sh` is run for each mammal listed in `input.txt` and outputs a json file containing all available genome information for the mammal and associated information such as:
	- date of release
	- whether the genome has gene annotations
	- assembly name
json files for each queried mammal are stored in `ncbi_queries/`. 

After creating the initial mammal's json file,`ncbi_queries/query.sh` then calls `ncbi_queries/json_parse.sh` which parses each mammal's json file to identify the current recommended reference genome for the mammal and whether the genome has a gene annotation file. This information is written to the output file. All queried genomes are in `ncbi_queries/output_all.csv`. Genomes that have a gene annotation can be prepped for Paipu - these genomes are written to `output.csv`

Output all has the following columns (these headers do not exist in the output.csv file at the moment):
`mammal_name, genome_accession, assembly_name, assembly_status, release_date, has_gene_annotation, gene_annotation_proivder, gene_annotation_release_date`

3. `output.csv` is what the `genome_prep.nf` nextflow pipeline expects as input.
`genome_prep.nf` prepares each mammalian genome by:	
	1) Downloading genome and gene annotation files for given mammal/genome.
	2) Creating `.fai` and `.dict` file needed for.
	3) Building the indexes needed for hisat2.
	4) Creating `.gff` files needed for DexSeq.
