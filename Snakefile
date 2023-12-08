## Imports
import csv

## Config file and directories
configfile: "config/config.json"
workdir: config["WORKFLOW"]["WORKDIR"]
log_dir   = "/orange/kgraim/panmammalian/Panmammalian/genomes/logs"

#Read in genome 
GENOME_ACCESSIONS = []
GENOME_PREFIX = []
ORGANISM_NAME = []
LAMEN_TERM = []
with open("initial_GCF_runs.csv", "r") as acc_file:
    csv_reader = csv.reader(acc_file)
    for row in csv_reader:
        if (len(row) > 0):
            GENOME_ACCESSIONS.append(row[5])
            ORGANISM_NAME.append(row[2].replace(" ", "_"))
            LAMEN_TERM.append(row[3].replace(" ", "_"))

print(GENOME_ACCESSIONS)
print(ORGANISM_NAME)
print(LAMEN_TERM)

os.getcwd()
## All rule

rule all: 
    input:
        expand(["{organism}" + "/" + "{genome}" + "/"+ "config.txt"], zip, organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "genomic.fna.fai"], zip, organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "genomic.dict"], zip, organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.1.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.2.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.3.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.4.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.5.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.6.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.7.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.8.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "DEXSeqGff.gff"],zip, organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS), 
        expand(["/orange/kgraim/panmammalian/Panmammalian/genomes/snpEff/data/{genome}/{organism}.txt"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS)
## DOWNLOADS RULE:
        # expand(["{organism}" + "/" + "{genome}" + "/" + "genomic.fna"], zip, organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        # expand(["{organism}" + "/" + "{genome}" + "/"+ "genomic.gtf"], zip, organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        #expand(["/orange/kgraim/panmammalian/Panmammalian/genomes/snpEff/data/{genome}/{organism}.txt"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS)
        # expand(["snpEff/data/{genome}/sequences.fa"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        # expand(["snpEff/data/{genome}/genes.gtf"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        # expand(["snpEff/data/{genome}/snpEffectPredictor.bin"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS) 
        # expand(["{organism}" + "/" + "{genome}" + "/" + "genomic.fna.gz"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        # expand(["{organism}" + "/" + "{genome}" + "/" + "genomic.gtf.gz"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        # expand(["snpEff" + "/" + "data" + "/" + "{genome}" + "/" + "{organism}.txt"],zip, organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS)
        #expand(["snpEff" + "/" + "data" + "/" + "{genome}" + "/" + "snpEffectPredictor.bin"],zip ,organism=ORGANISM_NAME,genome=GENOME_ACCESSIONS)
        
## Pipeline 
rule download_assemblies:
    envmodules:
        "edirect",
        "ncbi_cli"
    log:
        stdout=log_dir+ "/" + "{organism}" + "_" + "{genome}_DOWNLOAD.stdout", 
        stderr=log_dir + "/" + "{organism}" + "_" + "{genome}_DOWNLOAD.stderr"
    params:
        dir =  "{organism}",
        acc = "{genome}",
        master_dir = "/orange/kgraim/panmammalian/Panmammalian/genomes"


    output:
        out_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna",
        out_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf",
        config = "{organism}" + "/" + "{genome}" + "/" + "config.txt"

    shell:
        """
        set -o xtrace 
        mkdir -p {params.dir}
        datasets download genome accession {params.acc} --include genome,gtf --filename {params.dir}/{params.acc}.zip > {log.stdout} 2> {log.stderr}
        (cd {params.dir} && unzip {params.acc})
        mv {params.dir}/ncbi_dataset/data/{params.acc}/*genomic.fna {params.dir}/ncbi_dataset/data/{params.acc}/genomic.fna 
        mv {params.dir}/ncbi_dataset/data/assembly_data_report.jsonl {params.dir}/ncbi_dataset/data/{params.acc}/
        mv {params.dir}/ncbi_dataset/data/dataset_catalog.json {params.dir}/ncbi_dataset/data/{params.acc}/
        mv {params.dir}/ncbi_dataset/data/{params.acc}/ {params.dir}/
        mv {params.dir}/README.md {params.dir}/{params.acc}
        rm {params.dir}/{params.acc}.zip
        rm -r {params.dir}/ncbi_dataset
        echo "#FREYA PIPELINE CONFIG FOR {params.dir} CREATED ON $(date)" >> {output.config}
        echo "hisat2_version=2.2.1" >> {output.config}
        echo "fastqc_version=0.11.7" >> {output.config}
        echo "dexcount_version=1.42.0" >> {output.config}
        echo "picard_version=2.25.5" >> {output.config}
        echo "gatk_version=4.1.9.0" >> {output.config}
        echo "snpeff_version=5.0" >> {output.config}
        echo "samtools_version=1.15" >> {output.config}
        echo $'\n' >> {output.config}
        echo "#GENOME FILE: (must also have associated fai file)" >> {output.config}
        echo "CFFA={params.master_dir}/{output.out_fna}" >> {output.config}
        """

#can add to get more files: ,rna,cds,protein,genome,seq-report - ASK KILEY 
#possible values:

rule fai_build:
    envmodules:
        "samtools",
        "picard"
    log:
        stdout=log_dir+ "/" + "{organism}" + "_" + "{genome}_FAI.stdout", 
        stderr=log_dir + "/" + "{organism}" + "_" + "{genome}_FAI.stderr"
    input:
        fa = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna"
    output:
        fa = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna.fai",
        dict_in = "{organism}" + "/" + "{genome}" + "/" + "genomic.dict"
    shell:
        """
        set -o xtrace 
        samtools faidx {input.fa} > {log.stdout} 2> {log.stderr}
        picard CreateSequenceDictionary R={input.fa} > {log.stdout} 2> {log.stderr}
        """

rule hisat_build:
    envmodules:
        "hisat2/2.2.1"
    log:
        stdout=log_dir+ "/" + "{organism}" + "_" + "{genome}_HISAT.stdout", 
        stderr=log_dir + "/" + "{organism}" + "_" + "{genome}_HISAT.stderr"
    params:
        hisat_dir = "hisat2",
        target_dir = "{organism}" + "/" + "{genome}" + "/" + "hisat2",
        target_name = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic",
        config_file = "{organism}" + "/" + "{genome}" + "/" + "config.txt",
        master_dir = "/orange/kgraim/panmammalian/Panmammalian/genomes"
    input:
        in_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna" 
    output:
        target = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.1.ht2",
        target2 = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.2.ht2",
        target3 = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.3.ht2",
        target4 = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.4.ht2",
        target5 = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.5.ht2",
        target6 = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.6.ht2",
        target7 = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.7.ht2",
        target8 = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.8.ht2"
    shell:
        """
        mkdir -p {params.target_dir}
        hisat2-build {input.in_fna} {params.target_name} > {log.stdout} 2> {log.stderr}
        echo $'\n' >> {params.config_file}
        echo "#HISAT2 FILES:" >> {params.config_file}
        echo "HSX={params.master_dir}/{params.target_name}" >> {params.config_file}
        """
## add .fai 

rule dexseq:
    envmodules:
        "htseq/2.0.3"
    log:
        stdout=log_dir+ "/" + "{organism}" + "_" + "{genome}_DEXSEQ.stdout", 
        stderr=log_dir + "/" + "{organism}" + "_" + "{genome}_DEXSEQ.stderr"
    params:
        config_file = "{organism}" + "/" + "{genome}" + "/" + "config.txt",
        master_dir = "/orange/kgraim/panmammalian/Panmammalian/genomes"
    input:
        in_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf"
    output:
        target = "{organism}" + "/" + "{genome}" + "/" + "DEXSeqGff.gff",
    shell:
        """
        set -o xtrace 
        python ./scripts/dexseq_prepare_annotation.py {input.in_gtf} {output.target} > {log.stdout} 2> {log.stderr}
        echo $'\n' >> {params.config_file}
        echo "#DEXSEQ ANNOTATION GFF FILE:" >> {params.config_file}
        echo "DC_GFF={params.master_dir}/{output.target}" >> {params.config_file}
        """

# rule gzip:
#     input:
#         #in_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna",
#         in_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf",
#         halt_requirement1 = "{organism}" + "/" + "{genome}" + "/" + "DEXSeqGff.gff",
#         #halt_requirement2 = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.8.ht2"
#     output:
#         #out_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna.gz",
#         out_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf.gz"
#     shell:
#         """
#         gzip {input.in_gtf}
#         """

rule snpEff:
    envmodules:
        "snpeff/5.0"
    log:
        stdout=log_dir+ "/" + "{organism}" + "_" + "{genome}_SNPEFF.stdout", 
        stderr=log_dir + "/" + "{organism}" + "_" + "{genome}_SNPEFF.stderr"
    params:
        gen = "{genome}",
        org = "{organism}",
        config_file = "/orange/kgraim/panmammalian/Panmammalian/genomes/" + "{organism}" + "/" + "{genome}" + "/" + "config.txt",
        master_dir = "/orange/kgraim/panmammalian/Panmammalian/genomes"
        
    input:
        in_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna",
        in_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf"
        
    output:
        #out_dir = directory("snpEff/data/{genome}")
        #fa = "snpEff/data/{genome}/sequences.fa",
        #gtf = "snpEff/data/{genome}/genes.gtf",
        #target = "snpEff/data/{genome}/snpEffectPredictor.bin",
        out = "/orange/kgraim/panmammalian/Panmammalian/genomes/snpEff/data/{genome}/{organism}.txt"
    shell:
        """
        set -o xtrace 
        pwd
        echo "{params.gen}.genome : {params.org}" >> snpEff/snpEff.config
        mkdir -p snpEff/data/{params.gen}
        cd snpEff/data/{params.gen}
        ln -s /orange/kgraim/panmammalian/Panmammalian/genomes/{input.in_fna} sequences.fa
        ln -s /orange/kgraim/panmammalian/Panmammalian/genomes/{input.in_gtf} genes.gtf
        cd /orange/kgraim/panmammalian/Panmammalian/genomes/snpEff
        snpeff build -c snpEff.config -gtf22 -v {params.gen} > {log.stdout} 2> {log.stderr}
        echo $'\n' >> {params.config_file} 
        echo "#SNPEFF CONFIG:" >> {params.config_file}
        echo $'SnpEffConfigParam="/orange/kgraim/panmammalian/Panmammalian/genomes/snpEff/snpEff.config"' >> {params.config_file}
        echo "#SNPEFF GENOME:" >> {params.config_file}
        echo $'SnpEffGenomeParam="{params.gen}"' >> {params.config_file}
        echo "COMPLETED RUN ON $(date)" >> {output.out}
        echo "All log files in /orange/kgraim/panmammalian/Panmammalian/genomes/logs/" >> {output.out}
        """


