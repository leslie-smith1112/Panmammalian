## Imports
import csv
import pandas as pd
#for dir in ${dirs[@]}; do echo $dir; done
#newdir=$(echo $dir | awk '{split($0,a,"__"); print a[1]"__"a[3]}')
#mv $dir $newdir
## Config file and directories
configfile: "config/config.json"
workdir: config["WORKFLOW"]["WORKDIR"]
log_dir   = "/orange/kgraim/panmammalian/Panmammalian/genomes/logs"

#Read in assemb 
## Config file and directories
##configfile: "config/config.json"
#workdir: config["WORKFLOW"]["WORKDIR"]
#Read in assemb 
SPECIES = []
ASSEMBLY = []
GENOME = []
with open("TOGA/TOGA_Mammal_Alignments/tester_single.tsv", "r") as acc_file:
    for line in csv.reader(acc_file, delimiter='\t'):
        SPECIES.append(line[0].replace(" ", "_"))
        #COMMON.append(line[1].replace(" ", "_"))
        ASSEMBLY.append(line[4])
        GENOME.append(line[5])


## All rule

rule all: 
    input:
        ### DOWNLOAD ASSEMBLY AND GTF FILE ### 
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/"+ "genomic.fna"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "genomic.gtf"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
        # expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "genomic.gtf.gz"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
        # expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "config.txt"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME)
        
        ### FAI RULE ### 
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "genomic.fa.fai"],zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "genomic.dict"],zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),

        ### BUILD HISAT INDEXES ### 
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "genomic.1.ht2"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "genomic.2.ht2"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "genomic.3.ht2"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "genomic.4.ht2"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "genomic.5.ht2"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "genomic.6.ht2"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "genomic.7.ht2"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "genomic.8.ht2"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),
        
        ### BUILD DEXSEQ FILE ### 
         expand(["{spec}" + "/" + "{assemb}__{genome}" + "/" + "DEXSeqGff.gff"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME),

        ### SNPEFF ### 
         expand(["/orange/kgraim/panmammalian/Panmammalian/genomes/snpEff/data/{assemb}__{genome}/{spec}.txt"], zip, spec=SPECIES, assemb=ASSEMBLY, genome=GENOME)
## Pipeline 
rule download_assemblies:
    envmodules:
        "edirect",
        "ncbi_cli"
    log:
        stdout=log_dir+ "/" + "{spec}" + "__" + "{assemb}__{genome}__DOWNLOAD.stdout", 
        stderr=log_dir + "/" + "{spec}" + "__" + "{assemb}__{genome}__DOWNLOAD.stderr"
    params:
        species =  "{spec}",
        subspecies = "{assemb}__{genome}",
        genome_acc = "{genome}",
        acc = "{assemb}",
        temp_file = "TEMP.zip",
        master_dir = "/orange/kgraim/panmammalian/Panmammalian/genomes"


    output:
        out_fna = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "codonAlignments.fa",
        out_gtf = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "geneAnnotation.gtf",
        config = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "config.txt"

    shell:
        """
        set -o xtrace 
        if [ ! -d {param.species} ]; then
            mkdir -p {params.species}/{params.subspecies}
        mv TOGA/TOGA_Mammal_Alignments/{params.species}__{params.acc}/* {params.species}/{params.subspecies}/
        gzip -d {params.species}/{params.subspecies}/codonAlignments.fa.gz
        gzip -d {params.species}/{params.subspecies}/geneAnnotation.gtf.gz
        echo "#FREYA PIPELINE CONFIG FOR {params.species} {params.acc} {params.genome_acc} CREATED ON $(date)" >> {output.config}
        echo $'\n' >> {output.config}
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
        # datasets download genome accession {params.genome_acc} --include genome --filename {params.species}/{params.temp_file} > {log.stdout} 2> {log.stderr}
        # (cd {params.species} && unzip {params.temp_file})
        # mv {params.species}/ncbi_dataset/data/{params.genome_acc}/*genomic.fna {params.species}/ncbi_dataset/data/{params.genome_acc}/genomic.fna 
        # mv {params.species}/ncbi_dataset/data/assembly_data_report.jsonl {params.species}/ncbi_dataset/data/{params.genome_acc}/
        # mv {params.species}/ncbi_dataset/data/dataset_catalog.json {params.species}/ncbi_dataset/data/{params.genome_acc}/
        # mv {params.species}/ncbi_dataset/data/{params.genome_acc}/* {params.species}/{params.subspecies}
        # mv {params.species}/README.md {params.species}/{params.subspecies}
        # rm {params.species}/{params.temp_file}
        # rm -r {params.species}/ncbi_dataset
        # cp TOGA/TOGA_Mammal_Alignments/{params.species}__{params.acc}/geneAnnotation.gtf.gz {params.species}/{params.subspecies}/genomic.gtf.gz
#can add to get more files: ,rna,cds,protein,genome,seq-report - ASK KILEY 
#possible values:
rule fai_build:
    envmodules:
        "samtools",
        "picard"
    log:
        stdout=log_dir+ "/" + "{spec}" + "__" + "{assemb}__{genome}__FAI.stdout", 
        stderr=log_dir + "/" + "{spec}" + "__" + "{assemb}__{genome}__FAI.stderr"
    input:
        fa = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "codonAlignments.fa"
    output:
        fa = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "codonAlignments.fa.fai",
        dict_in = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "codonAlignments.dict"
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
        stdout=log_dir+ "/" + "{spec}" + "__" + "{assemb}__{genome}__HISAT.stdout", 
        stderr=log_dir + "/" + "{spec}" + "__" + "{assemb}__{genome}__HISAT.stderr"
    params:
        hisat_dir = "hisat2",
        target_dir = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2",
        target_name = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "codonAlignments",
        config_file = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "config.txt",
        master_dir = "/orange/kgraim/panmammalian/Panmammalian/genomes"
    input:
        in_fna = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "codonAlignments.fa" 
    output:
        target = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "codonAlignments.1.ht2",
        target2 = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "codonAlignments.2.ht2",
        target3 = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "codonAlignments.3.ht2",
        target4 = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "codonAlignments.4.ht2",
        target5 = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "codonAlignments.5.ht2",
        target6 = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "codonAlignments.6.ht2",
        target7 = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "hisat2" + "/" + "codonAlignments.7.ht2",
        target8 = "{spec}" + "/" + "{assemb}__{genome}"+ "/" + "hisat2" + "/" + "codonAlignments.8.ht2"
    shell:
        """
        mkdir -p {params.target_dir}
        hisat2-build {input.in_fna} {params.target_name} > {log.stdout} 2> {log.stderr}
        echo $'\n' >> {params.config_file}
        echo "#HISAT2 FILES:" >> {params.config_file}
        echo "HSX={params.master_dir}/{params.target_name}" >> {params.config_file}
        """


rule dexseq:
    envmodules:
        "htseq/2.0.3"
    log:
        stdout=log_dir+ "/" + "{spec}" + "__" + "{assemb}__{genome}__DEXSEQ.stdout", 
        stderr=log_dir + "/" + "{spec}" + "__" + "{assemb}__{genome}__DEXSEQ.stderr"
    params:
        config_file = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "config.txt",
        master_dir = "/orange/kgraim/panmammalian/Panmammalian/genomes"
    input:
        in_gtf = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "geneAnnotation.gtf"
    output:
        target = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "DEXSeqGff.gff",
    shell:
        """
        set -o xtrace 
        python ./scripts/dexseq_prepare_annotation.py {input.in_gtf} {output.target} > {log.stdout} 2> {log.stderr}
        echo $'\n' >> {params.config_file}
        echo "#DEXSEQ ANNOTATION GFF FILE:" >> {params.config_file}
        echo "DC_GFF={params.master_dir}/{output.target}" >> {params.config_file}
        """

# # rule gzip:
# #     input:
# #         #in_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna",
# #         in_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf",
# #         halt_requirement1 = "{organism}" + "/" + "{genome}" + "/" + "DEXSeqGff.gff",
# #         #halt_requirement2 = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.8.ht2"
# #     output:
# #         #out_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna.gz",
# #         out_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf.gz"
# #     shell:
# #         """
# #         gzip {input.in_gtf}
# #         """

rule snpEff:
    envmodules:
        "snpeff/5.0"
    log:
        stdout=log_dir+ "/" + "{spec}" + "__" + "{assemb}__{genome}__SNPEFF.stdout", 
        stderr=log_dir + "/" +"{spec}" + "__" + "{assemb}__{genome}__SNPEFF.stderr"
    params:
        gen = "{assemb}__{genome}",
        org = "{spec}",
        config_file = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "config.txt",
        master_dir = "/orange/kgraim/panmammalian/Panmammalian/genomes"
        
    input:
        in_fna = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "codonAlignments.fa",
        in_gtf = "{spec}" + "/" + "{assemb}__{genome}" + "/" + "geneAnnotation.gtf"
        
    output:
        #out_dir = directory("snpEff/data/{genome}")
        #fa = "snpEff/data/{genome}/sequences.fa",
        #gtf = "snpEff/data/{genome}/genes.gtf",
        #target = "snpEff/data/{genome}/snpEffectPredictor.bin",
        out = "/orange/kgraim/panmammalian/Panmammalian/genomes/snpEff/data/{assemb}__{genome}/{spec}.txt"
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