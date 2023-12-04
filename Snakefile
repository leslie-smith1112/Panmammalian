## Imports
import csv

## Config file and directories
configfile: "config/config.json"
workdir: config["WORKFLOW"]["WORKDIR"]
log_dir   = "/orange/kgraim/panmammalian/genomes/logs"

#Read in genome 
GENOME_ACCESSIONS = []
GENOME_PREFIX = []
ORGANISM_NAME = []
LAMEN_TERM = []
with open("het.csv", "r") as acc_file:
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
        # FASTQs
        #expand(fastq_dir + "/" + "{sample_sra_id}{pair_id}.fastq.gz", sample_sra_id=SRA_ACCESSIONS, pair_id=["", "_1","_2"])
        #expand([final_results_dir + "/" + "{sample_sra_id}" + "/" + "quant.sf", final_results_dir + "/" + "{sample_sra_id}" + "/" + "cmd_info.json", final_results_dir + "/" + "{sample_sra_id}" + "/" + "lib_format_counts.json"], sample_sra_id=SRA_ACCESSIONS, pair_id=[ "_1","_2"])
        #expand([final_results_dir + "/" + "{sample_sra_id}" + "/" + "quant.sf", final_results_dir + "/" + "{sample_sra_id}" + "/" + "cmd_info.json", final_results_dir + "/" + "{sample_sra_id}" + "/" + "lib_format_counts.json"], sample_sra_id=SRA_ACCESSIONS, pair_id=[ "_1","_2"])
        #expand("{organism}" + "/" + "{organism}"+ "_"+"{genome}_genomic.gff.gz", organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS)
        ## FIRST RULE BELOW
        # expand( "{organism}" + "/" + "{genome}" + "/" + "genomic.fna", organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        # expand("{organism}" + "/" + "{genome}" + "/"+ "genomic.gtf", organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS)
        #SECOND RULE 
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.1.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.2.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.3.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.4.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.5.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.6.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.7.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.8.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "DEXSeqGff.gff"],zip, organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        # zip_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf.gz",
        # zip_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna.gz"
        
## Pipeline 
#fna file 
rule download_assemblies:
    envmodules:
        "edirect",
        "ncbi_cli"
    #log:
     #   outlog = log_dir + "/" + "{organism}" + "_" + "{genome}.gz.log"
    params:
        dir =  "{organism}",
        acc = "{genome}"
    #input:
        #gen_accessions = "acc_list.txt"

    output:
        out_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna",
        out_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf"
        #fna_file =  "{organism}" + "/" + "{organism}" + "_" + "{genome}_genomic.gtf.gz"
        #gtf_file =  "{organism}" + "/" + "{organism}" + "_" + "{genome}_genomic.gff.gz"
        #fastq   = fastq_dir + "/" + "{sample_sra_id}.fastq.gz"
        #outfile= "{organism}/{genome}.zip"
       #fna_out="{params.dir}/ncbi_dataset/data/{params.acc}/{params.acc}{fna_name}genomic.fna",
       #jsonl_out="{params.dir}/ncbi_dataset/data/assembly_data_report.jsonl",
       #json_out="{params.dir}/ncbi_dataset/data/dataset_catalog.json"
       
    #wildcard_constraints:
     #   fna_name=".+"
    shell:
        """
        set -o xtrace
        mkdir -p {params.dir}
        datasets download genome accession {params.acc} --include genome,gtf --filename {params.dir}/{params.acc}.zip
        (cd {params.dir} && unzip {params.acc})
        mv {params.dir}/ncbi_dataset/data/{params.acc}/*genomic.fna {params.dir}/ncbi_dataset/data/{params.acc}/genomic.fna 
        mv {params.dir}/ncbi_dataset/data/assembly_data_report.jsonl {params.dir}/ncbi_dataset/data/{params.acc}/
        mv {params.dir}/ncbi_dataset/data/dataset_catalog.json {params.dir}/ncbi_dataset/data/{params.acc}/
        mv {params.dir}/ncbi_dataset/data/{params.acc}/ {params.dir}/
        mv {params.dir}/README.md {params.dir}/{params.acc}
        rm {params.dir}/{params.acc}.zip
        rm -r {params.dir}/ncbi_dataset
        """
        # set -o xtrace 
        # mkdir -p {params.dir}/{params.genome}
        
        # esearch -db assembly -query {params.genome} </dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_Genbank | 
        # while read -r url ; do
        #     fname=$(echo $url | grep -o 'GC[A,F]_.*' | sed 's/$/_genomic.fna.gz/') ;
        #     wget "$url/$fname" -O {params.dir}/{params.genome}/genomic.fna.gz; 
        # done ;
        # "


        # """
       

#can add to get more files: ,rna,cds,protein,genome,seq-report - ASK KILEY 
#possible values:

rule hisat_build:
    envmodules:
        "hisat2/2.2.1"
    params:
        hisat_dir = "hisat2",
        target_dir = "{organism}" + "/" + "{genome}" + "/" + "hisat2",
        target_name = "{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic"
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
        hisat2-build {input.in_fna} {params.target_name}
        """

rule dexseq:
    envmodules:
        "htseq/2.0.3"
    input:
        in_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf"
    output:
        target = "{organism}" + "/" + "{genome}" + "/" + "DEXSeqGff.gff",
        # zip_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf.gz",
        # zip_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna.gz"
    shell:
        """
        python ./scripts/dexseq_prepare_annotation.py {input.in_gtf} {output.target}
        """

rule snpEff:
    envmodules:
        "snpEff"
    input:
        in_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna",
        in_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf"
        
    output:
        target = "{organism}" + "/" + "{genome}" + "/" + "DEXSeqGff.gff",
        # zip_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf.gz",
        # zip_fna = "{organism}" + "/" + "{genome}" + "/" + "genomic.fna.gz"
    shell:
        """
        echo "{genome}.genome : {organism} >> snpEff/snpEff.config"
        mkdir snpEff/data/{genome}
        cd snpEff/data/{genome}
        ln -s /orange/kgraim/panmammalian/Panmammalian/genomes/{input.in_fna} sequences.fa
        ln -s /orange/kgraim/panmammalian/Panmammalian/genomes/{input.in_gtf} genes.gtf 
        cd /orange/kgraim/panmammalian/Panmammalian/genomes/snpEff
        snpeff build -c snpEff.config -gtf22 -v {organism}
        """


'''  
unzip {params.acc}.zip
    genome:     genomic sequence
                            * rna:        transcript
                            * protein:    amnio acid sequences
                            * cds:        nucleotide coding sequences
                            * gff3:       general feature file
                            * gtf:        gene transfer format
                            * gbff:       GenBank flat file
                            * seq-report: sequence report file
        '''
'''
mkdir -p {params.dir}
        esearch -db assembly -query {params.acc} </dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | while read -r url ; do  \
        fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ; wget --output-document {output.fna_file} "$url/$fname" ; done ;

rule fastqc:
    envmodules:
        "fastqc/0.11.7"
    log:
        log_dir + "/" + "{sample_sra_id}.fastqc.log"

    input:
        #fastq   = fastq_dir + "/" + "{sample_sra_id}.fastq.gz"
        fastq_1 = fastq_dir + "/" + "{sample_sra_id}_1.fastq.gz",
        fastq_2 = fastq_dir + "/" + "{sample_sra_id}_2.fastq.gz"

    output:
        #fastqc_html0 = fastqc_dir + "/" + "{sample_sra_id}.fastqc.html",
        #fastqc_zip0 = fastqc_dir + "/" + "{sample_sra_id}.fastqc.zip"
        fastqc_html1 = fastqc_dir + "/" + "{sample_sra_id}_1.fastqc.html",
        fastqc_zip1 = fastqc_dir + "/" + "{sample_sra_id}_1.fastqc.zip",
        fastqc_html2 = fastqc_dir + "/" + "{sample_sra_id}_2.fastqc.html",
        fastqc_zip2 = fastqc_dir + "/" + "{sample_sra_id}_2.fastqc.zip",
        #these go below: fastqc {input.fastq_1} --quiet --o {fastqc_dir}/{wildcards.sample_sra_id}
        #fastqc {input.fastq_2} --quiet --o {fastqc_dir}/{wildcards.sample_sra_id}
        #fastqc {input.fastq} --quiet --o {fastqc_dir}/{wildcards.sample_sra_id}

    shell:
        """
        fastqc {input.fastq_1} --quiet --o {fastqc_dir}/{wildcards.sample_sra_id}
        fastqc {input.fastq_2} --quiet --o {fastqc_dir}/{wildcards.sample_sra_id}
        """

"""


