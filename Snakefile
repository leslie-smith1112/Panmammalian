## Imports
import csv

## Config file and directories
configfile: "config/config.json"
workdir: config["WORKFLOW"]["WORKDIR"]
log_dir   = "/orange/kgraim/panmammalian/genomes/logs"

#Read in genome 
GENOME_ACCESSIONS = []
ORGANISM_NAME = []
LAMEN_TERM = []
with open("dog_test.csv", "r") as acc_file:
    csv_reader = csv.reader(acc_file)
    for row in csv_reader:
        if (len(row) > 0):
            GENOME_ACCESSIONS.append(row[5])
            ORGANISM_NAME.append(row[2].replace(" ", "_"))
            LAMEN_TERM.append(row[3].replace(" ", "_"))

print(GENOME_ACCESSIONS)
print(ORGANISM_NAME)
hisat_index=range(1,9)

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
        #expand( "{organism}" + "/" + "{genome}" + "/" + "genomic.fna", organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        #expand("{organism}" + "/" + "{genome}" + "/"+ "genomic.gtf", organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS)
        #SECOND RULE 
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.1.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.2.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.3.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.4.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.5.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.6.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.7.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "hisat2" + "/" + "genomic.8.ht2"],zip ,organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS),
        expand(["{organism}" + "/" + "{genome}" + "/" + "DEXSeqGff.gff"],zip, organism=ORGANISM_NAME, genome=GENOME_ACCESSIONS)
        
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
        mv {params.dir}/ncbi_dataset/data/{params.acc}/ {params.dir}/.
        mv {params.dir}/README.md {params.dir}/{params.acc}
        rm {params.dir}/{params.acc}.zip
        rm -r {params.dir}/ncbi_dataset
        """

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
        "htseq/2.0.3",
        "python"
    input:
        in_gtf = "{organism}" + "/" + "{genome}" + "/" + "genomic.gtf"
    output:
        target = "{organism}" + "/" + "{genome}" + "/" + "DEXSeqGff.gff"
    shell:
        """
        python ./scripts/dexseq_prepare_annotation.py {input.in_gtf} {output.target}
        """
  
'''  
rule snpeff_config_creation: 
    
rule snpeff:
    envmodules:
        "snpeff"
    input:
    output:
    shell:
    """
    echo "{genome}.genome : {organism} >> snpEff.config"
    mkdir snpEff/data/{genome}
    cd snpEff/data/{genome}
    ln -s ../../../genomes/{organism}/{genome}/{genomic.fna} sequences.fa
    ln -s ../../../genomes/{organism}/{genome}/{genomic.gtf} genes.gtf 

    """
        echo "sacCer.genome : Yeast" >> snpEff.config
    





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

rule fastp:
    envmodules:
        "fastqc/0.11.7",
        "fastp"
    log:
        log_dir + "/" + "{sample_sra_id}.fastp_trim.log"
    input:
        #fastq   = fastq_dir + "/" + "{sample_sra_id}.fastq.gz"
        fastq_1 = fastq_dir + "/" + "{sample_sra_id}_1.fastq.gz",
        fastq_2 = fastq_dir + "/" + "{sample_sra_id}_2.fastq.gz"
  
    output:
        #trim = fastq_trimmed_dir + "/" + "{sample_sra_id}" + "/" + "{sample_sra_id}_fastp.fastq.gz",
        trim_1 = fastq_trimmed_dir + "/" + "{sample_sra_id}" + "/" + "{sample_sra_id}_1_fastp.fastq.gz",
        trim_2 = fastq_trimmed_dir + "/" + "{sample_sra_id}" + "/" + "{sample_sra_id}_2_fastp.fastq.gz",
        json_out = fastp_dir + "/" + "{sample_sra_id}" + "/" + "{sample_sra_id}_fastp.json",
        html_out = fastp_dir + "/" + "{sample_sra_id}" + "/" + "{sample_sra_id}_fastp.html"       
        #below: fastp -i {input.fastq_1} -I {input.fastq_2} -o {output.trim_1} -O {output.trim_2} --qualified_quality_phred 15 --length_required 20 --report_title {wildcards.sample_sra_id} --json {output.json_out} --html {output.html_out} 
        #singleend: fastp -i {input.fastq} -o {output.trim}  --qualified_quality_phred 15 --length_required 20 --report_title {wildcards.sample_sra_id} --json {output.json_out} --html {output.html_out}

    shell:
        """
        mkdir -p {fastq_trimmed_dir}/{wildcards.sample_sra_id}
        fastp -i {input.fastq_1} -I {input.fastq_2} -o {output.trim_1} -O {output.trim_2} --qualified_quality_phred 15 --length_required 20 --report_title {wildcards.sample_sra_id} --json {output.json_out} --html {output.html_out} 
        """
rule indexprep:
    params:
        genome = genome_dir + "/" + config["GENOMES"]["PRIMARY_ASSEMBLY"],
        transcripts =  genome_dir  + "/" + config["GENOMES"]["TRANSCRIPTS"]
    output:
        target_decoys = "decoys.txt",
        reference_file = "gentrome.fa.gz"
    shell:
        """
        grep "^>" <(gunzip -c {params.genome}) | cut -d " " -f 1 > {output.target_decoys} 
        sed -i.bak -e 's/>//g' {output.target_decoys}
        cat {params.transcripts} {params.genome} > {output.reference_file}

        """

rule salmonindex:
    envmodules:
        "salmon"
    input:
        target_decoys = "decoys.txt",
        reference_file = "HOMO_SAPIENS_TRANSCRIPTOME_SHORT_1601578642.tar.gz"
        #reference_file = "gentrome.fa.gz"
    output:
        reflens = salmon_dir + "/" + "complete_ref_lens.bin",
        dupclus = salmon_dir + "/" + "duplicate_clusters.tsv",
        pos = salmon_dir + "/" + "pos.bin",
        refAccumLengths = salmon_dir + "/" + "refAccumLengths.bin",
        refseq = salmon_dir + "/" + "refseq.bin",
        ctable = salmon_dir + "/" + "ctable.bin",
        info = salmon_dir + "/" + "info.json",
        preindex = salmon_dir + "/" + "pre_indexing.log",
        refindex = salmon_dir + "/" + "ref_indexing.log",
        seq = salmon_dir + "/" + "seq.bin",
        mphf = salmon_dir + "/" + "mphf.bin",
        ctgoffsets = salmon_dir + "/" + "ctg_offsets.bin",
        rank = salmon_dir + "/" + "rank.bin",
        reflengths = salmon_dir + "/" + "reflengths.bin",
        version = salmon_dir + "/" + "versionInfo.json"

    shell:
        """
        salmon index -t {input.reference_file} -d {input.target_decoys} -p 12 -i {salmon_dir} --gencode
        """

rule salmonbuild:
    envmodules:
        "salmon/0.13.1"
    input:
        tsv = salmon_dir + "/" + "duplicate_clusters.tsv",
        txt = salmon_dir + "/" + "genes_to_transcripts.txt",
        hashbin = salmon_dir + "/" + "hash.bin",
        header = salmon_dir + "/" + "header.json",
        index = salmon_dir + "/" + "indexing.log",
        quasi_index = salmon_dir + "/" + "quasi_index.log",
        refInfo = salmon_dir + "/" + "refInfo.json",
        rsdbin = salmon_dir + "/" + "rsd.bin",
        sabin = salmon_dir + "/" + "sa.bin",
        txp = salmon_dir + "/" + "txpInfo.bin",
        version = salmon_dir + "/" + "versionInfo.json",
        #trim = fastq_trimmed_dir + "/" + "{sample_sra_id}" + "/" + "{sample_sra_id}_fastp.fastq.gz"
        #trim_1 =  fastq_trimmed_dir + "/" + "{sample_sra_id}" + "/" + "{sample_sra_id}_1_fastp.fastq.gz",
        #trim_2 = fastq_trimmed_dir + "/" + "{sample_sra_id}" + "/" + "{sample_sra_id}_2_fastp.fastq.gz"
        #trim = fastq_dir + "/" + "{sample_sra_id}.fastq.gz",
        trim_1 = fastq_dir + "/" + "{sample_sra_id}_1.fastq.gz",
        trim_2 = fastq_dir + "/" + "{sample_sra_id}_2.fastq.gz"
    params:
        result_dir = final_results_dir + "/" + "{sample_sra_id}"

    output:
        quant = final_results_dir + "/" + "{sample_sra_id}" +  "/" + "quant.sf",
        cmd = final_results_dir + "/" + "{sample_sra_id}" +  "/" + "cmd_info.json",
        libformat= final_results_dir + "/" + "{sample_sra_id}" +  "/" + "lib_format_counts.json"     
        #salmon quant -i {salmon_dir} -l A -1 {input.trim_1} -2 {input.trim_2} -o {params.result_dir} --gcBias --seqBias --threads 4 
        #salmon quant -i {salmon_dir} -l A -r {input.trim} -o {params.result_dir}  --seqBias --threads 4
    shell:
        """
        salmon quant -i {salmon_dir} -l A -1 {input.trim_1} -2 {input.trim_2} -o {params.result_dir} --gcBias --seqBias --threads 4 
        """

"""
rule expression:
    envmodules:
        "R"
    params:
        indir = "salmon_quant",
        accession_file = config["INPUT"]["ACCESSIONS_LIST"],
        SRP_study = config["INPUT"]["SRP_STUDY"]
    input:
        quant = final_results_dir + "/" + "{sample_sra_id}" +  "/" + "quant.sf"
        # input dir: final_results
    output:
        outfile = "expression.tsv"
    shell:
        
        Rscript /home/leslie.smith1/blue_kgraim/leslie.smith1/new_SRA_dowloads/get_gene_expression_tximport.R {params.accession_file} {params.indir} {Routdir} {params.SRP_study}
        """
        '''

