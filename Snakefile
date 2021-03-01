from datetime import datetime
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

# --- Settings --- #
# Current time in YYYY-MM-DD-hh-mm format:
timeStamp = datetime.now().strftime("%Y-%m-%d")

# Configuration is read from the config yaml:
configfile: 'configuration.yaml'

# schema version is now hardcoded, version will be read from the command line later:
schemaFile = config['global']['schema_file']
logFile = f"{config['global']['logDir']}/evidence_parser.{timeStamp}.log"

# --- All rules --- #
rule all:
    input:
        GS.remote(f"{config['ClinGen']['outputBucket']}/ClinGen-{timeStamp}.json.gz"),
        GS.remote(f"{config['PheWAS']['outputBucket']}/phewas_catalog-{timeStamp}.json.gz"),
        GS.remote(f"{config['SLAPEnrich']['outputBucket']}/slapenrich-{timeStamp}.json.gz"),
        GS.remote(f"{config['Gene2Phenotype']['outputBucket']}/gene2phenotype-{timeStamp}.json.gz"),
        GS.remote(f"{config['Gene2Phenotype']['outputBucket']}/gene2phenotype_unmapped_diseases-{timeStamp}.txt"),
        GS.remote(f"{config['CRISPR']['outputBucket']}/crispr-{timeStamp}.json.gz"),
        GS.remote(f"{config['SysBio']['outputBucket']}/sysbio-{timeStamp}.json.gz"),
        GS.remote(f"{config['Phenodigm']['outputBucket']}/phenodigm-{timeStamp}.json.gz")

# --- Help Rules --- #
## help               : prints help comments for Snakefile
rule help:
    input: "Snakefile"
    shell:
        "sed -n 's/^##//p' {input}"

## clean               : prints help comments for Snakefile
rule clean:
    shell:
        "rm -rf tmp"

# --- Data sources parsers --- #
## clingen          : processes the Gene Validity Curations table from ClinGen
rule ClinGen:
    input:
        f"tmp/ClinGen-Gene-Disease-Summary-{timeStamp}.csv"
    output:
        evidenceFile=GS.remote(f"{config['ClinGen']['outputBucket']}/ClinGen-{timeStamp}.json.gz"),
        unmappedFile=GS.remote(f"{config['ClinGen']['outputBucket']}/ClinGen-Gene-unmapped-{timeStamp}.lst")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/ClinGen.py -i {input} -o {output.evidenceFile} -u {output.unmappedFile}
        """

## geneticsPortal           : processes lead variants from the Open Targets Genetics portal
rule geneticsPortal:
    input:
        locus2gene=f"tmp/locus2gene-{timeStamp}/",
        toploci=f"tmp/toploci-{timeStamp}.parquet",
        study=f"tmp/studies-{timeStamp}.parquet",
        variantIndex=f"tmp/variantIndex-{timeStamp}.parquet",
        ecoCodes=f"tmp/ecoCodes-{timeStamp}.tsv"
    params:
        threshold=config['GeneticsPortal']['threshold']
    output:
        evidenceFile=GS.remote(f"{config['GeneticsPortal']['outputBucket']}/genetics_portal_evidences{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        '''
        python modules/GeneticsPortal.py \
        --locus2gene {input.locus2gene} \
        --toploci {input.toploci} \
        --study {input.study} \
        --variantIndex {input.variantIndex} \
        --ecoCodes {input.ecoCodes} \
        --threshold {params}
        --outputFile {output}
        '''

## phewas           : processes phenome-wide association studies data from PheWAS
rule phewas:
    input:
        inputFile=f"tmp/phewas_catalog-{timeStamp}.csv",
        consequencesFile=f"tmp/phewas_w_consequences-{timeStamp}.csv",
        diseaseMapping=f"tmp/phewascat_mappings-{timeStamp}.tsv"
    output:
        evidenceFile=GS.remote(f"{config['PheWAS']['outputBucket']}/phewas_catalog-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        '''
        python modules/PheWAS.py \
            --inputFile {input.inputFile} \
            --consequencesFile {input.consequencesFile} \
            --diseaseMapping {input.diseaseMapping} \
            --outputFile {output.evidenceFile}
        '''

## slapenrich           : processes cancer-target evidence strings derived from SLAPenrich
rule slapenrich:
    input:
        inputFile=f"tmp/slapenrich-{timeStamp}.csv",
        diseaseMapping=f"tmp/cancer2EFO_mappings-{timeStamp}.tsv"
    output:
        evidenceFile=GS.remote(f"{config['SLAPEnrich']['outputBucket']}/slapenrich-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/SLAPEnrich.py \
            --inputFile {input.inputFile} \
            --diseaseMapping {input.diseaseMapping} \
            --outputFile {output.evidenceFile} \
        """

## gene2Phenotype           : processes four gene panels from Gene2Phenotype
rule gene2Phenotype:
    input:
        ddPanel=f"tmp/DDG2P-{timeStamp}.csv.gz",
        eyePanel=f"tmp/EyeG2P-{timeStamp}.csv.gz",
        skinPanel=f"tmp/SkinG2P-{timeStamp}.csv.gz",
        cancerPanel=f"tmp/CancerG2P-{timeStamp}.csv.gz"
    output:
        evidenceFile=GS.remote(f"{config['Gene2Phenotype']['outputBucket']}/gene2phenotype-{timeStamp}.json.gz"),
        unmappedDiseases=GS.remote(f"{config['Gene2Phenotype']['outputBucket']}/gene2phenotype_unmapped_diseases-{timeStamp}.txt")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/Gene2Phenotype.py \
            --dd_panel {input.ddPanel} \
            --eye_panel {input.eyePanel} \
            --skin_panel {input.skinPanel} \
            --cancer_panel {input.cancerPanel} \
            --output_file {output.evidenceFile} \
            --unmapped_diseases_file {output.unmappedDiseases}
        """

## crispr           : processes cancer therapeutic targets using CRISPR–Cas9 screens
rule crispr:
    input:
        evidenceFile=f"tmp/crispr_evidence-{timeStamp}.csv",
        descriptionsFile=f"tmp/crispr_descriptions-{timeStamp}.tsv",
        cellTypesFile=f"tmp/crispr_cell_lines-{timeStamp}.tsv"
    output:
        evidenceFile=GS.remote(f"{config['CRISPR']['outputBucket']}/crispr-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/CRISPR.py \
        --evidence_file {input.evidenceFile} \
        --descriptions_file {input.descriptionsFile} \
        --cell_types_file {input.cellTypesFile} \
        --output_file {output.evidenceFile}
        """

## progeny           : processes gene expression data from TCGA derived from PROGENy
rule progeny:
    input:
        inputFile=f"tmp/progeny_normalVStumor_opentargets-{timeStamp}.txt",
        diseaseMapping=f"tmp/progeny_cancer2EFO_mappings-{timeStamp}.tsv",
        pathwayMapping=f"tmp/pathway2Reactome_mappings-{timeStamp}.tsv"
    output:
        evidenceFile=GS.remote(f"{config['PROGENy']['outputBucket']}/progeny-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/PROGENY.py \
        --inputFile {input.inputFile} \
        --diseaseMapping {input.diseaseMapping} \
        --pathwayMapping {input.pathwayMapping} \
        --outputFile {output.evidenceFile}
        """

## phenodigm           : processes target-disease evidence querying the IMPC SOLR API
rule phenodigm:
    output:
        evidenceFile=GS.remote(f"{config['Phenodigm']['outputBucket']}/phenodigm-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/MouseModels.py \
        --update-cache \
        --outputFile {output.evidenceFile}
        """

## sysbio           : processes key driver genes for specific diseases that have been curated from Systems Biology papers
rule sysbio:
    input:
        evidenceFile=f"tmp/sysbio_evidence-{timeStamp}.tsv",
        studyFile=f"tmp/sysbio_publication_info-{timeStamp}.tsv"
    output:
        evidenceFile=GS.remote(f"{config['SysBio']['outputBucket']}/sysbio-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/SystemsBiology.py \
        --evidenceFile {input.evidenceFile} \
        --studyFile {input.studyFile} \
        --outputFile {output.evidenceFile}
        """

# --- Fetching input data and uploading to GS --- #
## fetchClingen          : fetches the Gene Validity Curations table from ClinGen

rule fetchClingen:
    params:
        webSource=config['ClinGen']['webSource']
    output:
        bucket=GS.remote(f"{config['ClinGen']['inputBucket']}/ClinGen-Gene-Disease-Summary-{timeStamp}.csv"),
        local=f"tmp/ClinGen-Gene-Disease-Summary-{timeStamp}.csv"    
    log:
        GS.remote(logFile)
    shell:
        """
        curl  {params.webSource} > {output.local}
        gsutil cp {output.local} {output.bucket}
        """

## fetchGeneticsPortal           : fetches Open Targets Genetics data from their private bucket (not working)
rule fetchGeneticsPortal:
    input:
        #locus2gene=GS.remote(config['GeneticsPortal']['locus2gene']),
        toploci=GS.remote(config['GeneticsPortal']['toploci']),
        #study=GS.remote(config['GeneticsPortal']['study']),
        #variantIndex=GS.remote(config["GeneticsPortal"]["variantIndex"]),
        #ecoCodes=GS.remote(config['GeneticsPortal']['ecoCodes'])
    output:
        #locus2gene=f"tmp/locus2gene-{timeStamp}/*",
        toploci=f"tmp/toploci-{timeStamp}.parquet",
        #study=f"tmp/studies-{timeStamp}.parquet",
        #variantIndex=f"tmp/variantIndex-{timeStamp}.parquet",
        #ecoCodes=f"tmp/ecoCodes-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        gsutil cp {input.toploci} {output.toploci}
        """


## fetchPhewas           : fetches the PheWAS data and an enriched table from GS and a disease mapping look-up
rule fetchPhewas:
    input:
        inputFile=GS.remote(f"{config['PheWAS']['inputBucket']}/phewas-catalog-19-10-2018.csv"),
        consequencesFile=GS.remote(f"{config['PheWAS']['inputBucket']}/phewas_w_consequences.csv")
    params:
        diseaseMapping=config['PheWAS']['diseaseMapping'],
    output:
        inputFile=f"tmp/phewas_catalog-{timeStamp}.csv",
        consequencesFile=f"tmp/phewas_w_consequences-{timeStamp}.csv",
        diseaseMapping=f"tmp/phewascat_mappings-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        curl  {params.diseaseMapping} > {output.diseaseMapping}
        gsutil cp {input.inputFile} {output.inputFile}
        gsutil cp {input.consequencesFile} {output.consequencesFile}
        """

## fetchSlapenrich           : fetches SLAPenrich table from GS
rule fetchSlapenrich:
    input:
        inputFile=GS.remote(f"{config['SLAPEnrich']['inputBucket']}/slapenrich_opentargets-21-12-2017.tsv")
    params:
        diseaseMapping=config['SLAPEnrich']['diseaseMapping']
    output:
        inputFile=f"tmp/slapenrich-{timeStamp}.csv",
        diseaseMapping=f"tmp/cancer2EFO_mappings-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        cp  {params.diseaseMapping} {output.diseaseMapping}
        gsutil cp {input.inputFile} {output.inputFile}
        """

## fetchGene2Phenotype           : fetches four gene panels downloaded from Gene2Phenotype
rule fetchGene2Phenotype:
    params:
        webSource_dd_panel=config['Gene2Phenotype']['webSource_dd_panel'],
        webSource_eye_panel=config['Gene2Phenotype']['webSource_eye_panel'],
        webSource_skin_panel=config['Gene2Phenotype']['webSource_skin_panel'],
        webSource_cancer_panel=config['Gene2Phenotype']['webSource_cancer_panel']

    output:
        ddBucket=GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/DDG2P-{timeStamp}.csv.gz"),
        ddLocal=f"tmp/DDG2P-{timeStamp}.csv.gz",
        eyeBucket=GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/EyeG2P-{timeStamp}.csv.gz"),
        eyeLocal=f"tmp/EyeG2P-{timeStamp}.csv.gz",
        skinBucket=GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/SkinG2P-{timeStamp}.csv.gz"),
        skinLocal=f"tmp/SkinG2P-{timeStamp}.csv.gz",
        cancerBucket=GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/CancerG2P-{timeStamp}.csv.gz"),
        cancerLocal=f"tmp/CancerG2P-{timeStamp}.csv.gz"
    log:
        GS.remote(logFile)
    shell:
        """
        curl  {params.webSource_dd_panel} > {output.ddLocal}
        gsutil cp {output.ddLocal} {output.ddBucket}
        curl  {params.webSource_eye_panel} > {output.eyeLocal}
        gsutil cp {output.eyeLocal} {output.eyeBucket}
        curl  {params.webSource_skin_panel} > {output.skinLocal}
        gsutil cp {output.skinLocal} {output.skinBucket}
        curl  {params.webSource_cancer_panel} > {output.cancerLocal}
        gsutil cp {output.cancerLocal} {output.cancerBucket}
        """

## fetchCrispr           : fetches three tables from GS
rule fetchCrispr:
    input:
        evidenceFile=GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_evidence.tsv"),
        descriptionsFile=GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_descriptions.tsv"),
        cellTypesFile=GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_cell_lines.tsv")
    output:
        evidenceFile=f"tmp/crispr_evidence-{timeStamp}.csv",
        descriptionsFile=f"tmp/crispr_descriptions-{timeStamp}.tsv",
        cellTypesFile=f"tmp/crispr_cell_lines-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        gsutil cp {input.evidenceFile} {output.evidenceFile}
        gsutil cp {input.descriptionsFile} {output.descriptionsFile}
        gsutil cp {input.cellTypesFile} {output.cellTypesFile}
        """

## fetchProgeny           : fetches PROGENy table from GS, and two disease-mapping and pathway-mapping look-up tables
rule fetchProgeny:
    input:
        inputFile=GS.remote(f"{config['PROGENy']['inputBucket']}/progeny_normalVStumor_opentargets.txt"),
        diseaseMapping=config['PROGENy']['diseaseMapping'],
        pathwayMapping=config['PROGENy']['pathwayMapping']
    output:
        inputFile=f"tmp/progeny_normalVStumor_opentargets-{timeStamp}.txt",
        diseaseMapping=f"tmp/progeny_cancer2EFO_mappings-{timeStamp}.tsv", # solve ambiguity between progeny and slapenrich
        pathwayMapping=f"tmp/pathway2Reactome_mappings-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        cp  {input.diseaseMapping} {output.diseaseMapping}
        cp  {input.pathwayMapping} {output.pathwayMapping}
        gsutil cp {input.inputFile} {output.inputFile}
        """

## fetchSysbio           : fetches evidence data and study level information from GS
rule fetchSysbio:
    input:
        evidenceFile=GS.remote(f"{config['SysBio']['inputBucket']}/sysbio_evidence-31-01-2019.tsv"),
        studyFile=GS.remote(f"{config['SysBio']['inputBucket']}/sysbio_publication_info_nov2018.tsv")
    output:
        evidenceFile=f"tmp/sysbio_evidence-{timeStamp}.tsv",
        studyFile=f"tmp/sysbio_publication_info-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        gsutil cp {input.evidenceFile} {output.evidenceFile}
        gsutil cp {input.studyFile} {output.studyFile}
        """

'''
rule PanelApp:
    input:
        inputFile = GS.remote(f"{config['PanelApp']['inputBucket']}/All_genes_20200928-1959.tsv")
    output:
        evidenceFile=GS.remote(f"{config['PanelApp']['outputBucket']}/genomics_england-{timeStamp}.json.gz"),
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/GenomicsEnglandPanelApp.py -i {input.inputFile} -o {output.evidenceFile}
        """
##
## Running intOGen parser
##

rule fetchIntogen:
    input:
        evidenceFile=GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_evidence.tsv"),
        descriptionsFile=GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_descriptions.tsv"),
        cellTypesFile=GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_cell_lines.tsv")
    output:
        evidenceFile=f"tmp/crispr_evidence-{timeStamp}.csv",
        descriptionsFile=f"tmp/crispr_descriptions-{timeStamp}.tsv",
        cellTypesFile=f"tmp/crispr_cell_lines-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        gsutil cp {input.evidenceFile} {output.evidenceFile}
        gsutil cp {input.descriptionsFile} {output.descriptionsFile}
        gsutil cp {input.cellTypesFile} {output.cellTypesFile}
        """

rule intogen:
    input:
        inputGenes = GS.remote(f"{config['intOGen']['inputBucket']}/Compendium_Cancer_Genes.tsv")
        inputCohorts = GS.remote(f"{config['intOGen']['inputBucket']}/cohorts.tsv")
        diseaseMapping = "resources/cancer2EFO_mappings.tsv"
    output:
        evidenceFile=GS.remote(f"{config['intOGen']['outputBucket']}/intogen-{timeStamp}.json.gz"),
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/IntOGen.py -g {input.inputGenes} -c {input.inputCohorts} -d {input.diseaseMapping} -o {output.evidenceFile}
        """

'''