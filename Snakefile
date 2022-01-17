from datetime import datetime
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

GS = GSRemoteProvider()
HTTP = HTTPRemoteProvider()

# --- Settings --- #
# Current date in YYYY-MM-DD format:
timeStamp = datetime.now().strftime("%Y-%m-%d")

# Configuration is read from the config yaml:
configfile: 'configuration.yaml'
logFile = f"{config['global']['logDir']}/evidence_parser-{timeStamp}.log"

# --- All rules --- #
rule all:
    input:
        GS.remote(f"{config['cancerBiomarkers']['outputBucket']}/cancer_biomarkers-{timeStamp}.json.gz"),
        GS.remote(f"{config['ClinGen']['outputBucket']}/clingen-{timeStamp}.json.gz"),
        GS.remote(f"{config['CRISPR']['outputBucket']}/crispr-{timeStamp}.json.gz"),
        GS.remote(f"{config['EPMC']['outputBucket']}/epmc-{timeStamp}.json.gz"),
        GS.remote(f"{config['Gene2Phenotype']['outputBucket']}/gene2phenotype-{timeStamp}.json.gz"),
        GS.remote(f"{config['intOGen']['outputBucket']}/intogen-{timeStamp}.json.gz"),
        GS.remote(f"{config['Orphanet']['outputBucket']}/orphanet-{timeStamp}.json.gz"),
        GS.remote(f"{config['PanelApp']['outputBucket']}/genomics_england-{timeStamp}.json.gz"),
        GS.remote(f"{config['Phenodigm']['evidenceOutputBucket']}/phenodigm-{timeStamp}.json.gz"),
        GS.remote(f"{config['Phenodigm']['phenotypesOutputBucket']}/mouse_phenotypes-{timeStamp}.json.gz"),
        GS.remote(f"{config['PheWAS']['outputBucket']}/phewas_catalog-{timeStamp}.json.gz"),
        GS.remote(f"{config['PROGENy']['outputBucket']}/progeny-{timeStamp}.json.gz"),
        GS.remote(f"{config['SLAPEnrich']['outputBucket']}/slapenrich-{timeStamp}.json.gz"),
        GS.remote(f"{config['SysBio']['outputBucket']}/sysbio-{timeStamp}.json.gz"),
        GS.remote(f"{config['TEP']['outputBucket']}/tep-{timeStamp}.json.gz")
        GS.remote(f"{config['TargetSafety']['outputBucket']}/safetyLiabilities-{timeStamp}.json.gz")

# --- Auxiliary Rules --- #
## help                     : prints help comments for Snakefile
rule help:
    input: "Snakefile"
    shell:
        "sed -n 's/^##//p' {input}"

## clean                    : prints help comments for Snakefile
rule clean:
    shell:
        "rm -rf tmp"

# --- Data sources parsers --- #
## cancerBiomarkers          : processes the Cancers Biomarkers database from Cancer Genome Interpreter
rule cancerBiomarkers:
    input:
        biomarkers_table = GS.remote(f"{config['cancerBiomarkers']['inputBucket']}/cancerbiomarkers-2018-05-01.tsv"),
        source_table = GS.remote(f"{config['cancerBiomarkers']['inputBucket']}/cancer_biomarker_source.jsonl"),
        disease_table = GS.remote(f"{config['cancerBiomarkers']['inputBucket']}/cancer_biomarker_disease.jsonl"),
    params:
        # Downloaded separately using wget, because FTP.RemoteProvider cannot handle recursive directory downloads.
        drug_index = config['cancerBiomarkers']['drugIndex'],
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        GS.remote(f"{config['cancerBiomarkers']['outputBucket']}/cancer_biomarkers-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        wget -q -r ftp://{params.drug_index}
        python modules/cancerBiomarkers.py \
          --biomarkers_table {input.biomarkers_table} \
          --source_table {input.source_table} \
          --disease_table {input.disease_table} \
          --drug_index {params.drug_index} \
          --output_file {output}
        opentargets_validator --schema {params.schema} {output}
        """

## clingen                  : processes the Gene Validity Curations table from ClinGen
rule clingen:
    params:
        summaryTableWeb = config['ClinGen']['webSource'],
        cacheDir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        summaryTable = GS.remote(f"{config['ClinGen']['inputBucket']}/clingen-Gene-Disease-Summary-{timeStamp}.csv"),
        evidenceFile = GS.remote(f"{config['ClinGen']['outputBucket']}/clingen-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        # Download directly because HTTP RemoteProvider does not handle retries correctly.
        wget -q -O clingen_summary.csv {params.summaryTableWeb}
        # Retain the original summary table and store that in GCS.
        cp clingen_summary.csv {output.summaryTable}
        python modules/ClinGen.py \
          --input_file clingen_summary.csv \
          --output_file {output.evidenceFile} \
          --cache_dir {params.cacheDir} \
          --local
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

## crispr                   : processes cancer therapeutic targets using CRISPR–Cas9 screens
rule crispr:
    input:
        evidenceFile = GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_evidence.tsv"),
        descriptionsFile = GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_descriptions.tsv"),
        cellTypesFile = GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_cell_lines_enriched_2021-10-22.tsv")
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = GS.remote(f"{config['CRISPR']['outputBucket']}/crispr-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/CRISPR.py \
          --evidence_file {input.evidenceFile} \
          --descriptions_file {input.descriptionsFile} \
          --cell_types_file {input.cellTypesFile} \
          --output_file {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

## epmc                     : processes target/disease evidence strings from ePMC cooccurrence files
rule epmc:
    input:
        inputCooccurences = directory(GS.remote(config['EPMC']['inputBucket']))
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = GS.remote(f"{config['EPMC']['outputBucket']}/epmc-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/EPMC.py \
          --cooccurrenceFile {input.inputCooccurences} \
          --output {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

## gene2Phenotype           : processes four gene panels from Gene2Phenotype
rule gene2Phenotype:
    input:
        ddPanel = HTTP.remote(config['Gene2Phenotype']['webSource_dd_panel']),
        eyePanel = HTTP.remote(config['Gene2Phenotype']['webSource_eye_panel']),
        skinPanel = HTTP.remote(config['Gene2Phenotype']['webSource_skin_panel']),
        cancerPanel = HTTP.remote(config['Gene2Phenotype']['webSource_cancer_panel'])
    params:
        cacheDir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        ddBucket = GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/DDG2P-{timeStamp}.csv.gz"),
        eyeBucket = GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/EyeG2P-{timeStamp}.csv.gz"),
        skinBucket = GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/SkinG2P-{timeStamp}.csv.gz"),
        cancerBucket = GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/CancerG2P-{timeStamp}.csv.gz"),
        evidenceFile= GS.remote(f"{config['Gene2Phenotype']['outputBucket']}/gene2phenotype-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        # Retain the inputs and save to GCS.
        cp {input.ddPanel} {output.ddBucket}
        cp {input.eyePanel} {output.eyeBucket}
        cp {input.skinPanel} {output.skinBucket}
        cp {input.cancerPanel} {output.cancerBucket}
        python modules/Gene2Phenotype.py \
          --dd_panel {input.ddPanel} \
          --eye_panel {input.eyePanel} \
          --skin_panel {input.skinPanel} \
          --cancer_panel {input.cancerPanel} \
          --output_file {output.evidenceFile} \
          --cache_dir {params.cacheDir} \
          --local
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

## intogen                  : processes cohorts and driver genes data from intOGen
rule intogen:
    input:
        inputGenes = GS.remote(f"{config['intOGen']['inputBucket']}/Compendium_Cancer_Genes.tsv"),
        inputCohorts = GS.remote(f"{config['intOGen']['inputBucket']}/cohorts.tsv"),
        diseaseMapping = config['intOGen']['diseaseMapping']
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = GS.remote(f"{config['intOGen']['outputBucket']}/intogen-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/IntOGen.py \
          --inputGenes {input.inputGenes} \
          --inputCohorts {input.inputCohorts} \
          --diseaseMapping {input.diseaseMapping} \
          --outputFile {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

## Orphanet                 : Processing disease/target evidence from Orphanet
rule orphanet:
    input:
        HTTP.remote(config['Orphanet']['webSource'])
    params:
        cacheDir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        GS.remote(f"{config['Orphanet']['outputBucket']}/orphanet-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/Orphanet.py \
          --input_file {input} \
          --output_file {output} \
          --cache_dir {params.cacheDir} \
          --local
        opentargets_validator --schema {params.schema} {output}
        """

## panelApp                 : processes gene panels data curated by Genomics England
rule panelApp:
    input:
        inputFile = GS.remote(f"{config['PanelApp']['inputBucket']}/All_genes_20200928-1959.tsv")
    params:
        cacheDir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = GS.remote(f"{config['PanelApp']['outputBucket']}/genomics_england-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/PanelApp.py \
          --input-file {input.inputFile} \
          --output-file {output.evidenceFile} \
          --cache_dir {params.cacheDir}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

## phenodigm                : processes target-disease evidence and mouseModels dataset by querying the IMPC SOLR API
rule phenodigm:
    params:
        cacheDir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile=GS.remote(f"{config['Phenodigm']['evidenceOutputBucket']}/phenodigm-{timeStamp}.json.gz"),
        mousePhenotypes=GS.remote(f"{config['Phenodigm']['phenotypesOutputBucket']}/mouse_phenotypes-{timeStamp}."
                                  f"json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/PhenoDigm.py \
          --cache-dir {params.cacheDir} \
          --output-evidence {output.evidenceFile} \
          --output-mouse-phenotypes {output.mousePhenotypes}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

## phewas                   : processes phenome-wide association studies data from PheWAS
rule phewas:
    input:
        inputFile = GS.remote(f"{config['PheWAS']['inputBucket']}/phewas-catalog-19-10-2018.csv"),
        consequencesFile = GS.remote(f"{config['PheWAS']['inputBucket']}/phewas_w_consequences.csv"),
        diseaseMapping = HTTP.remote(config['PheWAS']['diseaseMapping'])
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = GS.remote(f"{config['PheWAS']['outputBucket']}/phewas_catalog-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/PheWAS.py \
          --inputFile {input.inputFile} \
          --consequencesFile {input.consequencesFile} \
          --diseaseMapping {input.diseaseMapping} \
          --outputFile {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

## progeny                  : processes gene expression data from TCGA derived from PROGENy
rule progeny:
    input:
        inputFile = GS.remote(f"{config['PROGENy']['inputBucket']}/progeny_normalVStumor_opentargets.txt"),
        diseaseMapping = config['PROGENy']['diseaseMapping'],
        pathwayMapping = config['PROGENy']['pathwayMapping']
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = GS.remote(f"{config['PROGENy']['outputBucket']}/progeny-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/PROGENY.py \
          --inputFile {input.inputFile} \
          --diseaseMapping {input.diseaseMapping} \
          --pathwayMapping {input.pathwayMapping} \
          --outputFile {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

## slapenrich               : processes cancer-target evidence strings derived from SLAPenrich
rule slapenrich:
    input:
        inputFile = GS.remote(f"{config['SLAPEnrich']['inputBucket']}/slapenrich_opentargets-21-12-2017.tsv"),
        diseaseMapping = config['SLAPEnrich']['diseaseMapping']
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = GS.remote(f"{config['SLAPEnrich']['outputBucket']}/slapenrich-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/SLAPEnrich.py \
          --inputFile {input.inputFile} \
          --diseaseMapping {input.diseaseMapping} \
          --outputFile {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

## sysbio                   : processes key driver genes for specific diseases that have been curated from Systems Biology papers
rule sysbio:
    input:
        evidenceFile = GS.remote(f"{config['SysBio']['inputBucket']}/sysbio_evidence-31-01-2019.tsv"),
        studyFile = GS.remote(f"{config['SysBio']['inputBucket']}/sysbio_publication_info_nov2018.tsv")
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = GS.remote(f"{config['SysBio']['outputBucket']}/sysbio-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/SystemsBiology.py \
          --evidenceFile {input.evidenceFile} \
          --studyFile {input.studyFile} \
          --outputFile {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

# --- Target annotation data sources parsers --- #
## TEP                    : Fetching Target Enabling Packages (TEP) data from Structural Genomics Consortium
rule TargetEnablingPackages:
    params:
        schema = f"{config['global']['schema']}/opentargets_tep.json"
    output:
        GS.remote(f"{config['TEP']['outputBucket']}/tep-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/TEP.py  \
          --output_file {output}
        opentargets_validator --schema {params.schema} {output}
        """

## Target Safety        : processes data from different sources that describe target safety liabilities
rule TargetSafety:
    input:
        toxcast = GS.remote(config['TargetSafety']['toxcast'])
    output:
        GS.remote(f"{config['TargetSafety']['outputBucket']}/safetyLiabilities-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/TargetSafety.py  \
            --toxcast {input.toxcast} \
            --output {output}
        """
