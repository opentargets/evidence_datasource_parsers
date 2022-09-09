import os
import tarfile
from datetime import datetime
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

GS = GSRemoteProvider()
HTTP = HTTPRemoteProvider()

# Settings.
timeStamp = datetime.now().strftime('%Y-%m-%d')  # YYYY-MM-DD.
configfile: 'configuration.yaml'

# Extracting EFO version from the configuration and export to the environment:
EFO_version = config['global']['EFOVersion']
os.environ['EFO_VERSION'] = EFO_version

# The master list of all files with their local and remote filenames to avoid code duplication. Only the files specified
# in this list will be generated and uploaded by the "all" and "local" rules.
ALL_FILES = [
    ('cancer_biomarkers.json.gz', GS.remote(f"{config['cancerBiomarkers']['outputBucket']}/cancer_biomarkers-{timeStamp}.json.gz")),
    ('chembl.json.gz', GS.remote(f"{config['ChEMBL']['outputBucket']}/chembl-{timeStamp}.json.gz")),
    ('clingen.json.gz', GS.remote(f"{config['ClinGen']['outputBucket']}/clingen-{timeStamp}.json.gz")),
    ('clingen-Gene-Disease-Summary.csv', GS.remote(f"{config['ClinGen']['inputBucket']}/clingen-Gene-Disease-Summary-{timeStamp}.csv")),
    ('crispr.json.gz', GS.remote(f"{config['CRISPR']['outputBucket']}/crispr-{timeStamp}.json.gz")),
    ('epmc.json.gz', GS.remote(f"{config['EPMC']['outputBucket']}/epmc-{timeStamp}.json.gz")),
    ('gene_burden.json.gz', GS.remote(f"{config['GeneBurden']['outputBucket']}/gene_burden-{timeStamp}.json.gz")),
    ('gene2phenotype.json.gz', GS.remote(f"{config['Gene2Phenotype']['outputBucket']}/gene2phenotype-{timeStamp}.json.gz")),
    ('DDG2P.csv.gz', GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/DDG2P-{timeStamp}.csv.gz")),
    ('EyeG2P.csv.gz', GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/EyeG2P-{timeStamp}.csv.gz")),
    ('SkinG2P.csv.gz', GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/SkinG2P-{timeStamp}.csv.gz")),
    ('CancerG2P.csv.gz', GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/CancerG2P-{timeStamp}.csv.gz")),
    ('intogen.json.gz', GS.remote(f"{config['intOGen']['outputBucket']}/intogen-{timeStamp}.json.gz")),
    ('orphanet.json.gz', GS.remote(f"{config['Orphanet']['outputBucket']}/orphanet-{timeStamp}.json.gz")),
    ('genomics_england.json.gz', GS.remote(f"{config['PanelApp']['outputBucket']}/genomics_england-{timeStamp}.json.gz")),
    # ('impc.json.gz', GS.remote(f"{config['IMPC']['evidenceOutputBucket']}/impc-{timeStamp}.json.gz")),
    # ('mouse_phenotypes.json.gz', GS.remote(f"{config['IMPC']['phenotypesOutputBucket']}/mouse_phenotypes-{timeStamp}.json.gz")),
    ('progeny.json.gz', GS.remote(f"{config['PROGENy']['outputBucket']}/progeny-{timeStamp}.json.gz")),
    ('slapenrich.json.gz', GS.remote(f"{config['SLAPEnrich']['outputBucket']}/slapenrich-{timeStamp}.json.gz")),
    ('sysbio.json.gz', GS.remote(f"{config['SysBio']['outputBucket']}/sysbio-{timeStamp}.json.gz")),
    ('tep.json.gz', GS.remote(f"{config['TEP']['outputBucket']}/tep-{timeStamp}.json.gz")),
    ('safetyLiabilities.json.gz', GS.remote(f"{config['TargetSafety']['outputBucket']}/safetyLiabilities-{timeStamp}.json.gz")),
]
LOCAL_FILENAMES = [f[0] for f in ALL_FILES]
REMOTE_FILENAMES = [f[1] for f in ALL_FILES]

# Auxiliary rules.
rule help:                    # Print help comments for Snakefile.
    input: "Snakefile"
    shell:
        "sed -n 's/^rule//p' {input}"

rule all:                     # Generate all files and upload them to Google Cloud Storage with the current datestamps.
    input:
        LOCAL_FILENAMES
    output:
        REMOTE_FILENAMES + [f"{config['global']['logDir']}/evidence_parser-{timeStamp}.tar.gz"]
    run:
        # The way this works is that remote filenames are actually represented by Snakemake as pseudo-local files.
        # Snakemake will catch the "os.rename" (the same way it would have caught a "mv" call) and proceed with
        # uploading the files.
        for local_filename, remote_filename in zip(input, output[:-1]):
            os.rename(local_filename, remote_filename)
        # Compress, timestamp, and upload the logs.
        with tarfile.open(output[-1], 'w:gz') as archive:
            for filename in os.listdir('log'):
                archive.add(os.path.join('log', filename))

rule local:                   # Generate all files, but do not upload them.
    input:
        LOCAL_FILENAMES

# Data source parsers.
rule cancerBiomarkers:        # Process the Cancers Biomarkers database from Cancer Genome Interpreter.
    input:
        biomarkers_table = GS.remote(config['cancerBiomarkers']['inputAssociationsTable']),
        source_table = GS.remote(config['cancerBiomarkers']['inputSourceTable']),
        disease_table = GS.remote(config['cancerBiomarkers']['inputBucket']),
        drug_index = GS.remote(config['cancerBiomarkers']['drugIndex'])
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        'cancer_biomarkers.json.gz'
    log:
        'log/cancer_biomarkers.log'
    shell:
        """
        # In this and the following rules, the exec call redirects the output of all subsequent commands (both STDOUT
        # and STDERR) to the specified log file.
        exec &> {log}
        python modules/cancerBiomarkers.py \
          --biomarkers_table {input.biomarkers_table} \
          --source_table {input.source_table} \
          --disease_table {input.disease_table} \
          --drug_index {input.drug_index} \
          --output_file {output}
        opentargets_validator --schema {params.schema} {output}
        """

rule chembl:                  # Add the category of why a clinical trial has stopped early to the ChEMBL evidence.
    input:
        evidenceFile = GS.remote(config['ChEMBL']['evidence']),
        stopReasonCategories = GS.remote(config['ChEMBL']['stopReasonCategories'])
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = 'chembl.json.gz'
    log:
        'log/chembl.log'
    shell:
        """
        exec &> {log}
        python modules/ChEMBL.py  \
            --chembl_evidence {input.evidenceFile} \
            --predictions {input.stopReasonCategories} \
            --output {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output}
        """

rule clingen:                 # Process the Gene Validity Curations table from ClinGen.
    params:
        summaryTableWeb = config['ClinGen']['webSource'],
        cacheDir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = 'clingen.json.gz',
        summaryTable = 'clingen-Gene-Disease-Summary.csv'
    log:
        'log/clingen.log'
    shell:
        """
        exec &> {log}
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

rule crispr:                  # Process cancer therapeutic targets using CRISPR–Cas9 screens.
    input:
        evidenceFile = GS.remote(config['CRISPR']['inputAssociationsTable']),
        descriptionsFile = GS.remote(config['CRISPR']['inputDescriptionsTable']),
        cellTypesFile = GS.remote(config['CRISPR']['inputCellTypesTable']")
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = 'crispr.json.gz'
    log:
        'log/crispr.log'
    shell:
        """
        exec &> {log}
        python modules/CRISPR.py \
          --evidence_file {input.evidenceFile} \
          --descriptions_file {input.descriptionsFile} \
          --cell_types_file {input.cellTypesFile} \
          --output_file {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule epmc:                    # Process target/disease evidence strings from ePMC cooccurrence files.
    input:
        inputCooccurences = GS.remote(config['EPMC']['inputBucket'])
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = 'epmc.json.gz'
    log:
        'log/epmc.log'
    shell:
        """
        exec &> {log}
        python modules/EPMC.py \
          --cooccurrences {input.inputCooccurences} \
          --output {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule geneBurden:              # Processes gene burden data from various burden analyses
    input:
        azPhewasBinary = GS.remote(config['GeneBurden']['azPhewasBinary']),
        azPhewasQuant = GS.remote(config['GeneBurden']['azPhewasQuantitative']),
        curation = HTTP.remote(config['GeneBurden']['curation']),
        genebass = GS.remote(config['GeneBurden']['genebass']),
    output:
        evidenceFile = "gene_burden.json.gz"
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    log:
        'log/gene_burden.log'
    shell:
        """
        exec &> {log}
        python modules/GeneBurden.py \
            --az_binary_data {input.azPhewasBinary} \
            --az_quant_data {input.azPhewasQuant} \
            --curated_data {input.curation} \
            --genebass_data {input.genebass} \
            --output {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule gene2Phenotype:          # Processes four gene panels from Gene2Phenotype
    input:
        ddPanel = HTTP.remote(config['Gene2Phenotype']['webSource_dd_panel']),
        eyePanel = HTTP.remote(config['Gene2Phenotype']['webSource_eye_panel']),
        skinPanel = HTTP.remote(config['Gene2Phenotype']['webSource_skin_panel']),
        cancerPanel = HTTP.remote(config['Gene2Phenotype']['webSource_cancer_panel']),
        cardiacPanel = HTTP.remote(config['Gene2Phenotype']['webSource_cardiac_panel'])
    params:
        cacheDir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        ddBucket = 'DDG2P.csv.gz',
        eyeBucket = 'EyeG2P.csv.gz',
        skinBucket = 'SkinG2P.csv.gz',
        cancerBucket = 'CancerG2P.csv.gz',
        cardiacBucket = 'CardiacG2P.csv.gz',
        evidenceFile = 'gene2phenotype.json.gz'
    log:
        'log/gene2Phenotype.log'
    shell:
        """
        exec &> {log}
        # Retain the inputs. They will be later uploaded to GCS.
        cp {input.ddPanel} {output.ddBucket}
        cp {input.eyePanel} {output.eyeBucket}
        cp {input.skinPanel} {output.skinBucket}
        cp {input.cancerPanel} {output.cancerBucket}
        cp {input.cardiacPanel} {output.cardiacBucket}
        python modules/Gene2Phenotype.py \
          --dd_panel {input.ddPanel} \
          --eye_panel {input.eyePanel} \
          --skin_panel {input.skinPanel} \
          --cancer_panel {input.cancerPanel} \
          --cardiac_panel {input.cardiacPanel} \
          --output_file {output.evidenceFile} \
          --cache_dir {params.cacheDir} \
          --local
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule intogen:                 # Process cohorts and driver genes data from intOGen.
    input:
        inputGenes = GS.remote(config['intOGen']['inputAssociationsTable']),
        inputCohorts = GS.remote(config['intOGen']['inputCohortsTable']),
        diseaseMapping = config['intOGen']['diseaseMapping']
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = 'intogen.json.gz'
    log:
        'log/intogen.log'
    shell:
        """
        exec &> {log}
        python modules/IntOGen.py \
          --inputGenes {input.inputGenes} \
          --inputCohorts {input.inputCohorts} \
          --diseaseMapping {input.diseaseMapping} \
          --outputFile {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule orphanet:                # Process disease/target evidence from Orphanet.
    input:
        HTTP.remote(config['Orphanet']['webSource'])
    params:
        cacheDir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        'orphanet.json.gz'
    log:
        'log/orphanet.log'
    shell:
        """
        exec &> {log}
        python modules/Orphanet.py \
          --input_file {input} \
          --output_file {output} \
          --cache_dir {params.cacheDir}
        opentargets_validator --schema {params.schema} {output}
        """

rule panelApp:                # Process gene panels data curated by Genomics England.
    input:
        inputFile = GS.remote(config['PanelApp']['inputAssociationsTable'])
    params:
        cacheDir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = 'genomics_england.json.gz'
    log:
        'log/genomics_england.log'
    shell:
        """
        exec &> {log}
        python modules/PanelApp.py \
          --input-file {input.inputFile} \
          --output-file {output.evidenceFile} \
          --cache_dir {params.cacheDir}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule impc:                    # Process target-disease evidence and mouseModels dataset by querying the IMPC SOLR API.
    params:
        cacheDir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile='impc.json.gz',
        mousePhenotypes='mouse_phenotypes.json.gz'
    log:
        'log/impc.log'
    shell:
        """
        exec &> {log}
        python modules/IMPC.py \
          --cache-dir {params.cacheDir} \
          --output-evidence {output.evidenceFile} \
          --output-mouse-phenotypes {output.mousePhenotypes} \
          --score-cutoff 41
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule progeny:                 # Process gene expression data from TCGA derived from PROGENy.
    input:
        inputFile = GS.remote(config['PROGENy']['inputAssociationsTable']),
        diseaseMapping = config['PROGENy']['inputDiseaseMapping'],
        pathwayMapping = config['PROGENy']['inputPathwayMapping']
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = 'progeny.json.gz'
    log:
        'log/progeny.log'
    shell:
        """
        exec &> {log}
        python modules/PROGENY.py \
          --inputFile {input.inputFile} \
          --diseaseMapping {input.diseaseMapping} \
          --pathwayMapping {input.pathwayMapping} \
          --outputFile {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule slapenrich:              # Process cancer-target evidence strings derived from SLAPenrich.
    input:
        inputFile = GS.remote(config['SLAPEnrich']['inputAssociationsTable']),
        diseaseMapping = config['SLAPEnrich']['inputDiseaseMapping']
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = 'slapenrich.json.gz'
    log:
        'log/slapenrich.log'
    shell:
        """
        exec &> {log}
        python modules/SLAPEnrich.py \
          --inputFile {input.inputFile} \
          --diseaseMapping {input.diseaseMapping} \
          --outputFile {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule sysbio:                  # Process key driver genes for specific diseases that have been curated from Systems Biology papers.
    input:
        evidenceFile = GS.remote(f"{config['SysBio']['inputBucket']}/sysbio_evidence-31-01-2019.tsv"),
        studyFile = GS.remote(f"{config['SysBio']['inputBucket']}/sysbio_publication_info_nov2018.tsv")
    params:
        schema = f"{config['global']['schema']}/opentargets.json"
    output:
        evidenceFile = 'sysbio.json.gz'
    log:
        'log/sysbio.log'
    shell:
        """
        exec &> {log}
        python modules/SystemsBiology.py \
          --evidenceFile {input.evidenceFile} \
          --studyFile {input.studyFile} \
          --outputFile {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

# --- Target annotation data sources parsers --- #
rule targetEnablingPackages:  # Fetching Target Enabling Packages (TEP) data from Structural Genomics Consortium
    params:
        schema = f"{config['global']['schema']}/opentargets_tep.json"
    output:
        'tep.json.gz'
    log:
        'log/tep.log'
    shell:
        """
        exec &> {log}
        python modules/TEP.py  \
          --output_file {output}
        opentargets_validator --schema {params.schema} {output}
        """

rule targetSafety:            # Process data from different sources that describe target safety liabilities.
    input:
        toxcast = GS.remote(config['TargetSafety']['toxcast'])
    params:
        ae = config['TargetSafety']['adverseEvents'],
        sr = config['TargetSafety']['safetyRisk'],
        schema = f"{config['global']['schema']}/opentargets_target_safety.json"
    output:
        'safetyLiabilities.json.gz'
    log:
        'log/safetyLiabilities.log'
    shell:
        """
        exec &> {log}
        python modules/TargetSafety.py  \
            --adverse_events {params.ae} \
            --safety_risk {params.sr} \
            --toxcast {input.toxcast} \
            --output {output}
        opentargets_validator --schema {params.schema} {output}
        """
