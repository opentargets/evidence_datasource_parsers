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
    ('baseline_expression.json.gz', GS.remote(f"{config['baselineExpression']['outputBucket']}/baseline_expression-{timeStamp}.json.gz")),
    ('cancer_biomarkers.json.gz', GS.remote(f"{config['cancerBiomarkers']['outputBucket']}/cancer_biomarkers-{timeStamp}.json.gz")),
    ('chembl.json.gz', GS.remote(f"{config['ChEMBL']['outputBucket']}/chembl-{timeStamp}.json.gz")),
    ('clingen.json.gz', GS.remote(f"{config['ClinGen']['outputBucket']}/clingen-{timeStamp}.json.gz")),
    ('clingen-Gene-Disease-Summary.csv', GS.remote(f"{config['ClinGen']['inputBucket']}/clingen-Gene-Disease-Summary-{timeStamp}.csv")),
    ('project_score.json.gz', GS.remote(f"{config['ProjectScore']['outputBucket']}/project_score-{timeStamp}.json.gz")),
    ('essentiality.json.gz', GS.remote(f"{config['Essentiality']['outputBucket']}/essentiality-{timeStamp}.json.gz")),
    ('gene_burden.json.gz', GS.remote(f"{config['GeneBurden']['outputBucket']}/gene_burden-{timeStamp}.json.gz")),
    ('gene2phenotype.json.gz', GS.remote(f"{config['Gene2Phenotype']['outputBucket']}/gene2phenotype-{timeStamp}.json.gz")),
    ('DDG2P.csv.gz', GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/DDG2P-{timeStamp}.csv.gz")),
    ('EyeG2P.csv.gz', GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/EyeG2P-{timeStamp}.csv.gz")),
    ('SkinG2P.csv.gz', GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/SkinG2P-{timeStamp}.csv.gz")),
    ('CancerG2P.csv.gz', GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/CancerG2P-{timeStamp}.csv.gz")),
    ('SkeletalG2P.csv.gz', GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/SkeletalG2P-{timeStamp}.csv.gz")),
    (config['ChemicalProbes']['probesExcelDump'], GS.remote(f"{config['ChemicalProbes']['inputBucket']}/config['ChemicalProbes']['probesExcelDump']-{timeStamp}.xlsx")),
    (config['ChemicalProbes']['probesMappingsTable'], GS.remote(f"{config['ChemicalProbes']['inputBucket']}/config['ChemicalProbes']['probesMappingsTable']-{timeStamp}.csv")),
    ('intogen.json.gz', GS.remote(f"{config['intOGen']['outputBucket']}/intogen-{timeStamp}.json.gz")),
    ('orphanet.json.gz', GS.remote(f"{config['Orphanet']['outputBucket']}/orphanet-{timeStamp}.json.gz")),
    ('genomics_england.json.gz', GS.remote(f"{config['PanelApp']['outputBucket']}/genomics_england-{timeStamp}.json.gz")),
    ('impc.json.gz', GS.remote(f"{config['IMPC']['evidenceOutputBucket']}/impc-{timeStamp}.json.gz")),
    ('mouse_phenotypes.json.gz', GS.remote(f"{config['IMPC']['phenotypesOutputBucket']}/mouse_phenotypes-{timeStamp}.json.gz")),
    ('progeny.json.gz', GS.remote(f"{config['PROGENy']['outputBucket']}/progeny-{timeStamp}.json.gz")),
    ('slapenrich.json.gz', GS.remote(f"{config['SLAPEnrich']['outputBucket']}/slapenrich-{timeStamp}.json.gz")),
    ('sysbio.json.gz', GS.remote(f"{config['SysBio']['outputBucket']}/sysbio-{timeStamp}.json.gz")),
    #('tep.json.gz', GS.remote(f"{config['TEP']['outputBucket']}/tep-{timeStamp}.json.gz")),
    ('safetyLiabilities.json.gz', GS.remote(f"{config['TargetSafety']['outputBucket']}/safetyLiabilities-{timeStamp}.json.gz")),
    ('chemicalProbes.json.gz', GS.remote(f"{config['ChemicalProbes']['outputBucket']}/chemicalProbes-{timeStamp}.json.gz")),
    ('crispr_screens.json.gz', GS.remote(f"{config['CrisprScreens']['outputBucket']}/crispr_screens-{timeStamp}.json.gz")),
    ('pharmacogenetics.json.gz', GS.remote(f"{config['Pharmacogenetics']['outputBucket']}/cttv012-{timeStamp}_pgkb.json.gz")),
    # PPP specific parsers:
    ('ot_crispr.json.gz', GS.remote(f"{config['OT_CRISPR']['outputBucket']}/ot_crispr-{timeStamp}.json.gz")),
    ('validation_lab.json.gz', GS.remote(f"{config['ValidationLab']['outputBucket']}/validation_lab-{timeStamp}.json.gz")),
    ('encore.json.gz', GS.remote(f"{config['Encore']['outputBucket']}/encore-{timeStamp}.json.gz"))
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
rule baselineExpression:      # Calculate baseline expression data from GTEx V8.
    params:
        gtex_source_data_path = config["baselineExpression"]["gtexSourceDataPath"],
        tissue_name_to_uberon_mapping_path = os.path.join(
            config['global']['curation_repo'],
            config['baselineExpression']['tissueNameToUberonMappingPath']
        ),
        schema = os.path.join(
            config['global']['schema'],
            config['baselineExpression']['schema']
        )
    output:
        "baseline_expression.json.gz"
    log:
        "log/baseline_expression.log"
    shell:
        """
        # In this and the following rules, the exec call redirects the output of all subsequent commands (both STDOUT
        # and STDERR) to the specified log file.
        exec &> {log}
        python modules/baseline_expression/baseline.py \
            --gtex-source-data-path {params.gtex_source_data_path} \
            --tissue-name-to-uberon-mapping-path {params.tissue_name_to_uberon_mapping_path} \
            --output-file-path {output}
        opentargets_validator --schema {params.schema} {output}
        """

rule essentiality:            # Process essentiality data from DepMap.
    params:
        tissue_mapping = f"{config['global']['curation_repo']}/{config['Essentiality']['depmap_tissue_mapping']}",
        schema = f"{config['global']['schema']}/schemas/gene-essentiality.json",
        input_folder = config['Essentiality']['inputBucket']
    output:
        'essentiality.json.gz'
    log:
        'log/essentiality.json.log'
    shell:
        """
        exec &> {log}
        # copy files from bucket to the home folder:
        gsutil -m cp -r "{params.input_folder}/*" ~/
        python modules/Essentiality.py \
            --depmap_input_folder  ~/ \
            --depmap_tissue_mapping {params.tissue_mapping} \
            --output_file {output}
        opentargets_validator --schema {params.schema} {output}
        """

rule cancerBiomarkers:        # Process the Cancers Biomarkers database from Cancer Genome Interpreter.
    input:
        biomarkers_table = GS.remote(config['cancerBiomarkers']['inputAssociationsTable']),
        source_table = GS.remote(config['cancerBiomarkers']['inputSourceTable']),
        disease_table = GS.remote(config['cancerBiomarkers']['inputDiseaseTable']),
        drug_index = GS.remote(config['cancerBiomarkers']['drugIndex'])
    params:
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
    output:
        'cancer_biomarkers.json.gz'
    log:
        'log/cancer_biomarkers.log'
    shell:
        """
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
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
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
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
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

rule projectScore:                  # Process cancer therapeutic targets using CRISPR–Cas9 screens.
    input:
        evidenceFile = GS.remote(config['ProjectScore']['geneScores']),
        cellTypesFile = GS.remote(config['ProjectScore']['cellLinesTable']),
        cell_passport_file = HTTP.remote({config['global']['cell_passport_file']})
    params:
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json",
        uberontoCellLineMapping = f"{config['global']['curation_repo']}/{config['Essentiality']['depmap_tissue_mapping']}"
    output:
        evidenceFile = 'project_score.json.gz'
    log:
        'log/project_score.log'
    shell:
        """
        exec &> {log}
        python modules/ProjectScore.py \
            --evidence_file {input.evidenceFile} \
            --cell_types_file {input.cellTypesFile} \
            --cell_passport_file {input.cell_passport_file} \
            --cell_line_to_uberon_mapping {params.uberontoCellLineMapping} \
            --output_file {output.evidenceFile} 
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule geneBurden:              # Processes gene burden data from various burden analyses
    input:
        azPhewasBinary = GS.remote(config['GeneBurden']['azPhewasBinary']),
        azPhewasQuant = GS.remote(config['GeneBurden']['azPhewasQuantitative']),
        azGenesLinks = GS.remote(config['GeneBurden']['azPhewasGenesLinks']),
        azPhenoLinks = GS.remote(config['GeneBurden']['azPhewasPhenotypesLinks']),
        curation = HTTP.remote(f"{config['global']['curation_repo']}/{config['GeneBurden']['curation']}"),
        genebass = GS.remote(config['GeneBurden']['genebass']),
        finngen = GS.remote(config['GeneBurden']['finngen']),
        finngen_manifest = HTTP.remote(f"{config['GeneBurden']['finngenManifest']}")
    output:
        evidenceFile = "gene_burden.json.gz"
    params:
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
    log:
        'log/gene_burden.log'
    shell:
        """
        exec &> {log}
        python modules/GeneBurden.py \
            --az_binary_data {input.azPhewasBinary} \
            --az_quant_data {input.azPhewasQuant} \
            --az_genes_links {input.azGenesLinks} \
            --az_phenotypes_links {input.azPhenoLinks} \
            --curated_data {input.curation} \
            --genebass_data {input.genebass} \
            --finngen_data {input.finngen} \
            --finngen_manifest {input.finngen_manifest} \
            --output {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule gene2Phenotype:          # Processes four gene panels from Gene2Phenotype
    input:
        ddPanel = HTTP.remote(config['Gene2Phenotype']['webSource_dd_panel']),
        eyePanel = HTTP.remote(config['Gene2Phenotype']['webSource_eye_panel']),
        skinPanel = HTTP.remote(config['Gene2Phenotype']['webSource_skin_panel']),
        cancerPanel = HTTP.remote(config['Gene2Phenotype']['webSource_cancer_panel']),
        cardiacPanel = HTTP.remote(config['Gene2Phenotype']['webSource_cardiac_panel']),
        skeletalPanel = HTTP.remote(config['Gene2Phenotype']['webSource_skeletal_panel'])
    params:
        cacheDir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
    output:
        ddBucket = 'DDG2P.csv.gz',
        eyeBucket = 'EyeG2P.csv.gz',
        skinBucket = 'SkinG2P.csv.gz',
        cancerBucket = 'CancerG2P.csv.gz',
        cardiacBucket = 'CardiacG2P.csv.gz',
        skeletalBucket = 'SkeletalG2P.csv.gz',
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
        cp {input.skeletalPanel} {output.skeletalBucket}
        python modules/Gene2Phenotype.py \
          --panels {input.ddPanel} {input.eyePanel} {input.skinPanel} {input.cancerPanel} {input.cardiacPanel} {input.skeletalPanel} \
          --output_file {output.evidenceFile} \
          --cache_dir {params.cacheDir}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule intogen:                 # Process cohorts and driver genes data from intOGen.
    input:
        inputGenes = GS.remote(config['intOGen']['inputAssociationsTable']),
        inputCohorts = GS.remote(config['intOGen']['inputCohortsTable']),
        diseaseMapping = config['intOGen']['diseaseMapping']
    params:
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
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
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
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
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
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
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
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
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
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
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
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
        evidenceFile = GS.remote(config['SysBio']['inputAssociationsTable']),
        studyFile = GS.remote(config['SysBio']['inputStudyTable'])
    params:
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
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
        schema = f"{config['global']['schema']}/schemas/TEP.json"
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

rule crisprScreens:           # Generating disease/target evidence based on various sources of CRISPR screens.
    params:
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json",
        crispr_brain_mapping = f"{config['global']['curation_repo']}/{config['CrisprScreens']['crispr_brain_mapping']}"
    output:
        'crispr_screens.json.gz'
    log:
        'log/crispr_screens.log'
    shell:
        """
        exec &> {log}
        python modules/crispr_screens.py  \
            --crispr_brain_mapping {params.crispr_brain_mapping} \
            --output {output}
        opentargets_validator --schema {params.schema} {output}
        """

rule chemicalProbes:          # Process data from the Probes&Drugs portal.
    input:
        rawProbesExcel = HTTP.remote(f"https://www.probes-drugs.org/media/download/probes/{config['ChemicalProbes']['probesExcelDump']}"),
        probesXrefsTable = HTTP.remote(f"https://www.probes-drugs.org/media/download/id_mapping/{config['ChemicalProbes']['probesMappingsTable']}")
    params:
        schema = f"{config['global']['schema']}/schemas/chemical_probes.json"
    output:
        rawProbesExcel = config['ChemicalProbes']['probesExcelDump'],
        probesXrefsTable = config['ChemicalProbes']['probesMappingsTable'],
        evidenceFile = 'chemicalProbes.json.gz'
    log:
        'log/chemicalProbes.log'
    shell:
        """
        exec &> {log}
        # Retain the inputs. They will be later uploaded to GCS.
        cp {input.rawProbesExcel} {output.rawProbesExcel}
        cp {input.probesXrefsTable} {output.probesXrefsTable}
        python modules/chemicalProbes.py \
            --probes_excel_path {input.rawProbesExcel} \
            --probes_mappings_path {input.probesXrefsTable} \
            --output {output.evidenceFile}
        opentargets_validator --schema {params.schema} {output.evidenceFile}
        """

rule ot_crispr:               # Generating PPP evidence for OTAR CRISPR screens
    input:
        study_table = GS.remote(config['OT_CRISPR']['config']),
    params:
        data_folder = config['OT_CRISPR']['data_directory'],
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
    output:
        'ot_crispr.json.gz'
    log:
        'log/ot_crispr.log'
    shell:
        """
        exec &> {log}
        
        # Fetching the data from the bucket to the home folder:
        mkdir -p ~/ot_crispr_data
        gsutil -m cp -r "{params.data_folder}/*" ~/ot_crispr_data/
        
        # Call parser script:
        python partner_preview_scripts/ot_crispr.py \
            --study_table {input.study_table} \
            --data_folder ~/ot_crispr_data \
            --output {output}
        
        # Validate schema:
        opentargets_validator --schema {params.schema} {output}
        """

rule encore:                  # Generating PPP evidence for ENCORE
    params:
        data_folder = config['Encore']['data_directory'],
        config = config['Encore']['config'],
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
    output:
        'encore.json.gz'
    input:
        cell_passport_table = HTTP.remote(config['global']['cell_passport_file'], keep_local=True),
    log:
        'log/encore.log'
    shell:
        """
        exec &> {log}
        python partner_preview_scripts/encore_parser.py \
            --output_file {output} \
            --parameter_file {params.config} \
            --data_folder {params.data_folder} \
            --cell_passport_file {input.cell_passport_table}
        opentargets_validator --schema {params.schema} {output}
        """

rule validation_lab:          # Generating PPP evidence for Validation Lab
    params:
        data_folder = config['ValidationLab']['data_directory'],
        config = config['ValidationLab']['config'],
        schema = f"{config['global']['schema']}/schemas/disease_target_evidence.json"
    input:
        cell_passport_table = HTTP.remote(config['global']['cell_passport_file'], keep_local=True),
    output:
        'validation_lab.json.gz'
    log:
        'log/validation_lab.log'
    shell:
        """
        exec &> {log}
        python partner_preview_scripts/ValidationLab.py \
            --parameter_file {params.config} \
            --data_folder {params.data_folder} \
            --cell_passport_file {input.cell_passport_table} \
            --output_file {output}
        opentargets_validator --schema {params.schema} {output}
        """

envvars:
    "OPENAI_API_KEY"

rule Pharmacogenetics:                     # Generating pharmacogenetics evidence
    input:
        evidenceFile = GS.remote(config['Pharmacogenetics']['evidence']),
    params:
        schema = f"{config['global']['schema']}/{config['Pharmacogenetics']['schema']}",
        phenotypes = f"{config['global']['curation_repo']}/{config['Pharmacogenetics']['phenotypes']}",
        phenotypes_optional_output = "pharmgkb_phenotypes.json",
        cache_dir = config['global']['cacheDir'],
        openai_api_key = os.environ['OPENAI_API_KEY']
    output:
        evidence = "pharmacogenetics.json.gz"
    log:
        'log/pharmacogenetics.log'
    shell:
        """
        exec &> {log}
        python modules/Pharmacogenetics.py \
            --pharmgkb_evidence_path {input.evidenceFile} \
            --extracted_phenotypes_path {params.phenotypes} \
            --openai_api_key {params.openai_api_key} \
            --output_evidence_path {output.evidence} \
            --output_phenotypes_path {params.phenotypes_optional_output} \
            --cache_dir {params.cache_dir}
        opentargets_validator --schema {params.schema} {output.evidence}
        """

rule targetSafety:            # Process data from different sources that describe target safety liabilities.
    input:
        toxcast = GS.remote(config['TargetSafety']['toxcast']),
        aopwiki = GS.remote(config['TargetSafety']['aopwiki']),
        pharmacogenetics = rules.Pharmacogenetics.output.evidence,
        brennan = GS.remote(config['TargetSafety']['brennan']),
    params:
        ae = f"{config['global']['curation_repo']}/{config['TargetSafety']['adverseEvents']}",
        sr = f"{config['global']['curation_repo']}/{config['TargetSafety']['safetyRisk']}",
        cache_dir = config['global']['cacheDir'],
        schema = f"{config['global']['schema']}/schemas/target_safety.json"
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
            --aopwiki {input.aopwiki} \
            --pharmacogenetics {input.pharmacogenetics} \
            --brennan {input.brennan} \
            --cache_dir {params.cache_dir} \
            --output {output}
        opentargets_validator --schema {params.schema} {output}
        """
