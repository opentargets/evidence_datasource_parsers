from datetime import datetime
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

# Current time in YYYY-MM-DD-hh-mm format:
timeStamp = datetime.now().strftime("%Y-%m-%d")

# Configuration is read from the config yaml:
configfile: 'configuration.yaml'

# schema version is now hardcoded, version will be read from the command line later:
schemaFile = config['global']['schema_file']
logFile = f"{config['global']['logDir']}/evidence_parser.{timeStamp}.log"


##
## Running all rules
##
rule all:
    input:
        GS.remote(f"{config['ClinGen']['inputBucket']}/ClinGen-Gene-Disease-Summary-{timeStamp}.csv"),
        evidenceFile=GS.remote(f"{config['ClinGen']['outputBucket']}/ClinGen-{timeStamp}.json")

##
## Running ClinGen parser
##
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

rule ClinGen:
    input:
        f"tmp/ClinGen-Gene-Disease-Summary-{timeStamp}.csv"
    output:
        evidenceFile=GS.remote(f"{config['ClinGen']['outputBucket']}/ClinGen-{timeStamp}.json"),
        unmappedFile=GS.remote(f"{config['ClinGen']['outputBucket']}/ClinGen-Gene-unmapped-{timeStamp}.lst")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/ClinGen.py -i {input} -o {output.evidenceFile} -u {output.unmappedFile}
        """

##
## Running Genetics Portal parser
##

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

rule GeneticsPortal:
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

##
## Running PheWAS parser
##

rule fetchPheWAS:
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

rule PheWAS:
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

##
## Running SLAPEnrich parser
##

rule fetchSLAPEnrich:
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

rule SLAPEnrich:
    input:
        inputFile=f"tmp/slapenrich-{timeStamp}.csv",
        diseaseMapping=f"tmp/cancer2EFO_mappings-{timeStamp}.tsv"
    output:
        evidenceFile=GS.remote(f"{config['SLAPEnrich']['outputBucket']}/slapenrich-{timeStamp}.json.gz"),
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/SLAPEnrich.py \
            --inputFile {input.inputFile} \
            --diseaseMapping {input.diseaseMapping} \
            --outputFile {output.evidenceFile} \
        """

##
## Running Gene2Phenotype parser
##
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

rule Gene2Phenotype:
    input:
        ddPanel=f"tmp/DDG2P-{timeStamp}.csv.gz",
        eyePanel=f"tmp/EyeG2P-{timeStamp}.csv.gz",
        skinPanel=f"tmp/SkinG2P-{timeStamp}.csv.gz",
        cancerPanel=f"tmp/CancerG2P-{timeStamp}.csv.gz"
    output:
        evidenceFile=GS.remote(f"{config['Gene2Phenotype']['outputBucket']}/gene2phenotype-{timeStamp}.json"),
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

'''

rule PROGENy:
    input:
        inputFile = GS.remote(f"{config['PROGENy']['inputBucket']}/progeny_normalVStumor_opentargets.txt"),
        diseaseMapping = "resources/cancer2EFO_mappings.tsv",
        pathwayMapping = "resources/pathway2Reactome_mappings.tsv"
    output:
        evidenceFile=GS.remote(f"{config['PROGENy']['outputBucket']}/progeny-{timeStamp}.json.gz"),
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/PROGENY.py -i {input.inputFile} -d {input.diseaseMapping} -p {input.pathwayMapping} -o {output.evidenceFile}
        """

rule intOGen:
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
'''