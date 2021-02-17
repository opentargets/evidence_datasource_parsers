from datetime import datetime
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

# Current time in YYYY-MM-DD-hh-mm format:
timeStamp = datetime.now().strftime("%Y-%m-%d")

# Configuration is read from the config yaml::
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
        webSource=config['ClinGen']['webSource'],
        tempFile=f"/tmp/ClinGen-Gene-Disease-Summary-{timeStamp}.csv"
    output:
        GS.remote(f"{config['ClinGen']['inputBucket']}/ClinGen-Gene-Disease-Summary-{timeStamp}.csv")
    log:
        GS.remote(logFile)
    shell:
        """
        curl  {params.webSource} > {params.tempFile}
        gsutil cp {params.tempFile} {output}
        """

rule ClinGen:
    input:
        GS.remote(f"{config['ClinGen']['inputBucket']}/ClinGen-Gene-Disease-Summary-{timeStamp}.csv")
    output:
        evidenceFile=GS.remote(f"{config['ClinGen']['outputBucket']}/ClinGen-{timeStamp}.json"),
        unmappedFile=GS.remote(f"{config['ClinGen']['outputBucket']}/ClinGen-Gene-unmapped-{timeStamp}.lst")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/ClinGen.py -i {input} -o {output.evidenceFile} -u {output.unmappedFile}
        """
'''
rule PheWAS:
    input:
        inputFile = GS.remote(f"{config['PheWAS']['inputBucket']}/phewas-catalog-19-10-2018.csv")
        consequencesFile = GS.remote(f"{config['PheWAS']['inputBucket']}/phewas_w_consequences.csv")
        diseaseMapping = GS.remote(config['PheWAS']['diseaseMapping'])
    output:
        evidenceFile=GS.remote(f"{config['PheWAS']['outputBucket']}/phewas_catalog-{timeStamp}.json.gz"),
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/PheWAS.py -i {input.inputFile} -c {input.consequencesFile} -d {input.diseaseMapping} -o {output.evidenceFile}
        """

rule SLAPEnrich:
    input:
        inputFile = GS.remote(f"{config['SLAPEnrich']['inputBucket']}/slapenrich_opentargets-21-12-2017.tsv")
        diseaseMapping = "resources/cancer2EFO_mappings.tsv"
    output:
        evidenceFile=GS.remote(f"{config['PheWAS']['outputBucket']}/slapenrich-{timeStamp}.json.gz"),
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/SLAPEnrich.py -i {input.inputFile} -d {input.diseaseMapping} -o {output.evidenceFile}
        """

rule PROGENy:
    input:
        inputFile = GS.remote(f"{config['PROGENy']['inputBucket']}/progeny_normalVStumor_opentargets.txt")
        diseaseMapping = "resources/cancer2EFO_mappings.tsv"
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
## Running genetics portal evidence generation. 
## * input files under opentargets-genetics project area
## * output is generated in opentargets-platform bucket.    
# rule GeneticsPortal:
#     params:
#         locus2gene='/temp/genetics/' + config['GeneticsPortal']['locus2gene'].split('/')[-1],
#         toploci='/temp/genetics/' + config['GeneticsPortal']['toploci'].split('/')[-1],
#         study='/temp/genetics/' + config['GeneticsPortal']['study'].split('/')[-1],
#         variantIndex='/temp/genetics/' + config['GeneticsPortal']['variantIndex'].split('/')[-1],
#         ecoCodes='/temp/genetics/' + config['GeneticsPortal']['ecoCodes'].split('/')[-1],
#         outputFile=config['GeneticsPortal']['outputFile'],
#         threshold=config['GeneticsPortal']['threshold']
#     output:
#         (config['GeneticsPortal']['outputFile'] % timeStamp)
#     threads: 32
#     log:
#         f'logs/GeneticsPortal.log'
#     shell:
#         '''
#         python modules/GeneticsPortal_prepare_data.py \
#             --locus2gene {params.locus2gene} \
#             --toploci {params.toploci} \
#             --study {params.study} \
#             --variantIndex {params.variantIndex} \
#             --ecoCodes {params.ecoCodes} \
#             --outputFile {output}
#         '''


