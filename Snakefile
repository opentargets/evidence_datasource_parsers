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


