from datetime import datetime

# Current date in YYYY-MM-DD format:
today = datetime.now().strftime("%Y-%m-%d")

# Configuration is read from the config yaml:
configfile: 'configuration.yaml'

## Running genetics portal evidence generation. 
## * input files under opentargets-genetics project area
## * output is generated in opentargets-platform bucket.    
rule GeneticsPortal:
    params:
        locus2gene=config['GeneticsPortal']['locus2gene'],
        toploci=config['GeneticsPortal']['toploci'],
        study=config['GeneticsPortal']['study'],
        variantIndex=config['GeneticsPortal']['variantIndex'],
        ecoCodes=config['GeneticsPortal']['ecoCodes'],
        outputFile=config['GeneticsPortal']['outputFile'],
        threshold=config['GeneticsPortal']['threshold']
    output:
        (config['GeneticsPortal']['outputFile'] % today)
    threads: 32
    log:
        f'logs/GeneticsPortal.log'
    shell:
        '''
        python modules/GeneticsPortal_prepare_data.py \
            --locus2gene {params.locus2gene} \
            --toploci {params.toploci} \
            --study {params.study} \
            --variantIndex {params.variantIndex} \
            --ecoCodes {params.ecoCodes} \
            --outputFile {output}
        '''

