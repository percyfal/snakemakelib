# -*- snakemake -*-
#
# CENTIPEDE: Transcription factor footprinting and binding site prediction
# install.packages("CENTIPEDE", repos="http://R-Forge.R-project.org") 
#  
# http://centipede.uchicago.edu/
#
include: '../settings.rules'

config_default = {
    'bio.ngs.motif.centipede' : {
        'options' : '',
    },
}

update_config(config_default, config)
config = config_default

