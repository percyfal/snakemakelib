# -*- snakemake -*- 
include: "settings.rules"

rule comp_utils_link:
    """Link source to target from srcpath to tgtpath"""
    input: os.path.join("{srcpath}", "{fn}")
    output: os.path.join("{tgtpath}", "{fn}")
    shell: "ln"
