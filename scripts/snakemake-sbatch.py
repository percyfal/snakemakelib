#!/usr/bin/env python
# Copyright (C) 2015 by Per Unneberg
"""
snakemake-sbatch.py

"""
import os
import sys
import argparse
import subprocess
from snakemake.utils import read_job_properties

# attributes = {
#     'dep_str': dependencies,
#     'job_name': 'snakemake_{0}'.format(job_properties['rule']),
#     'partition': rule_conf['partition'],
#     'cores': rule_conf['cores'],
#     'account': self.config['sbatch_general']['account'],
#     'log_file': self.ofiles[0] + '-slurm.out' if len(self.ofiles) > 0 else 'snakemake-{0}-slurm.out'.format(self.rule),
#                     'extra_parameters': rule_conf.get('extra_parameters', "")
#             }
#             sbatch_cmd = """sbatch --output={log_file} {dep_str} -A {account} -p {partition} -n {cores} -t {days}-{hours}:{minutes}:00 \
#                     -J {job_name} {extra_parameters} {sbatch_job_path} \
#                     '{script_name}'""".format(**attributes)
# do something useful with the threads
# threads = job_properties[threads]

# print (job_properties)
# os.system("qsub -t {threads} {script}".format(threads=threads, script=jobscript))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dependencies", nargs="*", help="{{dependencies}} string given by snakemake\n")
    parser.add_argument("snakescript", help="Snakemake generated shell script with commands to execute snakemake rule\n")
    args = parser.parse_args()
    try:
        jobscript = sys.argv[1]
        job_properties = read_job_properties(jobscript)
        print (sys.argv[1], file=sys.stderr)
        print ("snakescript: {}".format(args.snakescript), file=sys.stderr)
        print ("Dependencies: {}".format(args.dependencies), file=sys.stderr)
        threads = job_properties.get('threads', 1)
        sbatch_cmd = """sbatch -t 00:10:00 -p devel -n {threads} -A b2015134 {script}""".format(threads=threads, script=jobscript)
        print(sbatch_cmd, file=sys.stderr)
        popenrv = subprocess.Popen(sbatch_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).communicate()
    except:
        print (jobscript, job_properties)
