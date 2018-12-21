shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

localrules: all
# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs

# load cluster config file
CLUSTER = json.load(open(config['CLUSTER_JSON']))
FILES = json.load(open(config['SAMPLES_JSON']))

import csv
import os

SAMPLES = sorted(FILES.keys())
print ("SAMPLES")
print (SAMPLES)

## list all samples by sample_name and sample_type
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        MARK_SAMPLES.append(sample + "-" + sample_type)
print ("MARK_SAMPLES")
print (MARK_SAMPLES)

# which sample_type is used as control for calling peaks: e.g. Input, IgG...
CONTROL = config["control"]
CONTROLS = [sample for sample in MARK_SAMPLES if CONTROL in sample]
CASES = [sample for sample in MARK_SAMPLES if CONTROL not in sample]
print ("CONTROL")
print (CONTROL)

print ("CONTROLS")
print (CONTROLS)

print ("CASES")
print (CASES)

## multiple samples may use the same control input/IgG files
CONTROLS_UNIQUE = list(set(CONTROLS))
print (CONTROLS_UNIQUE)

## list BAM files
CONTROL_BAM = expand("01.bam/{sample}.bam", sample=CONTROLS_UNIQUE)
CASE_BAM = expand("01.bam/{sample}.bam", sample=CASES)

print ("CONTROL_BAM")
print (CONTROL_BAM)

print ("CASE_BAM")
print (CASE_BAM)

ALL_SAMPLES = CASES + CONTROLS_UNIQUE
print ("ALL_SAMPLES")
print (ALL_SAMPLES)

ALL_BAM = []
ALL_BAM = expand("01.bam/{sample}.bam", sample = ALL_SAMPLES)
print("ALL_BAM")
print(ALL_BAM)

ALL_VCF = []

for case in CASES:
    print("*************new")
    print("case:")
    print(case)  
#    print ("CASES:")
#    print (CASES)  
    sample = "-".join(case.split("-")[0:-1])
#    print("sample:")
#    print(sample)  
    control = sample + "-" + CONTROL
    print("control:")
    print(control)
    print ("CONTROLS:")
    print (CONTROLS)
    if control in CONTROLS:
        print("#########")
        print("case:")
        print(case) 
        print("control:")
        print(control)
        ALL_VCF.append("02.vcf/1c_somatic_m2_{}_vs_{}.vcf.gz".format(case, control))
        #ALL_VCF.append("02.vcf/2c_tumor_normal_m2_{}_vs_{}.bam".format(case, control))

print(ALL_VCF)   
TARGETS = []
TARGETS.extend(ALL_BAM)
TARGETS.extend(ALL_VCF)

localrules: all
rule all:
    input: TARGETS


## get a list of fastq.gz files for the same mark, same sample
def get_bam(wildcards):
    sample = "-".join(wildcards.sample.split("-")[0:-1])
    mark = wildcards.sample.split("-")[-1]
    return FILES[sample][mark]

## link,
rule ln_cp:
    input: get_bam
    output: "01.bam/{sample}.bam"
    log: "00.log/{sample}_cp"
    shell: 
       """
       ln -s {input} {output}
       """

rule all_m2:
    input : "01.bam/{control}.bam", "01.bam/{case}.bam" 
    output: "02.vcf/1c_somatic_m2_{case}_vs_{control}.vcf.gz", "02.vcf/2c_tumor_normal_m2_{case}_vs_{control}.bam"
    log: "00.log/{case}_vs_{control}_m2.log"
    params:
            tumor = "{case}",  normal = "{control}"
    shell:
        """
        gatk --java-options "-Xmx20g" Mutect2 \
        -R /data/exx/DNAseq/test/TEST2/vcf/ucsc.hg19.fasta \
        -I {input[1]} \
        -I {input[0]} \
        -tumor {params.tumor} \
        -normal {params.normal} \
        -pon /data/exx/DNAseq/test/TEST2/vcf/dbsnp_138.hg19.vcf \
        --germline-resource /data/exx/DNAseq/test/TEST2/vcf/dbsnp_138.hg19.vcf \
        --af-of-alleles-not-in-resource 0.0000025 \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        -O {output[0]} \
        -bamout {output[1]} 2> {log}
        """



