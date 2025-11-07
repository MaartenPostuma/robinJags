configfile: "config.yaml"
import pandas as pd
import os
import random
from snakemake.utils import Paramspace
paramspace = Paramspace(pd.read_csv("variablesToTest.tsv", sep="\t"))

rule all:
    input:
        scripts=expand("{path}/scripts/{params}.R",path=config["outputDir"],params=paramspace.instance_patterns),
        RDS=expand("{path}/output/{params}.Rds",path=config["outputDir"],params=paramspace.instance_patterns),
        png=expand("{path}/output/{params}.png",path=config["outputDir"],params=paramspace.instance_patterns),
        jag=expand("{path}/scripts/{params}.jag",path=config["outputDir"],params=paramspace.instance_patterns),
        report=expand("{path}/report.html",path=config["outputDir"])


#test2
rule rename:
    input:
        file= "src/defaultJagsScript.R"
    output:
        script=expand("{path}/scripts/{params}.R",path=config["outputDir"],params=paramspace.wildcard_pattern)
    resources:
        mem_mb= 1000,
        runtime= 10,
        cpus_per_task= 1
    shell:
        """
        cat src/defaultJagsScript.R | sed "s/Order/{wildcards.Order}/g" | sed "s/Variable/{wildcards.Variable}/g" > {output.script}
        """

rule runJags:
    input:
        script=expand("{path}/scripts/{params}.R",path=config["outputDir"],params=paramspace.wildcard_pattern)
    output:
        RDS=expand("{path}/output/{params}.Rds",path=config["outputDir"],params=paramspace.wildcard_pattern),
        png=expand("{path}/output/{params}.png",path=config["outputDir"],params=paramspace.wildcard_pattern),
        jag=expand("{path}/scripts/{params}.jag",path=config["outputDir"],params=paramspace.wildcard_pattern)
    conda:
        "src/jags.yaml"
    params:
        input1=config["countData"],
        input2=config["hourlyData"]
    resources:
        mem_mb= 10000,
        runtime= 480,
        cpus_per_task= 3
    shell:
        """
        Rscript {input.script} {params.input1} {params.input2} {output.RDS} {output.png} {output.jag}
        """

rule makeReport:
    input:
        RDS=expand("{path}/output/{params}.Rds",path=config["outputDir"],params=paramspace.instance_patterns)
    output:
        report=expand("{path}/report.html",path=config["outputDir"])
    params:
        outputDir=config["outputDir"]
    conda:
        "src/jags.yaml"
    resources:
        mem_mb= 10000,
        runtime= 480,
        cpus_per_task= 1
    shell:
        """
        R -e "rmarkdown::render('src/report.Rmd',output_file='report.html',params=list(args=c('{params.outputDir}')))"
        mv src/report.html {output.report}
        """