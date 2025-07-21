####this pipeline is based on the following link
####https://artic.network/ncov-2019
####https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
import os
rootdir = '/home/groupdirs/epibiocore/epishare/antarikshTyagi/altru/artic/G2021_65redo/G2021_65articredo/20210723_1942_X2_FAQ78961_83de785a/fastq_pass'
SAMPLES=os.listdir(rootdir)
print(SAMPLES)

sch_dir="/home/antariksh.tyagi/software/conda_envs/artic-ncov2019/primer_schemes"
seq_summ="/home/groupdirs/epibiocore/epishare/antarikshTyagi/altru/artic/G2021_65redo/G2021_65articredo/20210723_1942_X2_FAQ78961_83de785a/sequencing_summary_FAQ78961_4decf7a6.txt"
fast5_dir="/home/groupdirs/epibiocore/epishare/antarikshTyagi/altru/artic/G2021_65redo/G2021_65articredo/20210723_1942_X2_FAQ78961_83de785a/fast5_pass"
rule all:
    input:
        expand("minion_out/G2021_65_{sample}.minion.log.txt", sample = SAMPLES),
        expand("minion_out/G2021_65_{sample}.primertrimmed.rg.sorted.bam", sample = SAMPLES),
        expand("genomecov/{sample}.read_depth.tsv", sample = SAMPLES),
        expand("coverage_plots/{sample}_coverage_plot.png", sample = SAMPLES),
        "coverage_plots/combined_coverage_plot.png",
        "coverage_plots/combined_coverage_stats.txt"

rule guppyplex:
    input: 
        "/home/groupdirs/epibiocore/epishare/antarikshTyagi/altru/artic/G2021_65redo/G2021_65articredo/20210723_1942_X2_FAQ78961_83de785a/fastq_pass/{sample}"
    output:
        "guppy_fastq/G2021_65_{sample}.fastq"
    resources: time_min=600, mem_mb=5000, cpus=2
    params: partition = "talon-fat"
    shell:
        "artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory {input} --output {output}"

rule minion:
    input:
        fastq=rules.guppyplex.output,
    output:
        log="minion_out/G2021_65_{sample}.minion.log.txt",
        PTB="minion_out/G2021_65_{sample}.primertrimmed.rg.sorted.bam"
    resources: time_min=600, mem_mb=5000, cpus=2
    params: partition = "talon-fat", samplename = "minion_out/G2021_65_{sample}"
    shell:
        "artic minion --normalise 200 --threads {resources.cpus} --scheme-directory {sch_dir} "
        "--read-file {input.fastq} --fast5-directory {fast5_dir} "
        "--sequencing-summary {seq_summ} nCoV-2019/V3 {params.samplename}"

rule genomecov:
    input:
        rules.minion.output.PTB
    output:
        "genomecov/{sample}.read_depth.tsv"
    resources: time_min=600, mem_mb=5000, cpus=2
    params: partition = "talon-fat"
    shell:
        "bedtools genomecov -d -ibam {input} > {output}"

rule plotCoverage:
    input:
        rules.minion.output.PTB
    output:
        "coverage_plots/{sample}_coverage_plot.png"
    resources: time_min=600, mem_mb=5000, cpus=2
    params: partition = "talon-fat"
    shell:
        "plotCoverage -p {resources.cpus} -b {input} -o {output}"

rule plotCoverage_all:
    input:
        expand("minion_out/G2021_65_{sample}.primertrimmed.rg.sorted.bam", sample=SAMPLES)
    output:
        plot = "coverage_plots/combined_coverage_plot.png",
        stats = "coverage_plots/combined_coverage_stats.txt"
    resources: time_min=600, mem_mb=5000, cpus=24
    params: partition = "talon-fat"
    shell:
        "plotCoverage -p {resources.cpus} -b {input} -o {output.plot} > {output.stats} "
 

