dir_data='raw_data'
dir_ref, dir_reads, dir_control = 'reference', 'reads', 'control'
dir_fastqc='fastqc'

DIR_DICT={'reference': f'{dir_data}/{dir_ref}/',
          'reads': f'{dir_data}/{dir_reads}/',
          'control': f'{dir_data}/{dir_control}/'}
NAMES_DICT={'reference': ['new_dyn_soft_filtered_transcripts.fasta'],
          'reads': ['Primary_polyp_AZK_1_GAGTGGAT_L002_R1_001.fastq.bz2', 'Primary_polyp_AZK_2_ACTGATAT_L002_R1_001.fastq.bz2', 'Primary_polyp_AZK_3_ATTCCTTT_L002_R1_001.fastq.bz2'],
          'control':['Primary_polyp_DMSO_1_GTGGCCTT_L002_R1_001.fastq.bz2', 'Primary_polyp_DMSO_2_GTTTCGGA_L002_R1_001.fastq.bz2', 'Primary_polyp_DMSO_3_CGTACGTA_L002_R1_001.fastq.bz2']}

make_sorting=True

# rule reads_download:
#     output:
#         "{dir1}/{sample}-{n}.fastq.gz"
#     params:
#         URL=lambda wcs: URL_DICT[wcs.sample][int(wcs.n) - 1]  # https://stackoverflow.com/a/56043556
#     shell:
#         "wget -O {output} {params.URL}"

# rule reference_download:
#     output:
#         "{reference}.fasta"
#     shell:
#         "efetch -db nucleotide -id KF848938.1 -format fasta > {wildcards.reference}.fasta"

# rule reference_uzip:
#     input:
#         "{reference}.fasta.gz"
#     output:
#         "{reference}.fasta"
#     shell:
#         "gunzip -c {input} > {output}"

rule reads_unzip:
    output:
        f'{dir_data}/{dir_reads}/{{sample}}-{{n}}.fastq'
        # f'{dir_data}/{dir_control}/{{sample}}-{{n}}.fastq'
    params:
        DIR=lambda wcs: DIR_DICT[wcs.sample],
        FILE=lambda wcs: NAMES_DICT[wcs.sample][int(wcs.n) - 1]  # https://stackoverflow.com/a/56043556
    shell:
        'bzip2 -vcfdk {params.DIR}{params.FILE} > {output}'

rule inspect_quality:
    input:
        rules.reads_unzip.output
    output:
        "{dir_fastqc}/{sample}-{n}_fastqc.html",
        "{dir_fastqc}/{sample}-{n}_fastqc.zip"
    shell:
        "fastqc -o ./{wildcards.dir_fastqc} -t 16 {input}"

# rule bwa_index:
#     input:
#         "{reference}.fasta"
#     log:
#         "logs/bwa_logs/bwa_index_{reference}.log"
#     output:
#         "{reference}.fasta.amb",
#         "{reference}.fasta.ann",
#         "{reference}.fasta.bwt",
#         "{reference}.fasta.pac",
#         "{reference}.fasta.sa"
#     shell:
#         "bwa index {input} &> {log}"
        
# rule bwa_align:
#     input:        
#         rules.bwa_index.output,
#         ref="{reference}.fasta",
#         reads=lambda wcs: "{dir}/{{sample}}-{{n}}.fastq.gz".format(dir=DIR_DICT.get(wcs.sample, "reads"))  # https://stackoverflow.com/a/67944385
#     threads: 8
#     params:
#         readgroup=r"@RG\tID:1\tSM:{sample}-{n}\tPL:illumina"
#     log:
#         bwa_log="logs/bwa_logs/bwa_sort_{reference}_{sample}-{n}.log",
#         smt_log="logs/flagstat_{reference}_{sample}-{n}.log"
#     output:
#         "{reference}_{sample}-{n}_unsorted.bam"
#     shell:
#         "bwa mem -t {threads} -R '{params.readgroup}' {input.ref} {input.reads} 2> {log.bwa_log} | samtools view -S -b - | tee {output} | samtools flagstat - &> {log.smt_log}"
#         # 1. `samtools` нужно указывать `-` для входных данных, если они поступают из стандартного потока ввода, а не передаются файлом.
#         # 2. Чтобы перенаправить результат одновременно в и файл, и в стандартный поток вывода, нужно использовать команду `tee`.

# rule bam_sort:
#     input:
#         rules.bwa_align.output
#     output:
#         protected("{reference}_{sample}-{n}_sorted.bam")
#     threads: 8
#     shell:
#         "samtools sort --threads {threads} {input} > {output}"

# rule calculate_coverage:
#     input:
#         rules.bam_sort.output
#     output:
#         "logs/coverage_{reference}_{sample}-{n}.txt"
#     shell:
#         "samtools depth -a {input} > {output}"

# rule bam_index:
#     input:
#         rules.bam_sort.output
#     output:
#         "{reference}_{sample}-{n}_sorted.bam.bai"
#     threads: 8
#     shell:
#         "samtools index --threads {threads} {input} > {output}"

# rule variant_calling:
#     input:
#         bam=rules.bam_sort.output,
#         ref="{reference}.fasta"
#     log:
#         vr_log="logs/varscan_{reference}_{sample}-{n}_{freq}.log",
#         vr_res_parse="logs/varscan_{reference}_{sample}-{n}_{freq}_parse.txt"        
#     output:
#         "{reference}_{sample}-{n}_{freq}_VarScan_results.vcf"
#     shell:
#         """samtools mpileup -f {input.ref} {input.bam} -d 0 | varscan mpileup2snp --min-var-freq {wildcards.freq} --variants --output-vcf 1 1> >(tee {output}) 2> {log.vr_log} | awk 'NR>24 {{print $1, $2, $4, $5, $10}}' > {log.vr_res_parse}"""
#         # флаг -d 0 устанавливает количество обрабатываемых ридов для каждой позиции на максимальное значение 214748364: http://www.htslib.org/doc/samtools-mpileup.html#OPTIONS
#         # перенаправляем stdout и stderr в 2 разных файла: `1> {output} 2> {log}`
#         # используем [process substitution](https://stackoverflow.com/a/692407) чтобы передать в `tee` только stdout, а stderr перенаправить в лог: 1> >(tee {output}) 2> {log.vr_log}

# rule filter_awk_freq:
#     input:
#       "{reference}_{sample}-{n}_0.001_VarScan_results.vcf"
#     output:
#       "{reference}_{sample}-{n}_freq.csv"
#     shell:
#         """cat {input} | awk 'NR>24 {{print $10}}' | awk 'BEGIN {{FS=":";}}{{print $7}}' | awk -F'\t' -v OFS='\t' 'NR == 0 {{print $0; next}}{{print (NR), $0}}' > {output}"""
        
# rule filter_awk_01:
#     input:
#       "{reference}_{sample}-{n}_0.001_VarScan_results.vcf"
#     output:
#       "{reference}_{sample}-{n}_variants_01.csv"
#     shell:
#         """cat {input} | awk 'NR>24 {{print $1, $2, $4, $5}}' | awk -F'\t' -v OFS='\t' 'NR == 0 {{print $0; next}}{{print (NR), $0}}' > {output}"""

# rule merge:
#     input:
#       file1 = rules.filter_awk_01.output,
#       file2 = rules.filter_awk_freq.output
#     output:
#       "{reference}_{sample}-{n}_merge.csv"
#     shell:
#         """join {input.file1} {input.file2} | awk '{{sub("%", "", $6) ; print}}' > {output}"""
