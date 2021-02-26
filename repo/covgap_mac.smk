
#validate(config, "config.schema.yaml")
from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(config["read_dir"]) if isfile(join(config["read_dir"], f))]
samples= sorted(set([s.strip("R[1-2]_fastq.gz") for s in onlyfiles]))

READD= os.path.abspath(config["read_dir"])
#clone the necessary files in there
import os
if (os.path.exists(os.path.abspath('./config/'))):
    print("\n Configuration file and dependencies already present, proceeding.. \n")
else:
    print("\n Configuration file and dependencies missing, starting the download from github.. \n")
    os.system('git clone https://github.com/appliedmicrobiologyresearch/covgap/ temp_soft/')
    os.system('mkdir config/')
    os.system('cp -r temp_soft/repo/config/* config/')
    os.system('cp temp_soft/repo/settings.yaml .')
    os.system('rm -rf temp_soft')
    os.system('chmod 775 config/dependencies/*')
    os.system('tar -zxvf config/envs/VariantBam.tar.gz')
configfile: "settings.yaml"
rule all:
    input:
        expand("result/{sample}/Mapping/{sample}.alignment.sam", sample=samples),
        expand("result/{sample}/Mapping/{sample}.alignment.bam", sample=samples),
        expand("result/{sample}/Mapping/{sample}.alignment.removed.duplicates.bam", sample=samples),
        expand("result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.reads.only.sam", sample=samples),
        expand("result/{sample}/Mapping/{sample}.alignment.removed.duplicates.unmapped.reads.only.sam", sample=samples),
        expand("result/{sample}/Mapping/{sample}.alignment.stats.tab", sample=samples),
        expand("result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.reads.only.sorted.bam", sample=samples),
        expand("result/{sample}/Mapping/{sample}.alignment.removed.duplicates.unmapped.reads.only.sorted.bam", sample=samples),
        expand("result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.1000trimmed.sam", sample=samples),
        expand("result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.1000trimmed.sorted.bam", sample=samples),
        expand("result/{sample}/Mapping/{sample}.average.coverage.tab", sample=samples),
        expand("result/{sample}/variantcall/{sample}.vcf", sample=samples),
        expand("result/{sample}/variantcall/{sample}.changes", sample=samples),
        expand("result/{sample}/variantcall/{sample}.depth.filtered.vcf", sample=samples),
        expand("result/{sample}/variantcall/{sample}.final.vcf", sample=samples),
        expand("result/{sample}/variantcall/{sample}.afM0.vcf", sample=samples),
        expand("result/{sample}/variantcall/{sample}.minority_alleles.vcf", sample=samples),
        expand("result/{sample}/variantcall/{sample}.minority_alleles_report.tech", sample=samples),
        expand("result/{sample}/variantcall/{sample}.annotated.variants.vcf", sample=samples),
        expand("result/{sample}/variantcall/{sample}.annotated.minorityvariants.vcf", sample=samples),
        expand("result/{sample}/variantcall/{sample}.lessthan50.bed", sample=samples),
        expand("result/{sample}/variantcall/{sample}.finalmask.bed", sample=samples),
        expand("result/{sample}/variantcall/{sample}.old.consensus.fasta", sample=samples),
        expand("result/{sample}/variantcall/{sample}.Nstats.tab", sample=samples),
        expand("result/{sample}/variantcall/{sample}.classification.tab",sample=samples)

rule length_trim_1:
    threads: workflow.cores * 0.75
    params:
        sw=config["QC_sliding_window"],
        phred=config["QC_phred_score"],
        minlen=config["QC_min_len"]
    input:
        r1=os.path.join(READD,"{sample}_R1.fastq.gz"),
        r2=os.path.join(READD,"{sample}_R2.fastq.gz")
    output:
        firsttrim_R1="result/{sample}/trimmomatic/r1.firsttrimmed.fastq.gz",
        firsttrimNP_R1="result/{sample}/trimmomatic/r1.firsttrimmed.not-paired.fastq.gz",
        firsttrim_R2="result/{sample}/trimmomatic/r2.firsttrimmed.fastq.gz",
        firsttrimNP_R2="result/{sample}/trimmomatic/r2.firsttrimmed.not-paired.fastq.gz",
        firsttrimLOG="result/{sample}/trimmomatic/read_trimm1_info"
    conda:
        "config/envs/Java_related_mac.yaml" 
    shell:
        "trimmomatic PE -threads {threads} -phred33 {input.r1} {input.r2} {output.firsttrim_R1} {output.firsttrimNP_R1} {output.firsttrim_R2} {output.firsttrimNP_R2} SLIDINGWINDOW:{params.sw}:{params.phred} MINLEN:{params.minlen} 2> {output.firsttrimLOG}"

rule adaptor_trim:
    input:
        trimmed1_R1="result/{sample}/trimmomatic/r1.firsttrimmed.fastq.gz",
        trimmed1_R2="result/{sample}/trimmomatic/r2.firsttrimmed.fastq.gz"
    params:
        adapters=config["adapters"]
    threads: workflow.cores * 0.75
    output:
        adapttrim_R1="result/{sample}/trimmomatic/r1.no_adaptors.fastq.gz",
        adapttrimNP_R1="result/{sample}/trimmomatic/r1.no_adaptors.not-paired.fastq.gz",
        adapttrim_R2="result/{sample}/trimmomatic/r2.no_adaptors.fastq.gz",
        adapttrimNP_R2="result/{sample}/trimmomatic/r2.no_adaptors.not-paired.fastq.gz",
        adapttrimLOG="result/{sample}/trimmomatic/adaptor_read_trimm_info"
    conda:
        "config/envs/Java_related_mac.yaml"
    shell:
        "trimmomatic PE -threads {threads} -phred33 {input.trimmed1_R1} {input.trimmed1_R2} {output.adapttrim_R1} {output.adapttrimNP_R1} {output.adapttrim_R2} {output.adapttrimNP_R2} ILLUMINACLIP:{params.adapters}:2:30:10 2> {output.adapttrimLOG}"

rule primer_trim:
    input:
        adapter_trimmed_R1="result/{sample}/trimmomatic/r1.no_adaptors.fastq.gz",
        adapter_trimmed_R2="result/{sample}/trimmomatic/r2.no_adaptors.fastq.gz"
    params:
        primers=config["primers"]
    threads: workflow.cores * 0.75
    output:
        primer_trimmed_R1="result/{sample}/trimmomatic/r1.no_adaptors_no_primers.fastq.gz",
        primer_trimmedNP_R1="result/{sample}/trimmomatic/r1.no_adaptors_no_primers.not-paired.fastq.gz",
        primer_trimmed_R2="result/{sample}/trimmomatic/r2.no_adaptors_no_primers.fastq.gz",
        primer_trimmedNP_R2="result/{sample}/trimmomatic/r2.no_adaptors_no_primers.not-paired.fastq.gz",
        primer_trimmedLOG="result/{sample}/trimmomatic/primer_read_trimm_info"
    conda:
        "config/envs/Java_related_mac.yaml"
    shell:
        "trimmomatic PE -threads {threads} -phred33 {input.adapter_trimmed_R1} {input.adapter_trimmed_R2} {output.primer_trimmed_R1} {output.primer_trimmedNP_R1} {output.primer_trimmed_R2} {output.primer_trimmedNP_R2} ILLUMINACLIP:{params.primers}:2:30:10 2> {output.primer_trimmedLOG}"

rule second_length_trim:
    input:
        primer_trimmed_R1="result/{sample}/trimmomatic/r1.no_adaptors_no_primers.fastq.gz",
        primer_trimmed_R2="result/{sample}/trimmomatic/r2.no_adaptors_no_primers.fastq.gz"
    threads: workflow.cores * 0.75
    params:
        sw=config["QC_sliding_window"],
        phred=config["QC_phred_score"],
        minlen=config["QC_min_len"]
    output:
        end_trimmed_R1="result/{sample}/trimmomatic/r1.no_adaptors_no_primers_trimmed.fastq.gz",
        end_trimmedNP_R1="result/{sample}/trimmomatic/r1.no_adaptors_no_primers_trimmed_not-paired.fastq.gz",
        end_trimmed_R2="result/{sample}/trimmomatic/r2.no_adaptors_no_primers_trimmed.fastq.gz",
        end_trimmedNP_R2="result/{sample}/trimmomatic/r2.no_adaptors_no_primers_trimmed_not-paired.fastq.gz",
        end_trimmedLOG="result/{sample}/trimmomatic/no_adaptor_no_primer_quality_read_trimm_info"
    conda:
        "config/envs/Java_related_mac.yaml"
    shell:
        "trimmomatic PE -threads {threads} -phred33 {input.primer_trimmed_R1} {input.primer_trimmed_R2} {output.end_trimmed_R1} {output.end_trimmedNP_R1} {output.end_trimmed_R2} {output.end_trimmedNP_R2} SLIDINGWINDOW:{params.sw}:{params.phred} MINLEN:{params.minlen} 2> {output.end_trimmedLOG}"

rule align:
    input:
        QF_R1="result/{sample}/trimmomatic/r1.no_adaptors_no_primers_trimmed.fastq.gz",
        QF_R2="result/{sample}/trimmomatic/r2.no_adaptors_no_primers_trimmed.fastq.gz"
    threads: workflow.cores * 0.75
    params:
        ref=config["ref_genome"]
    conda:
        "config/envs/Java_related_mac.yaml"
    output:
        sam="result/{sample}/Mapping/{sample}.alignment.sam"
    shell:
        "bwa mem -t {threads} {params.ref} {input.QF_R1} {input.QF_R2} > {output.sam}"

rule sam_sort:
    input:
        sam1="result/{sample}/Mapping/{sample}.alignment.sam"
    threads: workflow.cores * 0.25
    conda: 
        "config/envs/Java_related_mac.yaml"
    output:
        bam="result/{sample}/Mapping/{sample}.alignment.bam",
        NDbam="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.bam"
    shell:
        """
        samtools sort -@ {threads} -T result/{wildcards.sample}/Mapping/temp_sort -o {output.bam} {input.sam1}
        samtools rmdup {output.bam} {output.NDbam}
        samtools index {output.NDbam}
        """

rule mapping:
    input:
        Fbam="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.bam",
    params:    
        Rmap=config["mapr"],
        Rumap=config["unmapr"]
    threads: workflow.cores * 0.50
    output:
        map_sam="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.reads.only.sam",
        unmap_sam="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.unmapped.reads.only.sam"
    shell:
        """
        config/envs/VariantBam/src/variant {input.Fbam} -r {params.Rmap} > {output.map_sam} -v
        config/envs/VariantBam/src/variant {input.Fbam} -r {params.Rumap} > {output.unmap_sam} -v
        """

rule map_stats:
    input:
        map="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.reads.only.sam",
        unmap="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.unmapped.reads.only.sam"
    threads: workflow.cores * 0.50
    conda:
        "config/envs/Java_related_mac.yaml"
    output:
        map_stats="result/{sample}/Mapping/{sample}.alignment.stats.tab"
    shell:
        """
        countmap="$(samtools view -c {input.map})"
        countunmap="$(samtools view -c {input.unmap})"
        tot=$((countmap+countunmap))
        percmap=$(bc <<< "scale=4;$countmap/$tot*100")
        percunmap=$(bc <<< "scale=4;$countunmap/$tot*100")
        echo {wildcards.sample} > {output.map_stats}
        echo "Total read count= $tot , 100 %" >> {output.map_stats}
        echo "Mapped read count= $countmap , $percmap %" >> {output.map_stats}
        echo "Unmapped read count= $countunmap , $percunmap %" >> {output.map_stats}
        """

rule sam_sort2:
    input:
        map="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.reads.only.sam",
        unmap="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.unmapped.reads.only.sam"
    threads: workflow.cores * 0.25
    conda:
        "config/envs/Java_related_mac.yaml"
    output:
        mapped_S="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.reads.only.sorted.bam",
        unmapped_S="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.unmapped.reads.only.sorted.bam"
    shell:
        """
        samtools sort -@ {threads} -T result/{wildcards.sample}/Mapping/temp_sort {input.map} -o {output.mapped_S}
        samtools index {output.mapped_S}
        samtools sort -@ {threads} -T result/{wildcards.sample}/Mapping/temp_sort {input.unmap} -o {output.unmapped_S}
        samtools index {output.unmapped_S}
        """

rule downtrim_DE:
    input:
        mapped_S="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.reads.only.sorted.bam"
    threads: workflow.cores * 0.25
    params: 
        up_tr=config["uptrim_threshold"]
    output:
        trimmed_1000="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.1000trimmed.sam"
    shell:
        "config/envs/VariantBam/src/variant {input.mapped_S} -m {params.up_tr} > {output.trimmed_1000} -v"

rule sam_sort_DE:
    input:
        trimmed_sam1="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.1000trimmed.sam"
    threads: workflow.cores * 0.25
    conda:
        "config/envs/Java_related_mac.yaml"
    output:
        trimmed_bam="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.1000trimmed.sorted.bam"
    shell:
        """
        samtools sort -@ {threads} -T result/{wildcards.sample}/Mapping/temp_sort -o {output.trimmed_bam} {input.trimmed_sam1}
        samtools index {output.trimmed_bam}
        """

rule coverage_stats:
    input:
        mapped_S="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.reads.only.sorted.bam" 
    threads: workflow.cores * 0.25
    conda:
        "config/envs/Java_related_mac.yaml"
    output:
        ave_cov="result/{sample}/Mapping/{sample}.average.coverage.tab"
    shell:
        """
        samtools depth -a {input.mapped_S} > result/{wildcards.sample}/Mapping/{wildcards.sample}.temp.cov 
        awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}}' result/{wildcards.sample}/Mapping/{wildcards.sample}.temp.cov >> {output.ave_cov}
        rm result/{wildcards.sample}/Mapping/{wildcards.sample}.temp.cov
        """

rule variant_call: 
    input:
        frags="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.reads.only.sorted.bam"
    threads: workflow.cores * 1
    params:
        ref=config["ref_genome"]
    conda:
        "config/envs/Java_related_mac.yaml"
    output:
        vcf="result/{sample}/variantcall/{sample}.vcf",
        changes="result/{sample}/variantcall/{sample}.changes"
    shell:
        "pilon --threads {threads} --genome {params.ref} --frags {input.frags} --changes --variant --outdir result/{wildcards.sample}/variantcall/ --output {wildcards.sample}"

rule variant_filter:
    input:
        vcf="result/{sample}/variantcall/{sample}.vcf"
    threads: workflow.cores * 0.5
    params:
        af=config["variant_freq"],
        dp=config["variant_depth"]
    conda:
        "config/envs/vcflib_mac.yaml"
    output:
        vcf_filt="result/{sample}/variantcall/{sample}.depth.filtered.vcf",
        final_vcf="result/{sample}/variantcall/{sample}.final.vcf"
    shell:
        """
        vcffilter -f "SVTYPE = DEL & ! IMPRECISE" --or -f "SVTYPE = INS & ! IMPRECISE" --or -f "AF > {params.af} & DP > {params.dp}" {input.vcf} > {output.vcf_filt}
        grep -v "<DUP>" {output.vcf_filt}  > {output.final_vcf}
        """

rule variant_formats:
    input:
        final_vcf="result/{sample}/variantcall/{sample}.final.vcf"
    conda:
        "config/envs/bcftools_mac.yaml"
    output:
        compressed="result/{sample}/variantcall/{sample}.consvariants.vcf.gz",
        indexed="result/{sample}/variantcall/{sample}.consvariants.vcf.gz.tbi"
    shell:
        """
        bgzip -c -i {input.final_vcf} > {output.compressed}
        tabix -p vcf {output.compressed} -f > {output.indexed}
        """

rule minority_variants:
    input:
        raw_vcf="result/{sample}/variantcall/{sample}.vcf"
    conda:
        "config/envs/vcflib_mac.yaml"
    output:
        highAF="result/{sample}/variantcall/{sample}.afM0.vcf",
        minority_vcf="result/{sample}/variantcall/{sample}.minority_alleles.vcf",
        minority_techrep="result/{sample}/variantcall/{sample}.minority_alleles_report.tech"
    shell:
        """
        vcffilter -f "AF > 0" {input.raw_vcf} > {output.highAF}
        config/dependencies/minority_revealer_1.sh {output.highAF} result/{wildcards.sample}/variantcall/{wildcards.sample}
        """

rule annotation:
    input:
        primary_vcf="result/{sample}/variantcall/{sample}.final.vcf",
        minority_vcf="result/{sample}/variantcall/{sample}.minority_alleles.vcf"
    conda: 
        "config/envs/snpeff_mac.yaml"
    output:
        annot_primary_vcf="result/{sample}/variantcall/{sample}.annotated.variants.vcf",
        annot_minority_vcf="result/{sample}/variantcall/{sample}.annotated.minorityvariants.vcf"
    shell:
        """
        snpEff -v NC_045512.2 {input.primary_vcf} > {output.annot_primary_vcf}
        snpEff -v NC_045512.2 {input.minority_vcf} > {output.annot_minority_vcf}
        """
#note, not included the passage of single coverage plot drafting, in Pipe_10.3 line 306-307
rule maskfinding:
    input:
        Mbam="result/{sample}/Mapping/{sample}.alignment.removed.duplicates.mapped.reads.only.sorted.bam",
        Fvcf="result/{sample}/variantcall/{sample}.final.vcf"
    output:
        lowcovreg="result/{sample}/variantcall/{sample}.lessthan50.bed",
        mask="result/{sample}/variantcall/{sample}.finalmask.bed"
    params: 
        dp=config["variant_depth"]
    conda:
        "config/envs/bedtools_mac.yaml"
    shell:
        """
        bedtools genomecov -bga -ibam {input.Mbam} > temp_{wildcards.sample}_out.bed
        awk '$4<{params.dp}' temp_{wildcards.sample}_out.bed > {output.lowcovreg}
        bedtools subtract -a {output.lowcovreg} -b {input.Fvcf} > {output.mask}
        rm temp_{wildcards.sample}_out.bed
        """

rule consensuscall:
    input:
        consensus_variants="result/{sample}/variantcall/{sample}.consvariants.vcf.gz",
        mask="result/{sample}/variantcall/{sample}.finalmask.bed"
    output:
        pre_consensus="result/{sample}/variantcall/{sample}.old.consensus.fasta"
    conda:
        "config/envs/bcftools_mac.yaml"
    params:
        ref=config["ref_genome"]
    shell:
        "bcftools consensus {input.consensus_variants} -f {params.ref} -m {input.mask} -o {output.pre_consensus}"

rule Nstats:
    input: 
        cons="result/{sample}/variantcall/{sample}.old.consensus.fasta"
    output:
        nstats="result/{sample}/variantcall/{sample}.Nstats.tab"
    conda:
        "config/envs/seqtk_mac.yaml"
    shell:
        """
        echo {wildcards.sample} > {output.nstats}
        seqtk comp {input.cons} > temp{wildcards.sample}_seqtk.temp 
        awk '{{Tot+=$2; N+=$9}} END {{ print "Total length = ",Tot; print "Ns = ",N; print "Ns Perc = ", N/Tot*100,"%"}}' temp{wildcards.sample}_seqtk.temp >> {output.nstats}
        rm temp{wildcards.sample}_seqtk.temp
        """

rule tagging:
    input:
        cons="result/{sample}/variantcall/{sample}.old.consensus.fasta"
    output:
        clas="result/{sample}/variantcall/{sample}.classification.tab"
        #finalcons="result/{sample}/variantcall/{sample}.consensus.{tag}.fasta"
    params:
        nt=config["n_threshold"]
    conda:
        "config/envs/seqtk_mac.yaml"
    shell:
        """
        config/dependencies/tagger.sh {wildcards.sample} {params.nt}
        """
