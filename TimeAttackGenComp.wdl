workflow TimeAttackGenComp {

#Copyright 2021 Jamie K. Teer
#
#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    String output_name
    String output_dir

    # inputFile format depends on content. Tab-delimited, one row per sample
    #FASTQ: sample_name full_path_to_fastq_1 full_path_to_fastq_2
    #BAM:   sample_name full_path_to_bam full_path_to_bai
    #VCF:   sample_name full_path_to_vcf.gz full_path_to_tbi 
    File inputFile  
    Array[Array[File]] inputSamples = read_tsv(inputFile)
    Boolean runAlign
    Boolean runGenotype

    File target_bed
    String ref
    File ref_fasta

    String compare    
    String R_heatmap
    String plot_dist

    scatter (sample in inputSamples) {
        if (runAlign) {
            call align_pair {
                input: output_base=sample[0],
                    fq1=sample[1],
                    fq2=sample[2],
                    ref=ref
            }
        }

        if (runGenotype) {
            call call_vars {
                input: bam=select_first([align_pair.bam, sample[1]]),
                    bai=select_first([align_pair.bai, sample[2]]),
                    target_bed=target_bed,
                    ref_fasta=ref_fasta,
                    output_base=sample[0]
            }
        }
                
        call convert_vcf {
            input: vcf=select_first([call_vars.vcf, sample[1]]),
                tbi=select_first([call_vars.tbi, sample[2]]),
                target_bed=target_bed,
                output_name=sample[0]
        }

        call extract_af {
            input: vcf=select_first([call_vars.vcf, sample[1]]),
                tbi=select_first([call_vars.tbi, sample[2]]),
                target_bed=target_bed,
                output_name=sample[0]
        }

        call plot_af {
            input: af=extract_af.af,
                plot_dist=plot_dist,
                output_name=sample[0]
        }       

    }

    call compare_snvs {
        input: compare=compare,
            snvs=convert_vcf.simple_outs,
            target_bed=target_bed,
            output_name=output_name
    }

    call heatmap {
        input: comparison=compare_snvs.snv,
            R_heatmap=R_heatmap,
            output_name=output_name
    }

    call copy_results {
        input: output_dir=output_dir,
            comp_input=compare_snvs.snv,
            heatmap=heatmap.pdf,
            af_inputs= flatten([ extract_af.af, plot_af.dist, plot_af.chr ])
    }

    if (runGenotype) {
        call copy_vcf {
            input: output_dir=output_dir,
                vcf_inputs=select_all( flatten([call_vars.vcf, call_vars.tbi]) )
        }
    }
}

task align_pair {
    File fq1
    File fq2
    String ref
    String output_base
    String output_type = "wes"
    String snap_path = "/snap-2.0.1"

    command {
        ${snap_path}/snap-aligner \
            paired \
            ${ref} \
            ${fq1} \
            ${fq2} \
            -t 16 \
            -xf 2.0 \
            -so \
            -R '@RG\tID:${output_base}\tSM:${output_base}\tPL:ILLUMINA\tLB:${output_base}_${output_type}' \
            -o ${output_base}.bam
    }
    output {
        File bam = "${output_base}.bam"
        File bai = "${output_base}.bam.bai"
    }
    runtime {
        memory: "155G"
        pbs_cpu: "16"
        pbs_walltime: "6:00:00"
        docker: "timeattackgencomp_0.2"
    }
}

task call_vars {
    File bam
    File bai
    File target_bed
    File ref_fasta
    String output_base
    String bcf_path = "/bcftools-1.15.1"

    command {
        ${bcf_path}/bcftools mpileup \
            --output-type u \
            --no-BAQ \
            --max-depth 1000 \
            --fasta-ref ${ref_fasta} \
            --regions-file ${target_bed} \
            --annotate FORMAT/AD,FORMAT/DP \
            ${bam} \
            | ${bcf_path}/bcftools call \
                --multiallelic-caller \
                --output-type z \
                --keep-alts \
                --skip-variants indels \
                --output ${output_base}.vcf.gz \
        && ${bcf_path}/bcftools index -t ${output_base}.vcf.gz
    }
    output {
        File vcf = "${output_base}.vcf.gz"
        File tbi = "${output_base}.vcf.gz.tbi"
    }
    runtime {
        memory: "6G"
        pbs_walltime: "6:00:00"
        docker: "timeattackgencomp_0.2"
    }
}

task convert_vcf {
    File vcf
    File tbi
    File target_bed
    String output_name
    String bcf_path = "/bcftools-1.15.1"

    command <<<
        ${bcf_path}/bcftools query \
            -R ${target_bed} \
            -f '%CHROM\t%POS\t[%TGT]\n' \
            ${vcf} \
            | perl -a -F"\t" -nle 'unless ($F[2] eq "./.") {$F[2] =~ s/\///; $F[2] = join "", sort (split //, $F[2]);} print (join "\t", @F);' \
            > ${output_name}.smp.snv
    >>>
    output {
        File simple_outs = "${output_name}.smp.snv"
    }
    runtime {
        memory: "4G"
        pbs_walltime: "12:00:00"
        docker: "timeattackgencomp_0.2"
    }
}

task extract_af {
    File vcf
    File tbi
    File target_bed
    String output_name
    String bcf_path = "/bcftools-1.15.1"

    command <<<
        ${bcf_path}/bcftools query \
            -i 'FORMAT/DP>0 & QUAL >20' \
            -R ${target_bed} \
            -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%DP]\t[%AD]\n' ${vcf} \
            | perl -a -F"\t" -nle '@ad = split /,/, $F[6]; $af = sprintf("%0.3f", ($ad[1] / $F[5])); push @F, $af; print(join "\t", @F)' \
            > ${output_name}.af
    >>>
    output {
        File af = "${output_name}.af"
    }
    runtime {
        memory: "2G"
        pbs_walltime: "1:00:00"
        docker: "timeattackgencomp_0.2"
    }
}

task plot_af {
    File af
    String plot_dist
    String output_name

    command {
        R --vanilla < ${plot_dist} --args ${af} ${output_name}
    }
    output {
        File dist = "${output_name}.af.pdf"
        File chr = "${output_name}.af.chr.pdf"
    }
    runtime {
        memory: "2G"
        pbs_walltime: "1:00:00"
        docker: "timeattackgencomp_0.2"
    }
}

task compare_snvs {
    Array[File] snvs
    File target_bed
    String compare
    String output_name

    command {
        perl ${compare} \
            --raw_output ${output_name}.raw.out \
            --bed ${target_bed} \
            --missing './.' \
            ${sep=" " snvs} \
            > ${output_name}.snv.out.txt
    }
    output {
        File raw = "${output_name}.raw.out"
        File snv = "${output_name}.snv.out.txt"
    }
    runtime {
        memory: "50G"
        pbs_walltime: "180:00:00"
        docker: "timeattackgencomp_0.2"
    }
}

task heatmap {
    File comparison
    String output_name
    String R_heatmap

    command {
        R --vanilla < ${R_heatmap} --args ${comparison} ${output_name}
    }
    output {
        File pdf = "${output_name}.pdf"
    }
    runtime {
        memory: "5G"
        pbs_walltime: "2:00:00"
        docker: "timeattackgencomp_0.2"
    }
}

task copy_results {
    File comp_input
    File heatmap
    Array[File] af_inputs
    String output_dir

    command <<<
        mkdir -p ${output_dir}/af; \
        cp ${comp_input} ${output_dir}; \
        cp ${heatmap} ${output_dir}; \
        for i in ${sep=' ' af_inputs}; do cp $i ${output_dir}/af/; done;
    >>>
    runtime {
        memory: "1G"
        pbs_walltime: "2:00:00"
    }
}

task copy_vcf {
    String output_dir
    Array[File] vcf_inputs

    command <<<
        mkdir -p ${output_dir}/vcf; \
        for i in ${sep=' ' vcf_inputs}; do cp $i ${output_dir}/vcf/; done;
    >>>
    runtime {
        memory: "1G"
        pbs_walltime: "2:00:00"
    }
}
