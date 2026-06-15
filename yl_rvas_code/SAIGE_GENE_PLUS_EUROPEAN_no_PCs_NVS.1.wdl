# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <noah_fields@dfci.harvard.edu>
# # Distributed under the terms of the GNU GPL v2.0

version 1.0

workflow SAIGE_GENE_EUROPEAN {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    File rare_variants = "gs://fc-73a5dbdd-51d6-4c28-a211-234c026bdb26/step_9_run_vep/STEP_9_RUN_VEP/yl_rare_variants.001.tsv"
    File ld_pruned_vcf = "gs://fc-73a5dbdd-51d6-4c28-a211-234c026bdb26/submissions/bce6e8be-1b76-4a41-a9a6-a534a85c4e9d/STEP_11_GENETIC_RELATEDNESS/f0fe51fb-b26c-4598-a69d-b6f895c374f1/call-ConcatVcfsGenome/young_lung_1000G_snps.vcf.gz"
    File covariates_data = "gs://fc-73a5dbdd-51d6-4c28-a211-234c026bdb26/saige_gene_plus/inputs/covariate_data.tsv"
    File european_pcs = "gs://fc-73a5dbdd-51d6-4c28-a211-234c026bdb26/saige_gene_plus/inputs/PLINK_ANCESTRY.EUR.ONLY.top10PC.020725.eigenvec"
    String output_file = "gs://fc-73a5dbdd-51d6-4c28-a211-234c026bdb26/saige_gene_plus/outputs/"
    File groupfile_001 = "gs://fc-73a5dbdd-51d6-4c28-a211-234c026bdb26/step_9_run_vep/STEP_9_RUN_VEP/v2/ufc_singleton.groupfile"
    File groupfile_0001 = "gs://fc-73a5dbdd-51d6-4c28-a211-234c026bdb26/step_9_run_vep/STEP_9_RUN_VEP/ufc_0001.groupfile"
  }
  scatter (i in range(length(vcfs))){
    call process_vcf_part1 {
      input:
        vcf = vcfs[i],
        rare_variants = rare_variants
    }
  }
  call ConcatVcfs {
    input:
      vcfs = process_vcf_part1.out_vcf,
      vcf_idxs = process_vcf_part1.out_vcf_idx
  }
  call get_covariates {
    input:
      data = covariates_data,
      european_pcs = european_pcs
  }

  call plink {
    input: 
      vcf = ld_pruned_vcf,
      output_prefix = "ld_pruned"
  }
  
  call saige_gene_step0 {
    input:
      bed = plink.bed,
      bim = plink.bim,
      fam = plink.fam
  }

  call saige_gene_step1 {
    input:
      bed = plink.bed,
      bim = plink.bim,
      fam = plink.fam,
      sample_data = get_covariates.covariate_data,
      sparseGRMFile = saige_gene_step0.sparseGRMFile,
      sparseGRMSampleIDFile = saige_gene_step0.sparseGRMSampleIDFile
  }
  
  

  call saige_gene_step2 as saige_gene_step2_001{
    input:
      pathogenic_criteria = "yl",
      vcf = ConcatVcfs.merged_vcf,
      vcf_idx = ConcatVcfs.merged_vcf_idx,
      sampleFile = get_covariates.subjects_list,
      rda = saige_gene_step1.rda,
      group_file = groupfile_001,
      varianceRatio = saige_gene_step1.varianceRatio,
      sparseGRMFile = saige_gene_step0.sparseGRMFile,
      sparseGRMSampleIDFile = saige_gene_step0.sparseGRMSampleIDFile
  }
  call saige_gene_step2_beta as saige_gene_step2_beta_001{
    input:
      pathogenic_criteria = "yl",
      vcf = ConcatVcfs.merged_vcf,
      vcf_idx = ConcatVcfs.merged_vcf_idx,
      sampleFile = get_covariates.subjects_list,
      rda = saige_gene_step1.rda,
      group_file = groupfile_001,
      varianceRatio = saige_gene_step1.varianceRatio,
      sparseGRMFile = saige_gene_step0.sparseGRMFile,
      sparseGRMSampleIDFile = saige_gene_step0.sparseGRMSampleIDFile
  }
  call saige_gene_step2 as saige_gene_step2_0001{
    input:
      pathogenic_criteria = "yl",
      vcf = ConcatVcfs.merged_vcf,
      vcf_idx = ConcatVcfs.merged_vcf_idx,
      sampleFile = get_covariates.subjects_list,
      rda = saige_gene_step1.rda,
      group_file = groupfile_0001,
      varianceRatio = saige_gene_step1.varianceRatio,
      sparseGRMFile = saige_gene_step0.sparseGRMFile,
      sparseGRMSampleIDFile = saige_gene_step0.sparseGRMSampleIDFile
  }
  call saige_gene_step2_beta as saige_gene_step2_beta_0001{
    input:
      pathogenic_criteria = "yl",
      vcf = ConcatVcfs.merged_vcf,
      vcf_idx = ConcatVcfs.merged_vcf_idx,
      sampleFile = get_covariates.subjects_list,
      rda = saige_gene_step1.rda,
      group_file = groupfile_0001,
      varianceRatio = saige_gene_step1.varianceRatio,
      sparseGRMFile = saige_gene_step0.sparseGRMFile,
      sparseGRMSampleIDFile = saige_gene_step0.sparseGRMSampleIDFile
  }

  call sortSAIGE_Output {
    input:
      saige_output_001 = saige_gene_step2_001.out1,
      saige_output_beta_001 = saige_gene_step2_beta_001.out1,
      saige_output_0001 = saige_gene_step2_0001.out1,
      saige_output_beta_0001 = saige_gene_step2_beta_0001.out1,
      output_name = "yl"
  }

  call RecalibrateSaigeStats {
    input:
      raw_stats_tsv = sortSAIGE_Output.out1,
      output_prefix = "yl"
  }
}


task bcftools_norm {
  input{
    File vcf
    Int bcftools_norm_disk
  }
  command <<<
  set -x pipefail
  bcftools annotate -x FORMAT/AD -O z -o tmp1.vcf.gz ~{vcf}
  rm ~{vcf}
  bcftools annotate -x INFO/dbNSFP_REVEL_score tmp1.vcf.gz -Oz -o tmp2.vcf.gz
  rm tmp1.vcf.gz
  bcftools view -m2 -M2 --include 'FILTER=="ExcessHet" || FILTER=="PASS"' tmp2.vcf.gz -Oz -o tmp3.vcf.gz
  rm tmp2.vcf.gz
  bcftools norm -m - tmp3.vcf.gz -Oz -o tmp4.vcf.gz
  rm tmp2.vcf.gz
  bcftools view -e 'ALT="*"' tmp4.vcf.gz -Oz -o output.vcf.gz
  rm tmp4.vcf.gz
  >>>
  output{
    File out1 = "output.vcf.gz"
  }
  runtime{
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
    disks: "local-disk ~{bcftools_norm_disk} HDD"
  }
}

task process_vcf_part1 {
  input {
    File vcf
    File rare_variants
  }
  Int default_disk_gb = ceil(size(vcf,"GB") * 2) + 10
  command <<<
    set -eu -o pipefail
    bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' ~{vcf} -Oz -o tmp.vcf.gz
    bcftools view --include ID==@~{rare_variants} tmp.vcf.gz -O z -o out.vcf.gz
    bcftools index -t out.vcf.gz
    
  >>>
  output {
    File out_vcf = "out.vcf.gz"
    File out_vcf_idx = "out.vcf.gz.tbi"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    disks: "local-disk ~{default_disk_gb} HDD"
    preemptible: 3
  }
}

task process_vcf_part2 {
  input {
    File vcf
    File subjects_list
  }
  Int default_disk_gb = ceil(size(vcf,"GB") * 2) + 10
  command <<<
    set -eu -o pipefail
    echo "Check 1"
    bcftools view -S ~{subjects_list} --force-samples ~{vcf} -O z -o tmp.vcf.gz
    bcftools view -e 'ALT="*"' tmp.vcf.gz -Oz -o output.vcf.gz
    echo "Check 2"
    bcftools index --csi output.vcf.gz
    echo "Check 3"
  >>>
  output {
    File output_vcf = "output.vcf.gz"
    File output_vcf_idx = "output.vcf.gz.csi"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    disks: "local-disk ~{default_disk_gb} HDD"
    preemptible: 3
  }
}

task make_group_file_part1 {
  input {
    File vcf
    File autosomal_gene_file
    Int make_group_file_part1_memory = 4
  }
  command <<<
  set -euxo pipefail
  
  echo "Checkpoint 1"
  bcftools +split-vep ~{vcf} -f '%SYMBOL\t%CHROM:%POS:%REF:%ALT\t%IMPACT\t%Consequence\n' -d > tmp.txt
  grep -Eiv 'MODIFIER|LOW|nmd_transcript_variant' tmp.txt > tmp1.txt
  rm tmp.txt
  cut -f1,2 tmp1.txt > tmp2.txt
  rm tmp1.txt
  grep -v '^$' tmp2.txt > tmp3.txt
  rm tmp2.txt
  sort -u tmp3.txt > vep.tmp.txt
  rm tmp3.txt
  awk 'NR==FNR {genes[$1]; next} $1 in genes' ~{autosomal_gene_file} vep.tmp.txt > vep.txt
  >>>
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
    memory: "~{make_group_file_part1_memory} GB"
    disks: "local-disk 20 HDD"
  }
  output {
    File vep_output = "vep.txt"
  }
}

task make_group_file_part2 {
  input {
    File split_vep_output
  }
  command <<<
  set -eu -o pipefail
  python3 <<CODE
  from collections import defaultdict

  # Read the vep.txt file
  vep_file = "~{split_vep_output}"
  output_file = "groupfile.txt"

  # Dictionary to store SYMBOL and their associated variants
  symbol_to_variants = defaultdict(list)

  # Parse vep.txt
  with open(vep_file, "r") as file:
      for line in file:
          line = line.strip()
          if line:  # Skip empty lines
              symbol, variant = line.split('\t')
              symbol_to_variants[symbol].append(variant)

  # Write the output file
  with open(output_file, "w") as out_file:
      for symbol, variants in symbol_to_variants.items():
          # Line 1: SYMBOL var CHROM:%POS:%REF:%ALT (for all variants)
          var_line = f"{symbol} var " + " ".join(variants)

          out_file.write(var_line + "\n")
        
          # Line 2: SYMBOL anno Variant_of_Interest (repeated for the count of variants)
          anno_line = f"{symbol} anno " + " ".join(["Variant_of_Interest"] * len(variants))
          out_file.write(anno_line + "\n")
  CODE
  >>>
  output {
    File group_file = "groupfile.txt"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task get_covariates {
  input {
    File data
    File european_pcs
  }
  command <<<
  set -eu -o pipefail
  python3 <<CODE
  import pandas as pd
  from scipy.stats import zscore

  covariates_data = pd.read_csv("~{data}",sep='\t')
  covariates_data = covariates_data[(covariates_data["Cohort"] == "NVS") | (covariates_data["Cohort"] == "TOPMED")]
  covariates_data = covariates_data[covariates_data['histology'] != "Carcinoid"]
  
  covariates_data["cancer_status"] = covariates_data["case_control"].apply(lambda x: 1 if x == "case" else 0)
  #covariates_data["NVS_COHORT"] = covariates_data["Cohort"].apply(lambda x: 1 if x == "NVS" else 0)

  #Remove continental PC data
  #covariates_data = covariates_data.drop(columns=[f"PC{i}" for i in range(1, 11)])

  # Add in European PC data
  #european_pcs = pd.read_csv("~{european_pcs}",sep='\t')
  #covariates_data = pd.merge(covariates_data, european_pcs, on="sample_id")

  # Clear samples w/ NA in any relevant covariate
  covariates_data = covariates_data[
    covariates_data['Cohort'].notna() &
    covariates_data['sex_binary'].notna() &
    covariates_data['PC1_EUR'].notna() &
    covariates_data['PC2_EUR'].notna() &
    covariates_data['PC3_EUR'].notna() &
    covariates_data['PC4_EUR'].notna() &
    covariates_data['PC5_EUR'].notna()
  ]
  # Filter to those truly European
  covariates = [f"PC{i}" for i in range(1, 6)]
  covariates_data[covariates] = covariates_data[covariates].apply(zscore)

  covariates_data = covariates_data[
    (abs(covariates_data['PC1']) < 1) &
    (abs(covariates_data['PC2']) < 1)
  ]

  # Standard normalize PCs
  covariates = [f"PC{i}_EUR" for i in range(1, 6)]
  #covariates = [f"PC{i}" for i in range(1, 6)]
  covariates_data[covariates] = covariates_data[covariates].apply(zscore)
  #covariates_data['age_dx'] = zscore(covariates_data['age_dx'])
  
  # Write to file
  covariates_data.to_csv("nvs_covariates.tsv",sep='\t', index=False)
  covariates_data['sample_id'].to_csv("nvs_subjects.list", index=False, header=False)
  print("File written to nvs_covariates.tsv")
  CODE
  >>>
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
  output {
    File covariate_data = "nvs_covariates.tsv"
    File subjects_list = "nvs_subjects.list"
  }
}

task saige_gene_step0 {
  input {
    File bed
    File bim
    File fam
  }
  command <<<
  set -eu -o pipefail

  # Ensure all files are in the same directory
  cp ~{bed} ufc.bed
  cp ~{bim} ufc.bim
  cp ~{fam} ufc.fam

  if ! createSparseGRM.R       \
       --plinkFile=./ufc \
       --nThreads=4  \
       --outputPrefix=sparseGRM       \
       --numRandomMarkerforSparseKin=2000      \
       --relatednessCutoff=0.125; then
    echo "Error: Sparse GRM creation failed" >&2
    exit 1
  fi

  echo "Files created by SparseGRM.R"
  ls *mtx*
  >>>
  runtime {
    docker: "wzhou88/saige:1.1.9"
    preemptible: 3
    memory: "8 GB"
  }
  output {
    File sparseGRMFile = "sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
    File sparseGRMSampleIDFile = "sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
  }
}

task saige_gene_step1 {
  input {
    File bed
    File bim
    File fam
    File sample_data
    File sparseGRMFile
    File sparseGRMSampleIDFile
  }
  command <<<
  set -eu -o pipefail

  # Ensure all files are in the same directory
  cp ~{bed} young_lung.bed
  cp ~{bim} young_lung.bim
  cp ~{fam} young_lung.fam

  echo "Checkpoint 1"

  #run step 1
  step1_fitNULLGLMM.R \
      --plinkFile=./young_lung \
      --sparseGRMFile=~{sparseGRMFile} \
      --sparseGRMSampleIDFile=~{sparseGRMSampleIDFile} \
      --useSparseGRMtoFitNULL=TRUE \
      --isCateVarianceRatio=FALSE \
      --phenoFile=~{sample_data} \
      --phenoCol=cancer_status \
      --traitType=binary \
      --covarColList=sex_binary \
      --sampleIDColinphenoFile=sample_id \
      --traitType=binary        \
      --outputPrefix=./test_run

  echo "Checkpoint 5"
  ls test_run*
  echo "Checkpoint 6"
  >>>
  output {
    File rda = "test_run.rda"
    File varianceRatio = "test_run.varianceRatio.txt"
  }
  runtime {
    docker: "wzhou88/saige:1.3.0"
    preemptible: 3
    memory: "8 GB"
  }
}


task saige_gene_step2 {
  input {
    File vcf
    File vcf_idx
    String pathogenic_criteria
    File sampleFile
    File rda
    File group_file
    File varianceRatio
    File sparseGRMFile
    File sparseGRMSampleIDFile
  }

  command <<<

  if ! step2_SPAtests.R \
    --vcfFile=~{vcf} \
    --vcfFileIndex=~{vcf_idx} \
    --vcfField=GT \
    --SAIGEOutputFile=saige_gene_output.~{pathogenic_criteria}.txt \
    --LOCO=FALSE \
    --minMAF=0 \
    --minMAC=0.5 \
    --sampleFile=~{sampleFile} \
    --GMMATmodelFile=~{rda} \
    --varianceRatioFile=~{varianceRatio} \
    --sparseGRMFile=~{sparseGRMFile} \
    --sparseGRMSampleIDFile=~{sparseGRMSampleIDFile} \
    --groupFile=~{group_file} \
    --annotation_in_groupTest=T1,T1:T2,T1:T2:T3,T1:T2:T3:T4,T1:T2:T3:T4:T5,T2,T3,T4,T5 \
    --maxMAF_in_groupTest=0.1; then

    echo "Error: Step2 failed" >&2
    exit 0
  fi

  echo "Locations of file"
  ls saige_gene_output*
  >>>
  output {
    File out1 = "saige_gene_output.~{pathogenic_criteria}.txt"
    File out3 = "saige_gene_output.~{pathogenic_criteria}.txt.singleAssoc.txt"
  }
  runtime{
    docker:"wzhou88/saige:1.3.0"
    memory: "8 GB"
  }
}

task saige_gene_step2_beta {
  input {
    File vcf
    File vcf_idx
    String pathogenic_criteria
    File sampleFile
    File rda
    File group_file
    File varianceRatio
    File sparseGRMFile
    File sparseGRMSampleIDFile
  }

  command <<<

  if ! step2_SPAtests.R \
    --vcfFile=~{vcf} \
    --vcfFileIndex=~{vcf_idx} \
    --vcfField=GT \
    --SAIGEOutputFile=saige_gene_output.~{pathogenic_criteria}.txt \
    --LOCO=FALSE \
    --minMAF=0 \
    --minMAC=0.5 \
    --sampleFile=~{sampleFile} \
    --GMMATmodelFile=~{rda} \
    --varianceRatioFile=~{varianceRatio} \
    --sparseGRMFile=~{sparseGRMFile} \
    --sparseGRMSampleIDFile=~{sparseGRMSampleIDFile} \
    --groupFile=~{group_file} \
    --is_no_weight_in_groupTest=TRUE \
    --annotation_in_groupTest=T1,T1:T2,T1:T2:T3,T1:T2:T3:T4,T1:T2:T3:T4:T5,T2,T3,T4,T5 \
    --maxMAF_in_groupTest=0.1; then

    echo "Error: Step2 failed" >&2
    exit 0
  fi

  echo "Locations of file"
  ls saige_gene_output*
  >>>
  output {
    File out1 = "saige_gene_output.~{pathogenic_criteria}.txt"
    File out3 = "saige_gene_output.~{pathogenic_criteria}.txt.singleAssoc.txt"
  }
  runtime{
    docker:"wzhou88/saige:1.3.0"
    memory: "8 GB"
  }
}

task sortSAIGE_Output {
  input {
    File saige_output_001
    File saige_output_beta_001
    File saige_output_0001
    File saige_output_beta_0001
    String output_name
  }
  String out_filename = output_name + ".saige.tsv"
  command <<<
  set -x
    # 1% AF
    head -n 1 ~{saige_output_001} > header.txt
    tail -n +2 ~{saige_output_001} | grep -v Cauchy > tmp_001.tsv
    tail -n +2 ~{saige_output_beta_001} | grep -v Cauchy > tmp.beta_001.tsv
    
    paste <(cut -f1-6 tmp_001.tsv) <(cut -f7,8 tmp.beta_001.tsv) <(cut -f9- tmp_001.tsv) > tmp.updated_001.tsv
    awk 'BEGIN{OFS="\t"} NR==1 {print; next} { $3 = 0.0001; print }' tmp.updated_001.tsv > tmp.updated_again_001.tsv
    sort -k4,4g tmp.updated_again_001.tsv > tmp.final_001.tsv

    # 0.1% AF
    tail -n +2 ~{saige_output_0001} | grep -v Cauchy > tmp_0001.tsv
    tail -n +2 ~{saige_output_beta_0001} | grep -v Cauchy > tmp.beta_0001.tsv

    paste <(cut -f1-6 tmp_0001.tsv) <(cut -f7,8 tmp.beta_0001.tsv) <(cut -f9- tmp_0001.tsv) > tmp.updated_0001.tsv
    awk 'BEGIN{OFS="\t"} NR==1 {print; next} { $3 = 0.001; print }' tmp.updated_0001.tsv > tmp.updated_again_0001.tsv
    sort -k4,4g tmp.updated_again_0001.tsv > tmp.final_0001.tsv

    # Put it all together
    cat header.txt > ~{out_filename}
    cat tmp.final_001.tsv tmp.final_0001.tsv >> ~{out_filename}
    #cat tmp.final_001.tsv tmp.final_0001.tsv | awk -F'\t' '$10 > 1' >> ~{out_filename}
  >>>
  runtime{
    docker:"ubuntu:latest"
    preemptible: 3
  }
  output{
    File out1 = "~{out_filename}"
  }
}

task plink {
  input {
    File vcf
    String output_prefix
  }

  Int default_disk_gb = ceil(2 * size(vcf,"GB")) + 8

  command <<<
  set -eu -o pipefail

  # Creating Plink Files from VCF; --double-id flag makes IID and FID the same value. Needs to be changed eventually
  plink --vcf ~{vcf} --make-bed --out ~{output_prefix} --double-id

  >>>
  output {
    File bed = "~{output_prefix}.bed"
    File bim = "~{output_prefix}.bim"
    File fam = "~{output_prefix}.fam"
  }
  runtime {
    docker: "elixircloud/plink:1.9-20210614"
    disks: "local-disk ~{default_disk_gb} HDD"
    preemptible: 3
  }
}

task RecalibrateSaigeStats {
  input {
    File raw_stats_tsv
    String output_prefix

    Int disk_gb = 20
    Float mem_gb = 7.5
    Int n_cpu = 4
    String docker = "vanallenlab/g2c_analysis:4c4511c"
  }

  String out_tsv = output_prefix + ".saige.recalibrated.stats.tsv"

  command <<<
    set -eu -o pipefail

    /opt/pancan_germline_wgs/scripts/association/recalibrate_saige_gene.R \
      -i "~{raw_stats_tsv}" \
      -o "~{out_tsv}"
    gzip -f "~{out_tsv}"
  >>>

  output {
    File out1 = "~{out_tsv}.gz"
  }

  runtime {
    docker: docker
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
    maxRetries: 1
  }
}

task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String callset_name = "tmp"

    String bcftools_concat_options = ""

    Float mem_gb = 3.5
    Int cpu_cores = 2
    Int? disk_gb

    String bcftools_docker = "vanallenlab/bcftools"
  }

  String out_filename = callset_name + ".vcf.gz"

  Int default_disk_gb = ceil(2.5 * size(vcfs, "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools concat \
      ~{bcftools_concat_options} \
      --file-list ~{write_lines(vcfs)} \
      -O z \
      -o ~{out_filename} \
      --threads ~{cpu_cores}

    bcftools index --csi ~{out_filename}
  >>>

  output {
    File merged_vcf = "~{out_filename}"
    File merged_vcf_idx = "~{out_filename}.csi"
  }

  runtime {
    docker: bcftools_docker
    memory: mem_gb + " GB"
    cpu: cpu_cores
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
  }
}