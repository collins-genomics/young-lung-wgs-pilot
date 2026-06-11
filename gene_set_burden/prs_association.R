# Process data and report the frequency of doubleton carriers and MAF<5% (in our population) carriers
# Correlate with clinical variables (somatic driver genes and related info)

library(data.table)

setwd("~/Desktop/LungData")
levels <- c("CV")
tfams <- lapply(levels, function(lvl) read.table(paste0("PilotData/PRS_variants/", lvl, ".tfam"), as.is=TRUE))
names(tfams) <- levels
tpeds <- lapply(levels, function(lvl) read.table(paste0("PilotData/PRS_variants/", lvl, ".tped"), as.is=TRUE))
names(tpeds) <- levels

sampData <- fread("PilotData/GenesetData/Samples/all_covariates.tsv")[!(sample_id %like% "TPMCCDG")]
# sampData <- fread("~/Desktop/LungData/PilotData/Misc/YL_NVS_metadata_metrics_EUR_LUAD_only_45_under.tsv")

catIds <- list(
  "YL" = sampData$sample_id[!(sampData$sample_id %like% "NSLC")],
  "NVS" = sampData$sample_id[sampData$sample_id %like% "NSLC"],
  "YL_NVS" = sampData$sample_id
)

resTbl <- data.table()

for (cat in c("YL")) {
  for (lvl in levels[1]) {
    mat = tpeds[[lvl]]
    fam = tfams[[lvl]]
    mat[, 2] = paste("chr", mat[, 1], "_", mat[, 4], "_", mat[, 5], "_", seq_len(nrow(mat)), sep = "")
    mat.calls = mat[,5:ncol(mat)]
    row.names(mat.calls) = mat[,2]
    mat.calls[ mat.calls == "0" ] = NA
    maxval = apply(mat.calls,1,function(x) names(which.max(table(x))))
    miss.drop = is.na(maxval)
    mat.calls = mat.calls[!miss.drop,]
    maxval = maxval[!miss.drop]
    
    print("binarizing data")
    
    mat.calls.bin = matrix(NA,nrow=nrow(mat.calls),ncol=ncol(mat.calls))
    for ( i in 1:ncol(mat.calls) ) {
      mat.calls.bin[ !is.na(mat.calls[,i]) & mat.calls[,i] == maxval , i ] = 1
      mat.calls.bin[ !is.na(mat.calls[,i]) & mat.calls[,i] != maxval , i ] = 0
    }
    mat.calls.sum = mat.calls.bin[ , rep( c(T,F) , nrow(fam)) , drop=F ]
    mat.calls.sum = mat.calls.sum + mat.calls.bin[ , rep( c(F,T) , nrow(fam)) , drop=F ]
    
    print("joining with sample data")
    
    # keep unique samples and link with case/control status
    uni = sampData[sample_id %in% catIds[[cat]]]
    m = match(fam[,2],uni$sample_id)
    uni = uni[m,]
    
    keep = !is.na( m )
    fam = fam[keep,]
    mat.calls.sum = mat.calls.sum[,keep,drop=F]
    uni = uni[keep,]
    
    print("getting MAC")
    # count minor allele frequency and count
    maf = 1 - apply(mat.calls.sum,1,mean,na.rm=T)/2
    mac = apply( 2 - mat.calls.sum,1,sum,na.rm=T)
    
    # count MAC <= 2 carriers
    doubleton_carriers = apply( 2 - mat.calls.sum[ mac <= 2 , , drop=F ] , 2 , sum ,na.rm=T )
    
    print("running test")
    #for continuous burden
    burden = doubleton_carriers
    tsv=uni
    reg_prs = summary(glm(burden ~ tsv$adeno_Byun_PRS_PV7_scaledscore + tsv$sex + tsv$PC1_EUR + tsv$PC2_EUR + tsv$PC3_EUR + tsv$PC4_EUR + tsv$PC5_EUR,
                          family = "poisson"))
    result_pvalue_prs <- reg_prs$coef[2, 4]  # Extract p-value for PRS assoc
    result_beta_prs <- reg_prs$coef[2, 1]    # Extract beta for PRS assoc
    print(reg_prs)
    
    resTbl <- rbind(resTbl, data.table(category = paste0(cat, "_", lvl), 
                                       beta = reg_prs$coefficients[2,1],
                                       stderr = reg_prs$coefficients[2,2],
                                       pval = reg_prs$coefficients[2,4]))
  }
}

  