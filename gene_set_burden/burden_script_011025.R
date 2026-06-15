# Process data and report the frequency of doubleton carriers; as binary or continuous 

#setwd("/Users/jaclynlopiccolo/Downloads/")

args = commandArgs(trailingOnly=TRUE)
pre = args[1]
uniq_samples_file = args[2]

mat = read.table(paste(pre,".tped",sep=''),as.is=T)

mat[, 2] = paste("chr", mat[, 1], "_", mat[, 4], "_", mat[, 5], "_", seq_len(nrow(mat)), sep = "")

fam = read.table(paste(pre,".tfam",sep=''),as.is=T)

mat.calls = mat[,5:ncol(mat)]
row.names(mat.calls) = mat[,2]

# get most common allele
mat.calls[ mat.calls == "0" ] = NA
maxval = apply(mat.calls,1,function(x) names(which.max(table(x))))

# remove fields with only missing values
miss.drop = is.na(maxval)

mat.calls = mat.calls[!miss.drop,]
maxval = maxval[!miss.drop]

# convert to binary
mat.calls.bin = matrix(NA,nrow=nrow(mat.calls),ncol=ncol(mat.calls))
for ( i in 1:ncol(mat.calls) ) {
  mat.calls.bin[ !is.na(mat.calls[,i]) & mat.calls[,i] == maxval , i ] = 1
  mat.calls.bin[ !is.na(mat.calls[,i]) & mat.calls[,i] != maxval , i ] = 0
}

# sum up the alternating columns to get per-individual counts
mat.calls.sum = mat.calls.bin[ , rep( c(T,F) , nrow(fam)) , drop=F ]
mat.calls.sum = mat.calls.sum + mat.calls.bin[ , rep( c(F,T) , nrow(fam)) , drop=F ]

# keep unique samples and link with case/control status
uni = read.table(uniq_samples_file,head=T,as.is=T)
m = match(fam[,2],uni$sample_ID)
uni = uni[m,]

keep = !is.na( m )
fam = fam[keep,]
mat.calls.sum = mat.calls.sum[,keep,drop=F]
uni = uni[keep,]

# count minor allele frequency and count
maf = 1 - apply(mat.calls.sum,1,mean,na.rm=T)/2
mac = apply( 2 - mat.calls.sum,1,sum,na.rm=T)

# count MAC <= 2 carriers
doubleton_carriers = apply( 2 - mat.calls.sum[ mac <= 2 , , drop=F ] , 2 , sum ,na.rm=T )

#to adjust for PCs, load them in
tsv <- read.csv("./YLvsTOPMED_metrics_PCs_SNV_SV_withEUR_PC_ASHK_ancestry_label_EURONLY_carcinoids_excluded_45under_cases.tsv", sep="\t", header=TRUE)

#if looking at counts only, remove !=0 and just use burden = doubleton_carriers; can also test singleton carriers
#for continuous burden use below
burden = doubleton_carriers
#burden = doubleton_carriers !=0
y = uni$case_control == "case"
#since y (case/control status) is binary variable, use logistic regression (glm) regardless of what type of variable burden is
reg = summary( glm( y ~ burden + tsv$sex + tsv$PC1_EUR + tsv$PC2_EUR + tsv$PC3_EUR + tsv$PC4_EUR + tsv$PC5_EUR , family="binomial" ) )

result_pvalue = reg$coef[2,4]
result_beta = reg$coef[2,1]
result_OR = exp(result_beta)

###(pre, result_pvalue, result_beta, result_OR, "\n")
cat(args[1], result_pvalue, result_beta, result_OR, "\n")

#mean(doubleton_carriers[uni$case_control == "case"]!=0)
#mean(doubleton_carriers[uni$case_control == "control"]!=0)

