library(cn.mops)

# Input
window_length = 1000
base_dir = "/lustre/scratch118/infgen/team81/jl11/cnvs/S_pneumoniae_new"
bam_files = "bam_file_list.txt"

samples <- read.table(file=paste(base_dir, bam_files, sep="/"),
                      header=T, sep="\t", stringsAsFactors = F)

# Run algorithm
bamDataRanges <- getReadCountsFromBAM(samples$bam_path, 
                    sampleNames=samples$sample_name,mode="paired",WL=window_length)
cnv_result <- haplocn.mops(bamDataRanges)
int_cnv_result <- calcIntegerCopyNumbers(cnv_result)

# Save R data frames for later use
saveRDS(cnv_result,file="cnv_result.RData")
saveRDS(int_cnv_result,file="int_cnv_result.RData")

# Process into GWAS compatible data
segm <- as.data.frame(segmentation(int_cnv_result))
CNVs <- as.data.frame(cnvs(int_cnv_result))
CNVRegions <- as.data.frame(cnvr(int_cnv_result))

# Write output
write.csv(segm, "segmentation.csv")
write.csv(CNVs,file="cnvs.csv")
write.csv(CNVRegions,"cnvr.csv")


