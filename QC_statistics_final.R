#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--level"), action="store", default="dataset", type='character',
              help="Level at which the QC filtering is peformed. By default, it will be performed at the whole dataset. You can chose to peform it at the level of a specific metadata variable (such as, lane/pool)."),
  make_option(c("--md_vars"), action="store", default=NULL, type='character',
              help="Metadata variables file (optional). By default, the QC statistics will be calculated at the whole dataset. You can choose to calcualte them by metadata variable."),
  make_option(c("--qc_mad"), action="store", default="/tools/wg1-qc_filtering/qc_mad_final.tab", type='character',
              help="QC metrics file. By default, the QC metrics-thresholds will be: lower nCount_RNA and upper percent.mt, from MAD 1 to 5."),
  make_option(c("--azimuth_pairing"), action="store", default="/tools/wg1-qc_filtering/azimuth_l1_l2.csv", type='character',
              help="Azimuth l1-l2 pairing file."),
  make_option(c("--downsampling"), action="store", default=NULL, type='character',
              help="Downsample file with metadata variable and number of cells (optional)."),
  make_option(c("--in_dir"), action="store", default=NULL, type='character',
              help="Input directory: WG2 output directory."),
  make_option(c("--out_dir"), action="store", default="QC_statistics", type='character',
              help="Output directory.")
)
opt = parse_args(OptionParser(option_list=option_list))

#################################### Set Variables ####################################
 # Input directory
in.dir <- paste0(opt$in_dir, '/')

# Loading functions
functions_fn <- '/tools/wg1-qc_filtering/scripts/QC_functions.R'
print(paste0('Loading functions from: ',functions_fn))
source(functions_fn)

# Read azimuth_l1_l2 pairing file
print(paste0('Reading azimuth_l1_l2 pairing file in: ',opt$azimuth_pairing))
pairs.fn <- opt$azimuth_pairing
pairs.df <- read.csv2(pairs.fn)
colnames(pairs.df) <- c('l1','l2')

# Read QC-bound-MAD file
print(paste0('Reading QC metrics and MAD combinations file in: ',opt$qc_mad))
qc_mad.fn <- paste0(opt$qc_mad)
qc.mad_df <- read.table(qc_mad.fn, header = T)
qc.mad.temp_list <- split(qc.mad_df, qc.mad_df$QC_metric)
qc.mad_list <- lapply(qc.mad.temp_list, function(i) qc_mad_format(i))

# Setting MAD combination order
mad_label.list <- lapply(qc.mad_list, function(x) paste0('MAD_',x$mad))
mad.order <- unlist(lapply(mad_label.list[[1]], function(i) paste0(i,':',mad_label.list[[2]])))

# Output directory
out.dir <- paste0(opt$out_dir, 'WG1_WG2_summary/')

print(paste0('Creating output directory in: ',out.dir))
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Print report
print(paste0('Considering ', nrow(qc.mad_df), ' QC metrics: ', paste(qc.mad_df$QC_metric, collapse = ', ')))
info <- lapply(names(qc.mad_list), function(qc){
  print(qc)
  print(paste0('MAD min: ', min(qc.mad_list[[qc]][["mad"]])))
  print(paste0('MAD max: ', max(qc.mad_list[[qc]][["mad"]])))
  cat('\n')
  return(NULL)
})

#################################### Summarize ####################################
# Read input data (metadata or seurat object in case metadata is not there)
md.fn <- paste0(in.dir, '/step4_reduce/metadata.reduced_data.RDS')
so.fn <- paste0(in.dir, '/step4_reduce/reduced_data.RDS')
if(file.exists(md.fn)){
  print(paste0('Reading metadata in: ', md.fn))
  system.time(so_metadata <- readRDS(md.fn))
}else if(file.exists(so.fn)){
  print(paste0('Reading seurat object in: ', so.fn))
  system.time(so <- readRDS(so.fn))
  print('Extracting metadata from the seurat object...')
  system.time(so_metadata <- so@meta.data)
  rm(so)
  gc()
}else{
  err_in <- paste0('Neither metadata object (metadata.reduced_data.RDS) nor seurat object (reduced_data.RDS) are found in the input directory: ',in.dir)
  stop(err_in) 
}

# Add Azimuth l1 classification to the metadata
print('Adding Azimuth l1 classification to the metadata...')
system.time(so_metadata <- add_azimuth_l1(so_md = so_metadata, pairing = pairs.df))

# Relabeling scPred to scPred_l2, and adding scPred_l1 classification to the metadata
print('Modifying scPred classification to the metadata...')
colnames(so_metadata)[which(colnames(so_metadata)=="scpred_prediction")]="scpred.l2.prediction"
so_metadata["scpred.l1.prediction"] = pairs.df$l1[match(so_metadata$scpred.l2.prediction,pairs.df$l2)]

# Calculate MADs and assigning Outlier/NotOutlier tags per cell barcode
print('Calculating MADs and assigning Outlier/NotOutlier tags per cell barcode...')
system.time(qc_tag.df <- qc_tag(so_md = so_metadata, 
                                qc_mad_list = qc.mad_list,
                                filter_level = opt$level))
# Save extended output
qc_tag_fn <- paste0(out.dir,'qc_tag.rds')
print(paste0('Saving extended output in: ', qc_tag_fn))
system.time(saveRDS(qc_tag.df, qc_tag_fn))
