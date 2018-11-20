library(data.table)
library(plyr)
library(ggplot2)
library(plotly)

{
  #Local directory of TCGA files. 
  tcga_tpm_files = list.files("/Volumes/ADATA/TCGA_isoform_auto", recursive = T, full.names = T)
  tcga_tpm_files = tcga_tpm_files[grep("rnaseqv2", tcga_tpm_files)]
}

{
  #Generation of TPM from 'scaled estimate'
  for (file in tcga_tpm_files){
    print(basename(file))
    subtype = strsplit(basename(file), "\\.")[[1]][1]
    if (!(file.exists(paste0(dirname(file), "/", subtype, "_TPM.csv")))){
      tcga_data = fread(file)
      tcga_data = tcga_data[-1, c(1, seq(3, ncol(tcga_data), 2)), with = F]
      tcga_data_genes = tcga_data$`Hybridization REF`
      tcga_data_patients = colnames(tcga_data)[-1]
      tcga_data = sapply(tcga_data[, -1], as.numeric)
      tcga_data = as.data.table(as.matrix(tcga_data) * matrix(1000000, nrow = nrow(tcga_data), ncol = ncol(tcga_data)))
      tcga_data$`Hybridization REF` = tcga_data_genes
      fwrite(tcga_data, paste0(dirname(file), "/", subtype, "_TPM.csv"))
    }
  }
}

{
  #Read TPM files.
  tcga_tpm_files = list.files("/Volumes/ADATA/TCGA_isoform_auto", recursive = T, full.names = T)
  tcga_tpm_files = tcga_tpm_files[grep("_TPM.csv", tcga_tpm_files)]
  tcga_subtypes_extended
}

{
  #Extract TRAIL and TRAILshort expression from TPM files, formatted for ggplot. 
  tcga_trail = list()
  tcga_trail_short = list()
  tcga_patient_id = list()
  tcga_subtype = list()
  tcga_subtype_extended = list()
  for (file in tcga_tpm_files){
    print(basename(file))
    tcga_data = fread(file)
    tcga_data = tcga_data[tcga_data$`Hybridization REF` %in% c("uc003fid.2", "uc003fie.2")]
    tcga_patient_id[[basename(file)]] = colnames(tcga_data)[-ncol(tcga_data)]
    tcga_subtype[[basename(file)]] = rep(gsub("_TPM.csv", "", basename(file)), ncol(tcga_data) - 1)
    tcga_subtype_extended[[basename(file)]] = rep(tcga_subtypes_extended[[gsub("_TPM.csv", "", basename(file))]], ncol(tcga_data) - 1)
    tcga_trail[[basename(file)]] = unlist(tcga_data[tcga_data$`Hybridization REF` == "uc003fid.2"][, -ncol(tcga_data), with = F])
    tcga_trail_short[[basename(file)]] = unlist(tcga_data[tcga_data$`Hybridization REF` == "uc003fie.2"][, -ncol(tcga_data), with = F])
  }
  tcga_trail = data.table(patient_id = unlist(tcga_patient_id), subtype = unlist(tcga_subtype), subtype_extended = unlist(tcga_subtype_extended), 
                          trail = unlist(tcga_trail), trail_short = unlist(tcga_trail_short))
  #Following file of TRAIL and TRAILshort expression in TCGA is included in github.
  fwrite(tcga_trail, "TCGA_TRAIL_TPM.csv")
  tcga_trail_median = aggregate(trail_short ~ subtype_extended, tcga_trail, median)
  tcga_trail_median = tcga_trail_median[order(tcga_trail_median$trail_short, decreasing = T), ]
  tcga_trail$subtype_extended = factor(tcga_trail$subtype_extended, levels = tcga_trail_median$subtype_extended)
}

{
  #Expression plot for Fig. 1
  p <- plot_ly(tcga_trail, y = ~trail_short, color = ~subtype_extended, type = "box")
  print(p)
}




