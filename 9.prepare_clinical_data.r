tumours <- list.files("/data1/wuguojia/data/IPA_QTM_tcga/data/", pattern = "ipause.txt$", full.names = FALSE)
tumours <- gsub("_ipause.txt$", "", tumours)
tumours <- paste0("TCGA-", tumours)
no_info_Biolinks <- c("LAML", "GAM", "LGG", "OV", "PCPG", "PRAD", "SARC", "THYM", "UCEC", "UCS")

all.info <- NULL
for (tumour in tumours) {
  print(tumour)
  
  j <- sub("^TCGA-", "", tumour)
  
  ipa_pt <- file.path("/data1/wuguojia/data/IPA_QTM_tcga/data/", paste0(j, "_ipause.txt"))
  ipa_mat <- read.table(ipa_pt,header=T,sep="\t",stringsAsFactors=F,check.names=FALSE)
  sample_ids <- colnames(ipa_mat)[-1]
  print(length(unique(sample_ids)))
  
  pt <- file.path("/data1/wangwenhui/pu-DNAm/TCGA-GDC-DNAm/05.INFO", paste0("TCGA-", j, ".Xena.info.txt"))
  info <- read.table(pt,header=T,sep="\t",stringsAsFactors=F,check.names=FALSE)
  info <- as.data.frame(t(info))
  missing_ids <- setdiff(sample_ids, rownames(info))
  if (length(missing_ids) > 0) {
    na_rows <- matrix(NA, nrow = length(missing_ids), ncol = ncol(info),
                      dimnames = list(missing_ids, colnames(info)))
    na_df <- as.data.frame(na_rows, stringsAsFactors = FALSE)
    info <- rbind(info, na_df)
  }
  info <- info[sample_ids, ]
  
  if (ncol(info) < 3) {
    info$pathologic_stage.Xena <- NA
  }
  
  colnames(info) <- c("gender.Xena", "age_at_initial_pathologic_diagnosis.Xena", "pathologic_stage.Xena")
  info <- info[, c("gender.Xena", "age_at_initial_pathologic_diagnosis.Xena", "pathologic_stage.Xena")]
  
  info$Tumour <- j
  info$sample_id <- rownames(info)
  info <- unique(info)
  print(length(unique(rownames(info))))
  
  info_file_path <- paste0("/data1/wangwenhui/pu-DNAm/TCGA-GDC-DNAm/", tumour, "/", tumour, ".sample.consistent.info")
  info_biolinks <- read.table(info_file_path,header=T,sep="\t",stringsAsFactors=F,check.names=FALSE)
  info_biolinks <- as.data.frame(info_biolinks)
  
  if (ncol(info_biolinks) < 4) {
    info_biolinks$pathologic_stage.Biolinks <- NA
  } else {
    colnames(info_biolinks)[colnames(info_biolinks) == "pathologic_stage"] <- "pathologic_stage.Biolinks"
  }
  
  colnames(info_biolinks) <- c("sample_id", "age_at_initial_pathologic_diagnosis.Biolinks", "gender.Biolinks", "pathologic_stage.Biolinks")
  info_biolinks <- unique(info_biolinks,by="sample_id")
  info <- merge(info, info_biolinks, by = "sample_id", all.x = TRUE)
  rownames(info) <- info$sample_id
  info <- as.data.frame(info)
  print(length(unique(rownames(info))))
  
  liu_pt <- "/data1/wangwenhui/pu-DNAm/TCGA-GDC-DNAm/03.sample.info.txt"
  liu_info <- read.table(liu_pt,header=T,sep="\t",stringsAsFactors=F,check.names=FALSE,fill = TRUE)
  liu_info <- liu_info[, c("patient_barcode", "gender", "age_at_diagnosis", "pathologic_stage")]
  liu_info <- unique(liu_info)
  liu_info <- as.data.frame(liu_info)
  colnames(liu_info) <- c("sample_id", "gender.Liu", "age_at_initial_pathologic_diagnosis.Liu", "pathologic_stage.Liu")
  
  info <- merge(info, liu_info, by = "sample_id", all.x = TRUE)
  rownames(info) <- info$sample_id
  info <- as.data.frame(info)
  print(length(unique(rownames(info))))
  
  pt <- readLines("/data1/wangwenhui/pu-DNAm/01.input/Covariate_matrix/00.Covariates.pt.list")
  tong_pt <- pt[grepl(j, pt)]
  tong_info <- read.table(tong_pt,header=T,sep="\t",stringsAsFactors=F,check.names=FALSE,fill = TRUE)
  tong_info <- tong_info[c(1,2),]
  tong_info <- as.data.frame(t(tong_info))
  tong_info$sample_id <- rownames(tong_info)
  
  if ("gender" %in% colnames(tong_info)) {
    colnames(tong_info)[colnames(tong_info) == "gender"] <- "gender.Tong"
  } else {
    tong_info$gender.Tong <- NA
  }
  
  age_cols <- grep("age", colnames(tong_info), value = TRUE)
  colnames(tong_info)[colnames(tong_info) %in% age_cols] <- "age_at_initial_pathologic_diagnosis.Tong"
  tong_info <- tong_info[, c("gender.Tong", "age_at_initial_pathologic_diagnosis.Tong", "sample_id")]
  info <- merge(info, tong_info, by = "sample_id", all.x = TRUE)
  info <- as.data.frame(info)
  rownames(info) <- info$sample_id
  info <- info[, -1]
  print(length(unique(rownames(info))))
  
  all.info <- rbind(all.info, info)
}

all.info$gender.Liu <- ifelse(all.info$gender.Liu == "MALE", 1, 
                              ifelse(all.info$gender.Liu == "FEMALE", 0, all.info$gender.Liu))
all.info$gender.Liu <- as.numeric(all.info$gender.Liu)
all.info$pathologic_stage.Liu[all.info$pathologic_stage.Liu == "Stage X"] <- NA
all.info$pathologic_stage.Liu[grepl("^Stage 0$", all.info$pathologic_stage.Liu)] <- "0"
all.info$pathologic_stage.Liu[grepl("^Stage I$|^Stage IA|^Stage IB", all.info$pathologic_stage.Liu)] <- "1"
all.info$pathologic_stage.Liu[grepl("^Stage II$|^Stage IIA|^Stage IIB|^Stage IIC", all.info$pathologic_stage.Liu)] <- "2"
all.info$pathologic_stage.Liu[grepl("^Stage III$|^Stage IIIA|^Stage IIIB|^Stage IIIC", all.info$pathologic_stage.Liu)] <- "3"
all.info$pathologic_stage.Liu[grepl("^Stage IV$|^Stage IVA|^Stage IVB|^Stage IVC", all.info$pathologic_stage.Liu)] <- "4"
all.info$pathologic_stage.Liu <- as.numeric(all.info$pathologic_stage.Liu)

all.info$gender.Biolinks <- ifelse(all.info$gender.Biolinks == "MALE", 1, 
                                   ifelse(all.info$gender.Biolinks == "FEMALE", 0, all.info$gender.Biolinks))
all.info$gender.Biolinks <- as.numeric(all.info$gender.Biolinks)
all.info$pathologic_stage.Biolinks[all.info$pathologic_stage.Biolinks == "Stage X"] <- NA
all.info$pathologic_stage.Biolinks[grepl("^Stage 0$", all.info$pathologic_stage.Biolinks)] <- "0"
all.info$pathologic_stage.Biolinks[grepl("^Stage I$|^Stage IA|^Stage IB", all.info$pathologic_stage.Biolinks)] <- "1"
all.info$pathologic_stage.Biolinks[grepl("^Stage II$|^Stage IIA|^Stage IIB|^Stage IIC", all.info$pathologic_stage.Biolinks)] <- "2"
all.info$pathologic_stage.Biolinks[grepl("^Stage III$|^Stage IIIA|^Stage IIIB|^Stage IIIC", all.info$pathologic_stage.Biolinks)] <- "3"
all.info$pathologic_stage.Biolinks[grepl("^Stage IV$|^Stage IVA|^Stage IVB|^Stage IVC", all.info$pathologic_stage.Biolinks)] <- "4"
all.info$pathologic_stage.Biolinks <- as.numeric(all.info$pathologic_stage.Biolinks)

no_info_Biolinks <- c("LAML", "GAM", "LGG", "OV", "PCPG", "PRAD", "SARC", "THYM", "UCEC", "UCS", "BRCA")
all.info$pathologic_stage.Biolinks[all.info$Tumour %in% no_info_Biolinks] <- NA
head(all.info)

info_file_path <- "/data1/wuguojia/data/IPA_QTM_tcga/00.all.info.txt"
write.table(all.info, file=info_file_path, row.names=TRUE, col.names=TRUE, quote=FALSE, sep='\t')

### sample info unique
info <- read.table("/data1/wuguojia/data/IPA_QTM_tcga/00.all.info.txt",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
gender_cols <- c("gender.Xena", "gender.Biolinks", "gender.Liu", "gender.Tong")
info$gender <- apply(info[, gender_cols], 1, function(row) {
  unique_vals <- unique(row[!is.na(row)])
  if (length(unique_vals) == 0) {
    return(NA)
  } else {
    return(unique_vals[1])
  }
})
head(info[, c(gender_cols, "gender")])

age_cols <- c("age_at_initial_pathologic_diagnosis.Xena", "age_at_initial_pathologic_diagnosis.Biolinks", 
              "age_at_initial_pathologic_diagnosis.Liu", "age_at_initial_pathologic_diagnosis.Tong")
info$age <- apply(info[, age_cols], 1, function(row) {
  unique_vals <- unique(row[!is.na(row)])
  if (length(unique_vals) == 0) {
    return(NA)
  } else {
    return(unique_vals[1])
  }
})
head(info[, c(age_cols, "age")])

info_file_path <- "/data1/wuguojia/data/IPA_QTM_tcga/00.all.info.txt"
write.table(info, file=info_file_path, row.names=TRUE, col.names=TRUE, quote=FALSE, sep='\t')

#examine if all sample is included in info
tumours <- list.files("/data1/wuguojia/data/IPA_QTM_tcga/data/", pattern = "ipause.txt$", full.names = TRUE)
y<- NULL
for(tumour in tumours){
  x<- read.table(tumour,header=T,sep="\t",stringsAsFactors=F,check.names=FALSE)
  y<- append(y,colnames(x)[-1])  
}
x <- read.table("/data1/wuguojia/data/IPA_QTM_tcga/00.all.info.txt",header=T,sep="\t",stringsAsFactors=F,check.names=FALSE)
x <- rownames(x)
setdiff(y,x)
setdiff(x,y)
