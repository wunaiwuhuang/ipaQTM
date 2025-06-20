library(data.table)

input_dir <- "/data1/wuguojia/data/IPA_QTM_tcga/spqtm/data"
chain_file <- "/data1/wuguojia/data/IPA_QTM_tcga/spqtm/src/00.hg19ToHg38.over.chain.gz"
liftover_path <- "/home/liuxiaochuan/apps/toolbox/liftOver" #use liuxiaochuan's liftover

file_list <- list.files(input_dir, pattern = "_sploc\\.txt$", full.names = TRUE)

for (file in file_list) {
    message("Processing: ", file)
    dat <- fread(file)
    colnames(dat) <- c("id", "chr", "start", "end")
    #bulid bed file
    bed <- data.table(
    chr = dat$chr,
    start = dat$start,
    end = dat$end,
    id = dat$id
    )
    bed <- na.omit(bed)
    #build temp files
    bed_file <- tempfile(fileext = ".bed")
    out_file <- tempfile(fileext = ".bed")
    unmapped_file <- tempfile(fileext = ".bed")
    fwrite(bed, bed_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

    cmd <- sprintf("%s %s %s %s %s", liftover_path, bed_file, chain_file, out_file, unmapped_file)
    system(cmd)

    lifted <- fread(out_file)
    if (nrow(lifted) == 0) {
        warning("No coordinates converted for: ", file)
        next
    }

    colnames(lifted) <- c("chr", "start", "end", "id")

    # 合并原始数据和lifted坐标，按 id 对齐，优先用lifted坐标
    merged <- merge(dat, lifted[, .(id, new_chr = chr, new_start = start, new_end = end)], by = "id", all.x = TRUE)
    
    # 如果有新坐标，就替换；否则用原始
    merged[, final_chr := ifelse(!is.na(new_chr), new_chr, chr)]
    merged[, final_start := ifelse(!is.na(new_start), new_start, start)]
    merged[, final_end := ifelse(!is.na(new_end), new_end, end)]
    
    # 构造输出，保留原始格式顺序
    out <- merged[, .(id, chr = final_chr, start = final_start, end = final_end)]
    
    # 覆盖写入
    fwrite(out, file = file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}

message("All files processed and overwritten with hg38 coordinates.")
