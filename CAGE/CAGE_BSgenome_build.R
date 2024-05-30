# make bs genome for CAGEr
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/pericles/TEfind/database/")
library(Biostrings)
library(BSgenome)
file <- system.file("extdata",package = "BSgenome")
seqnames <- c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY")

forgeSeqFiles("flybase","dm6",seqnames, prefix = "", suffix = ".fasta.gz",
              seqs_srcdir = file, seqs_destdir = tempdir(), ondisk_seq_format ="rds")
forgeSeqFiles("flybase","dm6",seqnames, prefix = "", suffix = ".fasta.gz",
              seqs_srcdir = file, seqs_destdir = tempdir(), ondisk_seq_format ="2bit")


dm6 <- import(file.path(tempdir(),"single_sequences.2bit"))
forgeBSgenomeDataPkg("/media/pericles/TEfind/BSgenome_seed")


