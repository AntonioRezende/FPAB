#!/usr/local/bin/Rscript

cat("Digite sua sequencia de DNA\n")
seqDNA <- readLines(con = "stdin", n = 1)
seqDNA <- toupper(seqDNA)

tamanho <- nchar(seqDNA)
revDNA <- c()
dnaLetters <- c("A", "T", "C", "G")
seq_chars <- strsplit(seqDNA, NULL)[[1]]

if (all(seq_chars %in% dnaLetters)) {
  cat(sprintf("Sequencia de DNA = %s\n", seqDNA))
  cat("Reverso Complementar = ")
  
  for (i in tamanho:1) {
    nt <- seq_chars[i]
    if (nt == "A") {
      revDNA <- c(revDNA, "T")
    } else if (nt == "G") {
      revDNA <- c(revDNA, "C")
    } else if (nt == "C") {
      revDNA <- c(revDNA, "G")
    } else {
      revDNA <- c(revDNA, "A")
    }
  }
  
  cat(paste0(revDNA, collapse = ""), "\n")
  cat("Sequencia de invertida com sucesso\n")
} else {
  cat("Sequencia possui caracteres invalidos para ser DNA\n")
}