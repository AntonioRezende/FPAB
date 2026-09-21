#!/usr/local/bin/Rscript

amino_acidos_massa <- c(
  'A' = 71.03711,
  'C' = 103.0092,
  'D' = 115.0269,
  'E' = 129.0426,
  'F' = 147.0684,
  'G' = 57.02146,
  'H' = 137.0589,
  'I' = 113.0841,
  'K' = 128.0950,
  'L' = 113.0841,
  'M' = 131.0405,
  'N' = 114.0429,
  'P' = 97.05276,
  'Q' = 128.0586,
  'R' = 156.1011,
  'S' = 87.03203,
  'T' = 101.0477,
  'V' = 99.06841,
  'W' = 186.0793,
  'Y' = 163.0633
)

cat("Digite sua sequencia de Aminoácidos\n")
seqPTN <- readLines(con = "stdin", n = 1)
seqPTN <- toupper(seqPTN)

tamanho <- nchar(seqPTN)
peso <- 0
AA <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
seq_chars <- strsplit(seqPTN, NULL)[[1]]

if (all(seq_chars %in% AA)) {
  cat(sprintf("Sequencia da Proteina = %s\n", seqPTN))
  for (i in 1:tamanho) {
    peso <- peso + amino_acidos_massa[seq_chars[i]]
  }
  cat(sprintf("O peso molecular da sequencia é igual a %.2f Da\n", peso))
} else {
  cat("Sequencia possui caracteres que não são aminoácidos\n")
  invalid_chars <- setdiff(seq_chars, AA)
  cat(paste(invalid_chars, collapse = ", "), "\n")
}