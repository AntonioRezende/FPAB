#!/usr/local/bin/Rscript

codigo_genetico_dna <- c(
  'TTT' = 'F', 'TTC' = 'F', 'TTA' = 'L', 'TTG' = 'L',
  'TCT' = 'S', 'TCC' = 'S', 'TCA' = 'S', 'TCG' = 'S',
  'TAT' = 'Y', 'TAC' = 'Y', 'TAA' = '*', 'TAG' = '*',
  'TGT' = 'C', 'TGC' = 'C', 'TGA' = '*', 'TGG' = 'W',
  
  'CTT' = 'L', 'CTC' = 'L', 'CTA' = 'L', 'CTG' = 'L',
  'CCT' = 'P', 'CCC' = 'P', 'CCA' = 'P', 'CCG' = 'P',
  'CAT' = 'H', 'CAC' = 'H', 'CAA' = 'Q', 'CAG' = 'Q',
  'CGT' = 'R', 'CGC' = 'R', 'CGA' = 'R', 'CGG' = 'R',
  
  'ATT' = 'I', 'ATC' = 'I', 'ATA' = 'I', 'ATG' = 'M',
  'ACT' = 'T', 'ACC' = 'T', 'ACA' = 'T', 'ACG' = 'T',
  'AAT' = 'N', 'AAC' = 'N', 'AAA' = 'K', 'AAG' = 'K',
  'AGT' = 'S', 'AGC' = 'S', 'AGA' = 'R', 'AGG' = 'R',
  
  'GTT' = 'V', 'GTC' = 'V', 'GTA' = 'V', 'GTG' = 'V',
  'GCT' = 'A', 'GCC' = 'A', 'GCA' = 'A', 'GCG' = 'A',
  'GAT' = 'D', 'GAC' = 'D', 'GAA' = 'E', 'GAG' = 'E',
  'GGT' = 'G', 'GGC' = 'G', 'GGA' = 'G', 'GGG' = 'G'
)

cat("Digite o nome do arquivo de entrada\n")
arquivoDNA <- readLines(con = "stdin", n = 1)

linhas <- readLines(arquivoDNA)
cabecalho <- trimws(linhas[1])
seq <- trimws(linhas[2])

tamanho <- nchar(seq)
proteina <- c()

for (i in seq(1, tamanho - 2, by = 3)) {
  nt1 <- substr(seq, i, i)
  nt2 <- substr(seq, i + 1, i + 1)
  nt3 <- substr(seq, i + 2, i + 2)
  codon <- paste0(nt1, nt2, nt3)
  proteina <- c(proteina, codigo_genetico_dna[codon])
}

cat("Digite o nome do arquivo de saida\n")
arquivo_saida <- readLines(con = "stdin", n = 1)

output_text <- c(cabecalho, paste0(proteina, collapse = ""))
writeLines(output_text, arquivo_saida)

cat(sprintf("%s\n", cabecalho))
cat(paste0(proteina, collapse = ""), "\n")