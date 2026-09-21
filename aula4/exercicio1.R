#!/usr/local/bin/Rscript

cat("Digite sua sequencia de DNA\n")
seqDNA <- readLines(con = "stdin", n = 1)
seqDNA <- toupper(seqDNA)

cat(sprintf("Sequencia de DNA = %s\n", seqDNA))
cat("Sequencia de RNA = ")

chars <- strsplit(seqDNA, NULL)[[1]]
for (i in chars) {
  if (i == "T") {
    cat("U")
  } else {
    cat(i)
  }
}

cat("\n")
cat("Sequencia de DNA convertida para RNA\n")
