#!/usr/local/bioinformatic/conda3/bin/python

seqDNA = input("Digite sua sequencia de DNA\n")
seqDNA = seqDNA.upper()

seqRNA = ()

print(f"Sequencia de DNA = {seqDNA}")
print("Sequencia de RNA = ", end="")
for i in seqDNA:
    if (i =="T"):
        print("U",end="")
    else:
        print(i, end="")

print("")
print("Sequencia de DNA convertida para RNA")