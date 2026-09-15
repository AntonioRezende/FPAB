#!/usr/local/bioinformatic/conda3/bin/python

seqDNA = input("Digite sua sequencia de DNA\n")
seqDNA = seqDNA.upper()

tamanho = len(seqDNA)
#print(tamanho)

revDNA=[]
#revDNA.append("A")
#print(revDNA)

dnaLetters=set("ATCG")

if(set(seqDNA).issubset(dnaLetters)):
    print(f"Sequencia de DNA = {seqDNA}")
    print("Reverso Complementar = ", end="")
    for i in range((tamanho-1),-1,-1):
        #print(seqDNA[i])
        if (seqDNA[i] =="A"):
            revDNA.append("T")
        elif(seqDNA[i] == "G"):
            revDNA.append("C")
        elif(seqDNA[i] == "C"):
            revDNA.append("G")
        else:
            revDNA.append("A")

    print("".join(revDNA))
    print("Sequencia de invertida com sucesso")
else:
    print("Sequencia possui caracteres invalidos para ser DNA")
