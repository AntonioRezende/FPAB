#!/usr/local/bioinformatic/conda3/bin/python

amino_acidos_massa = {
    'A': 71.03711,
    'C': 103.0092,
    'D': 115.0269,
    'E': 129.0426,
    'F': 147.0684,
    'G': 57.02146,
    'H': 137.0589,
    'I': 113.0841,
    'K': 128.095,
    'L': 113.0841,
    'M': 131.0405,
    'N': 114.0429,
    'P': 97.05276,
    'Q': 128.0586,
    'R': 156.1011,
    'S': 87.03203,
    'T': 101.0477,
    'V': 99.06841,
    'W': 186.0793,
    'Y': 163.0633
}

seqPTN = input("Digite sua sequencia de Aminoácidos\n")
seqPTN = seqPTN.upper()

tamanho = len(seqPTN)
#print(tamanho)
peso=0
AA = set("ACDEFGHIKLMNPQRSTVWY")

if(set(seqPTN).issubset(AA)):
    print(f"Sequencia da Proteina = {seqPTN}")
    for i in range(0,tamanho):
        peso = peso + amino_acidos_massa[seqPTN[i]]
   
    print(f"O peso molecular da sequencia é igual a {peso:.2f} Da")
else:
    print("Sequencia possui caracteres que não são aminoácidos")
    print(set(seqPTN).difference(AA))
