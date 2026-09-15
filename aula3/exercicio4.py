#!/usr/local/bioinformatic/conda3/bin/python

codigo_genetico_dna = {
    # Fenilalanina / Leucina / Sinais de Parada
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    
    # Leucina / Prolina / Histidina / Glutamina / Arginina
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    
    # Isoleucina / Metionina (Início) / Treonina / Asparagina / Lisina / Serina / Arginina
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    
    # Valina / Alanina / Ácido Aspártico / Ácido Glutâmico / Glicina
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

arquivoDNA = input("Digite o nome do arquivo de entrada\n")

contaLinhas=0
seq=""
cabecalho=""
with open(arquivoDNA,'r') as fasta:
    for linha in fasta:
        contaLinhas=contaLinhas+1
        if(contaLinhas==2):
            seq=linha.strip()
        elif(contaLinhas==1):
            cabecalho=linha.strip()

#len(linha)
tamanho = len(seq)
#print(tamanho)
proteina=[]

for i in range(0,tamanho-2,3):
    nt1=seq[i]
    nt2=seq[i+1]
    nt3=seq[i+2]
    codon=nt1+nt2+nt3
    #print(codon)
    proteina.append(codigo_genetico_dna[codon])

arquivo_saida = input("Digite o nome do arquivo de saida\n")
with open(arquivo_saida,"w") as o:
    o.write(cabecalho+"\n")
    o.write("".join(proteina))

print(f"{cabecalho}")
print("".join(proteina))