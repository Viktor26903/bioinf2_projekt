import sys

filename = sys.argv[1]

def match(seq1, seq2, index1, index2, direction):
    count = 0
    if direction == "lijevo":
        for i in range(index1, index2+1):
            if seq1[i] == seq2[i]:
                count -= 1
        for i in range(index1, index2+1):
            if seq1[i] == seq2[i-1]:
                count += 1
    else:
        for i in range(index1, index2+1):
            if seq1[i] == seq2[i]:
                count -= 1
        for i in range(index1, index2+1):
            if seq1[i] == seq2[i+1]:
                count += 1
    return count

def provjeri_susjede(seq1, seq2, index1, index2):
    if index1 > 0 and index2 < min(len(seq1), len(seq2)) - 1:
        if seq1[index1-1] == "-" and seq1[index2+1] != "-" and seq2[index1-1] != "-" and seq2[index2+1] == "-":
            if match(seq1, seq2, index1, index2, "lijevo") >= 0:
                del seq1[index1 - 1]
                del seq2[index2 + 1]
                provjeri_susjede(seq1, seq2, index1-1, index2)
        elif seq1[index1-1] != "-" and seq1[index2+1] == "-" and seq2[index1-1] == "-" and seq2[index2+1] != "-":
            if match(seq1, seq2, index1, index2, "desno") >= 0:
                del seq1[index2 + 1]
                del seq2[index1 - 1]
                provjeri_susjede(seq1, seq2, index1-1, index2)
    return seq1, seq2, index2

def poravnaj(seq1,seq2):
    i = 0
    while i < min(len(seq1),len(seq2)):
        if seq1[i] == seq2[i] == "-":
            del seq1[i]
            del seq2[i]
            continue
        i += 1
    i = 0
    while i < min(len(seq1),len(seq2)):
        if i == 0:
            i += 1
            continue
        if seq1[i] == "-" and seq2[i] != "-" and seq1[i-1] != "-" and seq2[i-1] == "-":
            del seq1[i]
            del seq2[i-1]
            i -= 1
            seq1, seq2, i = provjeri_susjede(seq1, seq2, i, i)
        elif seq2[i] == "-" and seq1[i] != "-" and seq2[i-1] != "-" and seq1[i-1] == "-":
            del seq2[i]
            del seq1[i-1]
            i -= 1
            seq1, seq2, i = provjeri_susjede(seq1, seq2, i, i)
        i += 1
    return seq1,seq2

sequences = {}
with open(filename) as file:
    for line in file:
        if line[0] == ">":
            name=line[1:].strip()
            sequences[name] = ""
        else:
            sequences[name] += line.strip()

output = sys.argv[2]
with open(output,"w") as out:
    for ind1, name1 in enumerate(sequences):
        for ind2, name2 in enumerate(sequences):
            if ind2 <= ind1:
                continue
            seq1 = list(sequences[name1])
            seq2 = list(sequences[name2])
            seq1,seq2 = poravnaj(seq1,seq2)
            out.write(f">{name1}\t{name2}\n{''.join(seq1)}\n{''.join(seq2)}\n")
        

