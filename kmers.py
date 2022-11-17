
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)


def stream_kmers(text, k):
    list_Kmer=[]
    kmer=0
    rkmer=0
    for i in range(k-1):#Initialisation du premier kmer
        kmer=kmer<<2
        kmer+=encode(text[i])
        rkmer>>2
        rkmer+=encode(text[i])
    for nucl in text[k-1:]:
        mask=(1<<(k-1)*2)-1#Masque à appliquer
        kmer=kmer & mask
        kmer=kmer<<2# décalage de 2 bit pour le prochain nucleotide
        kmer=kmer+encode(nucl)#ajout du dernier nucleotide
        rkmer>>2
        rkmer+=(encode((nucl))<<((k-1)*2))#creation du reverse Kmer

        list_Kmer.append(min(kmer,rkmer))# max entre le kmer et le reverse
    return list_Kmer# liste des kmer de la strinf


def encode(string):# return la valeur du nucleotide assocé 
    dico={"A":0,"C":1,"T":2,"G":3}
    if string not in dico:
        return dico["A"]
    return dico[string]
