from loading import load_directory
from kmers import stream_kmers, kmer2str,encode
from collections import Counter
import time

def similarity(C, inter):# calcul de la similarité
    sim = inter / (C + inter)
    return sim
def jaccard(A, inter, B):# calcul de jaccard
    union = A + B + inter
    index = inter / union
    return index
def my_method(file1, file2, k):
    #input_f1 = 'data//{}'.format(file1) # path and name of the file
    #input_f2 = 'data//{}'.format(file2)
    #content1 = readfile(input_f1) # all the sequences merged in the file
    #content2 = readfile(input_f2)
    Kmer_list1 = stream_kmers("".join(file1),k)# creation des listes de kmer
    Kmer_list2 = stream_kmers("".join(file2),k)# creation des listes de kmer
    inter = len(set(Kmer_list1).intersection(Kmer_list2)) 
    '''We use the intersection of sets method to ensure the fastest possible runtime. 
    With k = 21, it is almost impossible for two elements in a set to be identical, 
    and even if a few do occur, they will not have much effect on the similarity score 
    due to the large length of the data''' 
    A = len(Kmer_list1) - inter
    B = len(Kmer_list2) - inter  
    return A, inter, B
def my_methodv2(file1,file2,k):
    #input_f1 = 'data//{}'.format(file1) # path and name of the file
    #input_f2 = 'data//{}'.format(file2)
    #content1 = readfile(input_f1) # all the sequences merged in the file
    #content2 = readfile(input_f2)
    inter=[]
    A=stream_kmers("".join(file1),k)
    A_dict=Counter(A)
    B=[]
    kmer=0
    rkmer=0
    seq="".join(file2)
    mask=(1<<((k-1)*2))-1
    revmask=(1<<(k*2))-1-3
    scores={"A":0,"C":1,"T":2,"G":3}
    rscores={"A":2,"C":3,"T":0,"G":1}
    for i in range(k-1):
        kmer=kmer<<2
        kmer+=encode(seq[i])
        rkmer>>2
        rkmer+=encode(seq[i])
    for nucl in seq[k-1:]:
        mask=(1<<(k-1)*2)-1#Masque à appliquer
        kmer=kmer & mask
        kmer=kmer<<2# décalage de 2 bit pour le prochain nucleotide
        kmer=kmer+encode(nucl)#ajout du dernier nucleotide
        rkmer>>2
        rkmer+=(encode((nucl))<<((k-1)*2))#creation du reverse Kmer
        B.append(min(kmer,rkmer))
        if min(kmer,rkmer) in A_dict:
            inter.append(min(kmer,rkmer))
            if A_dict[min(kmer,rkmer)]==1:
                del A_dict[min(kmer,rkmer)]
            else:
                A_dict[min(kmer,rkmer)]=A_dict[min(kmer,rkmer)]-1
    #print(len(inter))
    #print(len(A))
    #print(len(B))
    A = len(A) - len(inter)
    B = len(B) - len(inter)
    #print(A)
    #print(B)
    return A, len(inter), B

if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    filenames = list(files.keys())
    print("premiere methode")
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            start_t=time.time()
            A, inter, B =my_method(files[filenames[i]], files[filenames[j]], k)
            end_t=time.time()

            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter), similarity(B, inter))

            print("Time spent:", end_t-start_t)
    print("deuxieme methode")

    for i in range(len(files)):
        for j in range(i + 1, len(files)):
            start_t=time.time()
            A, inter, B =my_methodv2(files[filenames[i]], files[filenames[j]], k)
            end_t=time.time()
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter), similarity(B, inter))

            print("Time spent:", end_t-start_t)


