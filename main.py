from loading import load_directory
from kmers import stream_kmers, kmer2str
from collections import Counter
import time
import heapq

def minhash(s,listesequence,k):
    L=[]
    for kmer,rkmer in stream_kmers(listesequence,k):
        Kmer=min(kmer,rkmer) 
        Kmer=xorshift(Kmer)
        #if len(L)<s:
            #.append(Kmer)
       #else:
            #.sort()
            #ax=L[-1]
            #f Kmer<Max and Kmer not in L:
                #max=L.index(Max)
                #[imax]=Kmer
        #print(L)
        if len(L) < s :
            heapq.heappush(L, Kmer)    
        else:
            heapq.heappushpop(L, Kmer)



    return sorted(L)

def xorshift(x):
    x ^= x << 13
    x ^= x >> 7
    x ^= x << 17
    return x
def partition_minhash(seq, k, s):
    sketch = []
    for i in range(s):
        sketch.append(float('inf'))
    for kmer,rkmer in stream_kmers(seq,k): 
        ##print(min(kmer,rkmer))
        v=min(kmer,rkmer)
            ##print(v)
        x = xorshift(v)
            ##print(x)
        index = x%s
        #print(index)
        if x <= sketch[index]:
            sketch[index] = x
    return sketch
def index_file(sequences, k):
    return Counter(list_kmers(sequences, k))



def list_file(sequences, k):
    lst = []

    # Construct the list of all the kmers of all the sequences
    for seq in sequences:
        lst.extend([min(x) for x in stream_kmers(seq, k)])

    # Sort the list
    lst.sort()

    return lst

def mymethod_partition_minhas( file1,file2, k,s):
    """ Create an index containing all the kmers of the input sequences
    """
    A=partition_minhash("".join(file1),k,s)
    B=[]
    L=[]
    inter=[]
    A_dict=Counter(A)
    for i in range(s):
            L.append(float('inf'))
    for kmer,rkmer in stream_kmers("".join(file2),k):
          # Initialise the sketch list
        
        x = xorshift(min(kmer,rkmer))
        #print(x)
        index = x%s
        if x <= L[index]:
            L[index] = x
            if x in A_dict:
                inter.append(x)
                if A_dict[x]==1:
                    del A_dict[x]
                else:
                    A_dict[x]=A_dict[x]-1
    
 
    A = len(A) - len(inter)
    L = len(L) - len(inter)
    #print(A)
    #print(B)
    return A, len(inter), L

def mymethod_minhash( file1,file2, k,s):
    """ Create an index containing all the kmers of the input sequences
    """
    A=minhash(s,"".join(file1),k)
    L=[]
    inter=[]
    A_dict=Counter(A)
    for kmer,rkmer in stream_kmers("".join(file2),k):
        Kmer=min(kmer,rkmer) 
        Kmer=xorshift(Kmer)
        if len(L) < s  :
            heapq.heappush(L, Kmer)    
        else:
            heapq.heappushpop(L, Kmer)
            if Kmer in A_dict:
                inter.append(Kmer)
                if A_dict[Kmer]==1:
                    del A_dict[Kmer]
                else:
                    A_dict[Kmer]=A_dict[Kmer]-1

   
    A = len(A) - len(inter)
    L = len(L) - len(inter)
    #print(A)
    #print(B)
    return A, len(inter), L

def intersect_index(index, sequences, k):
    """ Create an index containing all the kmers of the input sequences
    """
    index_uniq = index.copy()
    query_uniq = 0
    intersection = 0

    for seq in sequences:
        for kmer, rkmer in stream_kmers(seq, k):
            minmer = min(kmer, rkmer)

            # Query not in index
            if minmer not in index_uniq:
                query_uniq += 1
            # Query in index => intersection
            else:
                intersection += 1
                index_uniq[minmer] -= 1
                if index_uniq[minmer] == 0:
                    del index_uniq[minmer]

    return sum(index_uniq.values()), intersection, query_uniq


def intersect_sorted_lists(lst1, lst2):
    idx1 = idx2 = 0
    # Dataset specific or intersection kmer counts
    A = inter = B = 0

    while (idx1 < len(lst1) and idx2 < len(lst2)):
        kmer1 = lst1[idx1]
        kmer2 = lst2[idx2]

        # Same kmer => intersection
        if kmer1 == kmer2:
            inter += 1
            idx1 += 1
            idx2 += 1
        # first list specific
        elif kmer1 < kmer2:
            A += 1
            idx1 += 1
        # second list specific
        else:
            B += 1
            idx2 += 1

    # Add remaining kmers of the non empty list
    A += len(lst1) - idx1
    B += len(lst2) - idx2

    return A, inter, B


def list_kmers(sequences, k):
    kmers = []
    for seq in sequences:
        kmers.extend([min(kmer, rkmer) for kmer, rkmer in stream_kmers(seq, k)])
    return kmers


def similarity(A, inter, B):
    # +1 added for pseudocount. Avoid divisions by 0
    A_similarity = inter / (inter + A + 1)
    B_similarity = inter / (inter + B + 1)

    return A_similarity, B_similarity


def jaccard(A, inter, B):
    return inter / (A + inter + B)



if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")

    k = 21

    # Loading
    indexes = {f:index_file(files[f], k) for f in files}
    lists = {f:list_file(files[f], k) for f in files}
    
    filenames = list(files.keys())
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            # Method 1 using index
            print("index")
            start=time.time()
            A, inter, B = intersect_index(indexes[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter, B))
            end=time.time()
            print("Time spent:", end-start)

     # Method 2 using sorted lists
        #print("lists")
        #start=time.time()
        #A, inter, B = intersect_sorted_lists(lists[filenames[i]], lists[filenames[j]])
        #end=time.time()
        #print("Time spent:", end-start)

        #print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter, B))

    for i in range(len(files)):
        for j in range(i+1, len(files)):
            print(" hash de base")

            start_t=time.time()
            A, inter, B =mymethod_minhash(files[filenames[i]], files[filenames[j]], k,10000)
            end_t=time.time()
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter,B))
            print("Time spent:", end_t-start_t)
            print("hash casier")

            start_t=time.time()
            A, inter, B=mymethod_partition_minhas(files[filenames[i]], files[filenames[j]], k,10000)
            end_t=time.time()

            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter,B))
            print("Time spent:", end_t-start_t)
