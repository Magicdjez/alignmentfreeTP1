from os import listdir, path


def load_fna(filename):
    texts = []
    txt = []

    with open(filename) as fp:
        for line in fp:
            if line[0] == '>':
                if len(txt) > 0:
                    texts.append("".join(txt))
                txt = []
            else:
                txt.append(line.strip())

    if len(txt) > 0:
        texts.append("".join(txt))
    return texts



def load_directory(directory):
    files = {}
    for filename in listdir(directory):
        if filename[filename.rfind('.')+1:] in ["fa", "fasta", "fna"]:
            files[filename] = load_fna(path.join(directory, filename))
    
    return files


if __name__ == "__main__":
    files = load_directory("data")
    print(len(files))
