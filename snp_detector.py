import Bio.SeqIO
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse

FILENAME = "salmonella-enterica.reads.fna"
VARIANT_FILENAME = "salmonella-enterica-variant.reads.fna"
K = 12
LEVEL = 35

class ReadsDictionary:
    def __init__(self):
        self.dictionary = {} #dictionary whose keys are k-mer and values are number of occurences of the k-mer in the reads

    def __init__(self, file):
        """Creates a ReadsDictionary and fills it with all the k-mer present on the file"""
        directory = os.getcwd()
        reads = Bio.SeqIO.parse(f"{directory}/{file}","fasta")
        it_read = iter(reads)
        self.dictionary = {}
        for read in it_read:
            item = read.seq
            self.fill_dict(item)

    def draw_histogram(self, k, nb_seq):
        """Draws the histogram of values and gives the statistics of the dict of k-mer, with nb_seq sequences"""
        list_hist = list(self.dictionary.values())
        print(f"k={k}, number of sequences = {nb_seq}:")
        print(f"    mean = {np.mean(list_hist)}, median = {np.median(list_hist)}, var = {np.var(list_hist)}")
        plt.figure()
        plt.hist(list_hist, bins='auto')
        plt.title(f"{k}-mer histogram for {nb_seq} sequences")
        plt.show()

    def count_sequences(self, k):
        """Counts the occurences of the k-mer on the dict contained in the reads"""
        directory = os.getcwd()
        reads = Bio.SeqIO.parse(f"{directory}/{FILENAME}","fasta")
        it_read = iter(reads)
        for current_read in it_read:
            item = current_read.seq
            for i in range(len(item)-k+1):
                k_mer = get_k_mer(item, k, i)
                if k_mer in self.dictionary:
                    self.dictionary[k_mer] += 1

    def fill_dict(self, current_read):
        """Fills the ReadsDictionary with current_read k-mer"""
        for i in range(0, len(current_read)-K+1):
            k_m = get_k_mer(current_read, K, i)
            if k_m in self.dictionary:
                self.dictionary[k_m] += 1
            else:
                self.dictionary[k_m] = 1

    def clear_dict(self):
        """Deletes errors of reading of the ReadsDictionary"""
        to_delete = []
        for sequence, occurences in self.dictionary.items():
            if occurences <= LEVEL:
                to_delete.append(sequence)
        for sequence in to_delete:
            del self.dictionary[sequence]


def get_histogram(k, nb_seq):
    """Creates a dictionary containing nb_seq k-mer and their occurences,
    and draws the histogram corresponding"""
    directory = os.getcwd()
    reads = Bio.SeqIO.parse(f"{directory}/{FILENAME}","fasta")
    read = next(reads).seq
    dictionary = ReadsDictionary()
    for _ in range(nb_seq):
        sequence_aleatoire = read[:k]
        while sequence_aleatoire in dictionary.dictionary:
            read = next(reads).seq
            sequence_aleatoire = read[:k]
        dictionary.dictionary[sequence_aleatoire] = 0
        read = next(reads).seq
    dictionary.count_sequences(k)
    dictionary.draw_histogram(k, nb_seq)


def get_histograms(args):
    """Draws histograms for nb_k_mer k-mers between k1 and k2"""
    for i in range(args.k1, args.k2):
        get_histogram(i, args.nb_k_mer)


def poisson_distribution():
    """Shows the histograms for Poisson distributions with 100000 values"""
    nb_reads = 1993166
    lam = nb_reads*(250 - K)/(5000000 - 250)
    epsilon = 0.01    
    n = 100000
    x = np.random.poisson(lam,n)
    y = np.random.poisson(lam*epsilon/4, n)
    valeurs = [i for i in range (max(x))]
    plt.hist([x,y], bins=valeurs, density=True)
    plt.show()


def get_k_mer(read, k, start):
    """Gives the k-mer starting at start on the read"""
    k_mer = read[start:start+k]
    return k_mer


def init_dict():
    """Return two ReadsDictionary containing solid k-mers"""
    healthy_dict = ReadsDictionary(FILENAME)
    mutated_dict = ReadsDictionary(VARIANT_FILENAME)
    healthy_dict.clear_dict()
    mutated_dict.clear_dict()
    return healthy_dict, mutated_dict


def compare_dict(healthy_dict, mutated_dict):
    """Compare healthy and mutated ReadsDictionarys and gives the ReadsDictionarys
    with the sequences which are not common between both"""
    listSNPhealth = []
    for healthy_key in healthy_dict.dictionary.keys():
        if healthy_key not in mutated_dict.dictionary.keys():
            listSNPhealth.append(healthy_key)
        else :
            del mutated_dict.dictionary[healthy_key]
    return listSNPhealth, mutated_dict.dictionary


def get_snps(args):
    """Returns the different k-mer of the two fna files"""
    healthy_dict, mutated_dict = init_dict()
    health_not_found_kmer, mutated_not_found_kmer = compare_dict(healthy_dict, mutated_dict)
    list_mutated = get_list(mutated_not_found_kmer)
    write_fasta(list_mutated, "mutated.fasta")


def write_fasta(mutated_seq, output_file):
    """Writes the sequences into a fasta file"""
    fasta_out = open(output_file, 'w')
    for seq_i in range(len(mutated_seq)):
        fasta_out.write('>seq' +str(seq_i))
        fasta_out.write('\n')
        fasta_out.write(str(mutated_seq[seq_i]))
        fasta_out.write('\n')

def get_list(dict):
    """Returns the keys of a dictionary into a list"""
    liste_k_mers = dict.keys()
    return list(liste_k_mers)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
    help='Available commands', dest='subcommand')

    # Histogram
    subparser = subparsers.add_parser("histograms", help="Draws histograms for nb_k_mer k-mers between k1 and k2")
    subparser.add_argument("k1", type=int)
    subparser.add_argument("k2", type=int)
    subparser.add_argument("nb_k_mers", type=int)
    subparser.set_defaults(func=get_histograms)
    
    # Poisson distribution
    subparser = subparsers.add_parser("poisson", help="Shows the histograms for Poisson distributions with 100000 values")
    subparser.set_defaults(func=poisson_distribution)

    #Comparison
    subparser = subparsers.add_parser("comparison", help="Returns the different k-mer of the two fna files")
    subparser.set_defaults(func=get_snps)

    cmd_args = parser.parse_args()
    res = cmd_args.func(cmd_args)