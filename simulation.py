#import synteny_index
import random
import numpy as np
import math
import copy
import synteny_index

class gene:
    def __init__(self, id, sequence, mutation_rate):
        self.id = id
        self.sequence = sequence
        self.mutation_rate = mutation_rate

class genome:
    def __init__(self, hgt_multiplier):
        self.genes = []
        self.hgt_multipler = hgt_multiplier

    def initialize(self, k, len, mut_rate_mean, mut_rate_sd):
        #iniailize genome with k genes
        for i in range(k):
            g = gene(None, None, None)
            g.id = i
            nucleotides = ['A', 'T', 'C', 'G']
            g.sequence= (''.join(random.choice(nucleotides) for _ in range(len)))
            g.mutation_rate = np.random.normal(mut_rate_mean, mut_rate_sd)
            self.genes.append(g)

    #def go_through_one_gen(self):
    #    for seq_ind in range(len(self.genes)):
    #        mut = self.genes[seq_ind].mutation_rate
    #        for i in range(len(self.genes[seq_ind].sequence)):
    #            if random.random() < mut:
    #                choice = random.choice(['A', 'T', 'C', 'G'])
    #                self.genes[seq_ind].sequence = self.genes[seq_ind].sequence[:i] + choice + self.genes[seq_ind].sequence[i+1:]
    #        return
    def go_through_one_gen(self):
        new_genes = []
        for g in self.genes:
            new_gene = copy.deepcopy(g)
            mut = new_gene.mutation_rate
            new_sequence = list(new_gene.sequence)

            for i in range(len(new_sequence)):
                if random.random() < mut:
                    new_sequence[i] = random.choice(['A', 'T', 'C', 'G'])

            new_gene.sequence = ''.join(new_sequence)
            new_genes.append(new_gene)

        self.genes = new_genes

    
    #swap genes
    def genome_arrangement(self):
        ind = random.randrange(0,len(self.genes))
        k = math.floor(np.random.normal(5, 2))
        ind2 = random.randrange(0,len(self.genes))
        while ind2 + k >= ind - k and ind2 - k <= ind + k:
            ind2 = random.randrange(0,len(self.genes))
        
        region1 = self.genes[ind - k:ind + k]
        region2 = self.genes[ind2 - k:ind2 + k]

    # Swap the two regions
        self.genes[ind - k:ind + k] = region2
        self.genes[ind2 - k:ind2 + k] = region1

    def HGT(self, genome2, gene):
        ''' parameters: genome2 object, gene object from this gene'''
        for i in genome2.genes:
            if i.id == gene.id:
                genome2.genes.remove(i)
        c = random.randrange(0, len(genome2.genes))
        newgene = copy.deepcopy(gene)
        newgene.id = newgene.id
        genome2.genes.insert(c, newgene)
        newgene.mutation_rate = newgene.mutation_rate * 1000
        print(genome2.genes[c].mutation_rate)

    def write_fasta(self, filename):
        with open(filename, "w") as fasta:
            for i in self.genes:
                fasta.write(">"+ str(i.id) + "\n")
                fasta.write(i.sequence + "\n\n")

    def copy_genome(self):
        new_genome = genome(self.hgt_multipler)
        for g in self.genes:
            new_gene = gene(g.id, g.sequence[:], g.mutation_rate)
            new_genome.genes.append(new_gene)
        return new_genome



def main():
    genome1 = genome(3)
    genome1.initialize(1000, 1000, 0.000003, 0.000001)
    genome2 = genome1.copy_genome()
    genome3 = genome1.copy_genome()
    genome4 = genome1.copy_genome()

    for i in range(100): # 100 generations
        genome1.go_through_one_gen()
        genome2.go_through_one_gen()
        genome3.go_through_one_gen()
        genome4.go_through_one_gen()
    
    genome1.write_fasta("genome1_pre.fasta")
    genome2.write_fasta("genome2_pre.fasta")
    genome3.write_fasta("genome3_pre.fasta")
    genome4.write_fasta("genome4_pre.fasta")

    genome2.genome_arrangement()
    genome2.genome_arrangement()
    genome2.genome_arrangement()
    genome2.genome_arrangement()

    hgt_genes = random.sample(list(range(1000)), 100)
    for i in hgt_genes:
        genome1.HGT(genome2, genome1.genes[i])
    print(hgt_genes)

    for i in range(50):
        genome1.go_through_one_gen()
        genome2.go_through_one_gen()
        genome3.go_through_one_gen()
        genome4.go_through_one_gen()

    genome1.write_fasta("genome1_post.fasta")
    genome2.write_fasta("genome2_post.fasta")
    genome3.write_fasta("genome3_post.fasta")
    genome4.write_fasta("genome4_post.fasta")

    simulation = synteny_index.synteny("annotations", "ntd", "orthogroups/Orthogroups.tsv")
    simulation.init_from_genome_objects([genome1, genome2, genome3, genome4])
    simulation.comprehensive_test(1, 2, 3, 4, 1000, 0.05, 10, sample_size=10, use_suggested_n=True)



main()