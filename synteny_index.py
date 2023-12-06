from Bio import SeqIO
import os
import sys
import pandas as pd
import math
from Bio import pairwise2
import random
import numpy as np

'''
Requirements:
- path directory of all genomes
- annotation files for genomes
- nucleotide files for genes
- list of suspected genes (input is protein accession number)
- orthofinder Orthogroups.tsv output
'''


class gene:
    def __init__(self, id, sequence, pos, orthogroup):
        self.id = id
        self.sequence = sequence
        self.pos = pos
        self.orthogroup = orthogroup

    def set_sequence(self, sequence):
        self.sequence = sequence

    


class synteny:
    def __init__(self, genomes, ntd, orthogroups):

        self.genomes = os.listdir(genomes) # genomes folder
        self.genomes_index = {}
        self.genomes.remove(".DS_Store")
        self.genomes.sort()
        self.ntd_files = os.listdir(ntd)
        self.ntd_files.remove(".DS_Store")
        self.ntd_files.sort()
        self.gene_pos = [] # list of genes in order
        self.orthogroups = orthogroups # link to orthogroups file
        self.orthogroup_list = {}
        self.gene_data = [] # for each genome, have a dictionary. for each dictionary, each gene id corresponds to a sequence object that stores orthogroup info, pos info, etc.

    def init_genes(self):
        self.read_annotation_file()
        self.read_ntd_files()
        self.initialize_orthogroup()

    def init_from_genome_objects(self, genomes):
        ind = 0
        for genome in genomes:
            self.gene_pos.append([])
            self.gene_data.append({})
            self.genomes_index[ind+1] = ind
            i = 0
            for g in genome.genes:
                self.gene_pos[ind].append(g.id)
                self.gene_data[ind][g.id] = gene(g.id, g.sequence, i, g.id)
                if g.id > 1000:
                    self.gene_data[ind][g.id].orthogroup = g.id // 1000
                if self.orthogroup_list.get(g.id) == None and g.id <= 1000:
                    self.orthogroup_list[g.id] = {}
                elif g.id > 1000 and self.orthogroup_list.get(g.id // 1000) == None:
                    self.orthogroup_list[g.id // 1000] = {}
                if g.id <= 1000:
                    if self.orthogroup_list[g.id].get(ind) == None:
                        self.orthogroup_list[g.id][ind] = []
                else:
                    if self.orthogroup_list[g.id // 1000].get(ind) == None:
                        self.orthogroup_list[g.id // 1000][ind] = []
                if g.id <= 1000:
                    self.orthogroup_list[g.id][ind].append(g.id)
                else:
                    self.orthogroup_list[g.id//1000][ind].append(g.id)
                #self.orthogroup_list[g.id].append(g.id)
                i +=1
            ind+=1

    def read_annotation_file(self):
        print("Reading annotation files...")
        id = 0
        for genome in self.genomes:
            genes_1 = pd.read_csv("annotations/"+ genome, sep="\t", comment='#')
            self.gene_data.append({})
            self.genomes_index[genome[:-4]] = id
            for index, row in genes_1.iterrows(): # fill up gene position dictionary
                if row[2] == "CDS":
                    gene_id = row[8].split(";")[0].split("=")[1]
                    start = int(row[3])
                    end = int(row[4])
                    if (self.gene_data[id].get(gene_id[4:])) == None:
                        self.gene_data[id][gene_id[4:]] = gene(None, None, None, None)
                    #self.gene_pos[index][gene_id] = (start, end)
                    self.gene_data[id][gene_id[4:]].pos = (start, end)
            id += 1
        for genome in self.gene_data:
            sorted_genes = sorted(genome, key=lambda k: genome[k].pos[0])
            self.gene_pos.append(sorted_genes)

    def read_ntd_files(self):
        print("Reading nucleotide files...")
        
        index = 0
        for file in self.ntd_files:
            if file.endswith("fna") or file.endswith("faa") or file.endswith("fa") or file.endswith("fasta"):
                seqs = SeqIO.parse("ntd/"+ file, "fasta")
                for seq in seqs:
                    try:
                        gene_id = seq.description.split("protein_id=")[1].split("]")[0]
                    except IndexError:
                        gene_id = seq.id.split("_cds_")[1]
                    if self.gene_data[index].get(gene_id) == None:
                        self.gene_data[index][gene_id] = gene(None, None, None, None)
                    #self.gene_seqs[index][gene_id] = str(seq.seq)
                    self.gene_data[index][gene_id].sequence = str(seq.seq)
                index+=1 

    def si_sig_probability(self, threshold, genome1, genome2, k):
        return math.ceil((k/(len(genome1)+len(genome2))) * (2 * (len(genome1)+len(genome2)) - len(genome1.symmetric_difference(genome2)))-math.sqrt(-k*math.log(threshold)))

    def apply_JC_correction(self, hamming):
        try:
            if hamming > 3/4:
                raise Exception("hamming distance must be below 3/4")
        except ValueError as e:
            return None
        return -(0.75)*math.log(1-(4/3)*hamming)
    
    def calc_expected_hamming(self, jc_corrected):
        return 0.75*(1-math.e**(-(4/3)*jc_corrected))
    
    def calc_hamming_distance(self, gene1, gene2):
        alignments = pairwise2.align.globalxx(gene1, gene2)
        alignment = alignments[0]
        seq1 = alignment.seqA
        seq2 = alignment.seqB
        distance = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
        return distance / len(alignment.seqA)

    def initialize_orthogroup(self):
        print("Initializing orthogroups...")
        orthogroups = pd.read_csv(self.orthogroups, sep="\t")
        for index, row in orthogroups.iterrows(): # fill up orthogroup dictionary
            if self.orthogroup_list.get(row[0]) == None:
                self.orthogroup_list[row[0]] = {}
            for i in range(1, len(row)):
                if not (pd.isna(row[i])):
                    list_genes = row[i].split(", ")
                    
                    for gene in list_genes:
                        ind=0
                        for genome in self.gene_data:
                            
                            if genome.get(gene) != None:
                                if self.orthogroup_list[row[0]].get(ind) == None:
                                    self.orthogroup_list[row[0]][ind] = []
                                genome[gene].orthogroup = row[0]
                                self.orthogroup_list[row[0]][ind].append(gene)
                                ind += 1
                                break
                            ind +=1
                            

    def calc_expected_distance(self, hamming_distance_gr, hamming_distance_rw, hamming_distance_sw):
        return (self.apply_JC_correction(hamming_distance_gr) / self.apply_JC_correction(hamming_distance_rw)) * self.apply_JC_correction(hamming_distance_sw)

    def calc_egs(self, n, deltar):
        return math.sqrt(-math.log(deltar / 2) / (2 * n)) 
    
    def calc_distance_diff(self, hamming_gs, hamming_distance_gr, hamming_distance_sw, hamming_distance_rw):
        return abs(hamming_gs - self.calc_expected_hamming(self.calc_expected_distance(hamming_distance_gr, hamming_distance_rw, hamming_distance_sw)))

    def calc_SI(self, gene, genome1, genome2, k):
        '''
        input: genome name, gene suspected to be horizontally transferred, k neighborhood
        '''
        genome1_ind = self.genomes_index[genome1]
        genome2_ind = self.genomes_index[genome2]
        if self.gene_data[genome1_ind][gene].orthogroup == None:
            #print("Error: no homologs found for " +  gene + " - could not calculate SI")
            return
        list_homologs = self.orthogroup_list[self.gene_data[genome1_ind][gene].orthogroup] # dictionary
        genome1_list = self.gene_pos[genome1_ind]
        genome2_list = self.gene_pos[genome2_ind]
        ind1 = genome1_list.index(gene)
        ind2 = None
        try:
            homolog = list_homologs[genome2_ind][0]
            ind2 = genome2_list.index(homolog)
            if len(list_homologs[genome2_ind]) > 1:
                homolog = list_homologs[genome2_ind][1]
                ind2 = genome2_list.index(homolog)
        except:
            pass
        if ind2 == None:
            print("ERROR: synteny index could not be calculated for " + gene + " - no homologs found in genome 2")
            return
        neighborhood1 = genome1_list[max(0,ind1-math.floor(k)):min(ind1+math.ceil(k),len(genome1_list))]
        neighborhood2 = genome2_list[max(0,ind2-math.floor(k)):min(ind2+math.ceil(k),len(genome2_list))]
        ct = 0
        #converts gene names into orthogroups
        for g in neighborhood1:
            neighborhood1[ct] = self.gene_data[genome1_ind][g].orthogroup
            ct+=1
        ct = 0
        for g in neighborhood2:
            neighborhood2[ct] = self.gene_data[genome2_ind][g].orthogroup
            ct+=1
        

        return len(set(item for item in neighborhood1 if item is not None) & set(item for item in neighborhood2 if item is not None))
    
    def get_HGT_suspected_genes(self, genome1, genome2, k, max_prob=0.05):
        print("obtaining HGT-suspected genes...")
        genome1_ind = self.genomes_index[genome1]
        genome2_ind = self.genomes_index[genome2]
        ct = 0
        G1 = self.gene_pos[genome1_ind][:]
        G2 = self.gene_pos[genome2_ind][:]
        
        #converts gene names into orthogroups
        for g in G1:
            G1[ct] = self.gene_data[genome1_ind][g].orthogroup
            ct+=1
        ct = 0
        for g in G2:
            G2[ct] = self.gene_data[genome2_ind][g].orthogroup
            ct+=1
        id = 0
        for ind in range(len(G1)):
            if G1[ind] == None:
                G1[ind] = id
                id += 1
        for ind in range(len(G2)):
            if G2[ind] == None:
                G2[ind] = id
                id += 1
        #print(G1)
        rec_k = - (math.log(max_prob)) / (4 * (1-(len(set(G1).symmetric_difference(set(G2)))/(len(set(G1)) + len(set(G2))))))
        rec_k = math.ceil(k)
        if k < rec_k:
            print("ERROR: minimum k size needs to be " + str(rec_k))
            return
        thres = self.si_sig_probability(max_prob,set(G1),set(G2),k)
        print("Calculated maximum significant SI is " + str(thres))
        HGT_suspected = []
        ind = self.genomes_index[genome1]
        for gene in self.gene_pos[0]:
            SI = self.calc_SI(gene, genome1, genome2,k)
            if SI != None:
                if SI < thres:
                    HGT_suspected.append(gene)
                    #print(gene + " - SI index: " + str(SI))       
                    #
        return HGT_suspected 
    
    def test_sig_HGT(self, gene, witness, strain1, strain2, ref1, ref2, n, deltar):
        strain1_ind = self.genomes_index[strain1]
        strain2_ind = self.genomes_index[strain2]
        ref1_ind = self.genomes_index[ref1]
        ref2_ind = self.genomes_index[ref2]
        list_homologs = self.orthogroup_list[self.gene_data[strain1_ind][gene].orthogroup]
        list_witness_homologs = self.orthogroup_list[self.gene_data[strain1_ind][witness].orthogroup]
        if len(list_homologs.keys()) < 4:
            print("Error: HGT suspected gene not present in all reference organisms and strains")
            return
        if len(list_witness_homologs.keys()) < 4:
            print("Error: witness gene not found in all reference organisms/strains")
            return
        print("gene being tested: " + str(gene))
        if gene == 100:
            print(self.gene_data[strain1_ind][gene].sequence)
            print(self.gene_data[strain2_ind][gene].sequence)
        hamming_gs = self.calc_hamming_distance(self.gene_data[strain1_ind][gene].sequence, self.gene_data[strain2_ind][list_homologs[1][0]].sequence)
        hamming_gr = self.calc_hamming_distance(self.gene_data[ref1_ind][list_homologs[2][0]].sequence, self.gene_data[ref2_ind][list_homologs[3][0]].sequence)
        hamming_ws = self.calc_hamming_distance(self.gene_data[strain1_ind][witness].sequence, self.gene_data[strain2_ind][list_witness_homologs[1][0]].sequence)
        hamming_wr = self.calc_hamming_distance(self.gene_data[ref1_ind][list_witness_homologs[2][0]].sequence, self.gene_data[ref2_ind][list_witness_homologs[3][0]].sequence)
        if hamming_gs > 3/4 or hamming_gr > 3/4 or hamming_ws > 3/4 or hamming_wr > 3/4:
            print("Error: hamming distance is too high")
            return
        if hamming_gs == 0 or hamming_gr == 0 or hamming_ws == 0 or hamming_wr == 0:
            print("Error: hamming distance is too low")
            return
        print(hamming_gs)
        egs = self.calc_egs(n, deltar)
        corr = self.calc_distance_diff(hamming_gs, hamming_gr, hamming_ws, hamming_wr)
        print("difference: " + str(corr) + " | egs threshold: " + str(egs))
        return True if corr > egs else False
    
    def comprehensive_test(self, strain1, strain2, ref1, ref2, n, deltar, k, use_suggested_n=False, sample_size=100):
        witness_genes = self.get_all_witness_genes(strain1, strain2, ref1, ref2)
        random_witness_sample = random.sample(witness_genes, sample_size)
        suspected_HGT = self.get_HGT_suspected_genes(strain1, strain2, k)
        print(suspected_HGT)
        print("corrected threshold given number of witness genes: " + str(deltar / len(suspected_HGT)))
        print("Number of HGT tested: " + str(len(suspected_HGT)))
        if use_suggested_n:
            n = len(suspected_HGT)
        putative_HGT = []
        for gene in suspected_HGT:
            ct = 0
            for w in random_witness_sample:
                if w != gene:
                    outcome = self.test_sig_HGT(gene, w, strain1, strain2, ref1, ref2, n, deltar / len(random_witness_sample))
                    if outcome: 
                        print("Witness gene " + str(w) + " produced significant results for " + str(gene))
                        ct += 1
                        if ct > 1:
                            print(str(gene) + " has passed significant threshold " + str(len(random_witness_sample) //2 ) +  " times. is deemed putative HGT")
                            putative_HGT.append(gene)
                            break
                    if outcome == None:
                        pass
        print("putative HGT: " + str(putative_HGT))
                    
        return putative_HGT
    
    def get_all_witness_genes(self, genome1, genome2, genome3, genome4):
        '''This functions gets all of the possible witness genes of a genome
        (criteria: must be present in all genomes, and must not be the HGT suspected gene)
        '''
        genome1_ind = self.genomes_index[genome1]
        witness_genes = []
        for orthogroup in self.orthogroup_list:
            if len(self.orthogroup_list[orthogroup].keys()) >= 4:
                for g in self.orthogroup_list[orthogroup][genome1_ind]:
                    witness_genes.append(g)
        return witness_genes




if __name__ == "__main__":
    test = synteny("annotations", "ntd", "orthogroups/Orthogroups.tsv")
    test.init_genes()
    #print(test.gene_data[3]["QRP91555.1"].sequence)
    #print(test.gene_data[3]["QRP91555.1"].pos)
    #print(test.gene_data[3]["QRP91555.1"].orthogroup)


    #HGT_suspected = test.get_HGT_suspected_genes("GCA_000008865.2", "MG1655", 100)
    #HGT_candidates = []
    #print(len(HGT_suspected))
    #for gene in HGT_suspected:
    #    res = test.test_sig_HGT(gene, "BAB38645.1", "GCA_000008865.2", "MG1655", "Wolbachia", "fragilis", 1472, 0.05)
    #    if res != None:
    #        print(gene + " res: "+ str(res))


    h1 = 0.0237
    h2 = 0.583
    h3 = 0.541
    h4 = 0.008
    lala = test.calc_expected_distance(h2,h3,h4)
    print(test.calc_expected_distance(h2,h3,h4))
    print(test.calc_expected_hamming(lala))
    #suspected = test.get_HGT_suspected_genes("GCA_000008865.2", "MG1655", 10)
    #print(test.calc_SI("BAB36796.2", "GCA_000008865.2", "MG1655", 10))
    #print(test.test_sig_HGT("BAB36796.2", "BAB37946.1", "GCA_000008865.2", "MG1655", "Wolbachia", "fragilis", 1472, 0.05))
    test.comprehensive_test("1_isoflavoniconvertens", "2_equolifaciens", "3_Lactococcus", "4_Enteroscipio", 1472, 0.05, 10, sample_size=10)