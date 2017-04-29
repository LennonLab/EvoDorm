from __future__ import division
import numpy as np
import itertools, os, math, argparse, random, time
from collections import Counter
import pandas as pd


mydir = os.path.realpath(__file__).rsplit('/', 2)[0]

def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    x = Counter(x)
    y = Counter(y)
    z = x + y
    return z


def sweepSimDorm(N, M, c, s, reps = 1000):
    count = 0
    N_freqs_list = []
    M_freqs_list = []
    while count < reps:
    # B = selected allele
        Mb = M
        MB = 0
        NB = 1
        Nb = N-1
        N_freqs = [(NB/N)]
        M_freqs = [(MB/M)]
        gs = [1]
        g = 1
        K = N/M
        d = (c* K) / N
        r = (c  ) / M
        while (NB + MB) < (N + M) and (NB + MB) != 0:
            # transistion probs https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2206092/
            Nb_Nb =  (1-  (NB/N)) * (1-  (NB/N)) * (1-d)
            Nb_NB = (NB/N) * (1- (NB/N)) * (1-d)
            NB_Nb = (NB/N) * (1- (NB/N)) * (1-d) * (1-s)
            NB_NB = ( (NB/N) * (NB/N) * (1-d) ) + ((NB/N) * (1- (NB/N)) * s * (1-d) )

            Nb_Mb = (1 - (NB/N)) * d
            Nb_MB = 0
            NB_MB = (NB/N) * d
            NB_Mb = 0

            Mb_Nb = (1 - (MB/M)) * r
            Mb_NB = 0
            MB_NB = (MB/M) * r
            MB_Nb = 0

            Mb_Mb = (1- (MB/M)) * (1-r)
            Mb_MB = 0
            MB_Mb = 0
            MB_MB = (MB/M) * (1-r)

            # chose b and B
            #print Nb_Nb + NB_Nb + Mb_Nb + MB_Nb + Nb_NB + NB_NB + Mb_NB + MB_NB
            #print Nb_Mb + NB_Mb + Mb_Mb + MB_Mb + Nb_MB + NB_MB + Mb_MB + MB_MB
            pop_N = list(np.random.choice(2, N, p=[Nb_Nb + NB_Nb + Mb_Nb + MB_Nb, \
                Nb_NB + NB_NB + Mb_NB + MB_NB]))
            pop_M = list(np.random.choice(2, M, p=[Nb_Mb + NB_Mb + Mb_Mb + MB_Mb, \
                Nb_MB + NB_MB + Mb_MB + MB_MB]))

            Nb = pop_N.count(0)
            NB = pop_N.count(1)
            Mb = pop_M.count(0)
            MB = pop_M.count(1)
            N_freqs.append((NB / N) )
            M_freqs.append((MB / M) )
            g += 1
            gs.append(g)
            # if the advantagous allele is lost from the population,
            # restart the simulation
            if (NB + MB) == 0:
                break
        if (NB + MB)  != (N + M):
            continue
        else:
            print N, M, c, s, count
            N_freqs_list.append(N_freqs)
            M_freqs_list.append(M_freqs)
            count +=1
    N_df = pd.DataFrame(N_freqs_list)
    M_df = pd.DataFrame(M_freqs_list)
    N_df.to_csv(mydir + '/data/Fig3b/N_sweep_N_' + str(N) + '_M_' + str(M) + '_c_' + \
        str(c) + '_s_' + str(s) + '_r_' + str(reps) + '.txt', header=False, index = False)
    M_df.to_csv(mydir + '/data/Fig3b/M_sweep_N_' + str(N) + '_M_' + str(M) + '_c_' + \
        str(c) + '_s_' + str(s) + '_r_' + str(reps) + '.txt', header=False, index = False)



class WF_dorm:

    def __init__(self, N = 100, M = 10, u = 0.0001, G = 100, g = 100, c = 10000, \
        random = False):
        self.N = N
        self.M = M
        self.u = u
        self.G = G
        self.c = int(c)
        self.g = g
        self.alphabet = ['A', 'T', 'G', 'C']
        self.random = random

    def generate_base_haplotype(self):
        base_haplotype = ''.join(["A" for i in range(self.G)])
        return base_haplotype

    def random_dna_sequence(self):
        return ''.join(np.random.choice(self.alphabet) for _ in range(self.G))

    def generate_pop(self):
        pop = {}
        pop['Active'] = {}
        pop['Dormant'] = {}
        if self.random == True:
            for x in range(self.N):
                new_haplotype = self.random_dna_sequence()
                if new_haplotype in pop['Active']:
                    pop['Active'][new_haplotype] += 1
                else:
                    pop['Active'][new_haplotype] = 1
            for x in range(self.M):
                new_haplotype = self.random_dna_sequence()
                if new_haplotype in pop['Dormant']:
                    pop['Dormant'][new_haplotype] += 1
                else:
                    pop['Dormant'][new_haplotype] = 1

        else:
            base_haplotype = self.generate_base_haplotype()
            pop['Active'][base_haplotype] = self.N
            pop['Dormant'][base_haplotype] = self.M
        return pop

    def choose_by_weight(self, weights):
        weights = map(float, weights)
        rndm = random.random() * sum(weights)
        for i, j in enumerate(weights):
            rndm -= j
            if rndm < 0:
                return i

    # mutation
    def get_mutation_count(self):
        mean = self.u * self.N * self.G
        return np.random.poisson(mean)

    def get_random_haplotype(self, pop, sub_pop):
        # Need random haplotype for active only
        #except KeyError:
        haplotypes = pop[sub_pop].keys()
        size_step = sum(pop[sub_pop].values())
        frequencies = [x/float(size_step) for x in pop[sub_pop].values()]
        total = sum(frequencies)
        frequencies = [x / total for x in frequencies]
        return np.random.choice(haplotypes, p=frequencies)

    def get_mutant(self, haplotype):
        site = np.random.randint(self.G)
        possible_mutations = ['A', 'C', 'G', 'T']
        possible_mutations.remove(haplotype[site])
        mutation = np.random.choice(possible_mutations)
        new_haplotype = haplotype[:site] + mutation + haplotype[site+1:]
        return new_haplotype

    def mutation_event(self, pop):
        haplotype = self.get_random_haplotype(pop, 'Active')
        if pop['Active'][haplotype] > 1:
            pop['Active'][haplotype] -= 1
            new_haplotype = self.get_mutant(haplotype)
            if new_haplotype in pop['Active']:
                pop['Active'][new_haplotype] += 1
            else:
                pop['Active'][new_haplotype] = 1

    def mutation_step(self, pop):
        mutation_count = self.get_mutation_count()
        for i in range(mutation_count):
            self.mutation_event(pop)

    # reproduce active pop, drift acts here
    def get_offspring_counts(self, pop):
        haplotypes = pop['Active'].keys()
        frequencies = [x/float(self.N) for x in pop['Active'].values()]
        total = sum(frequencies)
        frequencies = [x / total for x in frequencies]
        return list(np.random.multinomial(self.N, frequencies))

    def offspring_step(self, pop):
        counts = self.get_offspring_counts(pop)
        for (haplotype, count) in zip(pop['Active'].keys(), counts):
            if (count > 0):
                pop['Active'][haplotype] = count
            else:
                del pop['Active'][haplotype]

    def dormancy_step(self, pop):
        if self.M <= 0:
            pass
        else:
            K = self.N / self.M
            for i in range(self.c):
                new_haplotype = self.get_random_haplotype(pop, 'Active')
                pop['Active'][new_haplotype] -= 1
                if new_haplotype in pop['Dormant']:
                    pop['Dormant'][new_haplotype] += 1
                else:
                    pop['Dormant'][new_haplotype] = 1
            for i in range(self.c):
                new_haplotype = self.get_random_haplotype(pop, 'Dormant')
                pop['Dormant'][new_haplotype] -= 1
                if new_haplotype in pop['Active']:
                    pop['Active'][new_haplotype] += 1
                else:
                    pop['Active'][new_haplotype] = 1

    def time_step(self, pop):
        if self.u != 0:
            self.mutation_step(pop)
        self.offspring_step(pop)
        self.dormancy_step(pop)



    def clean_empty(self, d):
        if not isinstance(d, (dict, list)):
            return d
        if isinstance(d, list):
            return [v for v in (self.clean_empty(v) for v in d) if v]
        return {k: v for k, v in ((k, self.clean_empty(v)) for k, v in d.items()) if v}

    def simulate(self):
        pop = self.generate_pop()
        for i in range(self.g):
            #print i, self.c
            self.time_step(pop)
            pop = self.clean_empty(pop)
            #print pop
            #pop = dict((k, v) for k, v in pop['Active'].iteritems() if v > 0)
            #print pop
            #pop = dict((k, v) for k, v in pop['Dormant'].iteritems() if v > 0)
        return pop

    def simulate_save_history(self):
        for i in range(generations):
            self.time_step(pop)
            return pop


def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print '%s function took %0.3f ms' % (f.func_name, (time2-time1)*1000.0)
        return ret
    return wrap
