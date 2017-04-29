from __future__ import division
import math

# Genetic diversity parameters

def a1(n):
    '''Given n, this function returns the (n-1)th harmonic number'''
    n = int(n)
    return sum((1.0/d) for d in range(1,n))

def a2(n):
    '''Given n, this function returns the (n-1)th squared harmonic number'''
    n = int(n)
    return sum((1.0/(d**2)) for d in range(1,n))

def b1(n):
    '''Creates b1 for the variance of Tajima's theta'''
    n = int(n)
    return ((n+1) /  (3*(n-1)))

def b2(n):
    '''Creates b2 for the variance of Tajima's theta'''
    n = int(n)
    num = ((n**2) + n + 3) * 2
    den = 9 * n * (n-1)
    return num / den

def c1(n):
    n = int(n)
    '''Creates c1 for the variance of Tajima's theta'''
    return b1(n) - (1 / a1(n))

def c2(n):
    n = int(n)
    '''Creates c2 for the variance of Tajima's theta'''
    return b2(n) - ((n+2) / (a1(n) * n)) + (a2(n) / (a1(n) ** 2 ))

def e1(n):
    n = int(n)
    return c1(n) / a1(n)

def e2(n):
    n = int(n)
    return c2(n) / ( (a1(n) ** 2) + a2(n) )

def cn(n):
    n = int(n)
    num = ( 2 * n * a1(n) ) - (4 * (n-1))
    den = (n-1) * (n-2)
    return num / den

def vD(n):
    n = int(n)
    chunk1 = (a1(n) ** 2) / ( a2(n) + (a1(n) ** 2) )
    chunk2 = cn(n) - ((n + 1) / (n - 1))
    return 1 + (chunk1 * chunk2)

def uD(n):
    n = int(n)
    return a1(n) - 1 - vD(n)

def flatten(foo):
    for x in foo:
        if hasattr(x, '__iter__'):
            for y in flatten(x):
                yield y
        else:
            yield x



# AKA # of segregating sites
def wattersons_theta(pop, seq_length):
    '''Accepts a single dictionary of a population where keys are haplotypes'''
    '''and values are the counts'''
    seq_number = sum(pop.values())
    list_of_seqs =  list(flatten(pop.keys()))
    site_freq_spec = zip(*list_of_seqs)
    set_sfs = [set(sublist) for sublist in site_freq_spec]
    K = 0
    for site in set_sfs:
        if len(site) > 1:
            K += 1
    theta = K / a1(seq_number)
    return theta

def fu_and_li_theta(pop, seq_length):
    seq_number = sum(pop.values())
    # a haplotype occuring once is a prerequisite for a singleton
    # Select haplotypes that occur only once
    singletons = []
    for key, value in pop.iteritems():
        if value == 1:
            singletons.append(key)
    site_freq_spect = zip(*singletons)
    what = zip(*singletons)
    '''Indicates that there's a singleton'''
    site_freq_spec = [x for x in zip(*singletons) if x.count(x[0]) == (len(x) -1) ]
    n_1 = len(site_freq_spec)
    return n_1

def get_distance(seq_a, seq_b):
    diffs = 0
    length = len(seq_a)
    assert len(seq_a) == len(seq_b)
    for chr_a, chr_b in zip(seq_a, seq_b):
        if chr_a != chr_b:
            diffs += 1
    #return diffs / float(length)
    return diffs

# AKA Pi
def tajimas_theta(population):
    '''Accepts a single dictionary of a population where keys are haplotypes \
        and values are the counts'''
    pop_size = sum(population.values())
    haplotypes = population.keys()
    haplotype_count = len(haplotypes)
    diversity = 0
    for i in range(haplotype_count):
        haplotype_a = haplotypes[i]
        frequency_a = population[haplotype_a] / float(pop_size)
        for j in range(0, i):
            #for j in range(haplotype_count):
            haplotype_b = haplotypes[j]
            frequency_b = population[haplotype_b] / float(pop_size)
            frequency_pair = frequency_a * frequency_b
            diversity += frequency_pair * get_distance(haplotype_a, haplotype_b)
    return diversity * 2


def tajimas_D(population, seq_length):
    pop_size = sum(population.values())
    theta = wattersons_theta(population, seq_length)
    pi = tajimas_theta(population)
    S_n = theta * a1(pop_size)
    num = pi - theta
    den = math.sqrt( (e1(pop_size) * S_n) + (e2(pop_size) * S_n * (S_n-1)) )
    if den != 0:
        return num / den
    else:
        return 0

def fu_and_li_D(population, seq_length):
    'Try it with folded SFS data?'
    pop_size = sum(population.values())
    theta = wattersons_theta(population, seq_length)
    a1_pop = a1(pop_size)
    S_n = theta * a1_pop
    num = S_n - ( a1_pop * fu_and_li_theta(population, seq_length))
    den = math.sqrt((uD(pop_size) * S_n) + (vD(pop_size) * (S_n ** 2)))
    if den != 0:
        return num / den
    else:
        return 0
