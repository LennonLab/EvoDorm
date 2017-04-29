from __future__ import division
import simulationCode as sc
#import makeFigs as mf
import os, math, argparse
from collections import Counter
import numpy as np
from subprocess import call


def gen_log_space(limit, n):
    result = [1]
    if n>1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result)<n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # problem! same integer. we need to find next_value by artificially incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining values will scale correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    return np.array(map(lambda x: round(x)-1, result), dtype=np.uint64)


class run_simulation_hgt:

    def __init__(self, N, M, u, G, g, reps = 100):
        self.N = N
        self.M = M
        self.u = u
        self.G = G
        self.g = g
        self.reps = reps

    def run_HGT_dorm(self, c, R):
        OUT = open(os.path.realpath(__file__).rsplit('/', 2)[0] + \
        '/data/Fig5/G' + str(self.G) + '_U' + str(self.u) + '_R' + str(R) + \
        '_N' + str(self.N) + '_M' + str(self.M) + '_c' + str(c) + '.txt', 'w+')
        print>> OUT, 'replicate', 'pan_genome', 'core_genome'
        for rep in range(self.reps):
            print c, R, rep
            simulation = sc.HGT_dorm(N = self.N, M = self.M, U = self.u, \
                G = self.G, R = R, gens = self.g, c = c).get_gene_counts()
            pan = simulation[0]
            core = simulation[1]
            print>> OUT, rep, pan, core
        OUT.close()

    def submit_mason_jobs(self):
        cs = np.linspace(1, 100, num = 20, endpoint=True)
        rs = np.logspace(-5, -1, num = 20, base = 10)
        cs = np.rint(cs)
        #cs = cs[15:]
        rs = rs[0:1]
        for c in cs:
            for r in rs:
                c = int(c)
                OUT_folder_path = os.path.realpath(__file__).rsplit('/', 2)[0] + '/bash/mason_hgt_c/'
                OUT_path = OUT_folder_path + 'mason_hgt_c' + str(c) + '_rec' +  str(r) + '.sh'
                if os.path.exists(OUT_path) == True:
                    os.remove(OUT_path)
                OUT = open(OUT_path, 'w')
                print>> OUT, '#!/bin/bash'
                print>> OUT, '#PBS -k o'
                print>> OUT, '#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=30:00:00'
                print>> OUT, '#PBS -M wrshoema@umail.iu.edu'
                print>> OUT, '#PBS -m abe'
                print>> OUT,  '#PBS -j oe'
                print>> OUT,  ''
                print>> OUT, 'module load python'
                print>> OUT,  ''
                print>> OUT,  'python ' + os.path.realpath(__file__) + ' -5 -s' + \
                    ' -N ' + str(self.N) + ' -M ' + str(self.M) + ' -u ' + str(self.u) \
                    + ' -G ' + str(self.G) + ' -g ' + str(self.g) + ' -R ' + str(r) \
                    + ' -c ' + str(c) + ' -r ' + str(self.reps)

                OUT.close()
                #call('sh ' + OUT_path, shell=True)
                call('qsub ' + OUT_path, shell=True)


class run_simulation_evol_dist:

    def __init__(self, N, M, u, G, g, reps = 1000):
        self.N = N
        self.M = M
        self.u = u
        self.G = G
        self.g = g
        self.reps = reps

    def calc_evol_dist(self, seqs):
        # jukes cantor corrected evolutionary distance
        ref = list(sc.WF_dorm(G = self.G).generate_base_haplotype())
        # takes a list of seqs and calculates the number of substitutions
        zipped_seqs = zip(*seqs)
        sites_keep = []
        k = 0
        for x, site in enumerate(zipped_seqs):
            #print len(set(site))
            if len(set(site)) == 1:
                sites_keep.append((x, site[0]))
        for x in sites_keep:
            #print x
            if ref[x[0]] != x[1]:
                k += 1
        d = (-3/4)* math.log( 1 - ((4/3) * (k/ len(ref))) )
        return (k / len(ref), abs(d))

    def run_WF_dorm(self, c):
        OUT_data = open(os.path.realpath(__file__).rsplit('/', 2)[0] + \
            '/data/Fig6/Fig6_sim/Fig6_sim_c_' + str(c)
            + '.txt', 'w')
        print>> OUT_data, 'active_size', 'dormant_size', 'mut_rate_per_base', \
            'genome_size', 'gens', 'dormancy_number', 'replicate', \
            'sub_rate', 'evol_distance', 'sub_rate_N', 'evol_distance_N', \
            'sub_rate_M', 'evol_distance_M'
        for rep in range(self.reps):
            pop = sc.WF_dorm(N = self.N, M = self.M, \
                u = self.u, G = self.G, g = self.g, \
                c = c).simulate()
            if args.dormant_size > 0:
                pop_merged = sc.merge_two_dicts(pop['Active'], pop['Dormant'])
            else:
                pop_merged = pop
            data = self.calc_evol_dist(pop_merged)
            data_N = self.calc_evol_dist(pop['Active'])
            if 'Dormant' in pop:
                data_M = self.calc_evol_dist(pop['Dormant'])
                print>> OUT_data, self.N, self.M, \
                    self.u, self.G, self.g, \
                    c, rep, data[0], data[1], \
                    data_N[0], data_N[1], data_M[0], data_M[1]
            else:
                print>> OUT_data, self.N, self.M, \
                    self.u, self.G, self.g, \
                    c, rep, data[0], data[1], \
                    data_N[0], data_N[1], float('nan'), float('nan')

        OUT_data.close()

    #  N, M, u, G, g, c_max = 0, r = 1000
    def submit_mason_jobs(self):
        # seed bank size  = 100, max number that can switch is 100
        # self.gen_log_space for logarithmically-spaced integers
        #cs = gen_log_space(100, 50)
        # we only go from 1 to 100, so we can use linspace and just round the numbers
        cs = np.linspace(1, self.M, num = 50, endpoint=True)
        #cs = [1, 2]
        #cs = np.linspace(1, 2, num = 1, endpoint=True)
        cs = np.rint(cs)
        for i, c in enumerate(cs):
            c = int(c)
            OUT_folder_path = os.path.realpath(__file__).rsplit('/', 2)[0] + '/bash/mason_subs_c/'
            OUT_path = OUT_folder_path + 'mason_sub_' + str(c) + '.sh'
            if os.path.exists(OUT_path) == True:
                os.remove(OUT_path)
            OUT = open(OUT_path, 'w')
            print>> OUT, '#!/bin/bash'
            print>> OUT, '#PBS -k o'
            print>> OUT, '#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=30:00:00'
            print>> OUT, '#PBS -M wrshoema@umail.iu.edu'
            print>> OUT, '#PBS -m abe'
            print>> OUT,  '#PBS -j oe'
            print>> OUT,  ''
            print>> OUT, 'module load python'
            print>> OUT,  ''
            print>> OUT,  'python ' + os.path.realpath(__file__) + ' -6 -s' + \
                ' -N ' + str(self.N) + ' -M ' + str(self.M) + ' -u ' + str(self.u) \
                + ' -G ' + str(self.G) + ' -g ' + str(self.g) + ' -c ' + str(c) \
                + ' -r ' + str(self.reps)

            OUT.close()
            #call('sh ' + OUT_path, shell=True)
            call('qsub ' + OUT_path, shell=True)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "run simulations")
    parser.add_argument('-m', '--run_mason', help='run simulation on Mason',
        action='store_true')
    parser.add_argument('-s', '--run_simulation', help='Print more data',
        action='store_true')
    parser.add_argument('-d', '--run_simulation_locally', help='Print more data',
        action='store_true')

    parser.add_argument('-N', '--active_size', type = int, default = 100)
    parser.add_argument('-M', '--dormant_size', type = int, default = 10)
    parser.add_argument('-u', '--mut_rate_per_base', type = float, default = 0.0001)
    parser.add_argument('-G', '--genome_size', type = int, default = 100)
    parser.add_argument('-g', '--generations', type = int, default = 100)
    parser.add_argument('-c', '--number_entering_dormancy', type = int, default = 0)
    parser.add_argument('-r', '--replicates', type = int, default = 1000)
    parser.add_argument('-R', '--recomb_rate_per_base', type = float, default = 0.01)

    parser.add_argument('-1', '--Fig1', action='store_true')
    parser.add_argument('-2', '--Fig2', action='store_true')
    parser.add_argument('-3', '--Fig3', action='store_true')
    parser.add_argument('-4', '--Fig4', action='store_true')
    #parser.add_argument('-5', '--Fig5', action='store_true')
    #parser.add_argument('-6', '--Fig6', action='store_true')

    args = parser.parse_args()
    if args.Fig1 == True:
        print "This is a conceptual figure and includes no simulated data"
    elif args.Fig2 == True:
        Ms = [10, 100, 1000]
        c = 10
        for M in Ms:
            df = pd.DataFrame()
            for x in range(args.replicates):
                OUT = open(os.path.realpath(__file__).rsplit('/', 2)[0] + \
                    '/data/Fig2b/G' + str(args.generations) + '_S' + \
                    str(self.genome_size) + '_N' + str(self.active_size) + \
                    '_M' + M +'_c' + str(c) + '_I' + str(x) + '.txt', 'w+')

                print>> OUT, 'iteration', 'WT', 'PI', 'FT', 'TD', 'FD'
                pop = sc.WF_dorm(N = args.active_size, M = M, u = args.mut_rate_per_base, \
                    G = args.genome_size, g = args.generations, c = c).generate_pop()
                for g in range(args.generations):
                    sc.WF_dorm(N = args.active_size, M = M, u = args.mut_rate_per_base, \
                        G = args.genome_size, g = args.generations, c = c).time_step(pop)
                    pop_merged = gp.merge_two_dicts(pop['Active'], pop['Dormant'])
                    FD = qp.fu_and_li_D(pop_merged, args.genome_size)
                    TD = qp.tajimas_D(pop_merged, args.genome_size)
                    WT = qp.wattersons_theta(pop_merged, args.genome_size)
                    FT = qp.fu_and_li_theta(pop_merged, args. genome_size)
                    PI = qp.tajimas_theta(pop_merged)
                    print>> OUT, g, WT, PI, FT, TD, FD
                print M, x
                OUT.close()
                IN = pd.read_csv(os.path.realpath(__file__).rsplit('/', 2)[0] + \
                    '/data/Fig2b/G' + str(args.generations) + '_S' + \
                    str(self.genome_size) + '_N' + str(self.active_size) + \
                    '_M' + M +'_c' + str(c) + '_I' + str(x) + '.txt', sep = ' ')
                df[x] = IN.PI
            out_path = os.path.realpath(__file__).rsplit('/', 2)[0] + \
                '/data/Fig2b/merged/G' + str(args.generations) + '_S' + \
                str(self.genome_size) + '_N' + str(self.active_size) + \
                '_M' + M +'_c' + str(c) + '_Pi.txt'
            df.to_csv(path_or_buf = out_path, sep = ' ', index_label = None, \
                index = False, columns = None)

    elif args.Fig3 == True:
        #mf.Fig4()
        pass
    elif args.Fig6 == True:
        if args.run_mason == True:
            run_simulation_Fig5(N = args.active_size, M = args.dormant_size, \
                u = args.mut_rate_per_base, G = args.genome_size, \
                g = args.generations, reps = args.replicates).submit_mason_jobs()
        elif args.run_simulation == True:
            run_simulation_Fig5(N = args.active_size, \
                M = args.dormant_size, u = args.mut_rate_per_base, \
                G = args.genome_size, g = args.generations, reps = args.replicates \
                ).run_HGT_dorm(c = args.number_entering_dormancy, R = args.recomb_rate_per_base)
        elif args.run_simulation_locally == True:
            cs = np.linspace(1, args.dormant_size, num = 50, endpoint=True)
            for c in cs:
                c = int(c)
                run_simulation_Fig5(N = args.active_size, \
                    M = args.dormant_size, u = args.mut_rate_per_base, \
                    G = args.genome_size, g = args.generations, reps = args.replicates \
                    ).run_HGT_dorm(c = c, R = args.recomb_rate_per_base)

    elif args.Fig4 == True:
        if args.run_mason == True:
            '''Makes a mason job for a given vaule of c'''
            '''A max of N individuals can enter the seed bank each gen'''
            run_simulation_evol_dist(N = args.active_size, M = args.dormant_size, \
                u = args.mut_rate_per_base, G = args.genome_size, g = args.generations, \
                reps = args.replicates).submit_mason_jobs()
        elif args.run_simulation == True:
            run_simulation_evol_dist(N = args.active_size, \
                M = args.dormant_size, u = args.mut_rate_per_base, \
                G = args.genome_size, g = args.generations, \
                reps = args.replicates).run_WF_dorm(args.number_entering_dormancy)
        elif args.run_simulation_locally == True:
            cs = np.linspace(1, self.M, num = 50, endpoint=True)
            cs = np.rint(cs)
            for c in cs:
                c = int(c)
                run_simulation_evol_dist(N = args.active_size, \
                    M = args.dormant_size, u = args.mut_rate_per_base, \
                    G = args.genome_size, g = args.generations, \
                    reps = args.replicates).run_WF_dorm(c)
