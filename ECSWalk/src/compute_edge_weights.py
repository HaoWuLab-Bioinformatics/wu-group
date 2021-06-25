import networkx as nx
import os
import argparse
import scipy.stats as st
import numpy as np

#compute mutex values
from src.utils import *


def compute_ed():
    if args.verbose:
        print ("computing ed values...")

    G = nx.Graph()

    for e in edge_list:
         G.add_edge(e[0], e[1])

    G.to_undirected()
    G.remove_edges_from(G.selfloop_edges())

    N = len(genes)
    num_samples = len(load_patients_to_indices())

    weight_out_file = ed_out_file
    fhout = open(weight_out_file, 'w+')

    count_1s = 0
    count_0s = 0

    for i in range(len(edge_list)):
        e = edge_list[i]
        gene1_index = e[0]
        gene2_index = e[1]
        gene1 = id_to_gene[gene1_index]
        gene2 = id_to_gene[gene2_index]

        gene1_neighbours = G[gene1_index]
        gene2_neighbours = G[gene2_index]

        union_set1 = set()
        union_set2 = set()
        count1 = 0
        count2 = 0

        if gene1 in data:
            union_set1 = set(data[gene1])
            count1 = len(data[gene1])
        if gene2 in data:
            union_set2 = set(data[gene2])
            count2 = len(data[gene2])

        for gene in gene1_neighbours:

            if id_to_gene[gene] in data and id_to_gene[gene] != gene2_index:
                union_set1 = union_set1.union(data[id_to_gene[gene]])
                count1 += len(data[id_to_gene[gene]])
        for gene in gene2_neighbours:
            if id_to_gene[gene] in data and id_to_gene[gene] != gene1_index:
                union_set2 = union_set2.union(data[id_to_gene[gene]])
                count2 += len(data[id_to_gene[gene]])

        union_set1 = len(union_set1)
        union_set2 = len(union_set2)

        m1 = 0
        m2 = 0
        if count1 != 0:
            m1 = float(union_set1) / count1
        if count2 != 0:
            m2 = float(union_set2) / count2

        ed = (m1+m2)/2

        if ed == 1:
           count_1s += 1
        fhout.write(id_to_gene[e[0]]+ "\t" + id_to_gene[e[1]] + "\t" + str(ed)+"\n")
        print (id_to_gene[e[0]]+ "\t" + id_to_gene[e[1]] + "\t" + str(ed))

    fhout.close()


def compute_cd():
    if args.verbose:
        print ("computing cd values...")

    G = nx.Graph()

    for e in edge_list:
         G.add_edge(e[0], e[1])

    G.to_undirected()
    G.remove_edges_from(G.selfloop_edges())

    N = len(genes)

    weight_out_file = cd_out_file
    fhout = open(weight_out_file, 'w+')


    for i in range(len(edge_list)):
        e = edge_list[i]
        gene1_index = e[0]
        gene2_index = e[1]
        gene1 = id_to_gene[gene1_index]
        gene2 = id_to_gene[gene2_index]

        gene1_patients = []
        gene2_patients = []

        if gene1 in data:
            gene1_patients = data[gene1]
        if gene2 in data:
            gene2_patients = data[gene2]

        gene1_count = len(gene1_patients)
        gene2_count = len(gene2_patients)

        if gene1_count + gene2_count != 0:
            gene1_cover = float(gene1_count) / num_samples
            gene2_cover = float(gene2_count) / num_samples

            cd = gene1_cover * gene2_cover
            fhout.write(id_to_gene[e[0]]+ "\t" + id_to_gene[e[1]] + "\t" + str(cd)+"\n")

            print (id_to_gene[e[0]] + "\t" + id_to_gene[e[1]] + "\t" + str(cd))

    fhout.close()

def find_degree(degree,node):
    for de in degree:
        if de[0] != node:
            continue
        else:
            return de[1]

def gen_pro_list(degree, neib):
    pro_list = []
    node_sum = 0
    for node in neib:
        node_sum += find_degree(degree,node)
    for n in neib:
        n_degree = find_degree(degree,n)
        res = n_degree/node_sum
        pro_list.append(res)
    return (neib[0],pro_list)

def nei_degree_sum(degree,nodes,edge_list, max_degree):
    degree_sum = []
    for node in nodes:
         neib = []
         neib.append(node)
         for i in range(len(edge_list)):
            e = edge_list[i]
            gene1_index = e[0]
            gene2_index = e[1]
            if gene1_index == node:
                neib.append(e[1])
            elif gene2_index == node:
                neib.append(e[0])
            else:
                continue
         print(node)
         degree_sum.append(gen_pro_list(degree,neib))
    degree_sum.sort()
    new_de_sum = []
    for no_de in degree_sum:
        if len(no_de[1]) < max_degree+1:
            for i in range(max_degree + 1 - len(no_de[1])):
                no_de[1].append(0)
        no_de[1].sort()
        no_de[1].reverse()

        new_de_sum.append(no_de)
    print(new_de_sum)
    return new_de_sum

def compute_sim():  # 相似值

    if args.verbose:
        print("computing coverage values...")
    G = nx.Graph()
    for e in edge_list:
        G.add_edge(e[0], e[1])
    weight_out_file = sim_out_file
    fhout = open(weight_out_file, 'w+')

    degree = G.degree()
    nodes = G.nodes

    max_degree=0
    for i in degree:
        if i[1] > max_degree:
            max_degree = i[1]
    new_de_sum = nei_degree_sum(degree, nodes, edge_list, max_degree)

    for i in range(len(edge_list)):
        e = edge_list[i]
        gene1_index = e[0]
        gene2_index = e[1]

        for j in new_de_sum:
            node_num = j[0]
            pro_set = j[1]
            if gene1_index == node_num:
                gene1_num = np.asarray(pro_set)
            if gene2_index == node_num:
                gene2_num = np.array(pro_set)

        M = (gene1_num + gene2_num )/2
        JS = 1- (0.5*st.entropy(gene1_num,M) + 0.5*st.entropy(gene1_num,M))

        fhout.write(id_to_gene[e[0]]  + "\t" + id_to_gene[e[1]]  + "\t" + str(JS) + "\n")
        print(id_to_gene[e[0]] + "\t" + id_to_gene[e[1]] + "\t" + str(JS))
    fhout.close()


def compute_edge_weights(key):
    if args.verbose:
        print ("Assigning edge weights...")

    N = len(genes)

    ed_scores = {}
    with open(path_pre + 'deg_ed.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            ed_scores[line[0]+" "+line[1]] = float(line[2])

    cd_scores = {}
    with open(path_pre + 'deg_cd.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            cd_scores[line[0]+" "+line[1]] = float(line[2])

    sim_scores = {}
    with open(path_pre + 'deg_sim.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            sim_scores[line[0] + " " + line[1]] = float(line[2])

    weight_out_file = edge_out_file
    fhout = open(weight_out_file, 'w+')

    for i in range(len(edge_list)):
        e = edge_list[i]
        gene1_index = e[0]
        gene2_index = e[1]
        gene1 = id_to_gene[gene1_index]
        gene2 = id_to_gene[gene2_index]
        gene1_patients = []
        gene2_patients = []

        if gene1 in data:
            gene1_patients = data[gene1]
        if gene2 in data:
            gene2_patients = data[gene2]

        gene1_count = len(gene1_patients)
        gene2_count = len(gene2_patients)

        if gene1_count + gene2_count != 0:

                if key == 'ecswalk':
                    ed = ed_scores[gene1 + " " + gene2]
                    cd = cd_scores[gene1 + " " + gene2]
                    sim = sim_scores[gene1 + " " + gene2]

                    if ed != 0 and cd != 0 and sim != 0:
                        res = st.hmean((ed , cd)) * sim
                    else: res = 0

                    fhout.write(id_to_gene[e[0]] +"\t" + id_to_gene[e[1]] + "\t" + str(res)+"\n")
                    print (id_to_gene[e[0]] + "\t" + id_to_gene[e[1]] + "\t" + str(res))

    fhout.close()


if __name__ == '__main__':

    # parse arguments
    description = "Compute ed and cd scores and set their product as edge weights (heat)"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-m', "--model", type=str, required=False, default='ecswalk', help="model to assign heat")
    parser.add_argument('-t', '--threshold', type=float, required=False, default=0.7, help="ed_threshold")
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose')

    args = parser.parse_args()

    model = args.model
    threshold = args.threshold

    path_pre = '../out/edge_weights/pre/'
    main_path = '../out/edge_weights/'

    if not os.path.exists(path_pre):
        os.makedirs(path_pre)

    ed_out_file = path_pre + 'deg_ed.txt'
    cd_out_file = path_pre + 'deg_cd.txt'

    sim_out_file = path_pre + 'deg_sim.txt'
    edge_out_file = main_path + model + '.txt'


    num_samples = len(load_patients_to_indices())
    data = load_gene_vs_patient_data()
    genes = load_unique_genes()
    id_to_gene = load_id_to_gene()
    gene_to_id = load_gene_to_id()
    edge_list = load_edge_list()

    compute_ed()
    compute_cd()
    compute_sim()
    compute_edge_weights(model)




