

import networkx as nx

from src.utils import *
import numpy as np
import scipy.sparse as sp
import operator
import os
from collections import deque
import glob
import argparse


def find_largest_comp(list_comp_graph):
    large = list_comp_graph[0]
    large_index = 0
    for i in range(len(list_comp_graph)):
        if len(large) < len(list_comp_graph[i]):
            large = list_comp_graph[i]
            large_index = i

    return (large,large_index)


def list_graph_to_list_comp(list_graph):
    list_comp = []
    for i in range (len(list_graph)):
        list_comp.extend([list_graph[i].nodes])

    return list_comp


def max_comp_length(list_graph):
    max_length = len(list_graph[0])

    for i in range (len(list_graph)):
        if len(list_graph[i])> max_length:
            max_length = len(list_graph[i])

    return max_length

def find_max_outdegree(comp):
    max_out_degree = 0
    largest_node_id = 0
    largest_gene_id = ''

    for n in comp.nodes():
        out_n = comp.out_degree(n)
        if max_out_degree < out_n:
            max_out_degree = out_n
            largest_node_id = n
            largest_gene_id = id_to_gene[n]

    return (max_out_degree, largest_node_id, largest_gene_id)


def star_construction(large_comp, largest_node_id):
    all_neighbor_of_largest = set()
    star_nodes_set = set([largest_node_id])

    for e in large_comp.edges:

        if largest_node_id in e:

            all_neighbor_of_largest.update(e)

    all_neighbor_of_largest.remove(largest_node_id)

    for u in all_neighbor_of_largest:
        in_comp = True
        for e in large_comp.edges:
            if e[0] == u:
                v = e[1]

                if v != largest_node_id and (v not in all_neighbor_of_largest):

                    in_comp = False

                    break
            elif e[1] == u:
                v = e[0]

                if v != largest_node_id and (v not in all_neighbor_of_largest):
                    in_comp = False

                    break

        if in_comp == True:
            star_nodes_set.update([u])

    return list(star_nodes_set)


def split_large_components(list_comp_graph, large_threshold = 0):

    list_comp = list_graph_to_list_comp(list_comp_graph)
    set_all_genes_before = set([n for comp in list_comp for n in comp])

    threshold_small_comp = args.k_threshold


    gene_set = set()
    for i in range(len(list_comp)):
        gene_set = gene_set.union(set(list_comp[i]))

    max_out_degree_all = 0

    for comp in list_comp_graph:
         max_out_degree_pre, largest_node_id_pre, largest_gene_id_pre = find_max_outdegree(comp)

         if max_out_degree_all < max_out_degree_pre:
               max_out_degree_all = max_out_degree_pre


    if large_threshold == 0:
        threshold_large_comp = max_out_degree_all
    else:
        threshold_large_comp = large_threshold


    list_graph_leftover = []
    list_large_components = []
    for comp in list_comp_graph:
        if len(comp)>= threshold_large_comp:
            list_large_components.append(comp)

        else:
            list_graph_leftover.append(comp)


    all_modified_component_list = []

    for lc in list_large_components:
        main_comp_list = []
        small_comp_list = []
        #large_graph_queue
        large_comp_queue = deque()
        large_comp_queue.append(lc)
        while len(large_comp_queue) >0:

            largest_comp = large_comp_queue.popleft()
            max_out_degree, largest_node_id, largest_gene_id = find_max_outdegree(largest_comp)
            reduced_comps = largest_comp.copy()

            removable_nodes_list = star_construction(largest_comp,largest_node_id)


            reduced_comps.remove_nodes_from(removable_nodes_list)

            largest_gene_graph = largest_comp.subgraph(removable_nodes_list).copy()
            if len(largest_gene_graph) < threshold_small_comp:
                small_comp_list.append(largest_gene_graph)
            else:
                main_comp_list.append(largest_gene_graph)

            scc = nx.strongly_connected_components(reduced_comps)

            for comp in scc:

                if len(comp)> threshold_large_comp:
                    g = reduced_comps.subgraph(comp).copy()
                    large_comp_queue.append(g)

                elif len(comp)< threshold_small_comp:
                    g = reduced_comps.subgraph(comp).copy()
                    small_comp_list.append(g)

                else:
                    g = reduced_comps.subgraph(comp).copy()
                    main_comp_list.append(g)


        temp_list = []
        for scm in range(len(main_comp_list)):

            if len(main_comp_list[scm]) < threshold_small_comp:
                small_comp_list.append(main_comp_list[scm])

            else:
                temp_list.append(main_comp_list[scm])

        main_comp_list = temp_list[:]


        with open("../out/edge_weights/ecswalk.txt")as ecswalk:
            lines = ecswalk.readlines()

        for scomp in range(len(small_comp_list)):

            for gene_s in small_comp_list[scomp]:
                max_weight_ave = 0
                max_weight_flag = 100000

                for comp_index in range(len(main_comp_list)):
                    big_weight_sum = 0
                    big_edge_sum = 0
                    for gene_i in main_comp_list[comp_index]:

                        for line in lines:
                            e1, e2, weight = line.split()
                            if (gene_s, gene_i) == (e1, e2):
                                big_weight_sum += weight
                                big_edge_sum += 1
                            if (gene_i, gene_s) == (e1, e2):
                                big_weight_sum += weight
                                big_edge_sum += 1
                    if big_edge_sum == 0:
                        continue
                    big_weight_ave = big_weight_sum / big_edge_sum

                    if big_weight_ave > max_weight_ave:
                        max_weight_ave = big_weight_ave
                        max_weight_flag = comp_index
                small_weight_sum = 0
                small_edge_sum = 0
                for scomp in range(len(small_comp_list)):
                    for gene_n in small_comp_list[scomp] :
                        if gene_s == gene_n:
                            continue
                        for line in lines:
                            e1, e2, weight = line.split()
                            if (gene_s, gene_n) == (e1, e2):
                                big_weight_sum += weight
                                small_edge_sum += 1
                            if (gene_n, gene_s) == (e1, e2):
                                big_weight_sum += weight
                                small_edge_sum += 1
                if small_edge_sum ==0:
                    small_edge_sum = 1
                small_weight_ave = small_weight_sum / small_edge_sum

                if max_weight_flag == 100000:
                    continue
                if max_weight_ave >= small_weight_ave:
                    main_comp_list[max_weight_flag].add_nodes_from(gene_s)
                    temp_subgraph_nodes = main_comp_list[max_weight_flag].nodes
                    main_comp_list[max_weight_flag] = lc.subgraph(temp_subgraph_nodes).copy()


        all_modified_component_list.extend(main_comp_list[:])

    return list_graph_leftover[:] + all_modified_component_list[:]


if __name__ == '__main__':


    description = "Find connected components then isolate and split large components"#查找连接的组件，然后隔离并拆分大型组件
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-m', "--model", type=str, required=False, default='ecswalk',
                        help="model to find connected components")
    parser.add_argument('-k', "--k_threshold", type=int, required=False, default=3,
                        help="min length of module")
    parser.add_argument('-ns', '--numstart', type=int, required=False, default=1000, help="max number of genes")#最大基因数量
    parser.add_argument('-ne', '--numend', type=int, required=False, default=100, help="min number of genes")#最小基因数量
    parser.add_argument('-step', '--stepsize', type=int, required=False, default=100, \
                        help="step size of decrement from num start to num end")
    parser.add_argument('-ts', '--start_threshold', type=float, required=False, default=0.0002, help="starting threshold")
    parser.add_argument('-cc', '--cc_original', action='store_true', \
                        default=False, help="generate connected component files pre-split")
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose',default=True)

    args = parser.parse_args()

    id_to_gene = load_id_to_gene()
    gene_to_id = load_gene_to_id()

    model = args.model

    filepath = "../out/connected_components_original/" + model + "/"
    cc_path = "../out/connected_components_isolarge/" + model + "/"

    try:

        threshold_start = float(glob.glob(filepath+'cc_n{}_*'.format(args.numstart))[0].split('/')[-1].split('d')[1][:-4])-0.000001
    except:
        threshold_start = args.start_threshold

    path = cc_path

    if not os.path.exists(cc_path):
        os.makedirs(cc_path)
    if not os.path.exists(filepath) and args.cc_original:
        os.makedirs(filepath)

    LARGE_NUM = 100000
    k = args.k_threshold
    our_E = sp.load_npz("../out/random_walk/"+ model+"_sparse_matrix_e.npz")
    E = our_E.toarray()

    num_start = args.numstart
    num_end = args.numend


    N = len(id_to_gene)
    G = nx.DiGraph()

    for i in range(N):
        G.add_node(i)

    for i in range(N):
        for j in range(N):
            if E[i][j] > threshold_start:
                G.add_edge(j, i)


    list_graphs = []
    num_nodes = 0
    count_comp = 0
    smallest_weight = LARGE_NUM
    smallest_edge_str = ''
    for comp in nx.strongly_connected_components(G):
        if len(comp) >= k:

            subG =G.subgraph(comp).copy()
            list_graphs.append(subG)
            num_nodes += len(comp)

            for e in subG.edges():
                u = e[0]#6636
                v = e[1]#3170
                key = str(count_comp) + '_' + str(u) + '_' + str(v)
                if E[v][u] < smallest_weight:
                    smallest_edge_str = key
                    smallest_weight = E[v][u]
            count_comp += 1
            print(count_comp)
            if args.verbose:
                print (smallest_edge_str)

    if args.verbose:
        print (path)
        print ('number of nodes in the beginning is : ', num_nodes)
        print (smallest_edge_str)

    while num_start >= num_end:
        print("a" + str(num_start))
        iteration = 0
        threshold_found = 0
        num_target_genes = num_start
        while num_nodes > num_target_genes:
            if args.verbose:
                print(str(num_target_genes) + " " + str(num_nodes) + " " + str(iteration))

            smallest_edge = smallest_edge_str.split('_')
            threshold_found = smallest_weight
            graph_index = int(smallest_edge[0])
            subG = list_graphs[graph_index]
            num_nodes -= len(subG.nodes())


            subG.remove_edge(int(smallest_edge[1]), int(smallest_edge[2]))
            subG_comps = nx.strongly_connected_components(subG)
            del list_graphs[graph_index]
            for comp in subG_comps:
                if len(comp) >= k:
                    list_graphs.append(subG.subgraph(comp).copy())
                    num_nodes += len(comp)

            num_comps = len(list_graphs)

            smallest_weight = LARGE_NUM
            smallest_edge_str = ''
            for i in range(num_comps):
                graph = list_graphs[i]

                if len(graph.nodes()) >= k:
                    for e in graph.edges():
                        u = e[0]
                        v = e[1]
                        if E[v][u] < smallest_weight:
                            key = str(i) + '_' + str(u) + '_' + str(v)
                            smallest_edge_str = key
                            smallest_weight = E[v][u]
                            print(smallest_edge_str)
            iteration += 1

        if args.cc_original:

            outfilename = filepath + "cc_n" + str(num_target_genes) + "_d" + str(threshold_found) + ".txt"

            with open(outfilename, "w") as f:
                for i in range(len(list_graphs)):
                    for node_index in list_graphs[i]:
                        f.write(id_to_gene[node_index] + " ")
                    f.write("\n")
                f.close()


        largest_comp, largest_comp_index = find_largest_comp(list_graphs)
        component_list = list_graphs[:]
        largeset = set([id_to_gene[idx] for idx in largest_comp])

        largest_gene_comp_list = []

        final_comp_list_graphs = split_large_components(component_list)


        outfilename = path + "cc_n" + str(num_target_genes) + "_"+ str(len(largest_comp))+ "_"+\
                      str(len(final_comp_list_graphs)) + "_d" + str(threshold_found) + ".txt"


        with open(outfilename, "w") as f:
            for i in range(len(final_comp_list_graphs)):
                for node_index in final_comp_list_graphs[i]:
                    f.write(id_to_gene[node_index] + " ")
                f.write("\n")
            f.close()

        num_start -= args.stepsize

