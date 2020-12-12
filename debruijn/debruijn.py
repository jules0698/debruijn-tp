#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""
import argparse
import os
import sys
import random
import statistics
import itertools
import networkx as nx
import matplotlib
from operator import itemgetter

random.seed(9001)




__author__ = "Jules Collat"
__copyright__ = ""
__credits__ = ["Jules Collat"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Jules Collat"
__email__ = "collatjule@eisti.eu"
__status__ = "student"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help = "Fastq file")
    parser.add_argument('-k', dest = 'kmer_size', type=int,
                        default = 21, help = "K-mer size (default 21)")
    parser.add_argument('-o', dest = 'output_file', type=str,
                        default = os.curdir + os.sep + "contigs.fasta",
                        help = "Output contigs in fasta file")
    return parser.parse_args()

def read_fastq(fastq_file):
    with open(fastq_file, "rt") as fastqfile:
        for line in fastqfile:
            yield next(fastqfile).replace("\n", "")
            next(fastqfile)
            next(fastqfile)

def cut_kmer(read, kmer_size):
    for c, valeur in enumerate(read):
        if len(read[c:kmer_size+c]) == kmer_size:
            yield(read[c:kmer_size+c])
        
def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = dict()
    sequences = read_fastq(fastq_file)
    for sequence in sequences:
        for kmer in cut_kmer(sequence, kmer_size):
            if kmer not in kmer_dict.keys():
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] += 1
    return kmer_dict
    
def build_graph(kmer_dict):
    tree = nx.DiGraph()
    for kmer, weight in kmer_dict.items():
        tree.add_edge(kmer[0:len(kmer)-1], kmer[1:len(kmer)], weight=weight)
    return tree

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        for k in range(len(path)):
            if (k != 0 and k != len(path)-1) or (k == 0 and delete_entry_node) or (k == len(path)-1 and delete_sink_node): 
                if path[k] in graph.nodes:
                    graph.remove_node(path[k])
    return graph
            
def std(data):
    return statistics.stdev(data)
    

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    weight_garde = []
    length_garde = []
    path_garde = []
    for k in range(len(weight_avg_list)):
        max_weight=max(weight_avg_list)
        if weight_avg_list[k] == max_weight :
            weight_garde.append(path_list[k])
            length_garde.append(path_length[k])

    for i in range(len(length_garde)):
        max_length = max(length_garde)
        if length_garde[i] == max_length:
            path_garde.append(weight_garde[i])
    bestpath = random.randint(0,len(path_garde)-1)
    path_list.remove(path_garde[bestpath])
    graph_final = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph_final

def path_average_weight(graph, path):
    liste_poids = []
    for k in range(len(path)-1):
        noeud1 = path[k]
        noeud2 = path[k+1]
        liste_poids.append(graph[noeud1][noeud2]['weight'])
    return statistics.mean(liste_poids)

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = []
    path_length = []
    weight_avg_list = []
    delete_entry_node = False
    delete_sink_node = False
    path_list = list(nx.all_simple_paths(graph,ancestor_node,descendant_node))
    for path in path_list:
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph,path))
    graph_final = select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node, delete_sink_node)
    return graph_final

def simplify_bubbles(graph):
    couples_noeuds = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) > 1:
            liste_noeud_predecesseur = list(graph.predecessors(node))
            combinaisons_noeud = itertools.combinations(liste_noeud_predecesseur, 2)
            for combinaison_noeud in combinaisons_noeud:
                noeud_commun = nx.lowest_common_ancestor(graph, combinaison_noeud[0], combinaison_noeud[1], None)
                if noeud_commun is None :
                    continue
                else:
                    couples_noeuds.append((noeud_commun, node))
    for couple_noeud in couples_noeuds:
        if couple_noeud[0] in graph.nodes and couple_noeud[1] in graph.nodes:
            graph = solve_bubble(graph, couple_noeud[0], couple_noeud[1])
    return graph
 
def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    noeud_entree = []
    noeuds = graph.nodes()
    for noeud in noeuds:
        predecessors = graph.predecessors(noeud) 
        if len(list(predecessors)) == 0 :
            noeud_entree.append(noeud)
    return noeud_entree       

    

def get_sink_nodes(graph):
    noeud_sortie = []
    noeuds = graph.nodes()
    for noeud in noeuds:
        successors = graph.successors(noeud) 
        if len(list(successors)) == 0 :
            noeud_sortie.append(noeud)
    return noeud_sortie 

    

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for noeud_entree in starting_nodes:
        for noeud_sortie in ending_nodes:
            chemin = list(nx.all_simple_paths(graph, noeud_entree, noeud_sortie))
            for noeud in chemin:
                chaine = noeud[0]
                for i in range(1,len(noeud)):
                    chaine += noeud[i][-1]
                tuplee = (chaine,len(chaine))
                contigs.append(tuplee)
    return contigs
                    
def save_contigs(contigs_list, output_file):

    def fill(text, width=80):
        return os.linesep.join(text[i:i+width] for i in range (0, len(text), width))   
         
    with open(output_file,"wt") as outputfile:
        nombre_contigs = len(contigs_list)
        for k in range (0, nombre_contigs):
            outputfile.write(">contig_" + str(k) + " " + "len=" + str(contigs_list[k][1]) + "\n" + fill(contigs_list[k][0], width=80) + "\n")
    


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    fastq_file = args.fastq_file
    kmer_size = args.kmer_size
    output_file = args.output_file
    dictionnaire = build_kmer_dict(fastq_file, kmer_size)

    graphe = build_graph(dictionnaire)
    graphe_sans_bulle = simplify_bubbles(graphe)

    starting_nodes = get_starting_nodes(graphe_sans_bulle)
    sink_nodes = get_sink_nodes(graphe_sans_bulle)

    contigs = get_contigs(graphe_sans_bulle, starting_nodes, sink_nodes)
    save_contigs(contigs, output_file)

if __name__ == '__main__':
    
    
    main()
