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
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Jules Collat"
__copyright__ = "Universite Paris Diderot"
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
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

def read_fastq(fastq_file):
    with open(fastq_file,"rt") as fastqfile:
        for line in fastqfile:
            yield next(fastqfile).replace("\n","")
            next(fastqfile)
            next(fastqfile)
            
def cut_kmer(read, kmer_size):
    for c, valeur in enumerate(read):
        if len(read[c:kmer_size+c]) == kmer_size:
            yield (read[c:kmer_size+c])
        
def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = dict()
    sequences = read_fastq(fastq_file)
    for sequence in sequences:
        for kmer in cut_kmer(sequence, kmer_size):
            if kmer not in kmer_dict.keys():
                kmer_dict[kmer]=1
            else:
                kmer_dict[kmer]+=1
    return kmer_dict
    
def build_graph(kmer_dict):
    tree = nx.DiGraph()
    for kmer, weight in kmer_dict.items():
        tree.add_edge(kmer[0:len(kmer)-1], kmer[1:len(kmer)], weight=weight)
    return tree

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    return statistics.stdev(data)
    

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    liste_poids = []
    for k in range(len(path)-1):
        noeud1 = path[k]
        noeud2 = path[k+1]
        liste_poids.append(graph[noeud1][noeud2]['weight'])
    return statistics.mean(liste_poids)


    

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    noeud_entree =[]
    noeuds = graph.nodes()
    for noeud in noeuds:
        predecessors=graph.predecessors(noeud) 
        if len(list(predecessors)) == 0 :
            noeud_entree.append(noeud)
    return noeud_entree       

    

def get_sink_nodes(graph):
    noeud_sortie =[]
    noeuds = graph.nodes()
    for noeud in noeuds:
        successors=graph.successors(noeud) 
        if len(list(successors)) == 0 :
            noeud_sortie.append(noeud)
    return noeud_sortie 

    

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for noeud_entree in starting_nodes:
        for noeud_sortie in ending_nodes:
            chemin=list(nx.all_simple_paths(graph,noeud_entree,noeud_sortie))
            for noeud in chemin:
                chaine=noeud[0]
                for i in range(1,len(noeud)):
                    chaine+= noeud[i][-1]
                tuplee= (chaine,len(chaine))
                contigs.append(tuplee)
    return contigs
                    
    



def save_contigs(contigs_list, output_file):

    def fill(text, width=80):
        return os.linesep.join(text[i:i+width] for i in range (0, len(text), width))   
         
    with open(output_file,"wt") as outputfile:
        nombre_contigs = len(contigs_list)
        for k in range (0,nombre_contigs):
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

if __name__ == '__main__':
    
    
    main()
