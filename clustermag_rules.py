_author__ = 'brian'
# This file handles taking the input from vasp and calculating the fit.
# The functions that are important for the actual calculating the fit
# have comments

import numpy as np
#import BEG
import clusters
import js
import monomers
import m_structure
import os
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from sklearn.linear_model import Ridge
from sklearn import linear_model


#def write_beg_rules(rule_file):
#    output = open(rule_file, 'w')
#    new_rule = True
#    while new_rule == True:
#        BEG_rule = BEG.BEGObj()
#        BEG_rule.create_rule()
#        output.write('# ' + str(BEG_rule.name) + '\n')
#        output.write(str(BEG_rule.neighbor_order) + '\n')
#        output.write(str(BEG_rule.neighbor_arrangement) + '\n')
#        output.write(str(BEG_rule.home_atom_list) + '\n')
#        output.write(str(BEG_rule.neighbor_atom_list) + '\n')
#        output.write(str(BEG_rule.phase) + '\n')
#        output.write(str(BEG_rule.plane) + '\n')
#        if input('Add another rule? (Y/N):  ') == 'N':
#            new_rule = False
#    output.close()


#def read_beg_rules(rule_file):
#    input_file = open(rule_file, 'r')
#    lines = input_file.readlines()
#    BEG_rule_list = []
#    for i in range(len(lines)):
#        if '#' in lines[i]:
#            BEG_rule = BEG.BEGObj()
#            name = lines[i]
#            BEG_rule.set_name(name.strip('# '))
#            BEG_rule.set_neighbor_order(int(lines[i + 1]))
#            BEG_rule.set_neighbor_arrangement(lines[i + 2].strip())
#            BEG_rule.set_home_atom_list(lines[i + 3].split())
#            BEG_rule.set_neighbor_atom_list(lines[i + 4].split())
#            BEG_rule.set_phase(lines[i + 5].strip())
#            BEG_rule.set_plane(lines[i + 6].strip())
#            BEG_rule.set_composition(lines[i + 7].split())
#            BEG_rule_list.append(BEG_rule)
#    input_file.close()
#    return BEG_rule_list

def write_monomer_rules(rule_file):
    output = open(rule_file, 'w')
    new_rule = True
    while new_rule == True:
        Monomer_rule = monomers.MonomerObj()
        Monomer_rule.create_rule()
        output.write('# ' + str(Monomer_rule.name) + '\n')
        output.write(str(Monomer_rule.neighbor_order) + '\n')
        output.write(str(Monomer_rule.neighbor_arrangement) + '\n')
        output.write(str(Monomer_rule.home_atom_list) + '\n')
        output.write(str(Monomer_rule.neighbor_atom_list) + '\n')
        output.write(str(Monomer_rule.phase) + '\n')
        output.write(str(Monomer_rule.plane) + '\n')
        output.write(str(Monomer_rule.composition) + '\n')
        if input('Add another rule? (Y/N):  ') == 'N':
            new_rule = False
    output.close()


def read_monomer_rules(rule_file):
    input_file = open(rule_file, 'r')
    lines = input_file.readlines()
    Monomer_rule_list = []
    for i in range(len(lines)):
        if '#' in lines[i]:
            rule_name = lines[i]
            rule_name = rule_name.strip('# ')
            rule_neighbor_order = int(lines[i + 1])
            rule_neighbor_arrangement = lines[i + 2].strip()
            rule_home_atom_list = lines[i + 3].split()
            rule_home_atom_list = [int(i) for i in rule_home_atom_list]###### still part of the ugly thing I don't like
            rule_neighbor_atom_list = lines[i + 4].split()
            rule_neighbor_atom_list = [int(i) for i in rule_neighbor_atom_list]###### still part of the ugly thing I don't like
            rule_phase = lines[i + 5].strip()
            rule_plane = lines[i + 6].strip()
            Monomer_rule = clusters.ClusterObj(rule_name,rule_neighbor_order,rule_neighbor_arrangement,rule_home_atom_list,rule_neighbor_atom_list,rule_phase,rule_plane)
            Monomer_rule_list.append(Monomer_rule)
    input_file.close()
    return Monomer_rule_list


def write_cluster_rules(rule_file):
    output = open(rule_file, 'w')
    new_rule = True
    while new_rule == True:
        Cluster_rule = clusters.ClusterObj()
        Cluster_rule.create_rule()
        output.write('# ' + str(Cluster_rule.name) + '\n')
        output.write(str(Cluster_rule.neighbor_order) + '\n')
        output.write(str(Cluster_rule.neighbor_arrangement) + '\n')
        output.write(str(Cluster_rule.home_atom_list) + '\n')
        output.write(str(Cluster_rule.neighbor_atom_list) + '\n')
        output.write(str(Cluster_rule.phase) + '\n')
        output.write(str(Cluster_rule.plane) + '\n')
        output.write(str(Cluster_rule.composition) + '\n')
        if input('Add another rule? (Y/N):  ') == 'N':
            new_rule = False
    output.close()


def read_cluster_rules(rule_file):
    input_file = open(rule_file, 'r')
    lines = input_file.readlines()
    Cluster_rule_list = []
    for i in range(len(lines)):
        if '#' in lines[i]:
            rule_name = lines[i]
            rule_name = rule_name.strip('# ')
            rule_neighbor_order = int(lines[i + 1])
            rule_neighbor_arrangement = lines[i + 2].strip()
            rule_home_atom_list = lines[i + 3].split()
            rule_home_atom_list = [int(i) for i in rule_home_atom_list]###### still part of the ugly thing I don't like
            rule_neighbor_atom_list = lines[i + 4].split()
            rule_neighbor_atom_list = [int(i) for i in rule_neighbor_atom_list]###### still part of the ugly thing I don't like
            rule_phase = lines[i + 5].strip()
            rule_plane = lines[i + 6].strip()
            Cluster_rule = clusters.ClusterObj(rule_name,rule_neighbor_order,rule_neighbor_arrangement,rule_home_atom_list,rule_neighbor_atom_list,rule_phase,rule_plane)
            Cluster_rule_list.append(Cluster_rule)
    input_file.close()
    return Cluster_rule_list


def write_j_rules(rule_file):
    output = open(rule_file, 'w')
    new_rule = True
    while new_rule == True:
        J_rule = js.JObj()
        J_rule.create_rule()
        output.write('# ' + str(J_rule.name) + '\n')
        output.write(str(J_rule.neighbor_order) + '\n')
        output.write(str(J_rule.neighbor_arrangement) + '\n')
        output.write(str(J_rule.home_atom_list) + '\n')
        output.write(str(J_rule.neighbor_atom_list) + '\n')
        output.write(str(J_rule.phase) + '\n')
        output.write(str(J_rule.plane) + '\n')
        if input('Add another rule? (Y/N):  ') == 'N':
            new_rule = False
    output.close()


def read_j_rules(rule_file):
    input_file = open(rule_file, 'r')
    lines = input_file.readlines()
    J_rule_list = []
    for i in range(len(lines)):
        if '#' in lines[i]:
            rule_name = lines[i]
            rule_name = rule_name.strip('# ')
            rule_neighbor_order = int(lines[i + 1])
            rule_neighbor_arrangement = lines[i + 2].strip()
            rule_home_atom_list = lines[i + 3].split()
            rule_home_atom_list = [int(i) for i in rule_home_atom_list]###### still part of the ugly thing I don't like
            rule_neighbor_atom_list = lines[i + 4].split()
            rule_neighbor_atom_list = [int(i) for i in rule_neighbor_atom_list]###### still part of the ugly thing I don't like
            rule_phase = lines[i + 5].strip()
            rule_plane = lines[i + 6].strip()
            J_rule = js.JObj(rule_name,rule_neighbor_order,rule_neighbor_arrangement,rule_home_atom_list,rule_neighbor_atom_list,rule_phase,rule_plane)
            J_rule_list.append(J_rule)
    input_file.close()
    return J_rule_list





