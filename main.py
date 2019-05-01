__author__ = 'brian'
import compile_vasp_structures as cvs
import postprocess_vasp as ppv
import clustermag_rules as cmr
import calc_fitting_params as cfp
# import mc_functions as mc
import sys

aust_tol = 0.02
spin_style = ['threshold', 'threshold',
              'threshold']  # options for spin_tol. Assuming [Ni Mn In]. choose 'threshold' or 'factor'
spin_tol = [0.05, 2, 0]  # .05                          # insert spin parameters here, this assumes [Ni Mn In ]
species = ['Ni', 'Mn', 'In']  # this is the order that the post-processed data is reported, NEEDS TO BE Heusler format Ni2MnIn, Ni2FeGa.

root_dir = './NiMnIn_Vasp_Data'  # where the VASP directories are  #root_dir = '/Users/brian/Desktop/folder'
vasp_data_file = './NiMnIn_Data'  # generated in compile_vasp_structures>import_vasp, summarizes output of all VASP calculations
vasp_data_file_pp = './NiMnIn_Data_pp'  # post-processed version of VASP results with spins, positions selected
cluster_file = "./Cluster_Rules"  # cluster expansion rules
j_file = './J_Rules'  # heisenberg rules
monomer_file = './Monomer_Rules'
# fitting_structures_file = './'

# Determine what needs to be generated from scratch
vasp_data_exists = True  # should be False if I want to regenrate files
vasp_pp_exists = True  # postprocessing of VASP results or no?
Cluster_rules_exist = True  # define cluster rules
J_rules_exist = True  # define heisenberg rules                                  # results of fitting model
Fitting_params_exist = False

# write Cluster and J rules file if doesn't exist, then read the rules
if Cluster_rules_exist is False:  # writes cluster rules if doesn't already exist
    cmr.write_cluster_rules(cluster_file)
if J_rules_exist is False:  # writes j_rules if doesn't already exist
    cmr.write_j_rules(j_file)
Cluster_rules = cmr.read_cluster_rules(cluster_file)
J_rules = cmr.read_j_rules(j_file)
Monomer_rules = cmr.read_monomer_rules(monomer_file)
if Fitting_params_exist == False:
    # summarize VASP data
    if vasp_data_exists is False:  # will make summary of VASP results if it doesn't already exist
        cvs.import_vasp(root_dir, vasp_data_file, species)

    # Read the VASP datafile and initialize a structure object for each
    # postprocess according to user selected parameters above and the cluster and j rules
    # calculation of sums and checking for duplicates occurs in here now
    # if a given structure is considered a duplicate then it is not added to the structure_list
    M_structures = ppv.generate_m_structure(vasp_data_file, len(Monomer_rules), len(Cluster_rules), len(J_rules), aust_tol, spin_style, spin_tol, Monomer_rules, Cluster_rules, J_rules)
    ppv.write_structures_processedvasp(M_structures, vasp_data_file_pp)
    ppv.summarize_classification(M_structures)
    warning_threshold = 0.5
    ppv.summarize_fitting_structures(M_structures,warning_threshold)

    ## Ridge Regression Fitting with Regularization
    ##ppv.scale_enrg(M_structures)
    ####################################################################################################
    # Single fit
    ####################################################################################################
    #ppv.set_KJ_ratio(M_structures, Monomer_rules, Cluster_rules, J_rules, .35)
    Js, intercept = cfp.single_fit(M_structures, False, Monomer_rules, Cluster_rules, J_rules)
    print(intercept)
    cfp.write_fitting_parameters(M_structures, Monomer_rules, Cluster_rules, J_rules, Js, intercept, 200)
    cfp.plot_data3(M_structures, Monomer_rules, Cluster_rules, J_rules, Js, intercept, 200, False)
    ####################################################################################################
    # Double fit
    ####################################################################################################
    # Js1, intercept1 = cfp.first_fit(M_structures, True, Monomer_rules)
    # print(intercept1)
    # cfp.linearize(M_structures,"MONOMER",Monomer_rules,Js1,intercept1)
    # Js2, intercept2 = cfp.second_fit(M_structures, False, Cluster_rules, J_rules)
    # Js = []
    # for i in range(len(Js1)):
    #     Js.append(Js1[i])
    # for i in range(len(Js2)):
    #     Js.append(Js2[i])
    # cfp.write_fitting_parameters(M_structures, Monomer_rules, Cluster_rules, J_rules, Js, intercept1, 200)
    # cfp.plot_data3(M_structures, Monomer_rules, Cluster_rules, J_rules, Js, intercept1, 200, False)
    # cfp.plot_data4(M_structures, Cluster_rules, J_rules, Js, intercept1, 200, False)
    ####################################################################################################
    # Triple fit
    ####################################################################################################
    # Js1, intercept1 = cfp.first_fit(M_structures, True, Monomer_rules)
    # print(intercept1)
    # cfp.linearize(M_structures, "MONOMER", Monomer_rules, Js1, intercept1)
    # Js2, intercept2 = cfp.beg_fit(M_structures, False)
    # cfp.linearize(M_structures, "CLUSTER", Cluster_rules, Js2, intercept2)
    # Js3, intercept3 = cfp.spin_fit(M_structures, False)
    # Js = []
    # for i in range(len(Js1)):
    #     Js.append(Js1[i])
    # for i in range(len(Js2)):
    #     Js.append(Js2[i])
    # for i in range(len(Js3)):
    #     Js.append(Js3[i])
    # cfp.write_fitting_parameters(M_structures, Monomer_rules, Cluster_rules, J_rules, Js, intercept1+intercept2+intercept3, 200)
    # cfp.plot_data3(M_structures, Monomer_rules, Cluster_rules, J_rules, Js, intercept1+intercept2+intercept3, 200, False)
    # cfp.plot_data4(M_structures, Cluster_rules, J_rules, Js, intercept1+intercept2+intercept3, 200, True)

else:
    paramiters_file = open('FittingParameters', 'r')
    lines = paramiters_file.readlines()
    reading_params = True
    line_index = 1
    Js = []
    while reading_params == True:
        line = lines[line_index]
        if 'Actual Energy' in line or line.split() == []:
            reading_params = False
        else:
            line = line.split()
            Js.append(float(line[2]))
            line_index += 1
