# -*- coding: utf-8 -*-
from cobra import Reaction
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
import itertools
import pandas as pd
import os
# reduction degree of elements
rd_element = {'C': 4, 'H': 1, 'O': -2, 'N': -3,
              'S': 6, 'P': 5, 'Fe': 3, 'Ca': 2, 'Na': 1, 'Mg': 2}


def model_preprocess(model, specialgroup, specialgroup_begin, compartments):
    # preprosess the model, modify exchange reactions, add demand reactions etc
    newmodel = model.copy()
    remove_input(newmodel)  # find and remove carbon input in the model
    # allow atp/ADP exchange without limit, pyr-PEP get pathway with one reaction
    newmodel.reactions.get_by_id('ATPM').bounds = (-1000.0, 1000.0)
    special_demand(newmodel, specialgroup, specialgroup_begin, compartments)
    add_demandall(newmodel)
    return newmodel


def remove_input(model):
    # check and remove carbon source input in the model
    # return the updated model
    for met in model.metabolites:
        me = met.elements
        metid = met.id
        if 'C' in me:
            if 'EX_' + metid in model.reactions:
                ex_reaction = model.reactions.get_by_id('EX_' + metid)
                model.reactions.get_by_id(
                    ex_reaction.id).bounds = (0.0, 1000.0)
            elif 'SK_' + metid in model.reactions:
                sk_reaction = model.reactions.get_by_id('SK_' + metid)
                model.reactions.get_by_id(
                    sk_reaction.id).bounds = (0.0, 1000.0)
            elif 'sink_' + metid in model.reactions:
                sk_reaction = model.reactions.get_by_id('sink_' + metid)
                model.reactions.get_by_id(
                    sk_reaction.id).bounds = (0.0, 1000.0)
            elif 'DM_' + metid in model.reactions:
                dm_reaction = model.reactions.get_by_id('DM_' + metid)
                model.reactions.get_by_id(
                    dm_reaction.id).bounds = (0.0, 1000.0)
    return model
    # r=model.reactions.get_by_id("EX_co2_e")# CO2 fixation may be dealed seperately
    #r.bounds=(0, 1000)


def get_compartments(model):
    comps = []
    for metbolite in model.metabolites:
        comp = metbolite.compartment
        comps.append(comp)
    compartments = list(set(comps))
    return compartments


def special_demand(model, specialgroup, specialgroup_begin, compartments):
    # add demand reaction for special groups such as CoA,ACP, but not those containing these groups, to make sure metabolites containing
    # such groups can be consumed
    #    specialgroup: a list of group names often at the end of the metabolite ID such as coa
    #    specialgroup_begin: a list of group names often at the beginning of the metabolite ID, such as UDP uac for UDP-N-acetyl
    #    compartments: a list of compartments to add demand reactions 'c','p','e'
    for comp in compartments:
        for m in specialgroup:
            cid = m+'_'+comp
            if cid in model.metabolites and 'DM_' + cid not in model.reactions:
                # add demand reaction for CoA so that acyl-CoA can be consumed
                g = model.metabolites.get_by_id(cid)
                model.add_boundary(g, type='demand')
        for m in specialgroup_begin:
            cid = m+'_'+comp
            if cid in model.metabolites and 'DM_' + cid not in model.reactions:
                # add demand reaction for CoA so that acyl-CoA can be consumed
                g = model.metabolites.get_by_id(cid)
                model.add_boundary(g, type='demand')
    return model


def reduction_degree(met, specialgroup_begin):
    # calculate the reduction degree for a metabolite object, special process for coa, ACP containing molecules
    me = met.elements
    metid = met.id
    rd = 0
    for element, count in me.items():
        if element in rd_element:
            rd += count * rd_element[element]
    rd -= met.charge  # need to correct with charge
    # for CoA containing metabolites only consider the main part
    if "coa_" in metid and not metid[:-2] in ["coa", 'dpcoa']:
        rd -= 88
    # for ACP containing metabolites only consider the main part
    elif "ACP_" in metid and metid[0:3] != "ACP":
        rd -= 57
    # for THF containing metabolites only consider the main part, may be not necessary as the
    elif "thf_" in metid and metid[0:3] != "thf" and metid[0:6] != "mththf":
        # the non THF part only have one carbon
        rd -= 66
    # startin with udp, cdp, but not udp_c udp_p etc
    elif metid[0:3] in specialgroup_begin and metid[3] != "_" and metid[0:5] != "gdptp":
        if metid[0] == 'g':  # for gdp
            rd -= 28
        else:
            rd -= 30
    return rd


def cnumber(met, specialgroup_begin):
    # calculate the carbon atom numbers for a metabolite object, special process for coa, ACP containing molecules
    me = met.elements
    metid = met.id
    if 'C' in me:
        nc = me["C"]
        # for CoA containing metabolites only consider the main part
        if "coa_" in metid and not metid[:-2] in ["coa", 'dpcoa']:
            nc -= 21
        # for ACP containing metabolites only consider the main part
        elif "ACP_" in metid and metid[0:3] != "ACP":
            nc -= 11
        # for THF containing metabolites only consider the main part, may be not necessary as the
        elif "thf_" in metid and metid[0:3] != "thf" and metid[0:6] != "mththf":
            # the non THF part only have one carbon
            nc -= 19
        # startin with udp, cdp, but not udp_c udp_p etc
        elif metid[0:3] in specialgroup_begin and metid[3] != "_" and metid[0:5] != "gdptp":
            if metid[0] == 'a' or metid[0] == 'g':
                nc -= 10
            else:
                nc -= 9
    return nc


def add_demand_old(model, met, specialgroup_begin):
    # add a demand reaction for a metabolite to model, met should be the metabolite object
    # return the demand reaction
    metid = met.id
    # for dpCoA still add normal demand reaction
    if "coa_" in metid and not metid[:-2] in ["coa", 'dpcoa']:
        demand = Reaction('DM_' + metid)
        model.add_reaction(demand)
        # metid[-1] for the same compartment
        demand.build_reaction_from_string(metid+'--> coa_c', fwd_arrow="-->")
    # for CoA containing metabolites only consider the main part
    elif "ACP_" in metid and metid[0:3] != "ACP":
        demand = Reaction('DM_' + metid)
        model.add_reaction(demand)
        demand.build_reaction_from_string(metid+'--> ACP_c', fwd_arrow="-->")
    # mththf is different
    elif "thf_" in metid and metid[0:3] != "thf" and metid[0:6] != "mththf":
        demand = Reaction('DM_' + metid)
        model.add_reaction(demand)
        demand.build_reaction_from_string(metid+'--> thf_c', fwd_arrow="-->")
    # for cdp containing metabolites but not cdp_c
    elif metid[0:3] in specialgroup_begin and metid[3] != "_" and metid[0:5] != "gdptp":
        demand = Reaction('DM_' + metid)
        model.add_reaction(demand)
        if metid[0:3] == 'uac':  # for uac reduce udp
            demand.build_reaction_from_string(
                metid+'--> udp_c', fwd_arrow="-->")
        else:
            demand.build_reaction_from_string(
                metid+'--> '+metid[0:3]+'_c', fwd_arrow="-->")
    else:
        #dm_met = newmodel.metabolites.get_by_id(met)
        demand = model.add_boundary(met, type='demand')
    return demand


def add_demandall(model):
    # add demand reactions for all carbon containing metabolites,to allow coproducts to be produced in a pathway
    for met in model.metabolites:
        me = met.elements
        if 'C' in me:
            if 'DM_' + met.id not in model.reactions:  # exclude those demand reaction already added
                # demand=add_demand(model,met) #not use this one as udpM-->udp then udp but not M can be used as substrate for pyr production
                demand = model.add_boundary(met, type='demand')
    return model


def path(m, s, p, pathfolder, specialgroup, specialgroup_begin):
    # calculate pathways between two metabolites in a model
    # m: model, s: substrate ID, p: product ID,pathfolder: the folder to save the calculated pathways
    with m as model:  # to avoid change of the original model
        substrate = model.metabolites.get_by_id(s)
        product = model.metabolites.get_by_id(p)
        if 'SK_' + s in model.reactions:  # add sink reaction
            model.reactions.get_by_id('SK_' + s).bounds = (-10, 0)
            #rea=model.reactions.get_by_id('SK_' + s)
            # rea.lb=-10.0
            # rea.ub=-10.0
        else:
            model.add_boundary(substrate, type='sink')
            model.reactions.get_by_id('SK_'+s).bounds = (-10, 0)
            #model.add_boundary(substrate, type='sink', lb=10)
        if 'DM_' + p in model.reactions:
            demand = model.reactions.get_by_id('DM_' + p)
            # for special metabolites, change to a new demand reaction to produce acCOA, acACP etc
            if p[-5:-2] in specialgroup or p[0:3] in specialgroup_begin:
                demand.remove_from_model()
                # model.remove_reactions([demand])
                demand = add_demand(model, product, specialgroup_begin)
        else:
            demand = add_demand(model, product, specialgroup_begin)
        model.objective = demand
        try:
            fluxes = pfba(model).fluxes
        except:
            metpair = [s, p, 'no solution']
        else:

            # solution = pfba(model,fraction_of_optimum=0.5) #using fraction_of_optimum to get the backbone pathway without carbon fixation
            # but it reduces the rate of the input reaction even it is fixed at -10
            # not work for AcCoA to pyr, as other compounds in reaction coa_c + 2.0 flxso_c + pyr_c <=> accoa_c + co2_c + 2.0 flxr_c
            # still need to be balanced
            # rate=solution.objective_value #not use this because it is the sum of the fluxes
            rate = model.reactions.get_by_id('DM_' + p).flux
            # sometimes the input flux is not equal to the upper bound
            influx = abs(model.reactions.get_by_id('SK_' + s).flux)
            # print(rate)
            if rate > 1e-6:
                # pathfile = open(pathfolder+ "/"+  s + "--" + p + ".txt", 'w')  # save the pathway for checking
                # for r,v in fluxes.iteritems():
                #    if abs(v)>1e-6:
                #        pathfile.write(r + '\t' + model.reactions.get_by_id(r).build_reaction_string() + '\t' + str(v) + '\n')
                # pathfile.close()
                sc = cnumber(substrate, specialgroup_begin)
                srd = reduction_degree(substrate, specialgroup_begin)
                pc = cnumber(product, specialgroup_begin)
                prd = reduction_degree(product, specialgroup_begin)
                cyield = rate*pc/sc/influx  # calculate the carbon yield
                if prd == 0:
                    cyieldrd = 0
                else:
                    cyieldrd = srd*pc/sc/prd
                metpair = [s, substrate.formula, sc, srd, srd/sc, p, product.formula,
                           pc, prd, prd/pc, rate, cyield, cyieldrd, cyieldrd-cyield]
            else:
                metpair = "no path"

    return metpair


def get_metpair(rea, pi_pairs1, h_pairs1, pi_pairs2, h_pairs2, nh4_pairs, other_pairs, currency_mets):
    # get the metabolite links for a reaction, excluding links through currency metabolites
    # processing in the order of P & H transfer, N transfer and other transfers
    sub_pro = []
    mark = 0  # mark if there is currency metabolite pairs for P/H transfer in the reaction
    c2mark = 0  # mark if there is currency metabolite pairs for P/H transfer in the reaction, second batch
    nmark = 0  # mark if there is nh4 transfer currency metabolite pairs in the reaction, deal seperately
    omark = 0  # mark if there is other group transfer currency metabolite pairs in the reaction, deal seperately
    cmet = []  # temorary currency metabolite list
    c2met = []  # temorary currency metabolite list, second batch
    ncmet = []  # temorary currency metabolite list for N transfer
    ocmet = []  # temorary currency metabolite list for other pair transfer
    phmet = []  # metabolite for P/H transfer, to be excluded again in the last step if >1 pairs still remaining, especially for UDP ADPsugars
    ex_pairs1 = pi_pairs1+h_pairs1  # excluded currency metabolite pairs, first batch
    ex_pairs2 = pi_pairs2+h_pairs2  # excluded currency metabolite pairs, second batch
    for sp in ex_pairs1:
        phmet.append(sp[0])
        phmet.append(sp[1])
    for sp in ex_pairs2:
        phmet.append(sp[0])
        phmet.append(sp[1])
    phmet = list(set(phmet))  # remove repeats
    # also check if CoA is a currency metabolite in the last step
    phmet += ['coa_c', 'coa_p', 'coa_e']
    subs = [
        m.id for m in rea.reactants if 'C' in m.elements if m.id not in currency_mets]
    pros = [m.id for m in rea.products if 'C' in m.elements if m.id not in currency_mets]
    for s, p in itertools.product(subs, pros):
        if (s, p) in ex_pairs1 or (p, s) in ex_pairs1:  # need to consider direction
            mark = 1
            cmet.append(s)
            cmet.append(p)
        if (s, p) in ex_pairs2 or (p, s) in ex_pairs2:  # need to consider direction
            c2mark = 1
            c2met.append(s)
            c2met.append(p)
        if (s, p) in nh4_pairs or (p, s) in nh4_pairs:  # for nh4 transfer
            nmark = 1
            ncmet.append(s)
            ncmet.append(p)
        if (s, p) in other_pairs or (p, s) in other_pairs:  # for other pair transfer
            omark = 1
            ocmet.append(s)
            ocmet.append(p)
    # if len(sub_pro)>1 and mark==1:
    if mark == 1:  # process in order
        subs = [m for m in subs if m not in cmet]
        pros = [m for m in pros if m not in cmet]
    if c2mark == 1:  # process in order
        subsn = [m for m in subs if m not in c2met]
        prosn = [m for m in pros if m not in c2met]
        if subsn and prosn:  # proceed if only there are still other metabolites in the reactant and product list
            subs = subsn
            pros = prosn
    if nmark == 1:
        subsn = [m for m in subs if m not in ncmet]
        prosn = [m for m in pros if m not in ncmet]
        if subsn and prosn:  # proceed if only there are still other metabolites in the reactant and product list
            subs = subsn
            pros = prosn
    if omark == 1:
        subsn = [m for m in subs if m not in ocmet]
        prosn = [m for m in pros if m not in ocmet]
        if subsn and prosn:  # proceed if only there are still other metabolites in the reactant and product list
            subs = subsn
            pros = prosn
    if len(subs) > 1:  # to remove UTP, UDP etc.
        subsn = subs
        for m in subsn:
            if m in phmet:
                subs.remove(m)
                if len(subs) == 1:
                    break
    if len(pros) > 1:  # to remove UTP, UDP etc.
        prosn = pros
        for m in prosn:
            if m in phmet:
                pros.remove(m)
                if len(pros) == 1:
                    break
    for s, p in itertools.product(subs, pros):
        sub_pro.append((s, p))
    return sub_pro


def bow_tie_structure_for_model(modelfile, specialgroup, specialgroup_begin, euk_model_list, result_save_path):
    # calculate pathways between a precusor and other metabolites and obtain the bow tie structure
    # process the model to get a newmodel
    # model_name = eachmodelfile.split('\\')[1].split('.')[0] #for windows
    model_name = modelfile.split('/')[-1].split('.')[0]
    print(model_name, ' ')
    model = read_sbml_model(modelfile)
    compartments = get_compartments(model)
    newmodel = model_preprocess(
        model, specialgroup, specialgroup_begin, compartments)
    if not os.path.exists(result_save_path):
        os.makedirs(result_save_path)
    pathout = []  # list of list to save pathway information from the precusor
    pathin = []  # list of list to save pathway information to the precusor
    metall = []  # list for all C containing metabolites to calculate the IS subset
    bowtie = []
    pathexcept = []

    if model_name in euk_model_list:
        precusor = 'pep_c'  # the seed metabolite
    else:
        precusor = 'pyr_c'  # the seed metabolite

    for met in newmodel.metabolites:
        me = met.elements
        metid = met.id
        if 'C' in me and metid[:-2] != "co2" and metid[:-2] != "hco3" and metid[:-2] != "ACP":
            # print(metid)
            metall.append(metid)
            metpairout = path(newmodel, precusor, metid,
                              result_save_path, specialgroup, specialgroup_begin)
            metpairin = path(newmodel, metid, precusor,
                             result_save_path, specialgroup, specialgroup_begin)

            if metpairout == "no path" or metpairout[2] == 'nosolution':
                if metpairout[2] == 'nosolution':
                    pathexcept.append(metpairout)
                if metpairin == "no path" or metpairin[2] == 'nosolution':
                    bowtie.append([metid, "IS", metid[-1]])
                    if metpairin[2] == 'nosolution':
                        pathexcept.append(metpairin)
                else:
                    pathin.append(metpairin)
                    bowtie.append([metid, "IN", metid[-1]])
            elif metpairin == "no path" or metpairin[2] == 'nosolution':
                if metpairin[2] == 'nosolution':
                    pathexcept.append(metpairin)
                pathout.append(metpairout)
                bowtie.append([metid, "OUT", metid[-1]])
            else:
                pathin.append(metpairin)
                pathout.append(metpairout)
                bowtie.append([metid, "GSC", metid[-1]])

    dfbowtie = pd.DataFrame(
        bowtie, columns=['metabolite', 'bowtie', 'compartment'])
    dfout = pd.DataFrame(pathout, columns=['substrtae', 'formula', 'Cs', 'reduction degree', 'Crd', 'product',
                                           'formula', 'Cp', 'reduction degree', 'Crd', 'rate', 'path yield', 'rdyield', 'difference'])
    dfin = pd.DataFrame(pathin, columns=['substrtae', 'formula', 'Cs', 'reduction degree', 'Crd', 'product',
                                         'formula', 'Cp', 'reduction degree', 'Crd', 'rate', 'path yield', 'rdyield', 'difference'])
    save_path = result_save_path + model_name
    writer = pd.ExcelWriter(save_path+'_results.xlsx')
    dfout.to_excel(writer, 'out', index=False)
    dfin.to_excel(writer, 'in', index=False)
    dfbowtie.to_excel(writer, 'bowtie', index=False)
    if len(pathexcept) != 0:
        dfexcept = pd.DataFrame(pathexcept, columns=[
                                'substrat', 'product', 'state'])
        dfexcept.to_excel(writer, 'except', index=False)
    writer.save()


def add_demand(model, met, specialgroup_begin):
    # add a demand reaction for a metabolite to model, met should be the metabolite object
    # return the demand reaction
    metid = met.id
    # for dpCoA still add normal demand reaction
    if "coa_" in metid and not metid[:-2] in ["coa", 'dpcoa']:
        demand = Reaction('DM_' + metid)
        model.add_reaction(demand)
        # judge whether coa+metid[-2:]exist in the model, if not, use coa_c uniformly
        if 'coa'+metid[-2:] in model.metabolites:
            # metid[-1] for the same compartment
            demand.build_reaction_from_string(
                metid+'--> coa'+metid[-2:], fwd_arrow="-->")
        else:
            demand.build_reaction_from_string(
                metid+'--> coa_c', fwd_arrow="-->")
    # for CoA containing metabolites only consider the main part
    elif "ACP_" in metid and metid[0:3] != "ACP":
        demand = Reaction('DM_' + metid)
        model.add_reaction(demand)
        if 'ACP'+metid[-2:] in model.metabolites:
            demand.build_reaction_from_string(
                metid+'--> ACP'+metid[-2:], fwd_arrow="-->")
        else:
            demand.build_reaction_from_string(
                metid+'--> ACP_c', fwd_arrow="-->")
    # mththf is different
    elif "thf_" in metid and metid[0:3] != "thf" and metid[0:6] != "mththf":
        demand = Reaction('DM_' + metid)
        model.add_reaction(demand)
        if 'thf'+metid[-2:] in model.metabolites:
            demand.build_reaction_from_string(
                metid+'--> thf'+metid[-2:], fwd_arrow="-->")
        else:
            demand.build_reaction_from_string(
                metid+'--> thf_c', fwd_arrow="-->")
    # for cdp containing metabolites but not cdp_c
    elif metid[0:3] in specialgroup_begin and metid[3] != "_" and metid[0:5] != "gdptp":
        demand = Reaction('DM_' + metid)
        model.add_reaction(demand)
        if metid[0:3] == 'uac':  # for uac reduce udp
            if 'udp'+metid[-2:] in model.metabolites:
                demand.build_reaction_from_string(
                    metid+'--> udp'+metid[-2:], fwd_arrow="-->")
            else:
                demand.build_reaction_from_string(
                    metid+'--> udp_c', fwd_arrow="-->")
        else:
            if 'coa'+metid[-2:] in model.metabolites:
                demand.build_reaction_from_string(
                    metid+'--> '+metid[0:3]+metid[-2:], fwd_arrow="-->")
            else:
                demand.build_reaction_from_string(
                    metid+'--> '+metid[0:3]+'_c', fwd_arrow="-->")
    else:
        #dm_met = newmodel.metabolites.get_by_id(met)
        demand = model.add_boundary(met, type='demand')
    return demand
