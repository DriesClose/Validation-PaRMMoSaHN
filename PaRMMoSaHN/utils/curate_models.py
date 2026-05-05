## TODO: replace Linux paths to more general ones using Pathlib

import cobra as cb
import pandas as pd
import sys
import os

source_model_filename = sys.argv[1]
curated_model_filename = sys.argv[2]
curations_database = sys.argv[3]

def read_curations(curations_database):
    imbalances = pd.read_excel(curations_database, sheet_name = 'corrected_balances').fillna('')
    duplicates = pd.read_excel(curations_database, sheet_name = 'duplicate_reactions').fillna('')
    return (imbalances, duplicates)

def correct_imbalances(model, imbalances):
    # Loops over curation DB
    for rxn_name,form in zip(imbalances['reaction_id'], imbalances['new_formula']):
        try:
            rxn = model.reactions.get_by_id(rxn_name)
        except KeyError:
            continue
        if len(form) > 0:
            print('Correcting imbalanced ' + rxn_name)
            rxn.build_reaction_from_string(form)
        else:
            print('Removing imbalanced ' + rxn_name)
            rxn.remove_from_model(remove_orphans = True)
    return model

def remove_duplicates(model, duplicates):
    # Loops over curation DB
    for rxn_pair_names0, decision in zip(duplicates['Pair'], duplicates['Decision']):
        rxn_pair_names = eval('(' + rxn_pair_names0 + ')')
        try:
            rxn1 = model.reactions.get_by_id(rxn_pair_names[0])
            rxn2 = model.reactions.get_by_id(rxn_pair_names[1])
            print('Correcting duplicate pair ' + rxn_pair_names0)
            match decision:
                case 1:
                    print('Removing ' + rxn2.id)
                    rxn2.remove_from_model(remove_orphans = True)
                case -1:
                    print('Removing ' + rxn1.id)
                    rxn1.remove_from_model(remove_orphans = True)
                case 0:
                    print('Removing both ' + rxn1.id + ' & ' + rxn2.id)
                    rxn1.remove_from_model(remove_orphans = True)
                    rxn2.remove_from_model(remove_orphans = True)
                case 2:
                    print('Keeping both reactions ' + rxn1.id  + ' & ' + rxn2.id)
                case _:
                    raise ValueError("Unknown treatment case! Check your curation!")
        except KeyError:
            continue
    return model
        
model = cb.io.read_sbml_model(source_model_filename)
(imbalances, duplicates) = read_curations(curations_database)
model = correct_imbalances(model, imbalances)
model = remove_duplicates(model, duplicates)
cb.io.write_sbml_model(model, curated_model_filename)

