#! /bin/python

import sys
import pandas as pd

pan_reaction_table_file = sys.argv[1]
strain_blast_hits_file = sys.argv[2]
strain_reaction_table_file = sys.argv[3]

# pan_reaction_table_file = "../02-pan_model/clostridiales-all-Reactions.tbl"
# strain_reaction_table_file = "../03-strain_models/models/test.tbl"
# strain_blast_hits_file = "../03-strain_models/matches/C_acetobutylicum.tsv"

def format_number(num):
    if num.is_integer():
        return "{:.0f}".format(num)
    else:
        return str(num)

# Read homology match table of new strain
strain_blast_hits = pd.read_table(strain_blast_hits_file, header = None, 
                                  usecols=[0,1,11], names = ['Query', 'Target', 'Repr. bitscore'])

# Take the best hit
strain_blast_hits.drop_duplicates(subset = ['Query'], keep='first', inplace=True)

# Read panmodel reaction table, omitting the comment header
pan_reaction_table = pd.read_table(pan_reaction_table_file, comment = '#')

# Join match table and panmodel reaction table by representative gene IDs
strain_reaction_table = pd.merge(strain_blast_hits, pan_reaction_table, 
                                 how = 'inner', left_on='Target', right_on='stitle')
# Drop duplicate entries, retaining the highest ranked record for each query-representative-rxn-pathway combination
# (assuming first has highest bitscore)
strain_reaction_table.drop_duplicates(inplace = True, ignore_index = True,
                                      subset = ['Query', 'Target', 'rxn', 'pathway'])
# Sort by bitscore overall
strain_reaction_table.sort_values(by = "Repr. bitscore", ascending = False, inplace = True)
# Drop duplicate entries, retaining the highest ranked record (just ranked by bitscore) for each query-pathway pair
strain_reaction_table.drop_duplicates(subset = ['Query', 'pathway', 'dbhit'], keep = "first", inplace = True)
# Sort by index again, thus by order of the genes
strain_reaction_table.sort_index(inplace = True)
# Overwrite representative ID by query ID
strain_reaction_table['stitle'] = strain_reaction_table['Query']
# Drop columns from the homology match table
strain_reaction_table.drop(columns = ['Query', 'Target', 'Repr. bitscore'], inplace = True)
# Write filtered reaction table
strain_reaction_table.to_csv(strain_reaction_table_file,
                              sep = '\t', index = False, na_rep = "NA", float_format = format_number)
