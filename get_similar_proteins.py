'''
2024.02.20
'''

import pandas as pd

import pandas as pd
# import requests, json

# unq_data = pd.read_excel('final dataset/maize processing data/Seeds_maize.xlsx', sheet_name= 'Seeds')
unq_data = pd.read_excel('final dataset/rice processing data/Seeds_rice.xlsx', sheet_name= 'Seeds')


print(unq_data)
unq_data['PROTEIN NAME'] =  unq_data['PROTEIN NAME'].str.rstrip('.')
unq_data.drop(columns=['STRING ID', 'GO TERM', 'GO NAME', 'GO ASPECT', 'GO EVIDENCE CODE', 'REFERENCE', 'SOURCE'])

def concatenate_values(series):
    unique_values = set()
    for value in series:
        unique_values.update(value.split(', '))
    return ', '.join(sorted(unique_values))

filer_proteins = unq_data.astype(str)

aggregation_functions = {
    'PREFERRED NAME': concatenate_values,
    'GENE PRODUCT ID': concatenate_values,
    'FUNCTION': concatenate_values
}

filer_proteins = filer_proteins.groupby("PROTEIN NAME").agg(aggregation_functions).reset_index()
print(filer_proteins)

# with pd.ExcelWriter("AlignNemo/results-new/m_proteins.xlsx", 'openpyxl') as writer:
with pd.ExcelWriter("AlignNemo/results-new/r_proteins.xlsx", 'openpyxl') as writer:
    filer_proteins.to_excel(writer, sheet_name="Protein info", index=False)

