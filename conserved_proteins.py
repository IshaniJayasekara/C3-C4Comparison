'''
2024.02.20
'''

import pandas as pd

# cons_df = pd.read_excel('AlignNemo/results-new/clsuters.xlsx', sheet_name='maize clusters')
# seeds = pd.read_excel('final dataset/maize processing data/Seeds_maize.xlsx', sheet_name='Seeds')

cons_df = pd.read_excel('AlignNemo/results-new/clsuters.xlsx', sheet_name='rice clusters')
seeds = pd.read_excel('final dataset/rice processing data/Seeds_rice.xlsx', sheet_name='Seeds')

with pd.ExcelWriter("AlignNemo/results-new/r_cons_proteins.xlsx", engine="openpyxl") as w:

    for col in cons_df.columns:
        info_df = pd.DataFrame()
        for id in cons_df[col]:
            if id != None:
                info = seeds[seeds['GENE PRODUCT ID'] == id]
                info_df = pd.concat([info_df, info], ignore_index=True)
        
        info_df.drop(columns=['STRING ID', 'GO TERM', 'GO NAME', 'GO ASPECT', 'GO EVIDENCE CODE', 'REFERENCE', 'SOURCE'])

        def concatenate_values(series):
            unique_values = set()
            for value in series:
                unique_values.update(value.split(', '))
            return ', '.join(sorted(unique_values))

        filer_proteins = info_df.astype(str)

        aggregation_functions = {
            'PREFERRED NAME': concatenate_values,
            'GENE PRODUCT ID': concatenate_values,
            'FUNCTION': concatenate_values
        }

        filer_proteins = filer_proteins.groupby("PROTEIN NAME").agg(aggregation_functions).reset_index()
        print(filer_proteins)

        filer_proteins.to_excel(w, sheet_name=f'cluster_{col}', index=False)
        # print(info_df)
