'''
2023.11.25
Author - Ishani Jayasekara

Input- Raw datasets downloaded from GO-amigo api
Output- Modified dataset with uniprot ids

Retrieve uniprot ids for ensemble ids in Planteome dataset
'''

import pandas as pd
import requests, json

taxon_id = 4577

# rice taxon id = 39947
# maize taxon id = 4577

planteome_file = "maize processing data/maize_raw_data.xlsx"

# Load the Excel file
xls = pd.ExcelFile(planteome_file)

concatenated_data = pd.DataFrame() # Initialize an empty DataFrame to store the concatenated data

# Concatenate the Planteome data from each sheet
for sheet_name in xls.sheet_names:
    if sheet_name.lower().startswith("planteome"):
        sheet_data = pd.read_excel(planteome_file, sheet_name=sheet_name, header=None)
        concatenated_data = pd.concat([concatenated_data, sheet_data], ignore_index=True)

concatenated_data = concatenated_data.astype(str)

column_names = ['Object', 'Object name', 'Object type', 'GO TERM', 'GO NAME', 'GO ASPECT', 'Annotation extension', 'Taxon id', 'GO EVIDENCE CODE', 
                'With/from', 'REFERENCE', 'Assigned by']
concatenated_data.columns = column_names

drop_columns = ['Object', 'Object type', 'Annotation extension', 'Taxon id', 'With/from', 'Assigned by']
planteome_data = concatenated_data.drop(columns=drop_columns)


def concatenate_values(series):
    """
    Concatenate unique values from a pandas Series.
    Each value in the Series is expected to be a string with comma-separated items.
    """
    unique_values = set()
    for value in series:
        unique_values.update(value.split(', '))
    return ', '.join(sorted(unique_values))

aggregation_functions = {
        'GO TERM': concatenate_values,
        'GO NAME': concatenate_values,
        'GO ASPECT': concatenate_values,
        'GO EVIDENCE CODE': concatenate_values,
        'REFERENCE': concatenate_values
    }

# Group by 'Object name' and apply aggregation functions
planteome_data = planteome_data.groupby('Object name').agg(aggregation_functions).reset_index()

# Define the base URL for the UniProt API
url = "https://rest.uniprot.org/uniparc/search?query="

uniprot_id_list = []

for id in planteome_data.iloc[:, 0]:
    
    request = requests.get(f"{url}" + id + "&size=1&&format=json") # Send a request to the UniProt API to get UniParc cross-references for the given ID
    record = json.loads(request.text) # Load the response JSON
    results = record["results"] # Extract the results from the JSON response

    print(id)
    
    if (results != []):
        
        cross_references = results[0]["uniParcCrossReferences"] # Extract cross-references from the first result

        uniprot_id = []
        for ref in cross_references:
            # Check if the reference is from UniProtKB/TrEMBL or UniProtKB/Swiss-Prot and if it matches the taxon ID
            if (ref["database"] == "UniProtKB/TrEMBL" or ref['database'] == "UniProtKB/Swiss-Prot") and "organism" in ref:
                if ref["organism"]["taxonId"] == taxon_id:    
                    uniprot_id.append(ref["id"])
        print(uniprot_id)

        # Append the list of UniProt IDs to the main list
        uniprot_id_list.append(uniprot_id)
    
    else:
        # Append None if no results were found
        uniprot_id_list.append(None)

# Add the UniProt IDs as a new column in the DataFrame                
planteome_data['uniprot_ids'] = uniprot_id_list

# Save the processed data to a new Excel file
with pd.ExcelWriter("maize processing data/Planteome_uniprot_data.xlsx", 'openpyxl', mode='w') as writer:
    planteome_data.to_excel(writer, sheet_name= "Uniprot ids")