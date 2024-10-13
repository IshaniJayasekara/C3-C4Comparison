'''
2024.02.04
Author - Ishani Jayasekara

Input - Total list of proteins from all the databases, List of literature based data
Output - All information related to the protein from STRING database 

Get the information from STRING database

'''

import pandas as pd
import gzip
import requests, json

aliases_data = pd.read_table(gzip.open("final dataset/4577.protein.aliases.v12.0.txt.gz", mode='rb'))

uniprot_data = aliases_data[aliases_data['source'] == 'UniProt_ID']

# Create a mapping dictionary from STRING protein IDs to UniProt IDs
mapping_dict = dict(zip(uniprot_data["#string_protein_id"].str.split('.').str[1], uniprot_data["#string_protein_id"]))

# Load literature data from an Excel sheet and convert all data to string type
lit_data = pd.read_excel('final dataset/literature data.xlsx', sheet_name='Sheet2')
lit_data = lit_data.astype('str')


lit_ids_list = []
lit_df = pd.DataFrame()
for i, row in lit_data.iterrows():
    ids = [str(x).strip() for x in row['Maize'].split(',')] if row['Maize'] != 'nan' else []
    # ids = [str(x).strip() for x in row['Rice'].split(',')] if row['Rice'] != 'nan' else []
    lit_ids_list.extend(ids)

lit_df["GENE PRODUCT ID"] = lit_ids_list
lit_df['SOURCE'] = 'Literature'
lit_df['GO EVIDENCE CODE'] = 'lit'


ref_file = "final dataset/maize processing data/maize_ref_data.xlsx"
# ref_file = "final dataset/rice processing data/rice_ref_data.xlsx"
total_dataset = pd.read_excel(ref_file, sheet_name= "Total_data")

total_dataset = pd.concat([total_dataset, lit_df])

for id, row in total_dataset.iterrows():
    if row["GENE PRODUCT ID"] in mapping_dict.keys():
        total_dataset.loc[id, "STRING ID"] = mapping_dict[row["GENE PRODUCT ID"]]
    
    else:
        total_dataset.loc[id, "STRING ID"] = ""

total_dataset.replace([None], '', inplace= True)

not_in_string = total_dataset[total_dataset["STRING ID"].isna() | (total_dataset["STRING ID"] == "")]
string_data = total_dataset[(total_dataset["STRING ID"] != "")]

def concatenate_values(series):
    unique_values = set()
    for value in series:
        if value != "":
            unique_values.update(value.split(', '))
    return ', '.join(sorted(unique_values))

aggregation_functions = {
        'GENE PRODUCT ID':concatenate_values,
        'GO TERM': concatenate_values,
        'GO NAME': concatenate_values,
        'GO ASPECT': concatenate_values,
        'GO EVIDENCE CODE': concatenate_values,
        'REFERENCE': concatenate_values,
        'SOURCE':concatenate_values
    }


string_data = string_data.astype(str)
string_data = string_data.groupby('STRING ID').agg(aggregation_functions).reset_index()
  
for i, row in string_data.iterrows():
    string_data.loc[i, "GENE PRODUCT ID"] = row['STRING ID'].split('.')[1]


url = "https://string-db.org/api"


preferred_name_list = []
protein_name_list = []
functions_list = []

for i, row in string_data.iterrows():

    id = row["STRING ID"]
    print(id)

    r = requests.get(f"{url}/json/get_string_ids?identifiers=" + id)
    record = json.loads(r.text)
    if (record != []):
        results = record[0]

        string_data.loc[i, "PREFERRED NAME"] = results["preferredName"]
        

        if (";" in results["annotation"]):
            string_data.loc[i, "PROTEIN NAME"] = results["annotation"].split(';')[0].split(',')[0]
            string_data.loc[i, "FUNCTION"] = results["annotation"].split(';')[1].split(',')[0]

            # functions_list.append(results["annotation"].split(';')[1].split(',')[0])
        else:
            string_data.loc[i, "PROTEIN NAME"] = results["annotation"]
            string_data.loc[i, "FUNCTION"] = ""

    else:
        string_data.loc[i, "PREFERRED NAME"] = ""
        string_data.loc[i, "PROTEIN NAME"] = ""
        string_data.loc[i, "FUNCTION"] = ""


names = ('Uncharacterized protein.', 'Uncharacterized mitochondrial protein.')  # unambiguous   
string_data = string_data[~string_data["PROTEIN NAME"].isin(names)]

string_data = string_data.drop_duplicates()

  
evi_codes = ["IEA"]
non_curated_data =  string_data[string_data["GO EVIDENCE CODE"].isin(evi_codes)]
evidenced_data = string_data[~string_data["GO EVIDENCE CODE"].isin(evi_codes)]

with pd.ExcelWriter("final dataset/maize processing data/Seeds_maize.xlsx", engine='openpyxl', mode='w') as writer:
# with pd.ExcelWriter("final dataset/rice processing data/Seeds_rice_test_2.xlsx", engine='openpyxl', mode='w') as writer:

    total_dataset.to_excel(writer, sheet_name= "Total data", index=False)
    lit_df.to_excel(writer, sheet_name='Lit based seeds', index=False)
    not_in_string.to_excel(writer, sheet_name='Not in string', index=False)
    non_curated_data.to_excel(writer, sheet_name="DB predicted data", index=False)
    evidenced_data.to_excel(writer, sheet_name='Seeds', index=False)


# Process literature data to extract Maize IDs and create a DataFrame
lit_ids_list = []
lit_df = pd.DataFrame()
for i, row in lit_data.iterrows():
    ids = [str(x).strip() for x in row['Maize'].split(',')] if row['Maize'] != 'nan' else []
    # Alternatively, to extract Rice IDs, uncomment the following line:
    # ids = [str(x).strip() for x in row['Rice'].split(',')] if row['Rice'] != 'nan' else []
    lit_ids_list.extend(ids)

lit_df["GENE PRODUCT ID"] = lit_ids_list
lit_df['SOURCE'] = 'Literature'
lit_df['GO EVIDENCE CODE'] = 'lit'

# Load the total dataset from an Excel file
ref_file = "final dataset/maize processing data/maize_ref_data.xlsx"
# To process Rice data instead, uncomment the following line:
# ref_file = "final dataset/rice processing data/rice_ref_data.xlsx"
total_dataset = pd.read_excel(ref_file, sheet_name="Total_data")

# Combine the total dataset with the literature-based data
total_dataset = pd.concat([total_dataset, lit_df])

# Map STRING IDs using the UniProt ID mapping dictionary
for id, row in total_dataset.iterrows():
    total_dataset.loc[id, "STRING ID"] = mapping_dict.get(row["GENE PRODUCT ID"], "")

# Replace None values with empty strings for consistency
total_dataset.replace([None], '', inplace=True)

# Separate data into those with and without STRING IDs
not_in_string = total_dataset[total_dataset["STRING ID"].isna() | (total_dataset["STRING ID"] == "")]
string_data = total_dataset[total_dataset["STRING ID"] != ""]

# Define a function to concatenate unique, sorted values
def concatenate_values(series):
    unique_values = set()
    for value in series:
        if value != "":
            unique_values.update(value.split(', '))
    return ', '.join(sorted(unique_values))

# Define aggregation functions for each relevant column
aggregation_functions = {
    'GENE PRODUCT ID': concatenate_values,
    'GO TERM': concatenate_values,
    'GO NAME': concatenate_values,
    'GO ASPECT': concatenate_values,
    'GO EVIDENCE CODE': concatenate_values,
    'REFERENCE': concatenate_values,
    'SOURCE': concatenate_values
}

# Group data by STRING ID and apply the aggregation functions
string_data = string_data.astype(str)
string_data = string_data.groupby('STRING ID').agg(aggregation_functions).reset_index()

# Adjust 'GENE PRODUCT ID' based on the STRING ID
for i, row in string_data.iterrows():
    string_data.loc[i, "GENE PRODUCT ID"] = row['STRING ID'].split('.')[1]

# API URL for fetching additional protein information from the STRING database
url = "https://string-db.org/api"

# Initialize lists to store additional protein information
preferred_name_list = []
protein_name_list = []
functions_list = []

# Fetch protein information from the STRING API for each STRING ID
for i, row in string_data.iterrows():
    id = row["STRING ID"]
    print(id)

    # Request data from the STRING API
    r = requests.get(f"{url}/json/get_string_ids?identifiers=" + id)
    record = json.loads(r.text)
    if record:
        results = record[0]

        # Store the preferred name and other annotations
        string_data.loc[i, "PREFERRED NAME"] = results.get("preferredName", "")
        annotation = results.get("annotation", "")
        if ";" in annotation:
            string_data.loc[i, "PROTEIN NAME"] = annotation.split(';')[0].split(',')[0]
            string_data.loc[i, "FUNCTION"] = annotation.split(';')[1].split(',')[0]
        else:
            string_data.loc[i, "PROTEIN NAME"] = annotation
            string_data.loc[i, "FUNCTION"] = ""

    else:
        # Handle cases where no data is returned
        string_data.loc[i, "PREFERRED NAME"] = ""
        string_data.loc[i, "PROTEIN NAME"] = ""
        string_data.loc[i, "FUNCTION"] = ""

# Filter out rows with uncharacterized proteins
names = ('Uncharacterized protein.', 'Uncharacterized mitochondrial protein.')
string_data = string_data[~string_data["PROTEIN NAME"].isin(names)]

# Remove duplicate rows
string_data = string_data.drop_duplicates()

# Separate data into non-curated and curated sets based on evidence codes
evi_codes = ["IEA"]
non_curated_data = string_data[string_data["GO EVIDENCE CODE"].isin(evi_codes)]
evidenced_data = string_data[~string_data["GO EVIDENCE CODE"].isin(evi_codes)]

# Save the processed data to an Excel file with multiple sheets
with pd.ExcelWriter("final dataset/maize processing data/Seeds_maize.xlsx", engine='openpyxl', mode='w') as writer:
# To save Rice data instead, uncomment the following line:
# with pd.ExcelWriter("final dataset/rice processing data/Seeds_rice_test_2.xlsx", engine='openpyxl', mode='w') as writer:
    total_dataset.to_excel(writer, sheet_name="Total data", index=False)
    lit_df.to_excel(writer, sheet_name='Lit based seeds', index=False)
    not_in_string.to_excel(writer, sheet_name='Not in string', index=False)
    non_curated_data.to_excel(writer, sheet_name="DB predicted data", index=False)
    evidenced_data.to_excel(writer, sheet_name='Seeds', index=False)

