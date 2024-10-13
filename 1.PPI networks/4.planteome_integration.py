'''
2023.11.25
Author - Ishani Jayasekara

Input - Modified dataset of planteome data with its uniprot ids
Output - All the planteome data related to photosynthesis in an excel file with sheet name "Planteome"

Get the total list of planteome annotations


'''
import pandas as pd

planteome_file = "maize processing data/Planteome_uniprot_data.xlsx"

# Load Planteome data from the specified sheet
planteome_data = pd.read_excel(planteome_file, sheet_name="Uniprot ids")

# Filter out rows where 'uniprot_ids' column is empty (i.e., '[]')
uniprot_planteome_data = planteome_data[planteome_data["uniprot_ids"]!= "[]"]

uniprot_planteome_data['uniprot_ids'] = uniprot_planteome_data['uniprot_ids'].str.strip('[]').str.replace("'", "") # Clean and split 'uniprot_ids' into separate rows
uniprot_planteome_data = uniprot_planteome_data.assign(uniprot_ids=uniprot_planteome_data['uniprot_ids'].str.split(', ')).explode('uniprot_ids')


def concatenate_values(series):
    """
    Concatenates unique, sorted values from a pandas Series, separating them by commas.
    :param series: pandas Series with values to concatenate
    :return: Concatenated string of unique, sorted values
    """

    unique_values = set()
    for value in series:
        unique_values.update(value.split(', '))
    return ', '.join(sorted(unique_values))

# Group by 'uniprot_ids' and aggregate other columns
uniprot_planteome_data = uniprot_planteome_data.groupby('uniprot_ids').agg({
    'Object name': concatenate_values,
    'GO TERM': concatenate_values,
    'GO NAME': concatenate_values,
    'GO ASPECT': concatenate_values,
    'GO EVIDENCE CODE':concatenate_values,  
    'REFERENCE': concatenate_values  
}).reset_index()

uniprot_planteome_data = uniprot_planteome_data.drop(columns=["Object name"])
uniprot_planteome_data = uniprot_planteome_data.rename(columns={'uniprot_ids': 'GENE PRODUCT ID'})

uniprot_planteome_data['SOURCE'] = 'Planteome'

# Output processed data to the specified Excel file, appending to a new sheet named "Planteome"
ref_file = "maize processing data/maize_ref_data.xlsx"

with pd.ExcelWriter(ref_file, "openpyxl", mode="a") as writer:
    uniprot_planteome_data.to_excel(writer, sheet_name="Planteome", index=False)
    



