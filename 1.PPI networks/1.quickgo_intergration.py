'''
2023.11.25
Author - Ishani Jayasekara

Input - Raw datasets downloaded from GO-quickgo api
Output - All the quickgo data related to photosynthesis in an excel file with sheet name "Quickgo"

Get the total list of quickgo annotations

'''

import pandas as pd

def concatenate_values(series):
    """
    Concatenate unique values from a pandas Series.
    Each value in the Series is expected to be a string with comma-separated items.
    """
    unique_values = set()
    for value in series:
        unique_values.update(value.split(', '))
    return ', '.join(sorted(unique_values))

def process_data(sheet_data):
    """
    Process and clean the data from a given DataFrame (sheet_data).
    - Drops unnecessary columns
    - Aggregates values based on 'GENE PRODUCT ID'
    - Applies custom concatenation function to specific columns
    """
    
    df = sheet_data.astype(str)

    drop_columns = ['ANNOTATION EXTENSION','SYMBOL','QUALIFIER', 'WITH/FROM', 'TAXON ID', 'ASSIGNED BY', 'GENE PRODUCT DB', 'ECO ID']
    df = df.drop(columns=drop_columns)
    
    # Define aggregation functions for specific columns
    aggregation_functions = {
        'GO TERM': concatenate_values,
        'GO NAME': concatenate_values,
        'GO ASPECT': concatenate_values,
        'GO EVIDENCE CODE': concatenate_values,
        'REFERENCE': concatenate_values
    }

    # Group by 'GENE PRODUCT ID' and apply aggregation functions
    df = df.groupby("GENE PRODUCT ID").agg(aggregation_functions).reset_index()

    # Add a new column 'SOURCE' with a constant value 'Quickgo' 
    df['SOURCE'] = "Quickgo"
    
    return df


file_path = "maize processing data/maize_raw_data.xlsx"

# Load the Excel file
xls = pd.ExcelFile(file_path)

concatenated_data = pd.DataFrame() 

# Concatenate the data from each sheet into the main DataFrame
for sheet_name in xls.sheet_names:
    if sheet_name.lower().startswith("quickgo"):
        sheet_data = pd.read_excel(file_path, sheet_name=sheet_name)
        concatenated_data = pd.concat([concatenated_data, sheet_data], ignore_index=True)

quickgo_dataset = process_data(concatenated_data) 

# Save the processed data to a new Excel file
with pd.ExcelWriter('maize processing data/maize_ref_data.xlsx', 'openpyxl', 'w')as writer:
    quickgo_dataset.to_excel(writer, sheet_name='Quickgo', index=False)
