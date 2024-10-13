'''
2023.11.26
Author - Ishani Jayasekara

Input - Total lists of proteins from quickgo, amigo and planteome databases in excel 
Output - All proteins retrieved from databases in excel file

Concatonate data from quickgo, amigo, and planteome

'''
import pandas as pd

ref_file = "maize processing data/maize_ref_data.xlsx"

# Load data from QuickGO, Amigo, and Planteome sheets
quickgo = pd.read_excel(ref_file, sheet_name="Quickgo")
amigo = pd.read_excel(ref_file, sheet_name="Amigo")
planteome = pd.read_excel(ref_file, sheet_name="Planteome")

concat_data = pd.concat([quickgo, amigo, planteome]) # Concatenate the data from the three sources into one DataFrame

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

concat_data = concat_data.astype(str) # Convert all data to string type to ensure consistency during aggregation

# Define aggregation functions for each column to handle duplicate entries
aggregation_functions = {
    'GO TERM': concatenate_values,
    'GO NAME': concatenate_values,
    'GO ASPECT': concatenate_values,
    'GO EVIDENCE CODE': concatenate_values,
    'REFERENCE': concatenate_values, 
    'SOURCE': concatenate_values
}

# Group the concatenated data by 'GENE PRODUCT ID' and apply aggregation functions
concat_data = concat_data.groupby("GENE PRODUCT ID").agg(aggregation_functions).reset_index()

# Output the processed data to the specified Excel file, appending to a new sheet named "Total_data"
with pd.ExcelWriter(ref_file, 'openpyxl', mode='a') as writer:
    concat_data.to_excel(writer, sheet_name="Total_data", index=False)