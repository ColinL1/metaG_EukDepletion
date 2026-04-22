#!/usr/bin/python
import os
import json
import pandas as pd

## start actuall script 
def extract_info_from_json(json_file):
    with open(json_file, 'r') as file:
        data = json.load(file)

    summary = data['summary']['after_filtering']
    total_reads = summary['total_reads']
    total_bases = summary['total_bases']

    return {
        'file_name': os.path.basename(json_file),
        'total_reads': total_reads,
        'total_bases': total_bases
    }

# def list_full_paths(directory):
#     return [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith('.json')]

def process_dataframe(data):
    # Convert the data to a DataFrame
    df = pd.DataFrame(data)

    # Extract the common part of the file_name
    df['file_name_og'] = df['file_name']
    df['file_name'] = df['file_name'].str.replace(r'\..*$', '', regex=True)

    # Create 'Kingdom' column based on patterns in 'file_name_og' in a single row
    df['Kingdom'] = df['file_name_og'].replace({
        r'^(.*?)\.bacteria\.json$': 'bacteria mapped',
        r'^(.*?)\.non-bacteria\.corals\.mapped\.json$': 'coral mapped',
        r'^(.*?)\.non-bacteria\.corals\.unmapped\.symbiodiniaceae\.mapped\.json$': 'symbiodiniaceae mapped',
        r'^(.*?)\.non-bacteria\.corals\.unmapped\.symbiodiniaceae\.unmapped\.json$': 'unmapped'
    }, regex=True)

    # Replace "Kingdom" with "TOTAL" where "file_name_og" is equal to "file_name" + ".json"
    df.loc[df['file_name_og'] == df['file_name'] + '.json', 'Kingdom'] = 'TOTAL'

    # Filter out rows where 'Kingdom' and 'file_name_og' are the same
    df_filtered = df[df['Kingdom'] != df['file_name_og']]

    # Pivot the filtered DataFrame
    df_pivoted = df_filtered.pivot_table(index='file_name', columns=['Kingdom'], values=['total_reads', 'total_bases'], aggfunc='sum')

    # split and merge DataFrame (temporary fix)
    df_pivoted_split_1 = df_pivoted.loc[:,"total_bases"]
    df_pivoted_split_1['Type'] = 'total_bases'

    df_pivoted_split_2 = df_pivoted.loc[:,"total_reads"]
    df_pivoted_split_2['Type'] = 'total_reads'
    final_df = pd.concat([df_pivoted_split_1, df_pivoted_split_2])

    return final_df

# json_files = list_full_paths(args.json)

json_files = [file for file in os.listdir('.') if file.endswith('.json')]

data_list = []
for json_file in json_files:
    data_list.append(extract_info_from_json(json_file))
    
result_df = process_dataframe(data_list)
result_df.to_csv("report.csv")