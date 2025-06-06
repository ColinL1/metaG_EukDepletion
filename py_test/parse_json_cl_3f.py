import sys, os
import json
import pandas as pd
import argparse

# use argparse to grab command line arguments
parser = argparse.ArgumentParser("parse fastp json file into a clean csv")

parser.add_argument('-j', '--json', type = str,
                    help = "path to folder containg the json report file from fastp")
parser.add_argument('-o', '--output', type = str,
                    help = "csv output name")

# if no args given, print help and exit
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided
if args.json is None or \
        args.output is None:
    print("\n** a required input is missing\n"
            "** a json file and output name are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)


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

def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith('.json')]

def process_dataframe(data):
    # Convert the data to a DataFrame
    df = pd.DataFrame(data)

    # Extract the common part of the file_name
    df['file_name_og'] = df['file_name']
    df['file_name'] = df['file_name'].str.replace(r'\..*$', '', regex=True)

    # Create 'Kingdom' column based on patterns in 'file_name_og' in a single row
    df['Kingdom'] = df['file_name_og'].replace({
        
        # Por_02_BT_DESS_L3.json
        # r'^(.*?).json$': 'Bacteria',

        # Por_02_BT_DESS_L3.TRIM.bacteria.json
        r'^(.*?)\.TRIM\.bacteria\.json$': 'Bacteria',

        # Por_02_BT_DESS_L3.TRIM.json
        r'^(.*?)\.TRIM\.json$': 'Trim',

        # Por_02_BT_DESS_L3.TRIM.other.json
        r'^(.*?)\.TRIM\.other\.json$': 'Other',

        # Por_02_BT_DESS_L3.TRIM.scleractina.json
        r'^(.*?)\.TRIM\.scleractina\.json$': 'Scleractina',

        # Por_02_BT_DESS_L3.TRIM.symbiodiniaceae.json
        r'^(.*?)\.TRIM\.symbiodiniaceae\.json$': 'Symbiodiniaceae',
        
        # Por_02_BT_DESS_L3.TRIM.scleractina.json
        r'^(.*?)\.TRIM\.aiptasiidae\.json$': 'Aiptasiidae'
        
        # r'^(.*?)\.bacteria\.json$': 'Bacteria',
        # r'^(.*?)\.non-bacteria\.corals\.mapped\.json$': 'Corals',
        # r'^(.*?)\.non-bacteria\.all_scleractina\.mapped\.json$': 'Corals',

        # r'^(.*?)\.non-bacteria\.corals\.unmapped\.symbiodiniaceae\.mapped\.json$': 'Symbiodiniaceae',
        # r'^(.*?)\.non-bacteria\.all_scleractina\.unmapped\.all_symbiodiniaceae\.mapped\.json$': 'Symbiodiniaceae',

        # r'^(.*?)\.non-bacteria\.corals\.unmapped\.symbiodiniaceae\.unmapped\.json$': 'Other',
        # r'^(.*?)\.non-bacteria\.all_scleractina\.unmapped\.all_symbiodiniaceae\.unmapped\.json$': 'Other',

        # r'^(.*?)\.non-bacteria\.aiptasiidae\.mapped\.json$': 'aiptasiidae',
        # r'^(.*?)\.non-bacteria\.all_aiptasiiadae\.mapped\.json$': 'aiptasiidae',

        # r'^(.*?)\.non-bacteria\.aiptasiidae\.unmapped\.symbiodiniaceae\.mapped\.json$': 'Symbiodiniaceae',
        # r'^(.*?)\.non-bacteria\.all_aiptasiiadae\.unmapped\.all_symbiodiniaceae\.mapped\.json$': 'Symbiodiniaceae',

        # r'^(.*?)\.non-bacteria\.aiptasiidae\.unmapped\.symbiodiniaceae\.unmapped\.json$': 'Other',
        # r'^(.*?)\.non-bacteria\.all_aiptasiiadae\.unmapped\.all_symbiodiniaceae\.unmapped\.json$': 'Other'
    }, regex=True)

    # Replace "Kingdom" with "TOTAL" where "file_name_og" is equal to "file_name" + ".json"
    df.loc[df['file_name_og'] == df['file_name'] + '.json', 'Kingdom'] = 'TOTAL'

    # Filter out rows where 'Kingdom' and 'file_name_og' are the same
    df_filtered = df[df['Kingdom'] != df['file_name_og']]

    # Pivot the filtered DataFrame
    df_pivoted = df_filtered.pivot_table(index='file_name', columns=['Kingdom'], values=['total_reads', 'total_bases'], aggfunc = 'sum')

    # split and merge DataFrame (temporary fix)
    # df_pivoted_split_1 = df_filtered.loc[:,"total_bases"]
    # df_pivoted_split_1['Type'] = 'total_bases'

    # df_pivoted_split_2 = df_filtered.loc[:,"total_reads"]
    # df_pivoted_split_2['Type'] = 'total_reads'
    df_filtered_split_1 = df_filtered[["file_name","Kingdom", "total_bases"]]
    df_filtered_split_1['Type'] = 'total_bases'
    df_filtered_split_1 = df_filtered_split_1.rename(columns={'total_bases': 'total'})
    
    df_filtered_split_2 = df_filtered[["file_name","Kingdom", "total_reads"]]
    df_filtered_split_2['Type'] = 'total_reads'
    df_filtered_split_2 = df_filtered_split_2.rename(columns={'total_reads': 'total'})
    
    final_df = pd.concat([df_filtered_split_1, df_filtered_split_2])

    return final_df

json_files = list_full_paths(args.json)

# json_files = [file for file in os.listdir('.') if file.endswith('.json')]

data_list = []
for json_file in json_files:
    data_list.append(extract_info_from_json(json_file))
    
result_df = process_dataframe(data_list)
result_df.to_csv(args.output)