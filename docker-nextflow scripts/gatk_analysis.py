import pandas as pd

#Get bcftools processed vcf file results
tsv_file_path = '/path/to/bcftools/processed/VCF/file/results'

# Read the TSV file into a DataFrame
df1150 = pd.read_csv(tsv_file_path, delimiter='\t')
df1150['Sample'] = 'df1150'

#run previous code with multiple samples, merge them into one dataframe for following analysis
frames = [df1150, df1174,df1342,...]
result = pd.concat(frames)

# filtering coding-regioon based on column 'Consequence'
patterns = ['missense', 'splice_[ad]', 'stop_', 'frameshift', 'inframe_[id]']
regex = '|'.join(patterns)
result = result[result['Consequence'].str.contains(regex, regex=True, na=False)]

#replace the na in Gnomad_AF as 0
result['Gnomad_AF'] = result['Gnomad_AF'].fillna(0)
result['Gnomad_AF'] = pd.to_numeric(result['Gnomad_AF'], errors='coerce')
result = result[result['Gnomad_AF'] < 0.0025]

#calculate AF
def calculate_AF(row):
    # Split the AD string into a list of strings
    AD_list = row['AD'].split(',')

    # Convert the strings to integers
    a = int(AD_list[0])
    b = int(AD_list[1])
    DP = int(row['DP'])

    # Calculate AF
    AF = b / DP

    return AF

result['AF'] = result.apply(calculate_AF, axis=1)

    
#filtering results 
result['DP'] = pd.to_numeric(result['DP'], errors='coerce')
result = result[result['DP'] > 10]
result = result[result['AF'] > 0.01]

result

#import the driver_gene list
file_path = '/path/to/driver/gene/list/'
driver_genes = pd.read_csv(file_path)

#label gatk variants with drivergene
result['drivergene'] = result['SYMBOL'].notna() & result['SYMBOL'].isin(driver_genes['Gene'])
result['drivergene'] = result['drivergene'].map({True: 'yes', False: ''})

#filtering with drvier gene
driver = result[result['drivergene'] == 'yes']
