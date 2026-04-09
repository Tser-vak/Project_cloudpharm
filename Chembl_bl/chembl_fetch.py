import sqlite3
import pandas as pd
import sys
import os

# Add the directory containing Similar_labels.py to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Similar_labels import exc_extract
from Comp_old_new import comp_data

def read_sq(path,query):
    conn = sqlite3.connect(path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df

#def type_classification(df):
    if df == 4:
        return 'Approved'
    elif df == 3:
        return 'Phase 3'
    elif df == 2:
        return 'Phase 2'
    elif df == 1:
        return 'Phase 1'
    elif df == 0.5:
        return 'Early Phase 1'
    else:
        return 'Unknown'


#def classifie_phase(df):
    df['max_phase'] = df['max_phase'].astype(int)
    df['max_phase'] = df['max_phase'].apply(type_classification)
    return df


if __name__ == "__main__":
    path=r"C:\Users\Max\Desktop\New folder\chembl_36\chembl_36_sqlite\chembl_36.db"
    #keeping the apprioved with the proteins that are null 
    query= """SELECT 
    md.chembl_id,
    md.pref_name,
    md.max_phase,
    cs.canonical_smiles AS smiles,
    cs.standard_inchi_key AS inchi_key
FROM 
    molecule_dictionary md
INNER JOIN 
    compound_structures cs ON md.molregno = cs.molregno
INNER JOIN 
    compound_properties cp ON md.molregno = cp.molregno
WHERE 
    -- Keep phases 1 through 4
    md.max_phase IN (1, 2, 3, 4)
    
    -- 1. Exclude the standalone word "ION" (but keep CATION, etc.)
    AND (md.pref_name IS NULL OR NOT (
        md.pref_name = 'ION' OR
        md.pref_name LIKE 'ION %' OR
        md.pref_name LIKE '% ION %' OR
        md.pref_name LIKE '% ION'
    ))
    
    -- 2. Keep only valid SMILES strings that are 7 or more characters long
    AND cs.canonical_smiles IS NOT NULL
    AND LENGTH(cs.canonical_smiles) >= 7
    
    -- 3. Keep only compounds with a molecular weight strictly greater than 120
    AND cp.full_mwt > 120
ORDER BY 
    md.chembl_id;
"""
#Here we run the query to extract the data from the chembl database and store it in a dataframe called df,
# then we can use the classifie_phase function to classify the phases of the drugs if needed.
df=read_sq(path,query)
data_num=df['max_phase'].value_counts()
exel_map={ 
        'Full Data':(df,False),
        'Data Summary':(data_num,True) 
    }# excel map to gain the full data and the number of classes for each drug status

exc_extract(r'C:\Users\Max\Desktop\Project\Cloudscreen_data_engineering\chembl_bulk\chembl_results','Chembl_data_summary.xlsx',exel_map)

#===Comapiring the data if all phases from server excel with th chembl results
#path_1=r'C:\Users\Max\Desktop\Project\Cloudscreen_data_engineering\chembl_bulk\chembl_results\Chembl_data_summary.xlsx'
#path_2=r'C:\Users\Max\Desktop\Project\Cloudscreen_data_engineering\results\Data_Chmbl_pub.xlsx'
#(merged_df,not_merged_df)=comp_data(path_1,path_2)
#path_out=r'C:\Users\Max\Desktop\Project\Cloudscreen_data_engineering\chembl_bulk\chembl_results\Merged_Chembl_pub.csv'
#merged_df.to_csv(path_out, index=False)
#runs the def classifie_phase to classify the phases of the drugs if needed
#df1=df.copy()
#df1=classifie_phase(df1)