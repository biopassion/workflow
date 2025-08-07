import sys
sys.path.append('../')

## streamlit 
import streamlit as st
import streamlit.components.v1 as components

import os 
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs


# the current dir is where the main app.py is 
main_dir= "./"
# output dir
output_dir=os.path.join(main_dir,"output")
data_dir = os.path.join(main_dir,"data")


# Form title
st.title("Molecule Similarity")

st.markdown("### Input (Smiles)")

# esomeprazole, lansoprazole
input_1 = st.text_input("Molecule 1", 'COc1ccc2nc([nH]c2c1)[S@](=O)Cc1ncc(C)c(OC)c1C')
input_2 = st.text_input("Molecule 2", 'FC(F)(F)COc1ccnc(c1C)CS(=O)c2[nH]c3ccccc3n2')

if input_1 and input_2 and st.button("Submit"):
    
    smiles_list = [input_1, input_2]
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

    mol_1 = mol_list[0]
    mol_2 = mol_list[1]

    img1 = Draw.MolToImage(mol_1)
    img2 =  Draw.MolToImage(mol_2)



    col1, col2 = st.columns([2,2]) 
    with col1:
        st.image(img1)

    # Use the second column for the submit button
    with col2:
        st.image(img2)
        
    #Draw.MolsToGridImage([mol_1, mol_2],  legends=['mol_1','mol_2'], molsPerRow=2, subImgSize=(300, 250))

    rdk_fpg = rdFingerprintGenerator.GetRDKitFPGenerator()
    mfp_fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2)

    rdk_fps = [rdk_fpg.GetFingerprint(mol) for mol in [mol_1, mol_2]]
    mfp_fps = [mfp_fpg.GetFingerprint(mol) for mol in [mol_1, mol_2]]

    # output 
    rd_ts = str(DataStructs.TanimotoSimilarity(rdk_fps[0], rdk_fps[1]))
    mfp_ts = str(DataStructs.TanimotoSimilarity(mfp_fps[0], mfp_fps[1]))

    rd_dice = str(DataStructs.DiceSimilarity(rdk_fps[0], rdk_fps[1]))
    mfp_dice = str(DataStructs.DiceSimilarity(mfp_fps[0], mfp_fps[1]))

    # markdown
    st.markdown(f'RDKit fingerprint similarity using TanimotoSimilarity: {rd_ts}')
    st.markdown(f'RDKit fingerprint similarity using DiceSimilarity: {rd_dice}')
    
    #st.markdown(f'Morgan fingerprint similarity using TanimotoSimilarity: {mfp_ts}')
    #st.markdown(f'Morgan fingerprint similarity using DiceSimilarity: {mfp_dice}')
    
    