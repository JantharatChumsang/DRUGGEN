# from email.policy import default
# from tkinter import CENTER
# import Tkinter as tk
import streamlit as st 
from streamlit_option_menu import option_menu
from streamlit_lottie import st_lottie
import streamlit.components.v1 as stc

import codecs
# import pickle
import joblib
import imblearn
import requests
import random

#----------------------------------------------#

import json
import pandas as pd
import numpy as np

# import seaborn as sn
# import matplotlib.pyplot as plt
# from statistics import mean, stdev

#----------------------------------------------#
# import sklearn
# from sklearn import preprocessing
# from sklearn.model_selection import StratifiedKFold
# from sklearn import linear_model
# from sklearn.model_selection import train_test_split
# from sklearn.datasets import make_classification
# from sklearn.metrics import f1_score,accuracy_score
# from sklearn.metrics import plot_confusion_matrix
# from sklearn.metrics import balanced_accuracy_score
# from imblearn.ensemble import BalancedRandomForestClassifier
# from sklearn.pipeline import Pipeline
# from sklearn.impute import SimpleImputer
# from sklearn.ensemble import RandomForestRegressor
# from sklearn.metrics import r2_score,mean_absolute_error, mean_squared_error
# from sklearn.preprocessing import StandardScaler

#------------------------------------------------#

from rdkit.Chem import Descriptors, Lipinski,Draw
# from rdkit.Chem import Draw
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Avalon import pyAvalonTools
from rdkit import Chem, DataStructs
from rdkit.Chem.Lipinski import RotatableBondSmarts
# from rdkit.Chem.Draw import SimilarityMaps, IPythonConsole
# from rdkit.Chem.Draw import rdMolDraw2D
from chembl_webresource_client.new_client import new_client
# from pikachu.general import draw_smiles
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

# import os
# from os import path
# import zipfile
# import glob
# import random
import keras
from keras.utils import np_utils
from keras.models import load_model

import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras import layers
from tensorflow.keras import initializers
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping

#--------------------------------------------------------------------------------------------------------#
#start#
st.set_page_config(page_title='Drug Discovery - Handover Delays',  layout='wide', page_icon=':pill:')

t1, t2 = st.columns((0.15,1)) 
t1.image('images/index4.png', width = 170)
web_title = '<p style="text-align:; color:#3D0E04; font-size: 22px;">Web applications for Breast Cancer Novel Drug Discovery \n Using the ChEMBL Database and Deep Learning approach ChEMBL</p>'
t2.markdown(web_title, unsafe_allow_html=True)
web2_title = '<p style="text-align:; color:#3D0E04; font-size: 18px;">การพัฒนาเว็บแอพพลิเคชั่นสำหรับการค้นคว้ายาใหม่สำหรับมะเร็งเต้านมด้วยฐานข้อมูล</p>'
t2.markdown(web2_title, unsafe_allow_html=True)

### tab bar ####
selected = option_menu(
    menu_title=None, 
    options=["Home", "About us","Check your SMILES molecule", "Predict new SMILES molecule"], 
    icons=["house","book","check2-all","search"],
    menu_icon="cast",
    default_index=0,
    orientation="horizontal", #แนวนอน
    styles={
        "container": {"padding": "0!important","background-color": "#24A4AC"},
        "nav-link": {"font-size": "16px", "text-align": "center", "margin":"8px", "--hover-color": "#eee"},
        "nav-link-selected": {"background-color": "#24A4AC"},
    }
)

#### sticker image ####
def load_lottiefile(filepath: str):
    with open (filepath,"r") as f:
        return json.load(f)


def load_lottieurl(url: str):
    r = requests.get(url)
    if r.status_code != 200:
        return None
    return r.json()

#### import html ####
import streamlit.components.v1 as stc 
def st_webpage(page_html,width=1370,height=1550):
    page_file = codecs.open(page_html,'r')
    page =page_file.read()
    stc.html(page,width=width, height=height , scrolling = False)


#### selected tab bar Home ####
if selected =="Home":
    
    # ---- LOAD ASSETS ----
    st.write("##")

    # lottie_coding = load_lottieurl("https://assets7.lottiefiles.com/packages/lf20_nw19osms.json")
    # st_lottie(lottie_coding, height=450, key="coding")
    with st.expander("Dataset"):
        df1 = pd.read_csv("databi.csv")
        st.write(df1)
            

    # with st.expander("Power BI"):
    # t1,t2 = st.columns((0.15,1))
    # t1.image('images/qrcode_app.powerbi.com.png', width = 175)
    st_webpage('powerBI.html')
      
#-------------------------------------------------------#
if selected =="About us":
    with st.container():
        st.title("About us 👥")
        Welcome_title = '<p style="text-align:left; font-family: Raleway, sans-serif; color:#06BBCC; font-size: 20px; ">  Web applications for Breast Cancer Novel Drug Discovery Using the ChEMBL Database and Deep Learning approach ChEMBL</p>'
        st.markdown(Welcome_title, unsafe_allow_html=True)
        # t1, t2 = st.columns((0.08,1)) 
        st.image('images/chem.png', width = 90)
        st.write("website: https://www.ebi.ac.uk/chembl/")

            
    with st.container():
        st.write("---")
        st.header("Goal of the projects")
        """
        ```
        The user will bring the SMILES molecule. Come to check the Lipinski principle in drug making, see pIC50, active/inactive values, approve/nonapprove values. 
        of SMILES molecules and a new line of SMILES molecule genrates to further develop drug development in the future.
        ```
        """
        st.info('This is the introduction of a drug molecule from an old SMILES molecule into a new one from the eight target protein targets mTOR, HER2, aromatase, CDK4/6, Trop-2, Estrogen Receptor, PI3K and Akt of breast cancer. To researchers or individuals who wish to discover drugs or produce drugs in the drug discovery process to explore the possibilities of molecules before studying further research into the production of future drugs.')
        
        left_column, right_column = st.columns(2)
        with left_column:
            st.write("##")
            Ideal_title = '<p style="font-family: Poppins, sans-serif; color:#06BBCC; font-size: 20px; "> 🧬 Ideal.</p>'
            st.markdown(Ideal_title, unsafe_allow_html=True)
            st.info(
                """ 
                - It is a process of drug discovery that can Commonly discovered by predicting protein target molecules on web applications
                  And the drug molecules obtained from the prediction can be further studied before producing or developing drugs in the future.
             """)
            st.write("##")
            Reality_title = '<p style="font-family: Poppins, sans-serif; color:#06BBCC; font-size:20px; "> 🧬 Reality.</p>'
            st.markdown(Reality_title, unsafe_allow_html=True)
            st.info(
                """ 
                - In reality, the field of medicine is more complex than we think, starting with observing, experimenting and researching the properties of the natural surroundings. 
                  development of medicinal substances with the synthesis of chemical compounds or compounds that imitate important substances in nature, which the discovery and manufacture of each drug knowledge 
                  required in many disciplines Important substances that have medicinal properties and are available for sale. must be extracted synthesis or compound analysis number of more than ten thousand species 
                  To be selected to study the potency and toxicity of the drug in vitro.
            """)

            st.write("##")
            Consequences_title = '<p style="font-family: Poppins, sans-serif; color:#06BBCC; font-size: 20px; "> 🧬 Consequences.</p>'
            st.markdown(Consequences_title, unsafe_allow_html=True)
            st.info(
                """ 
                - To discover and produce the desired drug If there is no application or technology to help at all It takes an average period of up to 15 years and costs a minimum of approximately $800 million. 
                  Therefore, each discovery of a drug takes a long time and a large budget. The group therefore chose to discover new drugs. together with a machine learning model to help mitigate this problem for future drug development circles.
             """)
            st.write("##")
        with right_column:
            st.write("##")
            st.write("##")
            st.write("##")
            st.write("##")
            st.write("##")
            
            lottie2_coding = load_lottieurl("https://assets2.lottiefiles.com/packages/lf20_pk5mpw6j.json")
            st_lottie(lottie2_coding, height=400,  key="coding")
            st.write("##")
        
        st.header("Citation")

        """
        ```
        Chanin Nantasenamat. (2021). Python for Bioinformatics - Drug Discovery Using Machine Learning and Data Analysis. 
        link. https://www.youtube.com/watch?v=jBlTQjcKuaY&list=LL&index=44&ab_channel=freeCodeCamp.org

        Gashawmg. (Aug 3, 2020).Random Forest Model to Predict pIC50 of SARS-CoV 3CL protease inhibitors
        link. https://github.com/gashawmg/pIC50_prediction/blob/master/SARS-CoV%203CLpro_pIC50_prediction%20uisng_Morgan_Fingerprints.ipynb

        TanushGoel. (Aug 22, 2020).PharmaceuticAI
        link. https://github.com/TanushGoel/PharmaceuticAI
        ```
        """

        st.header(":mailbox: Get In Touch With Us!")
        contact_form = """
        <form action="https://formsubmit.co/jantharat.june@mail.kmutt.ac.th" method="POST">
            <input type="hidden" name="_captcha" value="false">
            <input type="text" name="name" placeholder="Your name" required>
            <input type="email" name="email" placeholder="Your email" required>
            <textarea name="message" placeholder="Your message here"></textarea>
            <button type="submit">Send</button>
        </form>
        """
        st.markdown(contact_form, unsafe_allow_html=True)

        # Use Local CSS File
        def local_css(file_name):
            with open(file_name) as f:
                st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
        local_css("style2.css")

        st.header("Feedback")
        t1,t2 = st.columns((0.05,1))
        t1.image('images/sodsri-logo.png', width = 45)
        t2.write("If you have any feedback, please reach out to us at https://www.facebook.com/thaissf.org/, https://thaissf.org/")
        
#---------------------------------------------------------#

if selected =="Check your SMILES molecule":
    st.title(f"Check your SMILES molecule")
    st.write(""" SMILES = Simplified Molecular Input Line Entry Specification """)
    
       
    canonical_smiles = st.text_input("1.Enter your SMILES molecules string")  

    if st.button("Predict"):
        try:
            if canonical_smiles=="" :
                st.write(f"Don't have SMILES molecules")
            
            else:
            
                model3 = joblib.load('pIC50_predictor1.joblib')
                model4 = joblib.load('active-inactive_predictor3.joblib')
                model5 = joblib.load('BalancedRandomForestClassifier_model6.joblib')

                def draw_compound(canonical_smiles):
                    pic = Chem.MolFromSmiles(canonical_smiles)
                    weight = Descriptors.MolWt(pic)
                    return Draw.MolToImage(pic, size=(400,400))
                picim = draw_compound(canonical_smiles)

                t1, t2 = st.columns(2)
                t1.write('')
                t1.write("""<style>.font-family: Poppins, sans-serif; {font-size:15px !important;}</style>""", unsafe_allow_html=True)
                t1.write('<p class="font-family: Poppins, sans-serif;">This is your smile molecule image</p>', unsafe_allow_html=True)
                # mol = Chem.MolFromSmiles(canonical_smiles)
                # col1.image(mol)
                t1.image(picim)
                
            
                def analyze_compound(canonical_smiles):
                    m = Chem.MolFromSmiles(canonical_smiles)
                    t2.success("The Lipinski's Rule stated the following: Molecular weight < 500 Dalton, Octanol-water partition coefficient (LogP) < 5, Hydrogen bond donors < 5, Hydrogen bond acceptors < 10 ")
                    t2.write('<p class="font-family: Poppins, sans-serif;">Molecule Weight: A molecular mass less than 500 daltons </p>', unsafe_allow_html=True)
                    t2.code(Descriptors.MolWt(m))
                    t2.write('<p class="font-family: Poppins, sans-serif;">LogP: An octanol-water partition coefficient (log P) that does not exceed 5</p>', unsafe_allow_html=True)
                    t2.code(Descriptors.MolLogP(m))
                    t2.write('<p class="font-family: Poppins, sans-serif;">Hydrogen bond donors: No more than 5 hydrogen bond donors (the total number of nitrogen–hydrogen and oxygen–hydrogen bonds)</p>', unsafe_allow_html=True)
                    t2.code(Lipinski.NumHDonors(m))
                    t2.write('<p class="font-family: Poppins, sans-serif;">Hydrogen bond acceptors: No more than 10 hydrogen bond acceptors (all nitrogen or oxygen atoms)</p>', unsafe_allow_html=True)
                    t2.code(Lipinski.NumHAcceptors(m))

                    if Descriptors.MolWt(m) <= np.array(500): 
                        if Descriptors.MolLogP(m) <= np.array(5):
                            if Lipinski.NumHDonors(m) <= np.array(5):
                                if Lipinski.NumHAcceptors(m) <= np.array(10):
                                    str = "your smile is well ✔️"
                                    return str
                                else:
                                    str = "Warning!! your SMILES molecule don't pass Lipinski's Rule ❌"
                                    return str
                            else:
                                str = "Warning!! your SMILES molecule don't pass Lipinski's Rule ❌"
                                return str
                        else:
                            str = "Warning!! your SMILES molecule don't pass Lipinski's Rule ❌"
                            return str
                    else:
                        str = "Warning!! your SMILES molecule don't pass Lipinski's Rule ❌"
                        return str
               
                t2.warning(analyze_compound(canonical_smiles))
            

                def prediction_pIC50(canonical_smiles):
                    test_morgan_fps = []
                    mol = Chem.MolFromSmiles(canonical_smiles) 
                    info = {}
                    temp = AllChem.GetMorganFingerprintAsBitVect(mol,2,2048,bitInfo=info)
                    test_morgan_fps.append(temp)
                    prediction = model3.predict(test_morgan_fps)
                    return prediction


                def get_h_bond_donors(mol):
                    idx = 0
                    donors = 0
                    while idx < len(mol)-1:
                        if mol[idx].lower() == "o" or mol[idx].lower() == "n":
                            if mol[idx+1].lower() == "h":
                                donors+=1
                        idx+=1
                    return donors
                def get_h_bond_acceptors(mol):
                    acceptors = 0
                    for i in mol:
                        if i.lower() == "n" or i.lower() == "o":
                            acceptors+=1
                    return acceptors

                m = Chem.MolFromSmiles(canonical_smiles)
                MW = Descriptors.MolWt(m)
                NA = m.GetNumAtoms()
                LP =  Descriptors.MolLogP(m)
                SA =  Descriptors.TPSA(m)
                mdataf = {'Molecule Weight': MW , 'ALogP': LP , 'HBD' : NA , 'HBA': SA}
                dfm  = pd.DataFrame([mdataf])
                my_array = np.array(dfm)

  
                predict_pIC50 = prediction_pIC50(canonical_smiles)
                prediction3 = ' '.join(map(str, predict_pIC50))
               
                
                prediction4 = model4.predict(my_array)
                prediction4_2 = ' '.join(map(str, prediction4))
                predictionprob4 = model4.predict_proba(my_array)
                
                prediction5 = model5.predict(my_array)
                predictionprob5 = model5.predict_proba(my_array)
                prediction5_2 = ' '.join(map(str, prediction5))
                # st.write(prediction5)
                
                # predictionprob44 = ' '.join(map(str, predictionprob4[:,1]))
                predictionprob55 = ' '.join(map(str, predictionprob5[:,1]))

                with open('style.css') as f:
                    st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)

                col1, col2, col3 = st.columns(3)
                col1.write("""<style>.font-family: Poppins, sans-serif; {font-size:15px !important;}</style>""", unsafe_allow_html=True)
                col1.write('<p class="font-family: Poppins, sans-serif;">Predicted your pIC50 from SMILES molecule 👇</p>', unsafe_allow_html=True)
                col1.code(prediction3) 
                
                col2.write("""<style>.font-family: Poppins, sans-serif; {font-size:15px !important;}</style>""", unsafe_allow_html=True)
                col2.write('<p class="font-family: Poppins, sans-serif;">Predicted your active/inactive Drug 👇</p>', unsafe_allow_html=True)
                col2.code(prediction4_2)
                # col2.write('<p class="font-family: Poppins, sans-serif;">Probability value predicted your active/inactive Drug👇</p>', unsafe_allow_html=True)
                # col2.code(predictionprob44)

                col3.write("""<style>.font-family: Poppins, sans-serif; {font-size:15px !important;}</style>""", unsafe_allow_html=True)
                col3.write('<p class="font-family: Poppins, sans-serif;">Predicted your approve/non-approve Drug👇</p>', unsafe_allow_html=True)
                col3.code(prediction5_2)
                col3.write('<p class="font-family: Poppins, sans-serif;">Probability value predicted your approve/non-approve Drug</p>', unsafe_allow_html=True)
                col3.code(predictionprob55)
        except:
             st.error(f"Your SMILES does not meet the principles of the Lipinski Rules!! ❌")

#------------------------------------------------------------#
# if selected =="Predict new SMILES molecule":
#     Welcome_title = '<p style="font-family: Poppins, sans-serif; color:#06BBCC; font-size: 20px; "> Web applications for Breast Cancer Novel Drug Discovery Using the ChEMBL Database and Deep Learning approach ChEMBL</p>'
#     st.markdown(Welcome_title, unsafe_allow_html=True)
#     st.title(f"Check your SMILES molecule")
#     st.write(""" SMILES = Simplified Molecular Input Line Entry Specification """)
#     predict_nsmiles = st.text_input("1.Enter your SMILES molecules string")
    

#     if st.button("Predict"):
#         # try:
#             if predict_nsmiles=="" :
#                 st.write(f"Don't have SMILES molecules")
            
#             else:
#                 df = pd.read_csv("pharmaceuticAI_all_compounds.smiles")
#                 model = load_model('model_final.h5')
#                 original = predict_nsmiles

#                 # double letters for one element turned into single letters that are not in the dataset
#                 double_to_single = {'Si':'q', 'Se':'w', 'Cn':'t', 'Sc':'y', 'Cl':'u', 'Sn':'z', 'Br':'x'} 
#                 single_to_double = {'q':'Si', 'w':'Se', 't':'Cn', 'y':'Sc', 'u':'Cl', 'z':'Sn', 'x':'Br'}
#                 elements_with_double_letters = list(double_to_single)

#                 element_set = ['C', '(', '=', 'O', ')','/','\\','.','@', 'N', 'c', '1', '$', '2', '3', '4', '#', 'n', 'F', 'u', '-', '[', 'H', ']', 's', 'o', 'S', 't', '5', '6', '+', 'P', 'I', 'x', 'y', 'q', 'B', 'w', '7', '8', 'e', '9', 'b', 'p', '%', '0', 'z','a']
#                 n_vocab = len(element_set)
#                 element_to_int = dict(zip(element_set, range(0, n_vocab)))
#                 int_to_element = {v: k for k, v in element_to_int.items()}
#                 sequence_length = 100 

#                 filey = open('pharmaceuticAI_all_compounds.smiles')
#                 structures = [line[:-1] for line in filey]
                
#                 data = structures 
#                 del structures

#                 def gen_structs(data):

#                     structs = []
#                     for structure in data:
#                         i = 0
#                         while i < len(structure):
#                             try:
#                                 if structure[i] + structure[i+1] in elements_with_double_letters:
#                                     structs.append(double_to_single[structure[i] + structure[i+1]])
#                                     i+=2
#                                 else:
#                                     structs.append(structure[i])
#                                     i+=1
#                             except:
#                                     structs.append(structure[i])
#                                     i+=1
#                         structs.append("$") # end token

#                     return structs

#                 def gen_data(structs):

#                     network_inp = []

#                     # create input sequences
#                     for i in range(0, len(structs) - sequence_length):
#                         sequence_in = structs[i:i + sequence_length]
#                         network_inp.append([element_to_int[char] for char in sequence_in])
                        
#                     n_patterns = len(network_inp)

#                     # reshape the input into a format compatible with CuDNNLSTM layers
#                     network_inp = np.reshape(network_inp, (n_patterns, sequence_length))

#                     return network_inp

#                 network_input = gen_data(gen_structs(data))

#                 def get_h_bond_donors(mol):
#                     idx = 0
#                     donors = 0
#                     while idx < len(mol)-1:
#                         if mol[idx].lower() == "o" or mol[idx].lower() == "n":
#                             if mol[idx+1].lower() == "h":
#                                 donors+=1
#                         idx+=1
#                     return donors

#                 def get_h_bond_acceptors(mol):
#                     acceptors = 0
#                     for i in mol:
#                         if i.lower() == "n" or i.lower() == "o":
#                             acceptors+=1
#                     return acceptors
                
#                 # Lipinski's “Rule of Five” - Constraints in order to maintain drug-like character within the compounds
#                 def rule_of_five(molecule):
#                     m = Chem.MolFromSmiles(molecule)
#                     if get_h_bond_donors(molecule) <= 5 and get_h_bond_acceptors(molecule) <= 10 and Descriptors.MolWt(m) <= 500 and Descriptors.MolLogP(m) <= 5:
#                         return True
#                     else:
#                         return False

#                 def similarity(molecule1, molecule2): # fraction of fingerprints the set of two molecules have in common
#                     m1 = Chem.MolFromSmiles(molecule1)
#                     bi1 = {}
#                     fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, radius=2, nBits=2048, bitInfo=bi1)
#                     fp1_bits1 = fp1.GetOnBits()

#                     m2 = Chem.MolFromSmiles(molecule2)
#                     bi2 = {}
#                     fp2 = AllChem.GetMorganFingerprintAsBitVect(m2, radius=2, nBits=2048, bitInfo=bi2)
#                     fp2_bits2 = fp2.GetOnBits()

#                     common = set(fp1_bits1) & set(fp2_bits2)
#                     combined = set(fp1_bits1) | set(fp2_bits2)

#                     return len(common)/len(combined) # recreation of DataStructs.TanimotoSimilarity

#                 def test_molecule(molecule):
#                     if molecule == None or len(molecule) <= 3:
#                         return False
#                     mol = Chem.MolFromSmiles(molecule)
#                     if mol == None:
#                         return False
#                     else:
#                         try:
#                             Draw.MolToImage(mol) # if molecule is not drawable, the molecule is not valid
#                             return True
#                         except:
#                             return False
                
#                 def complete(inp): # helper method for augment - returns a prediction in place of the removed element

#                     for word, initial in double_to_single.items(): 
#                         inp = inp.replace(word, initial) # replace double-letter elements to single-letter for model input

#                     net_inp = []

#                     for i in range(0, len(inp)):
#                         seq_in = inp[i]
#                         net_inp.append([element_to_int[char] for char in seq_in])

#                     found = False
#                     while not found: # find sequence with the end token "$" at the end
#                         start = np.random.randint(0, len(network_input)-1)
#                         pattern = network_input[start]
#                         if int_to_element[np.round(pattern[-1])] == "$":
#                             found = True

#                     net_in = []
#                     for i in pattern[len(inp):sequence_length]:
#                         net_in.append(i)
#                     for j in net_inp:
#                         net_in.append(j[0])

#                     net_in = np.reshape(net_in, (1, sequence_length))

#                     prediction = model.predict(net_in, verbose=0) # make prediction
#                     index = np.argmax(prediction)
#                     result = int_to_element[index]

#                     return result
                   
#                 def augment(compound, num_changes): # could enhance the pharmacokinetics and bioactivity of the compound

#                     for word, initial in double_to_single.items(): 
#                         compound = compound.replace(word, initial) # replace double-letter elements to single-letter for model input

#                     changes = np.random.randint(1, num_changes+1)
#                     for i in range(0, changes): # randomly removes certain amount of random elements in SMILES string compound and replaces them with prediction
#                         ind = np.random.randint(0, len(compound))
#                         changed = compound[ind]
#                         new_compound = compound[:ind]
                        
#                         result = complete(new_compound)
#                         if result == "$":
#                             return compound[:ind] # if an end token is predicted, return the part of the compound up to the changed index
#                         else:
#                             compound = compound[:ind] + result + compound[ind+1:] # add the prediction in place of the removed element

#                     return compound # return the augmented compound after all the changes have been made 

#                 def augment_repeat(inp, sim, changes, max_try):

#                     if len(inp) > sequence_length:
#                         inp = inp[:sequence_length]

#                     tries = 0
#                     while tries < max_try: # keep trying to make augmented molecules until the model has exceeded the max number of tries (max_try)

#                         augmented = augment(inp, changes)
#                         tries += 1 
#                         for word, initial in single_to_double.items(): 
#                             augmented = augmented.replace(word.lower(), initial) # replace single-letter elements back to original double-letter SMILES notation elements
#                         try:
#                             if test_molecule(augmented) and augmented != inp and rule_of_five(augmented): # make sure that the molecules are valid and drug-like
#                                 s = similarity(inp, augmented) # calculate similarity between the original molecule and the augmented molecule
#                                 if sim < s < 1: # make sure the molecule follows the similarity threshold
#                                     st.write("augmented", inp, "-->", augmented, "with similarity", s)
#                                     return augmented
#                         except:
#                             continue

#                     st.write("could not augment", inp, "within", tries, "tries") 

#                 def aug_list(inp_list, similarity=0.2, max_changes=20, max_tries=100):

#                     molecules = [augment_repeat(compound.replace('/', "").replace('@', "").replace('\\', "").replace('.', ""), similarity, max_changes, max_tries) for compound in inp_list]
#                     molecules = list(set(molecules)) # remove duplicates
#                     molecules = [i for i in molecules if i] 
#                     return molecules

                

#                 #encoder.save('final_model.h5')
#                 # model.load_weights("SMILES-best(2).hdf5")
                
#                 original = str(predict_nsmiles)
#                 augmented = aug_list([original])
#                 st.write(original) 
#                 st.write(augmented) 
                            
        # except:
        #      st.error(f"Your SMILES does not meet the principles of the Lipinski Rules!! ❌")


