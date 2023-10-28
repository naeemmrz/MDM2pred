# Streamlit functions for frontend and input parsing
import os
import PIL
import glob
import pickle
from math import e
import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


col1, col2, col3 = st.columns([1,4,1])
with col1:
	st.write("")
with col2:
	image = PIL.Image.open('logo.png')
	st.image(image, caption='Cartoon representation of MDM2 from PDB ID: 6I3S', use_column_width=False, width=500)
with col3:
	st.write("")

st.write("""Developed by **Naeem Abdul Ghafoor¹** @[Yildiz Neuro Lab](http://ynlab.mu.edu.tr/en/mdm2pred-6997)""")
st.write("""¹Department of Molecular Biology and Genetics, Graduate School of Natural and Applied Sciences, Mugla Sitki Kocman University, 48000 Mugla, Turkey.\n """)
st.write("\n")
st.markdown("<div style='text-align: justify;'><strong>MDM2pred</strong> is a powerful machine learning tool for predicting the inhibitory potency of compounds against the human E3 ubiquitin ligase MDM2, a key regulator of the tumor suppressor p53. Based on the KNeighbors Regressor algorithm, MDM2pred has been trained on a comprehensive dataset of 1647 known MDM2 inhibitors, achieving an impressive R² value of ~0.74 and an RMSE of ~0.70 (in pIC50 units) over a 10-fold cross-validation. By simply inputting the SMILE notation of any compound, MDM2pred  predicts its pIC50 value against MDM2 and returns the result as IC50. MDM2pred can be a valuable resource for researchers and drug developers looking to accelerate their early screening steps.</div>", unsafe_allow_html=True)
st.write("\n")
st.write("\n")


# Processeing input and generating the results
st.subheader("Please enter the SMILES for your compound:")
user_input = st.text_input("", "CC(=O)NC1=CC=C(C=C1)O")

## Input Control
if user_input is None:
	st.write(f"Waiting user input")
else:
	smile = user_input

## Input Conversion
try:
	m = Chem.MolFromSmiles(smile)
	csmi = Chem.rdmolfiles.MolToSmiles(m)
except:
	st.write(f"Please provide a valid SMILE")
	st.stop()
	
## Input Featurization
model_300dim = glob.glob(f"*model_300dim.pkl")
if len(model_300dim) == 0:
	os.system(f"curl -O https://raw.githubusercontent.com/samoturk/mol2vec/master/examples/models/model_300dim.pkl")
else:
	pass
with open('molecule.smi', 'w') as f:
	f.write(f"{csmi}\tid")
os.system('mol2vec featurize -i molecule.smi -o m2v_output.csv -m model_300dim.pkl -r 1 --uncommon UNK')
_ = pd.read_csv('m2v_output.csv')
features = _.drop(['Unnamed: 0', 'Smiles', 'ID'], axis=1)
	
## Prediction pIC50 and IC50 conversion
MDM2_KNN = pickle.load(open('MDM2_M2V_KNN_UA.sav', 'rb'))
pIC50 = round(MDM2_KNN.predict(features)[0], 3)
IC50_M = 10 ** (- pIC50)
IC50_uM = round(float(IC50_M) * (10 ** 9), 3)

## Compound image and name
def smile2png(smile):
	smi = Chem.MolFromSmiles(f'{smile}')
	d = rdMolDraw2D.MolDraw2DCairo(1500, 1500)
	d.DrawMolecule(smi)
	d.FinishDrawing()
	d.WriteDrawingText("input.png")
	return
smile2png(smile)
def smiles_to_iupac(smile):
	import requests
	CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"
	rep = "iupac_name"
	url = CACTUS.format(smile, rep)
	response = requests.get(url)
	response.raise_for_status()
	return response.text
try:
	smile_name = smiles_to_iupac(smile)
except:
	pass
# Displaying the result
st.write(f"\n")
st.write(f"The predict IC50 for the following compound is **{IC50_uM} nM** (pIC50 = {pIC50}).")


col1, col2, col3 = st.columns([1,4,1])
with col1:
	st.write("")
with col2:
	input_mol = PIL.Image.open('input.png')
	st.image(input_mol, use_column_width=False, width=450)
with col3:
	st.write("")

try:
	smile_name = smiles_to_iupac(smile)
	st.write(f"Compound IUPAC name: **{smile_name}**")
except:
	pass

st.write(f"The models benchmarks are:")
results = pd.read_csv('MDM2_M2V_KNN_UP_CV10_Results.tsv', sep=';')
results
