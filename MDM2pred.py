# Streamlit functions for frontend and input parsing
import streamlit as st
import PIL

col1, col2, col3 = st.columns([1,4,1])
with col1:
	st.write("")
with col2:
	image = PIL.Image.open('logo.png')
	st.image(image, use_column_width=False, width=500)
with col3:
	st.write("")

st.write("""Developed by **Naeem Abdul Ghafoor¹** & **Ayşegül Yildiz¹** @[Yildiz Neuro Lab](https://ynlab.mu.edu.tr/en)""")
st.write("""¹Mugla Sitki Kocman University, Faculty of Science, Dept. of Molecular Biology & Genetics, Mugla, Turkey.\n """)
st.write("\n")
st.markdown("<div style='text-align: justify;'><strong>MDM2pred</strong> is a machine learning application based on the KNNRegressor algorithm, it's trained on 1647 known inhibitors of the human E3 ubiquitin ligase, the primary negative regulator of the well-known suppressor p53. The KNN model backing MDM2pred achieves ~0.74 R² on test compounds (cross-validated) and has an RMSE of ~0.70 (pIC50 unit), the application takes the SMILE of any compound and predicts its pIC50 against MDM2, returning the result as IC50.</div>", unsafe_allow_html=True)
st.write("\n")
st.write("\n")
st.subheader("Please enter the SMILE for your compound:")
user_input = st.text_input("", "CC(=O)NC1=CC=C(C=C1)O")

# Functions for running the backend
def smile2png(smile):
	from rdkit import Chem
	from rdkit.Chem.Draw import rdMolDraw2D
	smi = Chem.MolFromSmiles(f'{smile}')
	d = rdMolDraw2D.MolDraw2DCairo(1500, 1500)
	d.DrawMolecule(smi)
	d.FinishDrawing()
	d.WriteDrawingText("input.png")
	return

def smi2canon(smile):
	from rdkit import Chem
	m = Chem.MolFromSmiles(smile)
	csmi = Chem.rdmolfiles.MolToSmiles(m)
	return csmi

def get_m2v(csmi):
	import os
	import pandas as pd
	with open('molecule.smi', 'w') as f:
		f.write(f"{csmi}\tid")
	os.system('mol2vec featurize -i molecule.smi -o m2v_output.csv -m model_300dim.pkl -r 1 --uncommon UNK')
	_ = pd.read_csv('m2v_output.csv')
	features = _.drop(['Unnamed: 0', 'Smiles', 'ID'], axis=1)
	return features

def get_prediction(features):
	import pickle
	MDM2_KNN = pickle.load(open('MDM2_M2V_KNN.sav', 'rb'))
	prediction = MDM2_KNN.predict(features)
	return prediction

def pIC50_2_IC50(pIC50):
	from math import e
	IC50_M = e ** (-pIC50)
	IC50_uM = round(float(IC50_M) * (10 ** 6), 3)
	return IC50_uM

def LOOResultsReproduce():
	"""
	Function to reproduce the reported results for the MDM2pred Model 
	"""
	import pickle
	import pandas as pd
	from sklearn.model_selection import KFold, cross_val_score

	def average(lst):
		return sum(lst) / len(lst)
	KNN_Model = pickle.load(open('MDM2_M2V_KNN.sav', 'rb'))
	X = pd.read_csv('MDM2_M2VX.csv')
	y = pd.read_csv('MDM2.csv')['pIC50']

	cv = KFold(n_splits=10, shuffle=True, random_state=48)
	R2 = cross_val_score(KNN_Model, X, y, scoring='r2', cv=cv, n_jobs=None)      
	RMSE = cross_val_score(KNN_Model, X, y, scoring='neg_root_mean_squared_error', cv=cv, n_jobs=None) * -1
	MSE = cross_val_score(KNN_Model, X, y, scoring='neg_mean_squared_error', cv=cv, n_jobs=None) * -1
	MAE = cross_val_score(KNN_Model, X, y, scoring='neg_mean_absolute_error', cv=cv, n_jobs=None) * -1
	
	results = {'Number of Compounds': X.shape[0], 'LOOCV R2': average(R2), 'LOOCV RMSE': average(RMSE),
	 'LOOCV MSE': average(MSE), 'LOOCV MAE': average(MAE)}
	df = pd.DataFrame(results, index=[0])
	return df


# Processeing input and generating the results
if user_input is None:
	st.write(f"Waiting user input")
else:
	smile = user_input
	try:
		csmi = smi2canon(smile)
		features = get_m2v(csmi)
		pIC50 = get_prediction(features)
		IC50_uM = pIC50_2_IC50(pIC50)
		# Displaying the result
		st.write(f"\n")
		st.write(f"The predict IC50 for the following compound is **{IC50_uM} μM** (pIC50 = {round(pIC50[0], 3)}).")
		smile2png(smile)

		col1, col2, col3 = st.columns([1,4,1])
		with col1:
			st.write("")
		with col2:
			input_mol = PIL.Image.open('input.png')
			st.image(input_mol, use_column_width=False, width=450)
		with col3:
			st.write("")
		
		st.write(f"The models benchmarks are:")
		st.write(LOOResultsReproduce())
	except:
		st.write(f"Please enter a valid SMILE :)")	
