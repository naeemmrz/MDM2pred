# MDM2pred

![Logo](https://github.com/naeemmrz/mdm2pred/blob/main/logo.png?raw=true)

<div style="text-align: justify"> 
MDM2pred is a machine learning application based on the KNNRegressor algorithm, it's trained on 1647 known inhibitors of the human E3 ubiquitin ligase (Mouse Double Minute 2; MDM2), the primary negative regulator of the well-known tumor suppressor p53. The KNN model backing MDM2pred achieves ~0.74 R² on test compounds (cross-validated) and has an RMSE of ~0.70 (pIC50 unit), the application takes the SMILE of any compound and predicts its pIC50 against MDM2, returning the result as IC50.
</div>


## Authors
**Naeem Abdul Ghafoor¹²** & **Ayşegül Yildiz²³** @[Yildiz Neuro Lab](https://ynlab.mu.edu.tr/en)
###### ¹UnivLyon, Université Claude Bernard Lyon 1, 69100 Villeurbanne, France
###### ²Department of Molecular Biology and Genetics, Graduate School of Natural and Applied Sciences, Mugla Sitki Kocman University, 48000 Mugla, Turkey
###### ³Department of Molecular Biology and Genetics, Faculty of Science, Mugla Sitki Kocman University, Mugla, Turkey.
###### [Click Here For The Corresponding Article](https://)


## Usage
## [![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://ynlab-mdm2pred.streamlit.app/)
- Click on the "Open in Streamlit" badge above. 
- Enter the SMILES for the compound of your interest.
- A 2D depiction of the SMILES and its predicted IC50 will be printed out within a few seconds.
  
  
## Reproduce the results
Download the project

```bash
  wget https://github.com/naeemmrz/mdm2pred.git
```

Unzip and enter the project directory

```bash
  unzip mdm2pred-main.zip
  cd mdm2pred-main
```

Install dependencies

```bash
  pip install -r requirements.txt
```

Run the Reproduce.py

```bash
  python Reproduce.py
```

The results will be printed out in the terminal.


## Run Locally
Download the project

```bash
  wget https://github.com/naeemmrz/mdm2pred.git
```

Unzip and enter the project directory

```bash
  unzip mdm2pred-main.zip
  cd mdm2pred-main
```

Install dependencies

```bash
  pip install -r requirements.txt
```

Start the application

```bash
  streamlit run MDM2pred.py
```

The Application will open in your default browser with the same interface as the online version.


## Acknowledgements
The development of MDM2pred was funded by the Research Support and Funding Office (BAP) of Mugla Sitki Kocman University under the project number 22/138/01/3/4.
