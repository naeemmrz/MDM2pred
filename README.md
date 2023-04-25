# MDM2pred

![Logo](https://github.com/naeemmrz/mdm2pred/blob/main/logo.png?raw=true)

<div style="text-align: justify"> 
MDM2pred is a powerful machine learning tool for predicting the inhibitory potency of compounds against the human E3 ubiquitin ligase MDM2, a key regulator of the tumor suppressor p53. Based on the KNeighbors Regressor algorithm, MDM2pred has been trained on a comprehensive dataset of 1647 known MDM2 inhibitors, achieving an impressive R² value of ~0.74 and an RMSE of ~0.70 (in pIC50 units) over a 10-fold cross-validation. By simply inputting the SMILE notation of any compound, MDM2pred  predicts its pIC50 value against MDM2 and returns the result as IC50. MDM2pred can be a valuable resource for researchers and drug developers looking to accelerate their early screening steps.
</div>


## Authors
**Naeem Abdul Ghafoor¹²** & **Ayşegül Yildiz¹²** @[Yildiz Neuro Lab](https://ynlab.mu.edu.tr/en)
###### ¹Department of Molecular Biology and Genetics, Graduate School of Natural and Applied Sciences, Mugla Sitki Kocman University, 48000 Mugla, Turkey
###### ²UnivLyon, Université Claude Bernard Lyon 1, 69100 Villeurbanne, France
###### ³Department of Molecular Biology and Genetics, Faculty of Science, Mugla Sitki Kocman University, Mugla, Turkey.
###### [Click Here For The Corresponding Article](https://)


## Usage
## [![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://ynlab-mdm2pred.streamlit.app/)
- To use MDM2pred, simply click on the "Open in Streamlit" badge above.
- Enter the SMILES notation for your compound of interest.
- Within seconds, MDM2pred will generate a 2D depiction of the compound and its predicted IC50 value against MDM2.
- Use this information to gain valuable insights into the inhibitory potential of your compound and accelerate your drug development efforts.
  
  
## To reproduce the results
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
