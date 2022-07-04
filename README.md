# MDM2pred
<div style="text-align: justify"> 
**MDM2pred** is a machine learning application based on the KNNRegressor algorithm, it's trained on 1647 known inhibitors of the human E3 ubiquitin ligase, the primary negative regulator of the well-known suppressor p53. The KNN model backing MDM2pred achieves ~0.74 R² on test compounds (cross-validated) and has an RMSE of ~0.70 (pIC50 unit), the application takes the SMILE of any compound and predicts its pIC50 against MDM2, returning the result as IC50.
</div>

## Authors
**Naeem Abdul Ghafoor¹** & **Ayşegül Yildiz¹** @[Yildiz Neuro Lab](https://ynlab.mu.edu.tr/en)
###### ¹Mugla Sitki Kocman University, Faculty of Science, Department of Molecular Biology and Genetics, Mugla, Turkey.
###### [Click Here For The Corresponding Article](https://)

## Usage
## [![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/naeemmrz/maspa.py/main/MasPA.py)
- Click on the "Open in Streamlit" badge above. 
- Enter the SMILE for the compound of your interest.
- A 2D dipiction of the SMILE and its predicted IC50 will be printed out within few seconds.
  
## Run Locally
Download the project

```bash
  wget https://github.com/naeemmrz/mdm2pred.git
```

Unzip and go to the project directory

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

#### The Application will open in your default browser with the same interface as the online version.
  
## Acknowledgements
The development of MDM2pred was funded by the Research Support and Funding Office (BAP) of Mugla Sitki Kocman University under the project number 22/138/01/3/4.
