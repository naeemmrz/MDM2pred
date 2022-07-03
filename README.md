# MDM2pred

Developed by **Naeem Abdul Ghafoor¹** & **Ayşegül Yildiz¹** @[Yildiz Neuro Lab](https://ynlab.mu.edu.tr/en)
¹Mugla Sitki Kocman University, Faculty of Science, Dept. of Molecular Biology & Genetics, Mugla, Turkey.

**MDM2pred** is a machine learning application based on the KNNRegressor algorithm, it's trained on 1647 known inhibitors of the human E3 ubiquitin ligase, the primary negative regulator of the well-known suppressor p53. The KNN model backing MDM2pred achieves ~0.74 R² on test compounds (cross-validated) and has an RMSE of ~0.70 (pIC50 unit), the application takes the SMILE of any compound and predicts its pIC50 against MDM2, returning the result as IC50.
