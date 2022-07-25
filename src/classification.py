#!/usr/bin/env python
# coding: utf-8

# # Gene Classification using ML methods

from sklearn.preprocessing import RobustScaler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.svm import SVC

import os
from datetime import datetime
from tqdm import tqdm

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers


# ## Split Genes into Classes
# P, LP, WN, LN, RN

DISEASE_NAMES = ["C0006142_Malignant_neoplasm_of_breast", "C0009402_Colorectal_Carcinoma", "C0023893_Liver_Cirrhosis_Experimental", 
    "C0036341_Schizophrenia", "C0376358_Malignant_neoplasm_of_prostate"]

SAVE_SCORES = True
SAVE_MODELS = True

train_ratio = 0.70
now = datetime.now()
current_time = now.strftime("%Y_%m_%d-%H_%M_%S")


for DISEASE_NAME in tqdm(DISEASE_NAMES):
    print("Disease: ", DISEASE_NAME)

    scores_RF = []
    scores_SVM = []
    scores_MLP = []

    
    GENE_APU_SCORES_PATH = "../data/APU_scores/" + DISEASE_NAME + "_features_Score.csv"
    GENE_FEATURES_PATH = "../data/NeDBIT_features/" + DISEASE_NAME + "_features"
    SAVE_PATH = "../results/classification/" + DISEASE_NAME + "/"
    TRAIN_SEEDS_PATH = "../data/curated_seed_genes/" + DISEASE_NAME + "_seed_genes.txt"  
    MODEL_PATH = "../models/" + DISEASE_NAME + "_"

    if not os.path.exists(SAVE_PATH):
            os.makedirs(SAVE_PATH)

    if not os.path.exists(MODEL_PATH):
            os.makedirs(MODEL_PATH)

    #load APU scores
    APU_scores_df = pd.read_csv(GENE_APU_SCORES_PATH) #diffusion scores for ALL THE GENES in the PPI

    train_seeds_df = pd.read_csv(TRAIN_SEEDS_PATH, header = None, sep = " ")
    train_seeds_df.columns = ["name", "GDA Score"]
    train_seeds_list = train_seeds_df["name"].values.tolist()

    gene_features_df = pd.read_csv(GENE_FEATURES_PATH)
    gene_features_df['name']=gene_features_df['name'].astype(str)
    
    APU_scores_df_not_seeds = APU_scores_df[~APU_scores_df['name'].isin(train_seeds_list)]
    APU_scores_df_not_seeds.shape

    APU_scores_df_not_seeds = APU_scores_df_not_seeds.sort_values(by = "out", ascending = False)
    


    # ### Pseudo-label assignment
    # We assign genes to pseudo-label according to their APU score
    pseudo_labels = pd.qcut(x = APU_scores_df_not_seeds["out"], q = 4, labels = ["RN", "LN", "WN", "LP"])
    APU_scores_df_not_seeds["label"] = pseudo_labels
    APU_scores_seeds = APU_scores_df[APU_scores_df['name'].isin(train_seeds_list)]
    APU_scores_seeds["label"] = "P"
    APU_scores_seeds

    genes_labels_df = pd.concat([APU_scores_seeds, APU_scores_df_not_seeds])
    genes_labels_df['name']=genes_labels_df['name'].astype(str)
    genes_labels_df

    features_and_labels_df = gene_features_df.merge(genes_labels_df, on = "name")
    features_and_labels_df["label"] = features_and_labels_df["label"].map({"P": 0, "LP": 1, "WN": 2, "LN":3, "RN":4})
    features_and_labels_df["label_binary"] = features_and_labels_df["label"].map({0: 0, 1: 1, 2: 1, 3:1, 4:1})
    features_and_labels_df.pop("class") #removing class attribute
    features_and_labels_df.pop("out")
    features_and_labels_df


    # ### Train-test split

    yData = features_and_labels_df.pop("label").values
    yData_binary = features_and_labels_df.pop("label_binary").values
    geneNames = features_and_labels_df.pop("name").values
 
    xData = features_and_labels_df.values

    num_features = xData.shape[1]

    # ### scale values using robustscaler

    transformer = RobustScaler().fit(xData)
    xData = transformer.transform(xData)

    import numpy as np
    np.random.seed(42)

    
    

    # train is now 75% of the entire data set
    X_train_multi, X_test_multi, y_train_multi, y_test_multi = train_test_split(xData, yData, test_size= 1 - train_ratio, random_state = 42, shuffle = True, stratify = yData)
    X_train_bin, X_test_bin, y_train_bin, y_test_bin = train_test_split(xData, yData_binary, test_size= 1 - train_ratio, random_state = 42, shuffle = True, stratify = yData_binary)

    # ## Model Training

    METHODS = ["RF", "SVM", "MLP"]
    labels_multi =  ["P", "LP", "WN", "LN", "RN"]
    labels_bin = ["P", "U"]

    labellings = ["binary", "multiclass"]

    for METHOD in METHODS:

        for labelling in labellings:
            if labelling == "binary":
                X_train = X_train_bin
                X_test = X_test_bin
                y_train = y_train_bin
                y_test = y_test_bin

                labels = labels_bin
            else:
                X_train = X_train_multi
                X_test = X_test_multi
                y_train = y_train_multi
                y_test = y_test_multi

                labels = labels_multi

            if METHOD == "RF":    
                rf = RandomForestClassifier(random_state=42)

                rf.fit(X_train, y_train)
                y_pred = rf.predict(X_test)

                acc = rf.score(X_test, y_test)
                
                if SAVE_MODELS:
                    filename_model = MODEL_PATH + METHOD + "_" + labelling
                    pickle.dump(rf, open(filename_model, 'wb+'))
                


                #print(class_report_df)
            if METHOD == "SVM":   
                
                clf = SVC(kernel = "rbf", decision_function_shape = "ovo")
                clf.fit(X_train, y_train)
                y_pred = clf.predict(X_test)
                classification_report(y_test, y_pred)
                acc = clf.score(X_test, y_test)

                if SAVE_MODELS:
                    filename_model = MODEL_PATH + METHOD + "_" + labelling
                    pickle.dump(clf, open(filename_model, 'wb+'))
                
            if METHOD == "MLP":
            
                INPUT_SIZE = num_features
                LEARNING_RATE = 1e-4 
                BATCH_SIZE = 32
                EPOCHS = 100
                num_classes = 5 if labelling == "multiclass" else 2

                X_train = np.array(X_train)
                y_train = np.array(y_train)

                def build_model():

                    model = None

                    model = keras.Sequential([
                        layers.Dense(64, activation=tf.nn.relu, input_shape=[INPUT_SIZE]),
                            layers.Dense(32, activation=tf.nn.relu),
                            layers.Dropout(0.3, noise_shape=None, seed=None),
                            layers.Dense(num_classes, activation='softmax')
                        ])

                    optimizer = tf.keras.optimizers.Adam(lr=LEARNING_RATE)

                    # classification
                    model.compile(loss='sparse_categorical_crossentropy',
                                    optimizer=optimizer,
                                    metrics=['accuracy'])

                    return model


                model = build_model()

                es = keras.callbacks.EarlyStopping(
                        monitor='val_loss', mode='min', verbose=1, patience=5, restore_best_weights = True)


                history = model.fit(X_train, y_train, epochs=EPOCHS, validation_data = (np.array(X_test), np.array(y_test)), verbose=1, batch_size=BATCH_SIZE, callbacks = [es])

                y_pred = model.predict_classes(np.array(X_test), verbose=1)
                
                acc = model.evaluate(np.array(X_test), np.array(y_test))[1]

                if SAVE_MODELS:
                    filename_model = MODEL_PATH + METHOD + "_" + labelling
                    model.save(filename_model)
                    model.save(filename_model + ".hdf5")
            
            
            class_report = classification_report(y_test, y_pred, output_dict=True)

            class_report_df = pd.DataFrame(class_report).transpose()
            if SAVE_SCORES:
                class_report_df.to_csv(SAVE_PATH + "/" + DISEASE_NAME + "_" + METHOD + "_" + labelling + ".csv")
                

            #create CM    
            norms = [None, "true"]
            for norm in norms:
                cm = confusion_matrix(y_test, y_pred, normalize=norm)

                plt.figure(figsize=(9,9))
                
                if norm == "true":
                    sns.heatmap(cm, annot=True, fmt=".3f", linewidths=.5, square = True, cmap = 'Blues', xticklabels = labels, yticklabels = labels)
                else:
                    sns.heatmap(cm, annot=True, fmt=".0f", linewidths=.5, square = True, cmap = 'Blues', xticklabels = labels, yticklabels = labels)
                plt.ylabel('Actual label')
                plt.xlabel('Predicted label')
                all_sample_title = 'Accuracy Score: {0}'.format(acc)
                plt.title(all_sample_title, size = 15)

                if SAVE_SCORES:
                    if norm == "true":
                        plt.savefig(SAVE_PATH + "/" + DISEASE_NAME + "_" + METHOD + "_normalized_" + labelling + ".png")
                    else:
                        plt.savefig(SAVE_PATH + "/" + DISEASE_NAME + "_" + METHOD + "_" + labelling  + ".png")

