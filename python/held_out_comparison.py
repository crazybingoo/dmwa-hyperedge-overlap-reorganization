"""Training-fold scaling and fixed logistic models; accepts authorized observations."""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
BASE = ["pair_density", "degree_gini", "clustering", "triplet_density"]
DMWA = ["S", "source_sink_gini_gap", "sink_minus_source_entropy_gap",
        "sink_minus_source_pr_gap", "reciprocal_fraction_mean", "dense_fraction_mean"]

def evaluate(obs):
    predictions=[]
    for patient in sorted(obs.patient_id.unique()):
        train = obs[obs.patient_id!=patient]; test = obs[obs.patient_id==patient]
        scaler = StandardScaler().fit(train[BASE+DMWA])
        def features(frame):
            z = scaler.transform(frame[BASE+DMWA])
            # Average the separately standardized role gaps before rescaling the block.
            return np.column_stack((z[:,:5],z[:,5:8].mean(axis=1),z[:,8:10]))
        tr,te = features(train),features(test)
        block_scale = StandardScaler().fit(tr[:,[5]])
        for m in [tr,te]:m[:,5:6]=block_scale.transform(m[:,[5]])
        pred = test[["patient_id","label"]].copy()
        for name,size in [("base",4),("extended",8)]:
            model = LogisticRegression(C=1.,solver="lbfgs",tol=1e-10,max_iter=1000)
            model.fit(tr[:,:size],train.label)
            pred[name] = model.predict_proba(te[:,:size])[:,1]
        predictions.append(pred)
    preds=pd.concat(predictions,ignore_index=True)
    for model in ["base","extended"]:preds[model+"_brier"]=(preds[model]-preds.label)**2
    scores=preds.groupby("patient_id")[["base_brier","extended_brier"]].mean()
    scores["delta_brier"]=scores.base_brier-scores.extended_brier
    return preds,scores
