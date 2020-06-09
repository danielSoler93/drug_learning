#!/usr/bin/env python
# coding: utf-8
import os
import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt
from itertools import product
from sklearn.metrics import confusion_matrix
from tensorboard.plugins.hparams import api as hp
from sklearn.model_selection import train_test_split

PATH_DATA = "../datasets/CYP/"


def construct_optimizer(hparams):
    if hparams[HP_OPTIMIZER] == "adam":
        return tf.keras.optimizers.Adam(learning_rate = hparams[HP_LR])
    elif hparams[HP_OPTIMIZER] == "sgd":
        return tf.keras.optimizers.SGD(learning_rate = hparams[HP_LR])
    elif hparams[HP_OPTIMIZER] == "RMSprop":
        return tf.keras.optimizers.RMSprop(learning_rate = hparams[HP_LR])

def train_test_model(hparams):
    internal_layers = [tf.keras.layers.Dropout(hparams[HP_DROPOUT])]+[tf.keras.layers.Dense(hparams[HP_NEURONS], kernel_regularizer=tf.keras.regularizers.l2(hparams[HP_L2]), activation='relu') for _ in range(hparams[HP_HIDDEN_LAYERS])]
    model = tf.keras.models.Sequential([
        tf.keras.layers.Dense(1024, activation='relu', input_shape=(1024,))]+
        internal_layers+[tf.keras.layers.Dense(1, activation="sigmoid")]
    )
    model.compile(optimizer=construct_optimizer(hparams), loss="binary_crossentropy", metrics=[tf.keras.metrics.BinaryAccuracy(name="binary_accuracy")])
    model.fit(train_data, train_labels, epochs=10, verbose=2)
    _, results = model.evaluate(test_data, test_labels, verbose=0)
    _, results_val = model.evaluate(features_only_2c9, labels_validation_2c9, verbose=0)
    return results, results_val


def run(run_dir, hparams):
  with tf.summary.create_file_writer(run_dir).as_default():
    hp.hparams(hparams)  # record the values used in this trial
    accuracy, accuracy_val = train_test_model(hparams)
    tf.summary.scalar("accuracy", accuracy, step=1)
    tf.summary.scalar("accuracy_val", accuracy_val, step=1)


shared_data = pd.read_csv(os.path.join(PATH_DATA, "shared_set_cyp.csv"))
labels_2c9 = (shared_data["p450-cyp2c9 Activity Outcome"] == "Active").values.astype(int)
labels_3a4 = (shared_data["p450-cyp3a4 Activity Outcome"] == "Active").values.astype(int)
validation_2c9_data = pd.read_csv(os.path.join(PATH_DATA, "only_2c9_set_cyp.csv"))
labels_validation_2c9 = (validation_2c9_data["p450-cyp2c9 Activity Outcome"] == "Active").values.astype(int)
validation_3a4_data = pd.read_csv(os.path.join(PATH_DATA, "only_3a4_set_cyp.csv"))
labels_validation_3a4 = (validation_3a4_data["p450-cyp3a4 Activity Outcome"] == "Active").values.astype(int)

features_shared = np.load("shared_set_features.npy")
features_only_2c9 = np.load("only_2c9_set_features.npy")
features_only_3a4 = np.load("only_3a4_set_features.npy")

train_data, test_data, train_labels, test_labels = train_test_split(features_shared, labels_2c9, stratify=labels_2c9)
set_params = "extended"

if set_params == "extended":
    HP_HIDDEN_LAYERS = hp.HParam("hidden_layers", hp.Discrete(list(range(2, 10))))
    HP_NEURONS = hp.HParam("neurons", hp.Discrete([2**i for i in range(5,10)]))
    HP_DROPOUT = hp.HParam("dropout", hp.Discrete([0.2, 0.5]))
    HP_OPTIMIZER = hp.HParam('optimizer', hp.Discrete(['adam', 'sgd','RMSprop']))
    HP_L2 = hp.HParam('l2 regularizer', hp.Discrete([.001,.01]))
    HP_LR = hp.HParam("learning_rate", hp.Discrete([0.001, 0.01, 0.1, 1.0, 10.0]))
    logs_path = 'logs/hparam_tuning'

elif set_params == "reduced":
    HP_HIDDEN_LAYERS = hp.HParam("hidden_layers", hp.Discrete(list(range(3, 10))))
    HP_NEURONS = hp.HParam("neurons", hp.Discrete([2**i for i in range(7,10)]))
    HP_DROPOUT = hp.HParam("dropout", hp.Discrete([0.2]))
    HP_OPTIMIZER = hp.HParam('optimizer', hp.Discrete(['sgd']))
    HP_L2 = hp.HParam('l2 regularizer', hp.Discrete([.001]))
    HP_LR = hp.HParam("learning_rate", hp.Discrete([0.001, 0.01, 0.1, 1.0, 10.0]))
    logs_path = "logs_2/hparam_tuning"

with tf.summary.create_file_writer(logs_path).as_default():
    hp.hparams_config(hparams=[HP_HIDDEN_LAYERS,HP_NEURONS, HP_DROPOUT, HP_OPTIMIZER, HP_L2, HP_LR],
                      metrics=[hp.Metric("accuracy", display_name='Accuracy'), hp.Metric("accuracy_val", display_name="Validation_accuracy")])

session_num = 0
looping = product(HP_NEURONS.domain.values, HP_HIDDEN_LAYERS.domain.values, HP_DROPOUT.domain.values, HP_OPTIMIZER.domain.values, HP_L2.domain.values, HP_LR.domain.values)
for neurons, hidden_lay, dropout, opt, l2, lr in looping:
    hp_params = {HP_NEURONS: neurons, HP_HIDDEN_LAYERS: hidden_lay, HP_DROPOUT: dropout, HP_OPTIMIZER: opt, HP_L2: l2, HP_LR: lr}
    if session_num % 10 == 0:
        tf.keras.backend.clear_session()
    run_name = f"run_{session_num}"
    print(f"---Starting trial: {run_name}")
    print({h.name: hp_params[h] for h in hp_params})
    run(logs_path + '/' + run_name, hp_params)
    session_num += 1
