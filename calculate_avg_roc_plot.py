import csv
import os
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import numpy as np
from scipy import interp
import sklearn.metrics as metrics
from rdkit.ML.Scoring import Scoring
from sklearn import svm, datasets
from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
#from ccdc.descriptors import StatisticalDescriptors
import operator
def read_csv(filename):
    '''
    ead csv file
    :param filename:
    :return:
    '''
    with open(filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for rows in csv_reader:
            yield rows


def main():


    data_folder = '/home/amukhopadhyay/ligand_screener_data/test_results_data'
    plt.figure(figsize=(5, 5))
    tprs = []
    mean_fpr = np.linspace(0, 1, 100)
    for file in os.listdir(data_folder):
        rows = read_csv(os.path.join(data_folder, file))
        scores = []
        for row in rows:
            scores.append([row[0], int(row[1])])
        fpr, tpr = Scoring.CalcROC(scores, 1)
        tpr = np.array(tpr)
        tprs.append(interp(mean_fpr, fpr, tpr))

        roc_auc = metrics.auc(fpr, tpr)


        plt.plot(fpr, tpr, 'b')
        plt.plot([0, 1], [0, 1],'r--')
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.ylabel('True Positive Rate')
        plt.xlabel('False Positive Rate')

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0

    plt.plot(mean_fpr, mean_tpr, color='green', alpha=.8)



    plt.savefig('avg_roc.png')















if __name__ == "__main__":
    main()