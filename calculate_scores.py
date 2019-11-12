from rdkit import Chem
from rdkit.ML.Scoring import Scoring
import csv
import matplotlib.pyplot as plt
import numpy as np
import sklearn.metrics as metrics

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
    rows = read_csv('screening_scores.csv')
    scores = []
    for row in rows:
        scores.append([row[0], int(row[1])])

    #print(scores)
    fractions = [0.01, 0.05, 0.1]
    print(Scoring.CalcAUC(scores, 1))
    print(Scoring.CalcBEDROC(scores, 1, 20))
    print(Scoring.CalcEnrichment(scores, 1, fractions))
    print(Scoring.CalcRIE(scores, 1, 20))
    print((Scoring.CalcAUC(scores, 1)))
    print((Scoring.CalcROC(scores, 1)))

    fpr, tpr = Scoring.CalcROC(scores, 1)
    roc_auc = metrics.auc(fpr, tpr)


    plt.title('Receiver Operating Characteristic')
    plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
    plt.legend(loc = 'lower right')
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.savefig('test_roc.png')

if __name__ == "__main__":
    main()