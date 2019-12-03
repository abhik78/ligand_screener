from rdkit import Chem
from rdkit.ML.Scoring import Scoring
import csv
import matplotlib.pyplot as plt
import numpy as np
import sklearn.metrics as metrics
from ccdc.descriptors import StatisticalDescriptors
import operator
import collections
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
    rows = read_csv('/home/amukhopadhyay/ligand_screener_testing/screening_scores.csv')
    scores = []
    for row in rows:
        scores.append([row[0], int(row[1])])

    #print(scores) rdkit methods
    #fractions = [0.01, 0.05, 0.1]
    #print(Scoring.CalcAUC(scores, 1))
    #print(Scoring.CalcBEDROC(scores, 1, 20))
    #print(Scoring.CalcEnrichment(scores, 1, fractions))
    #print(Scoring.CalcRIE(scores, 1, 20))
    #print((Scoring.CalcAUC(scores, 1)))
    #print((Scoring.CalcROC(scores, 1)))

    rank_stats = StatisticalDescriptors.RankStatistics(scores, activity_column=operator.itemgetter(1))
    print(round(rank_stats.EF(0.01), 1))
    print(round(rank_stats.EF(0.02), 1))
    print(round(rank_stats.EF(0.05), 1))
    print(round(rank_stats.EF(0.1), 1))
    print(round(rank_stats.AUC(), 1))
    print(round(rank_stats.BEDROC(alpha=20), 1))
    print(round(rank_stats.RIE(alpha=20), 1))


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