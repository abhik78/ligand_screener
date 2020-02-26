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
import argparse
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
        next(csv_reader)
        for rows in csv_reader:
            yield rows

def plot_curve(fpr, tpr, color):
    plt.plot(fpr, tpr, color, label = 'test')
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--result_dir', '-d', help='sdf file to generate conformers')
    parser.add_argument('--target_ids', '-f', help='single column csv file with target uniprot ids, header uniprot_id')
    parser.add_argument('--output_dir', '-o', help='directory where the output file will be written')
    #parser.add_argument('--label', '-l', help= 'method name as label of the plot')
    args = parser.parse_args()
    return args



def main():
    args = parse_arguments()
    unp_id_list = [row[0] for row in read_csv(args.target_ids)]

    category_list = [row[1] for row in read_csv(args.target_ids)]
    plt.figure(figsize=(5, 5))
    for unp_id, category in list(zip(unp_id_list, category_list)):
        print(unp_id)
        sub_dir_list = next(os.walk(os.path.join(args.result_dir)))[1]

        tprs = []
        mean_fpr = np.linspace(0, 1, 100)

        for sub_dir in sub_dir_list:

            scores_dir = os.path.join(args.result_dir, sub_dir, '{}'.format(unp_id))
            print(scores_dir)

            if os.path.isdir(scores_dir):
                os.chdir(scores_dir)

                if os.path.isdir(scores_dir):
                    for filename in os.listdir(scores_dir):

                        if filename.endswith('.csv'):
                            print(filename)

                            rows = read_csv(filename)
                            scores = []
                            for row in rows:
                                scores.append([row[0], int(row[1])])
                            fpr, tpr = Scoring.CalcROC(scores, 1)
                            tpr = np.array(tpr)
                            tprs.append(interp(mean_fpr, fpr, tpr))
        if tprs:


            mean_tpr = np.mean(tprs, axis=0)

            mean_tpr[-1] = 1.0

            if (category) == 'easy':
                plot_curve(fpr=fpr, tpr=tpr, color='blue')


            elif (category) == 'moderate':
                plot_curve(fpr=fpr, tpr=tpr, color='orange')

            elif (category) == 'hard':
                plot_curve(fpr=fpr, tpr=tpr, color='green')

            elif (category) == 'unfeasible':
                plot_curve(fpr=fpr, tpr=tpr, color='magenta')
        else:
            print("no info for target {}".format(unp_id))


    plt.savefig(os.path.join(args.output_dir, 'avg_roc.png'))

if __name__ == "__main__":
    main()