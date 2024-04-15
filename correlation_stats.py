import pickle
from os import listdir
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

def mean_confidence_interval(data, confidence=0.95):
    # Calculates mean and CI for a list of data
    a = 1.0 * np.array(data)
    n = len(a)
    mean, se = np.mean(a), scipy.stats.sem(a)
    ME_CI = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return mean, mean-ME_CI, mean+ME_CI

# Specify the file path where the pickle files are located
folder_path = "C:\\Users\\david\\OneDrive\\Desktop\\Imperial College London\\Year 3\\Project\\Results_2\\"

#cell_line = ['231','453','468','BT474','Cal51','MCF10A','MCF10A_TGFB','SUM159','T47D','wm']
cell_line = ['468']
    
log_space_numbers = np.logspace(0, 2, num=21, base=10)
bin_sizes = log_space_numbers.tolist()
#print(bin_sizes)

for line in cell_line:
    # Iterate over each cell line
    print(str(line))

    # Obtain the files for each bin size
    file_list = listdir(
        folder_path + str(line))
    # Reorganise the bin sizes into order
    for i in range(len(file_list)):
        file_list[i] = float(file_list[i][10:])
    file_list = sorted(file_list)
    for i in range(len(file_list)):
        file_list[i] = 'bin_width_' + str(file_list[i])

    CI_list = []
    corr_list_list = [[] for i in range(len(bin_sizes))]
    its = 0

    for folder in file_list:
        # Iterate over each bin size

        # Obtain the files for each FOV
        plk_list = listdir(
            folder_path + str(line) + '\\' + folder)

        corr_list = []

        for i in range(len(plk_list)):
            # For each FOV
            file_path = folder_path + str(line) + '\\' + folder + '\\' + plk_list[i]
            with open(file_path, 'rb') as file:
                # Deserialize and load the object from the file
                file_data = pickle.load(file)
                # Only use FOVs with more than one cell
                if len(file_data) != 0:
                    corr_list.append(
                        sum(file_data)/len(file_data))
        
        #print(corr_list)
        corr_list_list[its] = corr_list
        CI_list.append(mean_confidence_interval(corr_list, 0.95))
        its = its + 1

    # Calculate mean and CI for observed data
    mean_list = []
    lowerCI_list = []
    upperCI_list = []

    for i in range(len(CI_list)):
        mean_list.append(CI_list[i][0])
        lowerCI_list.append(CI_list[i][1])
        upperCI_list.append(CI_list[i][2])    

    # Calculate mean and CI for shuffled data
    shuffle_file_list = listdir(
        folder_path + str(line) + '_shuffled')
    shuffle_CI_list = []
    shuffle_corr_list = [[] for i in range(len(bin_sizes))]

    for i in range(len(shuffle_file_list)):
        shuffle_path = folder_path + str(line) + '_shuffled\\' + shuffle_file_list[i]
        with open(shuffle_path, 'rb') as file:
            # Deserialize and load the object from the file
            shuffle_data = pickle.load(file)[-21:]
            for j in range(len(shuffle_data)):
                shuffle_corr_list[j].append(shuffle_data[j])

    #print(shuffle_corr_list)
    for i in range(len(shuffle_corr_list)):
        shuffle_CI_list.append(mean_confidence_interval(shuffle_corr_list[i], 0.95))
        
    # Organise mean and CI for shuffled data
    shuffle_mean_list = []
    shuffle_lowerCI_list = []
    shuffle_upperCI_list = []

    for i in range(len(CI_list)):
        shuffle_mean_list.append(shuffle_CI_list[i][0])
        shuffle_lowerCI_list.append(shuffle_CI_list[i][1])
        shuffle_upperCI_list.append(shuffle_CI_list[i][2])

    #print('==============================================================')
    
    # Perform a t-test:
    alpha = 0.05

    ttest_list = []
    for i in range(len(bin_sizes)):
        t_stat, p_value = scipy.stats.ttest_ind(corr_list_list[i], shuffle_corr_list[i], equal_var='False')
        ttest_list.append(p_value)

    for i in range(len(bin_sizes)):
        if ttest_list[i] < alpha:
            ttest_list[i] = 0
        else:
            ttest_list[i] = 1

    # Plot observed data and CI
    plt.plot(bin_sizes, mean_list, color='black', marker='o', label='observed')
    for i in range(len(bin_sizes)):
        plt.plot([bin_sizes[i], bin_sizes[i]], [lowerCI_list[i], upperCI_list[i]], color='black', marker='_')

    # Plot randomised data and CI
    plt.plot(bin_sizes, shuffle_mean_list, color='blue', linestyle=':', marker='o', label='Shuffled')
    for i in range(len(bin_sizes)):
        plt.plot([bin_sizes[i], bin_sizes[i]], [shuffle_lowerCI_list[i], shuffle_upperCI_list[i]], color='blue', marker='_')

    # Add labels and a title
    plt.xlabel('log10 Bin size')
    plt.ylabel('log10 Correlation coefficient')
    plt.title('Temporal correlations in ' + str(line))
    plt.legend()

    plt.xscale('log')  # Set x-axis to log scale
    plt.yscale('log')  # Set y-axis to log scale
    plt.ylim([0.001, 1])
    plt.grid(True)

    # Make a graph and save as a png
    #plt.savefig(folder_path+'fig_'+str(line)+'_temporal_correlations.png', transparent=False, dpi='figure', format='png')
    #plt.savefig(folder_path+'equal_axis_'+str(line)+'_temporal_correlations_.png', transparent=False, dpi='figure', format='png')
    #plt.clf()

    # Show graph without saving
    plt.show()
