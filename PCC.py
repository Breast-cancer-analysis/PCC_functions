import cancer_functions as canf
import correlation_functions as corrf
from FOV_cell import get_all_filenames
from FOV_cell import FOV
from FOV_cell import cell
import pickle
import numpy as np
import os
import random

log_space_numbers = np.logspace(0, 2, num=21, base=10)
bin_sizes = log_space_numbers.tolist()
print(bin_sizes)

shuffle_mode = True

temp_mean_list = []
mean_list_bin = []
path_name = "C:\\Users\\david\\OneDrive\\Desktop\\Imperial College London\\Year 3\\Project\\Results_2\\"
#cell_line = ['231','453','468','BT474','Cal51','MCF10A','MCF10A_TGFB','SUM159','T47D','wm']
cell_line = ['231','453','468','BT474','Cal51','MCF10A','MCF10A_TGFB','SUM159','T47D','wm']


# cell_line = [231] # used to test the if the result is the same as Foust
for i in cell_line:
    if shuffle_mode == False:
        folder_name = str(i)
        iterations = 1
    elif shuffle_mode == True:
        # Perform shuffled PCC 10 times
        # Use Average_shuffled_data.py to average the 10 shuffled files to remove noise
        folder_name = str(i) + '_shuffled'
        iterations = 10
    else:
        exit()

    # Combine the folder name with the desired path
    folder_path = os.path.join(path_name, folder_name)
    # Use os.makedirs() to create the folder and any necessary parent directories
    os.makedirs(folder_path)

    for its in range(iterations):
        print(i)
        print(str(its+1) + '/' + str(iterations))

        csv_list = get_all_filenames(
            "C:\\Users\\david\\OneDrive\\Desktop\\Imperial College London\\Year 3\\Project\\" + str(i))# + "\\" + str(i))
        print(csv_list)
        for d in bin_sizes:
            if shuffle_mode == False:
                # Only create folder to store detailed PCC information when we are not in shuffle mode
                folder_name = 'bin_width_'+str(d)
                # Combine the folder name with the desired path
                folder_path = os.path.join(path_name+str(i), folder_name)
                # Use os.makedirs() to create the folder and any necessary parent directories
                os.makedirs(folder_path)

            for j in csv_list:
                # iterate each file(each FOV)
                view = FOV(j)
                print(view.name)
                all_spike = []
                for z in view.cells:
                    # Iterate over each cell
                    # There are two outputs for the get_spike_train function
                    temp_list = z.get_spike_train(bin_width=d)[0]

                    if shuffle_mode == True:
                        random.shuffle(temp_list)

                    # Some cells have no events at all, need to exclude them
                    if z.no_event == False:
                        all_spike.append(temp_list)
                    elif z.no_event == True:
                        view.cells_no_event.append(z)

                correlation_list = corrf.calculate_pairwise_corrs(
                    np.array(all_spike))
                # print(correlation_list)

                if len(correlation_list) != 0:
                    # Some FOVs just contain one cell, no way to calculate pairwise PCC of that
                    # Include the mean of FOVs with multiple cells
                    temp_mean_list.append(
                        sum(correlation_list)/len(correlation_list))

                if shuffle_mode == False:
                    # Only save details when we are not in the shuffle mode
                    file_name = 'PCC_' + view.name + '.plk'
                    # Saving each PCC for every pair of cells in every FOV, bin_width, and cell line
                    with open(path_name+str(i)+'\\bin_width_'+str(d)+'\\' + file_name, 'wb') as file:
                        pickle.dump(correlation_list, file)
            mean_list_bin.append(sum(temp_mean_list)/len(temp_mean_list))

            temp_mean_list = []

        # The PCC means which will be plotted
        if shuffle_mode == False:
            with open(path_name + 'all_bin_width_' + str(i) + '.plk', 'wb') as file:
                pickle.dump(mean_list_bin, file)
        elif shuffle_mode == True:
            with open(path_name + str(i) + '_shuffled\\all_bin_width_shuffle' + str(i) + '_' + str(int(its) + 1) + '.plk', 'wb') as file:
                pickle.dump(mean_list_bin, file)

        print(mean_list_bin)

        mean_list_bin = []
