import pickle

# Specify the file path where the pickle files are located
file_path = "C:\\Users\\david\\OneDrive\\Desktop\\Imperial College London\\Year 3\\Project\\Results_2\\"
#cell_line = ['231','453','468','BT474','Cal51','MCF10A','MCF10A_TGFB','SUM159','T47D','wm']
cell_line = ['231','453','468','BT474','Cal51','MCF10A','MCF10A_TGFB','SUM159','T47D','wm']

for line in cell_line:
    new_data = [0]*21
    # Loop over the 10 shuffled PCC files labelled 1 to 10
    for i in range(1,11):
        shuffle_path = file_path + str(line) + '_shuffled\\all_bin_width_shuffle' + str(line) + '_' + str(i) + '.plk'
        with open(shuffle_path, 'rb') as file:
            # Deserialize and load the object from the file
            shuffle_data = pickle.load(file)[-21:]
            for j in range(len(shuffle_data)):
                new_data[j] = new_data[j] + shuffle_data[j]

    for j in range(len(new_data)):
        # Average the shuffled PCC values
        new_data[j] = new_data[j] / 10

    # Save averaged PCC values in new file
    with open(file_path + 'all_bin_width_' + str(line) + '_shuffled.plk', 'wb') as file:
        pickle.dump(new_data, file)
