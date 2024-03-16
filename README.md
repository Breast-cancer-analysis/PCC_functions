# PCC_functions
David's various PCC codes

### PCC.py
Contains the code to calculate the pairwise PCC values for every pair of cells in each POV, and saves them as pickle files.

**shuffle_mode == False:** Actual obtained data
  
**shuffle_mode == True:** Shuffled PCC values (iterates 10 times)

### Average_shuffled_data
  Averages the PCC values from the 10 shuffled PCC files to smooth out the data and saves them as pickle files

### correlation_stats
  Calculate mean and confidence interval for the PCC values, plot the data, and perform t-tests to determine statistical difference from randomness
