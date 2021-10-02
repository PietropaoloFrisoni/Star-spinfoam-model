# folder where draws, ampls and operators will be stored
data_folder_path = "/home/frisus95/Scrivania"

# verbosity of execution flux, 0 is minimum, 2 is maximum
verbosity_flux = 2

# verbosity of random walk, 0 is minimum, 2 is maximum
verbosity_random_walk = 1

# combine every computed chain, assemble all data in a csv file (combining all computed configurations) and save it in main folder
assemble_and_save_final_data = true

# print final combined data on terminal after execution  
print_final_data_on_terminal = true

# usage example
Configurations = [[0.5, "BF", 10^6, 10^3, 0.40, [false, false], [true, false], [false, false, 1, 2]],
                  [1.0, "BF", 10^6, 10^3, 0.39, [false, false], [true, false], [false, false, 1, 2]],
                  [1.5, "BF", 10^6, 10^3, 0.37, [false, false], [true, false], [false, false, 1, 2]],
                 ]  

