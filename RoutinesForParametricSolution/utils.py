from os import listdir
from os.path import isfile, join
import os

import re 
import numpy as np
import time

from sklearn.utils.extmath import randomized_svd
from tqdm import tqdm

def ritual():
    os.system("source ~/.bashrc")
    os.system("./Allclean")
    os.system("source ~/.bashrc")
    os.system("./Allclean")
    os.system("source ~/.bashrc")

    # Actual run of the code
    # os.system("mpirun -n 4 msfrPimpleStructureFoam > log.msfrPimple")
    os.system("./Allrun")

def modify_pump_velocity(file_path, tau, omega = 0):
    # Read all lines from the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    for i, line in enumerate(lines):
        if "scalar v3 = -26.6*exp(-x/2);" in line:
            lines[i] = f"                 scalar v3 = -26.6*exp(-x/({tau:.4f})) * (1 + 0.2 * sin({omega:.4f} * x)); \n"
            break

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

def solve_msfr(case_i, ):
    
    print('------------------------------------------------------------')
    print('                       Case  '+str(case_i+1))
    print('------------------------------------------------------------')

    # Start time for wall-clock and process time
    start_wall = time.time()  # Wall-clock time
    start_process = time.process_time()  # Process time

    ritual()

    # End time for wall-clock and process time
    end_wall = time.time()
    end_process = time.process_time()

    # Calculate elapsed times
    wall_clock_time = end_wall - start_wall
    process_time = end_process - start_process

    # Print results
    print(f"Wall-clock time: {wall_clock_time:.6f} seconds")
    print(f"Process time: {process_time:.6f} seconds")

    # Clean Repository
    os.system("rm -r processors4")
    # os.system("rm log.*")

def remove_fields(snaps_fold, var_names):

    only_folders = [f for f in listdir(snaps_fold) if not isfile(join(snaps_fold, f))]
    sorted_folders = sorted(only_folders)

    pattern = re.compile(r'\d')

    # Filter folders containing numbers
    numeric_folders = [folder for folder in sorted_folders if pattern.search(folder)]

    # Remove from the candidate folders 0 and 0.orig
    numeric_folders = [item for item in numeric_folders if item not in ('0', '0.orig')]

    for folder in numeric_folders:

        if folder != '0.orig':
            only_files = [f for f in listdir(snaps_fold+folder+'/') if isfile(join(snaps_fold+folder+'/', f))]

            for file in only_files:    
                if not sum([file == field for field in var_names]) == 1:
                    print(snaps_fold+folder+'/'+file)
                    os.remove(snaps_fold+folder+'/'+file)

def import_field(reader, var_name: str):
        """
        Importing all time instances (**skipping zero folder**) from OpenFOAM directory.

        Parameters
        ----------
        var_name : str
            Name of the field to import.
            
        Returns
        -------
        field : list
            Imported list of functions (each element is a `numpy.ndarray`), sorted in time.
        time_instants : list
            Sorted list of time.
        """
        
        field = list()
        time_instants = list()
        for idx_t in tqdm(range(1, len(reader.time_values)), 'Importing '+var_name):
            reader.set_active_time_value(reader.time_values[idx_t])  
            
            field.append(reader.read()['internalMesh'].point_data[var_name])
            time_instants.append(reader.time_values[idx_t])
                
        return field, time_instants

def svd_compression(reader, var_names, is_vectors, filename):

    svd_compression = {
        'u': dict(),
        's': dict(),
        'vh': dict()
    }

    for field_i, field in enumerate(var_names):

        _snap, fom_times = import_field(reader, field)

        snaps_matrix = np.asarray([snap.flatten() for snap in _snap]).T

        svd_output = randomized_svd(snaps_matrix,  n_components = 100, n_iter = 'auto')

        svd_compression['u'][field] = svd_output[0]
        svd_compression['s'][field] = svd_output[1]
        svd_compression['vh'][field] = svd_output[2]
    svd_compression['time'] = fom_times
    svd_compression['is_vector'] = is_vectors

    np.savez_compressed(filename+'.npz', svd_compression)