import os
import glob

import numpy as np
import pyvista as pv
from utils import modify_pump_velocity, remove_fields, solve_msfr, svd_compression

print('Solving Parametric MSFR: oscillatory ULOFF scenario')
print(' ')

filename = 'MSFR_oscillULOFF_'

# Define Parameter
tau = 5.
omegas = np.linspace(0.5, 5, 25)

# Define Variables to Keep
energy_groups = 6
var_names = ['flux'+str(ii) for ii in range(1, energy_groups+1)]

prec_groups = 8
var_names.extend(['prec'+str(ii) for ii in range(1, prec_groups+1)])

var_names.extend(['p_rgh', 'q', 'T', 'U'])

# Define if the variabile are vectors or scalars
is_vectors = [False]* len(var_names)
is_vectors[-1] = True

for case_i, omega in enumerate(omegas):

    # Preparing directory
    os.system("cp -r BaseCase "+filename+str(case_i+1))
    os.chdir('./'+filename+str(case_i+1))

    # Modify the pump velocity
    modify_pump_velocity('./constant/fvOptions', tau = tau, omega=omega)

    # Solving MSFR
    solve_msfr(case_i)

    # Removing fields
    remove_fields('./', var_names)

    # Import OpenFOAM data
    os.system("touch case.foam")

    # Create pyvista reader
    foam_files = glob.glob(os.path.join('./', '*.foam'))[0]
    reader = pv.POpenFOAMReader(foam_files)

    # Compress and save data using SVD
    print(' ')
    print('SVD compression')
    svd_compression(reader, var_names, is_vectors, 
                    filename=filename+str(case_i+1))

    # Move results in the parent directory
    os.system("mv "+filename+str(case_i+1)+".npz ../"+filename+str(case_i+1)+".npz")

    # Change directory
    os.chdir('../')
    
    # # Clean Repository
    os.system("rm -r "+filename+str(case_i+1))