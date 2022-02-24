import os
import time


'''Inputs'''
path_to_gsas2 = "/home/danielmurphy/gsas2/"
save_directory = "/home/danielmurphy/Downloads/GSASdata/new_outputs/"
data_directory = "/home/danielmurphy/Downloads/GSASdata/"
refinement_method = "Pawley"
input_data_files = ["PBSO4.XRA", "PBSO4.CWN"]
instrument_files = ["INST_XRY.PRM","inst_d1a.prm"]
phase_file = "PbSO4-Wyckoff.cif"
project_name = "mantid_test"

x_min = [16.0,19.0]
x_max = [158.4,153.0]


'''Validation'''
if len(input_data_files) != len(instrument_files):
    raise ValueError(f"The number of input_data_files (histograms) ({len(input_data_files)}) must equal the"
                     + f"number of instrument_files ({len(instrument_files)})")

if len(x_min) != len(x_max):
    raise ValueError(f"The number of x_min values ({len(x_min)}) must equal the"
                     + f"number of x_max values ({len(x_max)})")

if len(x_min) != len(input_data_files):
    raise ValueError(f"The number of x_min values ({len(x_min)}) must equal the"
                     + f"number of input_data_files (histograms) ({len(input_data_files)})")


'''exec'''
start = time.time()
os.system(path_to_gsas2 + "bin/python "
          + "/home/danielmurphy/mantid/qt/python/mantidqtinterfaces/mantidqtinterfaces/Engineering/gui/"
          + "engineering_diffraction/tabs/gsas2/call_G2sc.py "
          + path_to_gsas2 + " "
          + save_directory + " "
          + data_directory + " "
          + refinement_method + " "
          + input_data_files[0] + " "
          + input_data_files[1] + " "
          + instrument_files[0] + " "
          + instrument_files[1] + " "
          + phase_file + " "
          + project_name + " "
          + str(x_min[0]) + " "
          + str(x_min[1]) + " "
          + str(x_max[0]) + " "
          + str(x_max[1])
          )
gsas_runtime = time.time() - start
print(f"\nGSASII Complete in {gsas_runtime} seconds.\n")


# project_path = save_directory + project_name + '.gpx'
result_filepath = save_directory + project_name + '.lst'

last_modified_time = os.path.getmtime(result_filepath)

if time.time() > (last_modified_time + 2):
    # if GSASII result file not modified in the last 2 seconds
    m_ti = time.ctime(last_modified_time)
    # Using the timestamp string to create a time object/structure
    t_obj = time.strptime(m_ti)
    # Transforming the time object to a timestamp of ISO 8601 format
    T_stamp = time.strftime("%Y-%m-%d %H:%M:%S", t_obj)
    print(f"The file located at the path {result_filepath} was last modified at {T_stamp}")
    raise ValueError('This GSASII operation failed as output files were from an old operation.')

print(f"GSASII result file found. Opening {result_filepath}")
with open(result_filepath, 'r') as file:
    result_string = file.read().replace('\n', '')
    for loop_histogram in input_data_files:
        where_loop_histogram = result_string.rfind(loop_histogram)
        if where_loop_histogram != -1:
            where_loop_histogram_wR = result_string.find('Final refinement wR =', where_loop_histogram)
            if where_loop_histogram_wR != -1:
                where_loop_histogram_wR_end = result_string.find('%', where_loop_histogram_wR)
                print(loop_histogram, result_string[ where_loop_histogram_wR : where_loop_histogram_wR_end + 1])
