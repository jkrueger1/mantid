import os
import time


'''Inputs Tutorial'''
# path_to_gsas2 = "/home/danielmurphy/gsas2/"
# save_directory = "/home/danielmurphy/Downloads/GSASdata/new_outputs/"
# data_directory = "/home/danielmurphy/Downloads/GSASdata/"
# refinement_method = "Pawley"
# input_data_files = ["PBSO4.XRA", "PBSO4.CWN"]
# histogram_indexing = []
# instrument_files = ["INST_XRY.PRM","inst_d1a.prm"]
# phase_files = ["PbSO4-Wyckoff.cif"]
# project_name = "mantid_test"
#
# x_min = [16.0,19.0]
# x_max = [158.4,153.0]
#

'''Inputs Mantid'''
path_to_gsas2 = "/home/danielmurphy/gsas2/"
save_directory = "/home/danielmurphy/Downloads/GSASdata/new_outputs/"
data_directory = "/home/danielmurphy/Desktop/GSASMantiddata_030322/"
refinement_method = "Pawley"
input_data_files = ["ENGINX_305761_307521_all_banks_TOF.gss"]  # ["ENGINX_305761_307521_bank_1_TOF.gss"]
histogram_indexing = [1]  # assume only indexing when using 1 histogram file
instrument_files = ["ENGINX_305738_bank_1.prm"]
phase_files = ["7217887.cif"]
project_name = "mantid_enginx1"

x_min = []
x_max = []

'''Matching Dictionary'''
number_of_inputs = {'data_files': len(input_data_files), 'histogram_indices': len(histogram_indexing),
                    'phases': len(phase_files), 'instruments': len(instrument_files), 'limits': len(x_min)}


'''Validation'''
if len(x_min) != len(x_max):
    raise ValueError(f"The number of x_min values ({number_of_inputs['x_min']}) must equal the"
                     + f"number of x_max values ({number_of_inputs['x_max']})")

number_histograms = len(input_data_files)
if histogram_indexing and len(input_data_files) == 1:
    number_histograms = len(histogram_indexing)
if number_of_inputs['limits'] != 0 and number_of_inputs['limits'] != number_histograms:
    raise ValueError(f"The number of x_min values ({number_of_inputs['x_min']}) must equal the"
                     + f"number of histograms ({number_histograms}))")

if histogram_indexing and len(input_data_files) > 1:
    raise ValueError(f"Histogram indexing can is currently only supported, when the "
                     + f"number of input_data_files ({number_of_inputs['histograms']}) == 1")


'''exec'''
main_call = (path_to_gsas2 + "bin/python "
             + "/home/danielmurphy/mantid/qt/python/mantidqtinterfaces/mantidqtinterfaces/Engineering/gui/"
             + "engineering_diffraction/tabs/gsas2/call_G2sc.py "
             + path_to_gsas2 + " "
             + save_directory + " "
             + data_directory + " "
             + refinement_method + " "
             + project_name + " "
             + str(number_of_inputs['data_files']) + " "
             + str(number_of_inputs['histogram_indices']) + " "
             + str(number_of_inputs['phases']) + " "
             + str(number_of_inputs['instruments']) + " "
             + str(number_of_inputs['limits']) + " ")

for input_data_file in input_data_files:
    main_call += (input_data_file + " ")
for histogram_index in histogram_indexing:
    main_call += (str(histogram_index) + " ")
for phase in phase_files:
    main_call += (phase + " ")
for instrument in instrument_files:
    main_call += (instrument + " ")

if x_min and x_max:
    for value in x_min:
        main_call += (str(value) + " ")
    for value in x_max:
        main_call += (str(value) + " ")

start = time.time()
os.system(main_call)
gsas_runtime = time.time() - start
print(f"\nGSASII Complete in {gsas_runtime} seconds.\n")


# project_path = save_directory + project_name + '.gpx'
result_filepath = save_directory + project_name + '.lst'

try:
    last_modified_time = os.path.getmtime(result_filepath)
except FileNotFoundError:
    raise FileNotFoundError('This GSASII operation must have failed as the output files were not found.')

if time.time() > (last_modified_time + 2):
    # if GSASII result file not modified in the last 2 seconds
    m_ti = time.ctime(last_modified_time)
    # Using the timestamp string to create a time object/structure
    t_obj = time.strptime(m_ti)
    # Transforming the time object to a timestamp of ISO 8601 format
    T_stamp = time.strftime("%Y-%m-%d %H:%M:%S", t_obj)
    print(f"The file located at the path {result_filepath} was last modified at {T_stamp}")
    raise ValueError('This GSASII operation must have failed as the output files were from an old operation.')

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
