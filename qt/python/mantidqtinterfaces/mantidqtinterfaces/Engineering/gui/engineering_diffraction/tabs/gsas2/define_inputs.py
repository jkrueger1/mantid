import subprocess
import time
import numpy as np
import matplotlib.pyplot as plt
from mantid.simpleapi import *
from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter

TIME_DELAY = 10  # seconds within gsas files must have been generated

'''Inputs Tutorial'''
# path_to_gsas2 = "/home/danielmurphy/gsas2/"
# save_directory = "/home/danielmurphy/Downloads/GSASdata/new_outputs/"
# data_directory = "/home/danielmurphy/Downloads/GSASdata/"
# refinement_method = "Rietveld"
# input_data_files = ["PBSO4.XRA", "PBSO4.CWN"]
# histogram_indexing = []
# instrument_files = ["INST_XRY.PRM","inst_d1a.prm"]
# phase_files = ["PbSO4-Wyckoff.cif"]
# project_name = "mantid_test"
#
# x_min = [16.0, 19.0]
# x_max = [158.4, 153.0]
#
# user_override_cell_length = [3.65, 3.65, 3.65] # in presenter force to empty or len3 list of floats
# refine_background = True
# refine_microstrain = True
# refine_sigma_one = False
# refine_gamma = False
#
# refine_histogram_scale_factor = True  # True by default

'''Inputs Mantid'''
path_to_gsas2 = "/home/danielmurphy/gsas2/"
save_directory = "/home/danielmurphy/Downloads/GSASdata/new_outputs/"
data_directory = "/home/danielmurphy/Desktop/GSASMantiddata_030322/"
refinement_method = "Pawley"
input_data_files = ["Save_gss_305761_307521_bank_1_bgsub.gsa"]  # ["ENGINX_305761_307521_all_banks_TOF.gss"]
histogram_indexing = [1]  # assume only indexing when using 1 histogram file
instrument_files = ["ENGINX_305738_bank_1.prm"]
phase_files = ["FE_GAMMA.cif"]  # ["7217887.cif"]
project_name = "220321script3"

x_min = [18401.18]
x_max = [50000.0]

user_override_cell_length = [3.65, 3.65, 3.65] # in presenter force to empty or len3 list of floats
refine_background = False
refine_microstrain = True
refine_sigma_one = False
refine_gamma = False

refine_histogram_scale_factor = True  # True by default

'''Generate Pawley Reflections'''
dmin = 1.0


def find_in_file(file_path, marker_string, start_of_value, end_of_value, strip_separator=None):
    value_string = None
    with open(file_path, 'r') as file:
        full_file_string = file.read().replace('\n', '')
        where_marker = full_file_string.rfind(marker_string)
        if where_marker != -1:
            where_value_start = full_file_string.find(start_of_value, where_marker)
            if where_value_start != -1:
                where_value_end = full_file_string.find(end_of_value, where_value_start + 1)
                value_string = full_file_string[where_value_start: where_value_end]
                if strip_separator:
                    value_string = value_string.strip(strip_separator + " ")
                else:
                    value_string = value_string.strip(" ")
    return value_string


def find_basis_block_in_file(file_path, marker_string, start_of_value, end_of_value):

    with open(file_path, 'r') as file:
        full_file_string = file.read()
        where_marker = full_file_string.find(marker_string)
        value_string = None
        if where_marker != -1:
            where_first_digit = -1
            index = int(where_marker)
            while index < len(full_file_string):
                if full_file_string[index].isdecimal():
                    where_first_digit = index
                    break
                index += 1
            if where_first_digit != -1:
                where_start_of_line = full_file_string.rfind(start_of_value, 0, where_first_digit-1)
                if where_start_of_line != -1:
                    where_end_of_block = full_file_string.find(end_of_value, where_start_of_line)
                    # if "loop" not found then assume the end of the file is the end of this block
                    value_string = full_file_string[where_start_of_line: where_end_of_block]
                    for remove_string in ["Biso", "Uiso"]:
                        value_string = value_string.replace(remove_string, "")
    return value_string


compressed_reflections = []
if refinement_method == 'Pawley':
    phase_filepath = os.path.join(data_directory, phase_files[0])

    basis = find_basis_block_in_file(phase_filepath, "atom", "\n", "loop")
    split_string_basis = basis.split()
    split_string_basis[0] = ''.join([i for i in split_string_basis[0] if not i.isdigit()])
    if len(split_string_basis[0]) == 2:
        split_string_basis[0] = ''.join([split_string_basis[0][0].upper(), split_string_basis[0][1].lower()])

    basis = " ".join(split_string_basis[0:6]) # assuming only one line in the basis block!!
    print("Basis of Scatterers:\n", "-"*30, "\n", basis, "\n", "-"*30, "\n")

    space_group = find_in_file(phase_filepath, "_symmetry_space_group_name_H-M", '"', '"', strip_separator='"')
    # Insert "-" to make the space group work with Mantid
    for index, character in enumerate(space_group):
        if character.isdigit():
            split_string = list(space_group)
            split_string.insert(index, "-")
            space_group = "".join(split_string)
    print("Space Group: ", space_group)

    if user_override_cell_length:
        cell_lengths = " ".join([str(user_override_cell_length[0]),
                                 str(user_override_cell_length[1]),
                                 str(user_override_cell_length[2])])
    else:
        cell_length_a = find_in_file(phase_filepath, '_cell_length_a', ' ', '_cell_length_b')
        cell_length_b = find_in_file(phase_filepath, '_cell_length_b', ' ', '_cell_length_c')
        cell_length_c = find_in_file(phase_filepath, '_cell_length_c', ' ', '_cell')
        cell_lengths = " ".join([cell_length_a, cell_length_b, cell_length_c])
    print("Cell Lengths:", cell_lengths)

    structure = CrystalStructure(cell_lengths, space_group, basis)
    generator = ReflectionGenerator(structure)

    hkls = generator.getUniqueHKLsUsingFilter(dmin, 3.0, ReflectionConditionFilter.StructureFactor)
    dValues = generator.getDValues(hkls)
    fSquared = generator.getFsSquared(hkls)
    pg = structure.getSpaceGroup().getPointGroup()
    # Make list of tuples and sort by d-values, descending, include point group for multiplicity.
    reflections = sorted([(hkl, d, fsq, len(pg.getEquivalents(hkl))) for hkl, d, fsq in zip(hkls, dValues, fSquared)],
                                    key=lambda x: x[1] - x[0][0]*1e-6, reverse=True)
    # 'HKL', 'd', 'F^2', 'M'

    for reflection in reflections:
        reflection = list(reflection)
        for index, elem in enumerate(reflection):
            reflection[index] = str(elem)
        compressed_reflections.append("#".join(reflection))


'''Matching Dictionary'''
number_of_inputs = {'data_files': len(input_data_files), 'histogram_indices': len(histogram_indexing),
                    'phases': len(phase_files), 'instruments': len(instrument_files), 'limits': len(x_min),
                    'Pawley Reflections': len(compressed_reflections),
                    "Override Cell Length": len(user_override_cell_length)}


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

if refinement_method == 'Pawley' and not compressed_reflections:
    raise ValueError(f"No Pawley Reflections were generated for the phases provided. Not calling GSASII.")


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
             + str(number_of_inputs['limits']) + " "
             + str(number_of_inputs['Pawley Reflections']) + " "
             + str(number_of_inputs['Override Cell Length']) + " "
             + str(1 if refine_background else 0) + " "
             + str(1 if refine_microstrain else 0) + " "
             + str(1 if refine_sigma_one else 0) + " "
             + str(1 if refine_gamma else 0) + " "
             + str(1 if refine_histogram_scale_factor else 0) + " ")

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

if compressed_reflections:
    for reflection in compressed_reflections:
        main_call += (reflection + " ")

if user_override_cell_length:
    for cell_length in user_override_cell_length:
        main_call += (str(cell_length) + " ")

start = time.time()

p = subprocess.Popen(
        [main_call],
        shell=True,
        stdin=None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        close_fds=True)
out, err = p.communicate()
double_line ="-"*33 + "\n" + "-"*33
# if err:
#     print("\n"*5 + double_line + "\n Error from GSAS-II \n" + double_line + "\n" )
#     raise ValueError(err.decode(), out.decode())
out = out.decode()
print("\n"*5 + double_line + "\n Commandline output from GSAS-II \n" + double_line + "\n" + out + double_line + "\n"*5)

gsas_runtime = time.time() - start
print(f"\nGSASII Complete in {gsas_runtime} seconds.\n")


# project_path = save_directory + project_name + '.gpx'
result_filepath = save_directory + project_name + '.lst'

try:
    last_modified_time = os.path.getmtime(result_filepath)
except FileNotFoundError:
    raise FileNotFoundError('This GSASII operation must have failed as the output .lst file was not found.')
for index in range(number_histograms):
    result_csv = save_directory + project_name + f"_{index}.csv"
    try:
        last_modified_time_csv = os.path.getmtime(result_csv)
    except FileNotFoundError:
        raise FileNotFoundError('This GSASII operation must have failed as the exported histogram .csv file/s was not found.')

if time.time() > (last_modified_time + TIME_DELAY) or time.time() > (last_modified_time_csv + TIME_DELAY):
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


def chop_to_limits(input_array, xmin, xmax):
    input_array[x <= xmin] = np.nan
    input_array[x >= xmax] = np.nan
    return input_array


for index in range(number_histograms):
    result_csv = save_directory + project_name + f"_{index}.csv"
    my_data = np.transpose(np.genfromtxt(result_csv, delimiter=",", skip_header=39))
    # x  y_obs	weight	y_calc	y_bkg	Q
    x = my_data[0]
    y_obs = chop_to_limits(np.array(my_data[1]), x_min[index], x_max[index])
    y_calc = chop_to_limits(np.array(my_data[3]), x_min[index], x_max[index])
    y_diff = y_obs - y_calc
    y_diff -= np.max(np.ma.masked_invalid(y_diff))
    y_bkg = chop_to_limits(np.array(my_data[4]), x_min[index], x_max[index])
    y_data = np.concatenate((y_obs, y_calc, y_diff, y_bkg))

    result_reflections_txt = os.path.join(save_directory, project_name + f"_reflections_{index}.txt")
    reflection_positions = np.loadtxt(result_reflections_txt)

    gsas_histogram = CreateWorkspace(DataX=np.tile(my_data[0], 4), DataY=y_data, NSpec=4)

    fig, axes = plt.subplots(num=project_name + f'_{index} GSASII Refinement', subplot_kw={'projection': 'mantid'})
    axes.plot(gsas_histogram, color='#1105f0', label='observed', linestyle='None', marker='+', wkspIndex=0)
    axes.plot(gsas_histogram, color='#246b01', label='calculated', wkspIndex=1)
    axes.plot(gsas_histogram, color='#09acb8', label='difference', wkspIndex=2)
    axes.plot(gsas_histogram, color='#ff0000', label='background', wkspIndex=3)
    axes.set_title(project_name + ' GSAS Refinement')
    _, y_max = axes.get_ylim()
    axes.plot(reflection_positions, [-0.10*y_max]*len(reflection_positions),
              color='#1105f0', label='reflections', linestyle='None', marker='|', mew=1.5, ms=8)
    axes.axvline(x_min[index], color='#246b01', linestyle='--')
    axes.axvline(x_max[index], color='#a80000', linestyle='--')
    legend = axes.legend(fontsize=8.0).draggable().legend
    plt.show()


open_project_call = (path_to_gsas2 + "bin/python " + path_to_gsas2 + "GSASII/GSASII.py "
                     + save_directory + project_name + ".gpx")
p = subprocess.Popen([open_project_call],
                     shell=True,
                     stdin=None,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     close_fds=True)
out, err = p.communicate()
