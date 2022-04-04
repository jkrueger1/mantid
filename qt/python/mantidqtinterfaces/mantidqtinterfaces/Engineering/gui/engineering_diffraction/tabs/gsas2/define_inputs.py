# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
import os
import subprocess
import time
import datetime

import matplotlib.pyplot as plt
import numpy as np
from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter
from mantid.simpleapi import CreateWorkspace
from mantidqtinterfaces.Engineering.gui.engineering_diffraction.tabs.gsas2 import parse_inputs


'''Inputs Tutorial'''
# path_to_gsas2 = "/home/danielmurphy/gsas2/"
# save_directory = "/home/danielmurphy/Downloads/GSASdata/new_outputs/"
# data_directory = "/home/danielmurphy/Downloads/GSASdata/"
# refinement_method = "Rietveld"
# data_files = ["PBSO4.XRA", "PBSO4.CWN"]
# histogram_indexing = []
# instrument_files = ["INST_XRY.PRM","inst_d1a.prm"]
# phase_files = ["PbSO4-Wyckoff.cif"]
# project_name = "mantid_test"
#
# x_min = [16.0, 19.0]
# x_max = [158.4, 153.0]
#
# override_cell_lengths = [3.65, 3.65, 3.65] # in presenter force to empty or len3 list of floats
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
data_files = ["Save_gss_305761_307521_bank_1_bgsub.gsa"]  # ["ENGINX_305761_307521_all_banks_TOF.gss"]
histogram_indexing = [1]  # assume only indexing when using 1 histogram file
instrument_files = ["ENGINX_305738_bank_1.prm"]
phase_files = ["FE_GAMMA.cif"]
project_name = "220321script3"

x_min = [18401.18]
x_max = [50000.0]

override_cell_lengths = [3.65, 3.65, 3.65] # in presenter force to empty or len3 list of floats
refine_background = True
refine_microstrain = True
refine_sigma_one = False
refine_gamma = False


def call_subprocess(command_string):
    shell_output = subprocess.Popen([command_string.replace('"','\\"')],
                                    shell=True,
                                    stdin=None,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    close_fds=True)
    return shell_output.communicate()


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


def read_phase_and_create_reflections(data_directory, phase_files, override_cell_lengths):

    # ASSUMING phase files has only 1 entry
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

    if override_cell_lengths:
        cell_lengths = " ".join([str(override_cell_lengths[0]),
                                 str(override_cell_lengths[1]),
                                 str(override_cell_lengths[2])])
    else:
        cell_length_a = find_in_file(phase_filepath, '_cell_length_a', ' ', '_cell_length_b')
        cell_length_b = find_in_file(phase_filepath, '_cell_length_b', ' ', '_cell_length_c')
        cell_length_c = find_in_file(phase_filepath, '_cell_length_c', ' ', '_cell')
        cell_lengths = " ".join([cell_length_a, cell_length_b, cell_length_c])
    print("Cell Lengths:", cell_lengths)

    structure = CrystalStructure(cell_lengths, space_group, basis)
    generator = ReflectionGenerator(structure)

    dmin = 1.0
    hkls = generator.getUniqueHKLsUsingFilter(dmin, 3.0, ReflectionConditionFilter.StructureFactor)
    dValues = generator.getDValues(hkls)
    pg = structure.getSpaceGroup().getPointGroup()
    # Make list of tuples and sort by d-values, descending, include point group for multiplicity.
    reflections = sorted([[list(hkl), d, len(pg.getEquivalents(hkl))] for hkl, d in zip(hkls, dValues)],
                                    key=lambda x: x[1] - x[0][0]*1e-6, reverse=True)
    print(type([reflections[0][0]]))
    return reflections


def chop_to_limits(input_array, x, xmin, xmax):
    input_array[x <= xmin] = np.nan
    input_array[x >= xmax] = np.nan
    return input_array


def read_gsas_lst_and_print_wR(result_filepath):
    with open(result_filepath, 'r') as file:
        result_string = file.read().replace('\n', '')
        for loop_histogram in data_files:
            where_loop_histogram = result_string.rfind(loop_histogram)
            if where_loop_histogram != -1:
                where_loop_histogram_wR = result_string.find('Final refinement wR =', where_loop_histogram)
                if where_loop_histogram_wR != -1:
                    where_loop_histogram_wR_end = result_string.find('%', where_loop_histogram_wR)
                    print(loop_histogram, result_string[where_loop_histogram_wR: where_loop_histogram_wR_end + 1])


def load_and_plot_gsas_histograms(save_directory, project_name, index, x_min, x_max):
    result_csv = save_directory + project_name + f"_{index}.csv"
    my_data = np.transpose(np.genfromtxt(result_csv, delimiter=",", skip_header=39))
    # x  y_obs	weight	y_calc	y_bkg	Q
    x_values = my_data[0]
    y_obs = chop_to_limits(np.array(my_data[1]), x_values,  x_min[index], x_max[index])
    y_calc = chop_to_limits(np.array(my_data[3]), x_values, x_min[index], x_max[index])
    y_diff = y_obs - y_calc
    y_diff -= np.max(np.ma.masked_invalid(y_diff))
    y_bkg = chop_to_limits(np.array(my_data[4]), x_values, x_min[index], x_max[index])
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
    axes.legend(fontsize=8.0).draggable().legend
    plt.show()


def format_shell_output(title, shell_output_string):
    double_line = "-" * (len(title)+2) + "\n" + "-" * (len(title)+2)
    return ("\n"*3 + double_line + "\n " + title + " \n" + double_line + "\n" + shell_output_string.decode() + double_line + "\n"*3)


''' Pre exec calculations '''
refine_histogram_scale_factor = True  # True by default

user_save_directory = os.path.join(save_directory, project_name)
temporary_save_directory = os.path.join(save_directory,
                                        datetime.datetime.now().strftime('tmp_EngDiff_GSASII_%Y-%m-%d_%H-%M-%S'))
make_temporary_save_directory = ("mkdir -p " + temporary_save_directory)
out_make_temporary_save_directory, err_make_temporary_save_directory = call_subprocess(make_temporary_save_directory)

if refinement_method == 'Pawley':
    mantid_pawley_reflections = read_phase_and_create_reflections(data_directory, phase_files, override_cell_lengths)


'''Validation'''
number_histograms = len(data_files)
if histogram_indexing and len(data_files) == 1:
    number_histograms = len(histogram_indexing)
if len(x_min) != 0 and len(x_min) != number_histograms:
    raise ValueError(f"The number of x_min values ({len(x_min)}) must equal the"
                     + f"number of histograms ({number_histograms}))")

if histogram_indexing and len(data_files) > 1:
    raise ValueError(f"Histogram indexing is currently only supported, when the "
                     + f"number of data_files ({len(data_files)}) == 1")

if refinement_method == 'Pawley' and not mantid_pawley_reflections:
    raise ValueError(f"No Pawley Reflections were generated for the phases provided. Not calling GSASII.")


'''exec'''

gsas2_inputs = parse_inputs.Gsas2Inputs(
    path_to_gsas2=path_to_gsas2,
    temporary_save_directory=temporary_save_directory,
    data_directory=data_directory,
    project_name=project_name,
    refinement_method=refinement_method,
    refine_background=refine_background,
    refine_microstrain=refine_microstrain,
    refine_sigma_one=refine_sigma_one,
    refine_gamma=refine_gamma,
    refine_histogram_scale_factor=refine_histogram_scale_factor,
    data_files=data_files,
    histogram_indexing=histogram_indexing,
    phase_files=phase_files,
    instrument_files=instrument_files,
    limits=[x_min, x_max],
    mantid_pawley_reflections=mantid_pawley_reflections,
    override_cell_lengths=override_cell_lengths
)

call_gsas2 = (path_to_gsas2 + "bin/python "
              + "/home/danielmurphy/mantid/qt/python/mantidqtinterfaces/mantidqtinterfaces/Engineering/gui/"
              + "engineering_diffraction/tabs/gsas2/call_G2sc.py "
              + parse_inputs.convert_Gsas2Inputs_to_json(gsas2_inputs)
              )

start = time.time()
out_call_gsas2, err_call_gsas2 = call_subprocess(call_gsas2)
gsas_runtime = time.time() - start
print(format_shell_output(title="Commandline output from GSAS-II", shell_output_string=out_call_gsas2))


gsas_project_filename = project_name + '.gpx'
gsas_project_filepath = os.path.join(temporary_save_directory, gsas_project_filename)
if gsas_project_filename not in os.listdir(temporary_save_directory):
    raise FileNotFoundError("GSAS-II call must have failed, as the output project file was not found.",
                            format_shell_output(title="Errors from GSAS-II", shell_output_string=err_call_gsas2))

gsas_log_filename = project_name + '.lst'
gsas_log_filepath = os.path.join(temporary_save_directory, gsas_log_filename)
if gsas_log_filename not in os.listdir(temporary_save_directory):
    raise FileNotFoundError("GSAS-II call must have failed, as the output log file was not found.",
                            format_shell_output(title="Errors from GSAS-II", shell_output_string=err_call_gsas2))

print(f"\nGSAS-II call complete in {gsas_runtime} seconds.\n")
print(f"GSAS-II .lst result file found. Opening {gsas_log_filename}")
read_gsas_lst_and_print_wR(gsas_log_filepath)

for index in range(number_histograms):
    load_and_plot_gsas_histograms(save_directory, project_name, index, x_min, x_max)

make_user_save_directory = ("mkdir -p " + user_save_directory)
out_make_user_save_directory, err_make_user_save_directory = call_subprocess(make_user_save_directory)
if err_make_user_save_directory:
    raise ValueError("Could not create the user save directory {user_save_directory}")

if os.listdir(temporary_save_directory): # currently just checking if there is something in the temp save directory
    for output_file in os.listdir(temporary_save_directory):
        os.rename(os.path.join(temporary_save_directory, output_file),
                  os.path.join(user_save_directory, output_file))
    os.rmdir(temporary_save_directory)
    print(f"\n\nOutput GSAS-II files saved in {user_save_directory}")

# # open GSASII project
# open_project_call = (path_to_gsas2 + "bin/python " + path_to_gsas2 + "GSASII/GSASII.py "
#                      + os.path.join(user_save_directory, project_name + ".gpx"))
# out_open_project_call, err_open_project_call = call_subprocess(open_project_call)
