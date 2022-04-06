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
from mantid.simpleapi import CreateWorkspace, logger
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
# d_spacing_min = 1.0
# refine_unit_cell = True
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
d_spacing_min = 1.0
refine_unit_cell = True


def call_subprocess(command_string):
    shell_output = subprocess.Popen([command_string.replace('"', '\\"')],
                                    shell=True,
                                    stdin=None,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    close_fds=True,
                                    universal_newlines=True)
    return shell_output.communicate()


def format_shell_output(title, shell_output_string):
    double_line = "-" * (len(title)+2) + "\n" + "-" * (len(title)+2)
    return "\n"*3 + double_line + "\n " + title + " \n" + double_line + "\n" + shell_output_string + double_line + "\n"*3


def find_in_file(file_path, marker_string, start_of_value, end_of_value, strip_separator=None):
    value_string = None
    with open(file_path, 'rt', encoding='utf-8') as file:
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
    with open(file_path, 'rt', encoding='utf-8') as file:
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


def read_basis(phase_file_path):
    basis_string = find_basis_block_in_file(phase_file_path, "atom", "\n", "loop")
    split_string_basis = basis_string.split()
    split_string_basis[0] = ''.join([i for i in split_string_basis[0] if not i.isdigit()])
    if len(split_string_basis[0]) == 2:
        split_string_basis[0] = ''.join([split_string_basis[0][0].upper(), split_string_basis[0][1].lower()])
    basis = " ".join(split_string_basis[0:6])  # assuming only one line in the basis block!!
    return basis


def read_space_group(phase_file_path):
    space_group = find_in_file(phase_file_path, "_symmetry_space_group_name_H-M", '"', '"', strip_separator='"')
    # Insert "-" to make the space group work
    for character_index, character in enumerate(space_group):
        if character.isdigit():
            split_string = list(space_group)
            split_string.insert(character_index, "-")
            space_group = "".join(split_string)
    return space_group


def choose_cell_lengths(overriding_cell_lengths, phase_file_path):
    if overriding_cell_lengths:
        cell_lengths = " ".join([str(overriding_cell_lengths[0]),
                                 str(overriding_cell_lengths[1]),
                                 str(overriding_cell_lengths[2])])
    else:
        cell_length_a = find_in_file(phase_file_path, '_cell_length_a', ' ', '_cell_length_b')
        cell_length_b = find_in_file(phase_file_path, '_cell_length_b', ' ', '_cell_length_c')
        cell_length_c = find_in_file(phase_file_path, '_cell_length_c', ' ', '_cell')
        cell_lengths = " ".join([cell_length_a, cell_length_b, cell_length_c])
    return cell_lengths


def create_pawley_reflections(cell_lengths, space_group, basis, dmin):
    structure = CrystalStructure(cell_lengths, space_group, basis)
    generator = ReflectionGenerator(structure)

    hkls = generator.getUniqueHKLsUsingFilter(dmin, 3.0, ReflectionConditionFilter.StructureFactor)
    dValues = generator.getDValues(hkls)
    pg = structure.getSpaceGroup().getPointGroup()
    # Make list of tuples and sort by d-values, descending, include point group for multiplicity.
    generated_reflections = sorted([[list(hkl), d, len(pg.getEquivalents(hkl))] for hkl, d in zip(hkls, dValues)],
                                   key=lambda x: x[1] - x[0][0]*1e-6, reverse=True)
    return generated_reflections


def check_for_output_file(temp_save_directory, name_of_project, file_extension, file_descriptor, gsas_error_string):
    gsas_output_filename = name_of_project + file_extension
    if gsas_output_filename not in os.listdir(temporary_save_directory):
        raise FileNotFoundError(f"GSAS-II call must have failed, as the output {file_descriptor} file was not found.",
                                format_shell_output(title="Errors from GSAS-II", shell_output_string=gsas_error_string))
    return os.path.join(temp_save_directory, gsas_output_filename)


def chop_to_limits(input_array, x, min_x, max_x):
    input_array[x <= min_x] = np.nan
    input_array[x >= max_x] = np.nan
    return input_array


def read_gsas_lst_and_print_wR(result_filepath):
    with open(result_filepath, 'rt', encoding='utf-8') as file:
        result_string = file.read().replace('\n', '')
        for loop_histogram in data_files:
            where_loop_histogram = result_string.rfind(loop_histogram)
            if where_loop_histogram != -1:
                where_loop_histogram_wR = result_string.find('Final refinement wR =', where_loop_histogram)
                if where_loop_histogram_wR != -1:
                    where_loop_histogram_wR_end = result_string.find('%', where_loop_histogram_wR)
                    logger.notice(loop_histogram)
                    logger.notice(result_string[where_loop_histogram_wR: where_loop_histogram_wR_end + 1])


def load_gsas_histogram(temp_save_directory, name_of_project, histogram_index, min_x, max_x):
    result_csv = os.path.join(temp_save_directory, name_of_project + f"_{histogram_index}.csv")
    my_data = np.transpose(np.genfromtxt(result_csv, delimiter=",", skip_header=39))
    # x  y_obs	weight	y_calc	y_bkg	Q
    x_values = my_data[0]
    y_obs = chop_to_limits(np.array(my_data[1]), x_values,  min_x[histogram_index], max_x[histogram_index])
    y_calc = chop_to_limits(np.array(my_data[3]), x_values, min_x[histogram_index], max_x[histogram_index])
    y_diff = y_obs - y_calc
    y_diff -= np.max(np.ma.masked_invalid(y_diff))
    y_bkg = chop_to_limits(np.array(my_data[4]), x_values, min_x[histogram_index], max_x[histogram_index])
    y_data = np.concatenate((y_obs, y_calc, y_diff, y_bkg))

    gsas_histogram = CreateWorkspace(OutputWorkspace=f"gsas_histogram_{histogram_index}",
                                     DataX=np.tile(my_data[0], 4), DataY=y_data, NSpec=4)
    return gsas_histogram


def load_gsas_reflections(temp_save_directory, name_of_project, histogram_index):
    result_reflections_txt = os.path.join(temp_save_directory, name_of_project + f"_reflections_{histogram_index}_Fe_gamma.txt")
    return np.loadtxt(result_reflections_txt)


def plot_gsas_histogram(gsas_histogram, reflection_positions, name_of_project, histogram_index, min_x, max_x):
    fig, axes = plt.subplots(num=name_of_project + f'_{histogram_index} GSAS-II Refinement',
                             subplot_kw={'projection': 'mantid'})
    axes.plot(gsas_histogram, color='#1105f0', label='observed', linestyle='None', marker='+', wkspIndex=0)
    axes.plot(gsas_histogram, color='#246b01', label='calculated', wkspIndex=1)
    axes.plot(gsas_histogram, color='#09acb8', label='difference', wkspIndex=2)
    axes.plot(gsas_histogram, color='#ff0000', label='background', wkspIndex=3)
    axes.set_title(project_name + ' GSAS Refinement')
    _, y_max = axes.get_ylim()
    axes.plot(reflection_positions, [-0.10*y_max]*len(reflection_positions),
              color='#1105f0', label='reflections', linestyle='None', marker='|', mew=1.5, ms=8)
    axes.axvline(min_x[histogram_index], color='#246b01', linestyle='--')
    axes.axvline(max_x[histogram_index], color='#a80000', linestyle='--')
    axes.legend(fontsize=8.0).draggable().legend
    plt.show()


''' Pre exec calculations '''
refine_histogram_scale_factor = True  # True by default

user_save_directory = os.path.join(save_directory, project_name)
temporary_save_directory = os.path.join(save_directory,
                                        datetime.datetime.now().strftime('tmp_EngDiff_GSASII_%Y-%m-%d_%H-%M-%S'))
make_temporary_save_directory = ("mkdir -p " + temporary_save_directory)
out_make_temporary_save_directory, err_make_temporary_save_directory = call_subprocess(make_temporary_save_directory)

if refinement_method == 'Pawley':
    phase_filepath = os.path.join(data_directory, phase_files[0])
    mantid_pawley_reflections = create_pawley_reflections(cell_lengths=choose_cell_lengths(override_cell_lengths,
                                                                                           phase_filepath),
                                                          space_group=read_space_group(phase_filepath),
                                                          basis=read_basis(phase_filepath),
                                                          dmin=d_spacing_min)

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
    raise ValueError(f"No Pawley Reflections were generated for the phases provided. Not calling GSAS-II.")

if len(instrument_files) != 1 or len(instrument_files) != number_histograms:
    raise ValueError(f'The number of instrument files ({len(instrument_files)}) must be 1 '
                     f'or equal to the number of input histograms {number_histograms}')


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
    override_cell_lengths=override_cell_lengths,
    refine_unit_cell=refine_unit_cell,
    d_spacing_min=d_spacing_min
)

call_gsas2 = (path_to_gsas2 + "bin/python "
              + "/home/danielmurphy/mantid/qt/python/mantidqtinterfaces/mantidqtinterfaces/Engineering/gui/"
              + "engineering_diffraction/tabs/gsas2/call_G2sc.py "
              + parse_inputs.Gsas2Inputs_to_json(gsas2_inputs)
              )

start = time.time()
out_call_gsas2, err_call_gsas2 = call_subprocess(call_gsas2)
gsas_runtime = time.time() - start
logger.notice(format_shell_output(title="Commandline output from GSAS-II", shell_output_string=out_call_gsas2))

gsas_project_filepath = check_for_output_file(temporary_save_directory, project_name, ".gpx", "project file", err_call_gsas2)
gsas_result_filepath = check_for_output_file(temporary_save_directory, project_name, ".lst", "result", err_call_gsas2)
logger.notice(f"\nGSAS-II call complete in {gsas_runtime} seconds.\n")

logger.notice(f"GSAS-II .lst result file found. Opening {project_name}.lst")
read_gsas_lst_and_print_wR(gsas_result_filepath)

for index_histograms in range(number_histograms):
    gsas_histogram_workspace = load_gsas_histogram(temporary_save_directory, project_name, index_histograms, x_min, x_max)
    reflections = load_gsas_reflections(temporary_save_directory, project_name, index_histograms)
    plot_gsas_histogram(gsas_histogram_workspace, reflections, project_name, index_histograms, x_min, x_max)

make_user_save_directory = ("mkdir -p " + user_save_directory)
out_make_user_save_directory, err_make_user_save_directory = call_subprocess(make_user_save_directory)
if err_make_user_save_directory:
    raise ValueError("Could not create the user save directory {user_save_directory}")

for output_file in os.listdir(temporary_save_directory):
    os.rename(os.path.join(temporary_save_directory, output_file),
              os.path.join(user_save_directory, output_file))
os.rmdir(temporary_save_directory)
logger.notice(f"\n\nOutput GSAS-II files saved in {user_save_directory}")

# # open GSAS-II project
# open_project_call = (path_to_gsas2 + "bin/python " + path_to_gsas2 + "GSASII/GSASII.py "
#                      + os.path.join(user_save_directory, project_name + ".gpx"))
# out_open_project_call, err_open_project_call = call_subprocess(open_project_call)
