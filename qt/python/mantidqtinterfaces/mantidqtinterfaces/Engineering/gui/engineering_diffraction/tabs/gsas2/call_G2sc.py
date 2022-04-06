# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
import os
import sys
import json


def print_histogram_R_factors(project):
    """prints profile rfactors for all histograms"""
    print(u"*** Weighted profile R-factor: Rwp " + os.path.split(project.filename)[1])
    for loop_histogram in project.histograms():
        print("\t{:20s}: {:.2f}".format(loop_histogram.name, loop_histogram.get_wR()))
    print("")


def calculate_number_histograms(number_data_files, number_histogram_indices):
    number_of_histograms = number_data_files
    if number_histogram_indices and number_data_files == 1:
        number_of_histograms = number_histogram_indices
    elif number_histogram_indices > 1 and number_data_files > 1:
        raise ValueError("Cannot use histogram indexing with multiple input data files.")
    return number_of_histograms


def assign_instruments_to_histograms(input_instrument_files, number_of_histograms):
    iparams_input = [input_instrument_files[0]] * number_histograms
    if len(input_instrument_files) > 1 and len(input_instrument_files) == number_of_histograms:
        iparams_input = input_instrument_files
    return iparams_input


def add_histograms(data_filenames, histogram_indices, project, input_data_directory, instruments):
    if not histogram_indexing:
        for loop_index, loop_data_filename in enumerate(data_filenames):
            project.add_powder_histogram(datafile=os.path.join(input_data_directory,
                                                               loop_data_filename),
                                         iparams=os.path.join(data_directory,
                                                              instruments[loop_index]),
                                         phases=project.phases()
                                         )
    else:
        for loop_index, loop_histogram_index in enumerate(histogram_indices):
            project.add_powder_histogram(datafile=os.path.join(input_data_directory,
                                                               data_files[0]),
                                         iparams=os.path.join(data_directory,
                                                              instruments[loop_index]),
                                         phases=project.phases(),
                                         databank=loop_histogram_index  # indexing starts at 1
                                         )


def add_pawley_reflections(pawley_reflections, project, d_min):
    for loop_gsas_phase in project.phases():
        loop_gsas_phase.data['General']['doPawley'] = True
        gsas_reflections = []
        for reflection in pawley_reflections:
            [h, k, l], d, multiplicity = reflection
            gsas_reflections.append([int(h), int(k), int(l), int(multiplicity), float(d), True, 100.0, d_min])
        loop_gsas_phase.data["Pawley ref"] = gsas_reflections


def set_max_number_cycles(number_cycles):
    gsas_project.data['Controls']['data']['max cyc'] = number_cycles


def enable_background(refine, project):
    for loop_histogram in project.histograms():
        if refine:
            loop_histogram.set_refinements({"Background": {"no. coeffs": 3, "refine": True}})
        else:
            loop_histogram.set_refinements({"Background": {"no. coeffs": 0, "refine": False}})


def enable_histogram_scale_factor(refine, project):
    # GSAS-II default is enabled
    if not refine:
        for loop_histogram in project.histograms():
            loop_histogram.SampleParameters["Scale"] = [1.0, False]


def enable_unit_cell(refine, override_unit_cell_lengths, project):
    for loop_phase in project.phases():
        if refine:
            loop_phase.set_refinements({"Cell": True})
        if override_unit_cell_lengths:
            loop_phase.data['General']['Cell'][1:4] = tuple(override_unit_cell_lengths)


def enable_limits(x_limits, project):
    x_min, x_max = x_limits
    if x_min and x_max:
        for loop_index, loop_histogram in enumerate(project.histograms()):
            loop_histogram.set_refinements({'Limits': [x_min[loop_index], x_max[loop_index]]})


def run_microstrain_refinement(refine, project, path_to_project):
    if refine:
        for loop_phase in project.phases():
            loop_phase.set_HAP_refinements({'Mustrain': {'type': 'isotropic', 'refine': True}})
        project.do_refinements()
        project.save(path_to_project)
        print_histogram_R_factors(project)


def run_parameter_refinement(refine, instrument_parameter_string, project, path_to_project):
    if refine:
        for loop_histogram in project.histograms():
            loop_histogram.set_refinements({'Instrument Parameters': [instrument_parameter_string]})
        project.do_refinements()
        project.save(path_to_project)
        print_histogram_R_factors(project)


def export_refinement_to_csv(temp_save_directory, name_of_project, project):
    for histogram_index, loop_histogram in enumerate(project.histograms()):
        loop_histogram.Export(os.path.join(temp_save_directory, name_of_project + f"_{histogram_index}.csv"),
                              ".csv", "histogram CSV file")


def export_reflections(temp_save_directory, name_of_project, project):
    for histogram_index, loop_histogram in enumerate(project.histograms()):
        phase_names = [list(loop_histogram.reflections().keys())[0]]
        for phase_name in phase_names:
            reflection_positions = loop_histogram.reflections()[phase_name]['RefList'][:, 5]
            with open(os.path.join(temp_save_directory,
                                   name_of_project + f"_reflections_{histogram_index}_{phase_name}.txt"),
                      'wt', encoding='utf-8') as file:
                file.write(reflection_positions)


'''Parse Inputs from Mantid'''
inputs_dict = json.loads(sys.argv[1])

path_to_gsas2 = inputs_dict['path_to_gsas2']
temporary_save_directory = inputs_dict['temporary_save_directory']
data_directory = inputs_dict['data_directory']
project_name = inputs_dict['project_name']
refinement_method = inputs_dict['refinement_method']
refine_background = inputs_dict['refine_background']
refine_microstrain = inputs_dict['refine_microstrain']
refine_sigma_one = inputs_dict['refine_sigma_one']
refine_gamma = inputs_dict['refine_gamma']
refine_histogram_scale_factor = inputs_dict['refine_histogram_scale_factor']
refine_unit_cell = inputs_dict['refine_unit_cell']
override_cell_lengths = inputs_dict['override_cell_lengths']
data_files = inputs_dict['data_files']
histogram_indexing = inputs_dict['histogram_indexing']
phase_files = inputs_dict['phase_files']
instrument_files = inputs_dict['instrument_files']
limits = inputs_dict['limits']
mantid_pawley_reflections = inputs_dict['mantid_pawley_reflections']
d_spacing_min = inputs_dict['d_spacing_min']


'''Call GSASIIscriptable'''

try:
    import_path = path_to_gsas2 + 'GSASII'
    sys.path.insert(0, import_path)
    import GSASIIscriptable as G2sc  # noqa: E402
except ModuleNotFoundError:
    raise ImportError(f"GSAS-II was not found at {import_path}")


project_path = os.path.join(temporary_save_directory, project_name + '.gpx')
gsas_project = G2sc.G2Project(filename=project_path)

for phase_file in phase_files:
    gsas_project.add_phase(os.path.join(data_directory, phase_file))

number_histograms = calculate_number_histograms(len(data_files), len(histogram_indexing))
assigned_instruments = assign_instruments_to_histograms(instrument_files, number_histograms)
add_histograms(data_files, histogram_indexing, gsas_project, data_directory, assigned_instruments)

if refinement_method == "Pawley" and mantid_pawley_reflections:
    add_pawley_reflections(mantid_pawley_reflections, gsas_project, d_spacing_min)

set_max_number_cycles(3)
enable_background(refine_background, gsas_project)
enable_histogram_scale_factor(refine_histogram_scale_factor, gsas_project)

enable_unit_cell(refine_unit_cell, override_cell_lengths, gsas_project)
enable_limits(limits, gsas_project)

gsas_project.save(project_path)
gsas_project.do_refinements()
gsas_project.save(project_path)
print_histogram_R_factors(gsas_project)

run_microstrain_refinement(refine_microstrain, gsas_project, project_path)
run_parameter_refinement(refine_sigma_one, 'sig-1', gsas_project, project_path)
run_parameter_refinement(refine_gamma, 'Y', gsas_project, project_path)

export_refinement_to_csv(temporary_save_directory, project_name, gsas_project)
if refinement_method == "Pawley":
    export_reflections(temporary_save_directory, project_name, gsas_project)
