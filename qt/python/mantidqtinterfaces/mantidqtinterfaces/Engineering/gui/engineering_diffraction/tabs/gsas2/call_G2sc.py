# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
import os
import sys
import json
import numpy as np


'''Parse Inputs from Mantid'''
inputs_json = sys.argv[1]
inputs_dict = json.loads(inputs_json)

path_to_gsas2 = inputs_dict['path_to_gsas2']
temporary_save_directory = inputs_dict['temporary_save_directory']
data_directory = inputs_dict['data_directory']
refinement_method = inputs_dict['refinement_method']
project_name = inputs_dict['project_name']

refine_background = inputs_dict['refine_background']
refine_microstrain = inputs_dict['refine_microstrain']
refine_sigma_one = inputs_dict['refine_sigma_one']
refine_gamma = inputs_dict['refine_gamma']
refine_histogram_scale_factor = inputs_dict['refine_histogram_scale_factor']

data_files = inputs_dict['data_files']
histogram_indexing = inputs_dict['histogram_indexing']

# Check the number of input histograms
number_histograms = len(data_files)
if histogram_indexing and len(data_files) == 1:
    number_histograms = len(histogram_indexing)
# if histogram_indexing and len(data_files) > 1 should be caught in Validation

phase_files = inputs_dict['phase_files']
instrument_files = inputs_dict['instrument_files']
x_min, x_max = inputs_dict['limits']
mantid_pawley_reflections = inputs_dict['mantid_pawley_reflections']
override_cell_lengths = inputs_dict['override_cell_lengths']


'''Call GSASIIscriptable'''
project_path = os.path.join(temporary_save_directory, project_name + '.gpx')

sys.path.insert(0, path_to_gsas2 + 'GSASII')
import GSASIIscriptable as G2sc  # noqa: E402
# Maybe add a try catch statement?


def HistStats(gpx):
    '''prints profile rfactors for all histograms'''
    print(u"*** profile Rwp, " + os.path.split(gpx.filename)[1])
    for hist in gpx.histograms():
        print("\t{:20s}: {:.2f}".format(hist.name, hist.get_wR()))
    print("")
    gpx.save()


gpx = G2sc.G2Project(filename=project_path)

gsas_phases = []
for phase_file in phase_files:
    gsas_phases.append(gpx.add_phase(os.path.join(data_directory, phase_file)))

# Assign instruments to histograms
if len(instrument_files) == 1:
    iparams_input = [instrument_files[0]] * number_histograms
elif len(instrument_files) > 1 and len(instrument_files) == number_histograms:
    iparams_input = instrument_files
else:
    raise ValueError(f'The number of instrument files ({len(instrument_files)}) must be 1 '
                     f'or equal to the number of input histograms {number_histograms}')

# Add histograms with instruments and phases
gsas_histograms = []
if not histogram_indexing:
    for data_file_index, input_data_file in enumerate(data_files):
        gsas_histograms.append(gpx.add_powder_histogram(datafile=os.path.join(data_directory,
                                                                              input_data_file),
                                                        iparams=os.path.join(data_directory,
                                                                             iparams_input[data_file_index]),
                                                        phases=gsas_phases
                                                        ))
else:
    for index_in_list, histogram_index in enumerate(histogram_indexing):
        gsas_histograms.append(gpx.add_powder_histogram(datafile=os.path.join(data_directory,
                                                                              data_files[0]),
                                                        iparams=os.path.join(data_directory,
                                                                             iparams_input[index_in_list]),
                                                        phases=gsas_phases,
                                                        databank=histogram_index,  # indexing starts at 1
                                                        ))

dmin = 1.0
peaks_to_add = set()
if refinement_method == "Pawley" and mantid_pawley_reflections:
    for gsas_phase in gsas_phases:
        gsas_phase.data['General']['doPawley'] = True

        gsas_reflections = []
        for reflection in mantid_pawley_reflections:
            [h, k, l], d, multiplicity = reflection
            gsas_reflections.append([int(h), int(k), int(l), int(multiplicity), float(d), True, 100.0, 1.0])
        gsas_phase.data["Pawley ref"] = gsas_reflections


# for i in G2sc.dictDive(phase.data['General'], 'paw'): print(i)

# increase # of cycles to improve convergence
gpx.data['Controls']['data']['max cyc'] = 3  # not in API

# tutorial step 4: turn on background refinement (Hist)
for index, histogram in enumerate(gsas_histograms):
    if refine_background:
        histogram.set_refinements({"Background": {"no. coeffs": 3, "refine": True}})
    else:
        histogram.set_refinements({"Background": {"no. coeffs": 0, "refine": False}})

if not refine_histogram_scale_factor:
    for gsas_histogram in gpx.histograms():
        gsas_histogram.SampleParameters["Scale"] = [1.0, False]

for gsas_phase in gpx.phases():
    gsas_phase.set_refinements({"Cell": True})
    if override_cell_lengths:
        gsas_phase.data['General']['Cell'][1:4] = tuple(override_cell_lengths)

if x_min and x_max:
    for index, histogram in enumerate(gsas_histograms):
        histogram.set_refinements({'Limits': [x_min[index], x_max[index]]})

gpx.save(project_path)
gpx.do_refinements()
gpx.save(project_path)
HistStats(gpx)

if refine_microstrain:
    for gsas_phase in gpx.phases():
        gsas_phase.set_HAP_refinements({'Mustrain': { 'type': 'isotropic', 'refine': True}})
    print("Refining Microstrain")
    gpx.do_refinements()
    gpx.save(project_path)
    HistStats(gpx)

if refine_sigma_one:
    for gsas_histogram in gpx.histograms():
        gsas_histogram.set_refinements({'Instrument Parameters': ['sig-1']})
    print("Refining Sigma-1")
    gpx.do_refinements()
    gpx.save(project_path)
    HistStats(gpx)

if refine_gamma:
    for gsas_histogram in gpx.histograms():
        gsas_histogram.set_refinements({'Instrument Parameters': ['Y']})
    print("Refining Gamma")
    gpx.do_refinements()
    gpx.save(project_path)
    HistStats(gpx)

for index, histogram in enumerate(gpx.histograms()):
    histogram.Export(os.path.join(temporary_save_directory, project_name + f"_{index}.csv"), ".csv", "histogram CSV file")

    # Assuming only one phase
if refinement_method == 'Pawley':
    phase_name = list(histogram.reflections().keys())[0]
    reflection_positions = np.transpose(np.array(histogram.reflections()[phase_name]['RefList']))[5]
    np.savetxt(os.path.join(temporary_save_directory, project_name + f"_reflections_{index}.txt"), reflection_positions)
