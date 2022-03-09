import os
import sys
count = 0


def counter():
    global count
    count += 1
    return count


path_to_gsas2 = sys.argv[counter()]
save_directory = sys.argv[counter()]
data_directory = sys.argv[counter()]
refinement_method = sys.argv[counter()]
project_name = sys.argv[counter()]

# number of each dynamic input
number_data_files = int(sys.argv[counter()])
number_histogram_indices = int(sys.argv[counter()])
number_phases = int(sys.argv[counter()])
number_instruments = int(sys.argv[counter()])
number_limits = int(sys.argv[counter()])

data_files = []
for i in range(number_data_files):
    data_files.append(sys.argv[counter()])

histogram_indexing = []
for i in range(number_histogram_indices):
    histogram_indexing.append(int(sys.argv[counter()]))

# Check the number of input histograms
number_histograms = len(data_files)
if histogram_indexing and len(data_files) == 1:
    number_histograms = len(histogram_indexing)
# if histogram_indexing and len(data_files) > 1 should be caught in Validation

# for now, all phases are applied to all histograms
phases = []
for i in range(number_phases):
    phases.append(sys.argv[counter()])

instruments = []
for i in range(number_instruments):
    instruments.append(sys.argv[counter()])

x_min = []
x_max = []
if number_limits != 0:
    for i in range(number_limits):
        x_min.append(float(sys.argv[counter()]))
    for i in range(number_limits):
        x_max.append(float(sys.argv[counter()]))


project_path = save_directory + project_name + '.gpx'

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
for phase_file in phases:
    gsas_phases.append(gpx.add_phase(os.path.join(data_directory, phase_file)))

# Assign instruments to histograms
if number_instruments == 1:
    iparams_input = [instruments[0]] * number_histograms
elif number_instruments > 1 and number_instruments == number_histograms:
    iparams_input = instruments
else:
    raise ValueError(f'The number of instrument files ({number_instruments}) must be 1 '
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

if refinement_method == "Pawley":
    for gsas_phase in gsas_phases:
        gsas_phase.data['General']['doPawley'] = True
# for i in G2sc.dictDive(phase.data['General'], 'paw'): print(i)

# increase # of cycles to improve convergence
gpx.data['Controls']['data']['max cyc'] = 8  # not in API

# tutorial step 4: turn on background refinement (Hist)
refdict0 = {"set": {"Background": {"no. coeffs": 3, "refine": True}}}

if x_min and x_max:
    for index, histogram in enumerate(gsas_histograms):
        histogram.set_refinements({'Limits': [x_min[index], x_max[index]]})

gpx.save(project_path)
gpx.do_refinements([refdict0])
gpx.save(project_path)
HistStats(gpx)

for index, histogram in enumerate(gpx.histograms()):
    histogram.Export(os.path.join(save_directory, project_name + f"_{index}.csv"), ".csv", "histogram CSV file")
