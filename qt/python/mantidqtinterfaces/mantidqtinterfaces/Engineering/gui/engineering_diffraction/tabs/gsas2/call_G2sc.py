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
number_histograms = int(sys.argv[counter()])
number_phases = int(sys.argv[counter()])
number_instruments = int(sys.argv[counter()])
number_limits = int(sys.argv[counter()])


histograms = []
for i in range(number_histograms):
    histograms.append(sys.argv[counter()])

# for now I use only the first phase file, even if there are multiple 040322
phases = []
for i in range(number_phases):
    phases.append(sys.argv[counter()])

instruments = []
for i in range(number_instruments):
    instruments.append(sys.argv[counter()])

print('LOOK here:', instruments)

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

gsas_histograms = []
for histogram_index, input_data_file in enumerate(histograms):
    gsas_histograms.append(gpx.add_powder_histogram(datafile=os.path.join(data_directory, input_data_file),
                                                    iparams=os.path.join(data_directory, iparams_input[histogram_index]),
                                                    databank=1,  # indexing starts at 1
                                                    phases=[gsas_phases[0]]
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
HistStats(gpx)
