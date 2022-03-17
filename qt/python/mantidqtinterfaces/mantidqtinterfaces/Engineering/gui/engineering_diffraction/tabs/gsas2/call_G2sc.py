import os
import sys
count = 0


def counter():
    global count
    count += 1
    return count


'''Parse Inputs from Mantid'''
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
number_reflections = int(sys.argv[counter()])

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

compressed_reflections = []
if number_reflections != 0:
    for i in range(number_reflections):
        compressed_reflections.append(sys.argv[counter()])


'''Call GSASIIscriptable'''
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
        gsas_histograms[data_file_index].SampleParameters["Scale"] = [1.0, False]
else:
    for index_in_list, histogram_index in enumerate(histogram_indexing):
        gsas_histograms.append(gpx.add_powder_histogram(datafile=os.path.join(data_directory,
                                                                              data_files[0]),
                                                        iparams=os.path.join(data_directory,
                                                                             iparams_input[index_in_list]),
                                                        phases=gsas_phases,
                                                        databank=histogram_index,  # indexing starts at 1
                                                        ))
        gsas_histograms[index_in_list].SampleParameters["Scale"] = [1.0, False]

dmin = 1.0
peaks_to_add = set()
if refinement_method == "Pawley":
    for gsas_phase in gsas_phases:
        gsas_phase.data['General']['doPawley'] = True

        # cell = gsas_phase.get_cell()
        # lattice_params = [cell['length_a'], cell['length_b'], cell['length_c'],
        #                   cell['angle_alpha'], cell['angle_beta'], cell['angle_gamma']]
        # reflections = G2sc.GenerateReflections(gsas_phase.data['General']['SGData']['SpGrp'],
        #                                        lattice_params, dmin=dmin)
        # for gsas_histogram in gsas_histograms:
        #     out = gsas_histogram.getHKLpeak(dmin, gsas_phase.data['General']['SGData'], lattice_params)
        #     print(out)

        pawley_reflections = []
        # 'HKL', 'd', 'F^2', 'M'
        for compressed_reflection in compressed_reflections:
            pawley_reflections.append(compressed_reflection.split("#"))

        gsas_reflections = []
        for reflection in pawley_reflections:
            h, k, l, = reflection[0][1:-1].split(",")
            d, F_sq, multiplicity = reflection[1:]
            gsas_reflections.append([int(h), int(k), int(l), int(multiplicity), float(d), True, 100.0, 1.0])

        gsas_phase.data["Pawley ref"] = gsas_reflections


# for i in G2sc.dictDive(phase.data['General'], 'paw'): print(i)

# increase # of cycles to improve convergence
gpx.data['Controls']['data']['max cyc'] = 8  # not in API

# tutorial step 4: turn on background refinement (Hist)
refdict0 = {"set": {"Background": {"no. coeffs": 3, "refine": True}}}
# refdict0.update({ 'Mustrain': { 'type': 'isotrpoic', 'refine': True}})

for p in gpx.phases():
    p.set_refinements({"Cell": True})
    gsas_phase.data['General']['Cell'][1:4] = 3.65, 3.65, 3.65
    print(gsas_phase.data['General']['Cell'])

if x_min and x_max:
    for index, histogram in enumerate(gsas_histograms):
        histogram.set_refinements({'Limits': [x_min[index], x_max[index]]})

for index, histogram in enumerate(gpx.histograms()):
    '''Manually add peaks'''
    # peak1 = histogram.add_peak(1, ttheta=38819.06646)  # this is cheeky, I'm using the TOF value in the ttheta input
    # peak2 = histogram.add_peak(1, ttheta=33619.962029999995)
    # peak3 = histogram.add_peak(1, ttheta=23767.413545)
    # peak4 = histogram.add_peak(1, ttheta=20277.14198)
    # peak5 = histogram.add_peak(1, ttheta=19409.27622)

    '''Add Generated reflections as peaks: plots not look reasonable'''
    # for peak_dspace_value in peaks_to_add:
    #     histogram.add_peak(1, dspace=peak_dspace_value)

    # histogram.set_peakFlags(area=True)
    # histogram.refine_peaks()
    # histogram.set_peakFlags(area=True, pos=True)
    # histogram.refine_peaks()
    # histogram.set_peakFlags(area=True, pos=True, sig=True, gam=True)
    # histogram.refine_peaks()

gpx.save(project_path)
gpx.do_refinements([refdict0])
gpx.save(project_path)

HistStats(gpx)

for index, histogram in enumerate(gpx.histograms()):
    histogram.Export(os.path.join(save_directory, project_name + f"_{index}.csv"), ".csv", "histogram CSV file")
