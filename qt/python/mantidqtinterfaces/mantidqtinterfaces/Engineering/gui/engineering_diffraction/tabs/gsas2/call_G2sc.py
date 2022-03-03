import os
import sys

path_to_gsas2 = sys.argv[1]
save_directory = sys.argv[2]
data_directory = sys.argv[3]
refinement_method = sys.argv[4]
input_data_file = sys.argv[5]
# input_data_file_2 = sys.argv[6]
instrument_file = sys.argv[6]
# instrument_file_2 = sys.argv[8]
phase_file = sys.argv[7]
# phase_files.append(sys.argv[11])

project_name = sys.argv[8]

x_min = [float(sys.argv[9])]
x_max = [float(sys.argv[10])]

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


# create a project with a default project name
gpx = G2sc.G2Project(filename=project_path)

print(os.path.join(data_directory, phase_file), "\n")
phase = gpx.add_phase(os.path.join(data_directory, phase_file))

histogram_1 = gpx.add_powder_histogram(datafile=os.path.join(data_directory, input_data_file),
                                       iparams=os.path.join(data_directory, instrument_file),
                                       databank=1,  # indexing starts at 1
                                       phases=[phase]
                                       )

phase.data['General']['doPawley'] = True
# for i in G2sc.dictDive(phase.data['General'], 'paw'): print(i)

# increase # of cycles to improve convergence
gpx.data['Controls']['data']['max cyc'] = 8  # not in API

# tutorial step 4: turn on background refinement (Hist)
refdict0 = {"set": {"Background": {"no. coeffs": 3, "refine": True}}}

# if x_min and x_max:
#     histogram_1.set_refinements({'Limits': [x_min[0], x_max[0]]})
#     # histogram_2.set_refinements({'Limits': [x_min[1], x_max[1]]})

gpx.save(project_path)
gpx.do_refinements([refdict0])
HistStats(gpx)
