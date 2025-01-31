# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
import json
import numpy as np

import abins
from abins.input import AbInitioLoader


class Tester(object):
    """Base class for testing Abins input loaders"""

    _loaders_extensions = {"CASTEPLoader": "phonon",
                           "CRYSTALLoader": "out",
                           "DMOL3Loader": "outmol",
                           "EuphonicLoader": "castep_bin",
                           "GAUSSIANLoader": "log",
                           "VASPLoader": "xml",
                           "VASPOUTCARLoader": "OUTCAR"}

    @staticmethod
    def _prepare_data(seedname):
        """Reads reference values from ASCII files

        :param seedname: Reference data will read from the file {seedname}_data.txt, except for the atomic displacements
            which are read from files {seedname}_atomic_displacements_data_{I}.txt, where {I} are k-point indices.
        :type seedname: str
        """

        with open(abins.test_helpers.find_file(seedname + "_data.txt")) as data_file:
            correct_data = json.loads(data_file.read().replace("\n", " "))

        num_k = len(correct_data["datasets"]["k_points_data"]["weights"])
        atoms = len(correct_data["datasets"]["atoms_data"])
        array = {}
        for k in range(num_k):

            temp = np.loadtxt(abins.test_helpers.find_file(
                    "{seedname}_atomic_displacements_data_{k}.txt".format(seedname=seedname, k=k))
                ).view(complex).reshape(-1)
            total_size = temp.size
            num_freq = int(total_size / (atoms * 3))
            array[str(k)] = temp.reshape(atoms, num_freq, 3)

            freq = correct_data["datasets"]["k_points_data"]["frequencies"][str(k)]
            correct_data["datasets"]["k_points_data"]["frequencies"][str(k)] = np.asarray(freq)

        correct_data["datasets"]["k_points_data"].update({"atomic_displacements": array})

        return correct_data

    @staticmethod
    def _cull_imaginary_modes(frequencies, displacements):
        from abins.constants import ACOUSTIC_PHONON_THRESHOLD
        finite_mode_indices = frequencies > ACOUSTIC_PHONON_THRESHOLD
        return (frequencies[finite_mode_indices], displacements[:, finite_mode_indices])

    def _check_reader_data(self, correct_data=None, data=None, filename=None, extension=None):
        # check data
        correct_k_points = correct_data["datasets"]["k_points_data"]
        items = data["datasets"]["k_points_data"]

        for k in correct_k_points["frequencies"]:
            correct_frequencies, correct_displacements = self._cull_imaginary_modes(
                correct_k_points["frequencies"][k], correct_k_points["atomic_displacements"][k])
            calc_frequencies, calc_displacements = self._cull_imaginary_modes(
                items["frequencies"][k], items["atomic_displacements"][k])
            self.assertEqual(True, np.allclose(correct_frequencies, calc_frequencies))
            self.assertEqual(True, np.allclose(correct_displacements, calc_displacements))
            self.assertEqual(True, np.allclose(correct_k_points["k_vectors"][k], items["k_vectors"][k]))
            self.assertEqual(correct_k_points["weights"][k], items["weights"][k])

        correct_atoms = correct_data["datasets"]["atoms_data"]
        atoms = data["datasets"]["atoms_data"]
        for item in range(len(correct_atoms)):
            self.assertEqual(correct_atoms["atom_%s" % item]["sort"], atoms["atom_%s" % item]["sort"])
            self.assertAlmostEqual(correct_atoms["atom_%s" % item]["mass"], atoms["atom_%s" % item]["mass"],
                                   delta=0.00001)  # delta in amu units
            self.assertEqual(correct_atoms["atom_%s" % item]["symbol"], atoms["atom_%s" % item]["symbol"])
            self.assertEqual(True, np.allclose(np.array(correct_atoms["atom_%s" % item]["coord"]),
                                               atoms["atom_%s" % item]["coord"]))

        # check attributes
        self.assertEqual(correct_data["attributes"]["hash"], data["attributes"]["hash"])
        self.assertEqual(correct_data["attributes"]["ab_initio_program"], data["attributes"]["ab_initio_program"])
        try:
            self.assertEqual(abins.test_helpers.find_file(filename + "." + extension),
                             data["attributes"]["filename"])
        except AssertionError:
            self.assertEqual(abins.test_helpers.find_file(filename + "." + extension.upper()),
                             data["attributes"]["filename"])

        # check datasets
        self.assertEqual(True, np.allclose(correct_data["datasets"]["unit_cell"], data["datasets"]["unit_cell"]))

    def _check_loader_data(self, correct_data=None, input_ab_initio_filename=None, extension=None, loader=None):

        try:
            read_filename = abins.test_helpers.find_file(input_ab_initio_filename + "." + extension)
            ab_initio_loader = loader(input_ab_initio_filename=read_filename)
        except ValueError:
            read_filename = abins.test_helpers.find_file(input_ab_initio_filename + "." + extension.upper())
            ab_initio_loader = loader(input_ab_initio_filename=read_filename)

        loaded_data = ab_initio_loader.load_formatted_data().extract()

        # k points
        correct_items = correct_data["datasets"]["k_points_data"]
        items = loaded_data["k_points_data"]

        for k in correct_items["frequencies"]:
            correct_frequencies, correct_displacements = self._cull_imaginary_modes(
                correct_items["frequencies"][k], correct_items["atomic_displacements"][k])
            calc_frequencies, calc_displacements = self._cull_imaginary_modes(
                items["frequencies"][k], items["atomic_displacements"][k])

            self.assertEqual(True, np.allclose(correct_frequencies, calc_frequencies))
            self.assertEqual(True, np.allclose(correct_displacements, calc_displacements))
            self.assertEqual(True, np.allclose(correct_items["k_vectors"][k], items["k_vectors"][k]))
            self.assertEqual(correct_items["weights"][k], items["weights"][k])

        # atoms
        correct_atoms = correct_data["datasets"]["atoms_data"]
        atoms = loaded_data["atoms_data"]

        for item in range(len(correct_atoms)):
            self.assertEqual(correct_atoms["atom_%s" % item]["sort"], atoms["atom_%s" % item]["sort"])
            self.assertAlmostEqual(correct_atoms["atom_%s" % item]["mass"], atoms["atom_%s" % item]["mass"],
                                   delta=0.00001)
            self.assertEqual(correct_atoms["atom_%s" % item]["symbol"], atoms["atom_%s" % item]["symbol"])
            self.assertEqual(True, np.allclose(np.array(correct_atoms["atom_%s" % item]["coord"]),
                                               atoms["atom_%s" % item]["coord"]))

    def check(self, *,
              name: str,
              loader: AbInitioLoader,
              extension: str = None):
        """Run loader and compare output with reference files

        Args:
            name: prefix for test files (e.g. 'ethane_LoadVASP')

        """
        if extension is None:
            extension = self._loaders_extensions[str(loader)]

        # get calculated data
        data = self._read_ab_initio(loader=loader, filename=name, extension=extension)

        # get correct data
        correct_data = self._prepare_data(name)

        # check read data
        self._check_reader_data(correct_data=correct_data, data=data, filename=name, extension=extension)

        # check loaded data
        self._check_loader_data(correct_data=correct_data,
                                input_ab_initio_filename=name,
                                extension=extension,
                                loader=loader)

    def _read_ab_initio(self, loader=None, filename=None, extension=None):
        """
        Reads data from .{extension} file.
        :param loader: ab initio loader
        :param filename: name of file with vibrational or phonon data (name + extension)
        :returns: vibrational or phonon data
        """
        # 1) Read data
        try:
            read_filename = abins.test_helpers.find_file(filename=filename + "." + extension)
            ab_initio_reader = loader(input_ab_initio_filename=read_filename)
        except ValueError:
            read_filename = abins.test_helpers.find_file(filename=filename + "." + extension.upper())
            ab_initio_reader = loader(input_ab_initio_filename=read_filename)

        data = self._get_reader_data(ab_initio_reader)

        # test validData method
        self.assertEqual(True, ab_initio_reader._clerk._valid_hash())

        return data

    @staticmethod
    def _get_reader_data(ab_initio_reader):
        """
        :param ab_initio_reader: object of type AbInitioLoader
        :returns: read data
        """
        abins_type_data = ab_initio_reader.read_vibrational_or_phonon_data()
        data = {"datasets": abins_type_data.extract(),
                "attributes": ab_initio_reader._clerk._attributes.copy()
                }
        data["datasets"].update({"unit_cell": ab_initio_reader._clerk._data["unit_cell"]})
        return data

    @classmethod
    def save_ab_initio_test_data(cls, ab_initio_reader, seedname):
        """
        Write ab initio calculation data to JSON file for use in test cases

        :param ab_initio_reader: Reader after import of external calculation
        :type ab_initio_reader: abins.AbInitioProgram
        :param filename: Seed for text files for JSON output. Data will be written to the file {seedname}_data.txt,
            except for the atomic displacements which are written to files {seedname}_atomic_displacements_data_{I}.txt,
            where {I} are k-point indices.
        :type filename: str

        """

        data = cls._get_reader_data(ab_initio_reader)

        # Discard advanced_parameters cache as this is not relevant to loader tests
        del data["attributes"]["advanced_parameters"]

        displacements = data["datasets"]["k_points_data"].pop("atomic_displacements")
        for i, eigenvector in displacements.items():
            with open('{seedname}_atomic_displacements_data_{i}.txt'.format(seedname=seedname, i=i), 'wt') as f:
                eigenvector.flatten().view(float).tofile(f, sep=' ')

        with open('{seedname}_data.txt'.format(seedname=seedname), 'wt') as f:
            json.dump(cls._arrays_to_lists(data), f, indent=4, sort_keys=True)

    @classmethod
    def _arrays_to_lists(cls, mydict):
        """Recursively convert numpy arrays in a nested dict to lists (i.e. valid JSON)

        Returns a processed *copy* of the input dictionary: in-place values will not be altered."""
        clean_dict = {}
        for key, value in mydict.items():
            if isinstance(value, np.ndarray):
                clean_dict[key] = value.tolist()
            elif isinstance(value, dict):
                clean_dict[key] = cls._arrays_to_lists(value)
            else:
                clean_dict[key] = value
        return clean_dict
