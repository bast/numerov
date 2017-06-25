import yaml
import sys
import numpy
from .solve import solve_numerov
from .constants import constants
from .__init__ import __version__
import math


def main():
    if len(sys.argv) == 1:
        # user has given no arguments: print help and exit
        sys.stderr.write("usage: python numerov.py <input.yml>\n")
        sys.exit(-1)

    input_file = sys.argv[-1]
    with open(input_file, 'r') as stream:
        try:
            input_data = yaml.load(stream)
        except yaml.YAMLError as e:
            sys.stderr.write(e)
            sys.stderr.write("\n")
            sys.exit(-1)

    displacements = []
    pot_energies = []
    exp_values = []
    for step in input_data['steps']:
        displacements.append(step['displacement'])
        pot_energies.append(step['pot_energy'])
        exp_values.append(step['exp_value'])

    displacements = numpy.array(displacements)
    pot_energies = numpy.array(pot_energies)
    exp_values = numpy.array(exp_values)

    # shift potential such that minimum is at zero
    pot_energies -= min(pot_energies)

    pot_energy_coefs = numpy.polyfit(displacements, pot_energies, input_data['degree_pot_energy'])
    exp_value_coefs = numpy.polyfit(displacements, exp_values, input_data['degree_exp_value'])

    transition_frequency_previous = sys.float_info.max
    displacement_range = (-0.5, 0.5)
    while True:
        qs, psi_squared, energies_hartree, averaged_exp_values_au = solve_numerov(pot_energy_coefs,
                                                                                  exp_value_coefs,
                                                                                  displacement_range,
                                                                                  input_data['num_steps'],
                                                                                  input_data['num_solutions'],
                                                                                  input_data['energy_precision_hartree'],
                                                                                  input_data['reduced_mass_amu'] * constants['amu_to_au'])
        transition_frequency = (energies_hartree[input_data['num_solutions'] - 1] - energies_hartree[0]) * constants['hartree_to_cm1']
        if abs(transition_frequency - transition_frequency_previous) < 1.0e-1:
            qs = numpy.around(qs, decimals=3).tolist()
            psi_squared = numpy.around(psi_squared, decimals=5).tolist()
            break
        displacement_range = (displacement_range[0] - 0.1, displacement_range[1] + 0.1)
        transition_frequency_previous = transition_frequency

    # get harmonic frequency from numerov
    # this can be used as a sanity check to verify whether this matches
    # the provided input frequency
    fourth_lowest_energy = sorted(pot_energies)[3]
    h_x = []
    h_y = []
    for (displacement, pot_energy) in zip(displacements, pot_energies):
        if pot_energy < fourth_lowest_energy:
            h_x.append(displacement)
            h_y.append(pot_energy)
    pot_energy_coefs_harmonic = numpy.polyfit(h_x, h_y, 2)
    _, _, energies_hartree_harmonic, _ = solve_numerov(pot_energy_coefs_harmonic,
                                                       exp_value_coefs,
                                                       displacement_range,
                                                       input_data['num_steps'],
                                                       input_data['num_solutions'],
                                                       input_data['energy_precision_hartree'],
                                                       input_data['reduced_mass_amu'] * constants['amu_to_au'])
    for i in range(len(exp_value_coefs)):
        exp_value_coefs[-i - 1] *= math.factorial(i)
    for i in range(len(pot_energy_coefs)):
        pot_energy_coefs[-i - 1] *= math.factorial(i)

    print(yaml.dump(input_data))
    print('version: {0}'.format(__version__))
    print('displacement_range:')
    print('  min: {0}'.format(displacement_range[0]))
    print('  max: {0}'.format(displacement_range[1]))
    print('exp_value_coefs: {0}'.format(list(exp_value_coefs)))
    print('pot_energy_coefs: {0}'.format(list(pot_energy_coefs)))
    print('pot_energy_coefs_harmonic: {0}'.format(list(pot_energy_coefs_harmonic)))
    print('energies_hartree: {0}'.format(list(energies_hartree)))
    print('energies_hartree_harmonic: {0}'.format(list(energies_hartree_harmonic)))
    print('averaged_exp_values_au: {0}'.format(list(averaged_exp_values_au)))
    print('qs: {0}'.format(qs))
    print('psi_squared: {0}'.format(psi_squared))
