# Copyright 2017 Radovan Bast. All rights reserved.
# Use of this source code is governed by a the Mozilla Public License v2.0 that
# can be found in the LICENSE file. If a copy of the license was not
# distributed with this file, you can obtain one at
# http://mozilla.org/MPL/2.0/.


import yaml
import sys
import numpy
from _numerov import solve_numerov, constants
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

    fourth_lowest_energy = sorted(pot_energies)[3]
    h_x = []
    h_y = []
    for (displacement, pot_energy) in zip(displacements, pot_energies):
        if pot_energy < fourth_lowest_energy:
            h_x.append(displacement)
            h_y.append(pot_energy)
    pot_energy_coefs_harmonic = numpy.polyfit(h_x, h_y, 2)

    transition_frequency_previous = sys.float_info.max
    displacement_range = (-0.5, 0.5)
    while True:
        qs, psi_squared, energies_hartree, averaged_exp_values_au = solve_numerov(pot_energy_coefs,
                                                                                  exp_value_coefs,
                                                                                  displacement_range,
                                                                                  input_data['num_steps'],
                                                                                  input_data['num_solutions'],
                                                                                  input_data['energy_precision'],
                                                                                  input_data['reduced_mass_amu'] * constants['amu_to_au'])
        transition_frequency = (energies_hartree[input_data['num_solutions'] - 1] - energies_hartree[0]) * constants['hartree_to_cm1']
        if abs(transition_frequency - transition_frequency_previous) < 1.0e-1:
            break
        displacement_range = (displacement_range[0] - 0.1, displacement_range[1] + 0.1)
        transition_frequency_previous = transition_frequency

    print(yaml.dump(input_data))
    print('displacement_range:')
    print('  min: {0}'.format(displacement_range[0]))
    print('  max: {0}'.format(displacement_range[1]))

    # get harmonic frequency from numerov
    _, _, energies_hartree_harmonic, _ = solve_numerov(pot_energy_coefs_harmonic,
                                                       exp_value_coefs,
                                                       displacement_range,
                                                       input_data['num_steps'],
                                                       input_data['num_solutions'],
                                                       input_data['energy_precision'],
                                                       input_data['reduced_mass_amu'] * constants['amu_to_au'])
    for i in range(len(exp_value_coefs)):
        exp_value_coefs[-i - 1] *= math.factorial(i)
    for i in range(len(pot_energy_coefs)):
        pot_energy_coefs[-i - 1] *= math.factorial(i)

    print('exp_value_coefs: {0}'.format(list(exp_value_coefs)))
    print('pot_energy_coefs: {0}'.format(list(pot_energy_coefs)))
    print('energies_hartree: {0}'.format(list(energies_hartree)))
    print('energies_hartree_harmonic: {0}'.format(list(energies_hartree_harmonic)))
    print('averaged_exp_values_au: {0}'.format(list(averaged_exp_values_au)))


if __name__ == '__main__':
    main()
