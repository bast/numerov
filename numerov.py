# Copyright 2017 Radovan Bast. All rights reserved.
# Use of this source code is governed by a the Mozilla Public License v2.0 that
# can be found in the LICENSE file. If a copy of the license was not
# distributed with this file, you can obtain one at
# http://mozilla.org/MPL/2.0/.


import yaml
import sys
import numpy
from _numerov import solve_numerov


def main():
    constants = {'amu_to_au': 1822.888479031408,
                 'hartree_to_hz': 6.579683920721e15,
                 'hartree_to_cm1': 2.194746313705e5}

    if len(sys.argv) == 1:
        # user has given no arguments: print help and exit
        sys.stderr.write("usage: python numerov.py <input.yml>\n")
        sys.exit(-1)

    input_file = sys.argv[-1]
    with open(input_file, 'r') as stream:
        try:
            input_data = yaml.load(stream)
        except yaml.YAMLError as e:
            print(e)

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

    print()
    print('reduced mass       =', input_data['reduced_mass_amu'])
    print('harmonic frequency =', input_data['harmonic_frequency_cm1'])
    print('degree pot_energy   =', input_data['degree_pot_energy'])
    print('degree exp_value    =', input_data['degree_exp_value'])
    print('energy precision   =', input_data['energy_precision'])
    print('nr solutions       =', input_data['num_solutions'])
    print('nr steps           =', input_data['num_steps'])

    transition_frequency_previous = sys.float_info.max
    displacement_range = (-0.5, 0.5)
    while True:
        q, psi_squared, energies_hartree, averaged_exp_values_au = solve_numerov(pot_energy_coefs,
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

    # factor two is hardcoded for PNC difference
    diff_au = 2.0 * (averaged_exp_values_au[input_data['num_solutions'] - 1] - averaged_exp_values_au[0])
    diff_hz = diff_au * constants['hartree_to_hz']

    # get harmonic frequency from numerov
    _, _, energies_hartree, _ = solve_numerov(pot_energy_coefs_harmonic,
                                              exp_value_coefs,
                                              displacement_range,
                                              input_data['num_steps'],
                                              input_data['num_solutions'],
                                              input_data['energy_precision'],
                                              input_data['reduced_mass_amu'] * constants['amu_to_au'])
    transition_frequency_harmonic = (energies_hartree[input_data['num_solutions'] - 1] - energies_hartree[0]) * constants['hartree_to_cm1']

    p0 = exp_value_coefs[-1] * constants['hartree_to_hz']
    p1 = exp_value_coefs[-2] * constants['hartree_to_hz']
    p2 = exp_value_coefs[-3] * constants['hartree_to_hz'] * 2.0
    v3 = 6.0 * pot_energy_coefs[-4]
    mass = input_data['reduced_mass_amu'] * constants['amu_to_au']
    frequency = input_data['harmonic_frequency_cm1'] / constants['hartree_to_cm1']

    delta_2 = p2
    delta_3 = delta_2 - p1 * v3 / (mass * frequency * frequency)

    n = input_data['to_level']
    delta_2 *= n / (mass * frequency)
    delta_3 *= n / (mass * frequency)

    print()
    print('%12s %10s %10s %10s %10s %10s %10s %7s %7s %4s' % ('p(0)', 'p(1)', 'p(2)', 'v(3)',
                                                              'delta(2)', 'delta(3)', 'numerov', 'freq', 'freq h', 'q'))
    print('%12.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %7.1f %7.1f %4.1f' % (p0, p1, p2, v3,
                                                                                  delta_2, delta_3,
                                                                                  diff_hz,
                                                                                  transition_frequency,
                                                                                  transition_frequency_harmonic,
                                                                                  displacement_range[1]))

    plot_file = input_file.replace('.in', '')


if __name__ == '__main__':
    main()
