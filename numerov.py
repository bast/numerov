# Copyright 2017 Radovan Bast. All rights reserved.
# Use of this source code is governed by a the Mozilla Public License v2.0 that
# can be found in the LICENSE file. If a copy of the license was not
# distributed with this file, you can obtain one at
# http://mozilla.org/MPL/2.0/.


import yaml
import sys
import numpy


def run_numerov(pot_energy_coefs,
                exp_value_coefs,
                displacement_range,
                num_steps,
                num_solutions,
                energy_precision,
                reduced_mass_au):

    n = num_steps + 1
    step = (displacement_range[1] - displacement_range[0]) / num_steps
    step2 = step**2.0

    q = numpy.zeros(n)
    pot_energies = numpy.zeros(n)
    exp_values = numpy.zeros(n)
    g = numpy.zeros(n)
    psi = numpy.zeros(n)
    psi_squared = numpy.zeros([num_solutions, n])
    energies_hartree = numpy.zeros(num_solutions)
    averaged_exp_values_au = numpy.zeros(num_solutions)

    for i in range(n):
        q[i] = displacement_range[0] + i * step
        pot_energies[i] = numpy.polyval(pot_energy_coefs, q[i])
        exp_values[i] = numpy.polyval(exp_value_coefs, q[i])

    energy_guess = 1e-4
    energy_step = 1e-4
    num_nodes_last = 0

    while num_nodes_last < num_solutions:

        g = reduced_mass_au * 2.0 * (pot_energies - energy_guess)

        psi[0] = 0.0
        psi[1] = 1.0e-6

        for i in range(n):
            if i > 1:
                t0 = g[i - 1] * psi[i - 1] * step2 * 5.0 / 6.0
                t1 = g[i - 2] * psi[i - 2] * step2 / 12.0
                t2 = 2.0 * psi[i - 1] - psi[i - 2] + t0 + t1
                t3 = 1.0 - g[i] * step2 / 12.0
                psi[i] = t2 / t3

        num_nodes = 0
        i_save = n
        for i in range(n):
            if i > 1:
                if psi[i - 1] != 0.0:
                    if (psi[i] / psi[i - 1]) < 0.0:
                        num_nodes = num_nodes + 1
                        i_save = i

        psi[i_save + 1:] = 0.0

        psi = psi / numpy.sqrt(numpy.dot(psi, psi))

        if (abs(energy_step) < energy_precision) and (num_nodes > num_nodes_last):

            norm = 0.0
            done = 0
            for i in range(len(psi)):
                norm = norm + psi[i] * psi[i]
                if (norm > 0.001) and not done:
                    i_left = i
                    done = 1

            norm = 0.0
            done = 0
            for i in range(n - 1, -1, -1):
                norm = norm + psi[i] * psi[i]
                if (norm > 0.001) and not done:
                    i_right = i
                    done = 1

            psi_squared[num_nodes - 1] = psi**2.0
            energies_hartree[num_nodes - 1] = energy_guess
            averaged_exp_values_au[num_nodes - 1] = numpy.dot(psi, exp_values * psi)

            energy_step = 1.0e-4
            num_nodes_last = num_nodes_last + 1

        if num_nodes > num_nodes_last:
            if energy_step > 0.0:
                energy_step = energy_step / (-10.0)
        else:
            if energy_step < 0.0:
                energy_step = energy_step / (-10.0)

        energy_guess = energy_guess + energy_step

    return q, psi_squared, energies_hartree, averaged_exp_values_au


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
        pot_energies.append(step['potential'])
        exp_values.append(step['property'])

    displacements = numpy.array(displacements)
    pot_energies = numpy.array(pot_energies)
    exp_values = numpy.array(exp_values)

    # shift potential such that minimum is at zero
    pot_energies -= min(pot_energies)

    pot_energy_coefs = numpy.polyfit(displacements, pot_energies, input_data['degree_potential'])
    exp_value_coefs = numpy.polyfit(displacements, exp_values, input_data['degree_property'])

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
    print('degree potential   =', input_data['degree_potential'])
    print('degree property    =', input_data['degree_property'])
    print('energy precision   =', input_data['energy_precision'])
    print('nr solutions       =', input_data['num_solutions'])
    print('nr steps           =', input_data['num_steps'])

    transition_frequency_previous = sys.float_info.max
    displacement_range = (-0.5, 0.5)
    while True:
        q, psi_squared, energies_hartree, averaged_exp_values_au = run_numerov(pot_energy_coefs,
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
    _, _, energies_hartree, _ = run_numerov(pot_energy_coefs_harmonic,
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


def test_run_numerov():
    pot_energy_coefs = numpy.array([-2.45560869e-01, -8.88252151e-03,
                                    1.24439946e-01, 1.93259856e-01, 2.78860663e-01, -5.62738650e-05,
                                    -5.78784571e-08])
    exp_value_coefs = numpy.array([-1.45680171e-14, -2.78094078e-16,
                                   1.62725432e-15, -1.99732822e-15, 3.08772558e-15, 2.06211298e-15,
                                   -7.30656049e-15])
    displacement_range = (-0.5, 0.5)
    num_steps = 21
    num_solutions = 3
    energy_precision = 1e-12
    reduced_mass_au = 26245.03
    q, psi_squared, energies_hartree, averaged_exp_values_au = run_numerov(pot_energy_coefs,
                                                                           exp_value_coefs,
                                                                           displacement_range,
                                                                           num_steps,
                                                                           num_solutions,
                                                                           energy_precision,
                                                                           reduced_mass_au)

    q_ref = [-0.5, -0.45238095238095238, -0.40476190476190477,
             -0.35714285714285715, -0.30952380952380953, -0.26190476190476192,
             -0.2142857142857143, -0.16666666666666669, -0.11904761904761907,
             -0.071428571428571452, -0.023809523809523836, 0.023809523809523725,
             0.071428571428571397, 0.11904761904761907, 0.16666666666666663,
             0.21428571428571419, 0.26190476190476186, 0.30952380952380953,
             0.3571428571428571, 0.40476190476190466, 0.45238095238095233, 0.5]
    assert numpy.allclose(q, q_ref)

    psi_0_squared_ref = [0.0, 5.8115604540099223e-11, 4.5928270821175377e-09,
                         2.2581144389236216e-07, 7.0160833333126289e-06, 0.00013741125594800236,
                         0.0016860686764846274, 0.012837351969918346, 0.059868399265284117,
                         0.16833949851844193, 0.28036406815649312, 0.27146736153858203,
                         0.15004107659776397, 0.046521109971706026, 0.0079563591949426953,
                         0.00073700740789436057, 3.6124036211418563e-05, 9.0535270296150204e-07,
                         1.0972855205998383e-08, 4.808188792726663e-11, 4.9096949798023428e-10, 0.0]

    assert numpy.allclose(psi_squared[0], psi_0_squared_ref)

    energies_hartree_ref = [0.0023035871039000012, 0.0068957376078000003, 0.011460087763699986]
    assert numpy.allclose(energies_hartree, energies_hartree_ref)

    averaged_exp_values_au_ref = [-7.3021260330311477e-15, -7.2926761553595848e-15, -7.2818315936417763e-15]
    assert numpy.allclose(averaged_exp_values_au, averaged_exp_values_au_ref, rtol=1e-05, atol=1e-30)

if __name__ == '__main__':
    main()
