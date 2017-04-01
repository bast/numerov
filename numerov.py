import yaml
import sys
import numpy


def run_numerov(input_data,
                constants,
                pot_energy_coefs,
                exp_value_coefs,
                displacement_range):

    n = input_data['num_steps'] + 1
    step = (displacement_range[1] - displacement_range[0]) / input_data['num_steps']
    step2 = step**2.0

    q = numpy.zeros(n)
    pot_energies = numpy.zeros(n)
    exp_values = numpy.zeros(n)
    g = numpy.zeros(n)
    psi = numpy.zeros(n)
    psi_squared = numpy.zeros([input_data['num_solutions'], n])
    energy = numpy.zeros(input_data['num_solutions'])

    for i in range(n):
        q[i] = displacement_range[0] + i * step
        pot_energies[i] = numpy.polyval(pot_energy_coefs, q[i])
        exp_values[i] = numpy.polyval(exp_value_coefs, q[i])

    energy_guess = 1e-4
    energy_step = 1e-4
    num_nodes_last = 0

    expectation_value_reference = 0.0
    zero_point_energy = 0.0

    while num_nodes_last < input_data['num_solutions']:

        g = input_data['reduced_mass_amu'] * constants['amu_to_au'] * 2.0 * (pot_energies - energy_guess)

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

        if (abs(energy_step) < input_data['energy_precision']) and (num_nodes > num_nodes_last):

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

            expectation_value = 2.0 * numpy.dot(psi, exp_values * psi)
            psi_squared[num_nodes - 1] = psi**2.0

            energy[num_nodes - 1] = energy_guess

            if num_nodes == 1:
                expectation_value_reference = expectation_value
                zero_point_energy = energy_guess

            transition_frequency = (energy_guess - zero_point_energy) * constants['hartree_to_cm1']

            diff_au = (expectation_value - expectation_value_reference)
            diff_hz = diff_au * constants['hartree_to_hz']

            energy_step = 1.0e-4
            num_nodes_last = num_nodes_last + 1

        if num_nodes > num_nodes_last:
            if energy_step > 0.0:
                energy_step = energy_step / (-10.0)
        else:
            if energy_step < 0.0:
                energy_step = energy_step / (-10.0)

        energy_guess = energy_guess + energy_step

    return q, psi_squared, energy, diff_hz, transition_frequency


def main():
    version = '1.1.0'

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
    print('version            =', version)
    print()
    print('reduced mass       =', input_data['reduced_mass_amu'])
    print('harmonic frequency =', input_data['harmonic_frequency_cm1'])
    print('degree potential   =', input_data['degree_potential'])
    print('degree property    =', input_data['degree_property'])
    print('energy precision   =', input_data['energy_precision'])
    print('nr solutions       =', input_data['num_solutions'])
    print('nr steps           =', input_data['num_steps'])

    transition_frequency = sys.float_info.max
    transition_frequency_previous = -sys.float_info.max
    displacement_range = (-0.5, 0.5)
    while abs(transition_frequency - transition_frequency_previous) > 1.0e-1:
        q, psi_squared, energy, diff_hz, transition_frequency = run_numerov(input_data,
                                                                            constants,
                                                                            pot_energy_coefs,
                                                                            exp_value_coefs,
                                                                            displacement_range)
        displacement_range = (displacement_range[0] - 0.1, displacement_range[1] + 0.1)
        transition_frequency_previous = transition_frequency

    # get harmonic frequency from numerov
    _, _, _, _, transition_frequency_harmonic = run_numerov(input_data,
                                                            constants,
                                                            pot_energy_coefs_harmonic,
                                                            exp_value_coefs,
                                                            displacement_range)

    p0 = exp_value_coefs[-1] * constants['hartree_to_hz']
    p1 = exp_value_coefs[-2] * constants['hartree_to_hz']
    p2 = 2.0 * exp_value_coefs[-3] * constants['hartree_to_hz']
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
