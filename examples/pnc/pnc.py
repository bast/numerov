import sys
import yaml


constants = {'amu_to_au': 1822.888479031408,
             'hartree_to_hz': 6.579683920721e15,
             'hartree_to_cm1': 2.194746313705e5}


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

    exp_value_coefs = input_data['exp_value_coefs']
    pot_energy_coefs = input_data['pot_energy_coefs']
    averaged_exp_values_au = input_data['averaged_exp_values_au']
    energies_hartree = input_data['energies_hartree']
    energies_hartree_harmonic = input_data['energies_hartree_harmonic']

    p0 = exp_value_coefs[-1] * constants['hartree_to_hz']
    p1 = exp_value_coefs[-2] * constants['hartree_to_hz']
    p2 = exp_value_coefs[-3] * constants['hartree_to_hz']
    v3 = pot_energy_coefs[-4]

    mass = input_data['reduced_mass_amu'] * constants['amu_to_au']
    frequency = input_data['harmonic_frequency_cm1'] / constants['hartree_to_cm1']

    delta_2 = p2
    delta_3 = delta_2 - p1 * v3 / (mass * frequency * frequency)

    to_level = 1
    delta_2 *= to_level / (mass * frequency)
    delta_3 *= to_level / (mass * frequency)

    diff_au = 2.0 * (averaged_exp_values_au[1] - averaged_exp_values_au[0])
    diff_hz = diff_au * constants['hartree_to_hz']

    transition_frequency = (energies_hartree[1] - energies_hartree[0]) * constants['hartree_to_cm1']
    transition_frequency_harmonic = (energies_hartree_harmonic[1] - energies_hartree_harmonic[0]) * constants['hartree_to_cm1']

    print('p0_hz: {0}'.format(p0))
    print('p1_hz: {0}'.format(p1))
    print('p2_hz: {0}'.format(p2))
    print('v3_hartree: {0}'.format(v3))
    print('delta2_hz: {0}'.format(delta_2))
    print('delta3_hz: {0}'.format(delta_3))
    print('numerov_hz: {0}'.format(diff_hz))
    print('freq_cm1: {0}'.format(transition_frequency))
    print('freq_harmonic_cm1: {0}'.format(transition_frequency_harmonic))


if __name__ == '__main__':
    main()
