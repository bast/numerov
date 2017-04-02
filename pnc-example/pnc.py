    p0 = exp_value_coefs[-1] * constants['hartree_to_hz']
    p1 = exp_value_coefs[-2] * constants['hartree_to_hz']
    p2 = exp_value_coefs[-3] * constants['hartree_to_hz']
    v3 = pot_energy_coefs[-4]

    mass = input_data['reduced_mass_amu'] * constants['amu_to_au']
    frequency = input_data['harmonic_frequency_cm1'] / constants['hartree_to_cm1']

    delta_2 = p2
    delta_3 = delta_2 - p1 * v3 / (mass * frequency * frequency)

    n = input_data['to_level']
    delta_2 *= n / (mass * frequency)
    delta_3 *= n / (mass * frequency)

    print(yaml.dump(input_data))
    print('displacement_range:')
    print('  min: {0}'.format(displacement_range[0]))
    print('  max: {0}'.format(displacement_range[1]))
    print()
    print('%12s %10s %10s %10s %10s %10s %10s %7s %7s %4s' % ('p(0)', 'p(1)', 'p(2)', 'v(3)',
                                                              'delta(2)', 'delta(3)', 'numerov', 'freq', 'freq h', 'q'))
    print('%12.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %7.1f %7.1f %4.1f' % (p0, p1, p2, v3,
                                                                                  delta_2, delta_3,
                                                                                  diff_hz,
                                                                                  transition_frequency,
                                                                                  transition_frequency_harmonic,
                                                                                  displacement_range[1]))


if __name__ == '__main__':
    main()
