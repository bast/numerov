import numpy


def solve_numerov(pot_energy_coefs,
                  exp_value_coefs,
                  displacement_range,
                  num_steps,
                  num_solutions,
                  energy_precision_hartree,
                  reduced_mass_au):

    n = num_steps + 1
    step = (displacement_range[1] - displacement_range[0]) / num_steps
    step2 = step**2.0

    qs = numpy.zeros(n)
    pot_energies = numpy.zeros(n)
    exp_values = numpy.zeros(n)
    g = numpy.zeros(n)
    psi = numpy.zeros(n)
    psi_squared = numpy.zeros([num_solutions, n])
    energies_hartree = numpy.zeros(num_solutions)
    averaged_exp_values_au = numpy.zeros(num_solutions)

    for i in range(n):
        qs[i] = displacement_range[0] + i * step
        pot_energies[i] = numpy.polyval(pot_energy_coefs, qs[i])
        exp_values[i] = numpy.polyval(exp_value_coefs, qs[i])

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

        # the wave function node will not match the last point
        # and right of this point there is some numerical leftover
        # here we cut the leftover away
        num_nodes = 0
        i_save = n
        for i in range(1, n):
            if psi[i - 1] != 0.0:
                if (psi[i] / psi[i - 1]) < 0.0:
                    num_nodes += 1
                    i_save = i
        psi[i_save + 1:] = 0.0

        # normalize the wave function
        psi = psi / numpy.sqrt(numpy.dot(psi, psi))

        if (abs(energy_step) < energy_precision_hartree) and (num_nodes > num_nodes_last):
            # we have found a solution and we move on to searching the next solution
            psi_squared[num_nodes - 1] = psi * psi
            energies_hartree[num_nodes - 1] = energy_guess
            averaged_exp_values_au[num_nodes - 1] = numpy.dot(psi, exp_values * psi)
            energy_step = 1.0e-4
            num_nodes_last += 1

        if num_nodes > num_nodes_last:
            if energy_step > 0.0:
                energy_step /= (-10.0)
        else:
            if energy_step < 0.0:
                energy_step /= (-10.0)

        energy_guess += energy_step

    return qs, psi_squared, energies_hartree, averaged_exp_values_au


def test_solve_numerov():
    pot_energy_coefs = numpy.array([-2.45560869e-01, -8.88252151e-03,
                                    1.24439946e-01, 1.93259856e-01, 2.78860663e-01, -5.62738650e-05,
                                    -5.78784571e-08])
    exp_value_coefs = numpy.array([-1.45680171e-14, -2.78094078e-16,
                                   1.62725432e-15, -1.99732822e-15, 3.08772558e-15, 2.06211298e-15,
                                   -7.30656049e-15])
    displacement_range = (-0.5, 0.5)
    num_steps = 21
    num_solutions = 3
    energy_precision_hartree = 1e-12
    reduced_mass_au = 26245.03
    q, psi_squared, energies_hartree, averaged_exp_values_au = solve_numerov(pot_energy_coefs,
                                                                             exp_value_coefs,
                                                                             displacement_range,
                                                                             num_steps,
                                                                             num_solutions,
                                                                             energy_precision_hartree,
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
