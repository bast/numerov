import numpy
import matplotlib.pyplot as plt
import yaml
import sys
import math
import click


def read_data(input_file):
    with open(input_file, 'r') as stream:
        try:
            input_data = yaml.load(stream)
            return input_data
        except yaml.YAMLError as e:
            sys.stderr.write(e)
            sys.stderr.write("\n")
            sys.exit(-1)


def get_smooth_curve(x, coefs):
    x_min = min(x)
    x_max = max(x)
    step = (x_max - x_min)/100.0
    xs = []
    ys = []
    q = x_min
    while q < x_max:
        q += step
        xs.append(q)
        ys.append(numpy.polyval(coefs, q))
    return xs, ys


@click.command()
@click.option('--force-field', help='File containing force-field data.')
@click.option('--exp-values', help='File containing expectation values along displacement.')
@click.option('--numerov-output', help='Output from numerov.')
@click.option('--img', help='Image file name.')
def main(force_field, exp_values, numerov_output, img):

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    numerov_data = read_data(numerov_output)

    qs = numerov_data['qs']
    psi_squared = numerov_data['psi_squared']
    for i in range(2):
        ax1.plot(qs, psi_squared[i], 'g-')

    ax1.set_xlabel('displacement')
    ax1.set_ylabel('pot', color='b')
    ax1.tick_params('y', colors='b')
    ax2.set_ylabel('exp', color='r')
    ax2.tick_params('y', colors='r')

    data = read_data(force_field)
    x = numpy.array(data['displacements_au'])
    y = numpy.array(data['energies_hartree'])
    # shift potential such that minimum is at zero
    y -= min(y)
    ax1.plot(x, y, 'bx')
    coefs = numerov_data['pot_energy_coefs']
    for i in range(len(coefs)):
        coefs[-i - 1] /= math.factorial(i)
    xs, ys = get_smooth_curve(x, coefs)
    ax1.plot(xs, ys, 'b.')

    data = read_data(exp_values)
    x = numpy.array(data['displacements_au'])
    y = numpy.array(data['exp_values_hartree'])
    ax2.plot(x, y, 'rx')
    coefs = numerov_data['exp_value_coefs']
    for i in range(len(coefs)):
        coefs[-i - 1] /= math.factorial(i)
    xs, ys = get_smooth_curve(x, coefs)
    ax2.plot(xs, ys, 'r.')

    fig.tight_layout()

    plt.savefig(img, dpi=300.0)


if __name__ == '__main__':
    main()
