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


@click.command()
@click.option('--numerov-output', help='Output from numerov.')
@click.option('--img', help='Image file name.')
def main(numerov_output, img):

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    numerov_data = read_data(numerov_output)
    qs = numerov_data['qs']
    psi_squared = numerov_data['psi_squared']


    colors = ['r', 'g', 'b', 'y']

    for i in range(len(psi_squared)):
        color = colors[i%len(colors)]
        ax1.plot(qs, [float(x) for x in psi_squared[i]], f'{color}-')

    plt.savefig(img, dpi=300.0)


if __name__ == '__main__':
    main()
