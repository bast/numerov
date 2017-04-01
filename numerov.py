version = '1.1.0'

amu_to_au      = 1822.888479031408
hartree_to_hz  = 6.579683920721e15
hartree_to_cm1 = 2.194746313705e5


import sys
import string
import os
import numpy
#import plot_numerov

from optparse import OptionParser

#-------------------------------------------------------------------------------

usage = '''
  example: ./%prog --file=CHFClBr_mode3_x2c_b3lyp'''

parser = OptionParser(usage)

parser.add_option('--file',
                  type='string',
                  action='store',
                  default=None,
                  help='File to parse',
                  metavar='FILE')
parser.add_option('--degree_potential',
                  type='int',
                  action='store',
                  default=6)
parser.add_option('--degree_property',
                  type='int',
                  action='store',
                  default=6)
parser.add_option('--energy_precision',
                  type='float',
                  action='store',
                  default=1.0e-12)
parser.add_option('--nr-solutions',
                  type='int',
                  action='store',
                  default=2)
parser.add_option('--nr-steps',
                  type='int',
                  action='store',
                  default=1000)
parser.add_option('--to-level',
                  type='int',
                  action='store',
                  default=1)

(options, args) = parser.parse_args()

if len(sys.argv) == 1:
    # user has given no arguments: print help and exit
    print parser.format_help().strip()
    sys.exit()


l = []


import yaml

with open(options.file, 'r') as stream:
    try:
        f = yaml.load(stream)
    except yaml.YAMLError as e:
        print(e)


harmonic_frequency_cm1 = f['harmonic_frequency_cm1']
reduced_mass_amu = f['reduced_mass_amu']

l = []
for step in f['steps']:
    l.append((step['displacement'], step['potential'], step['property']))
l.sort()

q_l = []
potential_l = []
property_l = []
for x in l:
    q_l.append(x[0])
    potential_l.append(x[1])
    property_l.append(x[2])

#-------------------------------------------------------------------------------

def run_numerov(potential_coef,
                property_coef,
                q_max):

    q_min = -q_max

    n     = options.nr_steps + 1
    step  = (q_max - q_min)/options.nr_steps
    step2 = step**2.0

    q           = numpy.zeros(n)
    potential   = numpy.zeros(n)
    property    = numpy.zeros(n)
    g           = numpy.zeros(n)
    psi         = numpy.zeros(n)
    psi_squared = numpy.zeros([options.nr_solutions, n])
    energy      = numpy.zeros(options.nr_solutions)

    for i in range(n):
        q[i]         = q_min + i*step
        potential[i] = numpy.polyval(potential_coef, q[i])
        property[i]  = numpy.polyval(property_coef,  q[i])

    energy_guess  = 1e-4
    energy_step   = 1e-4
    nr_nodes_last = 0

    expectation_value_reference = 0.0
    zero_point_energy           = 0.0

    while nr_nodes_last < options.nr_solutions:

          g = reduced_mass_amu*amu_to_au*2.0*(potential - energy_guess)

          psi[0]   = 0.0
          psi[1]   = 1.0e-6

          for i in range(n):
              if i > 1:
                 psi[i] = (2.0*psi[i-1] - psi[i-2] + g[i-1]*psi[i-1]*step2*5.0/6.0 \
                        + g[i-2]*psi[i-2]*step2/12.0)/(1.0 - g[i]*step2/12.0)

          nr_nodes = 0
          i_save   = n
          for i in range(n):
              if i > 1:
                  if psi[i-1] != 0.0:
                     if (psi[i]/psi[i-1]) < 0.0:
                        nr_nodes = nr_nodes + 1
                        i_save = i

          psi[i_save+1:] = 0.0

          psi = psi/numpy.sqrt(numpy.dot(psi, psi))

          if (abs(energy_step) < options.energy_precision) and (nr_nodes > nr_nodes_last):

             norm = 0.0
             done = 0
             for i in range(len(psi)):
                 norm = norm + psi[i]*psi[i]
                 if (norm > 0.001) and not done:
                    i_left = i
                    done   = 1

             norm = 0.0
             done = 0
             for i in range(n-1, -1, -1):
                 norm = norm + psi[i]*psi[i]
                 if (norm > 0.001) and not done:
                    i_right = i
                    done    = 1

             expectation_value = 2.0*numpy.dot(psi, property*psi)
             psi_squared[nr_nodes-1] = psi**2.0

             energy[nr_nodes-1] = energy_guess

             if nr_nodes == 1:
                expectation_value_reference = expectation_value
                zero_point_energy = energy_guess

             transition_frequency = (energy_guess - zero_point_energy)*hartree_to_cm1

             diff_au = (expectation_value - expectation_value_reference)
             diff_hz = diff_au*hartree_to_hz

             energy_step   = 1.0e-4
             nr_nodes_last = nr_nodes_last + 1

          if nr_nodes > nr_nodes_last:
             if energy_step > 0.0:
                energy_step = energy_step/(-10.0)
          else:
             if energy_step < 0.0:
                energy_step = energy_step/(-10.0)

          energy_guess = energy_guess + energy_step

    return q, psi_squared, energy, diff_hz, transition_frequency

#-------------------------------------------------------------------------------

potential_x  = numpy.array(q_l[:])
potential_y  = numpy.array(potential_l[:])
potential_y -= min(potential_y)

potential_coef = numpy.polyfit(potential_x, potential_y, options.degree_potential)

property_x = numpy.array(q_l[:])
property_y = numpy.array(property_l[:])

property_coef = numpy.polyfit(property_x, property_y, options.degree_property)

copy_y = []
for y in potential_y:
    copy_y.append(y)
copy_y.sort()
fourth_lowest_energy = copy_y[3]
h_x = []
h_y = []
for i in range(len(potential_x)):
    if potential_y[i] < fourth_lowest_energy:
        h_x.append(potential_x[i])
        h_y.append(potential_y[i])
potential_coef_harmonic = numpy.polyfit(h_x, h_y, 2)

#-------------------------------------------------------------------------------

print
print 'version            =', version
print
print 'reduced mass       =', reduced_mass_amu
print 'harmonic frequency =', harmonic_frequency_cm1
print 'degree potential   =', options.degree_potential
print 'degree property    =', options.degree_property
print 'energy precision   =', options.energy_precision
print 'nr solutions       =', options.nr_solutions
print 'nr steps           =', options.nr_steps

transition_frequency_previous = 0.0

q_max = 0.5
while True:
    q, psi_squared, energy, diff_hz, transition_frequency = run_numerov(potential_coef, property_coef, q_max)
    if abs(transition_frequency - transition_frequency_previous) < 1.0e-1:
        break
    q_max += 0.1
    transition_frequency_previous = transition_frequency

# get harmonic frequency from numerov
x1, x2, x3, x4, transition_frequency_harmonic = run_numerov(potential_coef_harmonic, property_coef, q_max)

p0        = property_coef[-1]*hartree_to_hz
p1        = property_coef[-2]*hartree_to_hz
p2        = 2.0*property_coef[-3]*hartree_to_hz
v3        = 6.0*potential_coef[-4]
mass      = reduced_mass_amu*amu_to_au
frequency = harmonic_frequency_cm1/hartree_to_cm1

n = options.to_level
delta_2 = p2
delta_3 = delta_2 - p1*v3/(mass*frequency*frequency)
delta_2 *= n/(mass*frequency)
delta_3 *= n/(mass*frequency)

print
print '%12s %10s %10s %10s %10s %10s %10s %7s %7s %4s' % ('p(0)', 'p(1)', 'p(2)', 'v(3)',
                                                     'delta(2)', 'delta(3)', 'numerov', 'freq', 'freq h', 'q')
print '%12.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %7.1f %7.1f %4.1f' % (p0, p1, p2, v3,
                                                                             delta_2, delta_3,
                                                                              diff_hz,
                                                                              transition_frequency,
                                                                              transition_frequency_harmonic,
                                                                              q_max)

plot_file = string.replace(options.file, '.in', '')

#plot_numerov.plot_numerov(potential_x,
#                          potential_y,
#                          potential_coef,
#                          property_x,
#                          property_y,
#                          property_coef,
#                          q,
#                          psi_squared,
#                          energy,
#                          'Q / a.u.',
#                          'V(Q) / E~B~h~N~',
#                          'E~B~PV~N~(Q) / E~B~h~N~',
#                          plot_file)
#
#os.system('convert -trim %s.ps %s.jpg' % (plot_file, plot_file))
