import numpy as np
import matplotlib.pyplot as plt
import yaml
import sys
import math


def read_data(input_file):
    with open(input_file, 'r') as stream:
        try:
            input_data = yaml.load(stream)
            return input_data
        except yaml.YAMLError as e:
            sys.stderr.write(e)
            sys.stderr.write("\n")
            sys.exit(-1)


input_file = sys.argv[-1]
data = read_data(input_file)

displacements = []
exp_values = []
pot_energies = []

for x in data['steps']:
    displacements.append(x['displacement'])
    exp_values.append(x['exp_value'])
    pot_energies.append(x['pot_energy'])

pot_energies = [p - min(pot_energies) for p in pot_energies]


exp_value_coefs = data['exp_value_coefs']
for i in range(len(exp_value_coefs)):
    exp_value_coefs[-i - 1] /= math.factorial(i)

pot_energy_coefs = data['pot_energy_coefs']
for i in range(len(pot_energy_coefs)):
    pot_energy_coefs[-i - 1] /= math.factorial(i)

pot_energy_coefs_harmonic = data['pot_energy_coefs_harmonic']
for i in range(len(pot_energy_coefs_harmonic)):
    pot_energy_coefs_harmonic[-i - 1] /= math.factorial(i)


q = min(displacements)
step = (max(displacements) - min(displacements))/100.0
qs = []
es = []
ps = []
ps_h = []
while q < max(displacements):
    q += step
    qs.append(q)
    es.append(np.polyval(exp_value_coefs, q))
    ps.append(np.polyval(pot_energy_coefs, q))
    ps_h.append(np.polyval(pot_energy_coefs_harmonic, q))


fig, ax1 = plt.subplots()

ax1.plot(displacements, pot_energies, 'bx')
ax1.plot(qs, ps, 'b.')
ax1.plot(qs, ps_h, 'b-')
ax1.set_xlabel('displacement')
ax1.set_ylabel('pot', color='b')
ax1.tick_params('y', colors='b')

for i in range(2):
    ax1.plot(data['qs'], data['psi_squared'][i], 'g-')

ax2 = ax1.twinx()
ax2.plot(displacements, exp_values, 'rx')
ax2.plot(qs, es, 'r.')
ax2.set_ylabel('exp', color='r')
ax2.tick_params('y', colors='r')

fig.tight_layout()
plt.show()
