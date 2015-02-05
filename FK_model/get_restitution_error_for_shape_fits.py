# Be sure to run this with Python 2.7 or later so that the subprocess.check_output function exists
import subprocess

voltage_solutions_filename = 'evaluatedSolutions_5p0_MAP_fixed_offset_21param_99bounds_500x200.txt'
restitution_GA_executable = './FK4V_5p0_cell_GA_restitution'
voltage_restitution_solutions_filename = 'evaluatedSolutions_5p0_MAP_fixed_offset_21param_99bounds_500x200_PLUS_restitution_error.txt'

with open(voltage_solutions_filename, 'r') as f_in:
  with open(voltage_restitution_solutions_filename,'a') as f_out: 
    individual = f_in.readline()
    for individual in f_in:
      individual = individual.split('\t')
      voltage_error = float(individual[21])
      parameters = individual[:21]
      parameters.append('delete.dat')
      parameters = ' '.join(parameters)
      if (voltage_error < 0.06):
        # print voltage_error
        # print (restitution_GA_executable + ' ' + parameters + ' delete.dat')
        restitution_error = subprocess.check_output([restitution_GA_executable, parameters])
        print restitution_error
        f_out.write('\t'.join(individual[:21]) + '\t' + str(voltage_error) + '\t' + restitution_error)
