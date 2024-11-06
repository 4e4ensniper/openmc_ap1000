import csv
import numpy as np
from math import pi
import openmc
import sys
sys.path.append('../')
from constants import split_number, core_height, r_fuel, r_hole, fr_number, n_fa, csv_path

#materials specifications
#we have 1 material complex for each split elemen

#function for determining the mean values ​​of quantities
#from their piecewise distributions
def average_value(n, filename, length):
    z_coordinates_in = []
    f_z_in = []
    z_coordinates_added = []
    f_z = []
    z_coordinates_ns = []
    z_coordinates = []
    eps = 1E-10
    # reading csv with water density
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            z_coordinates_in.append(float(row[0]))
            f_z_in.append(float(row[1]))
    d_z = length/n
    z_coordinates_added.append(float(-length/2 + d_z))
    for i in range(1, n-1):
        z_coordinates_added.append(z_coordinates_added[i-1]+d_z)
    z_coordinates_ns = np.array(z_coordinates_in + z_coordinates_added)
    z_coordinates = np.sort(z_coordinates_ns)
    f_z= np.interp(z_coordinates, z_coordinates_in, f_z_in)
    j = 1
    f_av = []
    z_coordinates_added.append(length/2)
    for i in range(0, n):
        sum_f_z = 0
        while z_coordinates[j] - z_coordinates_added[i] < eps:
            sum_f_z += (f_z[j-1] + f_z[j])*(z_coordinates[j] - z_coordinates[j-1])
            j += 1
            if j == len(z_coordinates):
                break
        f_av.append(sum_f_z/(2*d_z))
    return f_av

density_hc = []
fuel_temp_eff = []
helium_temp = []
shell_temp = []
hc_temp = []
hole_helium_temp = []
density_hc = average_value(split_number, csv_path + "density_hc_av(l).csv", core_height)
#np.savetxt(f"density_av.csv", density_hc, delimiter=",", fmt='%.15f')
fuel_temp_eff = average_value(split_number, csv_path + "fuel_temperature_eff(z).csv", core_height)
#np.savetxt(f"fuel_temp_eff_av.csv", fuel_temp_eff, delimiter=",", fmt='%.15f')
helium_temp = average_value(split_number, csv_path + "gaz_temperature_av(z).csv", core_height)
shell_temp = average_value(split_number, csv_path + "shell_temperature_av_av(z).csv", core_height)
hc_temp = average_value(split_number, csv_path + "T_hc_av(z).csv", core_height)
hole_helium_temp = average_value(split_number, csv_path + "fuel_temperature_av(z).csv", core_height)

fuel = []
gaz = []
shell = []
coolant = []
g_hole = []

#materials ID definition
#1??????? - first 1
#1__????? - material number
#1??___?? - fuel asssembly number
#1?????__ - split number from bottom

for i in range(0, n_fa):
    for j in range(0, split_number):

        #helium definition in central hole
        c_helium=openmc.Material(material_id = int(1E7 + 1E5 + i*1E2 + j), name = "He")
        c_helium.add_element('He', 1.0)
        c_helium.temperature = hole_helium_temp[j] + 273.15
        c_helium.set_density('g/cm3', 3.24e-3)
        g_hole.append(c_helium)

        #fuel definition
        fu = openmc.Material(material_id = int(1E7 + 2E5 + i*1E2 + j), name = "UO2")
        fu.add_element('U', 1.0, enrichment = 4.95)
        fu.add_element('O', 2.0)
        fu.set_density('g/cm3', 10.48)
        fu.temperature = fuel_temp_eff[j] + 273.15
        fu.volume = fr_number*pi*(r_fuel*r_fuel - r_hole*r_hole)*core_height*1E2/split_number#объем топлива в одном элементе разбиения ТВC по высоте
        fu.depletable = True
        fuel.append(fu)

        #helium definition in gap
        helium=openmc.Material(material_id = int(1E7 + 3E5 + i*1E2 + j), name = "He")
        helium.add_element('He', 1.0)
        helium.temperature = helium_temp[j] + 273.15
        helium.set_density('g/cm3', 3.24e-3)
        gaz.append(helium)

        #shell definition
        shell_alloy = openmc.Material(material_id = int(1E7 + 4E5 + i*1E2 + j), name = "110")
        shell_alloy.add_element('Zr', 0.99, percent_type='wo')
        shell_alloy.add_element('Nb',0.1, percent_type='wo')
        shell_alloy.temperature = shell_temp[j] + 273.15
        shell_alloy.set_density('g/cm3',6.5)
        shell.append(shell_alloy)

        #coolant definition
        water = openmc.Material(material_id = int(1E7 + 5E5 + i*1E2 + j), name = "H2O")
        water.add_element('H', 2.0)
        water.add_element('O', 1.0)
        water.set_density('g/cm3', density_hc[j]*1E-3)
        water.temperature = hc_temp[j] + 273.15
        water.add_s_alpha_beta('c_H_in_H2O')
        coolant.append(water)

#add water "fuel assemblys" for hexagonal form of core
for i in range(n_fa, n_fa + 18):
    for j in range (0, split_number):
        #coolant definition
        water = openmc.Material(material_id = int(1E7 + 5E5 + i*1E2 + j), name="H2O")
        water.add_element('H', 2.0)
        water.add_element('O', 1.0)
        water.set_density('g/cm3', density_hc[j]*1E-3)
        water.temperature = hc_temp[j] + 273.15
        water.add_s_alpha_beta('c_H_in_H2O')
        coolant.append(water)

#reactor vessel steel SA-508
steel_all = openmc.Material(name="SA508")
steel_all.add_element('C', 0.00206,'wo')
steel_all.add_element('Si', 0.00321, 'wo')
steel_all.add_element('Mn', 0.01280, 'wo')
steel_all.add_element('P', 0.00002, 'wo')
steel_all.add_element('S', 0.00014, 'wo')
steel_all.add_element('Cr', 0.00140, 'wo')
steel_all.add_element('Ni', 0.00644, 'wo')
steel_all.add_element('Fe', 0.969, 'wo')
steel_all.add_element('Mo', 0.00493, 'wo')
steel_all.temperature = 280.7 + 273.15
steel_all.set_density('g/cm3', 8.05)
shell.append(steel_all)

print("materials created!")
