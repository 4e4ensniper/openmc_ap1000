import csv
import numpy as np
import openmc
import sys
from math import log
from math import pi
from math import sqrt
from math import asin

import openmc.model
from fuel_assemblies import fa_types, find_name

sys.path.append('../')
from constants import split_number, core_height, csv_path
from constants import g1, g2, g3, g4, g5, h1, h2, h3, h4, h5, b_conc, dif_fu_cart
from constants import treton_input_files, treton_files_path, numbers, fr_number
from constants import r_fr, delta_shell, l_shell, r_fuel, l_g, l_fuel, r_hole

#materials specifications
#we have 1 material complex for each split elemen
#function for determining the mean values ​​of quantities
#from their piecewise distributions
eps = 0.0001
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
        csvreader = csv.reader(csvfile, delimiter='\t')
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

def cr_steel(i, j, num, temp):
    _08x18h10t = openmc.Material(material_id = int(1E7 + num*1E5 + i*1E2 + j), name = "08x18h10t")
    _08x18h10t.add_element('Fe', 67.665,'wo')
    _08x18h10t.add_element('Cr', 18.0,'wo')
    _08x18h10t.add_element('Ni', 10.0,'wo')
    _08x18h10t.add_element('Ti', 0.6,'wo')
    _08x18h10t.add_element('C', 0.08,'wo')
    _08x18h10t.add_element('Si', 0.8,'wo')
    _08x18h10t.add_element('Mn', 2.0,'wo')
    _08x18h10t.add_element('Mo', 0.5,'wo')
    _08x18h10t.add_element('S', 0.02,'wo')
    _08x18h10t.add_element('P', 0.035,'wo')
    _08x18h10t.add_element('Cu', 0.3,'wo')
    _08x18h10t.set_density('g/cm3', 7.85)
    _08x18h10t.temperature = temp
    return _08x18h10t

def b4c(i, j, num, temp):
    b4c_ = openmc.Material(material_id = int(1E7 + num * 1E5 + i*1E2 + j), name = "b4c_absorber")
    b4c_.set_density('g/cm3', 1.8)
    b4c_.add_element('B', 4.0, enrichment=80.0, enrichment_target='B10')
    b4c_.add_element('C', 1.0)
    b4c_.temperature = temp
    return b4c_

def helium_(i, j, num, temp):
    helium=openmc.Material(material_id = int(1E7 + num * 1E5 + i*1E2 + j), name = "He")
    helium.add_element('He', 1.0)
    helium.temperature = temp
    helium.set_density('g/cm3', 3.24e-3)
    return helium
def uo2_fuel(i, j, num, temp, enrich):
    fu = openmc.Material(material_id = int(1E7 + num * 1E5 + i*1E2 + j), name = f"UO2_{enrich}")
    fu.add_element('U', 1.0, enrichment = enrich, enrichment_type='wo')
    fu.add_element('O', 2.0)
    fu.set_density('g/cm3', 10.48)
    fu.temperature = temp
    return fu
def uo2_gdo2(i, j, num, temp, enrich, gdo2_pt):
    uo2 = uo2_fuel(i, j, num, temp, enrich)
    gdo2 = openmc.Material(material_id = int(1E7 + (num + 1) * 1E5 + i*1E2 + j), name = 'GdO2')
    gdo2.add_element('Gd', 2.0)
    gdo2.add_element('O', 3.0)
    gdo2.set_density('g/cm3', 7.407)
    gdo2_uo2 = openmc.Material.mix_materials([uo2, gdo2], [1-gdo2_pt*1E-2, gdo2_pt*1E-2], 'wo')
    gdo2_uo2.id = int(1E7 + (num + 2) * 1E5 + i*1E2 + j)
    gdo2_uo2.name = 'GdO2_UO2'
    gdo2_uo2.temperature = temp
    return gdo2_uo2
def dif_pos(array):
    dif_pos = []
    dict_pos = {}
    for i, number in enumerate(array):
        if number in dict_pos:
            continue
        else:
            dif_pos.append(i)
            dict_pos[number] = i
    return dif_pos


fuel = []
gaz = []
shell = []
coolant = []
g_hole = []
grey_rods = []

if treton_input_files == 0:
    density_hc = average_value(split_number, csv_path + "density_hc_av(l).csv", core_height)
    #np.savetxt(f"density_av.csv", density_hc, delimiter=",", fmt='%.15f')
    fuel_temp_eff = average_value(split_number, csv_path + "fuel_temperature_eff.csv", core_height)
    #np.savetxt(f"fuel_temp_eff_av.csv", fuel_temp_eff, delimiter=",", fmt='%.15f')
    helium_temp = average_value(split_number, csv_path + "gaz_temperature_av(z).csv", core_height)
    shell_temp = average_value(split_number, csv_path + "shell_temperature_av_av(z).csv", core_height)
    hc_temp_arr = average_value(split_number, csv_path + "T_hc_av(z).csv", core_height)
    hole_helium_temp = average_value(split_number, csv_path + "fuel_temperature_max_av(z).csv", core_height)
else:
    need_lines = dif_pos(numbers)
    shell_res = (log(r_fr/(r_fr - delta_shell)))/(2*pi*l_shell)
    gaz_res = (log((r_fr - delta_shell)/r_fuel))/(2*pi*l_g)
    delta_gap = r_fr + delta_shell - r_fuel
    fuel_res = 1/(4*pi*l_fuel)*(1-(2*r_hole*r_hole)/(r_fuel*r_fuel - r_hole*r_hole)*log(r_fuel/r_hole))
    with open(treton_files_path + 't_r_obl_naruj.dat', 'r') as file:
        shell_lines = file.readlines()
    with open(treton_files_path + 'ro_tepl.dat', 'r') as file:
        ro_lines = file.readlines()
    with open(treton_files_path + 't_tepl.dat', 'r') as file:
        hc_lines = file.readlines()
    with open(treton_files_path + 'Q6.txt', 'r') as file:
        q_lines = file.readlines()

    #materials ID definition
    #1??????? - first 1
    #1__????? - material number
    #1??___?? - fuel asssembly number
    #1?????__ - split number from bottom

for i in range(0, len(dif_fu_cart)):
    if treton_input_files == 1:
        out_shell_temp_line = []
        density_hc_line = []
        hc_temp_line = []
        q_line = []
        out_shell_temp_ = shell_lines[need_lines[i]].strip().split()
        out_shell_temp_line.extend([float(temp) for temp in out_shell_temp_[2:]])

        density_hc_ = ro_lines[need_lines[i]].strip().split()
        density_hc_line.extend([float(dens) for dens in density_hc_[2:]])

        hc_temp_arr_ = hc_lines[need_lines[i]].strip().split()
        hc_temp_line.extend([float(temp) for temp in hc_temp_arr_[2:]])

        q_line_ = q_lines[need_lines[i]].strip().split()
        q_line.extend(float(q) for q in q_line_)

    type = find_name(dif_fu_cart[i], fa_types)
    grey_rod = []
    for j in range(0, split_number):
        if treton_input_files == 1:
            q_l = q_line[j] * split_number/(fr_number * core_height)
            central_gaz_temp = out_shell_temp_line[j] + q_l * (shell_res + gaz_res + fuel_res)
            out_fuel_temp = out_shell_temp_line[j] + q_l * (shell_res + gaz_res)
            in_shell_temp = out_shell_temp_line[j] + q_l * shell_res
            shell_temp_av = 1/delta_shell*(in_shell_temp - out_shell_temp_line[j])/log(1-delta_shell/r_fr)*(r_fr-delta_shell)*(r_fr/(r_fr-delta_shell)*log(r_fr/(r_fr-delta_shell))-r_fr/(r_fr-delta_shell)+1)+in_shell_temp
            gaz_gap_temp_av = 1/delta_gap*(out_fuel_temp - in_shell_temp)/log(r_fuel/(r_fuel+delta_gap))*r_fuel*((r_fuel+delta_gap)/r_fuel*log((r_fuel+delta_gap)/r_fuel)-(r_fuel+delta_gap)/r_fuel+1)+out_fuel_temp

            b_ = q_l/(4*pi*l_fuel*(r_fuel*r_fuel - r_hole*r_hole))
            a_ = out_fuel_temp + b_*(r_fuel*r_fuel - 2*r_hole*r_hole*log(r_fuel/r_hole))
            t_eff_n = 1/sqrt(b_) * (sqrt(b_)*r_fuel/2*sqrt(a_ - b_*r_fuel*r_fuel) + a_/2 * asin(sqrt(b_/a_)*r_fuel) \
                                    - sqrt(b_)*r_hole/2*sqrt(a_ - b_*r_hole*r_hole) - a_/2 * asin(sqrt(b_/a_)*r_hole))
            t_eff_d = 1/sqrt(b_)*(asin(sqrt(b_/a_)*r_fuel) - asin(sqrt(b_/a_)*r_hole))
            t_eff = t_eff_n/t_eff_d

            hc_density = density_hc_line[j]
            hc_temp = hc_temp_line[j]

        else:
            central_gaz_temp = hole_helium_temp[j]
            shell_temp_av = shell_temp[j]
            gaz_gap_temp_av = helium_temp[j]
            t_eff = fuel_temp_eff[j]
            hc_density = density_hc[j]
            hc_temp = hc_temp_arr[j]

        #helium definition in central hole
        g_hole.append(helium_(i, j, 1, central_gaz_temp + 273.15))

        #fuel definition
        fu = uo2_fuel(i, j, 2, t_eff + 273.15, type["enrichment"])
        fuel.append(fu)
        if type["gdo2_wo"] != 0:
            grey = uo2_gdo2(i, j, 3, t_eff + 273.15, type["grey_enrichment"], type["gdo2_wo"])
            grey_rod.append(grey)

        #helium definition in gap
        gaz.append(helium_(i, j, 6, gaz_gap_temp_av + 273.15))
        #shell definition
        shell_alloy = openmc.Material(material_id = int(1E7 + 7E5 + i*1E2 + j), name = "110")
        shell_alloy.add_element('Zr', 0.99, percent_type='wo')
        shell_alloy.add_element('Nb',0.1, percent_type='wo')
        shell_alloy.temperature = shell_temp_av + 273.15
        shell_alloy.set_density('g/cm3',6.5)
        shell.append(shell_alloy)

        #coolant definition
        if b_conc > eps:
            b_ppm = 1/(1 + 61.83/18 * (1/(b_conc*1E-3)-1)) * 1E6
            water = openmc.model.borated_water(boron_ppm = b_ppm, density=hc_density*1E-3)
            water.id = int(1E7 + 8E5 + i*1E2 + j)
            water.temperature = hc_temp + 273.15
            water.name = 'H2O_b'
        else:
            water = openmc.Material(material_id = int(1E7 + 9E5 + i*1E2 + j), name = "H2O")
            water.add_element('H', 2.0)
            water.add_element('O', 1.0)
            water.set_density('g/cm3', hc_density*1E-3)
            water.temperature = hc_temp + 273.15
            water.add_s_alpha_beta('c_H_in_H2O')
        coolant.append(water)
    if type["gdo2_wo"] != 0:
        grey_rods.append(grey_rod)
cr_shell1 = []
boron_carbide1 = []

cr_shell2 = []
boron_carbide2 = []

cr_shell3 = []
boron_carbide3= []

cr_shell4 = []
boron_carbide4 = []

cr_shell5 = []
boron_carbide5 = []
for i in range(0, len(g1)):
    hc_temp_line_g1 = []
    if treton_input_files == 1:
        hc_temp_arr_ = hc_lines[need_lines[g1[i]-1]].strip().split()
        hc_temp_line_g1.extend([float(temp) for temp in hc_temp_arr_[2:]])
    else:
        hc_temp_line_g1 = hc_temp_arr
    for j in range(0, h1):
        cr_shell1.append(cr_steel(g1[i], j, 10, hc_temp_line_g1[split_number - h1 + j] + 273.15))
        boron_carbide1.append(b4c(g1[i], j, 11, hc_temp_line_g1[split_number - h1 + j] + 273.15))
for i in range(0, len(g2)):
    hc_temp_line_g2 = []
    if treton_input_files == 1:
        hc_temp_arr_ = hc_lines[need_lines[g2[i]-1]].strip().split()
        hc_temp_line_g2.extend([float(temp) for temp in hc_temp_arr_[2:]])
    else:
        hc_temp_line_g2 = hc_temp_arr
    for j in range(0, h2):
        cr_shell2.append(cr_steel(g2[i], j, 10, hc_temp_line_g2[split_number - h2 + j] + 273.15))
        boron_carbide2.append(b4c(g2[i], j, 11, hc_temp_line_g2[split_number - h2 + j] + 273.15))
for i in range(0, len(g3)):
    hc_temp_line_g3 = []
    if treton_input_files == 1:
        hc_temp_arr_ = hc_lines[need_lines[g3[i]-1]].strip().split()
        hc_temp_line_g3.extend([float(temp) for temp in hc_temp_arr_[2:]])
    else:
        hc_temp_line_g3 = hc_temp_arr
    for j in range(0, h3):
        cr_shell3.append(cr_steel(g3[i], j, 10, hc_temp_line_g3[split_number - h3 + j] + 273.15))
        boron_carbide3.append(b4c(g3[i], j, 11, hc_temp_line_g3[split_number - h3 + j] + 273.15))
for i in range(0, len(g4)):
    hc_temp_line_g4 = []
    if treton_input_files == 1:
        hc_temp_arr_ = hc_lines[need_lines[g4[i]-1]].strip().split()
        hc_temp_line_g4.extend([float(temp) for temp in hc_temp_arr_[2:]])
    else:
        hc_temp_line_g4 = hc_temp_arr
    for j in range(0, h4):
        cr_shell4.append(cr_steel(g4[i], j, 10, hc_temp_line_g4[split_number - h4 + j] + 273.15))
        boron_carbide4.append(b4c(g4[i], j, 11, hc_temp_line_g4[split_number - h4 + j] + 273.15))
for i in range(0, len(g5)):
    hc_temp_line_g5 = []
    if treton_input_files == 1:
        hc_temp_arr_ = hc_lines[need_lines[g5[i]-1]].strip().split()
        hc_temp_line_g5.extend([float(temp) for temp in hc_temp_arr_[2:]])
    else:
        hc_temp_line_g5 = hc_temp_arr
    for j in range(0, h5):
        cr_shell5.append(cr_steel(g5[i], j, 10, hc_temp_line_g5[split_number - h5 + j] + 273.15))
        boron_carbide5.append(b4c(g5[i], j, 11, hc_temp_line_g5[split_number - h5 + j] + 273.15))
if treton_input_files == 1:
    t_water_hc_outher = []
    t_water_hc_outher_ = hc_lines[0].strip().split()
    t_water_hc_outher.extend([float(temp) for temp in t_water_hc_outher_[2:]])
    av_temp_water = sum(t_water_hc_outher)/len(t_water_hc_outher)

    d_water_hc_outher = []
    d_water_hc_outher_ = ro_lines[0].strip().split()
    d_water_hc_outher.extend([float(temp) for temp in d_water_hc_outher_[2:]])
    av_d_water = sum(d_water_hc_outher)/len(d_water_hc_outher)
else:
    av_temp_water = sum(hc_temp_arr)/len(hc_temp_arr)
    av_d_water = sum(density_hc)/len(density_hc)
if b_conc > eps:
    b_ppm = 1/(1 + 61.83/18 * (1/(b_conc*1E-3)-1)) * 1E6
    water = openmc.model.borated_water(boron_ppm = b_ppm, density=av_d_water*1E-3)
    water.id = int(1E7 + 12E5 + len(dif_fu_cart)*1E2 + split_number)
    water.temperature = av_temp_water + 273.15
    water.name = 'H2O'
else:
    water = openmc.Material(material_id = int(1E7 + 12E5 + len(dif_fu_cart)*1E2 + split_number), name="H2O")
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.set_density('g/cm3', av_d_water*1E-3)
    water.temperature = av_temp_water + 273.15
    water.add_s_alpha_beta('c_H_in_H2O')
coolant.append(water)



#reactor vessel steel SA-508
steel_all = openmc.Material(material_id = int(1E7 + 13E5 + len(dif_fu_cart)*1E2 + split_number),name="SA508")
steel_all.add_element('C', 0.206,'wo')
steel_all.add_element('Si', 0.321, 'wo')
steel_all.add_element('Mn', 1.280, 'wo')
steel_all.add_element('P', 0.002, 'wo')
steel_all.add_element('S', 0.014, 'wo')
steel_all.add_element('Cr', 0.140, 'wo')
steel_all.add_element('Ni', 0.644, 'wo')
steel_all.add_element('Fe', 96.9, 'wo')
steel_all.add_element('Mo', 0.493, 'wo')
steel_all.temperature = 285.7 + 273.15
steel_all.set_density('g/cm3', 8.05)
shell.append(steel_all)
print("materials created!")
