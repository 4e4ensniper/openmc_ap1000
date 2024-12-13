import openmc
import openmc.deplete
import sys
import numpy as np
import os

from math import sqrt
from collections import Counter

#os.environ["OPENMC_CROSS_SECTIONS"] = "/home/adminsrv/projects/sections/endfb71_hdf5/cross_sections.xml"

sys.path.append('../')
from constants import n_fa, turnkey_size, core_barrel_in_r, core_barrel_out_r, core_height, split_number, line
from constants import batches, particles, inactive, numbers, g1, g2, g3, g4, g5, h1, h2, h3, h4, h5
from constants import dif_fu_cart, q_r

sys.path.append('../'+'fuel_assembly')
from fuel_assembly import full_fa, water_full_fa

sys.path.append('../'+'materials')
from materials import g_hole, fuel, gaz, shell, coolant, cr_shell1, boron_carbide1
from materials import cr_shell2, boron_carbide2, cr_shell3, boron_carbide3, cr_shell4, boron_carbide4, cr_shell5, boron_carbide5
from materials import grey_rods
from fuel_assemblies import find_name, fa_types

def write_floats_to_file(filename, float_array, elements_per_line):
    with open(filename, 'w') as file:
        for i in range(0, len(float_array), elements_per_line):
            line = float_array[i:i + elements_per_line]  # Получаем срез массива
            file.write('\t'.join(f'{num:.8f}' for num in line) + '\n')

if __name__ == '__main__':
    grey_array = [element for sublist in grey_rods for element in sublist]
    mats = openmc.Materials((*g_hole, *fuel, *gaz, *shell, *coolant,
                              *cr_shell1, *cr_shell2, *cr_shell3, *cr_shell4, *cr_shell5,
                             *boron_carbide1, *boron_carbide2, *boron_carbide3, *boron_carbide4, *boron_carbide5, *grey_array))
    mats.export_to_xml()

    fa_universe = []
    splits = []
    dif_fa_universe = []
    flag = False
    mat_num = 0
    fake_grey = []
    for i in range(0, split_number):
        fake_grey.append(None)
    for i in range(0, len(dif_fu_cart)):
        dict_fa = find_name(dif_fu_cart[i], fa_types)
        if dict_fa["grey_enrichment"] == 0:
            mat_grey = fake_grey
        else:
            mat_grey = grey_rods[mat_num]
            mat_num += 1
        for j in range(0, len(g1)):
            if i == g1[j]-1:
                fa_ = full_fa(i, g_hole, fuel, gaz, shell, coolant, boron_carbide1, cr_shell1, h1, mat_grey, dict_fa["grey_pos"])
                flag = True
                dif_fa_universe.append(fa_)
                splits += list(fa_.cells.values())
                break
        if flag:
            flag = False
            continue
        else:
            for j in range(0, len(g2)):
                if i == g2[j]-1:
                    fa_ = full_fa(i, g_hole, fuel, gaz, shell, coolant, boron_carbide2, cr_shell2, h2, mat_grey, dict_fa["grey_pos"])
                    flag = True
                    dif_fa_universe.append(fa_)
                    splits += list(fa_.cells.values())
                    break
            if flag:
                flag = False
                continue
            else:
                for j in range(0, len(g3)):
                    if i == g3[j]-1:
                        fa_ = full_fa(i, g_hole, fuel, gaz, shell, coolant, boron_carbide3, cr_shell3, h3, mat_grey, dict_fa["grey_pos"])
                        flag = True
                        dif_fa_universe.append(fa_)
                        splits += list(fa_.cells.values())
                        break
                if flag:
                    flag = False
                    continue
                else:
                    for j in range(0, len(g4)):
                        if i == g4[j]-1:
                            fa_ = full_fa(i, g_hole, fuel, gaz, shell, coolant, boron_carbide4, cr_shell4, h4, mat_grey, dict_fa["grey_pos"])
                            flag = True
                            dif_fa_universe.append(fa_)
                            splits += list(fa_.cells.values())
                            break
                    if flag:
                        flag = False
                        continue
                    else:
                        for j in range(0, len(g5)):
                            if i == g5[j]-1:
                                fa_ = full_fa(i, g_hole, fuel, gaz, shell, coolant, boron_carbide5, cr_shell5, h5, mat_grey, dict_fa["grey_pos"])
                                flag = True
                                dif_fa_universe.append(fa_)
                                splits += list(fa_.cells.values())
                                break
                        if flag:
                            flag = False
                            continue
                        else:
                            fa_ = full_fa(i, g_hole, fuel, gaz, shell, coolant, boron_carbide5, cr_shell5, 0, mat_grey, dict_fa["grey_pos"])
                            dif_fa_universe.append(fa_)
                            splits += list(fa_.cells.values())

    for i in range(0, n_fa):
        fa_universe.append(dif_fa_universe[numbers[i]-1])

    print("The core is divided into", len(splits), "elements.")
    fa_universe.append(water_full_fa(coolant[-1]))

    lat = openmc.HexLattice()
    lat.orientation = 'y'
    lat.center = (0., 0.)
    lat.pitch = (turnkey_size,)
    outer_coolant_cell = openmc.Cell(fill=coolant[-1])
    outer_universe = openmc.Universe(cells=(outer_coolant_cell,))
    lat.outer = outer_universe

    ring8 = [fa_universe[n_fa]]+[fa_universe[n_fa]]+[fa_universe[56]]+[fa_universe[43]]+[fa_universe[31]]+[fa_universe[20]]+[fa_universe[n_fa]]
    ring8 += [fa_universe[n_fa]]+[fa_universe[n_fa]]+[fa_universe[3]]+[fa_universe[2]]+[fa_universe[1]]+[fa_universe[0]] +[fa_universe[n_fa]]
    ring8 += [fa_universe[n_fa]]+[fa_universe[n_fa]]+[fa_universe[11]]+[fa_universe[21]]+[fa_universe[32]]+[fa_universe[44]] +[fa_universe[n_fa]]
    ring8 += [fa_universe[n_fa]]+[fa_universe[n_fa]]+[fa_universe[94]]+[fa_universe[107]]+[fa_universe[119]]+[fa_universe[130]] +[fa_universe[n_fa]]
    ring8 += [fa_universe[n_fa]]+[fa_universe[n_fa]]+[fa_universe[147]]+[fa_universe[148]]+[fa_universe[149]]+[fa_universe[150]] +[fa_universe[n_fa]]
    ring8 += [fa_universe[n_fa]]+[fa_universe[n_fa]]+[fa_universe[139]]+[fa_universe[129]]+[fa_universe[118]]+[fa_universe[106]] +[fa_universe[n_fa]]

    ring7 = [fa_universe[81]]+[fa_universe[68]]+[fa_universe[55]]+[fa_universe[42]]+[fa_universe[30]]+[fa_universe[19]]
    ring7 += [fa_universe[10]]+[fa_universe[9]]+[fa_universe[8]]+[fa_universe[7]]+[fa_universe[6]]+[fa_universe[5]]
    ring7 += [fa_universe[4]]+[fa_universe[12]]+[fa_universe[22]]+[fa_universe[33]]+[fa_universe[45]]+[fa_universe[57]]
    ring7 += [fa_universe[69]]+[fa_universe[82]]+[fa_universe[95]]+[fa_universe[108]]+[fa_universe[120]]+[fa_universe[131]]
    ring7 += [fa_universe[140]]+[fa_universe[141]]+[fa_universe[142]]+[fa_universe[143]]+[fa_universe[144]]+[fa_universe[145]]
    ring7 += [fa_universe[146]]+[fa_universe[138]]+[fa_universe[128]]+[fa_universe[117]]+[fa_universe[105]]+[fa_universe[93]]

    ring6 = [fa_universe[80]]+[fa_universe[67]]+[fa_universe[54]]+[fa_universe[41]]+[fa_universe[29]]
    ring6 += [fa_universe[18]]+[fa_universe[17]]+[fa_universe[16]]+[fa_universe[15]]+[fa_universe[14]]
    ring6 += [fa_universe[13]]+[fa_universe[23]]+[fa_universe[34]]+[fa_universe[46]]+[fa_universe[58]]
    ring6 += [fa_universe[70]]+[fa_universe[83]]+[fa_universe[96]]+[fa_universe[109]]+[fa_universe[121]]
    ring6 += [fa_universe[132]]+[fa_universe[133]]+[fa_universe[134]]+[fa_universe[135]]+[fa_universe[136]]
    ring6 += [fa_universe[137]]+[fa_universe[127]]+[fa_universe[116]]+[fa_universe[104]]+[fa_universe[92]]

    ring5 = [fa_universe[79]]+[fa_universe[66]]+[fa_universe[53]]+[fa_universe[40]]
    ring5 += [fa_universe[28]]+[fa_universe[27]]+[fa_universe[26]]+[fa_universe[25]]
    ring5 += [fa_universe[24]]+[fa_universe[35]]+[fa_universe[47]]+[fa_universe[59]]
    ring5 += [fa_universe[71]]+[fa_universe[84]]+[fa_universe[97]]+[fa_universe[110]]
    ring5 += [fa_universe[122]]+[fa_universe[123]]+[fa_universe[124]]+[fa_universe[125]]
    ring5 += [fa_universe[126]]+[fa_universe[115]]+[fa_universe[103]]+[fa_universe[91]]

    ring4 = [fa_universe[78]]+[fa_universe[65]]+[fa_universe[52]]
    ring4 += [fa_universe[39]]+[fa_universe[38]]+[fa_universe[37]]
    ring4 += [fa_universe[36]]+[fa_universe[48]]+[fa_universe[60]]
    ring4 += [fa_universe[72]]+[fa_universe[85]]+[fa_universe[98]]
    ring4 += [fa_universe[111]]+[fa_universe[112]]+[fa_universe[113]]
    ring4 += [fa_universe[114]]+[fa_universe[102]]+[fa_universe[90]]

    ring3 = [fa_universe[77]]+[fa_universe[64]]
    ring3 += [fa_universe[51]]+[fa_universe[50]]
    ring3 += [fa_universe[49]]+[fa_universe[61]]
    ring3 += [fa_universe[73]]+[fa_universe[86]]
    ring3 += [fa_universe[99]]+[fa_universe[100]]
    ring3 += [fa_universe[101]]+[fa_universe[89]]

    ring2 = [fa_universe[76]]+[fa_universe[63]]+[fa_universe[62]]+[fa_universe[74]]+[fa_universe[87]]+[fa_universe[88]]

    ring1 = [fa_universe[75]]

    lat.universes = [ring8, ring7, ring6, ring5, ring4, ring3, ring2, ring1]

    in_core_barrel = openmc.ZCylinder(surface_id=int(2E7 + 7E5 + 200*1E2 + 98), r=core_barrel_in_r)
    zdn_ = openmc.ZPlane(surface_id=int(2E7 + 9E5 + 200*1E2 + 97), z0= 0, boundary_type = 'vacuum')
    zup_ = openmc.ZPlane(surface_id=int(2E7 + 10E5 + 200*1E2 + 96), z0= core_height*1E2, boundary_type = 'vacuum')
    core_cell = openmc.Cell(cell_id = int(3E7 + 13E5 + 200*1E2 + 98), name="in_core_barrel", fill = lat, region = -in_core_barrel & -zup_ & +zdn_)


    out_core_barrel = openmc.ZCylinder(surface_id=int(2E7 + 7E5 + 200*1E2 + 94), r=core_barrel_out_r, boundary_type = "vacuum")
    barrel_cell = openmc.Cell(cell_id = int(3E7 + 13E5 + 200*1E2 + 97), name="core_barrel", fill = shell[len(shell)-1], region = -out_core_barrel & +in_core_barrel & -zup_ & +zdn_)
    reactor_universe = openmc.Universe(universe_id = int(4E7 + 7E5 + 201*1E2 + 60), name = "reactor", cells = [core_cell, barrel_cell])
    geom = openmc.Geometry(reactor_universe)
    geom.export_to_xml()


    p = openmc.Plot()

    #cross section drawing
    p.origin=(0, 0, core_height*1E2 - 5)
    p.filename = 'cluster_xy'
    p.basis = "xy"
    p.width = (350, 350)
    p.pixels = (6000, 6000)
    p.color_by = 'material'
    plots = openmc.Plots([p])
    plots.export_to_xml()
    openmc.plot_geometry()
    #------------
    # longitudinal section drawing
    length = 350
    p.origin=(0,0,length/2)
    p.filename = 'cluster_yz'
    p.basis = "yz"
    p.width = (500, 500)
    p.pixels = (5000, 5000)
    p.color_by = 'material'
    plots = openmc.Plots([p])
    plots.export_to_xml()
    openmc.plot_geometry()
    #-------------------



	#Computing settings
    set = openmc.Settings()
    set.batches = batches
    set.inactive = inactive
    set.particles = particles
    set.temperature = {'method': 'interpolation'}
    sourse_point = openmc.stats.Point(xyz=(0,0, core_height*1E2/2))
    set.source =openmc.Source(space = sourse_point)
    set.export_to_xml()

    tallies = openmc.Tallies()

    # Instantiate flux Tally in moderator and fuel
    tally = openmc.Tally(name='energy_rel')
    energy_filter = openmc.EnergyFilter([0., 20.0e6])
    tally.filters = [openmc.CellFilter(splits)]

    tally.filters.append(energy_filter)
    tally.scores = ['fission']
    tallies.append(tally)
    tallies.export_to_xml()

    openmc.run()

    sp = openmc.StatePoint(f'statepoint.{batches}.h5')#на каждый шаг расчета сделан файл с результатами
    energy_rel = sp.get_tally(name='energy_rel') #вернули энерговыделение

    #Get a pandas dataframe for the mesh tally data
    df = energy_rel.get_pandas_dataframe()#обработали текст

    # save results
    with open("results.txt", 'w') as f:
        f.write(df.to_string())
    values=energy_rel.get_values()
    values2=np.array(values.flatten())
    values2_stdev=np.array((energy_rel.get_values(value = 'std_dev')).flatten()) #stdev
    count = Counter(numbers)
    for i in range(0, len(dif_fu_cart)):
        divider = count[i+1]
        for j in range(i * split_number, (i + 1) * split_number):
            values2[j] /= divider
            values2_stdev[j] /= divider

    meanval=sum(values2)/len(values2)
    meanval_stdev = sqrt(sum(stddev**2 for stddev in values2_stdev))/len(values2_stdev)
    kv_=values2/meanval
    kv__stdev = []
    for i in range(0 , len(kv_)):
        kv__stdev.append(sqrt((values2_stdev[i]/meanval)**2 + (values2[i]/(meanval**2)*meanval_stdev)**2))

    kv_stdev = []
    kv = []
    q_element = []
    for i in range(0, n_fa):
        for j in range(0, split_number):
            kv.append(kv_[(numbers[i] - 1) * split_number + j])
            q_element.append(kv_[(numbers[i] - 1) * split_number + j] * (q_r / (n_fa * split_number)))
            kv_stdev.append(kv__stdev[(numbers[i] - 1) * split_number + j])
    write_floats_to_file("kv_stdev.txt", kv_stdev, split_number)
    write_floats_to_file("kv.txt", kv, split_number)
    write_floats_to_file("Q6.txt", q_element, split_number)

    #Kq calculation
    kq = []
    for i in range(0, n_fa):
        sum_ = 0
        for j in range(0, split_number):
            sum_ += kv[j+i*split_number]
        kq.append(sum_/split_number)
    xy =[]
    turnkey_size = turnkey_size/100
    y = (-(len(line)//2)//2 - len(line)//2 + 0.5)*turnkey_size/sqrt(3)
    for i in range(0,len(line)):
        x = -(line[i]/2 + 0.5)*turnkey_size
        for j in range(0,line[i]):
            x += turnkey_size
            xy.append((x,y))
        y += 1.5*turnkey_size/sqrt(3)
    combined_array = [(t[0], t[1], n) for t, n in zip(xy, kq)]
    np.savetxt(f"kq.txt", combined_array, delimiter="\t", fmt = "%.6f")

    #Kz and Kv in central fuel assembly calculation
    kz = []
    kv_central = []
    increment_z = core_height/split_number
    z = -core_height/2 + increment_z/2
    for i in range(0, split_number):
        sum_ = 0
        for j in range(0, n_fa):
            sum_ += kv[i+j*split_number]
        kz.append((z,sum_/n_fa))
        kv_central.append((z, kv_[(len(dif_fu_cart) - 1) * split_number + i], kv__stdev[(len(dif_fu_cart) - 1) * split_number + i]))
        z += increment_z
    np.savetxt(f"kz.txt", kz, delimiter="\t", fmt = "%.6f")
    np.savetxt(f"kv_central.txt", kv_central, delimiter="\t", fmt = "%.8f")
