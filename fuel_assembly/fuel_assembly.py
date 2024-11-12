import openmc
import openmc.deplete
import numpy as np
import math

import sys
sys.path.append('../')

from constants import split_number, core_height
sys.path.append('../'+'materials')
from materials import g_hole, fuel, gaz, shell, coolant
from assembly_element import fa_split, water_fa_split, turnkey_size

def full_fa(fa_num, c_gaz_list, fuel_list, gaz_list, shell_list, hc_list):
    fa_universe_return = openmc.Universe(universe_id = int(4E7 + 5E5 + fa_num*1E2 + 31), name=f'fa_universe_{fa_num}')
    for i in range (0, split_number):
        num = fa_num*split_number + i
        root_cell = fa_split(fa_num, i, c_gaz_list[num], fuel_list[num], gaz_list[num], shell_list[num], hc_list[num])
        fa_universe_return.add_cell(root_cell)

    return fa_universe_return

def water_full_fa(fa_num, hc_list):
    fa_universe_return = openmc.Universe(universe_id = int(4E7 + 6E5 + fa_num*1E2 + 32), name=f'fa_universe_{fa_num}')
    for i in range (0, split_number):
        num = fa_num*split_number + i
        root_cell = water_fa_split(fa_num, i, hc_list[num])
        fa_universe_return.add_cell(root_cell)

    return fa_universe_return

if __name__ == '__main__':

    mats = openmc.Materials((*g_hole, *fuel, *gaz, *shell, *coolant))
    mats.export_to_xml()
    zdn_ = openmc.ZPlane(surface_id=int(2E7 + 9E5 + 200*1E2 + 31), z0= 0, boundary_type = 'vacuum')
    zup_ = openmc.ZPlane(surface_id=int(2E7 + 10E5 + 200*1E2 + 32), z0= core_height*1E2, boundary_type = 'vacuum')
    assembly_prism = openmc.model.HexagonalPrism(edge_length = turnkey_size/math.sqrt(3), orientation = 'x', boundary_type = "reflective")
    fa_region = -assembly_prism & -zup_ & +zdn_
    full_fa_u = full_fa(1, g_hole, fuel, gaz, shell, coolant)
    fa_cell = openmc.Cell(cell_id =  int(3E7 + 13E5 + 200*1E2 + 31), name = "fa1", fill = full_fa_u)
    fa_cell.region = fa_region
    fa_universe = openmc.Universe(universe_id=int(4E7 + 7E5 + 201*1E2 + 33), cells=[fa_cell])
    geom = openmc.Geometry(fa_universe)
    geom.export_to_xml()


    z0 = 0
    z_us = 400


    p = openmc.Plot()

    #cross section drawing
    p.origin=(0,0,20)
    p.filename = 'fuel_assembly_xy'
    p.basis = "xy"
    p.width = (28, 28)
    p.pixels = (1000, 1000)
    p.color_by = 'material'
    plots = openmc.Plots([p])
    plots.export_to_xml()
    openmc.plot_geometry()
    #-----------------------

    # longitudinal section drawing
    p.origin=(0,0,core_height/2*1E2)
    p.filename = 'fuel_assembly_yz'
    p.basis = "xz"
    p.width = (360, 360)
    p.pixels = (2000, 2000)
    p.color_by = 'material'
    plots = openmc.Plots([p])
    plots.export_to_xml()
    openmc.plot_geometry()
    #-------------------


    #Computing settings
    batches = 2000
    inactive = 10
    particles = 10000

    settings = openmc.Settings()
    settings.batches = batches
    settings.inactive = inactive
    settings.particles = particles
    settings.output = {'tallies': True}
    settings.temperature = {'method': 'interpolation'}
    # instalation sorse
    sourse_point = openmc.stats.Point(xyz=(0,0, core_height*1E2/2))
    settings.source =openmc.Source(space = sourse_point)
    settings.export_to_xml()

    tallies = openmc.Tallies()
    splits = list(full_fa_u.cells.values())
    print(splits)

    # Instantiate flux Tally in moderator and fuel
    tally = openmc.Tally(name='energovidilenie')
    energy_filter = openmc.EnergyFilter([0., 20.0e6])
    #tally.filters = [ openmc.MaterialFilter(  fuel_array + gadalunm_array ) ]
    tally.filters = [openmc.CellFilter(splits)]

    tally.filters.append(energy_filter)
    tally.scores = ['fission']
    tallies.append(tally)
    tallies.export_to_xml()

    openmc.run()

    sp = openmc.StatePoint('statepoint.200.h5')#на каждый шаг расчета сделан файл с результатами
    energovidilenie = sp.get_tally(name='energovidilenie') #вернули энерговыделение

    #Get a pandas dataframe for the mesh tally data
    df = energovidilenie.get_pandas_dataframe()#обработали текст

    # сохраняем изначальные результаты в файл
    with open("results.txt", 'w') as f:
        f.write(df.to_string())
    values=energovidilenie.get_values()
    values2=np.array(values.flatten())
    meanval=sum(values2)/len(values2)
    kq=values2/meanval
    res = []
    z_incriment = core_height/split_number
    z_ = z_incriment/2
    for i in range(0, split_number):
        res.append((z_+z_incriment*i, kq[i]))
    np.savetxt(f"kq.csv", res, delimiter="\t", fmt = "%.6f")
