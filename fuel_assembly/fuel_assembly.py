import openmc
import openmc.deplete
import numpy as np
import math

import sys
sys.path.append('../')
from constants import split_number, core_height, n_fa, batches, particles, inactive, turnkey_size, h1
sys.path.append('../'+'materials')
from materials import g_hole, fuel, gaz, shell, coolant, cr_shell1, boron_carbide1
from assembly_element import fa_split

def full_fa(fa_num, c_gaz_list, fuel_list, gaz_list, shell_list, hc_list, b4c_list, cr_steel_list, cr_depth):
    fa_universe_return = openmc.Universe(universe_id = int(4E7 + 5E5 + fa_num*1E2 + split_number), name=f'fa_universe_{fa_num}')
    for i in range (0, split_number - cr_depth):
        num = fa_num*split_number + i
        root_cell = fa_split(fa_num, i, c_gaz_list[num], fuel_list[num], gaz_list[num], shell_list[num], hc_list[num], 0, 0, 0)
        fa_universe_return.add_cell(root_cell)
    for i in range (split_number - cr_depth, split_number):
        num = fa_num*split_number + i
        cr_num = i - (split_number - cr_depth)
        root_cell = fa_split(fa_num, i, c_gaz_list[num], fuel_list[num], gaz_list[num], shell_list[num], hc_list[num], b4c_list[cr_num], cr_steel_list[cr_num], 1)
        fa_universe_return.add_cell(root_cell)
    return fa_universe_return

def water_full_fa(coolant):
    fa_universe_return = openmc.Universe(universe_id = int(4E7 + 6E5 + (n_fa//6 + 1)*1E2 + split_number + 1), name=f'fa_universe_{n_fa}')
    zdn = openmc.ZPlane(surface_id=int(2E7 + 9E5 + (n_fa//6 + 1)*1E2 + split_number), z0= 0)
    zup = openmc.ZPlane(surface_id=int(2E7 + 10E5 + (n_fa//6 + 1)*1E2 + split_number+1), z0= core_height*1E2)
    zdn.boundary_type="vacuum"
    zup.boundary_type="vacuum"
    water_assembly_cell = openmc.Cell(cell_id = int(3E7 + 14E5 + (n_fa//6 + 1)*1E2 + split_number), name='water_cell_assembly')
    hex_prizm = openmc.model.HexagonalPrism(edge_length = turnkey_size/math.sqrt(3), orientation = 'x', boundary_type = "transmission")
    water_assembly_cell.region = -hex_prizm & -zup & +zdn
    water_assembly_cell.fill = coolant
    fa_universe_return.add_cell(water_assembly_cell)
    return fa_universe_return

if __name__ == '__main__':

    mats = openmc.Materials((*g_hole, *fuel, *gaz, *shell, *coolant, *cr_shell1, *boron_carbide1))
    mats.export_to_xml()
    zdn_ = openmc.ZPlane(surface_id=int(2E7 + 9E5 + 200*1E2 + 31), z0= 0, boundary_type = 'vacuum')
    zup_ = openmc.ZPlane(surface_id=int(2E7 + 10E5 + 200*1E2 + 32), z0= core_height*1E2, boundary_type = 'vacuum')
    assembly_prism = openmc.model.HexagonalPrism(edge_length = turnkey_size/math.sqrt(3), orientation = 'x', boundary_type = "reflective")
    fa_region = -assembly_prism & -zup_ & +zdn_
    full_fa_u = full_fa(1, g_hole, fuel, gaz, shell, coolant, boron_carbide1, cr_shell1, h1 )
    fa_cell = openmc.Cell(cell_id =  int(3E7 + 13E5 + 200*1E2 + 31), name = "fa1", fill = full_fa_u)
    fa_cell.region = fa_region
    fa_universe = openmc.Universe(universe_id=int(4E7 + 7E5 + 201*1E2 + 33), cells=[fa_cell])
    geom = openmc.Geometry(fa_universe)
    geom.export_to_xml()


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

    # Instantiate flux Tally in moderator and fuel
    tally = openmc.Tally(name='energovidilenie')
    energy_filter = openmc.EnergyFilter([0., 20.0e6])
    tally.filters = [openmc.CellFilter(splits)]

    tally.filters.append(energy_filter)
    tally.scores = ['fission']
    tallies.append(tally)
    tallies.export_to_xml()

    openmc.run()

    sp = openmc.StatePoint(f'statepoint.{batches}.h5')#на каждый шаг расчета сделан файл с результатами
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
    z_ = z_incriment/2-core_height/2
    for i in range(0, split_number):
        res.append((z_+z_incriment*i, kq[i]))
    np.savetxt(f"kq.csv", res, delimiter="\t", fmt = "%.6f")
