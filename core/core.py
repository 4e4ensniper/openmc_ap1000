import openmc
import openmc.deplete
import sys

sys.path.append('../')
from constants import n_fa, turnkey_size, core_barrel_in_r, core_barrel_out_r, core_height, split_number

sys.path.append('../'+'fuel_assembly')
from assembly_element import fa_split, water_fa_split
from fuel_assembly import full_fa, water_full_fa

sys.path.append('../'+'materials')
from materials import g_hole, fuel, gaz, shell, coolant


if __name__ == '__main__':

    mats = openmc.Materials((*g_hole, *fuel, *gaz, *shell, *coolant))
    mats.export_to_xml()

    fa_universe = []
    for i in range(0, n_fa):
        fa_universe.append(full_fa(i, g_hole, fuel, gaz, shell, coolant))
    for i in range(n_fa, n_fa + 18):
        fa_universe.append(water_full_fa(i, coolant))

    lat = openmc.HexLattice()
    lat.orientation = 'y'
    lat.center = (0., 0.)
    lat.pitch = (turnkey_size,)
    outer_coolant_cell = openmc.Cell(fill=coolant[-1])
    outer_universe = openmc.Universe(cells=(outer_coolant_cell,))
    lat.outer = outer_universe

    ring8 = [fa_universe[n_fa]]+[fa_universe[n_fa+1]]+[fa_universe[56]]+[fa_universe[43]]+[fa_universe[31]]+[fa_universe[20]]+[fa_universe[n_fa+2]]
    ring8 += [fa_universe[n_fa+3]]+[fa_universe[n_fa+4]]+[fa_universe[3]]+[fa_universe[2]]+[fa_universe[1]]+[fa_universe[0]] +[fa_universe[n_fa+5]]
    ring8 += [fa_universe[n_fa+6]]+[fa_universe[n_fa+7]]+[fa_universe[11]]+[fa_universe[21]]+[fa_universe[32]]+[fa_universe[44]] +[fa_universe[n_fa+8]]
    ring8 += [fa_universe[n_fa+9]]+[fa_universe[n_fa+10]]+[fa_universe[94]]+[fa_universe[107]]+[fa_universe[119]]+[fa_universe[130]] +[fa_universe[n_fa+11]]
    ring8 += [fa_universe[n_fa+12]]+[fa_universe[n_fa+13]]+[fa_universe[147]]+[fa_universe[148]]+[fa_universe[149]]+[fa_universe[150]] +[fa_universe[n_fa+14]]
    ring8 += [fa_universe[n_fa+15]]+[fa_universe[n_fa+16]]+[fa_universe[139]]+[fa_universe[129]]+[fa_universe[118]]+[fa_universe[106]] +[fa_universe[n_fa+17]]

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

    z0 = 0
    z_us = 400

    p = openmc.Plot()

    #cross section drawing
    p.origin=(0,0,100)
    p.filename = 'cluster_xy'
    p.basis = "xy"
    p.width = (350, 350)
    p.pixels = (2000, 2000)
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
    p.pixels = (1000, 1000)
    p.color_by = 'material'
    plots = openmc.Plots([p])
    plots.export_to_xml()
    openmc.plot_geometry()
    #-------------------



	#Computing settings
    batches = 20
    inactive = 10
    particles = 1000

    set = openmc.Settings()
    set.batches = batches
    set.inactive = inactive
    set.particles = particles
    set.temperature = {'method': 'interpolation'}
    sourse_point = openmc.stats.Point(xyz=(0,0, core_height*1E2/2))
    set.source =openmc.Source(space = sourse_point)
    set.export_to_xml()

    openmc.run()#, geometry_debug=True)
