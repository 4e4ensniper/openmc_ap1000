import openmc
import math
import sys
sys.path.append('../')
from constants import split_number, core_height, rod_pitch, turnkey_size

#regions ID definition
#2??????? - first 2
#2__????? - region number
#2??___?? - fuel asssembly number
#2?????__ - split number from bottom

def fa_split(fa_num, layer_num, c_gaz, fuel, gaz, shell, coolant, b4c, cr_steel, cr_key, g_fuel, grey_list):

    central_hole_cylinder = openmc.ZCylinder(surface_id=int(2E7 + 1E5 + fa_num*1E2 + layer_num), r=0.117)
    fuel_cylinder = openmc.ZCylinder(surface_id=int(2E7 + 2E5 + fa_num*1E2 + layer_num), r=0.378)
    gap_cylinder = openmc.ZCylinder(surface_id=int(2E7 + 3E5 + fa_num*1E2 + layer_num), r=0.386)
    clad_cylinder = openmc.ZCylinder(surface_id=int(2E7 + 4E5 + fa_num*1E2 + layer_num), r=0.455)
    #central tube
    central_tube_out_cylinder = openmc.ZCylinder(surface_id=int(2E7 + 5E5 + fa_num*1E2 + layer_num), r=0.65)
    central_tube_in_cylinder = openmc.ZCylinder(surface_id=int(2E7 + 6E5 + fa_num*1E2 + layer_num), r=0.55)
    #control sistem channel
    csc_out_cylinder = openmc.ZCylinder(surface_id=int(2E7 + 7E5 + fa_num*1E2 + layer_num), r=0.63)
    csc_in_cylinder = openmc.ZCylinder(surface_id=int(2E7 + 8E5 + fa_num*1E2 + layer_num), r=0.545)

    d_z = core_height*1E2/split_number
    #regions
    zdn = openmc.ZPlane(surface_id=int(2E7 + 9E5 + fa_num*1E2 + layer_num), z0= layer_num * d_z)
    zup = openmc.ZPlane(surface_id=int(2E7 + 10E5 + fa_num*1E2 + layer_num), z0= (layer_num + 1) * d_z)
    if layer_num == 0:
        zdn.boundary_type="vacuum"
        zup.boundary_type="transmission"
    elif layer_num == split_number - 1:
        zdn.boundary_type="transmission"
        zup.boundary_type="vacuum"
    else:
        zdn.boundary_type="transmission"
        zup.boundary_type="transmission"
    #fuel rod
    central_hole_region = -central_hole_cylinder
    fuel_region = +central_hole_cylinder & -fuel_cylinder
    gap_region = +fuel_cylinder & -gap_cylinder
    clad_region = +gap_cylinder & -clad_cylinder
    #central tube
    central_tube_region = -central_tube_out_cylinder & +central_tube_in_cylinder
    #control sistem channel
    csc_region = -csc_out_cylinder & +csc_in_cylinder

    #cell ID definition
    #3??????? - first 3
    #3__????? - cell number
    #3??___?? - fuel asssembly number
    #3?????__ - split number from bottom

    #universe ID definition
    #4??????? - first 3
    #4__????? - universe number
    #4??___?? - fuel asssembly number
    #4?????__ - split number from bottom

    #create fuel rod
    central_hole_cell = openmc.Cell(cell_id = int(3E7 + 1E5 + fa_num*1E2 + layer_num), fill = c_gaz, region = central_hole_region)
    fuel_cell = openmc.Cell(cell_id = int(3E7 + 2E5 + fa_num*1E2 + layer_num), fill = fuel, region = fuel_region)
    gap_cell = openmc.Cell(cell_id = int(3E7 + 3E5 + fa_num*1E2 + layer_num), fill = gaz, region = gap_region)
    clad_cell = openmc.Cell(cell_id = int(3E7 + 4E5 + fa_num*1E2 + layer_num), fill = shell, region = clad_region)
    water1_cell = openmc.Cell(cell_id = int(3E7 + 5E5 + fa_num*1E2 + layer_num), fill = coolant, region = +clad_cylinder)
    fr = openmc.Universe(universe_id = int(4E7 + 1E5 + fa_num*1E2 + layer_num), cells=[central_hole_cell, fuel_cell, gap_cell, clad_cell, water1_cell])

    #create grey rod
    central_hole_cell_g = openmc.Cell(cell_id = int(3E7 + 6E5 + fa_num*1E2 + layer_num), fill = c_gaz, region = central_hole_region)
    gap_cell_g = openmc.Cell(cell_id = int(3E7 + 7E5 + fa_num*1E2 + layer_num), fill = gaz, region = gap_region)
    clad_cell_g = openmc.Cell(cell_id = int(3E7 + 8E5 + fa_num*1E2 + layer_num), fill = shell, region = clad_region)
    water1_cell_g = openmc.Cell(cell_id = int(3E7 + 9E5 + fa_num*1E2 + layer_num), fill = coolant, region = +clad_cylinder)
    grey_cell = openmc.Cell(cell_id = int(3E7 + 10E5 + fa_num*1E2 + layer_num), fill = g_fuel, region = fuel_region)
    grey_f = openmc.Universe(universe_id = int(4E7 + 2E5 + fa_num*1E2 + layer_num), cells=[central_hole_cell_g, grey_cell, gap_cell_g, clad_cell_g, water1_cell_g])
    #create central tube
    central_tube_in_cell = openmc.Cell(cell_id = int(3E7 + 11E5 + fa_num*1E2 + layer_num), fill = coolant, region = -central_tube_in_cylinder)
    central_tube_cell = openmc.Cell(cell_id = int(3E7 + 12E5 + fa_num*1E2 + layer_num), fill = shell, region = central_tube_region)
    central_tube_cell.temperature = coolant.temperature
    water2_cell = openmc.Cell(cell_id = int(3E7 + 13E5 + fa_num*1E2 + layer_num), fill = coolant, region = +central_tube_out_cylinder)
    ct = openmc.Universe(universe_id = int(4E7 + 3E5 + fa_num*1E2 + layer_num), cells=[central_tube_in_cell, central_tube_cell, water2_cell])

    #create control sistem channel
    csc_cell = openmc.Cell(cell_id = int(3E7 + 14E5 + fa_num*1E2 + layer_num),fill = shell, region = csc_region)
    csc_cell .temperature = coolant.temperature
    water3_cell = openmc.Cell(cell_id = int(3E7 + 15E5 + fa_num*1E2 + layer_num), fill = coolant, region = +csc_out_cylinder)
    if cr_key == 0:
        csc_in_cell = openmc.Cell(cell_id = int(3E7 + 16E5 + fa_num*1E2 + layer_num), fill = coolant, region = -csc_in_cylinder)
        csc = openmc.Universe(universe_id = int(4E7 + 4E5 + fa_num*1E2 + layer_num), cells = [csc_in_cell, csc_cell, water3_cell])
    else:
        cr_in_cylinder = openmc.ZCylinder(surface_id=int(2E7 + 11E5 + fa_num*1E2 + layer_num), r=0.35)
        cr_out_cylinder = openmc.ZCylinder(surface_id=int(2E7 + 12E5 + fa_num*1E2 + layer_num), r=0.41)
        cr_shell = +cr_in_cylinder & -cr_out_cylinder
        cr_coolant = +cr_out_cylinder & -csc_in_cylinder
        cr_b4c_cell = openmc.Cell(cell_id = int(3E7 + 17E5 + fa_num*1E2 + layer_num), fill = b4c, region = -cr_in_cylinder)
        cr_shell_cell = openmc.Cell(cell_id = int(3E7 + 18E5 + fa_num*1E2 + layer_num), fill = cr_steel, region = cr_shell)
        cr_coolant_cell = openmc.Cell(cell_id = int(3E7 + 19E5 + fa_num*1E2 + layer_num), fill = coolant, region = cr_coolant)
        csc = openmc.Universe(universe_id = int(4E7 + 4E5 + fa_num*1E2 + layer_num), cells = [cr_b4c_cell, cr_shell_cell, cr_coolant_cell, csc_cell, water3_cell])



    #coolant universe
    all_water_out=openmc.Cell(cell_id = int(3E7 + 20E5 + fa_num*1E2 + layer_num), fill=coolant)
    all_water_out_u=openmc.Universe(universe_id = int(4E7 + 5E5 + fa_num*1E2 + layer_num), cells=[all_water_out])

    #lattice ID definition
    #5??????? - first 3
    #5__????? - lattice number
    #5??___?? - fuel asssembly number
    #5?????__ - split number from bottom

    #hex lattice of fuel assembly
    lat=openmc.HexLattice(lattice_id = int(5E7 + 1E5 + fa_num*1E2 + layer_num), name=f'assembly_{fa_num}_layer_{layer_num}')
    lat.center = (0.0, 0.0)
    lat.pitch = (rod_pitch,)
    lat.outer=all_water_out_u
    lat.orientation = 'x'
    ring0 = [ct]
    ring1 = [fr]*6
    ring2 = [fr]*12
    ring3 = [fr] + [csc, fr, fr]*5 + [csc, fr]
    ring4 = [fr]*24
    ring5 = ([csc] + [fr]*4)*6
    ring6 = ([fr]*3 + [csc] + [fr]*2)*6
    ring7 = [fr]*7*6
    ring8 = [fr]*8*6
    ring9 = [fr]*9*6
    ring10 = [fr]*10*6
    rings = [ring10, ring9, ring8, ring7, ring6, ring5, ring4, ring3, ring2, ring1, ring0]
    count = 0
    num = 0
    for i in range (0, len(rings)):
        for j in range(0, len(rings[i])):
            if num < len(grey_list):
                if count == grey_list[num]:
                    rings[i][j] = grey_f
                    num += 1
                count +=1
            else:
                break
    lat.universes = rings


    assembly_cell = openmc.Cell(cell_id = int(3E7 + 21E5 + fa_num*1E2 + layer_num), name=f'cell_assembly_{fa_num}_layer_{layer_num}' )
    hex_prizm = openmc.model.HexagonalPrism(edge_length = turnkey_size/math.sqrt(3), orientation = 'x', boundary_type = "transmission")
    assembly_cell.region = -hex_prizm & -zup & +zdn
    assembly_cell.fill = lat

    return assembly_cell
