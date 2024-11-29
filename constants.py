split_number = 30
core_height = 3.53 #m
n_fa = 151
r_fuel = 7.57/20 #sm
r_hole = 2.35/20 #sm
fr_number = 312
rod_pitch = 1.275 #sm
turnkey_size = 23.4 #sm
core_barrel_in_r = 339.72/2 #sm
core_barrel_out_r = 349.88/2 #sm
b_conc = 4 #g/kg
#The number of fuel assemblies per row, starting from the bottom of the core map.
line = [4, 7, 10, 11, 12, 13, 12, 13, 12, 13, 12, 11, 10, 7, 4]
#Path to the location of files with temperature distributions.
csv_path = "/home/ubuntu24/Desktop/openmc_ap1000/materials/temperature_distributions/"
#Calculation parameters
batches = 200
inactive = 10
particles = 10000
#Number of different fuel assemblies defined as different universes
dif_fu_cart = ['Z44B2', 'Z44B2',
               'Z40', 'Z24', 'Z33Z2', 'Z24',
               'Z24', 'Z13', 'Z13',
               'Z13', 'Z33Z9', 'Z24',
               'Z24', 'Z13',
               'Z24', 'Z33Z9',
               'Z13',
               'Z33Z2']
#An array that contains the numbers of the various universes of fuel assemblies on the cartogram.
half_numbers = [1, 2, 2, 1,
                3, 4, 5, 6, 5, 4, 3,
                1, 4, 7, 8, 9, 9, 8, 7, 4, 1,
                2, 5, 8, 10, 11, 12, 11, 10, 8, 5, 2,
                2, 6, 9, 11, 13, 14, 14, 13, 11, 9, 6, 2,
                1, 5, 9, 12, 14, 15, 16, 15, 14, 12, 9, 5, 1,
                4, 8, 11, 14, 16, 17, 17, 16, 14, 11, 8, 4,
                3, 7, 10, 13, 15, 17]

g1 = [9]
g2 = [5]
g3 = [12]
g4 = [7, 16]
g5 = [18, 13, 8]

h1 = 0
h2 = 0
h3 = 0
h4 = 0
h5 = 0

'''
half_numbers = [1,  2,  3,  4,
           5,  6,  7,  8,  9, 10,  5,
           4, 10, 11, 12, 13, 14, 15, 11,  6,  1,
           3,  9, 15, 16, 17, 18, 19, 16, 12,  7,  2,
           2,  8, 14, 19, 20, 21, 22, 20, 17, 13,  8,  3,
           1,  7, 13, 18, 22, 23, 24, 23, 21, 18, 14,  9, 4,
           6, 12, 17, 21, 24, 25, 25, 24, 22, 19, 15, 10,
           5, 11, 16, 20, 23, 25] #Half of the fuel assemblies number on the cartogram without the central cassette
'''
r_half_numbers = list(reversed(half_numbers))
half_numbers.append(len(dif_fu_cart))
numbers = half_numbers + r_half_numbers

