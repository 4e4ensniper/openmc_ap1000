split_number = 3
core_height = 3.53 #m
n_fa = 151
r_fuel = 7.57/20 #sm
r_hole = 2.35/20 #sm
fr_number = 312
rod_pitch = 1.275 #sm
turnkey_size = 23.4 #sm
core_barrel_in_r = 339.72/2 #sm
core_barrel_out_r = 349.88/2 #sm
#The number of fuel assemblies per row, starting from the bottom of the core map.
line = [4, 7, 10, 11, 12, 13, 12, 13, 12, 13, 12, 11, 10, 7, 4]
csv_path = "/home/ubuntu24/Desktop/openmc_ap1000/materials/temperature_distributions/"
batches = 200
inactive = 10
particles = 10000
n_dif = 26

half_numbers = [1,  2,  3,  4,
           5,  6,  7,  8,  9, 10,  5,
           4, 10, 11, 12, 13, 14, 15, 11,  6,  1,
           3,  9, 15, 16, 17, 18, 19, 16, 12,  7,  2,
           2,  8, 14, 19, 20, 21, 22, 20, 17, 13,  8,  3,
           1,  7, 13, 18, 22, 23, 24, 23, 21, 18, 14,  9, 4,
           6, 12, 17, 21, 24, 25, 25, 24, 22, 19, 15, 10,
           5, 11, 16, 20, 23, 25] #Half of the fuel assemblies number on the cartogram without the central cassette
r_half_numbers = list(reversed(half_numbers))
half_numbers.append(n_dif)
numbers = half_numbers + r_half_numbers
