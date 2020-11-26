

gap = 1e-3 # m
width_electrode = 10e-3 # m
width_cool = 10e-3 # m
ring_number = 6
massflow = 20 # sccm
pressure = 10 # torr


pi = 3.14159265358
r_out = 10e-3 # m
r_in = r_out - gap


torr2pa = 101325/760

s = pi * (r_out**2 - r_in**2) # m^2
volume = s * width_electrode * ring_number # m^3
print("\nvolume\n")
print(volume)

volflow = massflow * 760 / pressure # ml/min

vel = volflow * 1e-6 / 60 / s # m/s
print("\nvel\n")
print(vel)

time_electrode = width_electrode / vel # s
time_cool = width_cool / vel # s

print("\ntime_electrode\n")
print(time_electrode)
print("\ntime_cool\n")
print(time_cool)


