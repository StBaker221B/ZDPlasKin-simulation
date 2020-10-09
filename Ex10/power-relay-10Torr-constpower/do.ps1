
gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse01.f90 bolsig_x86_64.dll -o test01

gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse02.f90 bolsig_x86_64.dll -o test02

gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse03.f90 bolsig_x86_64.dll -o test03

gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse04.f90 bolsig_x86_64.dll -o test04

gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse05.f90 bolsig_x86_64.dll -o test05

gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse06.f90 bolsig_x86_64.dll -o test06

# gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse07.f90 bolsig_x86_64.so -o test07

# gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse08.f90 bolsig_x86_64.so -o test08

# gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse09.f90 bolsig_x86_64.so -o test09

# gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse10.f90 bolsig_x86_64.so -o test10

# export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH

.\test01
New-Item -path E:\Work\zdplaskin\Ex10\power-relay-10Torr -name 01 -type directory
move-item -path .\qt_*.txt -destination .\01

.\test02
New-Item -path E:\Work\zdplaskin\Ex10\power-relay-10Torr -name 02 -type directory
move-item -path .\qt_*.txt -destination .\02

.\test03
New-Item -path E:\Work\zdplaskin\Ex10\power-relay-10Torr -name 03 -type directory
move-item -path .\qt_*.txt -destination .\03

.\test04
New-Item -path E:\Work\zdplaskin\Ex10\power-relay-10Torr -name 04 -type directory
move-item -path .\qt_*.txt -destination .\04

.\test05
New-Item -path E:\Work\zdplaskin\Ex10\power-relay-10Torr -name 05 -type directory
move-item -path .\qt_*.txt -destination .\05

.\test06
New-Item -path E:\Work\zdplaskin\Ex10\power-relay-10Torr -name 06 -type directory
move-item -path .\qt_*.txt -destination .\06

