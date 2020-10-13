
# gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse01.f90 bolsig_x86_64.so -o test01

# gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse02.f90 bolsig_x86_64.so -o test02

# gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse03.f90 bolsig_x86_64.so -o test03

# gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse04.f90 bolsig_x86_64.so -o test04

gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse05.f90 bolsig_x86_64.so -o test05

gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse06.f90 bolsig_x86_64.so -o test06

gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse07.f90 bolsig_x86_64.so -o test07

gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse08.f90 bolsig_x86_64.so -o test08

gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse09.f90 bolsig_x86_64.so -o test09

gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse10.f90 bolsig_x86_64.so -o test10

export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH

# ./test01
# mkdir ./01
# mv ./qt_*.txt ./01

# ./test02
# mkdir ./02
# mv ./qt_*.txt ./02

# ./test03
# mkdir ./03
# mv ./qt_*.txt ./03

# ./test04
# mkdir ./04
# mv ./qt_*.txt ./04

./test05
mkdir ./05
mv ./qt_*.txt ./05

./test06
mkdir ./06
mv ./qt_*.txt ./06

./test07
mkdir ./07
mv ./qt_*.txt ./07

./test08
mkdir ./08
mv ./qt_*.txt ./08

./test09
mkdir ./09
mv ./qt_*.txt ./09

./test10
mkdir ./10
mv ./qt_*.txt ./10
