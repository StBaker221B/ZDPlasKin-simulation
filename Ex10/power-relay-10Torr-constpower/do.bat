
gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse01.f90 bolsig_x86_64_g.dll -o test01

REM gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse02.f90 bolsig_x86_64_g.dll -o test02

REM gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse03.f90 bolsig_x86_64_g.dll -o test03

REM gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse04.f90 bolsig_x86_64_g.dll -o test04

REM gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse05.f90 bolsig_x86_64_g.dll -o test05

REM gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse06.f90 bolsig_x86_64_g.dll -o test06

REM gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse07.f90 bolsig_x86_64.so -o test07

REM gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse08.f90 bolsig_x86_64.so -o test08

REM gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse09.f90 bolsig_x86_64.so -o test09

REM gfortran dvode_f90_m.f90 zdplaskin_m.F90 option_m.f90 input_MCSD.f90 pulse10.f90 bolsig_x86_64.so -o test10

REM export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH

.\test01
mkdir .\01
move .\qt_*.txt .\01

REM .\test02
REM mkdir .\02
REM move .\qt_*.txt .\02

REM .\test03
REM mkdir .\03
REM move .\qt_*.txt .\03

REM .\test04
REM mkdir .\04
REM move .\qt_*.txt .\04

REM .\test05
REM mkdir .\05
REM move .\qt_*.txt .\05

REM .\test06
REM mkdir .\06
REM move .\qt_*.txt .\06

