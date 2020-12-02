# LAMMPS input script for melting point of 3d FCC Cu
# Volume--Tempearature relation
units		metal
boundary	p p p
atom_style	atomic

lattice         fcc 3.615
region          box block 0 6 0 6 0 6
create_box      1 box
create_atoms    1 box

timestep	0.005
neighbor        2.0 bin
neigh_modify	delay 0 every 5 


pair_style      eam/alloy
pair_coeff      * * ../CuYM.eam.alloy Cu

# vars for V-T relation
#initial Temp and velocity is 10K 
variable	step equal step
variable	time equal step*dt
variable	temp equal temp
variable	vol equal vol
variable	Et equal etotal
variable	press equal press

variable	Tini equal 10
variable	Tdamp equal 1
velocity        all create ${Tini} 1234567

compute		pea all pe/atom
compute		kea all ke/atom

thermo          1000
thermo_style	custom step temp vol v_Et pe ke press

# equilibration run
dump		1 all custom 10000 rCu*.xyz id type x y z c_pea c_kea 
fix		1 all npt temp ${Tini} ${Tini} ${Tdamp} iso 0 0 10 
run	        50000
unfix		1
undump		1

# Heating from initial 10K to 2000K
reset_timestep  0
compute		csp all centro/atom fcc
compute		rdf all rdf 100
fix		rdf all ave/time 100 10 1000 c_rdf[*]  &
		ave running file CuMelt.rdf mode vector
dump		1 all custom 10000 Cu*.xyz id type x y z c_pea c_kea c_csp
fix             print all print 1000 "${step} ${temp} ${vol} ${Et} ${press}" file CuVT.dat

fix		1 all npt temp ${Tini} 2000 10 iso 0 0 10

thermo          1000
run	        140000

