# LAMMPS input script for diffusion of 3d FCC Cu
# mean-squared displacement via compute msd

# setup
variable	x equal 10
variable	y equal 10
variable	z equal 10

variable	a equal 3.615
variable        t equal 0.001
variable	T equal 1620
variable	Tdamp equal 100*$t
variable	p equal 1
variable	pdamp equal 1000*$t

units		metal
dimension	3
atom_style	atomic
timestep	$t
neighbor        2.0 bin
neigh_modify	delay 0 every 1

lattice         fcc $a
region          box block 0 $x 0 $y 0 $z
create_box      1 box
create_atoms    1 box

pair_style      eam/alloy
pair_coeff      * * ../CuYM.eam.alloy Cu

velocity        all create $T 1234567
fix		process all npt temp $T $T ${Tdamp} iso $p $p ${pdamp}

# equilibration run
thermo          1000
run	        50000

# msd = [dx, dy, dz,u^2=(dx^2+dy^2+dz^2)]
# summed and averaged displacements of atoms in the group
reset_timestep  0

variable	time equal step*$t

compute         msd all msd com yes
variable	msd equal c_msd[4] # c_msd[4] = u^2
fix             msd_tmp all vector 10 v_msd

variable        fitD equal slope(f_msd_tmp)/6/(10*dt) # in [A^2/ps] 
variable        D2ps equal (v_msd)/6/(step*dt+1.0e-6)
variable	aveD equal 0.5*(v_D2ps+v_fitD)*10 # in [nm*m/s]=1e-9m^2/s

thermo_style	custom step v_msd v_fitD v_D2ps temp press vol
fix             print all print 10 &
"${time} ${msd} ${fitD} ${D2ps} ${aveD}" file tMSD.dat

dump		myDump all atom 10000 dump.atom

thermo          10
run	        10000

