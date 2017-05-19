from __future__ import division
import os
import sys
import math
from hoomd_script import *

c = context.initialize()

buffer = option.get_user()
ChainLength = 18
SolvLength= 6
system=init.read_xml(filename='init.xml')

# pass --user=n as arguments, where n is the replication
# factor along the x-y dimensions
rep = int(buffer[0])

system.replicate(nx=rep, ny=rep, nz=1)
typeW1=group.type('W1')
typeW2=group.type('W2')
typeWall = group.union( name = "Wall" , a = typeW1 , b = typeW2)
typeB= group.type('B')
typeSolvBackBone = group.type('Solvent')
typeC= group.type('C')
typeCH3_W1=group.type('CH3_W1')
typeCH3_W2=group.type('CH3_W2')
typeSolvTail=group.type('SolvTail')
typeB_CH3W1=group.union(name="BCH3-particles", a = typeB , b = typeCH3_W1)
typeC_CH3W2=group.union(name="CCH3-particles", a = typeC , b = typeCH3_W2)
typeSolv= group.union(name="all-solvent" , a = typeSolvBackBone , b = typeSolvTail)
typeTethers= group.union(name="bc-particles" , a=typeB_CH3W1 , b = typeC_CH3W2)
typeAll= group.union(name="typeAll" , a=typeTethers , b = typeSolv)
typeNPT = group.union(name="typeNPT", a = typeAll  , b = typeWall)

#interaction between tethers

globals.msg.notice(1,"Number of polymeric particles %d\n" % (len(typeAll)))

tether_bulk = 46.0 / 300.0
tails = 98./300.0
bulk_tails =  67.0 / 300.0

#lj = pair.lj(r_cut=3)
lj = pair.lj(r_cut=3.3)
#nl = nlist.cell(r_buff=0.4)
#lj = pair.lj(r_cut=3.3, nlist=nl)
#lj = pair.lj(r_cut=4.0, nlist=nl)
#lj = pair.lj(r_cut=4.0)

#lj.set_params(mode="shift")

lj.pair_coeff.set('Solvent','W1' ,epsilon=0 , sigma=1.3144)
lj.pair_coeff.set('Solvent','W2' ,epsilon=0 , sigma=1.3144)
lj.pair_coeff.set('Solvent','Solvent' ,epsilon=1.0*tether_bulk , sigma=1.3144)
lj.pair_coeff.set('Solvent','B' ,epsilon=tether_bulk , sigma=1.3144)
lj.pair_coeff.set('Solvent','C' ,epsilon=tether_bulk , sigma=1.3144)
lj.pair_coeff.set('Solvent','CH3_W1' ,epsilon=bulk_tails , sigma=1.3144)
lj.pair_coeff.set('Solvent','CH3_W2' ,epsilon=bulk_tails , sigma=1.3144)
lj.pair_coeff.set('Solvent','SolvTail' ,epsilon=bulk_tails , sigma=1.3144)
lj.pair_coeff.set('B','B' ,epsilon=tether_bulk , sigma=1.3144)
lj.pair_coeff.set('C','C' ,epsilon=tether_bulk , sigma=1.3144)
lj.pair_coeff.set('B','C' ,epsilon=tether_bulk , sigma=1.3144)
lj.pair_coeff.set('B','CH3_W1' ,epsilon=bulk_tails , sigma=1.3144)
lj.pair_coeff.set('B','CH3_W2' ,epsilon=bulk_tails , sigma=1.3144)
lj.pair_coeff.set('C','CH3_W1' ,epsilon=bulk_tails , sigma=1.3144)
lj.pair_coeff.set('C','CH3_W2' ,epsilon=bulk_tails , sigma=1.3144)
lj.pair_coeff.set('CH3_W1','CH3_W1' ,epsilon=tails , sigma=1.3144)
lj.pair_coeff.set('CH3_W1','CH3_W2' ,epsilon=tails , sigma=1.3144)
lj.pair_coeff.set('CH3_W2' ,'CH3_W2' ,epsilon=tails , sigma=1.3144)

lj.pair_coeff.set('B','SolvTail' ,epsilon=bulk_tails , sigma=1.3144)
lj.pair_coeff.set('C','SolvTail' ,epsilon=bulk_tails , sigma=1.3144)
lj.pair_coeff.set('CH3_W1','SolvTail' ,epsilon=tails , sigma=1.3144)
lj.pair_coeff.set('CH3_W2','SolvTail' ,epsilon=tails , sigma=1.3144)
lj.pair_coeff.set('SolvTail' ,'SolvTail' ,epsilon=tails , sigma=1.3144)

lj.pair_coeff.set('B','C' ,epsilon=0.01 , sigma=2.62 , alpha=0.0)
lj.pair_coeff.set('W1','W2' ,epsilon=0.0 , sigma=1.0)
lj.pair_coeff.set('W1','W1' ,epsilon=0.0 , sigma=1.0)
lj.pair_coeff.set('W1','B' ,epsilon=0.0 , sigma = 1.0)
lj.pair_coeff.set('W1','C' ,epsilon=0.0 , sigma = 1.0)
lj.pair_coeff.set('W2','W2' ,epsilon=0.0 , sigma=1.0)
lj.pair_coeff.set('W2','B' ,epsilon=0.0 , sigma = 1.0)
lj.pair_coeff.set('W2','C' ,epsilon=0.0 , sigma = 1.0)
lj.pair_coeff.set('W1','CH3_W1' ,epsilon=0.0 , sigma = 1.0)
lj.pair_coeff.set('W1','CH3_W2' ,epsilon=0.0 , sigma = 1.0)
lj.pair_coeff.set('W1','SolvTail' ,epsilon=0.0 , sigma = 1.0)
lj.pair_coeff.set('W2','CH3_W1' ,epsilon=0.0 , sigma = 1.0)
lj.pair_coeff.set('W2','CH3_W2' ,epsilon=0.0 , sigma = 1.0)
lj.pair_coeff.set('W2','SolvTail' ,epsilon=0.0 , sigma = 1.0)
lj.pair_coeff.set('A','C' ,epsilon=0.0 , sigma=1.0 , alpha=0.0)

slj = pair.slj(r_cut= 1.3)
slj.set_params(mode="shift")
slj.pair_coeff.set('W2','Solvent', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('W1','Solvent', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('Solvent','Solvent', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('Solvent','B', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('Solvent','C', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('Solvent','CH3_W1', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('Solvent','CH3_W2', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('Solvent','SolvTail', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('W1','W2', epsilon = 0 , sigma =1.0)
slj.pair_coeff.set('W1','W1', epsilon = 0 , sigma =1.0)
slj.pair_coeff.set('W2','W2', epsilon = 0 , sigma =1.0)
slj.pair_coeff.set('W1','B', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('W1','C', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('W2','B', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('W2','C', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('W1','CH3_W1', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('W1','CH3_W2', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('W2','CH3_W1', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('W2','CH3_W2', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('W1','SolvTail', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('W2','SolvTail', epsilon = 1.0 , sigma =1.25 ,r_cut =1.25* 2**(1.0/6.0))
slj.pair_coeff.set('B','B', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('C','C', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('B','C', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('B','CH3_W1', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('B','CH3_W2', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('B','SolvTail', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('C','CH3_W1', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('C','CH3_W2', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('C','SolvTail', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('CH3_W1','CH3_W1', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('CH3_W1','CH3_W2', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('CH3_W2','CH3_W2', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('CH3_W1','SolvTail', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('CH3_W2','SolvTail', epsilon = 0, sigma = 1.0)
slj.pair_coeff.set('SolvTail','SolvTail', epsilon = 0, sigma = 1.0)

#Oleic-Oleic Bond
harmonicO_O= bond.harmonic(name = 'O-O-bond')
#harmonicO_O.bond_coeff.set('oleic', k=1063.47 , r0=0.38)
harmonicO_O.bond_coeff.set('oleic', k=9892.9362 , r0=0.5)

#Wall-Oleic Bond
harmonicO_O.bond_coeff.set('wall', k=9892.9362, r0=1.25)

#Angle Forces
angleWO = angle.harmonic()
angleWO.set_coeff('W-O' , k = 201.2878 , t0= 2.5307)
angleWO.set_coeff('O-O' , k = 201.2878 , t0= 1.92)
angleWO.set_coeff('kink' , k = 1000.2878 , t0= 2.8)

#Dihedral Forces
#def harmonic(theta , k):
#    return (0.5 * k * theta*theta, -k*theta)

dihedralField = dihedral.harmonic()
dihedralField.set_coeff('DihedralW1' , k = 2.6838 , d=1 , n= 1)
dihedralField.set_coeff('DihedralW2' , k = 1.4543 , d=-1 , n= 2)
dihedralField.set_coeff('DihedralW3' , k =5.4347 , d=1 , n= 3)
dihedralField.set_coeff('Dihedral1' , k =2.6838 , d=1 , n= 1)
dihedralField.set_coeff('Dihedral2' , k = 1.4543 , d=-1 , n= 2)
dihedralField.set_coeff('Dihedral3' , k =5.4347 , d=1 , n= 3)

integrate.mode_standard(dt= 0.005)

integrator=integrate.nve(group=typeAll , limit=0.01 )
zeroer = update.zero_momentum(period = 1)


run(200)
zeroer.disable()

integrator.disable()

integrator=integrate.npt(group=typeAll , tau=1.0 , T= 1.0 , tauP = 1.2 , P = 0.00078 , x=False , y = False , z = True , rescale_all=True )

#equilibrate NPT
run(1e4, profile=True, limit_hours=1)
#run(5e4)

#measure
#run(5e3)
#run(1e4)
