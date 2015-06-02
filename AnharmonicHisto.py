import numpy
import numpy as np
import math
import scipy.stats
import matplotlib.pyplot as plt
from KabschAlign import *
from IterativeMeansAlign import *
from MDAnalysis.core.Timeseries import *
from MDAnalysis import *
from numpy import *
from jade import *
import deps.timing as timing
import argparse

#============================================================
#a is mem-mapped array, b is array in RAM concatenated to a.       
def mmap_concat(a,b,res):
	c = np.memmap('dihedral_data.array', dtype='float32', mode='r+', shape=(4*res,a.shape[1]+b.shape[1]), order='F')
	c[:, a.shape[1]: ] = b
	return c

def get_phi(res):
#MDAnalysis' phi_selection requires a segid be present, this doesn't.  Use if MDanalysis fails.
	return res.universe.selectAtoms('resid %d and name C' %(res.id-1) ) + res.N + res.CA + res.C

def get_psi(res):
#MDAnalysis' psi_selection requires a segid be present, this doesn't.  Use if MDanalysis fails.
	return res.N + res.CA + res.C + res.universe.selectAtoms('resid %d and name N' %(res.id+1) )

#============================================================


#============================================================
#Main Code
def generate_ica(num_traj, opts):
	rad_gyr = []
	for i in range(num_traj):
	# This was to satisfy UBQ data formatting
	#	if i < 10:
	#		tm = '00';
	#	elif i < 100:
	#		tm = '0';
	#	else:
	#		tm = '';
		
		#Path to data
		#u = MDAnalysis.Universe("../data/kbh_pdb/1KBH%i_ww.pdb" %(i+1), "../data/lacie-kbhdata/1KBH_%i_50k.dcd" %(i+1), permissive=False);
		u = MDAnalysis.Universe("../../data/kbh_pdb/1KBH%i_ww.pdb" %(i+1), "lacie-kbhdata/1KBH_%i_50k.dcd" %(i+1), permissive=False);
		#u = MDAnalysis.Universe("lacie-ubqdata/protein.pdb", "lacie-ubqdata/pnas2013-native-1-protein-%s.dcd" %(tm+str(i)), permissive=False);
		atom = u.selectAtoms("backbone");

		phidat = TimeseriesCollection()
		psidat = TimeseriesCollection()

		timing.log('Processing Trajectory %d...'%(i+1))
		
		num_res = 73
		len_traj = len(u.trajectory)
		if opts.debug: print len_traj
		
		#Adjust to residue selection.
		for res in range(1,num_res+1):
			if opts.debug: print res
			#  selection of the atoms involved for the phi for resid res
			try: phi_sel = u.residues[res].phi_selection()
			except: phi_sel = get_phi(u.residues[res])

			#  selection of the atoms involved for the psi for resid res
			try: psi_sel = u.residues[res].psi_selection()
			except: psi_sel = get_psi(u.residues[res])

			#  adds selection to Timeseries object
			if len(phi_sel) == 4: phidat.addTimeseries(Timeseries.Dihedral(phi_sel)) 
			if len(psi_sel) == 4: psidat.addTimeseries(Timeseries.Dihedral(psi_sel))
		
		#Computes along timesteps
		phidat.compute(u.trajectory)
		psidat.compute(u.trajectory)
	
		#Converts to nd-array and 'flattens' (3d->2d)
		phidat =  array(phidat)
		phidat = phidat.reshape(phidat.shape[0],phidat.shape[2])
		psidat =  array(psidat)
		psidat = psidat.reshape(psidat.shape[0],psidat.shape[2])
		
		#np.save('psidat.npy', psidat)
		#np.save('phidat.npy', phidat)
			
		
		dihedral_dat = np.zeros((num_res,4,len_traj))
		#Data stored as | sin(phi) | cos(phi) | sin(psi) | cos(psi) |
		dihedral_dat[:,0,:] = np.sin(phidat)
		dihedral_dat[:,1,:] = np.cos(phidat)
		dihedral_dat[:,2,:] = np.sin(psidat)
		dihedral_dat[:,3,:] = np.cos(psidat)
		dihedral_dat = dihedral_dat.reshape(-1,len_traj)
		
		#Stores the array on disk rather than in memory.  Allows arbitrary data size
		if i == 0:
			fulldat = np.memmap('dihedral_data.array', dtype='float32', mode='w+', shape=(num_res*4, len_traj))
			fulldat[:,:] = dihedral_dat
		else:
			fulldat = mmap_concat(fulldat, dihedral_dat, num_res);
		for ts in u.trajectory:
			rad_gyr.append( atom.radiusOfGyration() )

	if opts.debug: print 'fulldat: ', fulldat.shape	
	
	# Some set up for running JADE
	Ncyc  = 1;
	subspace = 50;
	lastEig = subspace; # number of eigen-modes to be considered
	numOfIC = subspace; # number of independent components to be resolved
	
	# Prints the radius of Gyration of protein over time
	if opts.graph: 
		plt.plot(range(len(rad_gyr)),rad_gyr[:], 'r--', lw=2)
		plt.show()

	#icajade = np.memmap('icajade.array', dtype='float32', mode='w+', shape=(276, 10000))
	timing.log('All dihedral data stored, starting jade...')
	
	#JADE
	icajade = jadeR(fulldat, lastEig);
	if opts.debug: print 'icajade: ', numpy.shape(icajade);
	if opts.save: np.save('../save_data/icajade_PNAME.npy', icajade);

	#Scaling the data with ICA's cof
	icacoffs = icajade.dot(fulldat);
	icacoffs = numpy.asarray(icacoffs); 
	print 'icacoffs: ', numpy.shape(icacoffs);
	if opts.save: np.save('../save_data/icacoffs_PNAME.npy', icacoffs);
	
	if opts.graph:
		fig = plt.figure();
		ax = fig.add_subplot(111, projection='3d');
		ax.scatter(icacoffs[0,::100], icacoffs[1,::100], icacoffs[2,::100], marker='o', c=[0.6,0.6,0.6]); 
		print 'First 3-dimensions of ICA';
		plt.show();

	return icacoffs;
#============================================================

if __name__ == '__main__':
	timing.begin()
	parser = argparse.ArgumentParser();
	parser.add_argument('-g', action='store_true', dest='graph', default=False, help='Shows graphs that are produced, depending on flags');
	parser.add_argument('-v', action='store_true', dest='verbose', default=False, help='Runs program verbosely');
	parser.add_argument('-s', '--save', action='store_true', dest='save', default=False, help='Saves all data in ./save_data');
	parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='Prints debugging help.');

	values = parser.parse_args();
	global val;
	val = values;

	generate_ica(5,val);

