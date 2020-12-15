'''modifying the Jacobian inputs per integration time interval to account for changing particle phase fractions'''
# When several components contribute similar fractions to the particle phase of a size bin, the Jacobian 
# appears to work well (relatively low processing time) with only the diagonal elements complete.  However,
# when one or two components dominate the fraction, lower processing times can be achieved by stating
# the row elements of these components in the Jacobian

import numpy as np

def jac_up(Cp, rowvals, colptrs, num_asb, num_comp, jac_part_indx, H2Oi, Cw, jac_wall_indx, ser_H2O): # define function

	# inputs: --------------------------------------------------------------------
	# Cp - particle-phase concentrations (molecules/cc (air))
	# rowvals - rows of Jacobian
	# colptrs - index of rows in rows variable associated with each Jacobian column
	# num_asb - number of particle-phase size bins
	# num_comp - number of components
	# jac_part_indx - the original indices for the sparse Jacobian matrix representing the diagonal
	# 		elements for gas-particle partitioning
	# H2Oi - index for the water component
	# Cw - gas-phase concentration of water (molecules/cc (air))
	# jac_wall_indxn - the original indices for sparse Jacobian matrix representing the diagonal
	#	elements for gas-wall partitioning
	# ser_H2O - whether water gas-particle partitioning serialised
	# ------------------------------------------------------------------------------
	
	# if no size bins or no water or water gas-particle partitioning being solved separately, then nothing to change
	if ((num_asb == 0) or (Cw == 0.) or (ser_H2O == 1)):
		return(rowvals, colptrs, jac_part_indx, 0, [], np.ones((num_asb,1))*-1, jac_wall_indx, [])
		
	# identify rows of Jacobian to be affected - these indices relate to
	# the gas-phase rows of the Jacobian
	rw_indx = (np.array(H2Oi)).repeat(num_asb)
	
	# Number of elements in Jacobian to be modified as a result of high particle-phase 
	# fraction.  Multiply by two to account for both the gas and particle.  Multiply by
	# the number of component less itself since all columns per size bin will be affected
	# except for itself where the diagonal is already accounted for
	jac_mod_len = int(num_asb)*2.*(num_comp-1)
	
	# -------------------------------------------------------------------------------------------------------------
	# number of new Jacobian elements to be included from considering the particle-phase 
	# water effect, multiply by two to include effect on both the gas and particle phase
	jac_part_H2O_eff_num = (num_comp-1)*2*num_asb
	# ------------------------------------------------------------------------------------------------------------
	
	# number of jacobian particle effect indices allocated for gas-on-gas and 
	# gas-on-particle (starting index for jacobian particle effect indices 
	# (particle-on-gas and particle-on-particle effect))
	jpi_st = num_comp*(num_asb+1)
	
	jac_mod_len = 0 # count on the new number of elements in sparse Jacobian
	
	# new array for holding the indices of the sparse jacobian matrix that represent
	# the gas-particle partitioning diagonal elements 
	jac_part_indxn = np.zeros((len(jac_part_indx)))
	jac_part_indxn[:] = (jac_part_indx[:])
	jac_part_indxn = jac_part_indxn.astype('int')
	# same for gas-wall partitioning
	jac_wall_indxn = np.zeros((len(jac_wall_indx)))
	jac_wall_indxn[:] = (jac_wall_indx[:])
	jac_wall_indxn = jac_wall_indxn.astype('int')
	
	# new array for holding row indices of the sparse Jacobian
	rowvalsn = np.zeros((len(rowvals)))
	rowvalsn[:] = rowvals[:]
	rowvalsn = rowvalsn.astype('int')
	
	# new array for holding column indices of the sparse Jacobian
	colptrsn = np.zeros((len(colptrs)))
	colptrsn[:] = colptrs[:]
	colptrsn = colptrsn.astype('int')
	
	# array to hold sparse Jacobian indices representing the particle-on-gas and particle-on-particle
	# effect on particle-phase water
	jac_part_hmf_indx = np.zeros((0), dtype = int)
	
	# array to hold sparse Jacobian indices representing the particle-on-gas and particle-on-particle
	# effect of particle-phase water
	jac_part_H2O_indx = np.zeros((0), dtype = int)
	
	
	# modify Jacobian sparse matrix inputs to account for water in the particle phase
	for isb in range(int(num_asb)): # loop through size bins
			
		# starting indices for diagonal matrix for this size bin: first term described above, 
		# second is default inputs due to particle-on-gas and particle-on-particle diagonals 
		pi_sti = int(jpi_st+isb*(num_comp*2))
		# finishing indices for diagonal matrix for this size bin: terms described above
		pi_fii = int(jpi_st+(isb+1)*(num_comp*2))
			
		# new number of elements for Jacobian, multiply by 2 because of particle
		# effect on gas and particle, subtract 1 because diagonal already 
		# accounted for
		new_num_ele = (num_comp-1)*2
			
		# length of increment array to be applied to jacobian particle-phase diagonal part,
		# multiply by two to account for the gas and particle part
		inc_len = int(num_comp*2)
			
		# prepare array holding the increment per diagonal element
		indx_incr = np.arange(inc_len)
			
		# identify where diagonal occurs as no additional jacobian elements needed here
		comp_indx = rw_indx[isb]
		
		# the increments needed for the jacobian particle index
		indx_incr[comp_indx*2+1::] -= 1
			
		# repeat final increment over any larger size bins
		indx_incr = (np.append(indx_incr, np.ones(((num_asb-(isb+1))*num_comp*2))*indx_incr[-1])).astype('int')
			
		# new indices for jacobian particle-on-gas and particle-on-particle diagonals
		jac_part_indxn = np.concatenate((jac_part_indxn[0:pi_sti], (jac_part_indxn[pi_sti::]+indx_incr)))
		# push up the wall-on-gas and gas-on-wall elements of the gas-wall partitioning indices
		jac_wall_indxn[num_comp*2::] += (num_comp-1)*2
			
		# Jacobian indices for water's particle-on-gas and particle-on-particle rows
		ji_indx = ((jac_part_indxn[pi_sti+1:pi_fii])-(jac_part_indxn[pi_sti:pi_fii-1]) == 2)
			
		# new indices for the jacobian particle-on-gas and particle-on-particle
		# rows for water
		jac_part_hmf_indx = np.concatenate((jac_part_hmf_indx, ((jac_part_indxn[pi_sti+1:pi_fii][ji_indx])-1)))
			
		jac_mod_len += new_num_ele # count on the new number of elements
		# starting rowvals index for this size bin
		rv_indx_st = colptrsn[int(num_comp*(isb+1))]
		# finishing rowvals index for this size bin
		rv_indx_fi = colptrsn[int(num_comp*(isb+2))]
			
		rowvals_pre = rowvalsn[0:rv_indx_st]
		rowvals_pos = rowvalsn[rv_indx_fi::]
		rowvals_st = rowvalsn[rv_indx_st:rv_indx_fi]
			
		# reshape rowvals for this size bin so that components are in columns, effect on gas is in first row and
		# effect on particle is in third row
		rowvals_st = np.transpose(rowvals_st.reshape(int(num_comp), 2))
		rowvals_st = np.concatenate((rowvals_st[0, :].reshape(1, -1), np.ones((1, num_comp))*-1, rowvals_st[1, :].reshape(1, -1), np.ones((1, num_comp))*-1), axis=0)
			
		# shift down the necessary gas effect elements
		rowvals_st[1, rowvals_st[0, :]>comp_indx] = rowvals_st[0, rowvals_st[0, :]>comp_indx]
		rowvals_st[1, rowvals_st[0, :]<comp_indx] = comp_indx
		rowvals_st[0, rowvals_st[1, :]>comp_indx] = comp_indx
		# shift down the necessary particle effect elements
		rowvals_st[3, rowvals_st[2, :]>(comp_indx+(num_comp*(isb+1)))] = rowvals_st[2, rowvals_st[2, :]>(comp_indx+(num_comp*(isb+1)))]
		rowvals_st[3, rowvals_st[2, :]<(comp_indx+(num_comp*(isb+1)))] = (comp_indx+(num_comp*(isb+1)))
		rowvals_st[2, rowvals_st[3, :]>(comp_indx+(num_comp*(isb+1)))] = (comp_indx+(num_comp*(isb+1)))
			
		# ensure still integer
		rowvals_st = rowvals_st.astype('int')
			
		# ensure values to be used as indices are integer
		num_comp = int(num_comp)
			
		# get the new number of row indices per component, subtract the number prior to modification
		col_incr = np.sum(rowvals_st != -1, axis=0)-(colptrs[(num_comp*(isb+1))+1:(num_comp*(isb+2)+1)]-colptrs[(num_comp*(isb+1)):(num_comp*(isb+2))])
			
		# account for new Jacobian elements in colptrs
		colptrsn[(num_comp*(isb+1)+1):(num_comp*(isb+2)+1)] += np.cumsum(col_incr)
		# also account for in Jacobian elements further on
		colptrsn[(num_comp*(isb+2)+1)::] += (np.cumsum(col_incr))[-1]
			
		# flatten new rowvals into 1D array
		rowvals_st = rowvals_st.flatten(order = 'F')
		# keep only wanted elements
		rowvals_st = rowvals_st[rowvals_st != -1]
		# concatenate back to original
		rowvalsn = np.concatenate((rowvals_pre, rowvals_st, rowvals_pos))
		
		# ------------------------------------------------------------------------------------------------
		# include effect of particle-phase water on gas and particle phase of other components
		# (the particle-phase water column in the Jacobian)
		
		# first element of the particle effect diagonal indices array to be affected, first term is the 
		# elements affected by the gas phase, second term are elements affected in smaller size 
		# bin particle phase and final term is the relevant index in this size bin
		#jpi_indx = num_comp*(num_asb+1)+(num_comp*2*isb)+H2Oi*2
		# starting index for sparse Jacobian
		#jac_part_H2O_sti = jac_part_indxn[jpi_indx]
		
		# particle-on-gas
		# new indices before diagonal
		#ni = np.arange(jac_part_H2O_sti, jac_part_H2O_sti+H2Oi)
		#jac_part_H2O_indx = np.concatenate((jac_part_H2O_indx, ni))
		# push up the indices for the particle diagonal indices
		#jac_part_indxn[jpi_indx::] += len(ni)
		# push up the indices for the particle row indices
		#jac_part_hmf_indx[H2Oi*2+(((num_comp-1)*2)*isb)::] += len(ni)
		
		# new indices after diagonal
		#ni = np.arange(jac_part_H2O_sti+H2Oi+1, (jac_part_H2O_sti+H2Oi+1)+((num_comp-1)-H2Oi))
		#jac_part_H2O_indx = np.concatenate((jac_part_H2O_indx, ni))
		# push up the indices for the particle diagonal indices
		#jac_part_indxn[(jpi_indx+1)::] += len(ni)
		# push up the indices for the particle row indices
		#jac_part_hmf_indx[H2Oi*2+(((num_comp-1)*2)*isb)::] += len(ni)
		
		# push up the wall-on-gas and gas-on-wall elements of the gas-wall partitioning indices
		#jac_wall_indxn[num_comp*2::] += (num_comp-1)
		
		# ---------------------------------------
		# particle-on-particle
		# second element of the particle effect diagonal indices array to be affected
		#jpi_indx = num_comp*(num_asb+1)+(num_comp*2*isb)+H2Oi*2+1
		# starting index of sparse Jacobian
		#jac_part_H2O_sti = jac_part_indxn[jpi_indx]
		# new indices before diagonal
		#ni = np.arange(jac_part_H2O_sti, jac_part_H2O_sti+H2Oi)
		#jac_part_H2O_indx = np.concatenate((jac_part_H2O_indx, ni))
		
		# push up the indices for the particle diagonal indices
		#jac_part_indxn[jpi_indx::] += len(ni)
		# push up the indices for the particle row indices
		#jac_part_hmf_indx[H2Oi*2+(((num_comp-1)*2)*isb)::] += len(ni)
		
		# new indices after diagonal
		#ni = np.arange(jac_part_H2O_sti+H2Oi+1, (jac_part_H2O_sti+H2Oi+1)+((num_comp-1)-H2Oi))
		#jac_part_H2O_indx = np.concatenate((jac_part_H2O_indx, ni))
		# push up the indices for the particle diagonal indices
		#jac_part_indxn[(jpi_indx+1)::] += len(ni)
		# push up the indices for the particle row indices
		#jac_part_hmf_indx[H2Oi*2+(((num_comp-1)*2)*isb)::] += len(ni)
		
		# push up the wall-on-gas and gas-on-wall elements of the gas-wall partitioning indices
		#jac_wall_indxn[num_comp*2::] += (num_comp-1)
		
		# -----------------------------------------------------
		# increment rowvals and colptrs in response to new sparse Jacobian points
		#rowvals_pre = rowvalsn[0:colptrsn[num_comp+(num_comp*isb)+H2Oi]]
		#rowvals_pro = rowvalsn[colptrsn[num_comp+(num_comp*isb)+H2Oi+1]::]
		# note that with the two segments of the original rowvals above (rowvals_pre & rowvals_pro) 
		# the diagonal is omitted so need to account for this in the new rowals
		#rowvals_new = np.concatenate((np.arange(0, num_comp), np.arange(num_comp*(isb+1), num_comp*(isb+2)))) 
		#rowvalsn = np.concatenate((rowvals_pre, rowvals_new, rowvals_pro))
		#colptrsn[(num_comp+(num_comp*isb)+H2Oi+1)::] += (num_comp-1)*2
		 
		#jac_mod_len += (num_comp-1)*2 # count on the new number of elements in sparse Jacobian
	
	return(rowvalsn, colptrsn, jac_part_indxn, jac_mod_len, jac_part_hmf_indx, rw_indx, jac_wall_indxn, jac_part_H2O_indx)