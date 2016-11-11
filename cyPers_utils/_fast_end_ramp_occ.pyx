import numpy as np
import cython

@cython.cdivision(True) 
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)

def _fast_end_ramp_occ(rtime,rcts,cmin,a,b,t_rel):

  cdef int counts,i,ntimes,ntraps
  cdef double dt_o,dt_n

  ntimes = rtime.size
  ntraps = cmin.size
  
  totfill = np.zeros((ntimes-1),dtype=np.float_)
  above_cut_o   = np.zeros(ntraps,dtype=np.bool_)
  above_cut_n   = np.zeros(ntraps,dtype=np.bool_)
  below_cut_o   = np.zeros(ntraps,dtype=np.bool_)
  below_cut_n   = np.zeros(ntraps,dtype=np.bool_)
  states_o      = np.zeros(ntraps,dtype=np.bool_)
  states_n      = np.zeros(ntraps,dtype=np.bool_)
  dt_o          = 0.
  dt_n          = 0.
  occ_prob      = np.zeros(ntraps,dtype=np.float_)


  for i in range(ntimes-1):
    counts    = rcts[i]
    dt_n      = rtime[i] - rtime[i+1]
    
    above_cut_n = cmin <= counts
    below_cut_n = np.logical_not(above_cut_n)
    
    if (dt_o == dt_n):
      d_occ  = np.logical_xor(states_o,states_n)
      diff_a = np.logical_and(above_cut_n,np.logical_or(np.logical_xor(above_cut_n,above_cut_o),d_occ))
      diff_b = np.logical_and(below_cut_n,np.logical_or(np.logical_xor(below_cut_n,below_cut_o),d_occ))
      
    else:
      diff_a  = above_cut_n
      diff_b  = below_cut_n
      
    exp1 = np.exp(a[diff_a]*dt_n)
    occ_prob[diff_a] = states_n[diff_a] * exp1 + b[diff_a]*(1-exp1)
      
    exp2 = states_n[diff_b] * np.exp( dt_n / t_rel[diff_b] )
    occ_prob[diff_b] = exp2
    
    states_o = states_n
    states_n = np.random.uniform(size=ntraps)  < occ_prob
    totfill[i]  = np.sum(states_n)
    above_cut_o   = above_cut_n
    below_cut_o   = below_cut_o
    dt_o          = dt_n
    

  return totfill, states_n
