import numpy as np
from Pers_utils.Ramp import Ramp
from Pers_utils.MyExceptions import CustExc1M3V
import copy

class cPixTraps(object):

    '''
    Class to describe the ensemble of traps for a pixel
    in the persistence model

    :t_trp:
        Trapping times in sec

    :t_rel:
        Release times in sec

    :cmin:
        Number of counts in the diode
        that are necessary to expose
        each trap to free charge (i.e. if counts < cmin[i]
        the i-th trap cannot be filled)

    '''

    def __init__(self, t_trp, t_rel, cmin):
        self.t_trp = t_trp
        self.t_rel = t_rel
        self.cmin  = cmin

        if ( ( self.t_trp.shape[0] == self.t_rel.shape[0] == self.cmin.shape[0]) == False):
            raise CustExc1M3V('Trapping times, release times and min counts numpy.ndarrays must be of the same shape instead they are resp.',
                              self.t_trp.shape, self.t_rel.shape, self.cmin.shape)
        else:
            self.ntraps = self.t_trp.shape[0]

        '''
        Initialise the traps to empty
        '''
        self.states   = np.zeros(self.ntraps,dtype=np.bool_)
        self.occ_prob = np.zeros(self.ntraps,dtype=np.float_)

        '''
        Convenience variables for the diff. equations
        '''
        self.a = (self.t_trp + self.t_rel)/(self.t_trp*self.t_rel)
        self.b = self.t_rel/(self.t_trp + self.t_rel)
        self.above_cut   = np.zeros(self.ntraps,dtype=np.bool_)
        self.below_cut   = np.zeros(self.ntraps,dtype=np.bool_)
        self.prev_states = np.zeros(self.ntraps,dtype=np.bool_)
        self.dt          = 0.

        
    def end_ramp_occ(self,rmp):
        '''
        Method that computes the occupation of each trap in the pixel
        at the end of the rmp Ramp.  If a trap is filled at that
        time, then it will produce an electron with exponentially
        decaying probability. If it is empty it won't.

        :rmp:
            rmp is a Ramp object that describes the previous history
            of the pixel (in raw counts) up to the time
            when we start measuring persistence.
        '''

        cdef int counts,i,ntimes
        cdef float dt

        ntimes = len(rmp.rtime)
        self.totfill = np.zeros((ntimes-1),dtype=np.float_)
        rcts = rmp.rcts
        rtime = rmp.rtime

        for i in range(ntimes-1):
            counts    = rcts[i]
            dt        = rtime[i] - rtime[i+1]
            
            above_cut = self.cmin <= counts
            below_cut = np.logical_not(above_cut)

            if (dt == self.dt):
                d_occ  = np.logical_xor(self.prev_states,self.states)
                diff_a = np.logical_and(above_cut,np.logical_or(np.logical_xor(above_cut,self.above_cut),d_occ))
                diff_b = np.logical_and(below_cut,np.logical_or(np.logical_xor(below_cut,self.below_cut),d_occ))

            else:
                diff_a  = above_cut
                diff_b  = below_cut

            exp1 = np.exp(self.a[diff_a]*dt)
            self.occ_prob[diff_a] = self.states[diff_a].astype(np.float_) * exp1 + self.b[diff_a]*(1-exp1)

            exp2 = self.states[diff_b] * np.exp( dt / self.t_rel[diff_b] )
            self.occ_prob[diff_b] = exp2

            self.prev_states = self.states
            self.states      = np.random.random_sample(size=self.ntraps)  < self.occ_prob
            self.totfill[i]  = np.sum(self.states)
            self.above_cut   = above_cut
            self.below_cut   = below_cut
            self.dt = dt
                               
    def reset(self):
        '''
        Method to reset the occupancy of all traps to 0
        '''

        self.states[:]   = 0. 
        self.occ_prob[:] = 0.


    def get_trapped_charge(self,t_after):
        '''
        Method to measure the total charge still trapped at times t_after
        after the end of the ramp. It assumes the diode
        is no longer being exposed to light
        '''
        charge = np.zeros(len(t_after)-1)
        usert   = copy.deepcopy(self.t_rel[self.states])
 
        for i in range(len(t_after)-1):
            pmean     = (t_after[i+1]-t_after[i]) / usert
            exp       = np.exp( -1.*pmean )
            ikeep     = np.random.random_sample(size=exp.size) < exp
            charge[i] = np.sum(ikeep)
            usert     = usert[ikeep]
 
        return charge


    def follow_occ_prob(self,rmp,times):
        '''
        Method that computes the occupation probability
        of each trap in the pixel at given times.

        :rmp:
            rmp is a Ramp object that describes the  history
            of the pixel (in raw counts) up to the time
            when we start measuring persistence.

        :times:
            numpy array of times at which the occupancy is desired
        '''

        cdef int j,k,nshifts,ntimes
        cdef double f0,dt,exp1

        rcts   = rmp.rcts
        rtime  = rmp.rtime
        ntimes = times.size
        op     = np.zeros([self.ntraps,ntimes])

        for i in range(self.ntraps):

            '''
            Check when the trap switches between exposed and not-eposed
            '''
            above   = rcts > self.cmin[i]
            shifts  = np.logical_xor(above[1:],above[:-1])   # Indices of the last rtime before the shift
            nshifts = np.sum(shifts)

            '''
            Compute the occupancy probability up to each shift point
            '''
            shift_times = rtime[0:-1][shifts]
            states      = above[0:-1][shifts]
            shift_times = np.insert(shift_times,0,0.)
            states      = np.insert(states,0,above[0])

            f0 = 0.  # all traps are empty a time = 0
            j  = 0
            for k in range(nshifts):

                if (states[k+1] == True):
                    while (times[j] <= shift_times[k+1]):
                        dt  = shift_times[k]-times[j]
                        exp1 = np.exp(self.a[i]*dt)
                        op[i,j] = f0 * exp1 + self.b[i]*(1-exp1)
                        j = j +1
                    dt = shift_times[k]- shift_times[k+1]
                    exp1 = np.exp(self.a[i]*dt)
                    f0 = f0 * exp1 + self.b[i]*(1-exp1)
                else:
                    while (times[j] <= shift_times[k+1]):
                        dt  = shift_times[k]-times[j]
                        op[i,j] = f0 * np.exp( dt / self.t_rel[i] )
                        j = j +1
                    dt = shift_times[k]- shift_times[k+1]
                    f0 = f0 * np.exp( dt / self.t_rel[i] )

            while (j < ntimes):
                dt  = shift_times[-1]-times[j]
                op[i,j] = f0 * np.exp( dt / self.t_rel[i] )
                j = j +1
            
        return op


    def get_trapped_charge_diff(self,rmp,t_after):
        '''
        Method to measure the total charge that comes out of traps
        by times t_after after the end of the ramp rmp. It assumes the diode
        is no longer being exposed to light
        '''

        '''
        Get the occupancy probabilities at the end of the ramp and
        generate random numbers to get which traps are actually filled
        '''
        op0 = np.squeeze(self.follow_occ_prob(rmp,np.array([rmp.rtime[-1]])))
        self.states = np.random.random_sample(size=self.ntraps)  < op0

        '''
        Generate the instant of decay for the traps that are filled
        '''
        tdecay = -1 * self.t_rel[self.states]*np.log(np.random.random_sample(size=np.sum(self.states)))
        hist, edges = np.histogram(tdecay,bins=np.insert(t_after,0,-1.))
        charge = np.sum(self.states) - np.cumsum(hist)

        return charge

