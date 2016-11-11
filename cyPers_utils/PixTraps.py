import numpy as np
from Pers_utils.Ramp import Ramp
from Pers_utils.MyExceptions import CustExc1M3V
import copy
from numba import jit
from cyPers_utils._fast_end_ramp_occ import _fast_end_ramp_occ

class PixTraps(object):

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

        '''
        Convenience variables for the diff. equations
        '''
        self.a = (self.t_trp + self.t_rel)/(self.t_trp*self.t_rel)
        self.b = self.t_rel/(self.t_trp + self.t_rel)
        
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

        totfill, states = _fast_end_ramp_occ(rmp.rtime,rmp.rcts.astype(np.int),self.cmin,self.a,self.b,self.t_rel)
        self.totfill = np.array(totfill,dtype=np.float_)
        self.states  = np.array(states,dtype=np.bool_)

                               
    def reset(self):
        '''
        Method to reset the occupancy of all traps to 0
        '''

        self.states[:]   = 0.
        

    def get_acc_charge(self,t_after):
        '''
        Method to measure the accumulated charge at times t_after
        after the end of the ramp. It assumes the diode
        is no longer being exposed to light
        '''

        charge = np.zeros(len(t_after))

        usest   = copy.deepcopy(self.states[self.states])
        usert   = copy.deepcopy(self.t_rel[self.states])

        for i in range(len(t_after)):
            exp       = np.exp( -1*(t_after[i]) / usert )
            check     = np.random.random_sample(size=exp.size)
            idis      = check > exp 
            ikeep     = check < exp
            
            charge[i] = charge[i-1] + np.sum(idis)
            usest     = usest[ikeep]
            usert     = usert[ikeep]

        return charge
