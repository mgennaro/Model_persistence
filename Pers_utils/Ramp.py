import numpy as np
import matplotlib.pyplot as plt

class Ramp(object):

    '''
    Class that describes an observation (ramp) in
    raw counts vs. time.

    :rate:
        count rate of the source (counts/sec)

    :tint:
        integration time for a single exposure (sec)
        (i.e. time between the start of the integration
        and the diode reset)

    :nit:
        number of exposures

    :ttot:
        total time to simulate (sec)

    :deltatint ( optional ):
        time between the end of an exposure and the start
        of the next

    :sat ( optional ):
        saturation level in the pixel

    :bias ( optional ):
        counts at reset

    :tplateau ( optional ):
        time between the end of the integration and the reset
        (useful for some cal observations with special command)
    
    '''

    def __init__(self, rate,tint,nint,ttot,deltatint=10,sat=80000,bias=1000,tplateau=0):

        self.rate = rate
        self.tint = tint
        self.nint = nint
        self.ttot = ttot
        self.deltatint  = deltatint
        self.sat = sat
        self.bias = bias
        self.tplateau = tplateau

        self.rrgen()

    def rrgen(self):

        '''
        Generate a raw counts ramp as function of time
        This function is derived from ksl gen_ramps function
        However the generated ramp has one "time step"
        per each unit count increment at ramp-up while
        it has only the necessary time steps to account for
        reset, or stop of the illumination
        '''
        '''
        Seconds to get a dcount = 1
        '''
        inv_rate  = 1./self.rate

        '''
        Seconds to saturation
        '''
        sectosat = (self.sat-self.bias)*inv_rate

        self.rtime  = np.empty( shape=(0) )
        self.rcts   = np.empty( shape=(0) )

        ramps = 0
        while(ramps < self.nint):

            if (ramps == 0):
                offset = 0
            else:
                offset = self.rtime[-1]

            if (sectosat >= self.tint):
                timeh = offset + np.arange(0,self.tint,inv_rate)
                ctsh  = self.bias+np.arange(0,timeh.shape[0])
            else:
                timeh = offset + np.arange(0,sectosat,inv_rate)
                ctsh  = self.bias+np.arange(0,timeh.shape[0])
                timeh = np.append(timeh,offset+self.tint-1.e-5)
                ctsh  = np.append(ctsh,self.sat)


            timeh = np.append(timeh,offset+self.tint)
            ctsh  = np.append(ctsh,self.bias)
            timeh = np.append(timeh,timeh[-1]+self.deltatint)
            ctsh  = np.append(ctsh,self.bias)

            self.rtime = np.append(self.rtime,timeh)
            self.rcts  = np.append(self.rcts,ctsh)

            ramps += 1

        if(self.tplateau > 0):
            self.rtime = np.append(self.rtime,self.rtime[-1]+self.tplateau)
            self.rcts  = np.append(self.rcts,self.rcts[-1])

        self.rtime = np.append(self.rtime,self.ttot)
        self.rcts  = np.append(self.rcts,self.bias)


    def test_plot(self,**kwargs):
        '''
        Method to visually check the ramp
        '''

        plt.plot(self.rtime,self.rcts,**kwargs)
        plt.ylabel('Raw Counts in pixel')
        plt.xlabel('Time (s)')
        plt.ylim(0,1.05*self.sat)
