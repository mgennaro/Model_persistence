import numpy as np

def rrgen(rate,tint,nint,ttot,deltatint=10,sat=80000,bias=1000,tplateau=0):

    '''
    Generate a raw counts ramp as function of time
    This function is derived from ksl gen_ramps function
    However the generated ramp has one "time step"
    per each unit count increment at ramp-up while
    it has only the necessary time steps to account for
    reset, or stop of the illumination

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

    #Seconds to get a dcount = 1
    inv_rate  = 1./rate

    #Seconds to saturation
    sectosat = (sat-bias)*inv_rate

    ramps = 0
    time  = np.empty( shape=(0) )
    cts   = np.empty( shape=(0) )
    
    while(ramps < nint):

        if (ramps == 0):
            offset = 0
        else:
            offset = time[-1]

        if (sectosat >= tint):
            timeh = offset + np.arange(0,tint,inv_rate)
            ctsh  = bias+np.arange(0,timeh.shape[0])
        else:
            timeh = offset + np.arange(0,sectosat,inv_rate)
            ctsh  = bias+np.arange(0,timeh.shape[0])
            timeh = np.append(timeh,offset+tint-1.e-5)
            ctsh  = np.append(ctsh,sat)

        if(tplateau > 0):
            timeh = np.append(timeh,timeh[-1]+tplateau)
            ctsh  = np.append(ctsh,ctsh[-1])


        timeh = np.append(timeh,offset+tint)
        ctsh  = np.append(ctsh,bias)
        timeh = np.append(timeh,timeh[-1]+deltatint)
        ctsh  = np.append(ctsh,bias)


        time = np.append(time,timeh)
        cts  = np.append(cts,ctsh)

        ramps += 1

    time = np.append(time,ttot)
    cts  = np.append(cts,bias)


    return time,cts
