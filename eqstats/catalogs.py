import numpy as np

def omori_times(ncat, nevents, tmin, tmax, b, p=1., detectprob = None):
    """
    creates ncat synthetic realizations of an Omori decay in seismicity

    parameters are nevents (number of events), tmin (minimum catalog time, main shock is t=0)
    tmax (maximum catalog time), b (Omori time offset, R \propto 1/(b+t)^p)

    Inputs:
    ncat = number of realizations
    nevents = number of events per realization
    tmin = catalog start time (t=0 is main shock)
    tmax = catalog end time (t=0 is main shock)
    b, p = Omori parameters
    detectprob = function mapping event time to detection probability

    returns numpy array with shape (ncat, nevents)
    """

    assert(ncat > 0)
    assert(nevents > 0)
    assert(tmin > 0.)
    assert(tmax > tmin)
    assert(b > 0.)
    assert(p > 0.)

    if detectprob is None:
        detectprob = lambda x: 1.

    acceptedtimes = []

    for i in range(nevents*ncat):

        while True:

            times = np.random.random()

            if p == 1.:
                times = tmin + (b+tmin)*(((b+tmax)/(b+tmin))**times - 1.)
            else:
                times = -b + ((1.-times)/(b+tmin)**(p-1.)+times/(b+tmax)**(p-1.))**(-1./(p-1.))

            detect = detectprob(times)

            if detect >= np.random.random():
                acceptedtimes.append(times)
                break

    times = np.reshape(np.array(acceptedtimes), (ncat, nevents))

    times = np.sort(times)

    return times

def random_times(nevents, tmin = 0., tmax = 100.):
    "generates a random sequence of nevents events"

    times = tmin + (tmax-tmin)*np.random.random(nevents)
    times = np.sort(times)

    return times
