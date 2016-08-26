import numpy as np
from scipy.stats import kstest, binom_test
from scipy.misc import factorial
from scipy.special import erf
from statsmodels.stats.diagnostic import acorr_ljungbox

def ksexp_test(recurtimes, ttot = 100.):
    "Calculates statistic for ks exponential test"
    
    return kstest(recurtimes*(len(recurtimes)+1)/ttot, 'expon')

def ksunif_test(cumtimes, ttot = 100.):
    "Calculates p-value for ks uniform test"

    dstat, ksvalunif = kstest(cumtimes/ttot,'uniform')

    return ksvalunif

def acorr_test(recurtimes):
    "Calculates the first lag autocorrelation of a catalog"
    teststat, pval = acorr_ljungbox(recurtimes, lags=1)
    return pval[0]

def var_test(recurtimes):
    "Calculates variance of the recurrence times of a catalog"
    return np.var(recurtimes)/np.mean(recurtimes)**2

def condchi_test(cumtimes, ttot = 100., nbins = 100):
    "Calculates conditional chi-squared statistic"

    bincounts, bins = np.histogram(cumtimes, nbins, (0., ttot))

    obsrate = float(len(cumtimes))/ttot

    condchi = 0.

    for k in range(nbins):
        condchi += (bincounts[k]-obsrate)**2/obsrate

    return condchi

def multichi_test(cumtimes, ttot = 100., nbins = 100):
    "Calculates multinomial chi squared statistic"

    def poisson_cum(rate, k, ttot):
        try:
            total = 0.
            for i in range(0,k+1):
                total += rate**i/factorial(i)
            total *= np.exp(-rate)*ttot
        except OverflowError:
            total = ttot*0.5*(1.+erf((k+0.5-rate)/np.sqrt(2.*rate)))
        return total

    def poisson(rate, k, ttot):
        try:
            p = ttot*np.exp(-rate)*rate**k/factorial(k)
        except OverflowError:
            p = ttot*0.5*(erf((k+0.5-rate)/np.sqrt(2.*rate))-erf((k-0.5-rate)/np.sqrt(2.*rate)))
        return p

    bincounts, bins = np.histogram(cumtimes, nbins, (0., ttot))

    obsrate = float(len(cumtimes))/ttot

    kmin = 0
    tmpmax = int(obsrate)
    ekmin = poisson_cum(obsrate, kmin, ttot)
    etmpmax = poisson_cum(obsrate,tmpmax, ttot)

    while tmpmax-kmin > 1:
        newk = (kmin+tmpmax)//2
        eknew = poisson_cum(obsrate, newk, ttot)
        if eknew >= 5:
            tmpmax = newk
            etmpmax = eknew
        else:
            kmin = newk
            ekmin = eknew

    kmin = 0
    tmpmax = int(obsrate)
    ekmin = poisson_cum(obsrate, kmin, ttot)
    etmpmax = poisson_cum(obsrate,tmpmax, ttot)

    while tmpmax-kmin > 1:
        newk = (kmin+tmpmax)//2
        eknew = poisson_cum(obsrate, newk, ttot)
        if eknew >= 5:
            tmpmax = newk
            etmpmax = eknew
        else:
            kmin = newk
            ekmin = eknew

    kmin = tmpmax
    ekmin = etmpmax

    kmax = int(5*obsrate)
    ekmax = ttot-poisson_cum(obsrate,kmax-1, ttot)
    tmpmin = int(obsrate)
    etmpmin = poisson_cum(obsrate,tmpmin, ttot)

    while kmax - tmpmin > 1:
        newk = (kmax+tmpmin)//2
        eknew = ttot-poisson_cum(obsrate, newk, ttot)
        if eknew >= 5:
            tmpmin = newk
            etmpmin = eknew
        else:
            kmax = newk
            ekmax = eknew

    kmax = tmpmin
    ekmin = etmpmin

    binfloor = bincounts[:]
    binfloor[bincounts < kmin] = kmin
    binfloor[bincounts > kmax] = kmax
    
    nbins = kmax-kmin+1
    ec = np.empty(nbins)
    oc = np.empty(nbins)
    ec[0] = ekmin
    oc[0] = np.sum(binfloor == kmin)
    ec[nbins-1] = ekmax
    oc[nbins-1] = np.sum(binfloor == kmax)
    for k in range(1,nbins-1):
        ec[k] = poisson(obsrate,k+kmin, ttot)
        oc[k] = np.sum(binfloor == kmin+k)

    multichi = 0.

    for k in range(nbins):
        multichi += (oc[k]-ec[k])**2/ec[k]

    return multichi

def brownzhao_test(cumtimes, ttot = 100., nbins = 100):
    "Calculates Brown and Zhao chi-squared statistic"

    bincounts, bins = np.histogram(cumtimes, nbins, (0., ttot))

    y = np.empty(int(ttot))
    yhat = 0.

    for k in range(int(ttot)):
        y[k] = np.sqrt(bincounts[k]+3./8.)
        yhat += y[k]/ttot

    brownzhao = 0.

    for k in range(int(ttot)):
        brownzhao += 4.*(y[k]-yhat)**2

    return brownzhao

def bigtrig_test(cumtimes, mag, ttot = 100., mbig = 2.5):
    "Calculates test for triggering after big earthquakes"

    n = len(cumtimes)
    obsrate = float(n)/ttot

    nbig = int(np.sum(mag >= mbig))

    tbig = cumtimes[mag >= mbig]

    count = np.zeros(n)

    for tb in tbig:
        count[(cumtimes-tb > 0.)*(cumtimes-tb <= 1.)] = 1.

    tw = 0.
    for i in range(nbig):
        if i == 0:
            tw += 1.
        elif tbig[i]-tbig[i-1] > 1.:
            tw += 1.
        else:
            tw += tbig[i]-tbig[i-1]

    count[mag >= mbig] = 0.

    ntrig = np.sum(count)

    return binom_test(int(ntrig), n-nbig, tw*obsrate/float(n))

def get_pval(teststat, randvals):
    "convert test statistic to p value given array of test statistic for random catalog"

    return ((-(teststat > randvals)).astype(float)).sum()/float(len(randvals))

def get_power(pvals, thresh):
    "turns p-values into detection power given detection threshold"

    return np.sum(pvals < thresh)/float(len(pvals))
