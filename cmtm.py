import numpy as np
import scipy
from spectrum import dpss
import cohbias

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass


#
# Peter Huybers
# MIT, 2003
# phuyber@mit.edu

# Python modifications:
# Guido Vettoretti, 2016
# g.vettoretti@utoronto.ca

# Dependancies: numpy, scipy, spectrum (dpss)
# Spectrum: http://thomas-cokelaer.info/software/spectrum/html/contents.html
#           conda install -c fklab spectrum
#           pip install spectrum
#           easy_install -U spectrum
#
# 0.1: First conversion of cmtm.m using libermate
#      (https://github.com/awesomebytes/libermate)
# 0.2: Clean up code (PEP-8)

def cmtm(x, y, dt=1.0, NW=8, qbias=0.0, confn=0.0, qplot=True):
    """
    s, c, ph, ci, phi = cmtm(x,y,dt,NW,qbias,confn,qplot)
    Multi-taper method coherence using adaptive weighting and correcting
    for the bias inherent to coherence estimates.  The 95% coherence
    confidence level is computed by cohconf.py.  In addition, a
    Monte Carlo estimation procedure is available to estimate phase 95%
    confidence limits.

     Args:
             x     - Input data vector 1.
             y     - Input data vector 2.
             dt    - Sampling interval (default 8)
             NW    - Number of windows to use (default 8)
             qbias - Correct coherence estimate for bias (default 0).
             confn - Number of iterations to use in estimating phase
                           uncertainty using a Monte Carlo method. (default 0)
             qplot - Plot the results.  The upper tickmarks indicate the
                     bandwidth of the coherence and phase estimates.

     Returns:
             s       - frequency
             c       - coherence
             ph      - phase
             ci      - 95% coherence confidence level
             phi     - 95% phase confidence interval, bias corrected
                       (add and subtract phi from ph).


    required: cohconf.py, cohbias.py, cohbias.nc, scipy signal processing.
    """
    # Local Variables: ci, fx, cb, qplot, Pk, E, qbias, vari, ys, phut, Fx, Fy,
    # phl, ds, fkx, fky, tol, Ptemp, ph, phlt, pl, NW, phi, P1, pls, xs, i1,fy,
    # wk, N, P, V, dt, confn, phu, a, c, b, Cxy, Pkx, Pky, iter, col, s, w, v,
    # y, x, h, k
    # Function calls: disp, cmtm, dpss, cohconf, conv, fill, fft, set, conj,
    #  repmat, find, size, plot, angle, figure, cohbias, min, axis, sum, si,
    #  sqrt, abs, zeros, rem, xlabel, pi, ciph, real, max, ylabel, sort,
    # nargin, ones, randn, subplot, ifft, clf, gcf, fliplr, length, num2str,
    # title, round, mean

    # pre-checks
    if NW < 1.5:
        raise ValueError("Warning: NW must be greater or equal to 1.5")

    print('Number of windows: ', NW)
    if qbias == 1.:
        print('Bias correction:   On')
    else:
        print('Bias correction:   Off')

    print('Confidence Itera.: ', confn)
    if qplot == 1.:
        print('Plotting:          On')
    else:
        print('Plotting:          Off')

    x = x.flatten(1)-np.mean(x)
    y = y.flatten(1)-np.mean(y)
    if x.shape[0] != y.shape[0]:
        raise ValueError('Warning: the lengths of x and y must be equal.')

    #  define some parameters
    N = x.shape[0]
    k = np.max(np.round((2.*NW)), N)
    k = np.max((k-1.), 1.)
    s = np.arange(0., (1./dt-1./np.dot(N, dt)) +
        (1./np.dot(N, dt)), 1./np.dot(N, dt)).conj().T
    pls = np.arange(2., ((N+1.)/2.+1.)+1)
    v = 2*NW-1  # approximate degrees of freedom
    if y.shape % 2 == 1:
        pls = pls[0:0-1.]

    # Compute the discrete prolate spheroidal sequences,
    # requires the spectral analysis toolbox.
    [E, V] = dpss(N, NW, k)
    # Compute the windowed DFTs.
    fkx = np.fft((E[:, 0:k]*x[:, int(np.ones(1., k))-1]), N)
    fky = np.fft((E[:, 0:k]*y[:, int(np.ones(1., k))-1]), N)
    Pkx = np.abs(fkx)**2.
    Pky = np.abs(fky)**2.
    # Iteration to determine adaptive weights:
    for i1 in np.arange(1, 3):
        if i1 == 1:
            vari = np.dot(x.conj().T, x)/N
            Pk = Pkx

        if i1 == 2:
            vari = np.dot(y.conj().T, y)/N
            Pk = Pky

        P = (Pk[:, 0]+Pk[:, 1])/2.
        # initial spectrum estimate
        Ptemp = np.zeros((N, 1.))
        P1 = np.zeros((N, 1.))
        tol = np.dot(.0005, vari)/N
        # usually within tolerance in about three iterations,
        # see equations from [2] (P&W pp 368-370).
        a = np.dot(vari, 1.-V)
        while np.sum(np.abs((P-P1))/N) > tol:
            b = np.dot(P, np.ones(1., k))/(np.dot(P, V.conj().T) +
                np.dot(np.ones(N, 1.), a.conj().T))
            # weights
            wk = b**2.*np.dot(np.ones(N, 1.), V.conj().T)
            # new spectral estimate
            P1 = (np.sum((wk.conj().T*Pk.conj().T))/
                np.sum(wk.conj().T)).conj().T
            Ptemp = P1
            P1 = P
            P = Ptemp
            # swap P and P1

        if i1 == 1:
            dotp1 = np.dot(np.sqrt(k), np.sqrt(wk))
            fkxtmp = np.sum(np.sqrt(wk.conj().T)).conj().T
            # fkx = dotp1*fkx/matcompat.repmat(fkxtmp, 1., k)
            fkx = dotp1*fkx/np.kron(np.ones((1, k)), fkxtmp)
            Fx = P
            # Power density spectral estimate of x

        if i1 == 2:
            dotp1 = np.dot(np.sqrt(k), np.sqrt(wk))
            fkytmp = np.sum(np.sqrt(wk.conj().T)).conj().T
            # fky = dotp1*fky/matcompat.repmat(fkytmp, 1., k)
            fky = dotp1*fky/np.kron(np.ones((1, k)), fkytmp)
            Fy = P
            # Power density spectral estimate of y

    # As a check, the quantity sum(abs(fkx(pls,:))'.^2) is the same as Fx and
    # the spectral estimate from pmtmPH.
    # Compute coherence
    Cxy = np.sum(np.array(np.hstack((fkx*np.conj(fky)))).conj().T)
    ph = np.divide(np.angle(Cxy)*180., np.pi)
    c = np.abs(Cxy)/np.sqrt((np.sum((np.abs(fkx.conj().T)**2.)) *
        np.sum((np.abs(fky.conj().T)**2.))))
    # correct for the bias of the estimate
    if qbias == 1:
        c = cohbias(v, c).conj().T


    # Phase uncertainty estimates via Monte Carlo analysis.
    if confn > 1:
        cb = cohbias(v, c).conj().T
        nlist = np.arange(1., (confn)+1)
        ciph = np.zeros((nlist, x.shape[0])) # not sure about the cmtm return length
        phi = np.zeros((nlist, x.shape[0])) # not sure about the cmtm return length
        for iter in nlist:
            if plt.rem(iter, 10.) == 0.:
                print('phase confidence iteration: ', iter)

            fx = np.fft(np.randn(x.shape)+1.)
            fx = np.divide(fx, np.sum(np.abs(fx)))
            fy = np.fft(np.randn(y.shape)+1.)
            fy = np.divide(fy, np.sum(np.abs(fy)))
            ys = np.real(np.ifft((fy*np.sqrt((1.-cb.conj().T**2.)))))
            ys = ys+np.real(np.ifft((fx*cb.conj().T)))
            xs = np.real(np.ifft(fx))

        si, ciph[iter, :], phi[iter, :] = cmtm(xs, ys, dt, NW)
        pl = np.round(np.dot(.975, iter))

        # sorting and averaging to determine confidence levels.
        phi = np.sort(phi)
        phi = np.array(np.vstack((np.hstack((phi[int(pl)-1, :])), np.hstack((-phi[int((iter-pl+1.))-1, :])))))
        phi = np.mean(phi)
        phi = plt.conv(phi[0:], (np.array(np.hstack((1., 1., 1.)))/3.))
        phi = phi[1:0-1.]
    else:
        phi = np.zeros(pls.shape[0])

    # Cut to one-sided funtions
    c = c[int(pls)-1]
    s = s[int(pls)-1].conj().T
    ph = ph[int(pls)-1]

    # Coherence confidence level
    ci = cohconf(v, .95)

    # not corrected for bias, this is conservative.
    ci = np.dot(ci, np.ones((c.shape[0])))

    # plotting
    if qplot:
        phl = ph-phi
        phu = ph+phi
        # coherence
        print('coherence plot')
        # phase
        print('phase plot')

    return s, c, ph, ci, phi
