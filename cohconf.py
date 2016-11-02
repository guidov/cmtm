
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

#%function [cl]=cohconf(v,level,unbias,c);
#%
#%inputs:   v     - degrees of freedom
#%          level - confidence level, i.e. 0.95 gives the 95% c.l. (.95=default).
#%          bias  - are coherence estimates bias corrected (use cohbias.m), 0=no (default), 1=yes.
#%          c     - true coherence, 0 <= c < 1 (0=default).
#%                    
#%
#%outputs:  
#%          cl     - coherence value for selected confidence level 
#%
#%
#%
#%Peter Huybers
#%phuybers@mit.edu
#%MIT, 2003
def cohconf(n, level, unbias, c):

    # Local Variables: A, c, f, level, i1, k, cl, n, F, unbias, z, Fz, pl
    # Function calls: disp, cohconf, cohbias, fliplr, interp1, sum, sqrt, nargin, length, diff, pi, cumsum, find, gamma
    if nargin<1.:
        help(cohconf)
        return []
    
    
    if nargin<2.:
        level = .95
    
    
    if nargin<3.:
        unbias = 0.
    
    
    if nargin<4.:
        c = 0.
    
    
    if n<2.:
        np.disp('Warning: degress of freedom must be greater or equal to two.')
        return []
    
    
    if np.logical_or(level<=0., level >= 1.):
        np.disp('Warning: confidence level should be between zero and one.')
        return []
    
    
    #%Calculated according to: Amos and Koopmans, "Tables of the distribution of the
    #%coefficient of coherence for stationary bivariate Gaussian processes", Sandia
    #%Corporation, 1963
    #%Also see Priestly, 1981
    z = np.arange(0., 1.0005, .0005)
    for i1 in np.arange(1., (length(z))+1):
        A[0] = 1.
        for k in np.arange(1., (n-1.)+1):
            A[int((k+1.))-1] = np.dot(matdiv(np.dot(np.dot(A[int(k)-1], n-k), 2.*k-1.), np.dot(2.*n-2.*k+1., k)), matdiv(1.-np.dot(c, z[int(i1)-1]), 1.+np.dot(c, z[int(i1)-1]))**2.)
            
        f[int(i1)-1] = np.dot(matdiv(np.dot(matdiv(np.dot(np.dot(np.dot(2.*(n-1.), matixpower(1.-c**2., n)), z[int(i1)-1]), matixpower(1.-z[int(i1)-1]**2., n-2.)), np.dot(1.+np.dot(c, z[int(i1)-1]), matixpower(1.-np.dot(c, z[int(i1)-1]), 2.*n-1.))), plt.gamma((n-.5))), np.dot(np.sqrt(np.pi), plt.gamma(n))), np.sum(A))
        
    #%Use a quadratic Newton-Cotes methods to determine the cumulative sum
    for i1 in np.arange(2., (length(f)/2.)+1):
        F[int(i1)-1] = f[int((2.*(i0.)+1.))-1]+4.*f[int((2.*i1))-1]+f[int((2.*i1+1.))-1]
        
    F = matdiv(F, 6.*length(F))
    F = np.array(np.hstack((np.fliplr((1.-matdiv(np.cumsum(np.fliplr(F)), np.sum(F)))), 1.)))
    Fz = np.array(np.hstack((z[0:0-2.:2.], 1.)))
    pl = nonzero((np.diff(F) > 0.))
    pl = np.array(np.hstack((1., pl+1.)))
    cl = interp1(F[int(pl)-1], Fz[int(pl)-1], level)
    if unbias == 1.:
        cl = cohbias(n, cl)
    
    
    return [cl]