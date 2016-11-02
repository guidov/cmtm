
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

#%function [cu]=cohbias(v,cb);
#%
#%Corrects for the bias inherent to coherence estimates.  Note the Matlab
#%function cohere.m returns squared-coherence, and the square-root should
#%be used.  Coherence below the minimum expected value returns a zero. 
#%
#%Requires the file cohbias.mat.  If the file does not exist,
#%prompts whether it should be created -- note the calculation
#%takes roughly an hour on a 2 GHz machine (i.e. it should be 
#%easier to get the file from http://web.mit.edu/~phuybers/www/XCM/index.html.)    
#%
#%
#%inputs:   v   - degrees of freedom, single value or vector (2 <= n <= 50)
#%          cb - biased coherence estimate, single value of vector (0 <= c <= 1).
#%          
#%outputs:  cu  - unbiased cohernce estimate (always less than cb).
#%
#%
#%Peter Huybers
#%phuybes@mit.edu
#%MIT, 2003.
def cohbias(v, cb):

    # Local Variables: A, c, f, i1, cb, i3, i2, M, ec, n, ct, pl, expect, v, z, k, cu, qans
    # Function calls: disp, isnan, repmat, cohbias, interp1, sum, sqrt, nargin, length, strncmpi, input, exist, pi, find, gamma
    if nargin<2.:
        help(cohbias)
        return []
    
    
    if v<2.:
        np.disp('Warning: degress of freedom must be greater or equal to two.')
        return []
    
    
    if np.logical_or(cb<0., cb > 1.):
        np.disp('Warning: biased coherence should be between zero and one, inclusive.')
        return []
    
    
    if v > 50.:
        np.disp('using 50 degrees of freedom')
        v = 50.
    
    
    #%if nargin==0; help cohbias.m; return; end;
    if exist('cohbias.mat') == 0.:
        #%Cohbias.mat file should be down-loaded with cmtm.m and cohbias.m
    #%The routine is included primarily to show how it was created.
    n = np.arange(2., 51.0, 1.)
    c = np.arange(.1, 1.0001, .0001)
    np.disp('-- The file cohbias.mat does not exist within the path.')
    qans = input('-- To create this file now enter \'y\' or to skip \'n\'. \n--> ', 's')
    qans
    if strncmpi(qans, 'y', 1.):
        z = np.arange(0., 1.1, .1)
        for i3 in np.arange(1., (length(n))+1):
            np.disp(n[int(i3)-1])
            for i2 in np.arange(1., (length(c)-1.)+1):
                for i1 in np.arange(1., (length(z))+1):
                    A[0] = 1.
                    #%Calculated according to: Amos and Koopmans, "Tables of the distribution of the
                    #%coefficient of coherence for stationary bivariate Gaussian processes", Sandia
                    #%Corporation, 1963
                    #%
                    #%Also see the manuscript of Wunsch, C. "Time-Series Analysis.  A Heuristic Primer".  
                    for k in np.arange(1., (n[int(i3)-1]-1.)+1):
                        A[int((k+1.))-1] = np.dot(matdiv(np.dot(np.dot(A[int(k)-1], n[int(i3)-1]-k), 2.*k-1.), np.dot(2.*n[int(i3)-1]-2.*k+1., k)), matdiv(1.-np.dot(c[int(i2)-1], z[int(i1)-1]), 1.+np.dot(c[int(i2)-1], z[int(i1)-1]))**2.)
                        
                    f[int(i1)-1] = np.dot(matdiv(np.dot(matdiv(np.dot(np.dot(np.dot(2.*(n[int(i3)-1]-1.), matixpower(1.-c[int(i2)-1]**2., n[int(i3)-1])), z[int(i1)-1]), matixpower(1.-z[int(i1)-1]**2., n[int(i3)-1]-2.)), np.dot(1.+np.dot(c[int(i2)-1], z[int(i1)-1]), matixpower(1.-np.dot(c[int(i2)-1], z[int(i1)-1]), 2.*n[int(i3)-1]-1.))), plt.gamma((n[int(i3)-1]-.5))), np.dot(np.sqrt(np.pi), plt.gamma(n[int(i3)-1]))), np.sum(A))
                    
                #%Use a quadratic Newton-Cotes methods to determine the cumulative sum
                for i1 in np.arange(2., (length(f)/2.)+1):
                    M[int(i1)-1] = np.dot(np.array(np.hstack((f[int((2.*(i0.)+1.))-1], np.dot(+4., f[int((2.*i1))-1]), +f[int((2.*i1+1.))-1]))), z[int((2.*i1))-1])
                    
                expect[int(i3)-1,int(i2)-1] = matdiv(np.sum(np.array(np.hstack((M, 1.)))), 6.*length(M))
                
            expect[int(i3)-1,int((i2+1.))-1] = 1.
            
        #% save cohbias.mat expect n c; 
    else:
        #%if skip cohbias.mat calculation
        expect = matcompat.repmat(c, length(n), 1.)
        
    
    #%stop qans condition
    else:
        #%if cohbias.mat already exists
        #%load cohbias.mat expect n c;
        
        
    
    #%stop cohbias calculation
    cb = cb.flatten(1)
    c = c.flatten(1)
    n = n.flatten(1)
    v = v.flatten(1)
    for ct in np.arange(1., (length(c))+1):
        ec[int(ct)-1] = interp1(n, expect[:,int(ct)-1], v)
        
    for ct in np.arange(1., (length(cb))+1):
        cu[int(ct)-1] = interp1(ec, c, cb[int(ct)-1])
        
    cu = cu.flatten(1)
    pl = nonzero(np.logical_and(np.logical_and(np.isnan(cu) == 1., cb<1.), cb >= 0.))
    #%If cu is NaN while cb is between (0,1)
    cu[int(pl)-1] = 0.
    return [cu]