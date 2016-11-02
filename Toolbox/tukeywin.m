function w = tukeywin(n,r)
%TUKEYWIN Tukey window.
%   W = TUKEYWIN(N,R) returns an N-point Tukey window in a column vector.
%   The R parameter specifies the ratio of taper to constant sections.
%   This ratio is normalized to 1 (i.e., 0 < R < 1).  Note, at extreme
%   values of R, the Tukey window degenerates into other common windows.  
%   Thus when R = 1, it is equivalent to a Hanning window. Conversely,
%   for R = 0 the Tukey window assumes a constant value (i.e., boxcar).
%
%   EXAMPLE:
%      N = 64; 
%      w = tukeywin(N,0.5); 
%      plot(w); title('64-point Tukey window, Ratio = 0.5');
%
%   See also BARTLETT, BARTHANNWIN, BLACKMAN, BLACKMANHARRIS, BOHMANWIN, 
%            CHEBWIN, GAUSSWIN, HAMMING, HANN, KAISER, NUTTALLWIN, RECTWIN,
%            TRIANG, WINDOW.

%   Reference:
%     [1] fredric j. harris [sic], On the Use of Windows for Harmonic Analysis
%         with the Discrete Fourier Transform, Proceedings of the IEEE,
%         Vol. 66, No. 1, January 1978, Page 67, Equation 38.

%   Author(s): A. Dowd
%   Copyright 1988-2001 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2001/04/12 18:06:11 $

n=n+2; %get rid of zeros at the ends;

if r <= 0,
    w = ones(n,1);
elseif r >= 1,
    w = hann(n);
else
    t = linspace(0,1,n)';
    % Defines period of the taper as 1/2 period of a sine wave.
    per = r/2; 
    tl = floor(per*(n-1))+1;
    th = n-tl+1;
    % Window is defined in three sections: taper, constant, taper
    w = [ ((1+cos(pi/per*(t(1:tl) - per)))/2);  ones(th-tl-1,1); ((1+cos(pi/per*(t(th:end) - 1 + per)))/2)];
end

w=w(2:end-1);