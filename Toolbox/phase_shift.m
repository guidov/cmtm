%Phase shift a record
%
%function [S2]=phase_shift(S,phi);
%
%S2 is the phase-shifted record
%S is the original record
%phi is the angle to shift by in radians

function [S2]=phase_shift(S,phi);

S=S(:);
shift=exp(-i*phi*ones(1,floor(length(S)/2)));
if rem(length(S),2)>0,
  shift=[0 conj(shift)          fliplr(shift)]';
else;
  shift=[0 conj(shift(1:end-1)) fliplr(shift)]';
end;
S2=real(ifft(fft(S).*shift));


