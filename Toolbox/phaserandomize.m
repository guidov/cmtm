%Randomize the Fourier phases of a signal
%
%
%function [R]=PhaseRandomize(F);

function [R]=phaserandomize(F);

F=F(:);
rphase=exp(-i*rand(1,floor(length(F)/2))*2*pi);
if rem(length(F),2)>0,
  rphase=[0 conj(rphase)          fliplr(rphase)]';
else;
  rphase=[0 conj(rphase(1:end-1)) fliplr(rphase)]';
end;
R=real(ifft(fft(F).*rphase));


