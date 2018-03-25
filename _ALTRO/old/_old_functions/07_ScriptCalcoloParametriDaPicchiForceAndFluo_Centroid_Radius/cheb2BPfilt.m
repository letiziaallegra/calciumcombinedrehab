function sigfilt = cheb2BPfilt(sig,f_low,f_high,order,Fs,ripple)
% cheb2BPfilt band-pass filters the signal 'sig' between frequencies ['f_low' 'f_high']
% It uses a chebyshev type 2 filter of order 'order'
% 'Fs' is the sampling frequency of the signal 'sig'

if(nargin<6)
    sigLPfilt = cheb2LPfilt(sig,f_high,order,Fs);
    sigfilt = cheb2HPfilt(sigLPfilt,f_low,order,Fs);
else
    sigLPfilt = cheb2LPfilt(sig,f_high,order,Fs,ripple);
    sigfilt = cheb2HPfilt(sigLPfilt,f_low,order,Fs,ripple);
end