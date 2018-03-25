function sigfilt = cheb2LPfilt(sig,f,order,Fs,ripple)
% cheb2BPfilt band-pass filters the signal 'sig' between frequencies ['f_low' 'f_high']
% It uses a chebyshev type 2 filter of order 'order'
% 'Fs' is the sampling frequency of the signal 'sig'
if(nargin<5)
    ripple = 20;
end

Wn = f/(Fs/2); % Normalized stopband edge frequency
[B, A] = cheby2(order,ripple,Wn); 
sigfilt = filtfilt(B, A, sig);