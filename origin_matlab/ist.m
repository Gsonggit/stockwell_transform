function [ts] = inverse_st(st) 
% Returns the inverse of the Stockwell Transform of the REAL_VALUED timeseries. 
% Code by Robert Glenn Stockwell. 
% Reference is "Localization of the Complex Spectrum: The S Transform" 
% from IEEE Transactions on Signal Processing, vol. 44., number 4, April 1996, pages 998-1001. 
% 
%-------Inputs Needed------------------------------------------------ 
%   
%  S-Transform matrix created by the st.m function 
%-------Optional Inputs ------------------------------------------------ 
% none 
%-------Outputs Returned------------------------------------------------ 
% 
% ts     -a REAL-VALUED time series 
%--------Additional details----------------------- 
%   Copyright (c) by Bob Stockwell 
%   $Revision: 1.0 $  $Date: 2004/10/10  $ 
 
% sum over time to create the FFT spectrum for the positive frequencies 
stspe = sum(st,2); 
 
% get st matrix dimensions  
[nfreq,ntimes] = size(st); 
if rem(ntimes ,2) ~= 0 
    % odd number of points, so nyquist point is not aliased, so concatenate 
    % the reversed spectrum to create the negative frequencies 
    % drop the DC value 
    negspe = fliplr(stspe(2:nfreq)'); 
else 
     % even number of points 
     % therefore drop the first point (DC) and the last point (aliased nyqusit freq) 
     negspe = fliplr(stspe(2:nfreq-1)'); 
end 
 
% using symmetry of FFT spectrum of a real signal, recreate the negative frequencies from the positie frequencies 
 fullstspe = [conj(stspe')  negspe];     
 
% the time series is the inverse fft of this 
ts = ifft(fullstspe); 
% and take the real part, the imaginary part will be zero. 
ts = real(ts); 