function st=stquick(h);
% get from http://www.cora.nwra.com/~stockwel/index.php?module=phpwsbb&PHPWSBB_MAN_OP=view&PHPWS_MAN_ITEMS[]=37
H=fft(h);
n=length(h);
st=ifft([[1; zeros(n-1,1)], ...
[exp(-2*pi^2*([0:ceil(n/2)-1 -floor(n/2):-1]'* ...
(1./[1:ceil(n/2)-1 -floor(n/2):-1])).^2)]] .* ...
H(fliplr(toeplitz([n 1:n-1],[n:-1:1])))).';
return;
