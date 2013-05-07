function    predgaussian = getgaussian(pixels, xpeak, width, amp)
%  predgaussian = getgaussian(pixels,xpeak, width,amp); 
%
%  predgaussian = amp.*exp(-0.5 * ((pixels-xpeak)./width).^2 );
%
predgaussian = amp.*exp(-0.5 * ((pixels-xpeak)./width).^2 );
