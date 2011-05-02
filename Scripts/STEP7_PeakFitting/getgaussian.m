function    predlorentz = getgaussian(pixels,xpeak, width,amp); 

predlorentz = amp.*exp(-0.5 * ((pixels-xpeak)./width).^2 );
