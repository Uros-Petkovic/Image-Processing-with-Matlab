function H = lpfilter1(type, M, N, D0, n)
%LPFILTER Computes frequency domain lowpass filters
%   H = LPFILTER1(TYPE, M, N, D0, n) creates the transfer function of
%   a lowpass filter, H, of the specified TYPE and size (M-by-N).  To
%   view the filter as an image or mesh plot, it should be centered
%   using H = fftshift(H).
%
%   Valid values for TYPE, D0, and n are:
%
%   'ideal'    Ideal lowpass filter with cutoff frequency D0.  n need
%              not be supplied.  D0 must be positive
%
%   'btw'      Butterworth lowpass filter of order n, and cutoff D0.
%              The default value for n is 1.0.  D0 must be positive.
%
%   'gaussian' Gaussian lowpass filter with cutoff (standard deviation)
%              D0.  n need not be supplied.  D0 must be positive.

% Use function dftuv to set up the meshgrid arrays needed for 
% computing the required distances.
[U, V] = dftuv(M, N);

% Compute the distances D(U, V). Elipse
D = sqrt(2*U.^2 + 0.6*V.^2);

% Begin fiter computations.
switch type
case 'ideal'
   H = double(D <=D0);
case 'btw'
   if nargin == 4
      n = 1;
   end
   H = 1./(1 + (D./D0).^(2*n));
case 'gaussian'
   H = exp(-(D.^2)./(2*(D0^2)));
case 'trapez'
   if nargin == 4
       D1 = 1.2*D0;
   else
       D1 = n;
   end
   H = (D1 - D)./(D1 - D0);
   H(H<0) = 0;
   H(H>1) = 1;
case 'exp'
   if nargin == 4
      n = 1;
   end
   H = exp(-(D./D0).^n);
otherwise
   error('Unknown filter type.')
end