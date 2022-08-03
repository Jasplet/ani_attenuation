% MSAC_SPLITWAVE - Create a split shear-wave (natively) 
function [varargout]=msac_splitwave(fast,tlag,varargin)

noise = 0.0 ;
dfreq = 0.2 ;
delta = 0.05 ; % sampling rate
spol  = 45 ;
wavelet = 'dgauss' ;

% default to 5 times dominant wave period (can override)
nsamp = 1 + (1 / dfreq)*5 / delta ;

% do some checking
if nargout<2 | nargout >3
   error('MSAC_SPLITWAVE2: 2 or 3 outputs required.')
end

%  ** process the optional arguments
iarg = 1 ;
while iarg <= (length(varargin))
   switch lower(varargin{iarg})
      case 'noise'
         noise = varargin{iarg+1} ;
         iarg = iarg + 2 ;
      case 'dfreq'
         dfreq = varargin{iarg+1} ;
         iarg = iarg + 2 ;
      case 'delta'
         delta = varargin{iarg+1} ;
         iarg = iarg + 2 ;
      case 'nsamp'
         nsamp = varargin{iarg+1} ;
         iarg = iarg + 2 ;         
      case 'spol'
         spol = varargin{iarg+1} ;
         iarg = iarg + 2 ;         
      case 'wavelet'
         wavelet = varargin{iarg+1} ;
         iarg = iarg + 2 ; 
      otherwise 
         error(['Unknown option: ' varargin{iarg}]) ;   
   end   
end

% create a timebase
t = ([0:nsamp-1]-floor(nsamp/2)).*delta ;

% create the base wavelet and normalise
switch lower(wavelet)
   case 'dgauss'
      y=dgauss_wavelet(t,dfreq) ;
   case 'gabor'
      y=gabor_wavelet(t,dfreq,6,2*pi/5,0) ;
   otherwise
      error('MSAC_SPLITWAVE: Unknown wavelet')
end

% Make the source trace.
trS = msac_new(y,delta,'b',t(1),'e',t(end),'evla',0,'evlo',0,'evdp',500,...
               'stla',80,'stlo',0,'kstnm','SYN','reftime',[3001 1 12 00 00 00]) ;

% Set analysis window (twice period in each direction, plus tlag)
trS.a = -2*(1 / dfreq) ; trS.user0 = trS.a ;
trS.f = 2*(1 / dfreq)+tlag ; trS.user2 = trS.f ;

% generate N and E components
trN = trS ; trN.cmpaz = 0; trN.cmpinc = 90 ; trN.kcmpnm = 'N' ;
trE = trS ; trE.cmpaz = 90; trE.cmpinc = 90 ; trE.kcmpnm = 'E' ;
if nargout==3
   trZ = trS ; trZ.cmpaz = 00; trZ.cmpinc = 0 ; trZ.kcmpnm = 'Z' ;
end

spol

trN.x1 = trS.x1 .* cosd(spol);
trE.x1 = trS.x1 .* sind(spol);
if nargout==3
   trZ.x1 = zeros(1,trZ.npts);
end

% add noise if required. 
MaxAmp = max(abs(trS.x1)) ;

trN.x1 = trN.x1 + randn([1 trN.npts]).*MaxAmp.*noise ;
trE.x1 = trE.x1 + randn([1 trE.npts]).*MaxAmp.*noise ;
if nargout==3
   trZ.x1 = trZ.x1  + randn([1 trZ.npts]).*MaxAmp.*noise ;
end

% now apply the splitting
[trF,trS]=msac_rotate(trN,trE,fast) ;
trS = msac_tshift(trS,tlag,'int') ;
[trN,trE]=msac_rotate(trF,trS,-fast) ;

% finally, construct the output
if nargout == 2
   varargout = {trN trE} ;
elseif nargout == 3
   varargout = {trN trE trZ} ;
end

end

function y=gabor_wavelet(t,f0,gamma,v,t0)
   y = 2*cos(2.*pi.*f0.*(t-t0)+v).*exp((-4*pi.^2.*f0.^2.*(t-t0).^2)/gamma) ;
   y = y./max(y) ;
end

function y=dgauss_wavelet(t,dfreq)
   n=length(t);
   % make reference wavelet ; dfreq = 0.2 Hz
   [yref,tref]=gauswavf(-5,5,1001) ;

   % stretch / compress the baseline to match required dfreq
   tref = tref.*(0.2/dfreq) ;

   % interp / extrapolate onto the timebase
   y = interp1(tref,yref,t,'pchip',0) ;
   y = y./max(abs(y)) ;
end

