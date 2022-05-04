% ykw3_dtstar.m - Measures dt* for SKS phases recorded at Yellowknife array
% station YKW3.
fref=40; 
% Read sdb, parse list of events to loop over and measure
data = readtable('~/Projects/DeepMelt/CanadianShield/YKW3/YKW3_sks_splitting.sdb', 'FileType','delimitedtext');

% read traces
trE = msac_read('data/YKW3_splitex.BHE');
trN = msac_read('data/YKW3_splitex.BHN');
trZ = msac_read('data/YKW3_splitex.BHZ');
wbeg = trE.a;
wend = trE.f;
% Split example is YKW3_2001152_203657
fast = 44;
tlag = 0.95;
spol = 94.3;
% Rotate to fast direction to get fast/slow traces
[trF,trS]=msac_rotate(trN,trE,fast) ;
tstar_search = 0:0.05:0.5;
dtstar_meas = msac_measure_dtstar(trF,trS,trS.a,trS.f,fref, tstar_search) ;

fprintf('Measured delta t* is %4.3f \n',dtstar_meas)

% search over potential fast directions and measure dtstar

fast_search = [-90:1:90];

ifrF = zeros(1,length(fast_search));
ifrS = zeros(1,length(fast_search));
dtstar = zeros(1,length(fast_search));

for i=1:length(fast_search)
   [trF,trS]=msac_rotate(trN,trE,fast_search(i)) ;
   trS = msac_tshift(trS,-tlag,'int') ;
   ifrF(i)=msac_ifa_wwind(trF,trF.a,trF.f) ;
   ifrS(i)=msac_ifa_wwind(trS,trS.a,trS.f) ;
   dtstar(i) = msac_measure_dtstar(trF,trS,trS.a,trS.f,fref, tstar_search) ;
end

figure("Position",[0 0 600 2000]);
plot(fast_search,dtstar,'b-','LineWidth',1.5) ;
hold on 
aax=axis();
xline(fast,'k-','LineWidth',1.5) ;
xline(spol,'k--','LineWidth',1.5) ;
xlabel('Reference frame rotation (degrees)')
ylabel('dt* (s)')
title(['tlag = 1.5; dt* =' sprintf('%5.2f',dtstar_meas)])

% Search over potential tlag 