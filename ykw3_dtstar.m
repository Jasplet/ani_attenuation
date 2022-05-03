% ykw3_dtstar.m - Measures dt* for SKS phases recorded at Yellowknife array
% station YKW3.
fref=40; 
% Read sdb, parse list of events to loop over and measure
data = readtable('~/Projects/DeepMelt/CanadianShield/YKW3/YKW3_sks_splitting.sdb', 'FileType','delimitedtext');

% read traces
trE = msac_read('YKW3_splitex.BHE');
trN = msac_read('YKW3_splitex.BHN');
trZ = msac_read('YKW3_splitex.BHZ');
wbeg = trE.a;
wend = trE.f;
fast = 44;
tlag = 0.95;
% Rotate to fast direction to get fast/slow traces
[trF,trS]=msac_rotate(trN,trE,fast) ;
tstar_search = 0:0.01:2;
dtstar = msac_measure_dtstar(trF,trS,trS.a,trS.f,fref, tstar_search) ;

disp(dtstar)