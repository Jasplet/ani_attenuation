%% Script to measure t* for YKW3 data.
close all ; clear
ykw3_data = readtable('YKW3_sks_splitting.txt', 'Format','%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
% have changed.sdb to .txt to make matlab reasing easier.
path = '/Users/ja17375/Projects/DeepMelt/CanadianShield/YKW3/run/';
rows = height(ykw3_data);
fref = 10; 
rows = [1,2,3,4];
ykw3_dtstars = zeros(length(rows),1);
for row = 1:rows
   date = char(ykw3_data.DATE(row));
   time = char(ykw3_data.TIME(row));
   fast_true = ykw3_data.FAST(row);
   wbeg = ykw3_data.WBEG(row);
   wend = ykw3_data.WEND(row);
   fileid = sprintf('YKW3_%s_%s',date, time);
   
   trN = msac_read([path, fileid, '.BHN']);
   trE = msac_read([path, fileid, '.BHE']);
   trZ = msac_read([path, fileid, '.BHZ']);
   
   %Rotate to fast direction
   [trF,trS] = msac_rotate(trN,trE,fast_true) ;

   ifrF=msac_ifa_wwind(trF, wbeg, wend) ;
   ifrS=msac_ifa_wwind(trS, wbeg, wend) ;
    
   difr = ifrF - ifrS;

   [difrs, dtstars] =  msac_measure_dtstar_iter(trF, trS, wbeg, wend, fref);
   ykw3_dtstars(row) = dtstars(end);
   srange = [0:0.2:4];
   difr_range = zeros(1,length(srange));
   for i = 1:length(srange)
      tr_attenuated = msac_apply_tstar_operator(trF, fref, srange(i));
      ifr_trial = msac_ifa_wwind(tr_attenuated, tr_attenuated.a, tr_attenuated.f);
      difr = ifr_trial - ifrS;
      difr_range(i) = difr;
   end

% Plot iterative search for t* and search range
   figure() 
   plot(srange, difr_range, 'k-', 'Linewidth',1.5)
   hold on
   plot(dtstars, difrs, 'rx', 'Markersize',10)
   xlabel('\Delta t*')
   ylabel('\Delta IFr')
   title(['YKW3 ', date, '-',time])
end