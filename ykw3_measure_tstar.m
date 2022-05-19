%% Script to measure t* for YKW3 data.
close all ; clear
ykw3_data = readtable('YKW3_sks_splitting.txt', 'Format','%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
% have changed.sdb to .txt to make matlab reasing easier.
path = '/Users/ja17375/Projects/DeepMelt/CanadianShield/YKW3/run/';
rows = height(ykw3_data);
fref = 10; 
ykw3_dtstars = zeros(length(rows),1);
for row = 1:rows
   snr = ykw3_data.SNR(row);
   if snr >= 15
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
       
       %figure()
       
       ifrF=msac_ifa_wwind(trF, wbeg, wend) ;
       ifrS=msac_ifa_wwind(trS, wbeg, wend) ;

       difr = ifrF - ifrS;
       
       [difrs, dtstars] =  msac_measure_dtstar_iter(trF, trS, wbeg, wend, fref);
       ykw3_dtstars(row) = dtstars(end);
       srange = [0:0.1:4];
       difr_range = zeros(1,length(srange));
       if sign(difr) < 0
              % fast is more attenuated (has a smaller Ifr) than slow wave
              % so we want to measure dt* for the slow wave (as it is the
              % less attenuated
            tr_attn = trS;
            ifr_ref = ifrF;
          elseif sign(difr) > 0 
            tr_attn= trF;
            ifr_ref = ifrS;
          end
       for i = 1:length(srange)
          tr_attenuated = msac_apply_tstar_operator(tr_attn, fref, srange(i));
          ifr_trial = msac_ifa_wwind(tr_attenuated, wbeg, wend);
          difr = ifr_trial - ifr_ref;
          difr_range(i) = difr;
       end
   else
       fprintf('SNR = %4.2f which is less than threshold \n',snr)
       continue
   end

% Plot iterative search for t* and search range
   figure() 
   plot(srange, difr_range, 'k-', 'Linewidth',1.5)
   hold on
   % Plot absolute value of dt* here. This is becuase this figure (dIfr v
   % dt* is in the reference frame of the more attenuated trace (i.e dt* is
   % the t* operator added to the less attenuated trace to match the
   % instantaneous frequencies.
   % Only when we are thinking in terms of fast and slow waves does the
   % sign of dt* matter (as a negative dt* means fast is more attenuated
   % and vice-versa).
   plot(abs(dtstars), difrs, 'rx', 'Markersize',10)
   xlabel('\Delta t*')
   ylabel('\Delta IFr')
   xline(abs(ykw3_dtstars(row)), 'k--', 'LineWidth',1.5)
   yline(0, 'k:')
   title(['YKW3 ', date, '-',time])
end