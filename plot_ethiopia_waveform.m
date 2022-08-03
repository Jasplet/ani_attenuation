% Plot the waveforms for a station 
clc; clear; close all
stat = 'ODAS';
ifr_fast = -15;
ifr_dtstar = 1.55;
file = sprintf('/Users/ja17375/Projects/DeepMelt/Ethiopia/%s_sks_reprocessed.txt',stat);
data = readtable(file, 'Format','%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
% have changed.sdb to .txt to make matlab reasing easier.
path = sprintf('/Users/ja17375/Projects/DeepMelt/Ethiopia/data/%s/run/', stat);
snr_min = 5;
idx_test = find(data.SNR > snr_min);
%idx_test = 1:36;
fref=10;
n = length(idx_test);

for i=1:n
    idx = idx_test(i);
    date = char(data.DATE(idx));
    dtime = char(data.TIME(idx));
    if strcmp(date,'2016217')
        wbeg = data.WBEG(idx);
        wend = data.WEND(idx);
        spol = data.SPOL(idx);
        baz = data.BAZ(idx);
        sc_fast = data.FAST(idx);
        tlag = data.TLAG(idx);
        fileid = sprintf('%s_%s_%s',stat, date, dtime(1:4));
    
        trN = msac_read([path, fileid, '.BHN']);
        trE = msac_read([path, fileid, '.BHE']);
        samps = 1:length(trN.x1);
        time = trN.b + samps.*trN.delta;
        ind = find((time >=(wbeg-10)) & (time <=(wend+10)));
        fig = figure('Position',[500,500,1000,400]);
        [trF, trS] = msac_rotate(trN, trE, ifr_fast);
        % Assume a +ve dt* and attenuate fast trace
        trFA = msac_apply_tstar_operator(trF, fref, ifr_dtstar);
        plot(time(ind),trF.x1(ind), 'k-','LineWidth',2)
        hold on
        plot(time(ind), trS.x1(ind), 'r-','LineWidth',2)
        plot(time(ind), trFA.x1(ind), 'k--','LineWidth',2)
        %title(sprintf('%s = %4.2f',char(966), ifr_fast));
        xline(wbeg,'k-', 'LineWidth',2)
        xline(wend,'k-', 'LineWidth',2)
        title(['Station ',stat, '. Date: ',date, '-', dtime, '. \phi = ', ...
                num2str(ifr_fast), char(176), ', \Deltat* = ', num2str(ifr_dtstar), 's'],'FontSize',18)
        legend('"Fast"', '"Slow"', '"Fast" corrected','FontSize',18)
        xlabel('Time after event (s)','FontSize',18)
        xlim([wbeg-10, wend+10]);
        set(gca,'FontSize',22)
        exportgraphics(fig ,['/Users/ja17375/Projects/DeepMelt/Ethiopia/Figs/',stat, '_poster_2016217.eps'], ...
                                'BackgroundColor','none','Resolution',500)
    else
        disp('Skip')
    end
end