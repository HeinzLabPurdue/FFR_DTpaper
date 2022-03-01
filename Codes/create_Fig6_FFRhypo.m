clear;
clc;

figSize_cm= [5 5 8.3 7]; % [Xcorner Ycorner Xwidth Ywidth]
saveFig= 1;
%% Init
DirStruct.Codes= pwd;
DirStruct.Root= [fileparts(DirStruct.Codes) filesep];
DirStruct.Stimuli= [DirStruct.Root 'Stimuli' filesep];
DirStruct.Figures.eps= [DirStruct.Root 'Figures' filesep 'eps' filesep];
DirStruct.Figures.png= [DirStruct.Root 'Figures' filesep 'png' filesep];

figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

co= get(gca, 'ColorOrder');
co= co([1 2 3 3 5 4 4], :);
plt.col_gray= .4*ones(1,3);
plt.fontName='Arial';
plt.tick_len= [.025 .025];
plt.ytick_val= -80:40:0;
plt.lw= 1;

nh.tc= load([DirStruct.Root 'Data' filesep 'ANtc' filesep 'Q321_spikestimulusData.mat']);
nh.tc= nh.tc.spike_data;
nh.ind= [nh.tc.track]==2 & [nh.tc.unit]==7;
nh.tc= nh.tc(nh.ind).TC;

hi.tc= load([DirStruct.Root 'Data' filesep 'ANtc' filesep 'Q362_spikestimulusData.mat']);
hi.tc= hi.tc.spike_data;
hi.ind= [hi.tc.track]==3 & [hi.tc.unit]==3;
hi.tc= hi.tc(hi.ind).TC;

rng(0)
fs= 20e3;
dur= 5;
freqs= logspace(1, 4, 100);

amp= 10.^(-6*log2(freqs./freqs(1))/20);

[s1, t]= helper.create_sinusoid(freqs, fs, dur, amp);

plt.xticks_vals= [100 500 1e3 5e3 10e3]/1e3;
plt.xticks_labs= cellfun(@(x) num2str(x), num2cell(plt.xticks_vals), 'UniformOutput', false);
plt.fontSize= 9;
plt.nfft= 2^16;
plt.ax_lw= 1;

sp_ax(1)= subplot(211);
hold on
set(gca, 'ColorOrder', co, 'TickLength', plt.tick_len, 'LineWidth', plt.ax_lw, 'FontName', plt.fontName);
% fill([.5 11 11 .5 .5], [-150 -150 0 0 -150], [.8 .8 .9], 'FaceAlpha', .25);
[Pxx_dB, freq_psd, lHan]= helper.plot_dpss_psd(s1, fs, 'xunit', 'khz', 'nfft', plt.nfft);
set(lHan, 'linew', plt.lw);
xlim([95 1.05e4]/1e3);
% title('NH');
set(gca, 'XTick', plt.xticks_vals, 'XTickLabel', plt.xticks_labs, 'FontSize', plt.fontSize, 'YTick', plt.ytick_val, 'XTickLabel', '');
ylim([-120 0]);
grid off;
xlabel('');

sp_ax(2)= subplot(212);
set(gca, 'ColorOrder', co, 'TickLength', plt.tick_len, 'LineWidth', plt.ax_lw, 'FontName', plt.fontName);
hold on
% fill([.5 11 11 .5 .5], [-150 -150 0 0 -150], [.8 .8 .9], 'FaceAlpha', .25);
[~,~,lHan]= helper.plot_dpss_psd(s1, fs, 'xunit', 'khz', 'nfft', plt.nfft);
set(lHan, 'linew', plt.lw);
xlim([95 1.05e4]/1e3);
% title('HI');
set(gca, 'XTick', plt.xticks_vals, 'XTickLabel', plt.xticks_labs, 'FontSize', plt.fontSize, 'YTick', plt.ytick_val);
psd_freq_inds_below1k= freq_psd<0.5;
plot(freq_psd(psd_freq_inds_below1k), Pxx_dB(psd_freq_inds_below1k), 'color', plt.col_gray, 'linew', plt.lw);
set(gca, 'ColorOrder', co);
ylim([-120 0]);

nhFreqs= [.33 1 3];
nh_thresh= -75;
plt.lw_tc= 1.75;
freqs_kHz= freqs/1e3;
grid off;

freqSpread_around_cf= 1.1;
for freqVar= 1:length(nhFreqs)
    ind_cur_freq= dsearchn(freqs_kHz', nhFreqs(freqVar));
    tc_new= helper.get_shifted_tc(nh.tc, freqs_kHz(ind_cur_freq), db(amp(ind_cur_freq))-45, 'nh');
    
    psd_freq_ind1= dsearchn(freq_psd, freqs_kHz(ind_cur_freq)/freqSpread_around_cf);
    psd_freq_ind2= dsearchn(freq_psd, freqs_kHz(ind_cur_freq)*freqSpread_around_cf);
    
    subplot(211);
    ll= plot(tc_new.freqkHz, tc_new.TCfit, 'k-', 'linew', plt.lw_tc);
    plot(freq_psd(psd_freq_ind1:psd_freq_ind2), Pxx_dB(psd_freq_ind1:psd_freq_ind2), '-', 'color', plt.col_gray, 'linew', plt.lw);
    
    subplot(212);
    tc_new= helper.get_shifted_tc(hi.tc, freqs_kHz(ind_cur_freq), -60, 'hi');
    plot(tc_new.freqkHz, tc_new.TCfit, 'color', get(ll, 'color'), 'linew', plt.lw_tc);
    plot(tc_new.freqkHz, tc_new.TCfit, 'k-',  'linew', 2);
    
end

% set(gcf, 'Units', 'inches', 'Position', [1 1 8.5 7]);

%% work on labels


fig_SP_pan_A= subplot(211);
xlabel('');
ylabel('');
text(0, 1.1, '\bfA. NH\rm', 'FontSize', 11, 'Units', 'normalized', 'FontName', plt.fontName);
fig_SP_pan_B= subplot(212);
xlabel('Frequency (kHz)');
ylab_han= ylabel('PSD (dB/Hz)', 'FontSize', plt.fontSize);
text(0, 1.1, '\bfB. HI\rm', 'FontSize', 11, 'Units', 'normalized', 'FontName', plt.fontName);
ylab_han.Position(1)= .065;
ylab_han.Position(2)= 10;

%% Set axes placement/size
Xwidth=.85;
Ywidth=.35;
Xcorner=.125;
Yshift=.09;
Ycorner=.14;
% Yshift=0.06;

% A
set(sp_ax(2),'Position',[Xcorner Ycorner Xwidth Ywidth])
drawnow
% B
set(sp_ax(1),'Position',[Xcorner Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

%%

if saveFig
    saveas(gcf, [DirStruct.Figures.eps 'Fig6_ffr_hypo_nh_hi'], 'epsc');
    saveas(gcf, [DirStruct.Figures.png 'Fig6_ffr_hypo_nh_hi'], 'png');
end