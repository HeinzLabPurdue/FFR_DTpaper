clear;
clc;

figSize_cm= [5 5 12.5 4.5]; % [Xcorner Ycorner Xwidth Ywidth]
saveFigs= 0;

figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

%% Init 
DirStruct.Codes= pwd;
DirStruct.Root= [fileparts(DirStruct.Codes) filesep];
DirStruct.Stimuli= [DirStruct.Root 'Stimuli' filesep];
DirStruct.Figures.eps= [DirStruct.Root 'Figures' filesep 'eps' filesep];
DirStruct.Figures.png= [DirStruct.Root 'Figures' filesep 'png' filesep];

allChinData = load([DirStruct.Root 'Data' filesep 'Output' filesep 'all_chins_data.mat']);
allChinData= allChinData.allChinData;

nh_inds= strcmp({allChinData.group}', 'NH');
hi_inds= strcmp({allChinData.group}', 'PTS') & ([allChinData.chinID]' ~=369);

nSegs= numel(allChinData(1).env_power_ffr);
nh_tfs= reshape([allChinData(nh_inds).tfs_power_ffr], nSegs, sum(nh_inds));
nh_env= reshape([allChinData(nh_inds).env_power_ffr], nSegs, sum(nh_inds));
nh_ffr_tfs2env_dB= db(nh_tfs./nh_env);

hi_tfs= reshape([allChinData(hi_inds).tfs_power_ffr], nSegs, sum(hi_inds));
hi_env= reshape([allChinData(hi_inds).env_power_ffr], nSegs, sum(hi_inds));
hi_ffr_tfs2env_dB= db(hi_tfs./hi_env);

nh_mean_t2e= nanmean(nh_ffr_tfs2env_dB, 2); %t2e = TFS-to-ENV ratio in dB
hi_mean_t2e= nanmean(hi_ffr_tfs2env_dB, 2);
nh_std_t2e= nanstd(nh_ffr_tfs2env_dB, [], 2);
hi_std_t2e= nanstd(hi_ffr_tfs2env_dB, [], 2);

hi_minus_nh_t2e= (hi_mean_t2e-nh_mean_t2e);
hi_minus_nh_t2e_norm= (hi_mean_t2e-nh_mean_t2e)./sqrt(nh_std_t2e.*hi_std_t2e);

[~, maxSegInd]= max(hi_minus_nh_t2e);
[~, minSegInd]= min(hi_minus_nh_t2e);

[~, maxSegInd_norm]= max(hi_minus_nh_t2e_norm);
[~, minSegInd_norm]= min(hi_minus_nh_t2e_norm);

fprintf('Segment index for max TFS2ENV (HI-NH)= %d, min TFS2ENV= %d\n', maxSegInd_norm, minSegInd_norm);

%% Read stimulus 
[sig, fs]= audioread([DirStruct.Stimuli 'FLN_Stim_S_P.wav']);
sig= sig/rms(sig);
tWindow= 64e-3;
tStart= 0;


%% Plot PSDs
plt.fontName= 'Arial';
plt.fontSize= 9;
plt.ax_plt.lw= 1;
plt.panel_fSize= 11;
plt.tick_len= [.03 .03];
plt.nw= 1.5;
plt.yRange= 65;
plt.nfft= 2^nextpow2(2+tWindow*fs);
plt.xtick_vals_freq= [100 1e3 5e3]/1e3;
plt.xtick_labs_freq= cellfun(@(x) num2str(x), num2cell(plt.xtick_vals_freq), 'UniformOutput', false);
plt.mrkSize= 6;
plt.lw= 1.5;

%% Left subplot
% Max Diff T2E (HI-NH)
% subplot(1,5,1:3);
ax(1)= subplot(1,4,1:2);
hold on;

% title('Max TFS2ENV ');
sig_max= sig(round((maxSegInd_norm-1)*tWindow*fs+1):round(maxSegInd_norm*tWindow*fs));
[~,~, lHan1]= helper.plot_dpss_psd(sig_max, fs, 'nw', plt.nw, 'nfft', plt.nfft, 'yrange', plt.yRange, 'xunit', 'khz');
set(gca, 'XTick', plt.xtick_vals_freq, 'XTickLabel', plt.xtick_labs_freq, 'LineWidth', plt.ax_plt.lw, 'FontName', plt.fontName, 'FontSize', plt.fontSize);
set(lHan1, 'linew', plt.lw, 'color', helper.get_color('g'));

box off;
% title('Max T2E (HI-NH)');

% Min Diff T2E (HI-NH)
% subplot(212);
sig_min= sig(round((minSegInd_norm-1)*tWindow*fs+1):round(minSegInd_norm*tWindow*fs));
[~,~, lHan2]= helper.plot_dpss_psd(sig_min, fs, 'nw', plt.nw, 'nfft', plt.nfft, 'yrange', plt.yRange, 'xunit', 'khz');
set(lHan2, 'linew', plt.lw, 'color', helper.get_color('prp'));
set(gca, 'XTick', plt.xtick_vals_freq, 'XTickLabel', plt.xtick_labs_freq, 'TickLength', plt.tick_len);
box off;

[lgHan, icons]= legend('High FFR T2E', 'Low  FFR T2E','box' ,'off', 'location', 'southwest');
icons(3).XData= mean(icons(3).XData) + [0.05 +.15];
icons(5).XData= mean(icons(5).XData) + [0.05 +.15];
lgHan.Position(1)= .035;
lgHan.Position(2)= .285;
set(lgHan, 'fontsize', 9);

xlim([50 5.1e3]/1e3)
txt(1)= text(60/1e3, -15, 'A','FontName',plt.fontName,'FontSize',plt.panel_fSize);
set(txt(1),'Units','norm')
set(txt(1),'pos',[0.05    .95         0])

xlab_han_A= xlabel('Frequency (kHz)', 'Units', 'normalized');
ylabel('PSD (dB/Hz)');

%%
plt.axis_lims= [-40 0];

% ax(2)= subplot(2,5, [4 5]);
ax(2)= subplot(1,4,3);
db_offset= -db(max([nh_env(:);nh_tfs(:);hi_env(:);hi_tfs(:)]));
NH_max_env= db(nh_env(maxSegInd,:))+db_offset;
NH_max_tfs= db(nh_tfs(maxSegInd,:))+db_offset;
NH_min_env= db(nh_env(minSegInd,:))+db_offset;
NH_min_tfs= db(nh_tfs(minSegInd,:))+db_offset;


hold on;
plot(NH_max_env, NH_max_tfs, 'o', 'Color', helper.get_color('g'), 'markersize', plt.mrkSize, 'LineWidth', plt.lw);
plot(NH_min_env, NH_min_tfs, 'd', 'Color', helper.get_color('prp'), 'markersize', plt.mrkSize, 'LineWidth', plt.lw);
txt(2)= text(.05,.95,'B. NH', 'color', helper.get_color('b'),'FontName',plt.fontName,'FontSize',plt.panel_fSize,'Units','norm');
plot(plt.axis_lims, plt.axis_lims, '--', 'LineWidth', plt.lw, 'Color', helper.get_color('gray'));
set(gca, 'YAxisLocation','left', 'TickLength', plt.tick_len, 'XAxisLocation', 'bottom', 'LineWidth', plt.ax_plt.lw, 'FontName', plt.fontName, 'FontSize', plt.fontSize);
ylabel('TFS_{FFR}^{ power} (dB)', 'FontAngle', 'italic', 'Units', 'normalized','FontName',plt.fontName);
xlab_han_B= xlabel('ENV_{FFR}^{ power} (dB)', 'FontAngle', 'italic', 'Units', 'normalized','FontName',plt.fontName);
xlab_han_B.Position(1)= 1.1;
xlab_han_B.Position(2)= xlab_han_A.Position(2);
%%
% ax(3)= subplot(2,5, [9 10]);
ax(3)= subplot(1,4,4);
HI_max_env= db(hi_env(maxSegInd,:))+db_offset;
HI_max_tfs= db(hi_tfs(maxSegInd,:))+db_offset;
HI_min_env= db(hi_env(minSegInd,:))+db_offset;
HI_min_tfs= db(hi_tfs(minSegInd,:))+db_offset;

hold on;
plot(HI_max_env, HI_max_tfs, 'o', 'Color', helper.get_color('g'), 'markersize', plt.mrkSize, 'LineWidth', plt.lw);
plot(HI_min_env, HI_min_tfs, 'd', 'Color', helper.get_color('prp'), 'markersize', plt.mrkSize, 'LineWidth', plt.lw);
txt(3)= text(.05,.95, 'C. HI', 'color', helper.get_color('r'),'FontName',plt.fontName,'FontSize',plt.panel_fSize,'Units','norm');

plot(plt.axis_lims, plt.axis_lims, '--', 'LineWidth', plt.lw, 'Color', helper.get_color('gray'));
set(gca, 'YAxisLocation','left', 'TickLength', plt.tick_len, 'XAxisLocation', 'bottom', 'YTickLabel', '', 'LineWidth', plt.ax_plt.lw, 'FontName', plt.fontName, 'FontSize', plt.fontSize);

ylabel('');

linkaxes(ax(2:3));
ylim(plt.axis_lims)
xlim(plt.axis_lims)

%% define common y-axis


%% save 
% set(gcf, 'units', 'inches', 'position', [21 1 16 4.5]);


%% Set axes placement/size
Xwidth1=.385;
Xwidth23=.18;
Xcorner=.0875;
Ycorner=.25;
Xshift=.04;
Ywidth=.7;
% Yshift=0.06;

% A
set(ax(1),'Position',[Xcorner Ycorner Xwidth1 Ywidth])
drawnow
% B
set(ax(2),'Position',[Xcorner+Xwidth1+3*Xshift Ycorner Xwidth23 Ywidth])
drawnow

set(ax(3),'Position',[Xcorner+Xwidth1+Xwidth23+4*Xshift Ycorner Xwidth23 Ywidth])
drawnow

figName= 'Fig11_PSD_extreme_T2E_segments';
% saveas(gcf, figName,'png');

if saveFigs
    saveas(gcf, [DirStruct.Figures.eps figName], 'epsc');
    saveas(gcf, [DirStruct.Figures.png figName], 'png');
end