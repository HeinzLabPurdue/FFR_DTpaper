clear;
clc;

%%
figSize_cm= [1 2 12.5 5.2]; % [Xcorner Ycorner Xwidth Ywidth]
figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm}; 
figure(1); 
clf
set(gcf,figure_prop_name,figure_prop_val);


%% Init 
DirStruct.Codes= pwd;
DirStruct.Root= [fileparts(DirStruct.Codes) filesep];
DirStruct.Stimuli= [DirStruct.Root 'Stimuli' filesep];
DirStruct.Figures.eps= [DirStruct.Root 'Figures' filesep 'eps' filesep];
DirStruct.Figures.png= [DirStruct.Root 'Figures' filesep 'png' filesep];

CodesDirs= {[DirStruct.Codes filesep 'chronux_2_11' filesep 'helper'], ...
    [DirStruct.Codes filesep 'chronux_2_11' filesep filesep 'continuous']};
addpath(CodesDirs{:});

%%
plotVar= 0;
freq_output_spread= 0;

saveFig= 0;

[stim, fsOld]= audioread([DirStruct.Stimuli 'FLN_Stim_S_P.wav']);
fsStim= 10e3;
stim= helper.gen_resample(stim, fsOld, fsStim);

tStim= (1:length(stim))/fsStim;

nh.Dir= [DirStruct.Root 'Data' filesep 'ArtifactRemovedFFR' filesep 'SP-2019_03_19-Q371_SFR_artifact_NH' filesep];
nh.allfiles= dir([nh.Dir '*.mat']);

nh.art_file= nh.allfiles(contains({nh.allfiles.name}, '_artifact')).name;
nh.reg_file= nh.allfiles(contains({nh.allfiles.name}, 'snr_120') & ~contains({nh.allfiles.name}, '_artifact')).name;

temp= load([nh.Dir nh.art_file]);
fs_old= temp.data.Stimuli.RPsamprate_Hz;
fs_data= 5e3;

temp.data.AD_Data.AD_Avg_PO_V{1}= helper.gen_resample(temp.data.AD_Data.AD_Avg_PO_V{1}, fs_old, fs_data);
temp.data.AD_Data.AD_Avg_NP_V{1}= helper.gen_resample(temp.data.AD_Data.AD_Avg_NP_V{1}, fs_old, fs_data);

fMax= fs_data/2;

clf;
nh.art_data.pos= helper.remove_artifact_ffr(temp.data.AD_Data.AD_Avg_PO_V{1}, fs_data, plotVar, fMax, freq_output_spread);
clf;
nh.art_data.neg= helper.remove_artifact_ffr(temp.data.AD_Data.AD_Avg_NP_V{1}, fs_data, plotVar, fMax, freq_output_spread);
clf;

temp= load([nh.Dir nh.reg_file]);
fs_old= temp.data.Stimuli.RPsamprate_Hz;
temp.data.AD_Data.AD_Avg_PO_V{1}= helper.gen_resample(temp.data.AD_Data.AD_Avg_PO_V{1}, fs_old, fs_data);
temp.data.AD_Data.AD_Avg_NP_V{1}= helper.gen_resample(temp.data.AD_Data.AD_Avg_NP_V{1}, fs_old, fs_data);


nh.reg_data.pos= helper.remove_artifact_ffr(temp.data.AD_Data.AD_Avg_PO_V{1}, fs_data, plotVar, fMax, freq_output_spread);
clf;
nh.reg_data.neg= helper.remove_artifact_ffr(temp.data.AD_Data.AD_Avg_NP_V{1}, fs_data, plotVar, fMax, freq_output_spread);
clf;

t_data= (1:length(nh.reg_data.pos))/fs_data;

HalfPowerFrequency1= 70;
HalfPowerFrequency2= 2e3;
N_bp_half= 4;

bpFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);

%
nh.reg_data.pos= filtfilt(bpFilt, nh.reg_data.pos);
nh.reg_data.neg= filtfilt(bpFilt, nh.reg_data.neg);
nh.reg_data.env= (nh.reg_data.pos+nh.reg_data.neg)/2;
nh.reg_data.tfs= (nh.reg_data.pos-nh.reg_data.neg)/2;

%
nh.art_data.pos= filtfilt(bpFilt, nh.art_data.pos);
nh.art_data.neg= filtfilt(bpFilt, nh.art_data.neg);
nh.art_data.env= (nh.art_data.pos+nh.art_data.neg)/2;
nh.art_data.tfs= (nh.art_data.pos-nh.art_data.neg)/2;




lw1= .5;
fSize= 9;
fName='Arial';
lg_fSize= 9;
panel_label_fSize= 11;
gain = 20e3;
tickVals= [0 .5 1 1.5];

stim= .5*stim/rms(stim)*rms(nh.reg_data.pos/gain*1e6);

Amax=2;
ax(1)=subplot(121);
helper.set_colblind_order();


hold on;

tHan(1)= plot(t_data, Amax+nh.reg_data.pos/gain*1e6, 'LineWidth', lw1);
tHan(2)= plot(t_data, nh.art_data.pos/gain*1e6, 'LineWidth', lw1);
plot(tStim, -Amax+stim, 'k', 'LineWidth', lw1);
ylabel('Amplitude (\muV)');

tick_len= [.025 .025];
set(gca, 'FontSize', fSize,'fontname', fName, 'LineWidth', 1, 'TickLength', tick_len, 'XTick', tickVals)
xlabel('Time (sec)');

% l0= plot(nan, nan, 'k', 'linew', 1);

yRange= 42;
nw=5;
nfft= round(fs_data);
xtick_vals= [20 100 500 2e3];
xtick_labs= cellfun(@(x) num2str(x), num2cell(xtick_vals), 'UniformOutput', false);

stim= stim/400;

ax(2)=subplot(122);
helper.set_colblind_order();

hold on;
[~,~,l3]= helper.plot_dpss_psd(stim*25, fsStim, 'nfft', nfft, 'nw', nw, 'plotconf', false);
set(gca, 'ColorOrderIndex', 1);
[~,~,l1]= helper.plot_dpss_psd(nh.reg_data.pos, fs_data, 'nfft', nfft, 'nw', nw, 'plotconf', true);
[~,~,l2]= helper.plot_dpss_psd(nh.art_data.pos, fs_data, 'nfft', nfft, 'nw', nw, 'plotconf', true);
set(l3, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1);
set(l1, 'LineWidth', 1);
set(l2, 'LineWidth', 1);
axis tight;
yl= ylim;
ylimHard= [max(yl)-yRange+3 max(yl)+3];
ylim(ylimHard);
xlim([20 2e3]);

set(gca, 'FontSize', fSize, 'fontname', fName, 'XTick', xtick_vals, 'XTickLabel', xtick_labs, 'LineWidth', 1, 'TickLength', tick_len)
ylabel('PSD (dB/Hz)');
title('')
xlabel('')
grid off;

xVert= [60, 500];

matGreen= helper.get_color('g');
plot(xVert, [min(ylimHard) min(ylimHard)]+2, 'Color', matGreen, 'LineWidth', 3, 'LineStyle', '-');

[lg, icons]=legend([l1(1) l2(1) l3], '{\bf+}eartip', '{\bf-}eartip', 'stim', 'Location', 'northeast');


lg.FontSize= lg_fSize;
lg.Box= 'off';
lg.Position(1)= .77;
lg.Position(2)= .7;

%% legend
icons(4).XData= [0.25 .45];
icons(6).XData= [0.25 .45];
icons(8).XData= [0.25 .45];

xlabel('Frequency (Hz)')

txtHan= helper.add_subplot_letter(1, 2, panel_label_fSize);


%% Set axes placement/size
Xwidth=.385;
Xcorner=.075;
Xshift=.12;
Ywidth=.73;
Ycorner=.18;
% Yshift=0.06;

% A
set(ax(1),'Position',[Xcorner Ycorner Xwidth Ywidth])
set(txtHan(1),'pos',get(txtHan(1),'pos')+[0 0.02 0],'FontName',fName)
drawnow
% B
set(ax(2),'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
set(txtHan(2),'pos',get(txtHan(2),'pos')+[0 0.02 0],'FontName',fName)
drawnow

if saveFig
%     saveas(gcf, [DirStruct.Figures.eps 'Fig2_compare_no_eartip'], 'epsc');
%     saveas(gcf, [DirStruct.Figures.png 'Fig2_compare_no_eartip'], 'png');
    print([DirStruct.Figures.eps 'Fig2_compare_no_eartip.tiff'], '-dtiff' ,'-r600');
    print([DirStruct.Figures.eps 'Fig2_compare_no_eartip.png'], '-dpng' ,'-r600');
    print([DirStruct.Figures.png 'Fig2_compare_no_eartip.png'], '-dpng' ,'-r600');
end

rmpath(CodesDirs{:});