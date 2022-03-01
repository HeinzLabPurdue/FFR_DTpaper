clear;
clc;

saveFigure= 0;
figSize_cm= [1 1.5 8.3 8]; % [Xcorner Ycorner Xwidth Ywidth]
%% Init 
DirStruct.Codes= pwd;
DirStruct.Root= [fileparts(DirStruct.Codes) filesep];
DirStruct.Stimuli= [DirStruct.Root 'Stimuli' filesep];
DirStruct.Figures.eps= [DirStruct.Root 'Figures' filesep 'eps' filesep];
DirStruct.Figures.png= [DirStruct.Root 'Figures' filesep 'png' filesep];

anl.tMax= 1.3;
anl.fs= 20e3;
anl.winRes= 40e-3;
anl.nfft= 2^16;

[sig, anl.fsOrg]= audioread([DirStruct.Stimuli 'FLN_Stim_S_P.wav']);
sig= sig(1:ceil(anl.fsOrg*anl.tMax));

sig = helper.gen_resample(sig, anl.fsOrg, anl.fs);
t= (1:length(sig))/anl.fs;

voiced_times= helper.find_voicing_boundaries(sig, anl.fs);
voice_mask = sum(repmat(t, size(voiced_times,1), 1) > repmat(voiced_times(:,1), 1, numel(t)) &  ...
    repmat(t, size(voiced_times,1), 1) < repmat(voiced_times(:,2), 1, numel(t)), 1);

%% Plot parameters
plt.fontName='Arial';
plt.fontsize= 9;
plt.txSize= 11;
plt.ax_lw= 1;
plt.lw1= .75;
plt.lw2= 1.5;
plt.tick_freq= 0:4;
plt.tick_amp= [-.3 0 .3];
plt.tick_time= [0 .5 1 1.3];
plt.tick_mask= [0 1];
plt.tick_len= [.02 .02];

%% Plot 
figure(3);
clf;
figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1); clf
set(gcf,figure_prop_name,figure_prop_val);

ax(1)= subplot(5,1,1:4);
spectrogram(sig, blackman(round(anl.winRes*anl.fs)), round(anl.winRes*anl.fs*.9), anl.nfft, anl.fs, 'MinThreshold',-80, 'yaxis');
ylim([0 4e3]/1e3);
set(gca, 'fontsize', plt.fontsize, 'fontname', plt.fontName, 'xtick', []);

colorbar off
xlabel('');
set(gca, 'yTick', plt.tick_freq, 'Box', 'off', 'LineWidth', plt.ax_lw, 'TickLength', plt.tick_len);
text(0, 1.05, '\bfA\rm', 'FontSize', plt.txSize, 'fontname', plt.fontName, 'Units', 'normalized');

ax(2)= subplot(5,1,5);
plot(t, sig, 'linew', plt.lw1);
set(gca, 'XTick', plt.tick_time, 'YTick', plt.tick_amp)

yyaxis right;
plot(t, voice_mask, 'linew', plt.lw2);
ylim([-.5 1.1]);
hy_Voiced=ylabel('Voiced');
set(gca, 'YTick', plt.tick_mask);

yyaxis left;
xlabel('Time (sec)');
ylabel('Amp (Pa)');
set(gca, 'fontsize', plt.fontsize, 'fontname', plt.fontName, 'LineWidth', plt.ax_lw, 'TickLength', plt.tick_len);
linkaxes(ax, 'x');
box off;

axis tight;
text(0, 1.13, '\bfB\rm', 'FontSize', plt.txSize, 'fontname', plt.fontName, 'Units', 'normalized');

set(hy_Voiced,'units','norm')
get(hy_Voiced,'pos')
set(hy_Voiced,'pos',[1.05    0.5000  0])


%% Set axes placement/size
Xwidth=.82;
Xcorner=.092;
Yshift=.07;
Ycorner=.12;
Ywidth1=.56;
Ywidth2=.2;

% B
set(ax(2),'Position',[Xcorner Ycorner Xwidth Ywidth2])
drawnow
% A
set(ax(1),'Position',[Xcorner Ycorner+Ywidth2+Yshift Xwidth Ywidth1])
drawnow

%%


%% Save
fName_spectrogram= 'Fig1_stim_spectrogram';
if saveFigure
    saveas(gcf, [DirStruct.Figures.eps fName_spectrogram], 'epsc');
    saveas(gcf, [DirStruct.Figures.png fName_spectrogram], 'png');
end