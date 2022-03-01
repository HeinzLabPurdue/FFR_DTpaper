clear;
clc;

figSize_cm= [5 5 12.5 6.2]; % [Xcorner Ycorner Xwidth Ywidth]
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

chinID= 371;
saveFig= 0;

%% Init
DirStruct.Codes= pwd;
DirStruct.Root= [fileparts(DirStruct.Codes) filesep];
DirStruct.Stimuli= [DirStruct.Root 'Stimuli' filesep];
DirStruct.Figures.eps= [DirStruct.Root 'Figures' filesep 'eps' filesep];
DirStruct.Figures.png= [DirStruct.Root 'Figures' filesep 'png' filesep];


data_dir= [DirStruct.Root 'Data' filesep 'ArtifactRemovedFFR' filesep 'SP-2019_05_08-Q371_SFRpink500Hz_NH' filesep];
remove_artifact_here= 0; % because using cleaned data => already artifact removed
if ~remove_artifact_here
    fprintf('Assuming using artifact removed Data\n');
end
tStart= 0; tEnd= 1.3;

allfiles= dir([data_dir 'a*SFR*.mat']);
allfiles= allfiles(~(contains({allfiles.name}, 'latency') | contains({allfiles.name}, 'artifact')));

all_snrs= cell2mat(cellfun(@(x) str2double(strrep(x(regexp(x, 'snr_')+4 : regexp(x, '_atn')-1), 'm', '-')), {allfiles.name}, 'uniformoutput', false));
all_snrs(isnan(all_snrs))= [];
all_snrs= fliplr(unique(all_snrs));

[sig, fs_sig]= audioread([DirStruct.Stimuli 'FLN_Stim_S_P.wav']);
snrVar= 1;
[s_data_pos_filt, s_data_neg_filt, s_nf_filt, fs_data]= helper.get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs, remove_artifact_here);

snrVar= find(all_snrs==0);
[sn_data_pos_filt_p10, sn_data_neg_filt_p10, sn_nf_pos_filt_p10, fs_data_p10]= helper.get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs, remove_artifact_here);

snrVar= find(all_snrs==-20);
[sn_data_pos_filt_m20, sn_data_neg_filt_m20, sn_nf_pos_filt_m20, fs_data_m20]= helper.get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs, remove_artifact_here);

sig= helper.gen_resample(sig, fs_sig, fs_data_p10);
fs_sig= fs_data_p10;

plot_S_sig= 1;
legHan1= helper.create_panel_plot_s_vs_sn_tdt_at2snrs(1, fs_sig, sig, fs_data_p10, sn_data_pos_filt_p10, sn_data_neg_filt_p10, s_data_pos_filt, s_data_neg_filt, ...
    s_nf_filt, tStart, tEnd, plot_S_sig);

plot_S_sig= 0;
[legHan2, SPhandles]= helper.create_panel_plot_s_vs_sn_tdt_at2snrs(1, fs_sig, sig, fs_data_m20, sn_data_pos_filt_m20, sn_data_neg_filt_m20, s_data_pos_filt, s_data_neg_filt, ...
    s_nf_filt, tStart, tEnd, plot_S_sig);

sp_ax(1)= subplot(121); yyaxis left; leftHan(1)= gca;
xlab_han= xlabel('Frequency (Hz)', 'FontSize', 9);
sp_ax(2)= subplot(122); yyaxis left; leftHan(2)= gca;
linkaxes(leftHan, 'y');
ylim([-75 -35])

legHans= [legHan1([1 2 4]); legHan2; legHan1(3)];
[lg, icons]= legend(legHans, 'Stim', 'Quiet', '0 dB', '-20 dB', 'NF', 'box', 'off', 'location', 'northeast');
icons(6).XData= mean(icons(6).XData) + [0 +.2];
icons(8).XData= mean(icons(8).XData) + [0 +.2];
icons(10).XData= mean(icons(10).XData) + [0 +.2];
icons(12).XData= mean(icons(12).XData) + [0 +.2];
icons(14).XData= mean(icons(14).XData) + [0 +.2];

lg.FontSize= 9;
lg.Position(1)= .28;
lg.Position(2)= .6;
% set(gcf, 'units', 'centimeters', 'position', [50 3 30 16], 'Renderer','painters');

%% Set axes placement/size
Xwidth=.37;
Xcorner=.095;
Ycorner=.15;
Xshift=.07;
Ywidth=.77;
% Yshift=0.06;

% A
set(sp_ax(1),'Position',[Xcorner Ycorner Xwidth Ywidth])
drawnow
% B
set(sp_ax(2),'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
drawnow

% Xlabel
set(xlab_han, 'pos', [1e3 -79 -1]);
drawnow

fName_summary= sprintf('Fig4_nh_tdt_pink_example');
if saveFig
    print([DirStruct.Figures.eps fName_summary '.tiff'], '-dtiff' ,'-r600');
    print([DirStruct.Figures.eps fName_summary '.png'], '-dpng' ,'-r600');
%     saveas(gcf, [DirStruct.Figures.eps fName_summary], 'tiff');
%     saveas(gcf, [DirStruct.Figures.eps fName_summary], 'epsc');
%     saveas(gcf, [DirStruct.Figures.eps fName_summary], 'png');
    saveas(gcf, [DirStruct.Figures.png fName_summary], 'png');
end