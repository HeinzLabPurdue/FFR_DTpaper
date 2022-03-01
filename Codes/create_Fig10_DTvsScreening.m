clear;
clc;

figSize_cm= [5 5 8.3 4.1]; % [Xcorner Ycorner Xwidth Ywidth]
%% Init 
DirStruct.Codes= pwd;
DirStruct.Root= [fileparts(DirStruct.Codes) filesep];
DirStruct.Figures.eps= [DirStruct.Root 'Figures' filesep 'eps' filesep];
DirStruct.Figures.png= [DirStruct.Root 'Figures' filesep 'png' filesep];

addpath(DirStruct.Codes);

ABR_DataDir= [DirStruct.Root 'Data' filesep 'ABRfiles' filesep];
% DPoae_rootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/Baselines/'; % this should point to the ROOT folder that has all baseline EXPdata directories.
DPoae_rootDataDir= [fileparts(pwd) filesep 'BaselineFiles' filesep];

allChinData = load([DirStruct.Root 'Data' filesep 'Output' filesep 'all_chins_data.mat']);
allChinData = allChinData.allChinData;

freqs2use= [.5 1 2 4 8 0]*1e3;
midFreqs= [1 2 4 8]*1e3;
plt.fSize= 9;
saveFig= 0;

for iterVar= 1:length(allChinData)
    chinID= allChinData(iterVar).chinID;
    chinType= allChinData(iterVar).group;
    
    lf_power_audio= helper.dbspl(sqrt(allChinData(iterVar).lf_power_audio));
    hf_power_audio= helper.dbspl(sqrt(allChinData(iterVar).hf_power_audio));
    
    tfs_power_ffr= db(sqrt(allChinData(iterVar).tfs_power_ffr));
    env_power_ffr= db(sqrt(allChinData(iterVar).env_power_ffr));
    
    notNANinds= ~(isnan(lf_power_audio) | isnan(hf_power_audio) | isnan(tfs_power_ffr) | isnan(env_power_ffr));
    
    if any(notNANinds)
        
        mdl_raw= fitlm(lf_power_audio(notNANinds), tfs_power_ffr(notNANinds));
        allChinData(iterVar).DTslope_raw = mdl_raw.Coefficients.Estimate(2);
        
        mdl_norm= fitlm(hf_power_audio(notNANinds)- lf_power_audio(notNANinds), env_power_ffr(notNANinds) - tfs_power_ffr(notNANinds));
        allChinData(iterVar).DTslope_norm = mdl_norm.Coefficients.Estimate(2);
        
        if any(strcmpi(chinType, {'PTS', 'HI'}))
            chinType= 'HI';
        end
        
        %% Do ABR analyis
        curDir= dir(sprintf('%sQ%d_%s*', ABR_DataDir, chinID, chinType));
        if ~isempty(curDir)
            
            if numel(curDir)>1
                warning('Check that the right directory is selected');
                curDir= curDir(1);
            end
            
            fprintf('Using %s\n', curDir.name);
            
            curData_fName= dir(sprintf('%s%s%s*.mat', ABR_DataDir, curDir.name, filesep));
            curData= load(sprintf('%s%s%s%s', ABR_DataDir, curDir.name, filesep, curData_fName.name));
            curData= curData.abrs;
            [~, inds]= ismember(freqs2use, curData.thresholds(:,1));
            
            allChinData(iterVar).abr_thresh= nan(size(freqs2use));
            allChinData(iterVar).abr_thresh= curData.thresholds(find(inds),2); %#ok<FNDSB>
            
            allChinData(iterVar).midFreq_thresh= nanmean(allChinData(iterVar).abr_thresh(ismember(freqs2use, midFreqs)));
            
        else
            fprintf('%s \n ', '---');
            allChinData(iterVar).midFreq_thresh = nan;
        end
        
        %% Do DPOAE analysis

        allChinDirs= dir([DPoae_rootDataDir '*' num2str(chinID) '*']);
        if~isempty(allChinDirs)
            if strcmp(chinType, 'NH')
                dirNum= find(contains({allChinDirs.name}, {'nh', 'pre'},'IgnoreCase',true));
            elseif strcmp(chinType, 'PTS') | strcmp(chinType, 'HI')
                dirNum= find(contains({allChinDirs.name}, {'pts', 'post', 'follow', 'hi'}));
            else
                dirNum= nan;
            end
            if ~isempty(dirNum) && any(~isnan(dirNum))
                if numel(dirNum)>1
                    warning('Check that the right directory is selected');
                    dirNum= dirNum(1);
                end
                DataDir= allChinDirs(dirNum).name;
                dpFile= dir([DPoae_rootDataDir DataDir filesep '*dpoae*']);
                dpFile= [DPoae_rootDataDir DataDir filesep dpFile(1).name];
                calibFile= helper.get_lower_calibFile(dpFile);
                
                run(calibFile);
                calibData=ans;
                calibData=calibData.CalibData;
                
                out_DPOAE_data= helper.my_dpoae_analysis(dpFile);
                dpData=[out_DPOAE_data.dp_amp];
                dp_freqs= [out_DPOAE_data.freq2];
                
                freq_range2consider= [min(midFreqs) max(midFreqs)];
                inds2consider= dp_freqs>freq_range2consider(1) & dp_freqs<freq_range2consider(2);
                allChinData(iterVar).midFreq_DPamp= nanmean(dpData(inds2consider));
            else
                allChinData(iterVar).midFreq_DPamp= nan;
            end
        else
            allChinData(iterVar).midFreq_DPamp= nan;
        end
    else
        allChinData(iterVar).midFreq_thresh= nan;
        allChinData(iterVar).midFreq_DPamp= nan;
        allChinData(iterVar).DTslope_raw= nan;
        allChinData(iterVar).DTslope_norm= nan;
    end
    
end

abr_thresh_db= [allChinData.midFreq_thresh];
dp_amp_db= [allChinData.midFreq_DPamp];
dt_slope_raw= [allChinData.DTslope_raw];
dt_slope_norm= [allChinData.DTslope_norm];

%%
plt.fontName='Arial';
plt.plotYes= 1;
plt.verbose= 1;
plt.lw=1;
plt.lw2= 2;
plt.ax_lw= 1;
plt.mrkSize= 6;

figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

special_cond= [allChinData.chinID] == 369 & strcmp({allChinData.group}, 'PTS');

plt.tick_len= [.025 .025];
sp_ax(1)= subplot(121);
hold on;
plot(abr_thresh_db, dt_slope_norm, 'kd', 'linew', plt.lw, 'markersize', plt.mrkSize)
plot(abr_thresh_db(special_cond), dt_slope_norm(special_cond), 'd', 'color', helper.get_color('lr'), 'linew', plt.lw2, 'markersize', plt.mrkSize)
lHan= helper.plot_fitlm(abr_thresh_db, dt_slope_norm, plt.plotYes, plt.verbose);
set(lHan, 'linew', plt.lw2);

xlabel('ABR Thr. (dB SPL)');
ylabel('DT_{slope}', 'FontAngle', 'italic');
set(gca, 'TickLength', plt.tick_len, 'linew', plt.ax_lw, 'FontName', plt.fontName, 'FontSize', plt.fSize);
grid off;
xlim([12 40]);

sp_ax(2)= subplot(122);
hold on;
plot(dp_amp_db, dt_slope_norm, 'kd', 'linew', plt.lw, 'markersize', plt.mrkSize)
plot(dp_amp_db(special_cond), dt_slope_norm(special_cond), 'd', 'color', helper.get_color('lr'), 'linew', plt.lw2, 'markersize', plt.mrkSize)
lHan= helper.plot_fitlm(dp_amp_db, dt_slope_norm, plt.plotYes, plt.verbose);
set(lHan, 'linew', plt.lw2);

xlabel('DP Amp. (dB SPL) ');
grid off;
set(gca, 'TickLength', plt.tick_len, 'linew', plt.ax_lw, 'FontName', plt.fontName);

set(findall(gcf,'-property','FontSize'),'FontSize',plt.fSize)

% set(gcf, 'Units', 'inches', 'Position', [1 1 11 5]);
txtHan= helper.add_subplot_letter(1, 2, 11);
set(txtHan(1),'pos',get(txtHan(1),'pos')+[0 0.02 0],'FontName',plt.fontName)
set(txtHan(2),'pos',get(txtHan(2),'pos')+[0 0.02 0],'FontName',plt.fontName)

rmpath(DirStruct.Codes);

%% Set axes placement/size
Xwidth=.32;
Ywidth=.69;
Xcorner=.16;
Xshift=.14;
Ycorner=.22;
% Yshift=0.06;

% A
set(sp_ax(1),'Position',[Xcorner Ycorner Xwidth Ywidth])
drawnow
% B
set(sp_ax(2),'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
drawnow


if saveFig
    saveas(gcf, [DirStruct.Figures.png 'Fig10_baseline_vs_slope'], 'png');
    saveas(gcf, [DirStruct.Figures.eps 'Fig10_baseline_vs_slope'], 'epsc');
end
