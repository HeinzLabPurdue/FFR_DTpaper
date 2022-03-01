clear;
clc;
% close all;


figSize_cm= [5 5 12.5 5.2]; % [Xcorner Ycorner Xwidth Ywidth]
saveFigs= 0;
allChins= [358 360 366 367 370 371 373 374 379 369];

%% Init 
DirStruct.Codes= pwd;
DirStruct.Root= [fileparts(DirStruct.Codes) filesep];
DirStruct.Figures.eps= [DirStruct.Root 'Figures' filesep 'eps' filesep];
DirStruct.Figures.png= [DirStruct.Root 'Figures' filesep 'png' filesep];

addpath(DirStruct.Codes);

DataDir= [DirStruct.Root 'Data' filesep 'ABRfiles' filesep];

types= {'NH', 'HI'};
freqs2use_Hz= [.5 1 2 4 8]*1e3;
freqs2use_kHz= freqs2use_Hz/1e3;

thresh_data.nh.z= {};
thresh_data.hi.z= {};
thresh_data.nh.amp= {};
thresh_data.hi.amp= {};
nhChins= [];
hiChins= [];

for chinVar= 1:length(allChins)
    cur_ChinID= allChins(chinVar);
    for typeVar= 1:length(types)
        cur_type= types{typeVar};
        curDir= dir(sprintf('%sQ%d_%s*', DataDir, cur_ChinID, cur_type));
        if numel(curDir)==2 && cur_ChinID==369
            % Q369 has 2 post-exposure baselines. 
            curDir= curDir(strcmp({curDir.name}', 'Q369_HI_2019_01_16'));
        end
        if ~isempty(curDir)
            fprintf('Using %s\n', curDir.name);
            curData_fName= dir(sprintf('%s%s%s*.mat', DataDir, curDir.name, filesep));
            curData= load(sprintf('%s%s%s%s', DataDir, curDir.name, filesep, curData_fName.name));
            curData= curData.abrs;
            [~, inds]= ismember(freqs2use_Hz, curData.thresholds(:,1));
            if strcmp(cur_type, 'NH')
                thresh_data.nh.z= [thresh_data.nh.z, curData.thresholds(inds,2)];
                nhChins= [nhChins; cur_ChinID]; %#ok<*AGROW>
                %                 thresh_data.nh.amp= [thresh_data.nh.amp, curData.thresholds(:,3)];
            elseif strcmp(cur_type, 'HI')
                thresh_data.hi.z= [thresh_data.hi.z, curData.thresholds(inds,2)];
                hiChins= [hiChins; cur_ChinID];
                %                 thresh_data.hi.amp= [thresh_data.hi.amp, curData.thresholds(:,3)];
            end
            
        else
            warning('No %s data for Q%d\n', cur_type, cur_ChinID);
        end
        
    end
end


%% plot
plt.fontName='Arial';
plt.mrkSize= 6;
plt.lw1= 1;
plt.lw3= 2.5;
plt.fontSize= 9;
plt.tick_len= [.02 .02];
plt.ax_lw= 1;

nh_data= cell2mat(thresh_data.nh.z)';
hi_data= cell2mat(thresh_data.hi.z)';
outlier_chin= 369;

reg_ind_nh= find(ismember(nhChins, setxor(nhChins, outlier_chin)));
reg_ind_hi= find(ismember(hiChins, setxor(nhChins, outlier_chin)));
outlier_ind_nh= find(ismember(nhChins, outlier_chin));
outlier_ind_hi= find(ismember(hiChins, outlier_chin));


figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

sp_ax(1)= subplot(121);
[~, co_struct]= helper.set_colblind_order();
hold on;
plot(freqs2use_kHz, nh_data(reg_ind_nh, :)', '-o', 'color', co_struct.b, 'markersize', plt.mrkSize, 'linew', plt.lw1);
plot(freqs2use_kHz, nh_data(outlier_ind_nh, :)', '-d', 'color', co_struct.lb, 'markersize', plt.mrkSize, 'linew', plt.lw1);
plot(freqs2use_kHz, nanmean(nh_data(reg_ind_nh, :), 1), '-', 'color', 'b', 'markersize', plt.mrkSize, 'linew', plt.lw3);

plot(freqs2use_kHz, hi_data(reg_ind_hi, :)', '-o', 'color', co_struct.r,'markersize', plt.mrkSize, 'linew', plt.lw1);
plot(freqs2use_kHz, hi_data(outlier_ind_hi, :)', '-d', 'color', co_struct.lr,'markersize', plt.mrkSize, 'linew', plt.lw1);
plot(freqs2use_kHz, nanmean(hi_data(reg_ind_hi, :), 1), '-', 'color', 'r','markersize', plt.mrkSize, 'linew', plt.lw3);

set(gca, 'xscale', 'log', 'xtick', freqs2use_kHz, 'TickLength', plt.tick_len);
xlim([.4 10]);
xlab_han= xlabel('Frequency (kHz)');
ylabel('ABR Threshold (dB SPL)');
set(gca, 'fontsize', plt.fontSize, 'box', 'off', 'linew', plt.ax_lw, 'FontName', plt.fontName);


%%
sp_ax(2)= subplot(122);
[~, co_struct]= helper.set_colblind_order();
% warning('This should be the path to raw baseline/followup directories');
% DPoae_rootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/Baselines/'; 
DPoae_rootDataDir= [fileparts(pwd) filesep 'BaselineFiles' filesep];

xTicks= freqs2use_Hz;

% each column is for one animal
freq27= load([DirStruct.Root 'Data' filesep 'default_27freq_DPOAE.mat']);
freq27= freq27.default_27freq_DPOAE;
dp_data_nh= nan(27, length(allChins));
dp_data_hi= nan(27, length(allChins));

for chinVar=1:length(allChins)
    chinID=allChins(chinVar);
    allChinDirs= dir([DPoae_rootDataDir '*' num2str(chinID) '*']);
    
    NH.dirNum= find(contains(lower({allChinDirs.name}'), {'pre', 'nh'}));
    NH.DataDir= allChinDirs(NH.dirNum).name;
    NH.dpFile= dir([DPoae_rootDataDir NH.DataDir filesep '*dpoae*']);
    NH.dpFile= [DPoae_rootDataDir NH.DataDir filesep NH.dpFile(1).name];
    NH.calibFile= helper.get_lower_calibFile(NH.dpFile);
    
    HI.dirNum= find(contains(lower({allChinDirs.name}'), {'post', 'hi', 'pts', 'tts', 'follow'}));
    existHI= ~isempty(find(contains(lower({allChinDirs.name}'), {'post', 'hi', 'pts', 'tts', 'follow'}), 1));
    
    if existHI
        HI.DataDir= allChinDirs(HI.dirNum).name;
        HI.dpFile= dir([DPoae_rootDataDir HI.DataDir filesep '*dpoae*']);
        HI.dpFile= [DPoae_rootDataDir HI.DataDir filesep HI.dpFile(1).name];
        HI.calibFile= helper.get_lower_calibFile(HI.dpFile);
    end
    
    % Pre-exposure
    run(NH.calibFile);
    NH.calibData=ans;
    NH.calibData=NH.calibData.CalibData;
    %     run(NH.dpFile);
    %     NH.dpData= ans;
    out_DPOAE_data= helper.my_dpoae_analysis(NH.dpFile);
    NH.dpData=[out_DPOAE_data.dp_amp];
    NH.freqs_Hz= [out_DPOAE_data.freq2];
    calib_at_freqs=0*NH.freqs_Hz;
    for freqVar=1:length(NH.freqs_Hz)
        calib_at_freqs(freqVar)= helper.CalibInterp(NH.freqs_Hz(freqVar)/1e3, NH.calibData);
    end
    
    if existHI
        % post-exposure
        run(HI.calibFile);
        HI.calibData=ans;
        HI.calibData=HI.calibData.CalibData;
        %     run(HI.dpFile);
        %     HI.dpData= ans;
        
        out_DPOAE_data= helper.my_dpoae_analysis(HI.dpFile);
        HI.dpData=[out_DPOAE_data.dp_amp];
        HI.freqs_Hz= [out_DPOAE_data.freq2];
        
        calib_at_freqs=0*HI.freqs_Hz;
        for freqVar=1:length(HI.freqs_Hz)
            calib_at_freqs(freqVar)= helper.CalibInterp(HI.freqs_Hz(freqVar)/1e3, HI.calibData);
        end
    end
   
    %%
    hold on;
    
    if chinID~=369
        plot(NH.freqs_Hz/1e3, NH.dpData, '-o', 'color', co_struct.b, 'markersize', plt.mrkSize, 'linew', plt.lw1);
        plot(HI.freqs_Hz/1e3, HI.dpData, '-d', 'color', co_struct.r, 'markersize', plt.mrkSize, 'linew', plt.lw1);
    else
        plot(NH.freqs_Hz/1e3, NH.dpData, '-o', 'color', co_struct.lb, 'markersize', plt.mrkSize, 'linew', plt.lw1);
        if existHI
            plot(HI.freqs_Hz/1e3, HI.dpData, '-d', 'color', co_struct.lr, 'markersize', plt.mrkSize, 'linew', plt.lw1);
        end
    end
    if numel(NH.dpData)~=27
        fprintf('chin %d\n', chinID);
    end
    dp_data_nh(:, chinVar)= interp1(NH.freqs_Hz, NH.dpData, freq27);
    dp_data_hi(:, chinVar)= interp1(HI.freqs_Hz, HI.dpData, freq27);
end

plot(NH.freqs_Hz/1e3, nanmean(dp_data_nh, 2), '-', 'color', 'b', 'markersize', plt.mrkSize, 'linew', plt.lw3);
plot(HI.freqs_Hz/1e3, nanmean(dp_data_hi, 2), '-', 'color', 'r', 'markersize', plt.mrkSize, 'linew', plt.lw3);

set(gca, 'xscale', 'log', 'fontsize', plt.fontSize, 'xtick', xTicks/1e3, 'TickLength', plt.tick_len, 'linew', plt.ax_lw,'FontName',plt.fontName);
ylim([0 50]);
grid on;
ylabel('DP Amplitude (dB SPL)');
% xlabel('Frequency (kHz)');
xlim([450 10.1e3]/1e3);
box off;

lg(1)= plot(nan, nan, '-o', 'color', co_struct.b,'markersize', plt.mrkSize, 'linew', plt.lw1);
lg(2)= plot(nan, nan, '-d', 'color', co_struct.r,'markersize', plt.mrkSize, 'linew', plt.lw1);
[lg, icons]= legend(lg, 'NH', 'HI', 'Location', 'southeast', 'box', 'off');
grid off;

icons(3).XData= mean(icons(3).XData) + [-.15 +.15];
icons(5).XData= mean(icons(5).XData) + [-.15 +.15];

linkaxes(sp_ax, 'y');

txtHan= helper.add_subplot_letter(1, 2, 11);

rmpath(DirStruct.Codes);

%% Set axes placement/size
Xwidth=.39;
Xcorner=.085;
Xshift=.12;
Ywidth=.73;
Ycorner=.18;
% Yshift=0.06;

% A
set(sp_ax(1),'Position',[Xcorner Ycorner Xwidth Ywidth])
set(txtHan(1),'pos',get(txtHan(1),'pos')+[0 0.01 0],'FontName',plt.fontName)
drawnow
% B
set(sp_ax(2),'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
set(txtHan(2),'pos',get(txtHan(2),'pos')+[0 0.01 0],'FontName',plt.fontName)
drawnow

set(xlab_han, 'pos', [14 -6.736 -1]);
drawnow;

fName= 'Fig3_pooled_abr_dpoae';
if saveFigs
   saveas(gcf, [DirStruct.Figures.eps fName], 'epsc');
   saveas(gcf, [DirStruct.Figures.png fName], 'png');
end
