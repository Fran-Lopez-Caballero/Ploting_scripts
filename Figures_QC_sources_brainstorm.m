%% Define variables for QC of sources in old protocols

% Time series MAG, Time series GRAD, topos of each, brain left & right
% and source waveforms at A1 in left and right (8 figures, 4 columns, full
% page)

clear
root_dir = 'C:/Project'; % C:/Project or 'private/path/Project'
create_personalized_scouts = 1; % If that subject already has it, 0
Scout_to_check = 'A1'; % Only one at a time, and it will display left and right 
%{'52_L';'A1_L';'A4_L';'A5_L';'LBelt_L';'MBelt_L';'OP4_L';'PBelt_L';'RI_L';'STSdp_L';'STSvp_L';'52_R';'A1_R';'A4_R';'A5_R';'LBelt_R';'MBelt_R';'OP4_R';'PBelt_R';'RI_R';'STSdp_R';'STSvp_R'};
wave = {'LLR'}; % change to {'LLR', 'MLR'} if wanting to do both
modality_data = 'EEG'; % 'EEG' 'MEG' 'both_mod' % Except plot combined
% results_average_210208_2035_both_mod_Source_diag_LLR_Quietest
zoom_sources_6 = 4; zoom_sources = 2.5;
zoom_sources_all_modalities = 2; % For EEG + MEG + both
time_to_display = 0.056; % in seconds; N1 90, P2 150; Pa 30, P50 60 aprox
source_noise_tag = 'Regul'; % 'Regul' 'Eigen' 'diag'
WHAT_TO_PLOT = 6; 
% 1 = Subject & 93 only; 
% 2 = Subject & all; 
% 3 = Session & 93 only; 
% 4 = Session & all
% 5 = Same as 3 but MEG only and corregistration plot instead (to see
% effect of corregistration); 
% 6 = sources only of 6 sessions
% 7 = TO BE DONE: BUTTERFLY PLOT OF SOURCE WAVEFORM ACROSS SUBJECTS?
% 8 = EEG + MEG + BIMODAL average across sessions
% 9 = MAKE ANOTHER ONE IN THE PLOTS_EEG_MEG WITH SENSOR LOCAL MAXIMA (EEG MEG): ALL WAVEFORMS IN THE SAME PLOT
participant = {
    'xxxx' 
};
condition = {'11' '12' '13' '31' '32' '33' '51' '52' '53' '71' '72' '73' '91' '92' '93'};
load([root_dir '/Brainstorm_pipelines/session_block_array_sources.mat']) % Edited
plots_left_right_names = {'52_L','A1_L','A4_L','A5_L','LBelt_L','MBelt_L','OP4_L','PBelt_L','RI_L','STSdp_L','STSvp_L','52_R','A1_R','A4_R','A5_R','LBelt_R','MBelt_R','OP4_R','PBelt_R','RI_R','STSdp_R','STSvp_R'};
ScoutsArg = find(contains(plots_left_right_names,Scout_to_check));

%% Create personalized atlas for each subject (if not already there: Done for 26)

% If willing to create it in the high-defin cortex, select it as main in
% gui and run this again (for now useless since most subjects are looked at
% with the low resolution one)

if create_personalized_scouts
title_high_def = 'Black_color';
high_def_list = {'L_A1_ROI L';'L_LBelt_ROI L';'L_MBelt_ROI L';'R_A1_ROI R';'R_LBelt_ROI R';'R_MBelt_ROI R'};
high_def_new_names = {'A1_L';'LBelt_L';'MBelt_L';'A1_R';'LBelt_R';'MBelt_R'};
high_def_new_colors = {[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]}; % all black

title_plots_left_right = 'Plots_left_right';
plots_left_right_list = {'L_52_ROI L';'L_A1_ROI L';'L_A4_ROI L';'L_A5_ROI L';'L_LBelt_ROI L';'L_MBelt_ROI L';'L_OP4_ROI L';'L_PBelt_ROI L';'L_RI_ROI L';'L_STSdp_ROI L';'L_STSvp_ROI L';'R_52_ROI R';'R_A1_ROI R';'R_A4_ROI R';'R_A5_ROI R';'R_LBelt_ROI R';'R_MBelt_ROI R';'R_OP4_ROI R';'R_PBelt_ROI R';'R_RI_ROI R';'R_STSdp_ROI R';'R_STSvp_ROI R'};
plots_left_right_names = {'52_L';'A1_L';'A4_L';'A5_L';'LBelt_L';'MBelt_L';'OP4_L';'PBelt_L';'RI_L';'STSdp_L';'STSvp_L';'52_R';'A1_R';'A4_R';'A5_R';'LBelt_R';'MBelt_R';'OP4_R';'PBelt_R';'RI_R';'STSdp_R';'STSvp_R'};
plots_left_right_colors = {[0,51,102], [0,0,255], [102,0,204], [153,51,255],...
[0,0,204], [51,51,255], [0,153,153], [0,0,153], [204,0,204],...
[25,0,51], [51,0,51], [102,51,0], [255,0,0], [204,102,0], [255,153,51],...
[204,0,0], [255,51,51], [0,153,0], [153,0,0], [204,204,0],...
[51,25,0], [51,51,0]};

% Make the Left and right one 
for p = 1:length(participant)
    if ~exist([root_dir '/anat/' participant{p} '/brainstormsubject.mat'],'file');continue;end
    load([root_dir '/anat/' participant{p} '/brainstormsubject.mat'])
    load([root_dir '/anat/' Cortex]) % contains the atlas
    pos_new = find(strcmp({Atlas.Name},title_plots_left_right));
    if isempty(pos_new) % We have to create this atlas
        Pos_main_atlas = find(strcmp({Atlas.Name},'HCPMMP1'));
        % Ojo que puede ser hcpmmp1 también
        if isempty(Pos_main_atlas)
            Pos_main_atlas = find(strcmp({Atlas.Name},'hcpmmp1'));
        end
        if isempty(Pos_main_atlas)
            error(['no atlas under HCPMMP1 or hcpmmp1 names found for ' participant{p}])
        end
        % Crate new entry in Atlas with the name chosen
        Pos_new_atlas = size(Atlas,2)+1;
        Atlas(Pos_new_atlas).Name = title_plots_left_right; 
        % For every scout selected for this atlas, add it to the new entry
        for nal = 1:length(plots_left_right_list)
            Pos_orig = find(strcmp({Atlas(Pos_main_atlas).Scouts.Label},plots_left_right_list{nal}));
            Atlas(Pos_new_atlas).Scouts(nal) = Atlas(Pos_main_atlas).Scouts(Pos_orig);
            Atlas(Pos_new_atlas).Scouts(nal).Label = plots_left_right_names{nal};
            Atlas(Pos_new_atlas).Scouts(nal).Color = plots_left_right_colors{nal}/255;
        end
        % Make it the default atlas
        iAtlas = Pos_new_atlas;
        % Identify variables in file to be sure we save back the same
        variableInfo = who('-file',[root_dir '/anat/' Cortex]);
        save([root_dir '/anat/' Cortex],variableInfo{:})
    end
end

% Make the Black_color one (in this order so that Black_color one is the one defined by default)
for p = 1:length(participant)
    if ~exist([root_dir '/anat/' participant{p} '/brainstormsubject.mat'],'file');continue;end
    load([root_dir '/anat/' participant{p} '/brainstormsubject.mat'])    
    load([root_dir '/anat/' Cortex]) % contains the atlas
    pos_new = find(strcmp({Atlas.Name},title_high_def));
    if isempty(pos_new) % We have to create this atlas
        Pos_main_atlas = find(strcmp({Atlas.Name},'HCPMMP1'));
        % Ojo que puede ser hcpmmp1 también
        if isempty(Pos_main_atlas)
            Pos_main_atlas = find(strcmp({Atlas.Name},'hcpmmp1'));
        end
        if isempty(Pos_main_atlas)
            error(['no atlas under HCPMMP1 or hcpmmp1 names found for ' participant{p}])
        end
        % Crate new entry in Atlas with the name chosen
        Pos_new_atlas = size(Atlas,2)+1;
        Atlas(Pos_new_atlas).Name = title_high_def; 
        % For every scout selected for this atlas, add it to the new entry
        for nal = 1:length(high_def_list)
            Pos_orig = find(strcmp({Atlas(Pos_main_atlas).Scouts.Label},high_def_list{nal}));
            Atlas(Pos_new_atlas).Scouts(nal) = Atlas(Pos_main_atlas).Scouts(Pos_orig);
            Atlas(Pos_new_atlas).Scouts(nal).Label = high_def_new_names{nal};
            Atlas(Pos_new_atlas).Scouts(nal).Color = high_def_new_colors{nal}/255;
        end
        % Make it the default atlas
        iAtlas = Pos_new_atlas;
        % Identify variables in file to be sure we save back the same
        variableInfo = who('-file',[root_dir '/anat/' Cortex]);
        save([root_dir '/anat/' Cortex],variableInfo{:})
    end
end

% Reload subject
for p = 1:length(participant)
    try
        prot_subs = bst_get('ProtocolSubjects');
        if strcmp(participant{p},'@default_subject') % default anatomy
            current_sub = 0;
        else % Any subject
            current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
        end
        db_reload_subjects(current_sub);
    catch
        disp(['Please reload anatomy for participant ' participant{p} ' to see the scouts'])
    end
end

end

%% 1) Plot sensor data, topo, source location and scouts for each SUBJECT (Only 93 condition)
if WHAT_TO_PLOT == 1
for w = 1:length(wave)
if strcmp(wave{w},'LLR')
    root_dir_bs = [root_dir '/brainstorm_db/LLRseg']; 
elseif strcmp(wave{w},'MLR')
    root_dir_bs = [root_dir '/brainstorm_db/MLRseg']; 
end

for p = 1:length(participant)
    if ~exist([root_dir '/anat/' participant{p} '/brainstormsubject.mat'],'file');continue;end
    load([root_dir '/anat/' participant{p} '/brainstormsubject.mat'])
    SurfaceFile = Cortex; % In case we want to use something else in the future, like a whole brain
    
    for c = 15
        sensor_file_EEG = dir([root_dir_bs '/data/' participant{p} '/@intra/data*_Normal_EEG_' wave{w} '_' condition{c} '.mat']);
        if isempty(sensor_file_EEG); continue; end
        sensor_file_EEG = sensor_file_EEG.name;
        sensor_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/data*_Normal_MEG_' wave{w} '_' condition{c} '.mat']);
        if isempty(sensor_file_MEG); continue; end
        sensor_file_MEG = sensor_file_MEG.name;
        source_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/results*_' modality_data '_Source_' source_noise_tag '_' wave{w} '_' condition{c} '.mat']);
        if isempty(source_file_MEG); continue; end
        source_file_MEG = source_file_MEG.name;

        sFilesEEG = [participant{p} '/@intra/' sensor_file_EEG];
        sFilesMEG = [participant{p} '/@intra/' sensor_file_MEG];

        % Plot sensor data and scalp topographies
        hFig = view_timeseries(sFilesEEG,'EEG',[],'NewFigure'); %#ok<*NASGU>
        hFig = view_topography(sFilesEEG,'EEG','2DDisc',[],[],'NewFigure',[]);
        hFig = view_timeseries(sFilesMEG,'MEG',[],'NewFigure');
        hFig = view_topography(sFilesMEG,'MEG','2DDisc',[],[],'NewFigure',[]);

        % Plot sources
        OverlayFile = [participant{p} '/@intra/' source_file_MEG];
        view_surface_data(SurfaceFile, OverlayFile, []);

        % Plot source waveforms
        ResultsFiles{1,1} = [participant{p} '/@intra/' source_file_MEG];
        hFig = view_scouts(ResultsFiles,ScoutsArg,[]);
        pause
        close all
    end
end

end
end

%% 2) Plot sensor data, topo, source location and scouts for each SUBJECT AND CONDITION

if WHAT_TO_PLOT == 2
for w = 1:length(wave)
if strcmp(wave{w},'LLR')
    root_dir_bs = [root_dir '/brainstorm_db/LLRseg']; 
elseif strcmp(wave{w},'MLR')
    root_dir_bs = [root_dir '/brainstorm_db/MLRseg']; 
end

for p = 1:length(participant)
      
    if ~exist([root_dir '/anat/' participant{p} '/brainstormsubject.mat'],'file');continue;end
    load([root_dir '/anat/' participant{p} '/brainstormsubject.mat'])
    SurfaceFile = Cortex; % In case we want to use something else in the future, like a whole brain
   
    for c = 1:length(condition)
        sensor_file_EEG = dir([root_dir_bs '/data/' participant{p} '/@intra/data*_Normal_EEG_' wave{w} '_' condition{c} '.mat']);
        if isempty(sensor_file_EEG); continue; end
        sensor_file_EEG = sensor_file_EEG.name;
        sensor_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/data*_Normal_MEG_' wave{w} '_' condition{c} '.mat']);
        if isempty(sensor_file_MEG); continue; end
        sensor_file_MEG = sensor_file_MEG.name;
        source_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/results*_' modality_data '_Source_' source_noise_tag '_' wave{w} '_' condition{c} '.mat']);
        if isempty(source_file_MEG); continue; end
        source_file_MEG = source_file_MEG.name;

        sFilesEEG = [participant{p} '/@intra/' sensor_file_EEG];
        sFilesMEG = [participant{p} '/@intra/' sensor_file_MEG];

        % Plot sensor data and scalp topographies
        hFig = view_timeseries(sFilesEEG,'EEG',[],'NewFigure'); %#ok<*NASGU>
        hFig = view_topography(sFilesEEG,'EEG','2DDisc',[],[],'NewFigure',[]);
        hFig = view_timeseries(sFilesMEG,'MEG',[],'NewFigure');
        hFig = view_topography(sFilesMEG,'MEG','2DDisc',[],[],'NewFigure',[]);

        % Plot sources
        OverlayFile = [participant{p} '/@intra/' source_file_MEG];
        view_surface_data(SurfaceFile, OverlayFile, []);

        % Plot source waveforms
        ResultsFiles{1,1} = [participant{p} '/@intra/' source_file_MEG];
        hFig = view_scouts(ResultsFiles,ScoutsArg,[]);
        pause
        close all
    end
end

end
end

%% 3) Plot sensor data, topo, source location and scouts for each SESSION (Only 93 condition)
if WHAT_TO_PLOT == 3
for w = 1:length(wave)
if strcmp(wave{w},'LLR')
    root_dir_bs = [root_dir '/brainstorm_db/LLRseg']; 
elseif strcmp(wave{w},'MLR')
    root_dir_bs = [root_dir '/brainstorm_db/MLRseg']; 
end

for p = 1:length(participant)
      
    if ~exist([root_dir '/anat/' participant{p} '/brainstormsubject.mat'],'file');continue;end
    load([root_dir '/anat/' participant{p} '/brainstormsubject.mat'])
    SurfaceFile = Cortex; % In case we want to use something else in the future, like a whole brain
   
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        for c = 15
            sensor_file_EEG = dir([root_dir_bs '/data/' participant{p} '/@intra/data*_Normal_EEG_' wave{w} '_' condition{c} session{s} '.mat']);
            if isempty(sensor_file_EEG); continue; end
            sensor_file_EEG = sensor_file_EEG.name;
            sensor_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/data*_Normal_MEG_' wave{w} '_' condition{c} session{s} '.mat']);
            if isempty(sensor_file_MEG); continue; end
            sensor_file_MEG = sensor_file_MEG.name;
            source_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/results*_' modality_data '_Source_' source_noise_tag '_' wave{w} '_' condition{c} session{s} '.mat']);
            if isempty(source_file_MEG); continue; end
            source_file_MEG = source_file_MEG.name;
            
            sFilesEEG = [participant{p} '/@intra/' sensor_file_EEG];
            sFilesMEG = [participant{p} '/@intra/' sensor_file_MEG];
            
            % Plot sensor data and scalp topographies
            hFig = view_timeseries(sFilesEEG,'EEG',[],'NewFigure'); %#ok<*NASGU>
            hFig = view_topography(sFilesEEG,'EEG','2DDisc',[],[],'NewFigure',[]);
            hFig = view_timeseries(sFilesMEG,'MEG',[],'NewFigure');
            hFig = view_topography(sFilesMEG,'MEG','2DDisc',[],[],'NewFigure',[]);
            
            % Plot sources
            OverlayFile = [participant{p} '/@intra/' source_file_MEG];
            view_surface_data(SurfaceFile, OverlayFile, []);
            
            % Plot source waveforms
            ResultsFiles{1,1} = [participant{p} '/@intra/' source_file_MEG];
            hFig = view_scouts(ResultsFiles,ScoutsArg,[]);
            pause
            close all
        end
    end
end

end
end

%% 4) Plot sensor data, topo, source location and scouts for each SESSION AND CONDITION
if WHAT_TO_PLOT == 4
for w = 1:length(wave)
if strcmp(wave{w},'LLR')
    root_dir_bs = [root_dir '/brainstorm_db/LLRseg']; 
elseif strcmp(wave{w},'MLR')
    root_dir_bs = [root_dir '/brainstorm_db/MLRseg']; 
end

for p = 1:length(participant)
      
    if ~exist([root_dir '/anat/' participant{p} '/brainstormsubject.mat'],'file');continue;end
    load([root_dir '/anat/' participant{p} '/brainstormsubject.mat'])
    SurfaceFile = Cortex; % In case we want to use something else in the future, like a whole brain
   
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        for c = 1:length(condition)
            sensor_file_EEG = dir([root_dir_bs '/data/' participant{p} '/@intra/data*_Normal_EEG_' wave{w} '_' condition{c} session{s} '.mat']);
            if isempty(sensor_file_EEG); continue; end
            sensor_file_EEG = sensor_file_EEG.name;
            sensor_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/data*_Normal_MEG_' wave{w} '_' condition{c} session{s} '.mat']);
            if isempty(sensor_file_MEG); continue; end
            sensor_file_MEG = sensor_file_MEG.name;
            source_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/results*_' modality_data '_Source_' source_noise_tag '_' wave{w} '_' condition{c} session{s} '.mat']);
            if isempty(source_file_MEG); continue; end
            source_file_MEG = source_file_MEG.name;
            
            sFilesEEG = [participant{p} '/@intra/' sensor_file_EEG];
            sFilesMEG = [participant{p} '/@intra/' sensor_file_MEG];
            
            % Plot sensor data and scalp topographies
            hFig = view_timeseries(sFilesEEG,'EEG',[],'NewFigure'); %#ok<*NASGU>
            hFig = view_topography(sFilesEEG,'EEG','2DDisc',[],[],'NewFigure',[]);
            hFig = view_timeseries(sFilesMEG,'MEG',[],'NewFigure');
            hFig = view_topography(sFilesMEG,'MEG','2DDisc',[],[],'NewFigure',[]);
            
            % Plot sources
            OverlayFile = [participant{p} '/@intra/' source_file_MEG];
            view_surface_data(SurfaceFile, OverlayFile, []);
            
            % Plot source waveforms
            ResultsFiles{1,1} = [participant{p} '/@intra/' source_file_MEG];
            hFig = view_scouts(ResultsFiles,ScoutsArg,[]);
            pause
            close all
        end
    end
end

end
end

%% 5) Plot sensor data MEG with corregistrattion (SESSION and only 93 condition)
if WHAT_TO_PLOT == 5
for w = 1:length(wave)
if strcmp(wave{w},'LLR')
    root_dir_bs = [root_dir '/brainstorm_db/LLRseg']; 
elseif strcmp(wave{w},'MLR')
    root_dir_bs = [root_dir '/brainstorm_db/MLRseg']; 
end

for p = 1:length(participant)
      
    if ~exist([root_dir '/anat/' participant{p} '/brainstormsubject.mat'],'file');continue;end
    load([root_dir '/anat/' participant{p} '/brainstormsubject.mat'])
    SurfaceFile = Cortex; % In case we want to use something else in the future, like a whole brain
    
    % Identify which head mask to use
    if ~exist([root_dir '/anat/' Scalp],'file');continue;end
    head_surface_file = Scalp;
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for c = 15
            sensor_file_EEG = dir([root_dir_bs '/data/' participant{p} '/@intra/data*_Normal_EEG_' wave{w} '_' condition{c} session{s} '.mat']);
            if isempty(sensor_file_EEG); continue; end
            sensor_file_EEG = sensor_file_EEG.name;
            sensor_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/data*_Normal_MEG_' wave{w} '_' condition{c} session{s} '.mat']);
            if isempty(sensor_file_MEG); continue; end
            sensor_file_MEG = sensor_file_MEG.name;
            source_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/results*_' modality_data '_Source_' source_noise_tag '_' wave{w} '_' condition{c} session{s} '.mat']);
            if isempty(source_file_MEG); continue; end
            source_file_MEG = source_file_MEG.name;
            channel_file = [root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{1} '_normal/channel_vectorview306_acc1.mat'];
            if ~exist(channel_file,'file')
                continue; end
            channel_file = [participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{1} '_normal/channel_vectorview306_acc1.mat'];
            sFilesEEG = [participant{p} '/@intra/' sensor_file_EEG];
            sFilesMEG = [participant{p} '/@intra/' sensor_file_MEG];
            
            % Plot sensor data and scalp topographies
%             hFig = view_timeseries(sFilesEEG,'EEG',[],'NewFigure'); %#ok<*NASGU>
%             hFig = view_topography(sFilesEEG,'EEG','2DDisc',[],[],'NewFigure',[]);
            hFig = view_timeseries(sFilesMEG,'MEG',[],'NewFigure');
            hFig = view_topography(sFilesMEG,'MEG','2DDisc',[],[],'NewFigure',[]);
            
            % Plot sources
            OverlayFile = [participant{p} '/@intra/' source_file_MEG];
            view_surface_data(SurfaceFile, OverlayFile, []);
%             
            % Plot source waveforms
            ResultsFiles{1,1} = [participant{p} '/@intra/' source_file_MEG];
            hFig = view_scouts(ResultsFiles,ScoutsArg,[]);
            
            % Plot corregistration EEG and MEG (plot Brian?)
            hFig = channel_align_manual(channel_file,'MEG',0,'scalp');
            hFig = channel_align_manual(channel_file,'EEG',0,'scalp');

            pause
            close all
        end
    end
end

end
end

%% 6) Plot average of all 15 conditions per session in one single plot per subject
if WHAT_TO_PLOT == 6
for w = 1:length(wave)
if strcmp(wave{w},'LLR')
    root_dir_bs = [root_dir '/brainstorm_db/LLRseg']; 
elseif strcmp(wave{w},'MLR')
    root_dir_bs = [root_dir '/brainstorm_db/MLRseg']; 
end

for p = 1:length(participant)
      
    if ~exist([root_dir '/anat/' participant{p} '/brainstormsubject.mat'],'file');continue;end
    load([root_dir '/anat/' participant{p} '/brainstormsubject.mat'])
    SurfaceFile = Cortex; % In case we want to use something else in the future, like a whole brain
    
    % Identify which head mask to use
    if ~exist([root_dir '/anat/' Scalp],'file');continue;end
    head_surface_file = Scalp;
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        source_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/results*_' modality_data '_Source_' source_noise_tag '_' wave{w} '_ALL_' session{s} '.mat']);
%         source_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/results*_MEG_Source_Regul_' wave{w} '_93'  session{s} '.mat']);
        if isempty(source_file_MEG); continue; end
        source_file_MEG = source_file_MEG.name;
        channel_file = [root_dir_bs '/data/' participant{p} '/' wave{w} '_93' session{s} '_' block{1} '_normal/channel_vectorview306_acc1.mat'];
        if ~exist(channel_file,'file'); continue; end
        channel_file = [participant{p} '/' wave{w} '_93' session{s} '_' block{1} '_normal/channel_vectorview306_acc1.mat'];

        % Plot sources
        OverlayFile = [participant{p} '/@intra/' source_file_MEG];
        hFig = view_surface_data(SurfaceFile, OverlayFile, []);
        zoom(hFig, zoom_sources_6);
        figure_3d('SetStandardView', hFig, {'left','right'});
        panel_time('SetCurrentTime', time_to_display);
        
        cd C:\Project\Brainstorm_pipelines
        panel_scout_Fran('UpdateScoutsDisplay', 'all'); 
%             
%         % Plot source waveforms
%         ResultsFiles{1,1} = [participant{p} '/@intra/' source_file_MEG];
%         hFig = view_scouts(ResultsFiles,ScoutsArg,[]);

        % Option to check for every session along with corregistration
        % Plot corregistration EEG and MEG (plot Brian?)
%         hFig = channel_align_manual(channel_file,'MEG',0,'scalp');
%         hFig = channel_align_manual(channel_file,'EEG',0,'scalp');
%         pause
    end
    pause
    close all
end

end
end

%% 7) Plot average of all 15 conditions per session (single plot per session + corregistration)
if WHAT_TO_PLOT == 7
for w = 1:length(wave)
if strcmp(wave{w},'LLR')
    root_dir_bs = [root_dir '/brainstorm_db/LLRseg']; 
elseif strcmp(wave{w},'MLR')
    root_dir_bs = [root_dir '/brainstorm_db/MLRseg']; 
end

for p = 1:length(participant)
      
    if strcmp(participant{p},'2288') % because this was computed with the 30K surface
        error('check if participant 2288 sources are computed with the right anatomy')
%         if ~exist([root_dir '/anat/' participant{p} '/tess_cortex_pial_high.mat'],'file');continue;end
%         SurfaceFile = [participant{p} '/tess_cortex_pial_high.mat'];
    else
    if ~exist([root_dir '/anat/' participant{p} '/tess_cortex_pial_low.mat'],'file');continue;end
    SurfaceFile = [participant{p} '/tess_cortex_pial_low.mat'];
    end
    
    % Identify which head mask to use
    if ~exist([root_dir '/anat/' Scalp],'file');continue;end
    head_surface_file = Scalp;
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        source_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/results*_' modality_data '_Source_' source_noise_tag '_' wave{w} '_ALL_' session{s} '.mat']);
%         source_file_MEG = dir([root_dir_bs '/data/' participant{p} '/@intra/results*_MEG_Source_Regul_' wave{w} '_93'  session{s} '.mat']);
        if isempty(source_file_MEG); continue; end
        source_file_MEG = source_file_MEG.name;
        channel_file = [root_dir_bs '/data/' participant{p} '/' wave{w} '_93' session{s} '_' block{1} '_normal/channel_vectorview306_acc1.mat'];
        if ~exist(channel_file,'file'); continue; end
        channel_file = [participant{p} '/' wave{w} '_93' session{s} '_' block{1} '_normal/channel_vectorview306_acc1.mat'];

        % Plot sources
        OverlayFile = [participant{p} '/@intra/' source_file_MEG];
        hFig = view_surface_data(SurfaceFile, OverlayFile, []);
        zoom(hFig, zoom_sources);
        figure_3d('SetStandardView', hFig, {'left','right'});
        panel_time('SetCurrentTime', time_to_display);
%             
%         % Plot source waveforms
%         ResultsFiles{1,1} = [participant{p} '/@intra/' source_file_MEG];
%         hFig = view_scouts(ResultsFiles,ScoutsArg,[]);

        % Option to check for every session along with corregistration
        % Plot corregistration EEG and MEG (plot Brian?)
        hFig = channel_align_manual(channel_file,'MEG',0,'scalp');
        hFig = channel_align_manual(channel_file,'EEG',0,'scalp');
        pause
        close all
    end
end

end
end

%% 8) Plot source location Average of all sessions per subject, separately for EEG, MEG and bimodal

if WHAT_TO_PLOT == 8
    
% Exclusively here, change modality data variable
modality_data = {'EEG', 'MEG', 'both_mod'}; % 'EEG' 'MEG' 'both_mod'    

for w = 1:length(wave)
    
if strcmp(wave{w},'LLR')
    root_dir_bs = [root_dir '/brainstorm_db/LLRseg']; 
elseif strcmp(wave{w},'MLR')
    root_dir_bs = [root_dir '/brainstorm_db/MLRseg']; 
end

    for p = 1:length(participant)

        if ~exist([root_dir '/anat/' participant{p} '/brainstormsubject.mat'],'file');continue;end
        
        load([root_dir '/anat/' participant{p} '/brainstormsubject.mat'])
        SurfaceFile = Cortex; % In case we want to use something else in the future, like a whole brain
    
        % Identify which head mask to use
        if ~exist([root_dir '/anat/' Scalp],'file');continue;end
        head_surface_file = Scalp;
        
        for mode = 1:length(modality_data)
        
        % Define variables
        source_file_plot = dir([root_dir_bs '/data/' participant{p} '/@intra/results*_' modality_data{mode} '_Source_' source_noise_tag '_' wave{w} '_ALL.mat']);
        if isempty(source_file_plot); continue; end
        source_file_plot = source_file_plot.name;
%         channel_file = [root_dir_bs '/data/' participant{p} '/' wave{w} '_93' session{1} '_' block{1} '_normal/channel_vectorview306_acc1.mat'];
%         if ~exist(channel_file,'file'); continue; end
%         channel_file = [participant{p} '/' wave{w} '_93' session{1} '_' block{1} '_normal/channel_vectorview306_acc1.mat'];

        % Plot sources
        OverlayFile = [participant{p} '/@intra/' source_file_plot];
        hFig = view_surface_data(SurfaceFile, OverlayFile, []);
        zoom(hFig, zoom_sources_all_modalities);
        figure_3d('SetStandardView', hFig, {'left','right'});
        panel_time('SetCurrentTime', time_to_display);
        
        cd C:\Project\Brainstorm_pipelines
        panel_scout_Fran('UpdateScoutsDisplay', 'all'); 
%             
%         % Plot source waveforms
%         ResultsFiles{1,1} = [participant{p} '/@intra/' source_file_MEG];
%         hFig = view_scouts(ResultsFiles,ScoutsArg,[]);

        % Option to check for every session along with corregistration
        % Plot corregistration EEG and MEG (plot Brian?)
%         hFig = channel_align_manual(channel_file,'MEG',0,'scalp');
%         hFig = channel_align_manual(channel_file,'EEG',0,'scalp');

        end
        pause
        close all
    end
end   
end
