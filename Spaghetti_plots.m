%% Create spaghetti plots with any variable in the Mega_variable (at choice)
clear;clc;close all
name_Mega_variable = 'Mega_variable_June_22'; % contains the subjects to plot
% Save options
save_fig = 0; % 1 = YES; 0 = NO
type_of_plot = 'all'; % '4timepoints', 'basvseach', 'basvs1followup', 'all'
if strcmp(type_of_plot,'basvs1followup') || strcmp(type_of_plot,'all')
    priority_follow_up = [3 4 2]; % e.g. 3 4 2 would mean '6mo' first, if not found, '1yr', if not '3mo'
    priority_labels = {'6mo','1yr','3mo'}; % Order based on previous variable (e.g. '6mo','1yr','3mo')
end
y_axis_increase_percent = 1.2; % To make y axis a % larger than larger value % 1 = same; 3
save_fig_path = 'C:/Project/User/Cross_studies_data/Figures/';
% Check list in variable
load('C:/Project/User/Cross_studies_data/Scripts/list_vars_conglomerate.mat');
load('C:/Project/User/Project_MMN_global/decrease_good_vars.mat');
% [indx,tf] = listdlg('ListString',list_vars);
% vars_to_plot = list_vars(indx);

if 1 == 1 % silly way to compress
matching_measures = {'AGE','SEX','PSES','YRSED','VOCAB_TS'};
neuropsych_measures = {'OVERALLTSCR','MATRIX_TS','FULL2IQ','SPEEDTSCR','ATT_VIGTSCR',...
    'WMTSCR', 'VERBTSCR','SANITM','SAPITM'};
functioning_measures = {'ROLECURR','ROLELOW','ROLEHIGH','SOCIALCURR','SOCIALLOW',...
    'SOCIALHIGH','SFS_WITHDRAW_RS','SFS_INTERACT_RS','SFS_RECREAT_RS','SFS_OCCUP_RS',...
    'SFS_IND_PERF_RS','SFS_IND_COMP_RS','SFS_PROSOC_RS'};
clinical_measures = {'PANSSP_RS','PANSSN_RS','PANSST_RS','PSYRATS_AUD_HALL','PSYRATS_DELUSIONS'};
composite_measures = {'SANSSAPS_Level7_AudHall','SANSSAPS_Level7_UnusPercBeh',...
    'SANSSAPS_Level7_Delusions','SANSSAPS_Level7_ThDis','SANSSAPS_Level7_Inattention',...
    'SANSSAPS_Level7_Inexpress','SANSSAPS_Level7_Apathy','SANSSAPS_Level4_RealityDis',...
    'SANSSAPS_Level4_ThDis','SANSSAPS_Level4_Inexpress','SANSSAPS_Level4_Apathy',...
    'PANSS_Affect','PANSS_Disorg','PANSS_Negative','PANSS_Positive','PANSS_Resistance',...
    'BPRS_Total','BPRS_Positive','BPRS_Negative','BPRS_DeprAnx','BPRS_ActMania','BPRS_HostSusp'};
auditory_threshold_measures = {'QT_L_1000Hz', 'QT_R_1000Hz', 'QT_L_1500Hz', 'QT_R_1500Hz',...
    'QT_L_2000Hz','QT_R_2000Hz', 'QT_L_3000Hz', 'QT_R_3000Hz', 'QT_L_4000Hz', 'QT_R_4000Hz',...
    'AUD_THRESH_L','AUD_THRESH_R'};
duration_of_illness_measures = {'DAYS_SINCE_PRDM','DAYS_SINCE_1STEP','DAYS_SINCE_DISOR','DUP','PSYCH2SCAN','MED2SCAN'};
medication_measures = {'CPZ_equivalent'};
Struct_measures_vol = {'vol_A1_L','vol_A1_R','vol_LBelt_L','vol_LBelt_R',...
    'vol_PBelt_L','vol_PBelt_R','vol_OFC_L_merged','vol_OFC_R_merged',...
    'vol_IFG_L','vol_IFG_R','vol_AUDCORTEX_L','vol_AUDCORTEX_R'};
Struct_measures_thickness = {'thck_A1_L','thck_A1_R','thck_LBelt_L','thck_LBelt_R',...
    'thck_PBelt_L','thck_PBelt_R'};
end

vars_to_plot = medication_measures;
% matching_measures, neuropsych_measures, functioning_measures,
% composite_measures, auditory_threshold_measures,
% duration_of_illness_measures, medication_measures
% Struct_measures_thickness, Struct_measures_vol

% vars_to_plot = {'vol_A1_L';'thck_A1_L';'vol_A1_R';'thck_A1_R';'vol_LBelt_L';'thck_LBelt_L';'vol_LBelt_R';'thck_LBelt_R';'vol_IFG_L';'vol_IFG_R';'vol_AUDCORTEX_L';'vol_AUDCORTEX_R'};
% vars_to_plot = {'pMMN_A1_L';'dMMN_A1_L';'pMMN_A1_R';'dMMN_A1_R';'pMMN_IFG_L';'dMMN_IFG_L';'pMMN_IFG_R';'dMMN_IFG_R';'pMMN_AUDCORTEX_L';'dMMN_AUDCORTEX_L';'pMMN_AUDCORTEX_R';'dMMN_AUDCORTEX_R'};

%% Longitudinal plots

if strcmp(type_of_plot,'4timepoints') || strcmp(type_of_plot,'all')
load(['C:/Project/User/Cross_studies_data/Mega_variable/Longitudinal/' name_Mega_variable '.mat']);
time_point_label = {'A','B','C','D'};
% To be line and marker colors are the same (so long as loop is < 42 i)
% load('C:/Project/User/Project_MMN_global/default_matlab_colors_42.mat')
def_colors = [];
for i = 1:500 % Create a color selection of 500
    def_colors(i,:) = rand(1,3);
end
% Replace all empty cells with nans
idx = cellfun('isempty',Base);
Base(idx) = {NaN};
for i = 2:size(Base,1) % All but header
    if Base{i,2} == 1 % C
        Base{i,2} = 'C';
    elseif Base{i,2} == 2 % FE
        Base{i,2} = 'FE';
    end
end
% FE subjects
subj_labels_FE = Base(strcmp(Base(:,2),'FE'),1);
% C subjects
subj_labels_C = Base(strcmp(Base(:,2),'C'),1);

% White background
set(0,'defaultfigurecolor',[1 1 1]); % I want white backgrounds

for vtp = 1:length(vars_to_plot)
    % New plot 
    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = subplot(1,2,1); h(2) = subplot(1,2,2);
    
    %%  Controls to the left
    spa_plot = []; % Just to count number of non NaNs later
    legend_C = {}; pos = 1;
    for p = 1:length(subj_labels_C)
        % Find row where this participant is in base variable
        row_pos = find(strcmp(Base(:,1),subj_labels_C{p}));    
        % Create vector
        subject_vector = [];
        for tp = 1:length(time_point_label)
            col_pos = find(strcmp(Base(1,:),[time_point_label{tp} '_' vars_to_plot{vtp}]));
            subject_vector(tp) = Base{row_pos,col_pos};
        end
        spa_plot(p,:) = subject_vector;
        % Find last position of NaN in vector for this subject
        last_nonnan = find(~isnan(subject_vector));
        if ~isempty(last_nonnan)
            last_nonnan = last_nonnan(end);
        elseif isempty(last_nonnan)
            warning('subject with all NaN values');
            continue;
        end
        % Build subject legend (to be sure allnan subjects are not included)
        legend_C{pos} = subj_labels_C{p};
        legend_C{pos+1} = subj_labels_C{p};
        pos = pos + 2;
        % Plot
        if strcmp(subj_labels_C{p},'2214')
            disp 'here'
        end
        hold (h(1),'on')
        idxs = ~isnan(subject_vector);
        x = 1:4;
        plot(h(1),x(idxs), subject_vector(idxs),'color',def_colors(p,:));
        plot(h(1),subject_vector,'-s','MarkerSize',10,'color',def_colors(p,:),'MarkerFaceColor',def_colors(p,:));
    end
   	xlim(h(1),[0.5,4.5]);
    xticks(h(1),[1,2,3,4]);
    legend(h(1),legend_C);
    pom=findobj('type','legend');
    delete(pom);
    set(h(1),'TickLength',[0 0])
    how_many = find(~isnan(spa_plot(:,1)));
    num_bas = length(how_many);
    how_many = find(~isnan(spa_plot(:,2)));
    num_3mo = length(how_many);
    how_many = find(~isnan(spa_plot(:,3)));
    num_6mo = length(how_many);
    how_many = find(~isnan(spa_plot(:,4)));
    num_1yr = length(how_many);
    xticklabels(h(1),{['Baseline (' num2str(num_bas) ')'],['3mo (' num2str(num_3mo) ')'],['6mo (' num2str(num_6mo) ')'],['1yr (' num2str(num_1yr) ')']});
    xtl = get(h(1),'XTickLabel');  
    set(h(1),'XTickLabel',xtl,'fontsize',14)
    % Get max and min from spa_plot
    max_c = max(spa_plot(:));
    min_c = min(spa_plot(:));
    % Var in y axis (only left plot)
    ylabel_title = strrep(vars_to_plot{vtp},'_',' ');
    ylabel(h(1),ylabel_title,'fontsize',14,'FontWeight','bold');
    % Title
    title(h(1),'Controls','FontSize',14);
    
    %%  FE to the right
    spa_plot = []; % Just to count number of non NaNs later
    legend_FE = {}; pos = 1;
    for p = 1:length(subj_labels_FE)
        % Find row where this participant is in base variable
        row_pos = find(strcmp(Base(:,1),subj_labels_FE{p}));    
        % Create vector
        subject_vector = [];
        for tp = 1:length(time_point_label)
            col_pos = find(strcmp(Base(1,:),[time_point_label{tp} '_' vars_to_plot{vtp}]));
            subject_vector(tp) = Base{row_pos,col_pos};
        end
        spa_plot(p,:) = subject_vector;
        % Find last position of NaN in vector for this subject
        last_nonnan = find(~isnan(subject_vector));
        if ~isempty(last_nonnan)
            last_nonnan = last_nonnan(end);
        elseif isempty(last_nonnan)
            warning('subject with all NaN values');
            continue;
        end
        % Build subject legend (to be sure allnan subjects are not included)
        legend_FE{pos} = subj_labels_FE{p};
        legend_FE{pos+1} = subj_labels_FE{p};
        pos = pos + 2;
        % Plot
        hold (h(2),'on')
        idxs = ~isnan(subject_vector);
        x = 1:4;
        plot(h(2),x(idxs), subject_vector(idxs),'color',def_colors(p,:));
        plot(h(2),subject_vector,'-s','MarkerSize',10,'color',def_colors(p,:),'MarkerFaceColor',def_colors(p,:));
    end
   	xlim(h(2),[0.5,4.5]);
    xticks(h(2),[1,2,3,4]);
    legend(h(2),legend_FE);
    pom=findobj('type','legend');
    delete(pom);
    set(h(1),'TickLength',[0 0])
    how_many = find(~isnan(spa_plot(:,1)));
    num_bas = length(how_many);
    how_many = find(~isnan(spa_plot(:,2)));
    num_3mo = length(how_many);
    how_many = find(~isnan(spa_plot(:,3)));
    num_6mo = length(how_many);
    how_many = find(~isnan(spa_plot(:,4)));
    num_1yr = length(how_many);
    xticklabels(h(2),{['Baseline (' num2str(num_bas) ')'],['3mo (' num2str(num_3mo) ')'],['6mo (' num2str(num_6mo) ')'],['1yr (' num2str(num_1yr) ')']});
    xtl = get(h(2),'XTickLabel');  
    set(h(2),'XTickLabel',xtl,'fontsize',14)
    % Get max and min from spa_plot
    max_fe = max(spa_plot(:));
    min_fe = min(spa_plot(:));
    % Title
    title(h(2),'First episodes','FontSize',14);
    
    %% No adjust y axis for both
    max_total = max(max_c,max_fe);
    min_total = min(min_c,min_fe);
    
    % Add a percentage of increase in y axis
    if min_total < 0
        min_y_axis = (abs(min_total)*y_axis_increase_percent) * -1;
    else
        min_y_axis = min_total/y_axis_increase_percent;
    end
    if max_total < 0
        max_y_axis = (abs(min_total)/y_axis_increase_percent) * -1;
    else
        max_y_axis = max_total*y_axis_increase_percent;
    end
    
    for i = 1:2 % Two plots
        ylim(h(i),[min_y_axis,max_y_axis]);
    end
    
    if save_fig == 1
        saveas(gcf,[save_fig_path 'Long_' vars_to_plot{vtp} '.emf']);
    end
end
end

%% Bas vs 3mo, Bas vs 6mo, Bas vs 1yr individually

if strcmp(type_of_plot,'basvseach') || strcmp(type_of_plot,'all')
time_point_label = {'A','B','C','D'};
time_comparisons = {'bas3mo','bas6mo','bas1yr'};
load(['C:/Project/User/Cross_studies_data/Mega_variable/Longitudinal/' name_Mega_variable '.mat']);
% To be line and marker colors are the same (so long as loop is < 42 i)
% load('C:/Project/User/Project_MMN_global/default_matlab_colors_42.mat') 
% Replace all empty cells with nans
idx = cellfun('isempty',Base);
Base(idx) = {NaN};
for i = 2:size(Base,1) % All but header
    if Base{i,2} == 1 % C
        Base{i,2} = 'C';
    elseif Base{i,2} == 2 % FE
        Base{i,2} = 'FE';
    end
end
% FE subjects
subj_labels_FE = Base(strcmp(Base(:,2),'FE'),1);
% C subjects
subj_labels_C = Base(strcmp(Base(:,2),'C'),1);

% White background
set(0,'defaultfigurecolor',[1 1 1]); % I want white backgrounds

for vtp = 1:length(vars_to_plot)
    % New plot 
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = subplot(1,3,1); h(2) = subplot(1,3,2); h(3) = subplot(1,3,3);
    
    %%  Controls to the left in every plot
    
    % To display Ns later
    Count_C_bas3mo = 0;
    Count_C_bas6mo = 0;
    Count_C_bas1yr = 0;
    All_C_vector_bas3mo = [];
    All_C_vector_bas6mo = [];
    All_C_vector_bas1yr = [];
    legend_C_bas3mo = {};
    legend_C_bas6mo = {};
    legend_C_bas1yr = {};
    for p = 1:length(subj_labels_C)
        % Find row where this participant is in base variable
        row_pos = find(strcmp(Base(:,1),subj_labels_C{p}));
        for tc = 1:length(time_comparisons)
            if strcmp(time_comparisons{tc},'bas3mo')
                time_points_to_search = [1,2]; % from time_point_label variable ('A','B','C','D')
                plot_pos = 1;
            elseif strcmp(time_comparisons{tc},'bas6mo')
                time_points_to_search = [1,3];
                plot_pos = 2;
            elseif strcmp(time_comparisons{tc},'bas1yr')
                time_points_to_search = [1,4];
                plot_pos = 3;
            end
            
            % Create vector
            subject_vector = [];pos_vector = 1;
            for tp = time_points_to_search % Baseline and 3mo
                col_pos = find(strcmp(Base(1,:),[time_point_label{tp} '_' vars_to_plot{vtp}]));
                subject_vector(pos_vector) = Base{row_pos,col_pos};
                pos_vector = pos_vector +1;
            end
            % Find last position of NaN in vector for this subject
            nan_value = find(isnan(subject_vector)); 
            if ~isempty(nan_value)
                continue; % If there is a NaN in a pair, don't plot it
            else
                % Count to display Ns later
                eval(['Count_C_' time_comparisons{tc} ' = Count_C_' time_comparisons{tc} '+1;'])
                % if size(All_C_vector_bas1yr,1) == 14
                %   disp 'here';
                % end
                % Store to determine y axis later
                eval(['pos_all_vector = size(All_C_vector_' time_comparisons{tc} ',1) + 1;'])
                eval(['All_C_vector_' time_comparisons{tc} '(pos_all_vector,:) = subject_vector;'])
                % Build subject legend (to be sure allnan subjects are not included)
                eval(['pos = size(legend_C_' time_comparisons{tc} ',2) + 1;'])
                eval(['legend_C_' time_comparisons{tc} '{pos} = subj_labels_C{p};'])
                eval(['legend_C_' time_comparisons{tc} '{pos+1} = subj_labels_C{p};'])
                % Plot
                hold (h(plot_pos),'on')
                plot(h(plot_pos),1:2,subject_vector,'k');
                plot(h(plot_pos),1:2,subject_vector,'-o','MarkerSize',5,'color',[0,0,0],'MarkerFaceColor',[0,0,0]);
            end
        end
    end
    
    %%  FE to the right in every plot
    
    % To display Ns later
    Count_FE_bas3mo = 0;
    Count_FE_bas6mo = 0;
    Count_FE_bas1yr = 0;
    All_FE_vector_bas3mo = [];
    All_FE_vector_bas6mo = [];
    All_FE_vector_bas1yr = [];
    legend_FE_bas3mo = {};
    legend_FE_bas6mo = {};
    legend_FE_bas1yr = {};
    for p = 1:length(subj_labels_FE)
        % Find row where this participant is in base variable
        row_pos = find(strcmp(Base(:,1),subj_labels_FE{p}));
        for tc = 1:length(time_comparisons)
            if strcmp(time_comparisons{tc},'bas3mo')
                time_points_to_search = [1,2]; % from time_point_label variable ('A','B','C','D')
                plot_pos = 1;
            elseif strcmp(time_comparisons{tc},'bas6mo')
                time_points_to_search = [1,3];
                plot_pos = 2;
            elseif strcmp(time_comparisons{tc},'bas1yr')
                time_points_to_search = [1,4];
                plot_pos = 3;
            end
            
            % Create vector
            subject_vector = [];pos_vector = 1;
            for tp = time_points_to_search % Baseline and 3mo
                col_pos = find(strcmp(Base(1,:),[time_point_label{tp} '_' vars_to_plot{vtp}]));
                subject_vector(pos_vector) = Base{row_pos,col_pos};
                pos_vector = pos_vector +1;
            end
            % Find last position of NaN in vector for this subject
            nan_value = find(isnan(subject_vector)); 
            if ~isempty(nan_value)
                continue; % If there is a NaN in a pair, don't plot it
            else
                % Count to display Ns later
                eval(['Count_FE_' time_comparisons{tc} ' = Count_FE_' time_comparisons{tc} '+1;'])
                % if size(All_C_vector_bas1yr,1) == 14
                %   disp 'here';
                % end
                % Store to determine y axis later
                eval(['pos_all_vector = size(All_FE_vector_' time_comparisons{tc} ',1) + 1;'])
                eval(['All_FE_vector_' time_comparisons{tc} '(pos_all_vector,:) = subject_vector;'])
                % Build subject legend (to be sure allnan subjects are not included)
                eval(['pos = size(legend_FE_' time_comparisons{tc} ',2) + 1;'])
                eval(['legend_FE_' time_comparisons{tc} '{pos} = subj_labels_FE{p};'])
                eval(['legend_FE_' time_comparisons{tc} '{pos+1} = subj_labels_FE{p};'])
                % Plot
                hold (h(plot_pos),'on')
                plot(h(plot_pos),3:4,subject_vector,'r');
                plot(h(plot_pos),3:4,subject_vector,'-o','MarkerSize',5,'color',[1,0,0],'MarkerFaceColor',[1,0,0]);
            end
        end
    end
    legend_3mo = [legend_C_bas3mo';legend_FE_bas3mo'];
    legend_6mo = [legend_C_bas6mo';legend_FE_bas6mo'];
    legend_1yr = [legend_C_bas1yr';legend_FE_bas1yr'];
    legend_cell = {'3mo','6mo','1yr'};
    % Get max and min
    Total_values = [All_C_vector_bas3mo; All_FE_vector_bas3mo;...
        All_C_vector_bas6mo; All_FE_vector_bas6mo;...
        All_C_vector_bas1yr; All_FE_vector_bas1yr];
    max_total = max(Total_values(:));
    min_total = min(Total_values(:));
    
    % Get mean values to plot
    mean_names = {'mean_C_bas_3mo_1','mean_C_bas_3mo_2','mean_FE_bas_3mo_1','mean_FE_bas_3mo_2',...
        'mean_C_bas_6mo_1','mean_C_bas_6mo_2','mean_FE_bas_6mo_1','mean_FE_bas_6mo_2',...
        'mean_C_bas_1yr_1','mean_C_bas_1yr_2','mean_FE_bas_1yr_1','mean_FE_bas_1yr_2'};
    attempts = {'mean(All_C_vector_bas3mo(:,1),1)','mean(All_C_vector_bas3mo(:,2),1)','mean(All_FE_vector_bas3mo(:,1),1)','mean(All_FE_vector_bas3mo(:,2),1)',...
        'mean(All_C_vector_bas6mo(:,1),1)','mean(All_C_vector_bas6mo(:,2),1)','mean(All_FE_vector_bas6mo(:,1),1)','mean(All_FE_vector_bas6mo(:,2),1)',...
        'mean(All_C_vector_bas1yr(:,1),1)','mean(All_C_vector_bas1yr(:,2),1)','mean(All_FE_vector_bas1yr(:,1),1)','mean(All_FE_vector_bas1yr(:,2),1)'};
    
    % Silly way to avoid it ocupying a large space (in case there are no
    % values e.g. for controls in PANSS)
    for mn = 1:length(mean_names)
        try 
            eval([mean_names{mn} ' = ' attempts{mn} ';'])
        catch
            eval([mean_names{mn} ' = NaN'])
        end
    end
    
    % Add a percentage of increase in y axis
    if min_total < 0
        min_y_axis = (abs(min_total)*y_axis_increase_percent) * -1;
    else
        min_y_axis = min_total/y_axis_increase_percent;
    end
    if max_total < 0
        max_y_axis = (abs(min_total)/y_axis_increase_percent) * -1;
    else
        max_y_axis = max_total*y_axis_increase_percent;
    end
        
    for i = 1:3 % Three subplots
        % Plot mean values
        if i == 1
              plot(h(1),0.5:1.5,[mean_C_bas_3mo_1,mean_C_bas_3mo_1],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
              plot(h(1),1.5:2.5,[mean_C_bas_3mo_2,mean_C_bas_3mo_2],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
              plot(h(1),2.5:3.5,[mean_FE_bas_3mo_1,mean_FE_bas_3mo_1],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
              plot(h(1),3.5:4.5,[mean_FE_bas_3mo_2,mean_FE_bas_3mo_2],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
        elseif i == 2
              plot(h(2),0.5:1.5,[mean_C_bas_6mo_1,mean_C_bas_6mo_1],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
              plot(h(2),1.5:2.5,[mean_C_bas_6mo_2,mean_C_bas_6mo_2],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
              plot(h(2),2.5:3.5,[mean_FE_bas_6mo_1,mean_FE_bas_6mo_1],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
              plot(h(2),3.5:4.5,[mean_FE_bas_6mo_2,mean_FE_bas_6mo_2],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
        elseif i == 3
              plot(h(3),0.5:1.5,[mean_C_bas_1yr_1,mean_C_bas_1yr_1],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
              plot(h(3),1.5:2.5,[mean_C_bas_1yr_2,mean_C_bas_1yr_2],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
              plot(h(3),2.5:3.5,[mean_FE_bas_1yr_1,mean_FE_bas_1yr_1],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
              plot(h(3),3.5:4.5,[mean_FE_bas_1yr_2,mean_FE_bas_1yr_2],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
        end
        
        xlim(h(i),[0.5,4.5]);
        xticks(h(i),[1.5,3.5]);
        eval(['count_1 = Count_C_' time_comparisons{i} ';'])
        eval(['count_2 = Count_FE_' time_comparisons{i} ';'])
        xticklabels(h(i),{['C (' num2str(count_1) ')'],['FE (' num2str(count_2) ')']});
        eval(['legend(h(i),legend_' legend_cell{i} ');'])
        pom=findobj('type','legend');
        delete(pom);
        set(h(i),'TickLength',[0 0])
        xtl = get(h(i),'XTickLabel');  
        set(h(i),'XTickLabel',xtl,'fontsize',14)
        ylim(h(i),[min_y_axis,max_y_axis]);
    end
    
    % Identify how many increase (C)
    decreasing_C_bas3mo = 0;
    for ar = 1:size(All_C_vector_bas3mo,1)
        if All_C_vector_bas3mo(ar,1) > All_C_vector_bas3mo(ar,2)
            decreasing_C_bas3mo = decreasing_C_bas3mo + 1;
        end
    end
    decreasing_C_bas6mo = 0;
    for ar = 1:size(All_C_vector_bas6mo,1)
        if All_C_vector_bas6mo(ar,1) > All_C_vector_bas6mo(ar,2)
            decreasing_C_bas6mo = decreasing_C_bas6mo + 1;
        end
    end
    decreasing_C_bas1yr = 0;
    for ar = 1:size(All_C_vector_bas1yr,1)
        if All_C_vector_bas1yr(ar,1) > All_C_vector_bas1yr(ar,2)
            decreasing_C_bas1yr = decreasing_C_bas1yr + 1;
        end
    end
    
    % Identify how many increase (FE)
    decreasing_FE_bas3mo = 0;
    for ar = 1:size(All_FE_vector_bas3mo,1)
        if All_FE_vector_bas3mo(ar,1) > All_FE_vector_bas3mo(ar,2)
            decreasing_FE_bas3mo = decreasing_FE_bas3mo + 1;
        end
    end
    decreasing_FE_bas6mo = 0;
    for ar = 1:size(All_FE_vector_bas6mo,1)
        if All_FE_vector_bas6mo(ar,1) > All_FE_vector_bas6mo(ar,2)
            decreasing_FE_bas6mo = decreasing_FE_bas6mo + 1;
        end
    end
    decreasing_FE_bas1yr = 0;
    for ar = 1:size(All_FE_vector_bas1yr,1)
        if All_FE_vector_bas1yr(ar,1) > All_FE_vector_bas1yr(ar,2)
            decreasing_FE_bas1yr = decreasing_FE_bas1yr + 1;
        end
    end
    
    % Get worse variable depends on which var_to_plot we are dealing with
    is_in_decrease_good = find(strcmp(decrease_good,vars_to_plot{vtp}));
    if ~isempty(is_in_decrease_good)
        get_worse_C_bas3mo = num2str(Count_C_bas3mo-decreasing_C_bas3mo); % That is, the worse ones are the ones that increase
        get_worse_FE_bas3mo = num2str(Count_FE_bas3mo-decreasing_FE_bas3mo);
        get_worse_C_bas6mo = num2str(Count_C_bas6mo-decreasing_C_bas6mo);
        get_worse_FE_bas6mo = num2str(Count_FE_bas6mo-decreasing_FE_bas6mo);
        get_worse_C_bas1yr = num2str(Count_C_bas1yr-decreasing_C_bas1yr);
        get_worse_FE_bas1yr = num2str(Count_FE_bas1yr-decreasing_FE_bas1yr);
    else
        get_worse_C_bas3mo = num2str(decreasing_C_bas3mo);
        get_worse_FE_bas3mo = num2str(decreasing_FE_bas3mo);
        get_worse_C_bas6mo = num2str(decreasing_C_bas6mo);
        get_worse_FE_bas6mo = num2str(decreasing_FE_bas6mo);
        get_worse_C_bas1yr = num2str(decreasing_C_bas1yr);
        get_worse_FE_bas1yr = num2str(decreasing_FE_bas1yr);
    end
    
    percentage_worse_C_bas3mo = round((str2double(get_worse_C_bas3mo)/Count_C_bas3mo)*100); % round((17/38)*100)
    percentage_worse_FE_bas3mo = round((str2double(get_worse_FE_bas3mo)/Count_FE_bas3mo)*100);
    percentage_worse_C_bas6mo = round((str2double(get_worse_C_bas6mo)/Count_C_bas6mo)*100); % round((17/38)*100)
    percentage_worse_FE_bas6mo = round((str2double(get_worse_FE_bas6mo)/Count_FE_bas6mo)*100);
    percentage_worse_C_bas1yr = round((str2double(get_worse_C_bas1yr)/Count_C_bas1yr)*100); % round((17/38)*100)
    percentage_worse_FE_bas1yr = round((str2double(get_worse_FE_bas1yr)/Count_FE_bas1yr)*100);
        
    % Annotate how many show decrease
    annotation('textbox', [0.132, 0.82, 0.1, 0.1], 'String', [get_worse_C_bas3mo '/' num2str(Count_C_bas3mo) ' worsen (' num2str(percentage_worse_C_bas3mo) '%)'],'LineStyle','none','FontSize',12);
    annotation('textbox', [0.255, 0.82, 0.1, 0.1], 'String', [get_worse_FE_bas3mo '/' num2str(Count_FE_bas3mo) ' (' num2str(percentage_worse_FE_bas3mo) '%)'],'LineStyle','none','FontSize',14);
    annotation('textbox', [0.425, 0.82, 0.1, 0.1], 'String', [get_worse_C_bas6mo '/' num2str(Count_C_bas6mo) ' (' num2str(percentage_worse_C_bas6mo) '%)'],'LineStyle','none','FontSize',14);
    annotation('textbox', [0.535, 0.82, 0.1, 0.1], 'String', [get_worse_FE_bas6mo '/' num2str(Count_FE_bas6mo) ' (' num2str(percentage_worse_FE_bas6mo) '%)'],'LineStyle','none','FontSize',14);
    annotation('textbox', [0.708, 0.82, 0.1, 0.1], 'String', [get_worse_C_bas1yr '/' num2str(Count_C_bas1yr) ' (' num2str(percentage_worse_C_bas1yr) '%)'],'LineStyle','none','FontSize',14);
    annotation('textbox', [0.813, 0.82, 0.1, 0.1], 'String', [get_worse_FE_bas1yr '/' num2str(Count_FE_bas1yr) ' (' num2str(percentage_worse_FE_bas1yr) '%)'],'LineStyle','none','FontSize',14); % 'FontSize',14
    
    % Var in y axis (only left plot)
    ylabel_title = strrep(vars_to_plot{vtp},'_',' ');
    ylabel(h(1),ylabel_title,'fontsize',14,'FontWeight','bold');
    % Title
    title(h(1),'Baseline vs 3mo','FontSize',14);
    title(h(2),'Baseline vs 6mo','FontSize',14);
    title(h(3),'Baseline vs 1yr','FontSize',14);
    
    if save_fig == 1
        saveas(gcf,[save_fig_path 'Each_' vars_to_plot{vtp} '.emf']);
    end
end
end

%% Baseline vs one follow up

if strcmp(type_of_plot,'basvs1followup') || strcmp(type_of_plot,'all')

time_point_label = {'A','B','C','D'};
time_comparisons = {'bas3mo','bas6mo','bas1yr'};
load(['C:/Project/User/Cross_studies_data/Mega_variable/Longitudinal/' name_Mega_variable '.mat']);
% To be line and marker colors are the same (so long as loop is < 42 i)
% load('C:/Project/User/Project_MMN_global/default_matlab_colors_42.mat') 
% Replace all empty cells with nans
idx = cellfun('isempty',Base);
Base(idx) = {NaN};
for i = 2:size(Base,1) % All but header
    if Base{i,2} == 1 % C
        Base{i,2} = 'C';
    elseif Base{i,2} == 2 % FE
        Base{i,2} = 'FE';
    end
end
% FE subjects
subj_labels_FE = Base(strcmp(Base(:,2),'FE'),1);
% C subjects
subj_labels_C = Base(strcmp(Base(:,2),'C'),1);

% White background
set(0,'defaultfigurecolor',[1 1 1]); % I want white backgrounds

for vtp = 1:length(vars_to_plot)
    % New plot 
%     fig = figure('units','normalized','outerposition',[0 0 1 1]);
    fig = figure;
    
    %%  Controls to the left in single plot
    
    % To display Ns later
    count_C_priority_1 = 0;
    count_C_priority_2 = 0;
    count_C_priority_3 = 0;
    All_C_vector = [];
    legend_C = {};
    for p = 1:length(subj_labels_C)
        % Find row where this participant is in base variable
        row_pos = find(strcmp(Base(:,1),subj_labels_C{p}));
        % Create vector
        subject_vector = [];pos_vector = 1;
        for tp = [1,priority_follow_up(1)] % Baseline and priority 1
            col_pos = find(strcmp(Base(1,:),[time_point_label{tp} '_' vars_to_plot{vtp}]));
            if tp ~= 1 % so if we are in the follow-up
               if isnan(Base{row_pos,col_pos}) % The priority follow up is empty
                   col_pos = find(strcmp(Base(1,:),[time_point_label{priority_follow_up(2)} '_' vars_to_plot{vtp}]));
                   if isnan(Base{row_pos,col_pos}) % If it's still empty, go to priority number 3
                       col_pos = find(strcmp(Base(1,:),[time_point_label{priority_follow_up(3)} '_' vars_to_plot{vtp}]));
                       if isnan(Base{row_pos,col_pos}) % If still empty, this subject has no follow up, so don't plot
                           continue;
                       else
                           count_C_priority_3 = count_C_priority_3 + 1;
                       end
                   else
                       count_C_priority_2 = count_C_priority_2 + 1;
                   end
               else
                   count_C_priority_1 = count_C_priority_1 + 1;
               end
            end
            subject_vector(pos_vector) = Base{row_pos,col_pos};
            pos_vector = pos_vector + 1;
        end
        % Find last position of NaN in vector for this subject
        nan_value = find(isnan(subject_vector)); 
        if ~isempty(nan_value) || length(subject_vector) == 1 % There are NaN values or there is only baseline (size 1)
            continue; % Shouldn't happen, if it does go to next subject
        else

            pos_all_vector = size(All_C_vector,1) + 1;
            All_C_vector(pos_all_vector,:) = subject_vector;
            % Build subject legend (to be sure allnan subjects are not included)
            pos = size(legend_C,2) + 1;
            legend_C{pos} = subj_labels_C{p};
            legend_C{pos+1} = subj_labels_C{p};
            % Plot
            hold on
            plot(1:2,subject_vector,'k');
            plot(1:2,subject_vector,'-o','MarkerSize',5,'color',[0,0,0],'MarkerFaceColor',[0,0,0]);
        end
    end
    
    %%  FE to the right in single plot
    
    % To display Ns later
    count_FE_priority_1 = 0;
    count_FE_priority_2 = 0;
    count_FE_priority_3 = 0;
    All_FE_vector = [];
    legend_FE = {};
    for p = 1:length(subj_labels_FE)
        % Find row where this participant is in base variable
        row_pos = find(strcmp(Base(:,1),subj_labels_FE{p}));
        % Create vector
        subject_vector = [];pos_vector = 1;
        for tp = [1,priority_follow_up(1)] % Baseline and priority 1
            col_pos = find(strcmp(Base(1,:),[time_point_label{tp} '_' vars_to_plot{vtp}]));
            if tp ~= 1 % so if we are in the follow-up
               if isnan(Base{row_pos,col_pos}) % The priority follow up is empty
                   col_pos = find(strcmp(Base(1,:),[time_point_label{priority_follow_up(2)} '_' vars_to_plot{vtp}]));
                   if isnan(Base{row_pos,col_pos}) % If it's still empty, go to priority number 3
                       col_pos = find(strcmp(Base(1,:),[time_point_label{priority_follow_up(3)} '_' vars_to_plot{vtp}]));
                       if isnan(Base{row_pos,col_pos}) % If still empty, this subject has no follow up, so don't plot
                           continue;
                       else
                           count_FE_priority_3 = count_FE_priority_3 + 1;
                       end
                   else
                       count_FE_priority_2 = count_FE_priority_2 + 1;
                   end
               else
                   count_FE_priority_1 = count_FE_priority_1 + 1;
               end
            end
            subject_vector(pos_vector) = Base{row_pos,col_pos};
            pos_vector = pos_vector + 1;
        end
        % Find last position of NaN in vector for this subject
        nan_value = find(isnan(subject_vector)); 
        if ~isempty(nan_value) || length(subject_vector) == 1 % There are NaN values or there is only baseline (size 1)
            continue; % Shouldn't happen, if it does go to next subject
        else

            pos_all_vector = size(All_FE_vector,1) + 1;
            All_FE_vector(pos_all_vector,:) = subject_vector;
            % Build subject legend (to be sure allnan subjects are not included)
            pos = size(legend_FE,2) + 1;
            legend_FE{pos} = subj_labels_FE{p};
            legend_FE{pos+1} = subj_labels_FE{p};
            % Plot
            hold on;
            plot(3:4,subject_vector,'r');
            plot(3:4,subject_vector,'-o','MarkerSize',5,'color',[1,0,0],'MarkerFaceColor',[1,0,0]);
        end
    end
    legend_final = [legend_C';legend_FE'];
    % Get max and min
    Total_values = [All_C_vector; All_FE_vector];
    max_total = max(Total_values(:));
    min_total = min(Total_values(:));
        
    % Get mean values to plot
    try 
        mean_C_1 = mean(All_C_vector(:,1),1);
    catch
        mean_C_1 = NaN;
    end
    try
        mean_C_2 = mean(All_C_vector(:,2),1);
    catch
        mean_C_2 = NaN;
    end
    try
        mean_FE_1 = mean(All_FE_vector(:,1),1);
    catch
        mean_FE_1 = NaN;
    end
    try
        mean_FE_2 = mean(All_FE_vector(:,2),1);
    catch
        mean_FE_2 = NaN;
    end
    
    
    % Add a percentage of increase in y axis
    if min_total < 0
        min_y_axis = (abs(min_total)*y_axis_increase_percent) * -1;
    else
        min_y_axis = min_total/y_axis_increase_percent;
    end
    if max_total < 0
        max_y_axis = (abs(min_total)/y_axis_increase_percent) * -1;
    else
        max_y_axis = max_total*y_axis_increase_percent;
    end
        
    % Plot mean values
    plot(0.5:1.5,[mean_C_1,mean_C_1],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
    plot(1.5:2.5,[mean_C_2,mean_C_2],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
    plot(2.5:3.5,[mean_FE_1,mean_FE_1],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
    plot(3.5:4.5,[mean_FE_2,mean_FE_2],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
  
    xlim([0.5,4.5]);
    xticks([1.5,3.5]);
    count_1 = count_C_priority_1 + count_C_priority_2 + count_C_priority_3;
    count_2 = count_FE_priority_1 + count_FE_priority_2 + count_FE_priority_3;
    xticklabels({['C (' num2str(count_1) ')'],['FE (' num2str(count_2) ')']});
    legend(legend_final);
    pom=findobj('type','legend');
    delete(pom);
    set(gca,'TickLength',[0 0])
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',14)
    ylim([min_y_axis,max_y_axis]);

    
    % Identify how many increase (C)
    decreasing_C = 0;
    for ar = 1:size(All_C_vector,1)
        if All_C_vector(ar,1) > All_C_vector(ar,2)
            decreasing_C = decreasing_C + 1;
        end
    end
    
    % Identify how many increase (FE)
    decreasing_FE = 0;
    for ar = 1:size(All_FE_vector,1)
        if All_FE_vector(ar,1) > All_FE_vector(ar,2)
            decreasing_FE = decreasing_FE + 1;
        end
    end
    
    % Get worse variable depends on which var_to_plot we are dealing with
    is_in_decrease_good = find(strcmp(decrease_good,vars_to_plot{vtp}));
    if ~isempty(is_in_decrease_good)
        get_worse_C = num2str(count_1-decreasing_C); % That is, the worse ones are the ones that increase
        get_worse_FE = num2str(count_2-decreasing_FE);
    else
        get_worse_C = num2str(decreasing_C);
        get_worse_FE = num2str(decreasing_FE);
    end
    
    percentage_worse_C = round((str2double(get_worse_C)/count_1)*100); % round((17/38)*100)
    percentage_worse_FE = round((str2double(get_worse_FE)/count_2)*100);
    
    % Annotate how many show decrease
    annotation('textbox', [0.16, 0.85, 0.1, 0.1], 'String', [get_worse_C '/' num2str(count_1) ' worsen (' num2str(percentage_worse_C) '%)'],'LineStyle','none','FontSize',14);
    annotation('textbox', [0.61, 0.85, 0.1, 0.1], 'String', [get_worse_FE '/' num2str(count_2) ' (' num2str(percentage_worse_FE) '%)'],'LineStyle','none','FontSize',14);
    
    % Annotate how many subjects of each follow up
    annotation('textbox', [0.43, 0.2, 0.1, 0.1], 'String', [priority_labels{1} ': ' num2str(count_C_priority_1)],'LineStyle','none','FontSize',12);
    annotation('textbox', [0.43, 0.16, 0.1, 0.1], 'String', [priority_labels{2} ': ' num2str(count_C_priority_2)],'LineStyle','none','FontSize',12);
    annotation('textbox', [0.43, 0.12, 0.1, 0.1], 'String', [priority_labels{3} ': ' num2str(count_C_priority_3)],'LineStyle','none','FontSize',12);
    annotation('textbox', [0.82, 0.2, 0.1, 0.1], 'String', [priority_labels{1} ': ' num2str(count_FE_priority_1)],'LineStyle','none','FontSize',12);
    annotation('textbox', [0.82, 0.16, 0.1, 0.1], 'String', [priority_labels{2} ': ' num2str(count_FE_priority_2)],'LineStyle','none','FontSize',12);
    annotation('textbox', [0.82, 0.12, 0.1, 0.1], 'String', [priority_labels{3} ': ' num2str(count_FE_priority_3)],'LineStyle','none','FontSize',12);
    
    
    % Var in y axis (only left plot)
    ylabel_title = strrep(vars_to_plot{vtp},'_',' ');
    ylabel(ylabel_title,'fontsize',14,'FontWeight','bold');
    % Title
    title('Baseline vs follow-up','FontSize',14);
    
    if save_fig == 1
        saveas(gcf,[save_fig_path 'Bas_vs_1_follow_up_' vars_to_plot{vtp} '.emf']);
    end
    
end
    
end

%% Calculate Spearman's correlations (Baseline or Followup)
clear;clc; close all
set(0,'defaultfigurecolor',[1 1 1]); % I want white backgrounds
plot_only_significant_ones = 1; % 0 = NO; 1 = YES;
% Modify here for convenience if running only this step
name_Mega_variable = 'Mega_variable_June_22'; % Mega_variable_matched
bas_or_follow = 'Baseline'; % 'Baseline', '3mo, '6mo', '1yr'
if strcmp(bas_or_follow,'Baseline')
    time_label = 'A';
elseif strcmp(bas_or_follow,'3mo')
    time_label = 'B';
elseif strcmp(bas_or_follow,'6mo')
    time_label = 'C';
elseif strcmp(bas_or_follow,'1yr')
    time_label = 'D';
end
group_to_plot = {'FE'}; % 'FE', 'C', 'ALL'

% Silly way to compress a large section of cell arrays below
if 1 == 1
matching_measures = {'AGE','SEX','PSES','YRSED','VOCAB_TS'};
neuropsych_measures = {'OVERALLTSCR','MATRIX_TS','FULL2IQ','SPEEDTSCR','ATT_VIGTSCR',...
    'WMTSCR', 'VERBTSCR','SANITM','SAPITM'};
functioning_measures = {'ROLECURR','ROLELOW','ROLEHIGH','SOCIALCURR','SOCIALLOW',...
    'SOCIALHIGH','SFS_WITHDRAW_RS','SFS_INTERACT_RS','SFS_RECREAT_RS','SFS_OCCUP_RS',...
    'SFS_IND_PERF_RS','SFS_IND_COMP_RS','SFS_PROSOC_RS'};
clinical_measures = {'PANSSP_RS','PANSSN_RS','PANSST_RS','PSYRATS_AUD_HALL','PSYRATS_DELUSIONS'};
composite_measures = {'SANSSAPS_Level7_AudHall','SANSSAPS_Level7_UnusPercBeh',...
    'SANSSAPS_Level7_Delusions','SANSSAPS_Level7_ThDis','SANSSAPS_Level7_Inattention',...
    'SANSSAPS_Level7_Inexpress','SANSSAPS_Level7_Apathy','SANSSAPS_Level4_RealityDis',...
    'SANSSAPS_Level4_ThDis','SANSSAPS_Level4_Inexpress','SANSSAPS_Level4_Apathy',...
    'PANSS_Affect','PANSS_Disorg','PANSS_Negative','PANSS_Positive','PANSS_Resistance',...
    'BPRS_Total','BPRS_Positive','BPRS_Negative','BPRS_DeprAnx','BPRS_ActMania','BPRS_HostSusp'};
auditory_threshold_measures = {'QT_L_1000Hz', 'QT_R_1000Hz', 'QT_L_1500Hz', 'QT_R_1500Hz',...
    'QT_L_2000Hz','QT_R_2000Hz', 'QT_L_3000Hz', 'QT_R_3000Hz', 'QT_L_4000Hz', 'QT_R_4000Hz',...
    'AUD_THRESH_L','AUD_THRESH_R'};
duration_of_illness_measures = {'DAYS_SINCE_PRDM','DAYS_SINCE_1STEP','DAYS_SINCE_DISOR','DUP','PSYCH2SCAN','MED2SCAN'};
medication_measures = {'CPZ_equivalent'};
Struct_measures_vol = {'vol_A1_L','vol_A1_R','vol_LBelt_L','vol_LBelt_R',...
    'vol_PBelt_L','vol_PBelt_R','vol_OFC_L_merged','vol_OFC_R_merged',...
    'vol_IFG_L','vol_IFG_R','vol_AUDCORTEX_L','vol_AUDCORTEX_R'};
Struct_measures_thickness = {'thck_A1_L','thck_A1_R','thck_LBelt_L','thck_LBelt_R',...
    'thck_PBelt_L','thck_PBelt_R'};
end

% Define which variable/groups of variables to correlate
CORR_1 = [composite_measures]; % e.g. [matching_measures, neuropsych_measures]
CORR_2 = [Struct_measures_vol]; % e.g. [source_measures, PTE_measures];
% matching_measures, neuropsych_measures, functioning_measures,
% composite_measures, auditory_threshold_measures,
% duration_of_illness_measures, medication_measures
% Struct_measures_thickness, Struct_measures_vol

% OR set specific ones
% CORR_1 = {'norm_vol_A1_L'};
% CORR_2 = {'norm_vol_A1_L'};
% {'vol_A1_L';'thck_A1_L';'vol_A1_R';'thck_A1_R';'vol_AUDCORTEX_L';'vol_AUDCORTEX_R'};


pvals = []; p_vals_pos = 1; header_pvals = {};
for gtp = 1:length(group_to_plot)
load(['C:/Project/User/Cross_studies_data/Mega_variable/Longitudinal/' name_Mega_variable '.mat']);
% Replace 1's and 2's (used for SPSS) by typical labels for convenience
Mega_variable = Base;
for i = 2:size(Mega_variable,1) % All but header
    if Mega_variable{i,2} == 1 % C
        Mega_variable{i,2} = 'C';
    elseif Mega_variable{i,2} == 2 % FE
        Mega_variable{i,2} = 'FE';
    end
end
if strcmp(group_to_plot,'FE')
    % Find FE subjects
    warning('Selecting FE only...');
    group_indices = find(strcmp(Mega_variable(:,2),'FE'));
    group_indices = [1 group_indices']; % Add header column
    Mega_variable = Mega_variable([group_indices],:);
    color_group = [255 0 0];
elseif strcmp(group_to_plot,'C')
    % Find FE subjects
    warning('Selecting C only...');
    group_indices = find(strcmp(Mega_variable(:,2),'C'));
    group_indices = [1 group_indices']; % Add header column
    Mega_variable = Mega_variable([group_indices],:);
    color_group = [0 0 0];
elseif strcmp(group_to_plot,'ALL')
    color_group = [50 50 50];
    warning('Selecting FE and Controls...');
end

for bm = 1:length(CORR_2)
% Manually replace any empty values with NaN
pos_measu = find(strcmp(Mega_variable(1,:),[time_label '_' CORR_2{bm}]));
for i = 1:length(CORR_1)
    pos_col = find(strcmp(Mega_variable(1,:),[time_label '_' CORR_1{i}]));
    Cell_corr = Mega_variable(:,[pos_measu,pos_col]);
    % Correct for "[]" cells
    empty_cells = find(cellfun(@isempty,Cell_corr(:,2)));
    Cell_corr(empty_cells,2) = num2cell(NaN); % Just in case
    Table_corr = cell2table(Cell_corr(2:end,:));
    Table_corr.Properties.VariableNames = Cell_corr(1,:);
    try
        [R, PValue] = corrplot(Table_corr,'type','Spearman','testR', 'on');
        if plot_only_significant_ones == 1
            if PValue(2,1) >= 0.05
                fig = gcf;
                close(fig);
                continue;
            end
        end
        pvals(p_vals_pos) = PValue(2,1);
        table_columns = Table_corr.Properties.VariableNames;
        header_pvals{p_vals_pos} = [table_columns{1} '_&_' table_columns{2}];
        p_vals_pos = p_vals_pos + 1;
        % Check n included in this correlation
        Test_NaN = Table_corr{:,:}; % Convert table to Matrix
        [NaNrows, ~] = find(isnan(Test_NaN)); % Check rows with NaN
        n_corr = (size(Test_NaN,1)) - length(NaNrows); % determine n based on previous info
        % Delete all but the convenient subplot
        h = get(gcf, 'children');
        % get true x labels first (weird issue with corrplot)
        true_x_axis_labels = get(h(1), 'YTick'); 
        delete(h(1:5));
        % Adjust labels
        current_title = [time_label '_' CORR_2{bm} ' (n = ' num2str(n_corr) ' ' group_to_plot{gtp} ')'];
        current_title = strrep(current_title,'_',' ');
        xlabel(current_title,'Color','k');
        current_title = [time_label '_' CORR_1{i}];
        current_title = strrep(current_title,'_',' ');
        ylabel(current_title,'Color','k');
        % Adjust x axis labels (weird issue with corrplot)
        xticklabels(true_x_axis_labels)
        % Adjust color and style
        hline = findobj(gcf, 'type', 'line');
        set(hline(1),'Color','k') % Line fit
        set(hline(2),'Color',color_group/256) % Actual dots
        set(hline(2),'LineWidt',3) % Actual dots
    catch
        p_vals_pos = p_vals_pos + 1;
        continue;
    end
end
end
end
% Define pvals variable to back-track the exact pvalues if needed
PVALS_STRING = {};
PVALS_STRING(1,:) = header_pvals;
PVALS_STRING(2,:) = num2cell(pvals);
PVALS_STRING = PVALS_STRING';

%% Scatter plots with single data values (C vs FE)

clear;clc; close all
set(0,'defaultfigurecolor',[1 1 1]); % I want white backgrounds in plots
% path_list_cell = regexp(path,pathsep,'Split');
% if ~any(ismember('C:/Project/User/Project_MMN_baseline/Scripts',path_list_cell))
%     addpath('C:/Project/User/Project_MMN_baseline/Scripts'); % For UnivarScatter
% end
save_fig_path = 'C:/Project/User/Cross_studies_data/Figures/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFY THIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_Mega_variable = 'Mega_variable_June_22'; % Mega_variable_matched
save_fig = 0; % 1 = YES; 0 = NO
bas_or_follow = '6mo'; % 'Baseline', '3mo, '6mo', '1yr'
if strcmp(bas_or_follow,'Baseline')
    time_label = 'A';
elseif strcmp(bas_or_follow,'3mo')
    time_label = 'B';
elseif strcmp(bas_or_follow,'6mo')
    time_label = 'C';
elseif strcmp(bas_or_follow,'1yr')
    time_label = 'D';
end
plot_only_significant_ones = 0;
group_to_plot = {'FE','C'};
color_group = [[255 0 0]/256;[0 0 0]/256]; % Specific for scatter

% Silly way to compress a large section of cell arrays below
if 1 == 1
matching_measures = {'AGE','SEX','PSES','YRSED','VOCAB_TS'};
neuropsych_measures = {'OVERALLTSCR','MATRIX_TS','FULL2IQ','SPEEDTSCR','ATT_VIGTSCR',...
    'WMTSCR', 'VERBTSCR','SANITM','SAPITM'};
functioning_measures = {'ROLECURR','ROLELOW','ROLEHIGH','SOCIALCURR','SOCIALLOW',...
    'SOCIALHIGH','SFS_WITHDRAW_RS','SFS_INTERACT_RS','SFS_RECREAT_RS','SFS_OCCUP_RS',...
    'SFS_IND_PERF_RS','SFS_IND_COMP_RS','SFS_PROSOC_RS'};
clinical_measures = {'PANSSP_RS','PANSSN_RS','PANSST_RS','PSYRATS_AUD_HALL','PSYRATS_DELUSIONS'};
composite_measures = {'SANSSAPS_Level7_AudHall','SANSSAPS_Level7_UnusPercBeh',...
    'SANSSAPS_Level7_Delusions','SANSSAPS_Level7_ThDis','SANSSAPS_Level7_Inattention',...
    'SANSSAPS_Level7_Inexpress','SANSSAPS_Level7_Apathy','SANSSAPS_Level4_RealityDis',...
    'SANSSAPS_Level4_ThDis','SANSSAPS_Level4_Inexpress','SANSSAPS_Level4_Apathy',...
    'PANSS_Affect','PANSS_Disorg','PANSS_Negative','PANSS_Positive','PANSS_Resistance',...
    'BPRS_Total','BPRS_Positive','BPRS_Negative','BPRS_DeprAnx','BPRS_ActMania','BPRS_HostSusp'};
auditory_threshold_measures = {'QT_L_1000Hz', 'QT_R_1000Hz', 'QT_L_1500Hz', 'QT_R_1500Hz',...
    'QT_L_2000Hz','QT_R_2000Hz', 'QT_L_3000Hz', 'QT_R_3000Hz', 'QT_L_4000Hz', 'QT_R_4000Hz',...
    'AUD_THRESH_L','AUD_THRESH_R'};
duration_of_illness_measures = {'DAYS_SINCE_PRDM','DAYS_SINCE_1STEP','DAYS_SINCE_DISOR','DUP','PSYCH2SCAN','MED2SCAN'};
medication_measures = {'CPZ_equivalent'};
Struct_measures_vol = {'vol_A1_L','vol_A1_R','vol_LBelt_L','vol_LBelt_R',...
    'vol_PBelt_L','vol_PBelt_R','vol_OFC_L_merged','vol_OFC_R_merged',...
    'vol_IFG_L','vol_IFG_R','vol_AUDCORTEX_L','vol_AUDCORTEX_R'};
Struct_measures_thickness = {'thck_A1_L','thck_A1_R','thck_LBelt_L','thck_LBelt_R',...
    'thck_PBelt_L','thck_PBelt_R'};
end

% Which scales to compare 
var_scatter = [Struct_measures_vol];
% matching_measures, neuropsych_measures, functioning_measures,
% composite_measures, auditory_threshold_measures,
% duration_of_illness_measures, medication_measures
% Struct_measures_thickness, Struct_measures_vol

specific_var_scatter = {}; % Empty ({}) by default: e.g. 'F0_SNR_low_Constant' (will cancel previous)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(specific_var_scatter)
    scatter_vars = var_scatter;
else
    scatter_vars = specific_var_scatter;
end

load(['C:/Project/User/Cross_studies_data/Mega_variable/Longitudinal/' name_Mega_variable '.mat']);
Mega_variable = Base;

% Now prepare tables
% Since this Mega variable will have subject groups as numbers, 
% change them to labels before using it
for i = 2:size(Mega_variable,1) % Except first row
    if Mega_variable{i,2} == 1 % Controls
        Mega_variable{i,2} = 'C';
    elseif Mega_variable{i,2} == 2 % FE
        Mega_variable{i,2} = 'FE';
    end
end
for sv = 1:length(scatter_vars) % Brain measure
pos_measu = find(strcmp(Mega_variable(1,:),[time_label '_' scatter_vars{sv}]));
table_scatter = [];
for pg = 1:length(group_to_plot)
    group_indices = find(strcmp(Mega_variable(:,2),group_to_plot{pg}));
    for i = 1:length(group_indices)
        if isempty(Mega_variable{group_indices(i),pos_measu})
            table_scatter(i,pg) = NaN;
        else
            table_scatter(i,pg) = Mega_variable{group_indices(i),pos_measu};
        end
    end
end

% If different numbers of C and FE, it adds zeros to complete tables, correct for that
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(group_to_plot) ~= 2 % Just to be sure in case we include chronics
    error('If using more than two groups, reprogram next lines');
end
for pg = 1:length(group_to_plot)
    eval(['group_indices_' num2str(pg) ' = find(strcmp(Mega_variable(:,2),group_to_plot{pg}));'])
end
% If they are the same size it won't do anything
if length(group_indices_1) > length(group_indices_2)
    difference =  length(group_indices_1) - length(group_indices_2);
    % So this means first column is 'difference' rows longer than second 
    table_scatter(end+1-difference:end,2) = NaN;
    % Therefore those 'extra' positions at the end of columnn 2 are NaN
elseif length(group_indices_2) > length(group_indices_1)
    difference = length(group_indices_2) - length(group_indices_1);
    % So this means second column is 'difference' rows longer than first 
    table_scatter(end+1-difference:end,1) = NaN;
    % Therefore those 'extra' positions at the end of columnn 1 are NaN
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now, if one of the colums of table_scatter are all NaNs, continue to next
nan_col_1 = find(isnan(table_scatter(:,1)));
if length(nan_col_1) == length(table_scatter(:,1))
    continue; % To next variable
end
nan_col_2 = find(isnan(table_scatter(:,2)));
if length(nan_col_2) == length(table_scatter(:,2))
    continue; % To next variable
end

% Compute an independent sample t test and retrieve p value
[~,p_value,~,t_stats] = ttest2(table_scatter(:,1),table_scatter(:,2));
% Move onto next variable if not significant
if plot_only_significant_ones == 1
    if p_value >= 0.05
        continue;
    end
end

% Get mean values
mean_plot_left = nanmean(table_scatter(:,1));
mean_plot_right = nanmean(table_scatter(:,2));

point_settings = {'MarkerFaceColor',color_group,'MarkerEdgeColor','white','PointSize',80,'LineWidth',1};
plot_settings = [point_settings]; % Something weird about additional wiskers (in case needed)
figure;
[xPositions, yPositions, Label, RangeCut, FigHandles] = UnivarScatter(table_scatter,plot_settings{:});

% set(gcf,'Position',[0,0,600,300])
set(gcf,'Position',[500,250,300,300])
y_title = [time_label '_' scatter_vars{sv}];
y_title = strrep(y_title,'_',' ');
ylabel(y_title); 
xticklabels(group_to_plot) 
h=gca; h.XAxis.TickLength = [0 0];
h.YGrid = 'on';
h.GridLineStyle = '--';

% Add p value text
if p_value < 0.05
    color_p_value = 'red';
else
    color_p_value = 'black';
end
if p_value < 0.001
    label_pvalue = 'p < 0.001';
else
    str_pvalue = num2str(p_value);
    label_pvalue = ['p = ' str_pvalue(1:5)];
end
title({['\color{' color_p_value '}' label_pvalue '']})

% Add longer mean lines
hold on;
x_values = xlim;
plot([x_values(1)+x_values(1)*0.25 x_values(2)-x_values(2)*0.45],[mean_plot_left,mean_plot_left],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
hold on;
plot([x_values(2)-x_values(2)*0.35 x_values(2)-x_values(2)*0.05],[mean_plot_right,mean_plot_right],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);

if save_fig == 1
    saveas(gcf,[save_fig_path 'Scatter_' time_label '_' scatter_vars{sv} '.emf']);
end
end
