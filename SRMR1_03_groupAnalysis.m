% Author: Birgit Nierula
% nierula@cbs.mpg.de

%% SRMR1_03_groupAnalysis
% group statistics and figures

clear; clc
delete(gcp('nocreate')) % clear parallel pool


%% variables that need to be changed
% step
step_number = 8; disp(['step_number = ' step_number])
% sampling rate
sampling_rate = 1000;
% subjects
subjects = 1:36;


%% define variables and paths
% experiment
srmr_nr = 1;
% conditions
conditions = 2:3;

% set paths
datadir = '/data/p_02068/SRMR1_experiment/analyzed_data/';
anadir = '/data/pt_02068/analysis/final/';
bidsdir = '/data/p_02068/SRMR1_experiment/bids/';
setenv('CFGDIR', '/data/pt_02068/analysis/manuscript_sep/scripts/cfg_srmr1/')

setenv('RAWDIR', bidsdir) % here is the raw data
setenv('RPKDIR', [datadir 'Rpeak_detected/']) % here R-peak detected data (holds only ECG channel and trigger info)
setenv('ANADIR', [anadir 'tmp_data/']) % analysis directory
setenv('ESGDIR', [datadir 'esg/']);
setenv('EEGDIR', [datadir 'prepro_eeg_icaclean/'])
setenv('BSDIR', [datadir 'bs/']);
setenv('OTHERDIR', [datadir 'other/']);
setenv('GADIR', [datadir 'ga/']);
setenv('ZIMDIR', '/data/pt_02068/doc/LabBook_SRMR1/SRMR1/EXPERIMENT/analysis_ga/');
setenv('FIGUREPATH', [datadir 'figures/']); if ~exist(getenv('FIGUREPATH'), 'dir'), mkdir(getenv('FIGUREPATH')); end

% settings for figures
set(0, 'DefaulttextInterpreter', 'none')


% add toolboxes and other sources for scripts
addpath('/data/pt_02068/toolboxes/eeglab14_1_2b/') % eeglab toolbox
eeglab  % start eeglab and close gui
close


% scripts 
functions_path = '/data/pt_02068/analysis/manuscript_sep/scripts/functions/';
addpath(genpath(functions_path)) % scripts


switch step_number
    
    case 1
        %% save latency of all potentials (EEG, ESG, ENG,...)
        for isubject = 1:length(subjects)
            subject = subjects(isubject);
            for nerve = 1:2
                close all
                prepro_extractSepParameters(subject, nerve, srmr_nr);
            end
        end
        
    case 2
        %% detection of SEPs at different levels
        for icondition = 1:length(conditions)
            condition = conditions(icondition);
            ga_amplitudeAndLatency(subjects, condition, srmr_nr) % with NaNs for subjects where no potential was visible in averaged over all trials
            ga_amplitudeAndLatency_allSubjects(subjects, condition, srmr_nr) %  subjets where no potential was visible were substituted with amplitude at grand average latency (over all subj)
        end
    case 3     
        for icondition = conditions
            stat_number = 1 % detection 
            has_allsubj = true; % subjets where no potential was visible were substituted with the amplitude at grand average latency (over all subj)
            ga_detection_stats(subjects, icondition, srmr_nr, stat_number, has_allsubj);
            stat_number = 2 % mean and std
            ga_detection_stats(subjects, icondition, srmr_nr, stat_number, has_allsubj);
        end
        % make table
        get_detection_statsTable(conditions, srmr_nr)
    
    case 4    
        %% spinal isopotentialplot
        %% ESG
        dat_level = 3; % TH6-referenced esg data
        is_raw = false;
        for condition = conditions
            out.(['c' num2str(condition)]) = ga_combineData(subjects, condition, srmr_nr, ...
                is_raw, dat_level);
        end
        is_au = false; 
        if dat_level == 2 || dat_level == 5
            is_au = true;
        end
        % all subjects
        for condition = conditions
            dat = out.(['c' num2str(condition)]);
            figure_spinal_isopotentialplot(dat, subjects, condition, srmr_nr, is_au)
        end
    case 5    
        %% EEG
        dat_level = 1; % EEG data
        is_raw = false;
        for condition = conditions
            out.(['c' num2str(condition)]) = ga_combineData(subjects, condition, srmr_nr, ...
                is_raw, dat_level);
        end
        is_au = false; 
        if dat_level == 2 || dat_level == 5
            is_au = true;
        end
        % all subjects
        for condition = conditions
            title_str = '';
            dat = out.(['c' num2str(condition)]);
            figure_cortical_isopotentialplot(dat, subjects, condition, srmr_nr, is_au)
        end
    case 6
        %% early SEPs at different levels
        for condition = conditions            
            figure_stackedplot_seps(condition, srmr_nr, subjects)
        end
        
    case 7    
        %% spinal SEPs comparing anterior referenced data at anatomical point with CCA data
        for condition = conditions
            isubject = subjects;
            figure_spinalSEP(condition, isubject)
        end
        
    case 8
        %% single trial plots
        for condition = conditions
            isubject = 11; 
            if condition == 2
                c_axis = [-0.5 0.5];
            elseif condition == 3
                c_axis = [-0.5 0.5];
            end
            iscolorbar = true;
            trial_number = 1000;
            is_norm = true;
            figure_singleTrial_cca(condition, srmr_nr, isubject, ...
                c_axis, iscolorbar, trial_number, is_norm)
        end
    
    case 9
        %% preparing varables for robustness analysis
        % data structure: dat = subject x trial x channel
        ga_data4robustness_allsubj(subjects, conditions, srmr_nr)          
            
    case 10
        %% robustness analysis 
        robustnessSEP_wrapper
        
    case 11
        %% late potentials
        % same prepro but filter from 5 to 400 Hz
        % then cluster based permutation test + plots 
        for condition = conditions
            for isubject = subjects
                is_restcontrol = false; new_ref = [];
                epo_esg{isubject} = srmr_prepro_latepotentials_esg(isubject, condition, is_restcontrol, new_ref, srmr_nr);
                is_restcontrol = true;
                epo_control{isubject} = srmr_prepro_latepotentials_esg(isubject, condition, is_restcontrol, new_ref, srmr_nr);
            end
            trials = [cellfun(@(x) x.trials, epo_esg)' cellfun(@(x) x.trials, epo_control)'];
            % figure_late_seps(epo_esg, epo_control, condition)
            stats_late_seps(epo_esg, epo_control, condition, srmr_nr) % plots also the figure!
        end
        clear epo*
        
    case 12
        %% loop 5: TF plots
        % TF plots over all subjects
        is_evoked = true;
        for condition = conditions
            figure_spinalTFplots(subjects, condition, is_evoked)
        end
                
end
