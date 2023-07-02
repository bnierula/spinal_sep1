% Author: Birgit Nierula
% nierula@cbs.mpg.de

function figure_singleTrial_cca(condition, srmr_nr, isubject, ...
    c_axis, iscolorbar, trial_number, is_norm)


% set path
savepath_ga = getenv('GADIR');

% condition info
cond_info = get_conditionInfo(condition, srmr_nr);
cond_name = cond_info.cond_name;
nerve = cond_info.nerve;
str_stimulation = cond_info.str_stimulation;

if nerve == 2
    knee_chans = {'KneeM' 'Knee1' 'Knee2' 'Knee3' 'Knee4'};
    fname = [savepath_ga 'amplitudeAndLatency_allSubjects.mat'];
    load(fname, [cond_name '_values'])
end


% get data
dat_level = [3 5]; % 3=esg-TH6 5=esg-cca
dat_str = {'esg-th6Ref' 'esg-cca'};
for idat = 1:length(dat_level)
    idatlevel = dat_level(idat);
    is_raw = false;
    out = ga_combineData(isubject, condition, srmr_nr, ...
        is_raw, idatlevel);
    
    if nerve == 1
        chan_names = 'SC6';
        title_str = ['cervical ' dat_str{idat}];
        x_lim = [-20 60];
    elseif nerve == 2
        chan_names = 'L1';
        title_str = ['lumbar ' dat_str{idat}];
        x_lim = [-20 80];
    end
    
    epo = out{1};
    if idatlevel == 3
        epo.is_au = false;
        chan_idx = find(ismember({epo.chanlocs.labels}, chan_names));
    elseif idatlevel == 5
        epo.is_au = true;
        chan_idx = 1;
        chan_names = NaN;
    end

    
    %% squeeze data in time
    time_idx = find(epo.times >= x_lim(1) & epo.times <=x_lim(2));
    epo = pop_select(epo, 'point', time_idx);
    
    %% normalize
    tmp_values = squeeze(epo.data(chan_idx,:,:))';
    if is_norm
        [norm_values] = normalize_2zscore(tmp_values);
    else
        norm_values = tmp_values;
    end

    %% select trials
    idx_trials = 1:trial_number; 
    trialvalues = norm_values(idx_trials, :);

    %% make plots
    %% ------------
    fset = myFigureSettings(); % input: size(1) = width, size(2) = hight
    font_size = fset.font_size;
    fig_size(1) = fset.fig_size(1) - 1; % width
    fig_size(2) = fset.fig_size(2) * 2; % hight

    % ==================== FIGURE SETTINGS ========================
    figure;
    set(gcf, 'units', 'centimeters', 'position', [1 1 fig_size(1) fig_size(2)])


    % ==================== COLOR SETTINGS ========================
    color_code = {fset.mixed};
    face_alpha = 0.3;


    % ==================== GRAPH SETTINGS ========================
    line_width = 1;
    font_name = 'Roboto';
    if nerve == 1
        x_tickpoints = -20:10:60;
        x_ticklabels = {'-20' '' '0' '' '20' '' '40' '' '60'}';
    elseif nerve == 2
        x_tickpoints = -20:10:80;
        x_ticklabels = {'-20' '' '0' '' '20' '' '40' '' '60' '' '80'}';
    end

    set(gca,'linewidth',1)
    set(gca,'FontSize', font_size)
    set(gca,'FontName', font_name)

    % ==================== PLOT FIGURE ============================

    % plot single trials
    imagesc(trialvalues), caxis(c_axis)
    set(gca,'YDir','normal') 
    if iscolorbar
        colorbar();
    end

    xticks(find(ismember(epo.times, x_tickpoints)))
    xticklabels(x_ticklabels)
    xlabel('Time [ms]')
    ylabel ('Trial number')
    title(title_str)

    %% save
    fname = ['esg_' cond_name sprintf(['_singelTrial_sub-%03i_' dat_str{idat}], isubject)];
    print([getenv('FIGUREPATH') fname], '-dpng', '-painters') 
    print([getenv('FIGUREPATH') fname], '-dsvg', '-painters') 
end
