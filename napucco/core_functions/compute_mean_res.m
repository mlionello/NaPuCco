function [power, fpr] = compute_mean_res(path)
% compute_mean_res:
%   This function computes the mean power and false positive rate (FPR)
%   from a given data file specified by the path. The argument must match
%   the path to a folder containing 'grouping/' results.
%   According to the experiment simulated, the function will generate some
%   results plot based on bootstrap confidence interval, in the case of
%   fpr, and statistics convergence while varying padding size, in the case
%   of clusterbased power.
%
% Syntax:
%   [power, fpr] = compute_mean_res(path)
%
% Input Arguments:
%   1. path - String specifying the file path to the grouped results.
%
% Output Arguments:
%   1. power - Mean power calculated from the simulated experiments.
%   2. fpr - False Positive Rate (FPR) computed from the simulated experiments.
%
% Example:
%   filePath = 'path/to/your/data/folder';
%   [meanPower, meanFPR] = compute_mean_res(filePath);
%

% Check clusterbased false-positive inflation:
%
% max_k = squeeze(max(results.K_fpr(:, 1, :), [], 'all'));
% num_vx_ranges = size(results.K_fpr, 3);
% a = nan(num_vx_ranges, max_k);
% subj_range = 1;
% for i = 1:num_vx_ranges; max_i = max(results.K_fpr(:,subj_range,i)); for j =1 :max_i; a(i,j)=sum(results.K_fpr(:,subj_range,i)==j); end; end
% vx_range_ind = 10;
% test_dist = [];
% num_K_available = sum(~isnan(a(vx_range_ind, :)));
% for j = 1 : num_K_available; test_dist = [test_dist ones(1, a(vx_range_ind, j))*j]; end
% for j = 1 : num_K_available; p_val = tiedrank(-[j test_dist])/(length(test_dist)+1); disp(p_val(1)); end
% tiedrank(-[3, ones(1, 288), ones(1,671)*2, ones(1,41)*3, ones(1,1)*4])/1002; ans(1)
% tiedrank(-[2, ones(1, 288), ones(1,671)*2, ones(1,41)*3, ones(1,1)*4])/1002; ans(1)
%

homepath = fullfile(path);
[fpr_cell, power_cell, params_fpr, params_power] = parse_result_folder(homepath);
power_cell = cellfun(@(x) reshape(x, [1 size(x)]), power_cell, 'UniformOutput', false);
fpr_cell = cellfun(@(x) reshape(x, [1 size(x)]), fpr_cell, 'UniformOutput', false);
power_cat = cell2mat(power_cell);
fpr_cat = cell2mat(fpr_cell);

%fpr_cat = reshape(fpr_cat, size(fpr_cat, 1), size(fpr_cat, 2), 8, 5, size(fpr_cat, 4), size(fpr_cat, 5));
fpr_cat = reshape(fpr_cat, size(fpr_cat, 1), 8, 5, size(fpr_cat, 3), size(fpr_cat, 4));

%power_cat = reshape(power_cat, size(power_cat, 1), 8, 5, size(power_cat, 3), size(power_cat, 4));
power_cat = reshape(power_cat, size(power_cat, 1), 8, 5, size(power_cat, 3), size(power_cat, 4));

parts = split(homepath, '/');
out_folder = fullfile('/', parts{1:end-2});
out_folder = fullfile(out_folder, "groupingavg");
if ~isfolder(out_folder)
    mkdir(out_folder);
end
% load(fullfile(fullfile(parts{1:end-1}), 'opts.mat'), 'opts')
% load(fullfile(fullfile(parts{1:end-1}), 'hyparams.mat'), 'hyparams')

if any(~isnan(power_cat), 'all')
    if contains(path, 'cluster') && size(power_cat, 4)>1
        nb_subj = size(power_cat, 2);
        nb_noneffsubj = size(power_cat, 3);
        noneff_subj = 5;
        nb_subj = 3;
        pit_thr=1;

        dims = size(power_cat);
        [mean_power_ext, perm_rep] = create_power_shuffled_exp(power_cat, dims);
        plot_and_save_power1(squeeze(mean_power_ext(:, :, nb_subj, noneff_subj, pit_thr)), ...
            nb_subj, out_folder, params_power, dims, perm_rep);
        plot_and_save_power2(squeeze(mean_power_ext(:, :, nb_subj, noneff_subj, :)), ...
            out_folder, params_power, dims, nb_subj, perm_rep);
    end

    power = squeeze(mean(power_cat > 0));
    save(fullfile(out_folder, "power.mat"), 'power');
else
    disp("NO POWER OUTPUT")
    power = nan;
end

if any(~isnan(fpr_cat), 'all')
    dims = size(fpr_cat);
    if contains(path, 'voxel')
        fpr_cat = squeeze(any(fpr_cat, 2));
    end
    fpr = squeeze(mean(fpr_cat > 0)); % fpr is 1 (experiment is failed) if K is larger than the threshold in the null distribution
    
    noneff_subj_ind = 1;
    pi_thr_ind = 1;
    if contains(path, 'voxel')
        plot_and_save_CI_bootstrap(squeeze(fpr(:, noneff_subj_ind, :, pi_thr_ind)), ...
            squeeze(fpr_cat(:, :, noneff_subj_ind, :, pi_thr_ind)), ...
            params_fpr, path, out_folder, noneff_subj_ind, pi_thr_ind)
    else
        plot_and_save_CI_bootstrap(squeeze(fpr(:, :, pi_thr_ind)), ...
            squeeze(fpr_cat(:, :, :, pi_thr_ind)), ...
            params_fpr, path, out_folder, noneff_subj_ind, pi_thr_ind)
        %         plot_and_save_CI_bootstrap(squeeze(fpr(:, :, pi_thr_ind)), ...
        %             squeeze(fpr_cat(:, :, :, pi_thr_ind)), ...
        %             params_fpr, path, out_folder, noneff_subj_ind, pi_thr_ind)
    end
    save(fullfile(out_folder, "fpr.mat"), 'fpr');
else
    disp("NO FPR OUTPUT")
    fpr = nan;
end
end

%%
%% PARSING RESULTS INTO FPR AND POWER

function [fpr, power, params_fpr, params_power] = parse_result_folder(homepath)
fold = dir(fullfile(homepath));
fpr = {}; power = {}; params_fpr = nan; params_power = nan;
m = 1; o = 1;
cr = "";
for i = 1: length(fold)
    msg = compose("parsing file %d out of %d", i, length(fold));
    fprintf(cr + msg);      
    cr = repmat('\b',1 , strlength(msg));
    if ~contains(fold(i).name, ".mat"); continue; end
    file2load = fullfile(fold(i).folder, fold(i).name);
    load(file2load, 'results');
    if isfield(results, "fpr")
        fpr{m, 1} = results.fpr; %(p, s, v)
        m = m + 1;
        params_fpr = results;
    end
    if isfield(results, "power")
        power{o, 1} = results.power; %(p, s, v, v)
        o = o + 1;
        params_power = results;
    end
end
fprintf(cr)
end

%%
%% PLOT GENERATORS

function plot_and_save_CI_bootstrap(fpr, fpr_cat, params_fpr, path, homepath, alpha_thr, noneff_subj)
vx_range = params_fpr.non_eff_range;
subj_range = params_fpr.subj_range;
if contains(path, 'cluster')
    ttl = "size: %d^3 voxels";
    corrmethod = "Clusterbased";
elseif contains(path, 'voxel')
    ttl = "size: %d voxels";
    corrmethod = "Voxelwise";
end
n_bootstrap = 10000;
n_exp = size(fpr_cat, 1);
for n = 1: n_bootstrap
    perm = randi(n_exp, 1, n_exp);
    btstrp(n, :, :) = squeeze(mean(fpr_cat(perm, :, :) > 0));
end
ci = quantile(btstrp, [0.05, 0.95]);
fig = figure('Visible', 'on', 'Position', [0, 0, 1000, 600]); hold on;
for j = 1 : numel(vx_range)
    subplot(ceil(numel(vx_range)/3), 3, j)
    title(compose(ttl, vx_range(j)))
    for i = 1: numel(subj_range)
        rectangle('Position',[i - 0.15, ...
            ci(1, i, j), ...
            0.3, ...
            ci(2, i, j) - ci(1, i, j)]);
        line([i - 0.15, i + 0.15], ...
            [fpr(i, j), fpr(i, j)], ...
            'LineWidth', 2, ...
            'Color', 'black');
    end
    ylabel('fpr', 'FontSize', 18, 'Interpreter', 'latex')
    yticks([0.0 : 0.025 : 0.075])
    xlabel('number of combined subjects', 'FontSize', 18, 'Interpreter', 'latex')
    xticks(1: size(ci, 2))
    xticklabels(subj_range)
    yline(0.05, 'LineStyle', '--')
    ylim([0, 0.09])
end
sgtitle(compose("%s correction. Confidence Intervals by %d bootstraps (5-95) from %d simulations", corrmethod, n_bootstrap, size(fpr_cat,1)),'FontSize',18,  'interpreter', 'latex')
figfilename_init = sprintf("fpr_ci_btstrp_%d_%dexp_alphalvl%d", n_bootstrap, n_exp, alpha_thr);
figfilename = figfilename_init;
fid = 1;
while exist(fullfile(homepath, figfilename + ".png"), 'file')
    figfilename =  sprintf("%s_%02d", figfilename_init, fid);
    fid = fid + 1;
end
saveas(fig, fullfile(fullfile(homepath, figfilename + ".png")))
end

function plot_and_save_power1(mean_power_ext, nb_subj, homepath, params_power, dims, perm_rep)
padding_m = params_power.padding_range;
pow_thr = 0.03;
for pad_i = 1:numel(padding_m)
    padding_ind = padding_m(pad_i);
    mean_power = mean_power_ext(:, :, padding_ind);
    below_thr = abs(mean_power-mean_power(end, :)) < pow_thr;
    from0to1 = diff(below_thr) == 1;
    for k = 1 : size(from0to1, 2)
        tmp = find(from0to1(:, k) == 1, 1, 'last') + 1;
        lastSwitchPosition(k) = 0;
        if ~isempty(tmp)
            lastSwitchPosition(k) = tmp;
        end
    end
    fig = figure('Visible', 'off', 'Position', [0 0 900 600]);
    yline_lgn = compose("power mean %.2f", mean_power(end,1)); hold on;
    hist_lgn = sprintf("mu %.2f, std %.2f", mean(lastSwitchPosition), std(lastSwitchPosition));
    yline(mean_power(end,1), 'k--', 'DisplayName',yline_lgn);
    histogram(lastSwitchPosition, 10, 'Normalization', 'probability', 'DisplayName', hist_lgn)
    legend('show')
    for nline = 1:min(size(mean_power,2), 100)
        plot(1:size(mean_power,1), smooth(mean_power(:,nline)), 'Color', [0.4 0.4 0.4 0.3],'HandleVisibility','off');
        line([lastSwitchPosition(nline), lastSwitchPosition(nline)], [mean_power(end,1),0],'HandleVisibility','off', 'Color', [0.4 0.4 0.4 0.3],'LineStyle','--');
        scatter(lastSwitchPosition(nline),0, 'xk', 'HandleVisibility','off')
    end
    ylim([0,1])
    ylabel('Power', 'Interpreter','latex')
    xlabel('Number of experiments', 'Interpreter','latex')
    title(compose("Number of experiments needed for convergence (thr %.2f) %d subj, %d non-eff padding", pow_thr, params_power.subj_range(nb_subj), padding_ind), 'Interpreter','latex')
    figfilename = compose("nbexp2conv_%dsubj_%dpadding_%drep_%dpop", params_power.subj_range(nb_subj), padding_ind, perm_rep, dims(1));
    fid = 1;
    while exist(fullfile(homepath, figfilename + ".png"), 'file')
        figfilename =  sprintf("%s_%02d", figfilename, fid);
        fid = fid +1;
    end
    saveas(fig, fullfile(homepath, figfilename + ".png"))
end
end

function plot_and_save_power2(mean_power_ext, homepath, params_power, dims, nb_subj, perm_rep)
stop_iterations = [100, 200, 300];
padding_m = params_power.padding_range;
if stop_iterations(end)>dims(1)
    stop_iterations = int32(1:dims(1)/4:dims(1)*3/4); end
stop_iterations = [stop_iterations, dims(1)];
for alpha_thr = 1:3
    fig = figure('Visible', 'off', 'Position', [0 0 1200 600]);
    maxy = max(mean_power_ext(stop_iterations,:,padding_m,alpha_thr), [], 'all')+0.01;
    miny = min(mean_power_ext(stop_iterations,:,padding_m,alpha_thr), [], 'all')-0.01;
    for subplot_i = 1:length(stop_iterations)
        subplot(2, 2, subplot_i);
        plot(squeeze(mean_power_ext(stop_iterations(subplot_i),:,padding_m,alpha_thr))', ...
            'Color', [0.4 0.4 0.4 0.8]);
        title(compose('stop at %dth run', stop_iterations(subplot_i)));
        xlabel('padding (half width)');
        xticklabels(string(padding_m));
        ylabel('power');
        ylim([miny, maxy]);
    end
    sgtitle(compose('%d subjs, varying non eff padding and stopping at ith experiments ( %d exps) alphalvl %d',params_power.subj_range(nb_subj), dims(1),alpha_thr))
    figfilename = sprintf("power_varying_padding_different_stops_%dsubj_%dexp_alphalvl%d",params_power.subj_range(nb_subj), dims(1), alpha_thr);
    fid = 1;
    while exist(fullfile(homepath, figfilename + ".png"), 'file')
        figfilename =  sprintf("%s_%02d", figfilename, fid);
        fid = fid + 1;
    end
    saveas(fig, fullfile(fullfile(homepath, figfilename + ".png")))
end
end

function [mean_power_ext, perm_rep] = create_power_shuffled_exp(power_cat, dims)
perm_rep = min(dims(1), 500);
mean_power_ext = cell(dims(1), perm_rep);
cr="";
for j = 1: perm_rep
    msg2print = compose("perm nb: %03d/%d", j, perm_rep);
    fprintf(cr + msg2print, j, perm_rep); pause(eps);
    cr = repmat('\b', 1, strlength(msg2print));

    perm_indices = randperm(dims(1));
    shuffled_data = power_cat(perm_indices, :, :, :, :);
    for i = 1 : dims(1)
        shuffled_mean = squeeze(mean(shuffled_data(1: i, :, :, :, :), 1));
        mean_power_ext{i, j} =  reshape(shuffled_mean, [1, 1, size(shuffled_mean)]);
    end
end
mean_power_ext = cell2mat(mean_power_ext);
fprintf("\n")
end