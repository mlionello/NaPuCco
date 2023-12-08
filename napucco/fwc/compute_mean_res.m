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

homepath = fullfile(path);
[fpr_cat, power_cat, params_fpr, params_power] = parse_result_folder(homepath);
parts = split(homepath, '/');
if ~isfolder(fullfile(homepath, "groupingavg"))
    mkdir(fullfile(homepath, "groupingavg"));
end
% load(fullfile(fullfile(parts{1:end-1}), 'opts.mat'), 'opts')
% load(fullfile(fullfile(parts{1:end-1}), 'hyparams.mat'), 'hyparams')

if size(power_cat, 1) > 1
    if contains(path, 'cluster')
        nb_subj = 4;
        dims = size(power_cat);
        [mean_power_ext, perm_rep] = create_power_shuffled_exp(power_cat, dims);
        plot_and_save_power1(mean_power_ext, nb_subj, homepath, params_power, dims, perm_rep);
        plot_and_save_power2(mean_power_ext, homepath, params_power, dims, nb_subj, perm_rep);
    end

    power = squeeze(mean(power_cat > 0));
    save(fullfile(homepath, "groupingavg", "power.mat"), 'power');
else
    disp("NO POWER OUTPUT")
    power = nan;
end

if size(fpr_cat, 1) > 1
    dims = size(fpr_cat);
    if contains(path, 'voxel')
        fpr_cat = squeeze(any(fpr_cat, 2));
    end
    fpr = squeeze(mean(fpr_cat > 0)); % fpr is 1 (experiment is failed) if K is larger than the threshold in the null distribution
    plot_and_save_CI_bootstrap(fpr, fpr_cat, params_fpr, path, homepath)
    save(fullfile(homepath, "groupingavg", "fpr.mat"), 'fpr');
else
    disp("NO FPR OUTPUT")
    fpr = nan;
end
end

%%
%% PARSING RESULTS INTO FPR AND POWER

function [fpr, power, params_fpr, params_power] = parse_result_folder(homepath)
fold = dir(fullfile(homepath, "grouping"));
fpr = []; power = []; params_fpr = nan; params_power = nan;
m = 1; o = 1;
for i = 1: length(fold)
    if ~contains(fold(i).name, ".mat"); continue; end
    file2load = fullfile(fold(i).folder, fold(i).name);
    result = load(file2load, 'results');
    result = result.results;
    if isfield(result, "fpr")
        fpr(m, :, :, :, :) = result.fpr(:,:,:,:); %(p, s, v)
        m = m + 1;
        params_power = result;
    end
    if isfield(result, "power")
        power(o, :, :, :, :) = result.power; %(p, s, v, v)
        o = o + 1;
        params_fpr = result;
    end
end
end

%%
%% PLOT GENERATORS

function plot_and_save_CI_bootstrap(fpr, fpr_cat, params_fpr, path, homepath)
alpha_thr = 1;
vx_indices = params_fpr.non_eff_range;
if contains(path, 'cluster')
    vx_indices = vx_indices(end - min(5, length(vx_indices)-1): end);
    vx_range = vx_indices;
    ttl = "size: %d^3 voxels";
    corrmethod = "Clusterbased";
elseif contains(path, 'voxel')
    vx_indices = 1: size(fpr, 2);
    vx_range = params_fpr.non_eff_range;
    ttl = "size: %d voxels";
    corrmethod = "Voxelwise";
end
n_bootstrap = 1000;
n_exp = size(fpr_cat, 1);
for n = 1: n_bootstrap
    perm = randi(n_exp, 1, n_exp);
    btstrp(n, :, :, :) = squeeze(mean(fpr_cat(perm, :, :, :)>0));
end
ci = quantile(btstrp, [0.05, 0.95]);
fig = figure('Visible', 'off', 'Position', [0, 0, 1000, 600]); hold on;
for j = 1:numel(vx_indices)
    subplot(ceil(numel(vx_indices)/3), 3, j)
    title(compose(ttl, vx_range(j)))
    for i = 1: size(ci, 2)
        rectangle('Position',[i-0.15, ...
            ci(1, i, vx_indices(j), alpha_thr), ...
            0.3, ...
            ci(2, i, vx_indices(j), alpha_thr) - ci(1, i, vx_indices(j), alpha_thr)]);
        line([i-0.15, i+0.15], ...
            [fpr(i, vx_indices(j), alpha_thr), fpr(i, vx_indices(j), alpha_thr)], ...
            'LineWidth', 2, ...
            'Color', 'black');
    end
    ylabel('fpr','FontSize',18, 'Interpreter','latex')
    yticks([0.0:0.025:0.075])
    xlabel('number of combined subjects','FontSize',18, 'Interpreter','latex')
    xticks(1: size(ci, 2))
    xticklabels(params_fpr.subj_range)
    yline(0.05, 'LineStyle', '--')
    ylim([0, 0.09])
end
sgtitle(compose("%s correction. Confidence Intervals by %d bootstraps (5%%-95%%) from %d simulations", corrmethod, n_bootstrap, size(fpr_cat,1)),'FontSize',18,  'interpreter', 'latex')
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
padding_m = 8: 14;
padding_m = unique(min(padding_m, dims(3)));
pow_thr = 0.03;
for pad_i = padding_m
    mean_power = squeeze(mean_power_ext(:, :, nb_subj, pad_i, 1));
    below_thr = abs(mean_power-mean_power(end,:))<pow_thr;
    from0to1 = diff(below_thr) == 1;
    for k=1: size(from0to1, 2)
        tmp = find(from0to1(:,k) == 1, 1, 'last') + 1;
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
    title(compose("Number of experiments needed for convergence (thr %.2f) %d subj, %d non-eff padding", pow_thr, params_power.subj_range(nb_subj), pad_i), 'Interpreter','latex')
    figfilename = compose("nbexp2conv_%dsubj_%dpadding_%drep_%dpop", params_power.subj_range(nb_subj), pad_i, perm_rep, dims(1));
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
if stop_iterations(end)>dims(1)
    stop_iterations = int32(1:dims(1)/4:dims(1)*3/4); end
stop_iterations = [stop_iterations, dims(1)];
for alpha_thr = 1:3
    fig = figure('Visible', 'off', 'Position', [0 0 1200 600]);
    maxy = max(mean_power_ext(stop_iterations,:,nb_subj,:,alpha_thr), [], 'all')+0.01;
    miny = min(mean_power_ext(stop_iterations,:,nb_subj,:,alpha_thr), [], 'all')-0.01;
    for subplot_i = 1:length(stop_iterations)
        subplot(2, 2, subplot_i);
        plot(squeeze(mean_power_ext(stop_iterations(subplot_i),:,nb_subj,:,alpha_thr))', ...
            'Color', [0.4 0.4 0.4 0.4]);
        title(compose('stop at %dth run', stop_iterations(subplot_i)));
        xlabel('padding (half width)');
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
mean_power_ext = nan([dims(1), perm_rep, dims(2), dims(3), dims(4)]);
cr="";
for j = 1:perm_rep
    msg2print = compose("perm nb: %03d/%d", j, perm_rep);
    fprintf(cr + msg2print, j, perm_rep); pause(eps);
    cr = repmat('\b', 1, strlength(msg2print));
    perm_indices = randperm(dims(1));
    shuffled_data = power_cat(perm_indices, :, :, :);
    for i = 1:dims(1)
        mean_power_ext(i, j, :, :, :) =  squeeze(mean(shuffled_data(1: i, :, :, :) > 0, 1));
    end
end
fprintf("\n")
end