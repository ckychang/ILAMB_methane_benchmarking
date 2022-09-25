clear all
close all
clc;
degree = sprintf('%c', char(176));
mon = [31,28,31,30,31,30,31,31,30,31,30,31];
mon_leap = [31,29,31,30,31,30,31,31,30,31,30,31];
mon_sum = cumsum(mon);
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
    "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
    "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];
mcolor = [.5 1.0 .83; .90 .17 .31; .0 .0039 1.00; 1.00 .57 .69; .2 .5 0;...
    .00 .75 1.00; .76 .60 .42; .45 .63 .76; .53 .33 .04; 0 0 0];
month = [];
month_leap = [];
for m = 1:12
    month = [month; m*ones(mon(m),1)];
    month_leap = [month_leap; m*ones(mon_leap(m),1)];
end
%% read in ilamb score from obs runs
clear model_ranking ilamb_score_obs threshold model_group
clf;
cdf_threshold = [.7 .8 .9];
var_name = ["Bias Score","RMSE Score","Seasonal Cycle Score","Overall Score"];
path = '/Volumes/LaCie/research/ILAMB_CH4/FCH4_benchmark/Sep2022/results/';
for region = 1:3
    if (region==1)
        loc = 'site_level_tier1';
    elseif (region==2)
        loc = 'site_level_tier1_north_of_30N';
    elseif (region==3)
        loc = 'site_level_tier1_south_of_30N';
    end
    % read in bias, rmse, sesonality, and overall scores
    file=(['scalar_database.csv']);
    T = readtable([path loc '/_build/' file]);
    name = unique(T{:,4});
    % do not include ELM-v0 and ML for the comparison
    idx = find(contains(name, "ELM_v0")==1);
    name(idx) = [];
    if (region==1)
        model_list = unique(T{:,4});
        model_list(idx) = [];
    end
    for jj = 1:length(model_list) % models
        idx = find(strcmpi(T{:,4},strtrim(model_list(jj))));
        if (length(idx)>0)
            tmp = T(idx,:);
            for kk = 1:4 % variables
                idx = find(strcmpi(tmp{:,5},var_name(kk)));
                ilamb_score_obs(jj,kk, region) = str2double(tmp{idx,10});
%                 ilamb_score_obs(jj,kk, region) = (tmp{idx,10});
            end
        else
            ilamb_score_obs(jj,:, region) = nan;
        end
    end
    if (region==1)
        model_name_full = name;
    end
end

%% correlation between ilamb score and FCH4
clc;
close all;
clc;
close all;
vio_color = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];
benchmark = 1;
model_id = model_name_full;
model_included = strrep(model_id,'_','-');
% read in ilamb score
clear model_ranking ilamb_score threshold model_group
clf;
cdf_threshold = [.5:.1:.8];
var_name = ["Bias Score","RMSE Score","Seasonal Cycle Score","Overall Score"];
path = '/Volumes/LaCie/research/ILAMB_CH4/FCH4_benchmark/Sep2022/results/';
for region = 1:8
    if (region==1)
        loc = 'site_level_tier1';
    elseif (region==2)
        loc = 'site_level_tier1_north_of_30N';
    elseif (region==3)
        loc = 'site_level_tier1_south_of_30N';
    elseif (region==4)
        loc = 'site_level_tier1_bog';
    elseif (region==5)
        loc = 'site_level_tier1_fen';
    elseif (region==6)
        loc = 'site_level_tier1_marsh';
    elseif (region==7)
        loc = 'site_level_tier1_swamp';
    elseif (region==8)
        loc = 'site_level_tier1_tundra';
    end
    % read in bias, rmse, sesonality, and overall scores
    file=(['scalar_database.csv']);
    T = readtable([path loc '/_build/' file]);
    name = unique(T{:,4});
    % do not include ELM-v0 and ML for the comparison
    idx = find(contains(name, "ELM_v0")==1|contains(name, "Stanford")==1);
    name(idx) = [];
    if (region==1)
        model_list = unique(T{:,4});
        model_list(idx) = [];
        idx_TD_models = contains(model_list, "SURF")|contains(model_list, "GOSAT");
        TD_models = model_list(idx_TD_models);
        idx_BU_models = ~idx_TD_models;
        BU_models = model_list(idx_BU_models);
    end
    for jj = 1:length(model_list) % models
        idx = find(strcmpi(T{:,4},strtrim(model_list(jj))));
        if (length(idx)>0)
            tmp = T(idx,:);
            for kk = 1:4 % variables
                idx = find(strcmpi(tmp{:,5},var_name(kk)));
                ilamb_score(jj,kk, region) = str2double(tmp{idx,10});
            end
        else
            ilamb_score(jj,:, region) = nan;
        end
    end
end
for ii = 1:1
    if (ii==1)
        path = '/Volumes/LaCie/research/ILAMB_CH4/FCH4_benchmark/Sep2022/results/ML';
    end
    % read in bias, rmse, sesonality, and overall scores
    file=(['/_build/scalar_database.csv']);
    T = readtable([path file]);
    name = unique(T{:,4});
    % do not include ELM-v0 and Stanford for the comparison
    idx_excl = find(contains(name, "ELM_v0")==1|contains(name, "Stanford")==1);
%     name(idx_excl) = [];
    if (ii==1)
        model_list = unique(T{:,4});
        model_list(idx_excl) = [];
    end
    for jj = 1:length(model_list) % models
        idx = find(strcmpi(T{:,4},strtrim(model_list(jj))));
        if (length(idx)>0)
            tmp = T(idx,:);
            for kk = 1:4 % variables
                idx = find(strcmpi(tmp{:,5},var_name(kk)));
                ilamb_score(jj,kk,ii+8) = str2double(tmp{idx,10});
            end
        else
            ilamb_score(jj,:,ii+8) = nan;
        end
    end
end

% read in global wetland FCH4
idx_BU = 1:2:9;
idx_TD = 2:2:10;
fch4_BU_top20 = [];
fch4_TD_top20 = [];
for ii = 1:1
    clf;
    if (ii==1)
        filename=['./GMB/FCH4_budget_WAD2M_025deg_site.csv'];
        fig_type = 'obs';
    end
    % read in texts and values as a table
    T = readtable(filename);
    fch4_BU = T{:,idx_BU};
    fch4_TD = T{:,idx_TD};
    fch4_BU_top20 = [fch4_BU_top20; fch4_BU(:,end)];
    fch4_TD_top20 = [fch4_TD_top20; fch4_TD(:,end)];
    if (ii==1)
        fch4_BU_all = fch4_BU(:,1);
        fch4_TD_all = fch4_TD(:,1);
    end
end
desp = {'Obs, global', ['Obs, > 30' degree ' N'], ['Obs, < 30' degree ' N'],...
    'Obs, bog','Obs, fen','Obs, marsh','Obs, swamp','Obs, tundra', 'Machine-learning'};
R2_BU = [];
R2_TD = [];
for ii = 1:9
    subplot(3, 3, ii)
    xtmp = squeeze(ilamb_score(idx_BU_models,4,ii));
    ytmp = fch4_BU_all(isfinite(fch4_BU_all));
    idx = isfinite(xtmp)&isfinite(ytmp);
    h_BU = scatter(xtmp(idx), ytmp(idx), 'filled',...
            'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerFaceAlpha', 0.8);
    hold on
    rg_coef = polyfit(xtmp(idx), ytmp(idx), 1);
    fff = polyval(rg_coef,[min(xtmp(idx)), max(xtmp(idx))]);
    R2_BU = [R2_BU; corr(ytmp(idx), polyval(rg_coef,xtmp(idx)))^2];
    p = plot([min(xtmp(idx)), max(xtmp(idx))], fff, ...
        'color', [0 0.4470 0.7410],'linewidth', 2);
    p.Color(4) = 0.3;
%     TD results
    xtmp = squeeze(ilamb_score(idx_TD_models,4,ii));
    ytmp = fch4_TD_all(isfinite(fch4_TD_all));
    idx = isfinite(xtmp)&isfinite(ytmp);
    h_TD = scatter(xtmp(idx), ytmp(idx), 'filled',...
            'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerFaceAlpha', 0.8);
    rg_coef = polyfit(xtmp(idx), ytmp(idx), 1);
    fff = polyval(rg_coef,[min(xtmp(idx)), max(xtmp(idx))]);
    R2_TD = [R2_TD; corr(ytmp(idx), polyval(rg_coef,xtmp(idx)))^2];
    p = plot([min(xtmp(idx)), max(xtmp(idx))], fff, ...
        'color', [0.8500 0.3250 0.0980],'linewidth', 2);
    p.Color(4) = 0.3;
    t = text(0.02,1.15,[char(erase(id_name(ii),'"')) ' ' desp{ii}],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top',...
            'FontSize',12);
    t = text(0.552,0.32,['R^2 = ' num2str(sprintf('%1.2f',R2_BU(ii)))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top',...
            'FontSize',12,'color', [0 0.4470 0.7410]);
    t = text(0.552,0.20,['R^2 = ' num2str(sprintf('%1.2f',R2_TD(ii)))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top',...
            'FontSize',12,'color', [0.8500 0.3250 0.0980]);
    ylim([120 210])
    xlim([0.2 1])
    box on
end
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1) % sets ax1 to current axes
% universal left Y-axis
t = text(0.02,0.23,['Global wetland CH_4 emissions (TgCH_4 yr^{-1})'],'Units', 'Normalized',...
        'VerticalAlignment', 'Top','FontSize',12);
t.Rotation = 90;
% universal X-axis
t = text(0.41,0.05,['ILAMB overall scores'],'Units', 'Normalized',...
        'VerticalAlignment', 'Top','FontSize',12);
lgd = legend([h_BU, h_TD], {'BU models', 'TD models'},...
    'orientation','horizontal', 'box','off');
lgd.FontSize = 12;
lgd.Position = [0.63 0.01 0.3643 0.0452];

print('-dpng','-r300',...
    ['./FCH4_ilamb_cor_dist.png']);

%% evaluate refined FCH4 inferred from different benchmarks with bootstrapping
clc;
close all;
vio_color = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];
idx_BU = 1:2:9;
idx_TD = 2:2:10;
fch4_BU_top20 = [];
fch4_TD_top20 = [];
fch4_BU_top30 = [];
fch4_TD_top30 = [];
for ii = 1:8
    clf;
    if (ii==1)
        filename=['./GMB/FCH4_budget_WAD2M_025deg_site.csv'];
    elseif (ii==2)
        filename=['./GMB/FCH4_budget_WAD2M_025deg_site_north_of_30N.csv'];
    elseif (ii==3)
        filename=['./GMB/FCH4_budget_WAD2M_025deg_site_south_of_30N.csv'];
    elseif (ii==4)
        filename=['./GMB/FCH4_budget_WAD2M_025deg_site_bog.csv'];
    elseif (ii==5)
        filename=['./GMB/FCH4_budget_WAD2M_025deg_site_fen.csv'];
    elseif (ii==6)
        filename=['./GMB/FCH4_budget_WAD2M_025deg_site_marsh.csv'];
    elseif (ii==7)
        filename=['./GMB/FCH4_budget_WAD2M_025deg_site_swamp.csv'];
    elseif (ii==8)
        filename=['./GMB/FCH4_budget_WAD2M_025deg_site_tundra.csv'];
    end
    % read in texts and values as a table
    T = readtable(filename);
    fch4_BU = T{:,idx_BU};
    fch4_TD = T{:,idx_TD};
    idx = isfinite(fch4_BU(:,end));
    fch4_BU_top20(1:sum(idx),ii) = fch4_BU(idx,end);
    idx = isfinite(fch4_TD(:,end));
    fch4_TD_top20(1:sum(idx),ii) = fch4_TD(idx,end);
    if (ii==1)
        fch4_BU_all = fch4_BU(:,1);
        fch4_TD_all = fch4_TD(:,1);
    end
end
% clear tmp
fch4_BU_top20(fch4_BU_top20==0) = nan;
fch4_TD_top20(fch4_TD_top20==0) = nan;
% fch4_BU = fch4_BU_top20(isfinite(fch4_BU_top20));
% fch4_TD = fch4_TD_top20(isfinite(fch4_TD_top20));
tmp = nan(length(fch4_TD_all(:)),18);
tmp(1:length(fch4_BU_all),1) = fch4_BU_all;
tmp(1:length(fch4_TD_all),10) = fch4_TD_all;
for ii = 2:9
    idx = isfinite(fch4_BU_top20(:,ii-1));
    tmp(1:sum(isfinite(fch4_BU_top20(:,ii-1))),ii) = fch4_BU_top20(idx,ii-1);
end
for ii = 11:18
    idx = isfinite(fch4_TD_top20(:,ii-10));
    tmp(1:sum(isfinite(fch4_TD_top20(:,ii-10))),ii) = fch4_TD_top20(idx,ii-10);
end
subplot(3, 1, 1)
% data points
idx_BU = 1:9;
idx_TD = 10:18;
x_BU = 1:3:25;
x_TD = 2:3:26;
xxx = 3:3:27;
for n = 1:9
    idx = isfinite(tmp(:,idx_BU(n)));
    scatter(ones(sum(idx),1)*x_BU(n), tmp(idx,idx_BU(n)), 'filled',...
        'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', 'none');
    hold on
    idx = isfinite(tmp(:,idx_TD(n)));
    scatter(ones(sum(idx),1)*x_TD(n), tmp(idx,idx_TD(n)), 'filled',...
        'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', 'none');
end
% BU range
for n = 1:9
    patch([x_BU(n)-.2 x_BU(n)+.2 x_BU(n)+.2 x_BU(n)-.2],...
    [min(tmp(:,idx_BU(n))) min(tmp(:,idx_BU(n))) max(tmp(:,idx_BU(n))) max(tmp(:,idx_BU(n)))],...
    [0, 0.4470, 0.7410],'FaceAlpha',.3,'LineStyle','none')
end
% TD range
for n = 1:9
    patch([x_TD(n)-.2 x_TD(n)+.2 x_TD(n)+.2 x_TD(n)-.2],...
    [min(tmp(:,idx_TD(n))) min(tmp(:,idx_TD(n))) max(tmp(:,idx_TD(n))) max(tmp(:,idx_TD(n)))],...
    [0.8500, 0.3250, 0.0980],'FaceAlpha',.3,'LineStyle','none')
end
box on
% broad line
for n = 1:9
    plot([xxx(n) xxx(n)], [120 210], 'k-')
end
plot(x_BU, mean(tmp(:,idx_BU), 'omitnan'), 'x', 'Color', [0, 0.4470, 0.7410],...
        'MarkerSize', 10,'linewidth', 1, 'MarkerEdgeColor', [0, 0.4470, 0.7410]);
plot(x_TD, mean(tmp(:,idx_TD), 'omitnan'), 'x', 'Color', [0.8500, 0.3250, 0.0980],...
        'MarkerSize', 10,'linewidth', 1, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980]);
ylabel({'Global wetland CH_4'; 'emissions (TgCH_4 yr^{-1})'}, 'FontSize', 12)
% xticks(sort([x_BU, x_TD]))
% xticklabels({'BU models', 'TD models', 'BU models', 'TD models', ...
%     'BU models', 'TD models', 'BU models', 'TD models', 'BU models', 'TD models'})
xticks([x_BU+x_TD]/2)
xticklabels({'No constraint','Obs, global', ['Obs, > 30' degree ' N'], ['Obs, < 30' degree ' N'],...
    'Obs, bog','Obs, fen','Obs, marsh','Obs, swamp','Obs, tundra'});
xtickangle(30)
ylim([120 210])
xlim([0 max(xxx)])

for ii = 1:1
    if (ii==1)
        fch4_BU = fch4_BU_top20(:);
        fch4_TD = fch4_TD_top20(:);
        desp = 'Top 20% models';
    elseif (ii==2)
        fch4_BU = fch4_BU_top30(:);
        fch4_TD = fch4_TD_top30(:);
        desp = 'Top 30% models';
    end
    fch4_BU = fch4_BU(isfinite(fch4_BU));
    fch4_TD = fch4_TD(isfinite(fch4_TD));
    fch4_BU_dem = fch4_BU_all(isfinite(fch4_BU_all));
    fch4_TD_dem = fch4_TD_all(isfinite(fch4_TD_all));
    subplot(3, 1, 2)
    yyaxis left
    h_BU = histogram(fch4_BU, [120:5:220]);
%     set(gca,'xtick',[0:.2:1])
%     xlabel('Global wetland CH_4 emissions (TgCH_4 yr^{-1})', 'FontSize', 12)
%     ylabel('Frequency', 'FontSize', 12)
%     hhh.BinLimits = [0 1];
    h_BU.FaceColor = [0 0.4470 0.7410];
    h_BU.EdgeColor = 'none';
    h_BU.FaceAlpha = 0.3;
    hold on
    % TD
    h_TD = histogram(fch4_TD, [120:5:220]);
%     set(gca,'xtick',[0:.2:1])
    xlabel('Global wetland CH_4 emissions from top 20% models (TgCH_4 yr^{-1})', 'FontSize', 12)
    ylabel('Frequency', 'FontSize', 12)
%     hhh.BinLimits = [0 1];
    h_TD.FaceColor = [0.8500 0.3250 0.0980];
    h_TD.EdgeColor = 'none';
    h_TD.FaceAlpha = 0.3;
    
    ax1 = gca;
%     ax1.XLim = [0 1];
    ax1.YColor = 'k';
    yyaxis right
    % model meritocracy
    [f,xi] = ksdensity(fch4_BU,...
                'Function','cdf','Support','positive');
    % PDF distribution
    p1 = plot(xi(2:end), diff(f), '-','color',[0 0.4470 0.7410],...
        'linewidth',3);
    hold on
    tmp = diff(f);
    idx = find(tmp==max(tmp));
    ymax = tmp(idx);
    plot([xi(idx), xi(idx)], [0 tmp(idx)], ':',...
        'color',[0 0.4470 0.7410],'linewidth',2)
    text(xi(idx)+1, 0.003, sprintf('%3.2f',xi(idx)), 'color', [0 0.4470 0.7410])
    [f,xi] = ksdensity(fch4_TD,...
                'Function','cdf','Support','positive');
    p2 = plot(xi(2:end), diff(f), '-','color',[0.8500 0.3250 0.0980],...
        'linewidth',3);
    tmp = diff(f);
    idx = find(tmp==max(tmp));
    ymax = max(ymax,tmp(idx));
    plot([xi(idx), xi(idx)], [0 tmp(idx)], ':',...
        'color',[0.8500 0.3250 0.0980],'linewidth',2)
    text(xi(idx)+1, 0.003, sprintf('%3.2f',xi(idx)), 'color', [0.8500 0.3250 0.0980])
%     axis tight
    xlim([120 220])
    ylim([0 ymax+0.005])
%     text(xi(idx)+2, 0.002, sprintf('%3.2f',xi(idx)), 'color', [0.8500 0.3250 0.0980])
    ylabel('Estimated PDF', 'FontSize', 12)
    ax1 = gca;
%     ax1.XLim = [0 1];
    ax1.YColor = 'k';
    id = char(erase(id_name(2),'"'));
    t = text(0.02,0.98,[id],'Units',...
            'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12);
end
% boostraping mean
tmp = fch4_BU_top20(:);
tmp(~isfinite(tmp)) = [];
rng default  % For reproducibility
stats_BU = bootstrp(1000,@(x)[mean(x) std(x)],tmp);
tmp = fch4_TD_top20(:);
tmp(~isfinite(tmp)) = [];
rng default  % For reproducibility
stats_TD = bootstrp(1000,@(x)[mean(x) std(x)],tmp);
subplot(3, 1, 3)
yyaxis left
h_BU = histogram(stats_BU(:,1), [120:5:220]);
h_BU.FaceColor = [0 0.4470 0.7410];
h_BU.EdgeColor = 'none';
h_BU.FaceAlpha = 0.3;
hold on
% TD
h_TD = histogram(stats_TD(:,1), [120:5:220]);
xlabel('Global wetland CH_4 emissions from top 20% models (TgCH_4 yr^{-1})', 'FontSize', 12)
ylabel('Frequency', 'FontSize', 12)
h_TD.FaceColor = [0.8500 0.3250 0.0980];
h_TD.EdgeColor = 'none';
h_TD.FaceAlpha = 0.3;
ax1 = gca;
%     ax1.XLim = [0 1];
ax1.YColor = 'k';
yyaxis right
% model meritocracy
[f,xi] = ksdensity(stats_BU(:,1),...
            'Function','cdf','Support','positive');
% PDF distribution
p1 = plot(xi(2:end), diff(f), '-','color',[0 0.4470 0.7410],...
    'linewidth',3);
hold on
tmp = diff(f);
idx = find(tmp==max(tmp));
ymax = tmp(idx);
plot([xi(idx), xi(idx)], [0 tmp(idx)], ':',...
    'color',[0 0.4470 0.7410],'linewidth',2)
text(xi(idx)+6, 0.02, sprintf('%3.2f',xi(idx)), 'color', [0 0.4470 0.7410])
[f,xi] = ksdensity(stats_TD(:,1),...
            'Function','cdf','Support','positive');
p2 = plot(xi(2:end), diff(f), '-','color',[0.8500 0.3250 0.0980],...
    'linewidth',3);
tmp = diff(f);
idx = find(tmp==max(tmp));
ymax = max(ymax, tmp(idx));
plot([xi(idx), xi(idx)], [0 tmp(idx)], ':',...
    'color',[0.8500 0.3250 0.0980],'linewidth',2)
text(xi(idx)+6, 0.02, sprintf('%3.2f',xi(idx)), 'color', [0.8500 0.3250 0.0980])
axis tight
xlim([120 220])
ylim([0 ymax+0.005])
%     text(xi(idx)+2, 0.002, sprintf('%3.2f',xi(idx)), 'color', [0.8500 0.3250 0.0980])
ylabel('Estimated PDF', 'FontSize', 12)
ax1 = gca;
%     ax1.XLim = [0 1];
ax1.YColor = 'k';
id = char(erase(id_name(3),'"'));
t = text(0.02,0.98,[id],'Units',...
        'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12);

% lgd = legend([h_BU, h_TD], {'BU models', 'TD models'},...
%     'orientation','horizontal', 'box','off');
% lgd.FontSize = 12;
% lgd.Position = [0.54 0.455 0.3643 0.0452];
% % make an invisible axis handle
% ax_copy = axes('Position',get(gca,'Position'),'Visible','Off');
lgd = legend([h_BU h_TD],...
            {'BU models', 'TD models'},...
            'orientation','horizontal', 'box','off');
lgd.FontSize = 12;
lgd.Position = [0.543 0.93 0.3643 0.0452];
        
print('-dpng','-r300',...
    ['./FCH4_refined.png']);
