%% this program sort model scores outputted by ilamb
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

% read in CRUJRA forcing
path = '/Volumes/LaCie/research/ILAMB_CH4/FCH4_benchmark/forcing/';
filename = [path 'CRUJRA_monthly_Tair.nc'];
CRUJRA_Tair = ncread(filename, 'Tair')-273.15;
CRUJRA_lat = ncread(filename, 'lat');
CRUJRA_lon = ncread(filename, 'lon');
CRUJRA_time = ncread(filename, 'time');
path = '/Volumes/LaCie/research/ILAMB_CH4/FCH4_benchmark/forcing/';
% CRUJRA time lable
idx_CRUJRA_year = [];
idx_CRUJRA_month = [];
for yr = 2008:2017
    idx_CRUJRA_year = [idx_CRUJRA_year; yr*ones(12,1)];
    tmp = 1:12;
    idx_CRUJRA_month = [idx_CRUJRA_month; tmp'];
end

% sites included
meta_data = readtable('/Volumes/LaCie/research/ILAMB_CH4/FCH4_benchmark/forcing/Powell_site_LATLON_V3.xlsx');
siteID_modeling = string(erase(meta_data{:,2},'-'));

% meta data for measurements
meta_data_obs = readtable('/Volumes/LaCie/research/ILAMB_CH4/FCH4_benchmark/forcing/Site list and metadata.xlsx');
siteID_obs = string(erase(meta_data_obs{:,2},'-'));
site_lat_lon = [meta_data_obs{:, 7}, meta_data_obs{:, 8}];

% obs time lable
idx_obs_year = [];
idx_obs_month = [];
for yr = 2006:2018
    idx_obs_year = [idx_obs_year; yr*ones(12,1)];
    tmp = 1:12;
    idx_obs_month = [idx_obs_month; tmp'];
end

% read in FLUXNET-CH4 dataset
path = '../wetlands_tier1/';
file = folderFiles([path],...
        '*V2.csv');
% variable names for analysis
var_wanted = ["Year","DOY","FCH4_F_ANN_mean","TA_F","P_F"];
[rows, cols] = size(file);
daily_aggregate = [];
wetland_lat_lon = [];
for site = 1:rows
    T = readtable([path file(site,:)]);
    % read in the data header
    fileID = fopen([path file(site,:)]);
    tline = fgets(fileID);
    fclose(fileID);
    tline = strsplit(tline,',');
    header = string(erase(tline,'"'));
    daily_tmp = [];
    for variables = 1:length(var_wanted) % loop through variables
        if (sum(contains(header,var_wanted(variables)))>=1) 
            % if the requested variable exists
            idx = find(strcmpi(header,var_wanted(variables)));
            tmp = T{:,idx};
            tmp(tmp==-9999) = nan;
            % variable will be stored as a cell array if it mixes texts &
            % numerics
            if (iscell(tmp))
                % find nan values in the array
                idx_NA=find(contains(tmp,'NA'));
                idx_full = logical(1:length(tmp))';
                idx_full(idx_NA) = 0;
                xxx = find(idx_full==1);
                yyy = find(idx_full==0);
                % assign the numeric and sting extracted from the cell
                % array to a double array with nan
%                 zzz = tmp(idx_full);
                dim = size(tmp);
                data = zeros(dim);
                data(xxx) = cellfun(@str2num,tmp(xxx));
                data(yyy) = nan;
            else
                data = tmp;
            end
            daily_tmp(1:length(T{:,idx}), variables) = data;
        else % the requested data is not available
            daily_tmp(1:length(T{:,idx}), variables) = nan;
        end
        
    end
    % include site number
    daily_tmp(:,length(var_wanted)+1) = site;
    % exclude leap days
    idx = find(daily_tmp(:,2)==366);
    daily_tmp(idx,:) = [];
    % include month
    years = unique(daily_tmp(:,1));
    for yr = 1:length(years)
        idx = find(daily_tmp(:,1)==years(yr));
        if (length(idx)==365)
            daily_tmp(idx,length(var_wanted)+2) = month;
%         elseif (length(idx)==366)
%             daily_tmp(idx,length(var_wanted)+2) = month_leap;
        else
            idx = find(daily_tmp(:,1)==years(yr)&daily_tmp(:,2)<=mon_sum(1));
            daily_tmp(idx,length(var_wanted)+2) = 1;
            for m = 2:12
                idx = find(daily_tmp(:,1)==years(yr)&daily_tmp(:,2)>mon_sum(m-1)&daily_tmp(:,2)<=mon_sum(m));
                daily_tmp(idx,length(var_wanted)+2) = m;
            end      
        end
    end
    % lat and lon
    idx_obs = find(strcmpi(siteID_obs,file(site,1:5)));
    tmp = [site_lat_lon(idx_obs, :)];
    wetland_lat_lon = [wetland_lat_lon; tmp];
    daily_aggregate = [daily_aggregate; daily_tmp];
end
% calculate the number of site-years
site_years = [];
for site = 1:rows
    idx = find(daily_aggregate(:,6)==site);
    years = unique(daily_aggregate(idx, 1));
    for yr = 1:length(years)
        if (length(find(daily_aggregate(idx,1)==years(yr)))>30)
            tmp = [years(yr), site];
            site_years = [site_years; tmp];
        end
    end
end

%% read in ilamb score
clear model_ranking ilamb_score threshold model_group
clf;
cdf_threshold = [.7 .8 .9];
var_name = ["Bias Score","RMSE Score","Seasonal Cycle Score","Overall Score"];
for ii = 1:2
    if (ii==1)
        path = '/Volumes/LaCie/research/ILAMB_CH4/FCH4_benchmark/Sep2022/results/site_level_tier1';
    elseif (ii==2)
        path = '/Volumes/LaCie/research/ILAMB_CH4/FCH4_benchmark/Sep2022/results/ML';
    end
    % read in bias, rmse, sesonality, and overall scores
    file=(['/_build/scalar_database.csv']);
    T = readtable([path file]);
    name = unique(T{:,4});
    % do not include ELM-v0 for the comparison
    idx_ELMv0 = find(contains(name, "ELM_v0")==1);
    name(idx_ELMv0) = [];
    model_list = unique(T{:,4});
    model_list(idx_ELMv0) = [];
    for jj = 1:length(model_list) % models
        idx = find(strcmpi(T{:,4},strtrim(model_list(jj))));
        if (length(idx)>0)
            tmp = T(idx,:);
            for kk = 1:4 % variables
                idx = find(strcmpi(tmp{:,5},var_name(kk)));
                ilamb_score(jj,kk, ii) = str2double(tmp{idx,10});
%                 ilamb_score(jj,kk,ii) = (tmp{idx,10});
            end
        else
            ilamb_score(jj,:,ii) = nan;
        end
    end
    
    model_name{ii} = name;
    
end

%% compare the overall score inferred from BU and TD models, using 2008-2017
% categorize BU and TD models
clc;
close all;
clear threshold_TD_BU model_group_TD_BU
fig_header = {'model_rank_dist_obs_', 'model_rank_dist_ML_'};
for benchmark = 1:2
    model_id = model_name{benchmark};
    idx_TD = contains(model_id, "SURF")|contains(model_id, "GOSAT");
    TD_models = model_id(idx_TD);
    idx_ML = contains(model_id, "Stanford");
    BU_models = model_id(~idx_TD&~idx_ML);

    % individual model score, 2008-2017
    clf;
    model_included = strrep(model_id,'_','-');
    % rename UPCH4
    if (benchmark~=2)
        model_included{27}= 'UPCH4-0.25degree';
        model_included{28}= 'UPCH4-1degree';
    end
    % color_scheme = [3, 2, 4];
    vio_color = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];
    
    for ii = 4:4
        clf;
        p_agg = [];
        subplot(1, 2, 1)
        if (ii==1)
            desp = 'Sbias';
            label = ['Model bias score from ILAMB'];
        elseif (ii==2)
            desp = 'Srmse';
            label = ['Model RMSE score from ILAMB'];
        elseif (ii==3)
            desp = 'Sseason';
            label = ['Model seasonality score from ILAMB'];
        elseif (ii==4)
            desp = 'Soverall';
            label = ['Model overall score from ILAMB'];
        end
        
        idx = isfinite(ilamb_score(1:length(model_id),ii, benchmark));
        score = squeeze(ilamb_score(idx,ii, benchmark));
        name = model_included(idx);
        pos = 1:length(score);    
        xrange = max(score);
        [tmp,idx]=sort(score,'ascend');
        % categorize BU and TD models
        idx_TD = contains(name, "SURF")|contains(name, "GOSAT");
        pos_tmp = pos(idx_TD);
        idx_TD_loc = [];
        for jj = 1:length(pos_tmp)
            idx_TD_loc = [idx_TD_loc; find(idx==pos_tmp(jj))];
        end
        idx_ML = contains(name, "UPCH4");
        pos_tmp = pos(idx_ML);
        idx_ML_loc = [];
        for jj = 1:length(pos_tmp)
            idx_ML_loc = [idx_ML_loc; find(idx==pos_tmp(jj))];
        end
        pos_tmp = pos(~idx_TD&~idx_ML);
        idx_BU_loc = [];
        for jj = 1:length(pos_tmp)
            idx_BU_loc = [idx_BU_loc; find(idx==pos_tmp(jj))];
        end

        h_BU = barh(pos(idx_BU_loc), tmp(idx_BU_loc), .5);
        h_BU.LineStyle = 'none';
        hold on
        h_TD = barh(pos(idx_TD_loc),tmp(idx_TD_loc), .5);
        h_TD.LineStyle = 'none';
        if (benchmark~=2)
            h_ML = barh(pos(idx_ML_loc),tmp(idx_ML_loc), .5);
            h_ML.LineStyle = 'none';
            h_ML.BarWidth = 0.5;
        end
%         if (benchmark==3)
%             h_ML.BarWidth = 0.25;
%         end
        xlabel(label, 'FontSize', 12)
        set(gca,'ytick',[1:length(tmp)])
        set(gca,'yticklabel',[name(idx)])
        ylim([.5 length(score)+.5])
        xlim([0 xrange])
        subplot(1, 2, 2)
        if (benchmark~=2)
            vio = nan(max(length(idx_BU_loc), length(idx_TD_loc)), 3);
            vio(1:length(idx_BU_loc), 1) = tmp(idx_BU_loc);
            vio(1:length(idx_TD_loc), 2) = tmp(idx_TD_loc);
            vio(1:length(idx_ML_loc), 3) = tmp(idx_ML_loc);
        else
            vio = nan(max(length(idx_BU_loc), length(idx_TD_loc)), 2);
            vio(1:length(idx_BU_loc), 1) = tmp(idx_BU_loc);
            vio(1:length(idx_TD_loc), 2) = tmp(idx_TD_loc);
        end
        median_value = median(vio,'omitnan');
        std_value = std(vio,'omitnan');
        ppp = violinplot(vio);
        for j = 1:length(ppp)
            ppp(j).BoxColor = 'k';
            ppp(j).BoxWidth = 0.04;
            ppp(j).ShowData = 1;
            ppp(j).ShowMean = 1;
            ppp(j).MeanPlot.LineWidth = 2;
            ppp(j).ViolinColor = vio_color(j,:);
            ppp(j).MedianPlot.SizeData = 16;
        end
        box on
        ylabel(label, 'FontSize', 12)
        if (benchmark~=5)
            set(gca,'xtick',[1:3])
            set(gca,'xticklabel',{'BU models', 'TD models', 'ML models'},...
                'FontSize',12)
        else
            set(gca,'xtick',[1:2])
            set(gca,'xticklabel',{'BU models', 'TD models'},...
                'FontSize',12)
        end
        set(gca,'xticklabelRotation',30)
        ylim([0 1])
        id = char(erase(id_name(2),'"'));
        t = text(0.02,0.98,[id],'Units',...
                'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12);    
        x_pos = linspace(0.01, 0.72, length(ppp));
        if (benchmark==1)
            ylim([0.3 0.7])
        end
        for j = 1:length(ppp)
            t = text(x_pos(j),0.78,[num2str(sprintf('%2.2f',median_value(j))) char(177) ...
                num2str(sprintf('%2.1f',std_value(j)))],...
                'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
        end
        ax1 = axes('Position',[0 0 1 1],'Visible','off');
        axes(ax1) % sets ax1 to current axes
        if (benchmark~=2)
            lgd = legend([h_BU, h_TD, h_ML],...
                        {'BU models', 'TD models', 'ML models'},...
                        'orientation','horizontal', 'box','off');
        else
            lgd = legend([h_BU, h_TD],...
                {'BU models', 'TD models'},...
                'orientation','horizontal', 'box','off');
        end
        lgd.FontSize = 12;
        lgd.Position = [0.14 0.93 0.3643 0.0452];
        print('-dpng','-r300',['./' char(fig_header{benchmark}) desp '.png']);
    end
end
%% synthesize model ranking inferred from different scoring system
% based on site-level measurements
clf;
benchmark = 1;
model_id = model_name{benchmark};
model_included = strrep(model_id,'_','-');
model_included{27}= 'UPCH4-0.25degree';
model_included{28}= 'UPCH4-1degree';
tmp = squeeze(ilamb_score(:,:,benchmark));
hhh = image(tmp,'CDataMapping','scaled');
% set nan as blank
set(hhh,'AlphaData',~isnan(tmp))
% grad=colorGradient([1 1 1], [0.0 0.27 0.13],64);
grad=colorGradient([1 1 1], [0 0.5 0],64);
colormap(grad);
caxis([0 1.0])
c = colorbar('NorthOutside', 'Position', [0.28, 0.93, 0.6, 0.03]);
% overlay grid lines
[rows, columns, numberOfColorChannels] = size(tmp);
hold on;
lineSpacing = 1;
for col = .5 : lineSpacing : columns+.5
    line([col, col], [.5, rows+.5], 'Color', 'k', 'LineWidth', 1);
end
% indicate the top 5 models
idx = [];
score = tmp;
for ii = 1:5
    idx_tmp = [];
    if (ii==1)
        desp = '1^{st}';
    elseif (ii==2)
        desp = '2^{nd}';
    elseif (ii==3)
        desp = '3^{rd}';
    elseif (ii==4)
        desp = '4^{th}';
    elseif (ii==5)
        desp = '5^{th}';
    end
    for col = 1:columns
        idx_tmp = [idx_tmp, find(score(:,col)==max(score(:,col)))];
        score(idx_tmp(end), col) = nan;
        text(col, idx_tmp(end), desp);
    end
    idx = [idx; idx_tmp];
end
ax = ancestor(hhh, 'axes');
% Change properties of the axes
ax.YTick = [1:rows];
ax.YTick = [1:rows];
ax.YTickLabel = model_included;
ax.XTick = [1:columns];
ax.XTick = [1:columns];
ax.XTickLabel = {'Bias score','RMSE score','Seasonality score',...
        'Overall score'};
ax.YAxis.FontSize = 11;
ax.XAxis.FontSize = 12;
ax.XAxis.TickLabelRotation = 20;
% set colorbar description
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1) % sets ax1 to current axes
% set colorbar style
text(0.89,0.97,['ILMAB score'],'Units', 'Normalized', 'VerticalAlignment', 'Top')
print('-dpng','-r300',['./model_rank_dist_sen_site.png']);

%% aggregate obs benchmarks with ML
% read in ilamb score from obs runs
clear model_ranking ilamb_score_obs threshold model_group
clf;
cdf_threshold = [.7 .8 .9];
var_name = ["Bias Score","RMSE Score","Seasonal Cycle Score","Overall Score"];
path = '/Volumes/LaCie/research/ILAMB_CH4/FCH4_benchmark/Sep2022/results/';
for subset = 1:8
    if (subset==1)
        loc = 'site_level_tier1';
    elseif (subset==2)
        loc = 'site_level_tier1_north_of_30N';
    elseif (subset==3)
        loc = 'site_level_tier1_south_of_30N';
    elseif (subset==4)
        loc = 'site_level_tier1_bog';
    elseif (subset==5)
        loc = 'site_level_tier1_fen';
    elseif (subset==6)
        loc = 'site_level_tier1_marsh';
    elseif (subset==7)
        loc = 'site_level_tier1_swamp';
    elseif (subset==8)
        loc = 'site_level_tier1_tundra';
    end
    % read in bias, rmse, sesonality, and overall scores
    file=(['scalar_database.csv']);
    T = readtable([path loc '/_build/' file]);
    name = unique(T{:,4});
    % do not include ELM-v0 and ML for the comparison
    idx = find(contains(name, "ELM_v0")==1);
    name(idx) = [];
    if (subset==1)
        model_list = unique(T{:,4});
        model_list(idx) = [];
    end
    for jj = 1:length(model_list) % models
        idx = find(strcmpi(T{:,4},strtrim(model_list(jj))));
        if (length(idx)>0)
            tmp = T(idx,:);
            for kk = 1:4 % variables
                idx = find(strcmpi(tmp{:,5},var_name(kk)));
                ilamb_score_obs(jj,kk, subset) = str2double(tmp{idx,10});
%                 if (subset>3)
%                     ilamb_score_obs(jj,kk, subset) = str2double(tmp{idx,10});
%                 else
%                     ilamb_score_obs(jj,kk, subset) = (tmp{idx,10});
%                 end
            end
        else
            ilamb_score_obs(jj,:, subset) = nan;
        end
    end
end
clf;
benchmark = 1;
model_id = model_name{benchmark};
model_included = strrep(model_id,'_','-');
model_included{27}= 'UPCH4-0.25degree';
model_included{28}= 'UPCH4-1degree';
tmp = [squeeze(ilamb_score_obs(:,4,:)), nan(length(model_id),1), squeeze(ilamb_score(:,4,end))];
hhh = image(tmp,'CDataMapping','scaled');
% set nan as blank
set(hhh,'AlphaData',~isnan(tmp))
% grad=colorGradient([.5 0 0], [0 0 .5], 64);
grad=colorGradient([1 1 1], [0 0.5 0],64);
colormap(grad);
% colormap summer
% colormap(jet(256))
caxis([0.1 0.70])
c = colorbar('Eastoutside', 'Position', [0.92, 0.17, 0.03, 0.7]);
% overlay grid lines
[rows, columns, numberOfColorChannels] = size(tmp);
hold on;
lineSpacing = 1;
for col = .5 : lineSpacing : columns+.5
    line([col, col], [.5, rows+.5], 'Color', 'k', 'LineWidth', 1);
end
% indicate the top 5 models
idx = [];
score = tmp;
% col_val = 1:10;
% col_val(col_val==9) = [];
for ii = 1:5
    idx_tmp = [];
    if (ii==1)
        desp = '1^{st}';
    elseif (ii==2)
        desp = '2^{nd}';
    elseif (ii==3)
        desp = '3^{rd}';
    elseif (ii==4)
        desp = '4^{th}';
    elseif (ii==5)
        desp = '5^{th}';
    end
    for col = 1:columns
        if (col~=9)
            idx_tmp = [idx_tmp, find(score(:,col)==max(score(:,col)))];
            score(idx_tmp(end), col) = nan;
            text(col, idx_tmp(end), desp);
        end
    end
    idx = [idx; idx_tmp];
end
ax = ancestor(hhh, 'axes');
% Change properties of the axes
% ax.YTick = [1:rows];
ax.YTick = [1:rows];
ax.YTickLabel = model_included;
% get the current tick labeks
ticklabels = get(gca,'YTickLabel');
% prepend a color for each tick label
ticklabels_new = cell(size(ticklabels));
% categorize BU and TD models
idx_full_TD = contains(model_included, "SURF")|contains(model_included, "GOSAT");
idx_full_ML = contains(model_included, "UPCH4");
idx_full_BU = ~idx_full_TD&~idx_full_ML;
for i = 1:length(ticklabels)
    if (idx_full_TD(i)==1)
        ticklabels_new{i} = ['\color{red} ' ticklabels{i}];
    elseif (idx_full_BU(i)==1)
        ticklabels_new{i} = ['\color{blue} ' ticklabels{i}];
    elseif (idx_full_ML(i)==1)
        ticklabels_new{i} = ['\color{black} ' ticklabels{i}];
    end
end
% set the tick labels
set(gca, 'YTickLabel', ticklabels_new);
ax.XTick = [1:columns];
ax.XTickLabel = {'Obs (global)', ['Obs (>30' degree 'N)'], ['Obs (<30' degree 'N)'], ...
    'Obs (bog)', 'Obs (fen)', 'Obs (marsh)', 'Obs (swamp)', 'Obs (tundra)', ' ',...
    'Machine-learning'};
ax.YAxis.FontSize = 10;
ax.XAxis.FontSize = 12;
ax.XAxis.TickLabelRotation = 20;
% set colorbar description
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1) % sets ax1 to current axes
% set colorbar style
text(0.94,0.93,{'ILAMB'; 'score'},'Units', 'Normalized', ...
    'VerticalAlignment', 'Top', 'HorizontalAlignment', 'center')
% text(0.89,0.99,['ILAMB score'],'Units', 'Normalized', 'VerticalAlignment', 'Top')
% assign data description
text(0.41,0.9805,'Site observations','Units', 'Normalized', ...
    'VerticalAlignment', 'Top','Fontsize',12)
text(0.81,0.9805,'Global-scale','Units', 'Normalized', ...
    'VerticalAlignment', 'Top','Fontsize',12)
% arrow 1
y1 = [.94 .94];
x1 = [.22 .77];
annotation('doublearrow',x1,y1,'Color','k');
% % arrow 2
y1 = [.94 .94];
x1 = [.83 .91];
annotation('doublearrow',x1,y1,'Color','k');
print('-dpng','-r300',['./model_rank_dist_sen_ben_color_v6.png']);

%% filtering out better performance models using ilamb overall score
% categorize BU and TD models
clc;
figname = {'constraint_site', 'constraint_ML'};
for benchmark = 1:2
    close all;
    clear threshold_TD_BU model_group_TD_BU
    model_id = model_name{benchmark};
    idx_TD = contains(model_id, "SURF")|contains(model_id, "GOSAT");
    TD_models = model_id(idx_TD);
    idx_ML = contains(model_id, "Stanford");
    BU_models = model_id(~idx_TD&~idx_ML);
    cdf_threshold = [.5:.1:.9];
    for ii = 1:2
        subplot(2, 1, ii)
        if (ii==1)
            desp = 'BU models';
            score = squeeze(ilamb_score(~idx_TD&~idx_ML,end,benchmark));
        elseif (ii==2)
            desp = 'TD models';
            score = squeeze(ilamb_score(idx_TD,end,benchmark));
        end
        yyaxis left
        hhh = histogram(score, [0:.05:1]);
        set(gca,'xtick',[0:.2:1])
        xlabel('Model overall score from ILAMB', 'FontSize', 12)
        ylabel('Frequency', 'FontSize', 12)
        hhh.BinLimits = [0 1];
        hhh.FaceColor = 'k';
        hhh.EdgeColor = 'none';
        ax1 = gca;
        ax1.XLim = [0 1];
        ax1.YColor = 'k';
        hold on
        [f,xi] = ksdensity(score,...
                    'Function','pdf','Support','positive');
        % PDF distribution
        p1 = plot(xi, f, 'k-',...
            'linewidth',2);
        p1.Color(4) = 0.2;
        yyaxis right
        [f,xi] = ksdensity(score,...
                    'Function','cdf','Support','positive');
        % CDF distribution
        p2 = plot(xi, f,...
            'linewidth',2, 'color',[5 128 0]/255);
        ylabel('Estimated CDF', 'FontSize', 12)
        ax2 = gca;
        ax2.YColor = [5 128 0]/255;
        hold on
        for jj = 1:length(cdf_threshold)
            tmp = abs(f-cdf_threshold(jj));
            idx = find(tmp==min(tmp));
            ppp = plot([xi(idx) xi(idx)], [0 f(idx)], ':', ...
                     'linewidth',2, 'color',[5 128 0]/255);
            ppp.Color(4) = 0.4;
            ppp = plot([xi(idx) 1], [cdf_threshold(jj) cdf_threshold(jj)], ':', ...
                     'linewidth',2, 'color',[5 128 0]/255);
            ppp.Color(4) = 0.4;
            threshold_TD_BU(ii, jj) = length(find(score>xi(idx)));
            model_group_TD_BU{ii,jj} = find(score>xi(idx));
        end
        yticks([0 .5:.1:1])
        id = char(erase(id_name(ii),'"'));
        t = text(0.02,1.1,[id ' ' desp],'Units',...
                'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12);
    end
    print('-dpng','-r300',['./' char(figname{benchmark}) '.png']);
end
