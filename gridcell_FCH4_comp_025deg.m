clear all
close all
clc;
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
    "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
    "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];
% read in wetland fraction from WAD2M ver2
file = '/Volumes/LaCie/research/ILAMB_CH4/datasets/regrid_720x1440_WAD2M_wetlands_2000-2020_025deg_Ver2.0.nc';
% WAD2M, 2000-2018
WAD2M = ncread(file, 'Fw');
% extract 2008-2017
WAD2M(:, :, 1:96) = [];
WAD2M(:, :, end-35:end) = [];
WAD2M(WAD2M<=0) = NaN;
% generate global area
WAD2M_lat = ncread(file, 'lat');
WAD2M_lon = ncread(file, 'lon');
[LAT, LON] = meshgrid(WAD2M_lat, WAD2M_lon);
land_area = cdtarea(LAT, LON); 

% read in benchmark dataset
mon = [31,28,31,30,31,30,31,31,30,31,30,31];
years = linspace(2008, 2017, 10);
idx_years = linspace(2008, 2018, 121);
idx_years(end) = [];
fch4_BU = [];
fch4_TD = [];
fch4_ML = [];
fch4_BU_idx = [];
fch4_TD_idx = [];
fch4_ML_idx = [];
fch4_all = [];
% read in saturated FCH4 from wetlands
path = '/Volumes/LaCie/research/ILAMB_CH4/datasets/Aug_2022/FCH4_per_gridcell_area_720x1440/';
filename = folderFiles([path],'*.nc');
% do not include ELM-v0 for the comparison
idx_ELMv0 = find(contains(string(filename), "ELM_CRUJRA_WAD2M_SWAMPS_GCP2019_MASK11")==1);
filename(idx_ELMv0,:) = [];
for ii = 1:length(filename(:,1))
    model_name{ii} = strtrim(filename(ii,:));
end
idx_TD = contains(model_name, "SURF")|contains(model_name, "GOSAT");
TD_models = model_name(idx_TD);
idx_ML = contains(model_name, "Stanford_ML");
ML_models = model_name(idx_ML);
idx_BU = ~idx_TD&~idx_ML;
BU_models = model_name(idx_BU);
idx_weighted = contains(model_name, "BRTSim-BAMS")|contains(model_name, "ELM_v1")|contains(model_name, "Stanford");
years = 2008:2017;
for jj = 1:length(filename(:,1))
    fch4_tmp = nan(10,1);
    fch4_lat_budget_tmp = nan(1440, 720, 10);
    file = char(model_name(jj));
    fch4 = ncread([path file], 'FCH4');
    dim = size(fch4);
    model_lat = ncread([path file], 'lat');
    model_lon = ncread([path file], 'lon');
    fch4(fch4>10^20) = NaN;
    % convert FCH4 per wetland area to FCH4 per gridcell area
    % weighted by WAD2M fraction
    if (idx_weighted(jj))
        fch4 = fch4.*WAD2M;
    else
        fch4 = fch4;
    end   
    fch4_map(:,:,jj) = mean(fch4, 3, 'omitnan');
    for yr = 1:length(years)
        idx = find(idx_years>=years(yr)&idx_years<years(yr)+1);
        tmp = squeeze(fch4(:,:,idx)); % g C d-1 for all grid cells
        % sum of monthly values
        for mm = 1:12
            agg(:,:,mm) = squeeze(tmp(:,:,mm))*mon(mm).*land_area; % g C for all grid cells
        end
        tmp = sum(agg, 3, 'omitnan'); % global annual FCH4
        fch4_tmp(yr) = sum(tmp(:), 'omitnan'); % global FCH4 in g C
        fch4_lat_budget_tmp(:,:,yr) = sum(agg, 3, 'omitnan'); % global annual FCH4
    end
    % filled out missing data
    fch4_tmp(fch4_tmp==0) = NaN;
    fch4_lat_budget_tmp(fch4_lat_budget_tmp==0) = NaN;
    if (idx_TD(jj))
        fch4_TD = [fch4_TD, fch4_tmp];
        fch4_TD_idx = [fch4_TD_idx; jj];
    elseif(idx_BU(jj))
        fch4_BU = [fch4_BU, fch4_tmp];
        fch4_BU_idx = [fch4_BU_idx; jj];
    elseif(idx_ML(jj))
        fch4_ML = [fch4_ML, fch4_tmp];
        fch4_ML_idx = [fch4_ML_idx; jj];
    end
    fch4_all = [fch4_all, fch4_tmp];
    fch4_lat_budget(:,:,jj) = mean(fch4_lat_budget_tmp, 3, 'omitnan');
end
%% include site location in global maps
clc;
clf;
clear fch4_lat_mean fch4_lat_std fch4_lat_median lat lon fch4_lat_max fch4_lat_min
degree = sprintf('%c', char(176));
% read in benchmark dataset
path = '/Volumes/LaCie/research/ILAMB_CH4/FCH4_benchmark/Aug2022/ilamb_run/';
filename = 'FCH4_F_ANN_monthly_wetland_tier1.nc';
obs_fch4 = ncread([path filename], 'FCH4');
obs_fch4(obs_fch4>10^20) = nan;
idx_site = find(sum(isfinite(obs_fch4),2)>1);
obs_lat = ncread([path filename], 'lat');
obs_lon = ncread([path filename], 'lon');
obs_lat = obs_lat(idx_site);
obs_lon = obs_lon(idx_site);
plot_id = 2:2:8;
for ii = 1:4
    if (ii==1)
        idx = [fch4_BU_idx; fch4_TD_idx];
        desp = 'Ensemble-BU&TD';
    elseif (ii==2)
        idx = [fch4_BU_idx];
        desp = 'Ensemble-BU';
    elseif (ii==3)
        idx = [fch4_TD_idx];
        desp = 'Ensemble-TD';
    elseif (ii==4)
        idx = [fch4_ML_idx];
        desp = 'Machine learning';
    end
    fch4_data = fch4_map(:,:,idx);
    fch4 = mean(fch4_data, 3, 'omitnan');
    % latitudinal distribution
    fch4_regrid = fch4_lat_budget(:,:,idx);
    fch4_regrid = squeeze(sum(fch4_regrid, 'omitnan')); % latitudinal FCH4 in g C
    lat(:,ii) = model_lat;
    lon(:,ii) = model_lon;
    if (ii~=4)
        fch4_lat_mean(:,ii) = mean(fch4_regrid, 2, 'omitnan'); % mean latitudinal FCH4 across models
        fch4_lat_median(:,ii) = median(fch4_regrid, 2, 'omitnan'); % median latitudinal FCH4 across models
        fch4_lat_std(:,ii) = std(fch4_regrid', 'omitnan'); % latitudinal FCH4 in g C
        fch4_lat_max(:,ii) = max(fch4_regrid'); % latitudinal FCH4 in g C
        fch4_lat_min(:,ii) = min(fch4_regrid'); % latitudinal FCH4 in g C
    else
        fch4_lat_mean(:,ii) = fch4_regrid; % latitudinal FCH4 in g C
        fch4_lat_median(:,ii) = fch4_regrid; % latitudinal FCH4 in g C
        fch4_lat_std(:,ii) = 0; % latitudinal FCH4 in g C
        fch4_lat_max(:,ii) = fch4_regrid; % latitudinal FCH4 in g C
        fch4_lat_min(:,ii) = fch4_regrid; % latitudinal FCH4 in g C
    end
    % make plot
    subplot(4, 2, plot_id(ii))
    ax1 = subplot(4, 2, plot_id(ii));
    
    ax = worldmap('World');
    load coastlines
    zoom(1.2)
    tmp = fch4/12*16*1000; % mgCH4/d
    tmp(tmp==0) = nan;
    gradient = {'Black-red-yellow-white, custom control points', 'multigradient([0 0 0; 1 0 0; 1 1 0; 1 1 1], ''pts'', [1 5 6 7])'};
    colormap(gca, flipud(eval(gradient{2})));
    pcolorm(lat(:,ii), lon(:,ii), tmp')
    geoshow(ax, coastlat,coastlon,'Color','k')
    ax.GridLineStyle = 'none';
    framem('k-')
    mlabel off; plabel off
    gridm('off');
    caxis([0 40]);
    if (ii==4)
        for i = 1:length(obs_lon)
            % site location
            ppp = plotm(obs_lat(i), obs_lon(i),'o','MarkerEdgeColor',...
                'b', 'Markersize',3);
            hold on
        end
    end
    if (ii==1)
        c1 = colorbar();
        c1.Position = c1.Position+[.09 -0.6 -.01 0.4];
    end
    id = char(erase(id_name(ii+2),'"'));
    t = text(0.02,1.30,[id ' ' desp],'Units',...
            'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12);
end
% latitudinal distribution
% g C to Tg CH4
fch4_lat_mean = fch4_lat_mean/12*16/(10^12)*4;
fch4_lat_median = fch4_lat_median/12*16/(10^12)*4;
fch4_lat_std = fch4_lat_std/12*16/(10^12)*4;
fch4_lat_max = fch4_lat_max/12*16/(10^12)*4;
fch4_lat_min = fch4_lat_min/12*16/(10^12)*4;
subplot(4, 2, [1:2:5])
% subplot(4, 3, [1:3:10])
mcolor = [0 0 0; 0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];
ppp = [];
for ii = 2:4
    % mean with max&min       
    hhh = shadedErrorBar_flip(lat(:,ii), fch4_lat_mean(:,ii), [squeeze(fch4_lat_max(:,ii))'-fch4_lat_mean(:,ii)'; fch4_lat_mean(:,ii)'-squeeze(fch4_lat_min(:,ii))'],...
                'lineprops', {'-','color',mcolor(ii,:)}, 'transparent', 1);
    hhh.mainLine.LineWidth = 2;
    if (ii==2|ii==3)
        hhh.patch.FaceAlpha = 0.25;
    else
        hhh.patch.FaceAlpha = 0.0;
    end
    hold on
    ppp = [ppp, hhh.mainLine];
end
axis([0 12 -90 90])
a=gca;
a.XRuler.TickLabelGapOffset = -4;
t = text(0.02,0.98,['(a)'],'Units',...
            'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12);
ylabel(['Latitude in degrees (' degree ')'], 'FontSize', 12)   
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1) % sets ax1 to current axes
lgd = legend(ppp, {'Ensemble-BU',...
    'Ensemble-TD', 'Machine-learning'},...
    'orientation','vertical', 'box','off');
lgd.FontSize = 12;
lgd.Position = [0.16 0.380 0.3643 0.0452];
% global CH4 budget from individual approaches
vio_color = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];
fch4_budget = nan(max(length(fch4_BU(:)), length(fch4_TD(:))),1);
fch4_budget(1:length(fch4_BU(:)),1) = fch4_BU(:);
fch4_budget(1:length(fch4_TD(:)),2) = fch4_TD(:);
fch4_budget(1:length(fch4_ML(:)),3) = fch4_ML(:);
fch4_budget(fch4_budget==0) = nan;
% g C to Tg CH4
fch4_budget = fch4_budget/12*16/(10^12);
subplot(4, 2, 7)
median_value = median(fch4_budget,'omitnan');
mean_value = median(fch4_budget,'omitnan');
std_value = std(fch4_budget,'omitnan');
ppp = violinplot(fch4_budget);
for j = 1:length(ppp)
    ppp(j).BoxColor = 'k';
    ppp(j).BoxWidth = 0.04;
    ppp(j).ShowData = 0;
    ppp(j).ShowMean = 1;
    ppp(j).MeanPlot.LineWidth = 2;
    ppp(j).ViolinColor = vio_color(j,:);
    ppp(j).MedianPlot.SizeData = 16;
end
box on
set(gca,'xtick',[1:3])
set(gca,'xticklabel',{'BU models', 'TD models', 'ML models'},...
    'FontSize',12)
ylim([115 230])
id = char(erase(id_name(2),'"'));
t = text(0.02,0.98,[id],'Units',...
        'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12);
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1) % sets ax1 to current axes
t = text(0.12,0.315,[{'Wetland CH_4 emissions (TgCH_4 degree^{-1} yr^{-1})'}],'Units',...
            'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 10);
t = text(0.87,0.81,[{'mgCH_4 m^{-2} d^{-1}'}],'Units',...
            'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 10);
t = text(0.88,0.77,[{'(per gridcell)'}],'Units',...
            'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 10);
t = text(0.01,0.015,[{'Global wetland CH_4 emissions'}],'Units',...
            'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12, 'rotation', 90);
t = text(0.04,0.10,[{'(TgCH_4 yr^{-1})'}],'Units',...
            'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12, 'rotation', 90);
x_pos = linspace(0.14, 0.390, length(ppp));
model_year = sum(isfinite(fch4_budget));
for j = 1:length(ppp)
    t = text(x_pos(j),0.07,[ '(n = ' num2str(sprintf('%3.0f',model_year(j))) ')'],...
        'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
    t = text(x_pos(j)-0.012,0.04,[num2str(sprintf('%2.1f',mean_value(j))) char(177) ...
        num2str(sprintf('%2.1f',std_value(j)))],...
        'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
end
print('-dpng','-r300',['./GMB_comp_WAD2M_025deg_w_sites.png']);