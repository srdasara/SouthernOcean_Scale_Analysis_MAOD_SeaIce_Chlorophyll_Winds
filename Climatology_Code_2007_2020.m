


%% Initial Figures for class and for Jeff
% Contour plot of all vars. Seasonal versus Monthly might be better. 
% Time Series: all lats
% Time Series: 50 - 60; 60 - 70; 70 - 80; 80 - 90; ?

% Times are from Jan 2007 to Dec 2019



%%
cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Variables/2017_2020_vars
load('Total_timetable_SO_MOD.mat')
load('Total_timetable_SO_Depol_Ratio.mat')
load('Southern_Ocean_Wind_Speed_Data_ONEDEGREERES.mat')

load('Aqua.mat')


%%
LAT_aqua_step_edit = (linspace(Latitude_Subset_SouthernOcean(1), Latitude_Subset_SouthernOcean(end), 40))' ; 
LON_aqua_step_edit = (linspace(Longitude_Subset_SouthernOcean(1), Longitude_Subset_SouthernOcean(end), 360))'; 
LAT_aqua_step_edit = flip(LAT_aqua_step_edit); % this is for plotting purposes because pcolor flips image on y axis. 


%%
Master_Wind_Speed = flipud(Master_Wind_Speed); 

t1 = datetime(2007,01,01);
t2 = datetime(2020,12,31);
times = t1:calmonths(1):t2; 

times = times';

[winter_x, winter_y] = find(times.Month >= 6 & times.Month <= 8);
% Can check if this worked by typing 'times(winter_x)' in command window
times(winter_x)
Chl_a_winter = Master_chl_a_full_res_SouthernOcean; 
Chl_a_winter = Chl_a_winter(:,:,winter_x);

Wind_winter = Master_Wind_Speed;
Wind_winter = Wind_winter(:,:,winter_x); 

% For Spring: 
[spring_x, spring_y] = find(times.Month >= 9 & times.Month <= 11);
times(spring_x)
Chl_a_spring = Master_chl_a_full_res_SouthernOcean; 
Chl_a_spring = Chl_a_spring(:,:, spring_x); 

Wind_spring = Master_Wind_Speed;
Wind_spring = Wind_spring(:,:,spring_x);


% For Summer: 

[summer_x, summer_y] = find(times.Month >= 1 & times.Month <=2 | times.Month == 12);

times(summer_x)
Chl_a_summer = Master_chl_a_full_res_SouthernOcean; 
Chl_a_summer = Chl_a_summer(:,:, summer_x); 

Wind_summer = Master_Wind_Speed;
Wind_summer = Wind_summer(:,:,summer_x);

% For Fall: 

[fall_x, fall_y] = find(times.Month >= 3 & times.Month <= 5); 
times(fall_x) 
Chl_a_fall = Master_chl_a_full_res_SouthernOcean; 
Chl_a_fall = Chl_a_fall(:,:, fall_x); 

Wind_fall = Master_Wind_Speed;
Wind_fall = Wind_fall(:,:,fall_x);

% So now for the raw and smoothn climatologies, you can average all of
% these and plot... (for raw) or average all of them and then smoothn on
% the 2D. 

Chl_a_winter_mean = mean(Chl_a_winter,3 ,'omitnan');
Chl_a_spring_mean = mean(Chl_a_spring, 3 , 'omitnan'); 
Chl_a_summer_mean = mean(Chl_a_summer, 3, 'omitnan'); 
Chl_a_fall_mean = mean(Chl_a_fall, 3, 'omitnan'); 

Wind_winter_mean = nanmean(Wind_winter,3);
Wind_spring_mean = nanmean(Wind_spring, 3); 
Wind_summer_mean = nanmean(Wind_summer, 3); 
Wind_fall_mean = nanmean(Wind_fall, 3); 

%% adjust resolution


Aqua_Lon = LON_aqua_step_edit; 
Aqua_Lat = LAT_aqua_step_edit;

fun = @(block_struct) nanmean(block_struct.data, [ 1 2]); 

Chl_a_winter_mean_one_degree_res = blockproc(Chl_a_winter_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

Chl_a_spring_mean_one_degree_res = blockproc(Chl_a_spring_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

Chl_a_summer_mean_one_degree_res = blockproc(Chl_a_summer_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

Chl_a_fall_mean_one_degree_res = blockproc(Chl_a_fall_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 
% 


Chl_a_summer_mean_rot = rot90(Chl_a_summer_mean_one_degree_res); 
Chl_a_spring_mean_rot = rot90(Chl_a_spring_mean_one_degree_res); 
Chl_a_fall_mean_rot   = rot90(Chl_a_fall_mean_one_degree_res); 
Chl_a_winter_mean_rot = rot90(Chl_a_winter_mean_one_degree_res); 



%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.06], [0.025 0.025], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

subplot(2,2,1)
plot_SO_figure_subplot(Chl_a_winter_mean_rot, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('algae'), ...
    [0 1]); 

subplot(2,2,2)
plot_SO_figure_subplot(Chl_a_spring_mean_rot, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('algae'), ...
    [0 1]); 

subplot(2,2,3)
plot_SO_figure_subplot(Chl_a_summer_mean_rot, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('algae'), ...
    [0 1]); 

subplot(2,2,4)
plot_SO_figure_subplot(Chl_a_fall_mean_rot, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('algae'), ...
    [0 1]); 


hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.01 hp4(2)+0.2  0.015  hp4(2)+hp4(3) * 1.5]);
h.FontWeight = 'bold';
h.FontSize = 15;

%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
hy = ylabel(h, sprintf('mg m^{-3}'), 'FontSize', 18);
% hy.FontSize = 18;

%%
cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'SO_CHLA_Climatology_1_1.png','-dpng','-r300');         %  *// 300 dpi



%%

fig = figure; clf;

subplot(2,2,1)
plot_SO_figure_subplot(Wind_winter_mean, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('amp'), ...
    [0 12]); 

subplot(2,2,2)
plot_SO_figure_subplot(Wind_spring_mean, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('amp'), ...
    [0 12]); 

subplot(2,2,3)
plot_SO_figure_subplot(Wind_summer_mean, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('amp'), ...
    [0 12]); 

subplot(2,2,4)
plot_SO_figure_subplot(Wind_fall_mean, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('amp'), ...
    [0 12]); 


hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.01 hp4(2)+0.2  0.015  hp4(2)+hp4(3) * 1.5]);
h.FontWeight = 'bold';
h.FontSize = 15;

%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
hy = ylabel(h, sprintf('m s^{-1}'), 'FontSize', 18);
% hy.FontSize = 18;


%%

cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'SO_Wind_Climatology_1_1.png','-dpng','-r300');       %  *// 300 dpi





%% CMOD & Ice Now 

   % CALIPSO didn't start having data until June of 2006. This winter season will be shorter. 
   
    TR_winter_2007 = timerange('2007-06-01', '2007-09-01'); % Winter: June 1st to September 1st
    TR_winter_2008 = timerange('2008-06-01', '2008-09-01');
    TR_winter_2009 = timerange('2009-06-01', '2009-09-01');
    TR_winter_2010 = timerange('2010-06-01', '2010-09-01'); 
    TR_winter_2011 = timerange('2011-06-01', '2011-09-01');
    TR_winter_2012 = timerange('2012-06-01', '2012-09-01');
    TR_winter_2013 = timerange('2013-06-01', '2013-09-01');
    TR_winter_2014 = timerange('2014-06-01', '2014-09-01');
    TR_winter_2015 = timerange('2015-06-01', '2015-09-01');
    TR_winter_2016 = timerange('2016-06-01', '2016-09-01');
    TR_winter_2017 = timerange('2017-06-01', '2017-09-01');
    TR_winter_2018 = timerange('2018-06-01', '2018-09-01'); 
    TR_winter_2019 = timerange('2019-06-01', '2019-09-01');
    TR_winter_2020 = timerange('2020-06-01', '2020-09-01');
    
    TR_spring_2007 = timerange('2007-09-01', '2007-12-01'); % Spring: Sept 1st to Dec 1st
    TR_spring_2008 = timerange('2008-09-01', '2008-12-01'); 
    TR_spring_2009 = timerange('2009-09-01', '2009-12-01'); 
    TR_spring_2010 = timerange('2010-09-01', '2010-12-01'); 
    TR_spring_2011 = timerange('2011-09-01', '2011-12-01');
    TR_spring_2012 = timerange('2012-09-01', '2012-12-01'); 
    TR_spring_2013 = timerange('2013-09-01', '2013-12-01'); 
    TR_spring_2014 = timerange('2014-09-01', '2014-12-01');
    TR_spring_2015 = timerange('2015-09-01', '2015-12-01'); 
    TR_spring_2016 = timerange('2016-09-01', '2016-12-01'); 
    TR_spring_2017 = timerange('2017-09-01', '2017-12-01');
    TR_spring_2018 = timerange('2018-09-01', '2018-12-01'); 
    TR_spring_2019 = timerange('2019-09-01', '2019-12-01');
    
    TR_summer_2007 = timerange('2007-12-01', '2008-03-01'); % Summer: Dec 1st to March 1st
    TR_summer_2008 = timerange('2008-12-01', '2009-03-01'); 
    TR_summer_2009 = timerange('2009-12-01', '2010-03-01'); 
    TR_summer_2010 = timerange('2010-12-01', '2011-03-01'); 
    TR_summer_2011 = timerange('2011-12-01', '2012-03-01');
    TR_summer_2012 = timerange('2012-12-01', '2013-03-01'); 
    TR_summer_2013 = timerange('2013-12-01', '2014-03-01'); 
    TR_summer_2014 = timerange('2014-12-01', '2015-03-01');
    TR_summer_2015 = timerange('2015-12-01', '2016-03-01'); 
    TR_summer_2016 = timerange('2016-12-01', '2017-03-01'); 
    TR_summer_2017 = timerange('2017-12-01', '2018-03-01');
    TR_summer_2018 = timerange('2018-12-01', '2019-03-01'); 
    TR_summer_2019 = timerange('2019-12-01', '2020-03-01');
    
    %     TR_fall_2006   = timerange('2006-03-01', '2006-06-01'); % there
    %     is no CALIPSO DATA in fall of 2006
    TR_fall_2007   = timerange('2007-03-01', '2007-06-01'); %% Fall: March 1st to June 1st
    TR_fall_2008   = timerange('2008-03-01', '2008-06-01');
    TR_fall_2009   = timerange('2009-03-01', '2009-06-01');
    TR_fall_2010   = timerange('2010-03-01', '2010-06-01');
    TR_fall_2011   = timerange('2011-03-01', '2011-06-01');
    TR_fall_2012   = timerange('2012-03-01', '2012-06-01');
    TR_fall_2013   = timerange('2013-03-01', '2013-06-01');
    TR_fall_2014   = timerange('2014-03-01', '2014-06-01');
    TR_fall_2015   = timerange('2015-03-01', '2015-06-01');
    TR_fall_2016   = timerange('2016-03-01', '2016-06-01');
    TR_fall_2017   = timerange('2017-03-01', '2017-06-01');
    TR_fall_2018   = timerange('2018-03-01', '2018-06-01');
    TR_fall_2019   = timerange('2019-03-01', '2019-06-01');
    TR_fall_2020   = timerange('2020-03-01', '2020-06-01');

%%

Total_CMOD_winter = [];
Total_Lat_CMOD_winter = [];
Total_Lon_CMOD_winter = [];

Total_CMOD_spring = [];
Total_Lat_CMOD_spring = [];
Total_Lon_CMOD_spring = [];

Total_CMOD_summer = [];
Total_Lat_CMOD_summer = [];
Total_Lon_CMOD_summer = [];

Total_CMOD_fall = [];
Total_Lat_CMOD_fall = [];
Total_Lon_CMOD_fall = [];


Total_Ice_winter = [];
Total_Lat_Ice_winter = [];
Total_Lon_Ice_winter = [];

Total_Ice_spring = [];
Total_Lat_Ice_spring = [];
Total_Lon_Ice_spring = [];

Total_Ice_summer = [];
Total_Lat_Ice_summer = [];
Total_Lon_Ice_summer = [];

Total_Ice_fall = [];
Total_Lat_Ice_fall = [];
Total_Lon_Ice_fall = [];



for i = 2007:2020
    
    disp(i)
    
    SEASON = {'winter', 'spring', 'summer', 'fall'};
    
    for j = 1:length(SEASON)
        
        disp(j)
        
        % [eval(sprintf('CMOD_Bin_%d_winter', i)), eval(sprintf('night_%d_winter_occ', i)),...
        %     eval(sprintf('night_%d_winter_std')), eval(sprintf('night_%d_winter_err', i))] = ...
        
        if exist(sprintf('TR_%s_%d', SEASON{j}, i),'var') % only keep looping if the variable actually exists
            % TR_fall_2006 does not exist, which is why I need this
            % statement in here
            
            eval(sprintf('Lat_CMOD = Total_timetable_SO_MOD(TR_%s_%d,:).Total_Latitude_Surface;', SEASON{j},i))
            eval(sprintf('Lon_CMOD = Total_timetable_SO_MOD(TR_%s_%d, :).Total_Longitude_Surface;',SEASON{j}, i))
            eval(sprintf('OD       = Total_timetable_SO_MOD(TR_%s_%d, :).CMOD_Surface;', SEASON{j}, i))
            
            eval(sprintf('Lat_Ice = Total_timetable_SO_Depol_Ratio(TR_%s_%d,:).Total_Latitude_Ice;', SEASON{j},i))
            eval(sprintf('Lon_Ice = Total_timetable_SO_Depol_Ratio(TR_%s_%d,:).Total_Longitude_Ice;', SEASON{j},i))
            eval(sprintf('Ice     = Total_timetable_SO_Depol_Ratio(TR_%s_%d,:).Total_Surface_532_Integrated_Depolarization_Ratio;', SEASON{j},i)) 
            
            
            
            
            bad_Ice_values = Ice <= -0.2 | Ice > 1.2; 
            Ice(bad_Ice_values) = NaN;
            
            nan_ice        = isnan(Ice(:,1)); 
            Ice      = Ice(~nan_ice) ;
            Lat_Ice  = Lat_Ice(~nan_ice); 
            Lon_Ice  = Lon_Ice(~nan_ice); 
            
            
            
            
%             eval(sprintf('Lat_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Latitude_Wind;', SEASON{j}, i))
%             eval(sprintf('Lon_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Longitude_Wind;', SEASON{j}, i))
%             eval(sprintf('Wind     = Total_timetable_amsrmf(TR_%s_%d,:).Total_windamsrMF;', SEASON{j}, i))
         
            
            if j == 1
                
                Total_CMOD_winter     = vertcat(Total_CMOD_winter, OD);                 
                Total_Lat_CMOD_winter = vertcat(Total_Lat_CMOD_winter, Lat_CMOD);
                Total_Lon_CMOD_winter = vertcat(Total_Lon_CMOD_winter, Lon_CMOD);
                
                Total_Ice_winter     = vertcat(Total_Ice_winter, Ice);
                Total_Lat_Ice_winter = vertcat(Total_Lat_Ice_winter, Lat_Ice); 
                Total_Lon_Ice_winter = vertcat(Total_Lon_Ice_winter, Lon_Ice); 
                                
%                 Total_Wind_winter = vertcat(Total_Wind_winter, Wind); 
%                 Total_Lat_Wind_winter = vertcat(Total_Lat_Wind_winter, Lat_Wind); 
%                 Total_Lon_Wind_winter = vertcat(Total_Lon_Wind_winter, Lon_Wind); 
                
            elseif j == 2
                
                Total_CMOD_spring = vertcat(Total_CMOD_spring, OD);
                Total_Lat_CMOD_spring = vertcat(Total_Lat_CMOD_spring, Lat_CMOD);
                Total_Lon_CMOD_spring = vertcat(Total_Lon_CMOD_spring, Lon_CMOD);
                
                Total_Ice_spring = vertcat(Total_Ice_spring, Ice);
                Total_Lat_Ice_spring = vertcat(Total_Lat_Ice_spring, Lat_Ice); 
                Total_Lon_Ice_spring = vertcat(Total_Lon_Ice_spring, Lon_Ice); 
                                
%                 Total_Wind_spring = vertcat(Total_Wind_spring, Wind); 
%                 Total_Lat_Wind_spring = vertcat(Total_Lat_Wind_spring, Lat_Wind); 
%                 Total_Lon_Wind_spring = vertcat(Total_Lon_Wind_spring, Lon_Wind); 
                
            elseif j == 3
                
                Total_CMOD_summer = vertcat(Total_CMOD_summer, OD);               
                Total_Lat_CMOD_summer = vertcat(Total_Lat_CMOD_summer, Lat_CMOD);
                Total_Lon_CMOD_summer = vertcat(Total_Lon_CMOD_summer, Lon_CMOD);
                
                Total_Ice_summer = vertcat(Total_Ice_summer, Ice);
                Total_Lat_Ice_summer = vertcat(Total_Lat_Ice_summer, Lat_Ice); 
                Total_Lon_Ice_summer = vertcat(Total_Lon_Ice_summer, Lon_Ice); 
                                
%                 Total_Wind_summer = vertcat(Total_Wind_summer, Wind); 
%                 Total_Lat_Wind_summer = vertcat(Total_Lat_Wind_summer, Lat_Wind); 
%                 Total_Lon_Wind_summer = vertcat(Total_Lon_Wind_summer, Lon_Wind); 
                
            elseif j == 4
                
                Total_CMOD_fall = vertcat(Total_CMOD_fall, OD);                
                Total_Lat_CMOD_fall = vertcat(Total_Lat_CMOD_fall, Lat_CMOD);
                Total_Lon_CMOD_fall = vertcat(Total_Lon_CMOD_fall, Lon_CMOD);
                
                Total_Ice_fall = vertcat(Total_Ice_fall, Ice);
                Total_Lat_Ice_fall = vertcat(Total_Lat_Ice_fall, Lat_Ice); 
                Total_Lon_Ice_fall = vertcat(Total_Lon_Ice_fall, Lon_Ice); 
                                
%                 Total_Wind_fall = vertcat(Total_Wind_fall, Wind); 
%                 Total_Lat_Wind_fall = vertcat(Total_Lat_Wind_fall, Lat_Wind); 
%                 Total_Lon_Wind_fall = vertcat(Total_Lon_Wind_fall, Lon_Wind); 
                
            end
     
            
        end
        
        clear CMOD CMOD_OCC CMOD_STD CMOD_ERR 
        
    end
end

%%

[CMOD_one_degree_winter, CMOD_OCC_winter, CMOD_STD_winter, CMOD_ERR_winter]  = hist_wt_occ_tot(Total_Lat_CMOD_winter, Total_Lon_CMOD_winter, Total_CMOD_winter, Aqua_Lat', Aqua_Lon');
[CMOD_one_degree_spring, CMOD_OCC_spring, CMOD_STD_spring, CMOD_ERR_spring]  = hist_wt_occ_tot(Total_Lat_CMOD_spring, Total_Lon_CMOD_spring, Total_CMOD_spring, Aqua_Lat', Aqua_Lon');
[CMOD_one_degree_summer, CMOD_OCC_summer, CMOD_STD_summer, CMOD_ERR_summer]  = hist_wt_occ_tot(Total_Lat_CMOD_summer, Total_Lon_CMOD_summer, Total_CMOD_summer, Aqua_Lat', Aqua_Lon');
[CMOD_one_degree_fall, CMOD_OCC_fall, CMOD_STD_fall, CMOD_ERR_fall]          = hist_wt_occ_tot(Total_Lat_CMOD_fall, Total_Lon_CMOD_fall, Total_CMOD_fall, Aqua_Lat', Aqua_Lon');



[Ice_one_degree_winter, Ice_OCC_winter, Ice_STD_winter, Ice_ERR_winter]  = hist_wt_occ_tot(Total_Lat_Ice_winter, Total_Lon_Ice_winter, Total_Ice_winter, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_spring, Ice_OCC_spring, Ice_STD_spring, Ice_ERR_spring]  = hist_wt_occ_tot(Total_Lat_Ice_spring, Total_Lon_Ice_spring, Total_Ice_spring, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_summer, Ice_OCC_summer, Ice_STD_summer, Ice_ERR_summer]  = hist_wt_occ_tot(Total_Lat_Ice_summer, Total_Lon_Ice_summer, Total_Ice_summer, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_fall, Ice_OCC_fall, Ice_STD_fall, Ice_ERR_fall]          = hist_wt_occ_tot(Total_Lat_Ice_fall, Total_Lon_Ice_fall, Total_Ice_fall, Aqua_Lat', Aqua_Lon');


%%

cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Variables/2017_2020_vars
load('land_Mask_SO.mat')


Ice_one_degree_winter = mask3(Ice_one_degree_winter, Land_Mask_SO);
Ice_one_degree_spring = mask3(Ice_one_degree_spring, Land_Mask_SO);
Ice_one_degree_summer = mask3(Ice_one_degree_summer, Land_Mask_SO);
Ice_one_degree_fall   = mask3(Ice_one_degree_fall, Land_Mask_SO);

CMOD_one_degree_winter = mask3(CMOD_one_degree_winter, Land_Mask_SO);
CMOD_one_degree_spring = mask3(CMOD_one_degree_spring, Land_Mask_SO);
CMOD_one_degree_summer = mask3(CMOD_one_degree_summer, Land_Mask_SO);
CMOD_one_degree_fall = mask3(CMOD_one_degree_fall, Land_Mask_SO);

cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars
%%


fig = figure; clf;

subplot(2,2,1)
plot_SO_figure_subplot(Ice_one_degree_winter, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('ice'), ...
    [0.15 0.65]); 

subplot(2,2,2)
plot_SO_figure_subplot(Ice_one_degree_spring, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('ice'), ...
    [0.15 0.65]); 

subplot(2,2,3)
plot_SO_figure_subplot(Ice_one_degree_summer, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('ice'), ...
    [0.15 0.65]); 

subplot(2,2,4)
plot_SO_figure_subplot(Ice_one_degree_fall, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('ice'), ...
    [0.15 0.65]); 

hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.01 hp4(2)+0.2  0.015  hp4(2)+hp4(3) * 1.5]);
h.FontWeight = 'bold';
h.FontSize = 15;

%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
hy = ylabel(h, 'Depolarization Ratio \delta', 'FontSize', 16);
% hy.FontSize = 18;


%%
cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'SO_DepolRatio_Climatology_1_1_colorsedited.png','-dpng','-r300');       %  *// 300 dpi

%%

fig = figure; clf;

subplot(2,2,1)
plot_SO_figure_subplot(CMOD_one_degree_winter, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('tempo'), ...
    [0 0.08]); 

subplot(2,2,2)
plot_SO_figure_subplot(CMOD_one_degree_spring, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('tempo'), ...
    [0 0.08]); 

subplot(2,2,3)
plot_SO_figure_subplot(CMOD_one_degree_summer, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('tempo'), ...
    [0 0.08]); 

subplot(2,2,4)
plot_SO_figure_subplot(CMOD_one_degree_fall, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('tempo'), ...
    [0 0.08]); 

hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.01 hp4(2)+0.2  0.015  hp4(2)+hp4(3) * 1.5]);
h.FontWeight = 'bold';
h.FontSize = 15;

%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
hy = ylabel(h, 'MOD', 'FontSize', 16);
% hy.FontSize = 18;


%%
cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'SO_MOD_Climatology_1_1_landmasked.png','-dpng','-r300');       %  *// 300 dpi


%%

% Monthly averaged time series code below

Total_chl_a_monthly = zeros(1, 156);
Total_wind_monthly = zeros(1, 156);

for i = 1:156
    Total_chl_a_monthly(i) = nanmean(Master_chl_a_full_res_SouthernOcean(:,:,i), [1 2]);
    Total_wind_monthly(i)  = nanmean(Master_Wind_Speed(:,:,i), [1 2]);
end

timetable_CMOD_monthly_avg = retime(Total_timetable_SO_CMOD, 'monthly', @nanmean);
CMOD_Monthly_avg           = timetable_CMOD_monthly_avg.Total_CMOD_Night_Cloud_Free;
CMOD_Time_Months           = timetable_CMOD_monthly_avg.Total_Profile_Time_New_Night_Cloud_Free;
CMOD_Lat_Months            = timetable_CMOD_monthly_avg.Total_Latitude_Night_Cloud_Free;
CMOD_Lon_Months            = timetable_CMOD_monthly_avg.Total_Longitude_Night_Cloud_Free;


Ice_only = Total_timetable_SO_Depol_Ratio.Total_Surface_532_Integrated_Depolarization_Ratio;
Ice_only(Ice_only <= 0.15) = NaN; 

Total_timetable_SO_Depol_Ratio = addvars(Total_timetable_SO_Depol_Ratio, Ice_only); 

timetable_Depol_Ratio_monthly_avg = retime(Total_timetable_SO_Depol_Ratio, 'monthly', @nanmean);
Depol_Ratio_Monthly_avg           = timetable_Depol_Ratio_monthly_avg.Total_Surface_532_Integrated_Depolarization_Ratio;
Ice_only_Monthly_avg              = timetable_Depol_Ratio_monthly_avg.Ice_only;
Depol_Ratio_Time_Months           = timetable_Depol_Ratio_monthly_avg.Total_Profile_Time_New_Ice;
Depol_Ratio_Lat_Months            = timetable_Depol_Ratio_monthly_avg.Total_Latitude_Ice;
Depol_Ratio_Lon_Months            = timetable_Depol_Ratio_monthly_avg.Total_Longitude_Ice;




%%
 %%%%%% ----- SUBPLOTTED? ------ %%%%%%%%%
 % doc addaxis If needed
 
 
 black = rgb('black');
 grass_green = rgb('true green');
 ice_blue = rgb('dodger blue');
 wind_blue = rgb('scarlet');
 
 %%
 
 make_it_tight = true;
 subplot = @(m,n,p) subtightplot (m, n, p, [0.035 0.05], [0.1 0.025], [0.1 0.05]);
 if ~make_it_tight,  clear subplot;  end
 x = times;
 
 figure(1); clf;
 
 ax(1) = subplot(4,1,1);
 aa_splot(x, CMOD_Monthly_avg, '-',...
     'linewidth', 1.5, ...
     'Color', black);
 
 xlim([x(1) x(end)]) 
 ylim([0 0.085])
 
 ylabel('MOD')
% l(1)= legend('CMOD', 'Position',[0.76 0.92 0.1 0.0349]);

legend('MOD')
%  legend('Position', [0.105 0.94 0.16 0.0349])
 
  grid on
 set(gca, 'XTick', (times(1) : calmonths(6) : times(end)) );
%  xtickformat('MM-yyyy')
 
 
 AX=findall(0,'type','axes');
 set(AX, 'FontSize', 16)
 
 % CHLOROPHYLL
 
 ax(2) = subplot(4,1,2);
 aa_splot(x,  Total_chl_a_monthly,'-',...
     'linewidth', 1.5,...
     'MarkerSize', 4,...
     'MarkerFaceColor', grass_green,...
     'MarkerEdgeColor', grass_green,...
     'Color', grass_green);
 xlim([x(1)  x(end)]) 
 ylim([0.05  1.25])
 
 
 ylabel(sprintf('mg m^{-3}')),...
 legend('Chl-{\ita}')

 
 grid on
 set(gca, 'XTick', (times(1) : calmonths(6) : times(end)) );
%  xtickformat('MM-yyyy')
 
 
 AX=findall(0,'type','axes');
 set(AX, 'FontSize', 16)
 
 % ICE
 ax(3) = subplot(4,1,3);
 aa_splot(x,  Ice_only_Monthly_avg,'-',...
     'linewidth', 1.5,...
     'MarkerSize', 4,...
     'MarkerFaceColor', ice_blue,...
     'MarkerEdgeColor', ice_blue,...
     'Color', ice_blue);
  xlim([x(1)  x(end)]) 

 ylim([0.72  0.9])
 
 ylabel('\delta',...
     'FontName','Helvetica Neue');%         xlabel('Months')
 legend( 'Ice')
 
 
 
 grid on
 set(gca, 'XTick', (times(1) : calmonths(6) : times(end)) );
%  xtickformat('MM-yyyy')
 
 
 AX=findall(0,'type','axes');
 set(AX, 'FontSize', 16)
 
 %
 % WIND
 %
 ax(4) = subplot(4,1,4);
 aa_splot(x, Total_wind_monthly, '-',...
     'linewidth', 1.5,...
     'MarkerSize', 4,...
     'MarkerFaceColor', wind_blue,...
     'MarkerEdgeColor', wind_blue,...
     'Color', wind_blue);
 ylabel('m s^{-1}');
 legend('Wind')
 
 xlim([x(1)  x(end)]) 
 ylim([7 13])
 
 grid on
 set(gca, 'XTick', (times(1) : calmonths(6) : times(end)) );
 xtickformat('MM-yyyy')
 
 
 AX=findall(0,'type','axes');
 set(AX, 'FontSize', 22)
 
 
 xtickangle(38)
 
  [ax(1:3).XTickLabel] = deal([]);


an(1) = annotation(gcf,'textbox',... % I drew this on in the figure 
  [0.105 0.94 0.16 0.0349],...
    'String','(a)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(2) = annotation(gcf,'textbox',...
    [an(1).Position(1),...
    (ax(2).Position(2) + (an(1).Position(2) - ax(1).Position(2))),... % so this is the second axes, y position  + (y position of text box - y position of first axes)
    0.16 0.0349],...
    'String',{'(b)'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(3) = annotation(gcf, 'textbox', ...
    [an(1).Position(1),...
    (ax(3).Position(2) + (an(2).Position(2) - ax(2).Position(2))),...
    0.16 0.0349],...
    'String',{'(c)'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(4) = annotation(gcf, 'textbox', ...
    [an(1).Position(1),...
    (ax(4).Position(2) + (an(3).Position(2) - ax(3).Position(2))),...
    0.16 0.0349],...
    'String',{'(d)'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');


%%


set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'TimeSeries.png','-dpng','-r96');       %  *// 300 dpi
































