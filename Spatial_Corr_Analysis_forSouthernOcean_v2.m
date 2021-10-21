
%% 


% Southern Ocean Scale Correlation Analysis VERSION 2 with 

% MOD Surface values, redone depol ratio, data until 2020


%%Now, taking my seasonal grids of CMOD, Ice, Winds, and CHL-A, I will run a correlation through all of them and try and produce a surf plot 

 %%
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PhD_Phase_Two_SouthernOcean/Variables/2017_2020_vars

load('Total_timetable_SO_MOD_NEW.mat')
load('Total_timetable_SO_Depol_Ratio_NEW.mat')
load('Southern_Ocean_Wind_Speed_Data_ONEDEGREERES.mat')
load('Aqua.mat')

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PhD_Phase_Two_SouthernOcean/Figures
 
 
 %%

Master_Wind_Speed(Master_Wind_Speed == -500) = nan;
% for full seasons across the time range.
Master_Wind_Speed = Master_Wind_Speed(:,:,3:164);
Master_chl_a_full_res_SouthernOcean = Master_chl_a_full_res_SouthernOcean(:,:,3:164);

%%
LAT_aqua_step_edit = (linspace(Latitude_Subset_SouthernOcean(1), Latitude_Subset_SouthernOcean(end), 40))' ; 
LON_aqua_step_edit = (linspace(Longitude_Subset_SouthernOcean(1), Longitude_Subset_SouthernOcean(end), 360))'; 
LAT_aqua_step_edit = flip(LAT_aqua_step_edit); % this is for plotting purposes because pcolor flips image on y axis. 

% for i = 1:151 
% Total_chl_monthly_2(i) = nanmean(Master_chl_a(:,:,i), [1 2]); 
% end
% 
% Total_chl_a_monthly = Total_chl_monthly_2; 


%%
step_Lat = (linspace(Latitude_subset(1), Latitude_subset(end), 15)) ; 
step_Lon = (linspace(Longitude_subset(1), Longitude_subset(end), 46)); 
step_Lat = flip(step_Lat); % this is for plotting purposes because pcolor flips image on y axis. 



% At the moment this is just to check how many steps they are and if they
% correspond...
nstep_Lat = Latitude_subset(1) : -1 : Latitude_subset(end) ; 
nstep_Lon = Longitude_subset(1) : 1 : Longitude_subset(end); 
nstep_Lat = flip(nstep_Lat); 

%%




t1 = datetime(2007,03,01);
t2 = datetime(2020,08,31);
times = t1:calmonths(1):t2; 

times = times';

[winter_x, winter_y] = find(times.Month >= 6 & times.Month <= 8);
% Can check if this worked by typing 'times(winter_x)' in command window
times(winter_x) % 42 months fall into this category
Chl_a_winter = Master_chl_a_full_res_SouthernOcean; 
Chl_a_winter = Chl_a_winter(:,:,winter_x);

Wind_winter = Master_Wind_Speed;
Wind_winter = Wind_winter(:,:,winter_x); 

% For Spring: 
[spring_x, spring_y] = find(times.Month >= 9 & times.Month <= 11);
times(spring_x) % 39 months fall into this category
Chl_a_spring = Master_chl_a_full_res_SouthernOcean; 
Chl_a_spring = Chl_a_spring(:,:, spring_x); 

Wind_spring = Master_Wind_Speed;
Wind_spring = Wind_spring(:,:,spring_x);


% For Summer: 

[summer_x, summer_y] = find(times.Month >= 1 & times.Month <=2 | times.Month == 12);
times(summer_x) % 39 months fall into this category
Chl_a_summer = Master_chl_a_full_res_SouthernOcean; 
Chl_a_summer = Chl_a_summer(:,:, summer_x); 

Wind_summer = Master_Wind_Speed;
Wind_summer = Wind_summer(:,:,summer_x);

% For Fall: 

[fall_x, fall_y] = find(times.Month >= 3 & times.Month <= 5); 
times(fall_x) % 42 months fall into this category
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
Chl_a_fall_mean   = mean(Chl_a_fall, 3, 'omitnan'); 

Wind_winter_mean = mean(Wind_winter,3, 'omitnan');
Wind_spring_mean = mean(Wind_spring, 3 , 'omitnan'); 
Wind_summer_mean = mean(Wind_summer, 3, 'omitnan'); 
Wind_fall_mean   = mean(Wind_fall, 3, 'omitnan'); 

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

%% Plotting of seasonal climatologies:

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.06], [0.025 0.025], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure;clf;
subplot(2,2,1) 
plot_SO_figure_subplot(Chl_a_winter_mean_rot,...
    LON_aqua_step_edit,...
    LAT_aqua_step_edit,...
    cmocean('algae'),...
    [0 1.2],...
    'winter')
    

    hold on;

    [C, h] = m_contour(LON_aqua_step_edit, LAT_aqua_step_edit, winter_contour, 'linewi', 3, 'LineColor', [1 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);


subplot(2,2,2) 
plot_SO_figure_subplot(Chl_a_spring_mean_rot,...
    LON_aqua_step_edit,...
    LAT_aqua_step_edit,...
    cmocean('algae'),...
    [0 1.2],...
    'spring')


    hold on;

    [C, h] = m_contour(LON_aqua_step_edit, LAT_aqua_step_edit, spring_contour, 'linewi', 3, 'LineColor', [1 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

subplot(2,2,3) 
plot_SO_figure_subplot(Chl_a_summer_mean_rot,...
    LON_aqua_step_edit,...
    LAT_aqua_step_edit,...
    cmocean('algae'),...
    [0 1.2],...
    'summer')


    hold on;

    [C, h] = m_contour(LON_aqua_step_edit, LAT_aqua_step_edit, summer_contour, 'linewi', 3, 'LineColor', [1 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

subplot(2,2,4) 
plot_SO_figure_subplot(Chl_a_fall_mean_rot,...
    LON_aqua_step_edit,...
    LAT_aqua_step_edit,...
    cmocean('algae'),...
    [0 1.2],...
    'fall')


    hold on;

    [C, h] = m_contour(LON_aqua_step_edit, LAT_aqua_step_edit, fall_contour, 'linewi', 3, 'LineColor', [1 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);
    
    %%
set(gcf,'PaperPositionMode','auto')

print(fig,'CHLA_Climatology_with_iceedge_SO_scale','-dpng','-r300');       %  *// 300 dpi


%%
fig = figure;clf;
subplot(2,2,1) 
plot_SO_figure_subplot(flipud(Wind_winter_mean),...
    LON_aqua_step_edit,...
    LAT_aqua_step_edit,...
    cmocean('amp'),...
    [0 12])
    
subplot(2,2,2) 
plot_SO_figure_subplot(flipud(Wind_spring_mean),...
    LON_aqua_step_edit,...
    LAT_aqua_step_edit,...
    cmocean('amp'),...
    [0 12])

subplot(2,2,3) 
plot_SO_figure_subplot(flipud(Wind_summer_mean),...
    LON_aqua_step_edit,...
    LAT_aqua_step_edit,...
    cmocean('amp'),...
    [0 12])

subplot(2,2,4) 
plot_SO_figure_subplot(flipud(Wind_fall_mean),...
    LON_aqua_step_edit,...
    LAT_aqua_step_edit,...
    cmocean('amp'),...
    [0 12])



%% Here I am averaging data from each individual season for the temporal range of record.
clear Chl_a_seasons
count = 0;
for i = 1:3:148
    count = count+1; 
    Chl_a_seasons(:,:,count) = mean(Master_chl_a_full_res_SouthernOcean(:,:, i : i+2), 3, 'omitnan'); 
    
end

Chl_a_seasons = cat(3,Chl_a_seasons, Master_chl_a(:,:,151)); % concatenating the last month since it's the remainder 

%%

clear Chl_a_seasons
clear WINDS_SEASONS
count = 0;
for i = 1:3:162
    count = count+1; 
    disp(count)
    Chl_a_seasons(:,:,count) = mean(Master_chl_a_full_res_SouthernOcean(:,:, i : i+2), 3, 'omitnan'); 
    WINDS_SEASONS(:,:,count) = mean(Master_Wind_Speed(:,:,i : i+2), 3, 'omitnan');
end

% Chl_a_seasons = rot90(Chl_a_seasons);
% times_seasons = times(3) : calmonths(3) : times(155); 


% Chl_a_seasons = cat(3,Chl_a_seasons, Master_chl_a(:,:,151)); % concatenating the last month since it's the remainder 



%% This is a better way to adjust resolution of chl-a, that way, you can use nanmean for the missing values.

fun = @(block_struct) nanmean(block_struct.data, [ 1 2]); 

CHL_SEASONS = blockproc(Chl_a_seasons,...
    [length(Chl_a_seasons(:,1, :)) ./ length(LON_aqua_step_edit) length(Chl_a_seasons(1,:, 1)) ./ length(LAT_aqua_step_edit)],...
     fun); 


CHL_SEASONS_correct   = rot90(CHL_SEASONS);
WINDS_SEASONS_correct = flipud(WINDS_SEASONS);

clear CHL_SEASONS WINDS_SEASONS

save('CHL_SEASONS_correct.mat','CHL_SEASONS_correct', '-v7.3')
save('WINDS_SEASONS_correct.mat', 'WINDS_SEASONS_correct', '-v7.3')

%% CMOD, Winds, & Ice Now 

   % CALIPSO didn't start having data until June of 2006. This winter season will be shorter. 
   
   % 14
    TR_winter_2007 = timerange('2007-06-01', '2007-09-01'); 
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
    
    % 13
    TR_spring_2007 = timerange('2007-09-01', '2007-12-01');% Spring: Sept 1st to Dec 1st
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
    
    % 13    
    TR_summer_2007 = timerange('2007-12-01', '2008-03-01');% Summer: Dec 1st to March 1st
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
    
    % 14
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


% 54 seasons in total 


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


% 
% Total_Wind_winter = [];
% Total_Lat_Wind_winter = [];
% Total_Lon_Wind_winter = [];
% 
% Total_Wind_spring = [];
% Total_Lat_Wind_spring = [];
% Total_Lon_Wind_spring = [];
% 
% Total_Wind_summer = [];
% Total_Lat_Wind_summer = [];
% Total_Lon_Wind_summer = [];
% 
% Total_Wind_fall = [];
% Total_Lat_Wind_fall = [];
% Total_Lon_Wind_fall = [];



clear CMOD_SEASONS SeaIce_SEASONS 
for i = 2007:2020
    
    disp(i)
    
    SEASON = {'fall','winter', 'spring', 'summer'};
    
    for j = 1:length(SEASON)
        
        disp(j)
        
        % [eval(sprintf('CMOD_Bin_%d_winter', i)), eval(sprintf('night_%d_winter_occ', i)),...
        %     eval(sprintf('night_%d_winter_std')), eval(sprintf('night_%d_winter_err', i))] = ...
        
        if exist(sprintf('TR_%s_%d', SEASON{j}, i),'var') % TR_season_year in accordance with sprintf function
            
            % only keep looping if the variable actually exists
            % TR_fall_2006 does not exist, which is why I need this
            % statement in here. It will skip fall 2006, which doesn't
            % exist, and start with winter 2006, proceeding onwards until
            % end of the timeframe.. 
            
            eval(sprintf('Lat_CMOD = Total_timetable_SO_MOD(TR_%s_%d,:).Total_Latitude_Surface;', SEASON{j},i))
            eval(sprintf('Lon_CMOD = Total_timetable_SO_MOD(TR_%s_%d, :).Total_Longitude_Surface;',SEASON{j}, i))
            eval(sprintf('OD       = Total_timetable_SO_MOD(TR_%s_%d, :).CMOD_Surface;', SEASON{j}, i))
            
            eval(sprintf('Lat_Ice = Total_timetable_SO_Depol_Ratio(TR_%s_%d,:).Total_Latitude_Ice;', SEASON{j},i))
            eval(sprintf('Lon_Ice = Total_timetable_SO_Depol_Ratio(TR_%s_%d,:).Total_Longitude_Ice;', SEASON{j},i))
            eval(sprintf('Ice     = Total_timetable_SO_Depol_Ratio(TR_%s_%d,:).Total_Surface_532_Integrated_Depolarization_Ratio;', SEASON{j},i)) 
            
%             eval(sprintf('Lat_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Latitude_Wind;', SEASON{j}, i))
%             eval(sprintf('Lon_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Longitude_Wind;', SEASON{j}, i))
%             eval(sprintf('Wind     = Total_timetable_amsrmf(TR_%s_%d,:).Total_windamsrMF;', SEASON{j}, i))
                                                
            [CMOD, CMOD_OCC, CMOD_STD, CMOD_ERR]         = hist_wt_occ_tot(Lat_CMOD, Lon_CMOD, OD, LAT_aqua_step_edit', LON_aqua_step_edit');
            [SeaIce, SeaIce_OCC, SeaIce_STD, SeaIce_ERR] = hist_wt_occ_tot(Lat_Ice, Lon_Ice, Ice, LAT_aqua_step_edit', LON_aqua_step_edit');
%             [Winds, Winds_OCC, Winds_STD, Winds_ERR]     = hist_wt_occ_tot(Lat_Wind, Lon_Wind, Wind, step_Lat', step_Lon');
            
            
            
            
            
            % CMOD 3D Seasons Now
            if ~exist('CMOD_SEASONS', 'var')
                
                CMOD_SEASONS     = CMOD;
                SeaIce_SEASONS   = SeaIce;
%                 Winds_SEASONS    = Winds;
                
                CMOD_OCC_SEASONS = CMOD_OCC;
                CMOD_STD_SEASONS = CMOD_STD;
                CMOD_ERR_SEASONS = CMOD_ERR;
                
                SeaIce_OCC_SEASONS = SeaIce_OCC;
                SeaIce_STD_SEASONS = SeaIce_STD;
                SeaIce_ERR_SEASONS = SeaIce_ERR;
                
%                 Winds_OCC_SEASONS = Winds_OCC;
%                 Winds_STD_SEASONS = Winds_STD;
%                 Winds_ERR_SEASONS = Winds_ERR;
                
            else
                
                CMOD_SEASONS(:,:, end + 1) = cat(3, CMOD);
                SeaIce_SEASONS(:,:, end + 1) = cat(3, SeaIce);
%                 Winds_SEASONS(:,:, end + 1) = cat(3, Winds);
                
                CMOD_OCC_SEASONS(:,:, end + 1) = cat(3, CMOD_OCC);
                CMOD_STD_SEASONS(:,:, end + 1) = cat(3, CMOD_STD);
                CMOD_ERR_SEASONS(:,:, end + 1) = cat(3, CMOD_ERR);
                
                SeaIce_OCC_SEASONS(:,:, end + 1) = cat(3, SeaIce_OCC);
                SeaIce_STD_SEASONS(:,:, end + 1) = cat(3, SeaIce_STD);
                SeaIce_ERR_SEASONS(:,:, end + 1) = cat(3, SeaIce_ERR);
                
%                 Winds_OCC_SEASONS(:,:, end + 1) = cat(3, Winds_OCC);
%                 Winds_STD_SEASONS(:,:, end + 1) = cat(3, Winds_STD);
%                 Winds_ERR_SEASONS(:,:, end + 1) = cat(3, Winds_ERR);
                
            end                
            
        end
        
        clear CMOD CMOD_OCC CMOD_STD CMOD_ERR SeaIce SeaIce_OCC SeaIce_STD SeaIce_ERR %Winds Winds_OCC Winds_STD Winds_ERR 
        
    end
end
%%
MOD_SEASONS = CMOD_SEASONS;

%%

save('allvars_SEASONS.mat',...
    'MOD_SEASONS',...
    'SeaIce_SEASONS',...
    'CHL_SEASONS_correct',...
    'WINDS_SEASONS_correct',...
    'LAT_aqua_step_edit',...
    'LON_aqua_step_edit',...
    '-v7.3');

%%
%%%%%% This was the procedure used to get seasonal ice contour lines %%%%%%

% I basically used climatological averages as the contour lines. 

% Initializing of my variables for the loop below. 


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


%
for i = 2007:2020
    
    disp(i)
    
    SEASON = {'fall','winter', 'spring', 'summer'};
    
    for j = 1:length(SEASON)
        
        disp(j)
        
        % [eval(sprintf('CMOD_Bin_%d_winter', i)), eval(sprintf('night_%d_winter_occ', i)),...
        %     eval(sprintf('night_%d_winter_std')), eval(sprintf('night_%d_winter_err', i))] = ...
        
        if exist(sprintf('TR_%s_%d', SEASON{j}, i),'var') % only keep looping if the variable actually exists
            % TR_fall_2006 does not exist, which is why I need this
            % statement in here
            
            eval(sprintf('Lat_Ice = Total_timetable_SO_Depol_Ratio_NEW(TR_%s_%d,:).Total_Latitude_Ice;', SEASON{j},i))
            eval(sprintf('Lon_Ice = Total_timetable_SO_Depol_Ratio_NEW(TR_%s_%d,:).Total_Longitude_Ice;', SEASON{j},i))
            eval(sprintf('Ice     = Total_timetable_SO_Depol_Ratio_NEW(TR_%s_%d,:).Total_Surface_532_Integrated_Depolarization_Ratio;', SEASON{j},i)) 
            
            
            % I dont need to filter this out anymore since the
            % 'Saving_out_CMOD_Ice_Wind.m' file took care of this. 
%             
%             bad_Ice_values = Ice <= -0.2 | Ice > 1.2; 
%             Ice(bad_Ice_values) = NaN;
%             
%             nan_ice        = isnan(Ice(:,1)); 
%             Ice      = Ice(~nan_ice) ;
%             Lat_Ice  = Lat_Ice(~nan_ice); 
%             Lon_Ice  = Lon_Ice(~nan_ice); 
%             
            
            
%             
%             eval(sprintf('Lat_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Latitude_Wind;', SEASON{j}, i))
%             eval(sprintf('Lon_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Longitude_Wind;', SEASON{j}, i))
%             eval(sprintf('Wind     = Total_timetable_amsrmf(TR_%s_%d,:).Total_windamsrMF;', SEASON{j}, i))
%          
            
            if j == 1
                
                Total_Ice_fall     = vertcat(Total_Ice_fall, Ice);
                Total_Lat_Ice_fall = vertcat(Total_Lat_Ice_fall, Lat_Ice); 
                Total_Lon_Ice_fall = vertcat(Total_Lon_Ice_fall, Lon_Ice); 
                                
            
            elseif j == 2
                
                Total_Ice_winter     = vertcat(Total_Ice_winter, Ice);
                Total_Lat_Ice_winter = vertcat(Total_Lat_Ice_winter, Lat_Ice); 
                Total_Lon_Ice_winter = vertcat(Total_Lon_Ice_winter, Lon_Ice); 
                            
                       
            elseif j == 3
                
                Total_Ice_spring     = vertcat(Total_Ice_spring, Ice);
                Total_Lat_Ice_spring = vertcat(Total_Lat_Ice_spring, Lat_Ice); 
                Total_Lon_Ice_spring = vertcat(Total_Lon_Ice_spring, Lon_Ice); 
                      
                
            elseif j == 4
                
                Total_Ice_summer     = vertcat(Total_Ice_summer, Ice);
                Total_Lat_Ice_summer = vertcat(Total_Lat_Ice_summer, Lat_Ice); 
                Total_Lon_Ice_summer = vertcat(Total_Lon_Ice_summer, Lon_Ice); 
                 
            end
     
            
        end
        
        clear CMOD CMOD_OCC CMOD_STD CMOD_ERR 
        
    end
end

%%


[Ice_one_degree_winter, Ice_OCC_winter, Ice_STD_winter, Ice_ERR_winter]  = hist_wt_occ_tot(Total_Lat_Ice_winter, Total_Lon_Ice_winter, Total_Ice_winter, LAT_aqua_step_edit', LON_aqua_step_edit');
[Ice_one_degree_spring, Ice_OCC_spring, Ice_STD_spring, Ice_ERR_spring]  = hist_wt_occ_tot(Total_Lat_Ice_spring, Total_Lon_Ice_spring, Total_Ice_spring, LAT_aqua_step_edit', LON_aqua_step_edit');
[Ice_one_degree_summer, Ice_OCC_summer, Ice_STD_summer, Ice_ERR_summer]  = hist_wt_occ_tot(Total_Lat_Ice_summer, Total_Lon_Ice_summer, Total_Ice_summer, LAT_aqua_step_edit', LON_aqua_step_edit');
[Ice_one_degree_fall, Ice_OCC_fall, Ice_STD_fall, Ice_ERR_fall]          = hist_wt_occ_tot(Total_Lat_Ice_fall, Total_Lon_Ice_fall, Total_Ice_fall, LAT_aqua_step_edit', LON_aqua_step_edit');

%

winter_contour = Ice_one_degree_winter; 
spring_contour = Ice_one_degree_spring; 
summer_contour = Ice_one_degree_summer;
fall_contour   = Ice_one_degree_fall;


%%

save('ice_contour_seasons.mat', ...
    'winter_contour',...
    'spring_contour',...
    'summer_contour',...
    'fall_contour',...
    '-v7.3'); 

%%
%%%%% Keeping track of my seasons here %%%%%%%%%%%%


t1 = datetime(2007,03,01);
t2 = datetime(2020,08,31);
times = t1:calmonths(1):t2; 

clear t1 t2

times_seasons = times(1) : calmonths(3) : times(end); 



%%

% You are not using this function in the updated version of this
% manuscript, so commenting it out: 
% 
% CMOD_SEASONS_SMOOTHN = zeros(size(CMOD_SEASONS)); 

for i = 1:50
    CMOD_SEASONS_SMOOTHN(:,:,i) = smoothn(CMOD_SEASONS(:,:,i), 'robust') ;
end


CMOD_SEASONS_SMOOTHN(:,:,39) = CMOD_SEASONS(:,:,39); 
CMOD_SEASONS_SMOOTHN(:,:,11) = CMOD_SEASONS(:,:,11); 
CMOD_SEASONS_SMOOTHN(:,:,31) = CMOD_SEASONS(:,:,31); 
CMOD_SEASONS_SMOOTHN(:,:,35) = CMOD_SEASONS(:,:,35); 
CMOD_SEASONS_SMOOTHN(:,:,19) = CMOD_SEASONS(:,:,19); 
CMOD_SEASONS_SMOOTHN(:,:,47) = CMOD_SEASONS(:,:,47);


%%

% I need to find the mean of each season and then plot out how it's changing from that mean across time. 


%% SeaIce Test with masking values below 0.15

SeaIce_SEASONS_masked = SeaIce_SEASONS; 

SeaIce_SEASONS_masked(SeaIce_SEASONS_masked < 0.15) = NaN;



 
%%

%%%%%%%%% Here I start calculating the spatial correlation between seasons %%%%%%%%% 

% Matrix preallocation 
CMOD_CHL_RHO_MATRIX  = zeros(40,360,4) ; 
CMOD_CHL_PVAL_MATRIX = zeros(40,360,4) ; 

CMOD_ICE_RHO_MATRIX  = zeros(40,360,4) ; 
CMOD_ICE_PVAL_MATRIX = zeros(40,360,4) ; 

CMOD_WIND_RHO_MATRIX  = zeros(40,360,4) ; 
CMOD_WIND_PVAL_MATRIX = zeros(40,360,4) ;  

CHL_ICE_RHO_MATRIX = zeros(40,360,4) ; 
CHL_ICE_PVAL_MATRIX = zeros(40,360,4) ;  

CHL_WIND_RHO_MATRIX = zeros(40,360,4) ; 
CHL_WIND_PVAL_MATRIX = zeros(40,360,4) ; 

WIND_ICE_RHO_MATRIX = zeros(40,360,4) ; 
WIND_ICE_PVAL_MATRIX = zeros(40,360,4) ; 
    
%%
clear x y z k 

for i = 1 : length(LAT_aqua_step_edit(:,1))
    % disp(i)
    
    for j = 1 : length(LON_aqua_step_edit(:,1))
        %  disp(j)
        
        for s = 1:4
          
%             disp(s)
            
            x = squeeze(MOD_SEASONS(i,j,s:4:51));
            y = squeeze(CHL_SEASONS_correct(i,j,s:4:51));
            z = squeeze(SeaIce_SEASONS_masked(i,j,s:4:51));
            k = squeeze(WINDS_SEASONS_correct(i,j,s:4:51));
            
            % CMOD & CHL
            if sum(isnan(x)) || sum(isnan(y)) <= 3
                [CMOD_CHL_RHO_MATRIX(i, j,s), CMOD_CHL_PVAL_MATRIX(i, j, s)]  = corr(x, y,'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(x)) || sum(isnan(y)) > 3
                isnan([CMOD_CHL_RHO_MATRIX(i, j,s)  , CMOD_CHL_PVAL_MATRIX(i, j, s)]);
            end
            
            
            % CMOD & Ice
            if sum(isnan(x)) || sum(isnan(z)) <= 3
                [CMOD_ICE_RHO_MATRIX(i, j,s), CMOD_ICE_PVAL_MATRIX(i, j,s)]   = corr(x, z, 'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(x)) || sum(isnan(z)) > 3
                isnan([CMOD_ICE_RHO_MATRIX(i, j,s)  , CMOD_ICE_PVAL_MATRIX(i, j,s)]); 
            end
             
            
            % CMOD & Wind
            if sum(isnan(x)) || sum(isnan(k)) <= 3
                [CMOD_WIND_RHO_MATRIX(i, j,s), CMOD_WIND_PVAL_MATRIX(i, j,s)] = corr(x, k, 'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(x)) || sum(isnan(k)) > 3
                isnan([CMOD_WIND_RHO_MATRIX(i, j,s) , CMOD_WIND_PVAL_MATRIX(i, j,s)]);
            end
            
            
            % CHL & Ice
            if sum(isnan(z)) || sum(isnan(y)) <= 3
                [CHL_ICE_RHO_MATRIX(i, j,s), CHL_ICE_PVAL_MATRIX(i, j,s)]     = corr(z, y, 'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(z)) || sum(isnan(y)) > 3
                isnan([CHL_ICE_RHO_MATRIX(i, j,s)   , CHL_ICE_PVAL_MATRIX(i, j,s)]);
            end
            
            
            % CHL & Wind
            if sum(isnan(k)) || sum(isnan(y)) <= 3
                [CHL_WIND_RHO_MATRIX(i,j,s), CHL_WIND_PVAL_MATRIX(i,j,s)]     = corr(k, y, 'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(k)) || sum(isnan(y)) > 3
                isnan([CHL_WIND_RHO_MATRIX(i,j,s)   , CHL_WIND_PVAL_MATRIX(i,j,s)]);
            end
            
            
            % Wind & Ice
            if sum(isnan(k)) || sum(isnan(z)) <= 3
                [WIND_ICE_RHO_MATRIX(i,j,s), WIND_ICE_PVAL_MATRIX(i,j,s)]     = corr(k, z, 'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(k)) || sum(isnan(z)) > 3
                isnan([WIND_ICE_RHO_MATRIX(i,j,s)   , WIND_ICE_PVAL_MATRIX(i,j,s)]);
            end
%             
        end
        
    end
    
end

%%

% %%
% 
% clear x y z k 
% 
% for i = 1 : length(nstep_Lat(1,:))
%     % disp(i)
%     
%     for j = 1 : length(nstep_Lon(1,:))
%         %  disp(j)
%         
%         for s = 1:4
%             
%             x = (squeeze(CMOD_SEASONS(i,j,s:4:end)));
%             y = squeeze(CHL_SEASONS(i,j,s:4:end)); 
%             z = squeeze(SeaIce_SEASONS_masked(i,j,s:4:end)); 
%             k = squeeze(WINDS_SEASONS(i,j,s:4:end)); 
%             
%             if sum(isnan(x)) || sum(isnan(y)) || sum(isnan(z)) || sum(isnan(k)) <= 5 
%             
%                 
%             
%             [CMOD_CHL_RHO_MATRIX(i, j,s), CMOD_CHL_PVAL_MATRIX(i, j, s)]  = corr(x, y,'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CMOD_ICE_RHO_MATRIX(i, j,s), CMOD_ICE_PVAL_MATRIX(i, j,s)]   = corr(x, z, 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CMOD_WIND_RHO_MATRIX(i, j,s), CMOD_WIND_PVAL_MATRIX(i, j,s)] = corr(x, k, 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CHL_ICE_RHO_MATRIX(i, j,s), CHL_ICE_PVAL_MATRIX(i, j,s)]     = corr(z, y, 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CHL_WIND_RHO_MATRIX(i,j,s), CHL_WIND_PVAL_MATRIX(i,j,s)]     = corr(k, y, 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [WIND_ICE_RHO_MATRIX(i,j,s), WIND_ICE_PVAL_MATRIX(i,j,s)]     = corr(k, z, 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             elseif sum(isnan(x)) || sum(isnan(y)) || sum(isnan(z)) || sum(isnan(k)) > 5 
%                 
%                 [CMOD_CHL_RHO_MATRIX(i, j,s)  , CMOD_CHL_PVAL_MATRIX(i, j, s)]   = NaN;
%                 [CMOD_ICE_RHO_MATRIX(i, j,s)  , CMOD_ICE_PVAL_MATRIX(i, j,s)]    = NaN;
%                 [CMOD_WIND_RHO_MATRIX(i, j,s) , CMOD_WIND_PVAL_MATRIX(i, j,s)]   = NaN;
%                 [CHL_ICE_RHO_MATRIX(i, j,s)   , CHL_ICE_PVAL_MATRIX(i, j,s)]     = NaN;
%                 [CHL_WIND_RHO_MATRIX(i,j,s)   , CHL_WIND_PVAL_MATRIX(i,j,s)]     = NaN;
%                 [WIND_ICE_RHO_MATRIX(i,j,s)   , WIND_ICE_PVAL_MATRIX(i,j,s)]     = NaN;
%                 
%             end
%             
%         end
%         
%     end
%     
% end
% 
% %%
%    
% for i = 1 : length(nstep_Lat(1,:))
%     % disp(i)
%     
%     for j = 1 : length(nstep_Lon(1,:))
%         %  disp(j)
%         
%         for s = 1:4
%             
%             
%             [CMOD_CHL_RHO_MATRIX(i, j,s), CMOD_CHL_PVAL_MATRIX(i, j, s)] = corr(squeeze(CMOD_SEASONS(i,j,s:4:end)),...
%                 squeeze(CHL_SEASONS(i,j,s:4:end)),'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CMOD_ICE_RHO_MATRIX(i, j,s), CMOD_ICE_PVAL_MATRIX(i, j,s)] = corr(squeeze(CMOD_SEASONS(i,j,s:4:end)),...
%                 squeeze(SeaIce_SEASONS_masked(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CMOD_WIND_RHO_MATRIX(i, j,s), CMOD_WIND_PVAL_MATRIX(i, j,s)] = corr(squeeze(CMOD_SEASONS(i,j,s:4:end)), ...
%                 squeeze(Winds_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CHL_ICE_RHO_MATRIX(i, j,s), CHL_ICE_PVAL_MATRIX(i, j,s)] = corr(squeeze(SeaIce_SEASONS_masked(i,j,s:4:end)),...
%                 squeeze(CHL_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CHL_WIND_RHO_MATRIX(i,j,s), CHL_WIND_PVAL_MATRIX(i,j,s)] = corr(squeeze(Winds_SEASONS(i,j,s:4:end)),...
%                 squeeze(CHL_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [WIND_ICE_RHO_MATRIX(i,j,s), WIND_ICE_PVAL_MATRIX(i,j,s)] = corr(squeeze(Winds_SEASONS(i,j,s:4:end)),...
%                 squeeze(SeaIce_SEASONS_masked(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             
%         end
%         
%     end
%     
% end
% 
% 
% 
% %%
% 
% for i = 1 : length(nstep_Lat(1,:))
%     % disp(i)
%     
%     for j = 1 : length(nstep_Lon(1,:))
%         %  disp(j)
%         
%         for s = 1:4
%             
%             [CMOD_CHL_RHO_MATRIX(i, j,s), CMOD_CHL_PVAL_MATRIX(i, j, s)] = corr(squeeze(CMOD_SEASONS_SMOOTHN(i,j,s:4:end)),...
%                 squeeze(CHL_SEASONS(i,j,s:4:end)),'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CMOD_ICE_RHO_MATRIX(i, j,s), CMOD_ICE_PVAL_MATRIX(i, j,s)] = corr(squeeze(CMOD_SEASONS_SMOOTHN(i,j,s:4:end)),...
%                 squeeze(SeaIce_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CMOD_WIND_RHO_MATRIX(i, j,s), CMOD_WIND_PVAL_MATRIX(i, j,s)] = corr(squeeze(CMOD_SEASONS_SMOOTHN(i,j,s:4:end)), ...
%                 squeeze(Winds_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CHL_ICE_RHO_MATRIX(i, j,s), CHL_ICE_PVAL_MATRIX(i, j,s)] = corr(squeeze(SeaIce_SEASONS(i,j,s:4:end)),...
%                 squeeze(CHL_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [CHL_WIND_RHO_MATRIX(i,j,s), CHL_WIND_PVAL_MATRIX(i,j,s)] = corr(squeeze(Winds_SEASONS(i,j,s:4:end)),...
%                 squeeze(CHL_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
%             
%             [WIND_ICE_RHO_MATRIX(i,j,s), WIND_ICE_PVAL_MATRIX(i,j,s)] = corr(squeeze(Winds_SEASONS(i,j,s:4:end)),...
%                 squeeze(SeaIce_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
%             
%         end
%         
%     end
%     
% end
%%


Lon = LON_aqua_step_edit;
Lat = LAT_aqua_step_edit;

%%%%%%%%% PLOTTING OF SPATIAL CORRELATION FIGURES BELOW %%%%%%%%%%%%%%%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.06], [0.025 0.025], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 

subplot(2,2,4) % fall
plot_SO_figure_subplot(SeaIce_SEASONS(:,:,1), ...
    Lon,... 
    Lat,...
    cmocean('ice'),...
    [0 .8]);
    hold on;


    [C, h] = m_contour(Lon, Lat, fall_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65 ; 
    clabel(C, h,'fontsize',13);
%     clegendm(C,h,'m')

subplot(2,2,1) % winter
plot_SO_figure_subplot(SeaIce_SEASONS(:,:,2), ...
    Lon,...
    Lat,...
    cmocean('ice'),...
    [0 .8]);

    hold on;

    [C, h] = m_contour(Lon, Lat, winter_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

subplot(2,2,2) % spring
plot_SO_figure_subplot(SeaIce_SEASONS(:,:,3), ...
    Lon,...
    Lat,...
     cmocean('ice'),...
    [0 .8]);
    hold on;

    [C, h] = m_contour(Lon, Lat, spring_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

subplot(2,2,3) % summer
plot_SO_figure_subplot(SeaIce_SEASONS(:,:,4), ...
    Lon,...
    Lat,...
    cmocean('ice'),...
    [0 .8]);

    hold on;
    
    [C, h] = m_contour(Lon, Lat, summer_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

 
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.01 hp4(2)+0.2  0.015  hp4(2)+hp4(3) * 1.5]);
h.FontWeight = 'bold';
h.FontSize = 15;

% h = colorbar('Location', 'westoutside');
% % ylabel(h, sprintf(ylabel_text))
% %
% h.FontSize = 18;
% h.Position = [0.1 0.1782 0.0161 0.6293] ; 
%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
 hy = ylabel(h, 'Ice test', 'FontSize', 18);

%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.02], [0.05 0.05], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 

subplot(2,2,4) % fall
plot_SO_figure_subplot(CMOD_ICE_RHO_MATRIX(:,:,1), ... 
    Lon,...
    Lat,...
    redblue,...
    [-1 1])
    hold on;


    [C, h] = m_contour(Lon, Lat, fall_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65 ; 
    clabel(C, h,'fontsize',13);
%     clegendm(C,h,'m')

subplot(2,2,1) % winter
plot_SO_figure_subplot(CMOD_ICE_RHO_MATRIX(:,:,2), ... 
    Lon,...
    Lat,...
    redblue,...
    [-1 1]); 

    hold on;

    [C, h] = m_contour(Lon, Lat, winter_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

subplot(2,2,2) % spring
plot_SO_figure_subplot(CMOD_ICE_RHO_MATRIX(:,:,3), ... 
    Lon,...
    Lat,...
    redblue,...    
    [-1 1]); 
    hold on;

    [C, h] = m_contour(Lon, Lat, spring_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

subplot(2,2,3) %summer
plot_SO_figure_subplot(CMOD_ICE_RHO_MATRIX(:,:,4), ... 
    Lon,...
    Lat,...
    redblue,...
    [-1 1]);
    hold on;
    
    [C, h] = m_contour(Lon, Lat, summer_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

 
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.01 hp4(2)+0.125  0.015  hp4(2)+hp4(3) * 1.5]);
h.FontWeight = 'bold';
h.FontSize = 15;

% h = colorbar('Location', 'westoutside');
% % ylabel(h, sprintf(ylabel_text))
% %
% h.FontSize = 18;
% h.Position = [0.1 0.1782 0.0161 0.6293] ; 
%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
 hy = ylabel(h, 'MOD & Ice seasonal correlation coefficient', 'FontSize', 18);

%%


cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'Updated_MOD_Ice_Correlation_1_1.png','-dpng','-r300');       %  *// 300 dpi
%  print(gcf, 'Ju

%%
% %%
% 
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
% if ~make_it_tight,  clear subplot;  end
% 
% fig = figure;clf;
% 
% % chl_season = rot90(Chl_a_seasons_one_degree_res); 
% 
% % test_2 = smoothn(test, 'robust'); 
% subplot(2,2,1)
% plot_BSea_figure_subplot(CMOD_ICE_PVAL_MATRIX(:,:,1), ...
%     step_Lon,...
%     step_Lat,...
%     [0 0.2],...
%     'Winter'); 
% 
%     hold on;
% 
% 
%     [C, h] = m_contour(step_Lon, step_Lat, winter_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]);  
%     h.LevelList = 0.65 ; 
% %     clabel(C, h,'fontsize',13);
% %     clegendm(C,h,'m')
% 
% subplot(2,2,2)
% plot_BSea_figure_subplot(CMOD_ICE_PVAL_MATRIX(:,:,2), ...
%     step_Lon,...
%     step_Lat,...
%     [0 0.2],...
%     'Spring'); 
% 
%     hold on;
% 
%     [C, h] = m_contour(step_Lon, step_Lat, spring_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
%     h.LevelList = 0.65  ; 
% %     clabel(C, h,'fontsize',13);
% 
% subplot(2,2,3)
% plot_BSea_figure_subplot(CMOD_ICE_PVAL_MATRIX(:,:,3), ...
%     step_Lon,...
%     step_Lat,...
%     [0 0.2],...
%     'Summer'); 
%     hold on;
% 
%     [C, h] = m_contour(step_Lon, step_Lat, summer_contour, 'linewi', 5.5, 'LineColor',[ 0 0 0]); 
%     h.LevelList = 0.65  ; 
% %     clabel(C, h,'fontsize',13);
% 
% subplot(2,2,4)
% plot_BSea_figure_subplot(CMOD_ICE_PVAL_MATRIX(:,:,4), ...
%     step_Lon,...
%     step_Lat,...
%     [0 0.2],...
%     'Fall'); 
% 
%     hold on;
%     [C, h] = m_contour(step_Lon, step_Lat, fall_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
%     h.LevelList = 0.65  ; 
% %     clabel(C, h,'fontsize',13);
% 
%  
%     
% hp4 = get(subplot(2,2,4),'Position');
% h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
% h.FontWeight = 'bold';
% h.FontSize = 15;
% %  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
%  hy = ylabel(h, 'CMOD & Ice seasonal p-values', 'FontSize', 18);

%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.02], [0.05 0.05], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 

subplot(3,1,3) % fall
plot_SO_figure_subplot(CHL_ICE_RHO_MATRIX(:,:,1), ... 
    Lon,...
    Lat,...
    redblue,...
    [-1 1])
    hold on;


    [C, h] = m_contour(Lon, Lat, fall_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65 ; 
    clabel(C, h,'fontsize',13);
%     clegendm(C,h,'m')

% subplot(2,2,1) % winter
% plot_SO_figure_subplot(CHL_ICE_RHO_MATRIX(:,:,2), ... 
%     Lon,...
%     Lat,...
%     redblue,...
%     [-1 1]); 
% 
%     hold on;
% 
%     [C, h] = m_contour(Lon, Lat, winter_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
%     h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(3,1,1) % spring
plot_SO_figure_subplot(CHL_ICE_RHO_MATRIX(:,:,3), ... 
    Lon,...
    Lat,...
    redblue,...    
    [-1 1]); 
    hold on;

    [C, h] = m_contour(Lon, Lat, spring_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

subplot(3,1,2) %summer
plot_SO_figure_subplot(CHL_ICE_RHO_MATRIX(:,:,4), ... 
    Lon,...
    Lat,...
    redblue,...
    [-1 1]);
    hold on;
    
    [C, h] = m_contour(Lon, Lat, summer_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

 
hp4 = get(subplot(3,1,3),'Position');
h = colorbar('Position', [0.65 0.2  0.02  0.6]);
h.FontWeight = 'bold';
h.FontSize = 15;

% h = colorbar('Location', 'westoutside');
% % ylabel(h, sprintf(ylabel_text))
% %
% h.FontSize = 18;
% h.Position = [0.1 0.1782 0.0161 0.6293] ; 
%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
 hy = ylabel(h, 'Ice & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);

%%


cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'Updated_CHL_Ice_Correlation_1_1_3figsubplot.png','-dpng','-r300');       %  *// 300 dpi
%  print(gcf, 'Ju

%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.02], [0.05 0.05], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 

subplot(2,2,4) % fall
plot_SO_figure_subplot(CMOD_CHL_RHO_MATRIX(:,:,1), ... 
    Lon,...
    Lat,...
    redblue,...
    [-1 1])
    hold on;


    [C, h] = m_contour(Lon, Lat, fall_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65 ; 
    clabel(C, h,'fontsize',13);
%     clegendm(C,h,'m')

subplot(2,2,1) % winter
plot_SO_figure_subplot(CMOD_CHL_RHO_MATRIX(:,:,2), ... 
    Lon,...
    Lat,...
    redblue,...
    [-1 1]); 

    hold on;

    [C, h] = m_contour(Lon, Lat, winter_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

subplot(2,2,2) % spring
plot_SO_figure_subplot(CMOD_CHL_RHO_MATRIX(:,:,3), ... 
    Lon,...
    Lat,...
    redblue,...    
    [-1 1]); 
    hold on;

    [C, h] = m_contour(Lon, Lat, spring_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

subplot(2,2,3) %summer
plot_SO_figure_subplot(CMOD_CHL_RHO_MATRIX(:,:,4), ... 
    Lon,...
    Lat,...
    redblue,...
    [-1 1]);
    hold on;
    
    [C, h] = m_contour(Lon, Lat, summer_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

 
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.01 hp4(2)+0.125  0.015  hp4(2)+hp4(3) * 1.5]);
h.FontWeight = 'bold';
h.FontSize = 15;

% h = colorbar('Location', 'westoutside');
% % ylabel(h, sprintf(ylabel_text))
% %
% h.FontSize = 18;
% h.Position = [0.1 0.1782 0.0161 0.6293] ; 
%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
 hy = ylabel(h, 'MOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);

 %%
cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'Updated_MOD_CHL_Correlation_1_1.png','-dpng','-r300');       %  *// 300 dpi
%  print(gcf, 'Ju


%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.02], [0.05 0.05], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 

subplot(2,2,4) % fall
plot_SO_figure_subplot(CMOD_WIND_RHO_MATRIX(:,:,1), ... 
    Lon,...
    Lat,...
    redblue,...
    [-1 1])
    hold on;


    [C, h] = m_contour(Lon, Lat, fall_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65 ; 
    clabel(C, h,'fontsize',13);
%     clegendm(C,h,'m')

subplot(2,2,1) % winter
plot_SO_figure_subplot(CMOD_WIND_RHO_MATRIX(:,:,2), ... 
    Lon,...
    Lat,...
    redblue,...
    [-1 1]); 

    hold on;

    [C, h] = m_contour(Lon, Lat, winter_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

subplot(2,2,2) % spring
plot_SO_figure_subplot(CMOD_WIND_RHO_MATRIX(:,:,3), ... 
    Lon,...
    Lat,...
    redblue,...    
    [-1 1]); 
    hold on;

    [C, h] = m_contour(Lon, Lat, spring_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

subplot(2,2,3) %summer
plot_SO_figure_subplot(CMOD_WIND_RHO_MATRIX(:,:,4), ... 
    Lon,...
    Lat,...
    redblue,...
    [-1 1]);
    hold on;
    
    [C, h] = m_contour(Lon, Lat, summer_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
    clabel(C, h,'fontsize',13);

 
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.01 hp4(2)+0.125  0.015  hp4(2)+hp4(3) * 1.5]);
h.FontWeight = 'bold';
h.FontSize = 15;

% h = colorbar('Location', 'westoutside');
% % ylabel(h, sprintf(ylabel_text))
% %
% h.FontSize = 18;
% h.Position = [0.1 0.1782 0.0161 0.6293] ; 
%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
 hy = ylabel(h, 'MOD & Wind Speed seasonal correlation coefficient', 'FontSize', 18);

 %%
cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'Updated_MOD_WIND_Correlation_1_1.png','-dpng','-r300');       %  *// 300 dpi
%  print(gcf, 'Ju

%% STATISTICAL CALCULATIONS %%


% Matrix preallocation 
randn_CMOD_CHL_RHO_MATRIX  = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 
randn_CMOD_CHL_PVAL_MATRIX = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 

randn_CMOD_ICE_RHO_MATRIX  = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 
randn_CMOD_ICE_PVAL_MATRIX = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 

randn_CMOD_WIND_RHO_MATRIX  = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 
randn_CMOD_WIND_PVAL_MATRIX = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 

randn_CHL_ICE_RHO_MATRIX = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 
randn_CHL_ICE_PVAL_MATRIX = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 

randn_CHL_WIND_RHO_MATRIX = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 
randn_CHL_WIND_PVAL_MATRIX = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 

randn_WIND_ICE_RHO_MATRIX = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 
randn_WIND_ICE_PVAL_MATRIX = zeros([length(MOD_SEASONS(:, 1, 1)), length(MOD_SEASONS(1, :, 1))]) ; 
    

%%

CMOD_randn_SEASONS = randn(size(MOD_SEASONS)); 
CHL_randn_SEASONS = randn(size(CHL_SEASONS_correct)); 
SeaIce_randn_SEASONS = randn(size(SeaIce_SEASONS_masked)); 
Winds_randn_SEASONS = randn(size(WINDS_SEASONS_correct)); 

% might not need the below:
CMOD_randn_SEASONS   = spatialPattern(size(CMOD_SEASONS),-1);
CHL_randn_SEASONS    = spatialPattern(size(CMOD_SEASONS),-1);
SeaIce_randn_SEASONS = spatialPattern(size(CMOD_SEASONS),-1);
Winds_randn_SEASONS  = spatialPattern(size(CMOD_SEASONS),-1);

%%
   
for i = 1 : length(LAT_aqua_step_edit(:,1))
    % disp(i)
    
    for j = 1 : length(LON_aqua_step_edit(:,1))
        %  disp(j)
        
        for s = 1:4
            
            [randn_CMOD_CHL_RHO_MATRIX(i, j,s), randn_CMOD_CHL_PVAL_MATRIX(i, j, s)] = corr(squeeze(CMOD_randn_SEASONS(i,j,s:4:end)),...
                squeeze(CHL_randn_SEASONS(i,j,s:4:54)),'Type', 'Pearson', 'rows', 'pairwise');
            
            [randn_CMOD_ICE_RHO_MATRIX(i, j,s), randn_CMOD_ICE_PVAL_MATRIX(i, j,s)] = corr(squeeze(CMOD_randn_SEASONS(i,j,s:4:end)),...
                squeeze(SeaIce_randn_SEASONS(i,j,s:4:54)), 'Type', 'Pearson', 'rows', 'pairwise');
            
            [randn_CMOD_WIND_RHO_MATRIX(i, j,s), randn_CMOD_WIND_PVAL_MATRIX(i, j,s)] = corr(squeeze(CMOD_randn_SEASONS(i,j,s:4:end)), ...
                squeeze(Winds_randn_SEASONS(i,j,s:4:54)), 'Type', 'Pearson', 'rows', 'pairwise');
            
            [randn_CHL_ICE_RHO_MATRIX(i, j,s), randn_CHL_ICE_PVAL_MATRIX(i, j,s)] = corr(squeeze(SeaIce_randn_SEASONS(i,j,s:4:end)),...
                squeeze(CHL_randn_SEASONS(i,j,s:4:54)), 'Type', 'Pearson', 'rows', 'pairwise');
            
            [randn_CHL_WIND_RHO_MATRIX(i,j,s), randn_CHL_WIND_PVAL_MATRIX(i,j,s)] = corr(squeeze(Winds_randn_SEASONS(i,j,s:4:end)),...
                squeeze(CHL_randn_SEASONS(i,j,s:4:54)), 'Type', 'Pearson', 'rows', 'pairwise');
            
            [randn_WIND_ICE_RHO_MATRIX(i,j,s), randn_WIND_ICE_PVAL_MATRIX(i,j,s)] = corr(squeeze(Winds_randn_SEASONS(i,j,s:4:end)),...
                squeeze(SeaIce_randn_SEASONS(i,j,s:4:54)), 'Type', 'Pearson', 'rows', 'pairwise');
            
        end
        
    end
    
end


%%

% 
% data_permuted = permute(CMOD_ICE_RHO_MATRIX,[3 1 2]);
% rand_permuted = permute(randn_CMOD_CHL_RHO_MATRIX, [ 3 1 2]);
% 
% % [ h, p ] = ttest(CMOD_ICE_RHO_MATRIX, randn_CMOD_CHL_RHO_MATRIX);
% [h,p] = ttest(data_permuted,rand_permuted);
% 
% testh = squeeze(h);
% testp = squeeze(p);
% % [ h, p ] = ttest2(CMOD_ICE_RHO_MATRIX, randn_CMOD_CHL_RHO_MATRIX);

% missing values in same place as for CMOD & Ice, Ice & chl-a rho matrices

% data_array = CMOD_ICE_RHO_MATRIX;
% random_array = randn_CMOD_ICE_RHO_MATRIX; 
% 
% save('data_array.mat', 'data_array'); 
% save('random_array.mat', 'random_array'); 
% 
% [row, col] = find(isnan(CMOD_ICE_RHO_MATRIX));




%%
% rand matrix creation with NaNs in same places as NaNs in data matrices

% CMOD & ICe
test_CMOD_ICE_RHO_MATRIX = CMOD_ICE_RHO_MATRIX; 
test_CMOD_ICE_RHO_MATRIX(test_CMOD_ICE_RHO_MATRIX==0) = NaN;

NAN_randn_CMOD_ICE_RHO_MATRIX = nan(size(CMOD_ICE_RHO_MATRIX));
NAN_randn_CMOD_ICE_RHO_MATRIX(~isnan(test_CMOD_ICE_RHO_MATRIX)) = randn_CMOD_ICE_RHO_MATRIX(~isnan(test_CMOD_ICE_RHO_MATRIX)); 


% CHL & ICE
test_CHL_ICE_RHO_MATRIX = CHL_ICE_RHO_MATRIX; 
test_CHL_ICE_RHO_MATRIX(test_CHL_ICE_RHO_MATRIX==0) = NaN; 

NAN_randn_CHL_ICE_RHO_MATRIX = nan(size(randn_CHL_ICE_RHO_MATRIX)); 
NAN_randn_CHL_ICE_RHO_MATRIX(~isnan(test_CHL_ICE_RHO_MATRIX)) = randn_CHL_ICE_RHO_MATRIX(~isnan(test_CHL_ICE_RHO_MATRIX)); 


% CMOD & CHL
test_CMOD_CHL_RHO_MATRIX = CMOD_CHL_RHO_MATRIX; 
test_CMOD_CHL_RHO_MATRIX(test_CMOD_CHL_RHO_MATRIX==0) = NaN;

NAN_randn_CMOD_CHL_RHO_MATRIX = nan(size(test_CMOD_CHL_RHO_MATRIX)); 
NAN_randn_CMOD_CHL_RHO_MATRIX(~isnan(test_CMOD_CHL_RHO_MATRIX)) = randn_CMOD_CHL_RHO_MATRIX(~isnan(test_CMOD_CHL_RHO_MATRIX)); 


% CMOD & WIND
test_CMOD_WIND_RHO_MATRIX = CMOD_WIND_RHO_MATRIX; 
test_CMOD_WIND_RHO_MATRIX(CMOD_WIND_RHO_MATRIX == 0) = NaN; 

NAN_randn_CMOD_WIND_RHO_MATRIX = nan(size(test_CMOD_WIND_RHO_MATRIX)); 
NAN_randn_CMOD_WIND_RHO_MATRIX(~isnan(test_CMOD_WIND_RHO_MATRIX)) = randn_CMOD_WIND_RHO_MATRIX(~isnan(test_CMOD_WIND_RHO_MATRIX));


%%
% trying to keep NaNs in same position as where ice was removed so that
% random matrices statistics can be computed
nx = length(CMOD_ICE_RHO_MATRIX(:,1,1));
ny = length(CMOD_ICE_RHO_MATRIX(1,:,1));
nz = length(CMOD_ICE_RHO_MATRIX(1,1,:));

CMOD_ICE_DATA_vector = reshape(CMOD_ICE_RHO_MATRIX, [nx*ny*nz 1]);
CMOD_CHLA_DATA_vector = reshape(CMOD_CHL_RHO_MATRIX, [nx*ny*nz 1]);
CMOD_WIND_DATA_vector = reshape(CMOD_WIND_RHO_MATRIX, [nx*ny*nz 1]); 
ICE_CHLA_DATA_vector = reshape(test_CHL_ICE_RHO_MATRIX, [nx*ny*nz 1]); 

% RAND_vector1 = reshape(randn_CMOD_CHL_RHO_MATRIX, [nx*ny*nz 1] ) ;
RAND_vector1 = reshape(NAN_randn_CMOD_ICE_RHO_MATRIX, [nx*ny*nz 1]); 
RAND_vector2 = reshape(NAN_randn_CHL_ICE_RHO_MATRIX, [nx*ny*nz 1]); 
RAND_vector3 = reshape(NAN_randn_CMOD_CHL_RHO_MATRIX, [nx*ny*nz 1]); 
RAND_vector4 = reshape(NAN_randn_CMOD_WIND_RHO_MATRIX, [nx*ny*nz 1]); 

% Check assumptions before conducting ANOVA / Kruskal - Wallis

% If h = 1, this indicates the rejection of the null hypothesis at the Alpha significance level.
% If h = 0, this indicates a failure to reject the null hypothesis at the Alpha significance level.

% returns a test decision for the null hypothesis that 
% the data in vector x is from a population with a normal distribution, 
% using the Anderson-Darling test

% For ANOVA, samples need to come from a normal distribution. if h = 0, you
% fail to reject the null that the samples come from a normal distribution,
% and can continue with the ANOVA. 

h_rand1     = adtest(RAND_vector1); 
h_rand2     = adtest(RAND_vector2); 
h_rand3     = adtest(RAND_vector3);
h_rand4     = adtest(RAND_vector4);
h_CMOD_wind = adtest(CMOD_WIND_DATA_vector); 
h_CMOD_ice  = adtest(CMOD_ICE_DATA_vector); 
h_CMOD_chla = adtest(CMOD_CHLA_DATA_vector); 
h_ice_chla  = adtest(ICE_CHLA_DATA_vector); 

%%


% save('DATA.mat', 'DATA'); % it's now saved for future work and saving these figures and stuff

% DATA 


DATA = cat(2,RAND_vector1, RAND_vector2, RAND_vector3, RAND_vector4, CMOD_ICE_DATA_vector,...
    CMOD_CHLA_DATA_vector, CMOD_WIND_DATA_vector, ICE_CHLA_DATA_vector);


groupnames = {'random MOD & ice', 'random ice & chl-a', 'random MOD & chl-a','random MOD & wind',...
    'MOD & ice', 'MOD & chl-a', 'MOD & wind', 'ice & chl-a'}; 


[p_anova,tbl_anova, stats_anova]                         = anova1(DATA, groupnames);

% p_kruskalwallis = kruskalwallis(DATA)

% 
% Fstat  = tbl{2,5};
% pvalue = tbl{2,6};

%%

fig = figure; clf; 
[c_anova, m_anova, h_anova, nms_anova]  = multcompare(stats_anova); 
set(gca, 'FontSize', 16);
title('Multiple Comparison Plot, ANOVA'); 
grid off

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 


%%
print(fig,'Multiple_Comparison_Plot_ANOVA.png','-dpng','-r96');       %  *// 300 dpi


writematrix(c_anova,'c_anova.xls')

%%

fig = figure;clf; 
[p_kruskalwallis, tbl_kruskalwallis, stats_kruskalwallis] = kruskalwallis(DATA, groupnames);
set(gca, 'FontSize', 22);
xtickangle(30)

fig = figure; clf;
[c_kruskalwallis, m_kruskalwallis,...
    h_kruskalwallis, nms_kruskalwallis] = multcompare(stats_kruskalwallis); 

set(gca, 'FontSize', 22);
title('Multiple Comparison Plot, Kruskal-Wallis'); 
xlabel('Testing whether mean ranks of groups are significantly different')
% ylabel('Groups')
% grid on

%%

writematrix(c_kruskalwallis,'c_kruskalwallis.xls')
%%
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'Multiple_Comparison_Plot_KruskalWallis.png','-dpng','-r300');       %  *// 300 dpi


%%
% writematrix(c_kruskalwallis,'c_kruskalwallis.xls')
% c returns a matrix of the pairwise comparison results from a multiple
% comparison test using information contained in stats structure


%%
% 
% 
% DATA = cat(2,RAND_vector1, RAND_vector2, RAND_vector3, RAND_vector4, CMOD_ICE_DATA_vector,...
%     CMOD_CHLA_DATA_vector, CMOD_WIND_DATA_vector, ICE_CHLA_DATA_vector);
% 
% groupnames = {'random MOD & ice', 'random ice & chl-a', 'random MOD & chl-a','random MOD & wind', 'MOD & ice', 'MOD & chl-a', 'MOD & wind', 'ice & chl-a'}; 


[p,tbl,stats] = anova1(DATA, groupnames);
[p,tbl,stats] = kruskalwallis(DATA, groupnames);

Fstat  = tbl{2,5};
pvalue = tbl{2,6};

[c, m, h, nms] = multcompare(stats); 





% [h, p] = ttest(CMOD_ICE_DATA_vector, RAND_vector) ;
% [h, p] = ttest(CMOD_CHLA_DATA_vector, RAND_vector) ;

%% this was from a help document regarding 3D arrays and seeing if
% significant differences arise. 
a = rand(52,21,49);
b = rand(52,21,49);
x = permute(a,[3 1 2]);
y = permute(b, [3 1 2]); 

[h,p] = ttest2(permute(a,[3 1 2]),permute(b,[3 1 2]));
testh = squeeze(h);
testp = squeeze(p);

[h, p] = ttest(rand(10,1), rand(10,1)) ;












%% % Plotting of spatial statistical distributions given random numbers 

Lon = LON_aqua_step_edit;
Lat = LAT_aqua_step_edit;

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.02], [0.05 0.05], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure;clf;

subplot(2,2,4) %fall
plot_SO_figure_subplot(NAN_randn_CMOD_ICE_RHO_MATRIX(:,:,1), ...
    Lon,...
    Lat,...
    redblue,...
    [-1 1]); 

    hold on;


    [C, h] = m_contour(Lon, Lat, fall_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65 ; 
%     clabel(C, h,'fontsize',13);
%     clegendm(C,h,'m')

subplot(2,2,1) % winter
plot_SO_figure_subplot(NAN_randn_CMOD_ICE_RHO_MATRIX(:,:,2), ...
    Lon,...
    Lat,...
    redblue,...
    [-1 1]); 

    hold on;

    [C, h] = m_contour(Lon, Lat, winter_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2,2,2) % spring
plot_SO_figure_subplot(NAN_randn_CMOD_ICE_RHO_MATRIX(:,:,3), ...
    Lon,...
    Lat,...
    redblue,...
    [-1 1]); 
    hold on;

    [C, h] = m_contour(Lon, Lat, spring_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2,2,3) %summer
plot_SO_figure_subplot(NAN_randn_CMOD_ICE_RHO_MATRIX(:,:,4), ...
    Lon,...
    Lat,...
    redblue,...
    [-1 1]); 

    hold on;
    [C, h] = m_contour(Lon, Lat, summer_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
 hy = ylabel(h, 'Random number MOD & Ice seasonal correlation coefficient', 'FontSize', 18);





%%
cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'RandomNumber_Seasonal_Spatial_Correlation_MOD_ICE_1_1_masked_newcolormap.png','-dpng','-r300');       %  *// 300 dpi

%%

fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 
subplot(3, 1, 1 ) % spring
plot_SO_figure_subplot(NAN_randn_CHL_ICE_RHO_MATRIX(:,:,3), ...
    Lon,...
    Lat,...
    redblue, ...
    [-1 1]); 

    hold on;

    [C, h] = m_contour(Lon, Lat, spring_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(3,1, 2) % summer
plot_SO_figure_subplot(NAN_randn_CHL_ICE_RHO_MATRIX(:,:,4), ...
    Lon,...
    Lat,...
    redblue, ...
    [-1 1]); 

    hold on;

    [C, h] = m_contour(Lon, Lat, summer_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(3,1, 3) % fall
plot_SO_figure_subplot(NAN_randn_CHL_ICE_RHO_MATRIX(:,:,1), ...
    Lon,...
    Lat,...
    redblue, ...
    [-1 1]); 

    hold on;
    [C, h] = m_contour(Lon, Lat, fall_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

 
hp4 = get(subplot(3,1,3),'Position');
h = colorbar('Position', [0.65 0.2  0.02  0.6]);
h.FontWeight = 'bold';
h.FontSize = 15;


 hy = ylabel(h, 'Random Number ice & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
% hy = ylabel(h, 'CMOD & Ice seasonal correlation coefficient', 'FontSize', 18);



%%
cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'RandomNumber_Seasonal_Spatial_Correlation_ice_chla_1_1_masked_newcolormap.png','-dpng','-r300');       %  *// 300 dpi
%  



%%

fig = figure;clf;

subplot(2,2,4) %fall
plot_SO_figure_subplot(NAN_randn_CMOD_CHL_RHO_MATRIX(:,:,1), ...
    Lon,...
    Lat,...
    redblue, ...
    [-1 1]);

hold on;

[C, h] = m_contour(Lon, Lat, fall_contour, 'linewi', 3, 'LineColor', [ 0 0 0]);
h.LevelList = 0.65  ;


subplot(2,2,1) % winter
plot_SO_figure_subplot(NAN_randn_CMOD_CHL_RHO_MATRIX(:,:,2), ...
    Lon,...
    Lat,...
    redblue, ...
    [-1 1]);

hold on;

[C, h] = m_contour(Lon, Lat, winter_contour, 'linewi', 3, 'LineColor', [ 0 0 0]);
h.LevelList = 0.65  ;

subplot(2,2,2) %  spring
plot_SO_figure_subplot(NAN_randn_CMOD_CHL_RHO_MATRIX(:,:,3), ...
    Lon,...
    Lat,...
    redblue, ...
    [-1 1]);

hold on;

[C, h] = m_contour(Lon, Lat, spring_contour, 'linewi', 3, 'LineColor', [ 0 0 0]);
h.LevelList = 0.65  ;


subplot(2,2,3) % summer
plot_SO_figure_subplot(NAN_randn_CMOD_CHL_RHO_MATRIX(:,:,4), ...
    Lon,...
    Lat,...
    redblue, ...
    [-1 1]);

hold on;

[C, h] = m_contour(Lon, Lat, summer_contour, 'linewi', 3, 'LineColor', [ 0 0 0]);
h.LevelList = 0.65  ;

hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;

hy = ylabel(h, 'Random Number MOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);


%%
cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'RandomNumber_Seasonal_Spatial_Correlation_MOD_chla_1_1_masked_newcolormap.png','-dpng','-r300');       %  *// 300 dpi
%  


%%

fig = figure;clf;

subplot(2,2,4) %fall
plot_SO_figure_subplot(NAN_randn_CMOD_WIND_RHO_MATRIX(:,:,1), ...
    Lon,...
    Lat,...
    redblue, ...
    [-1 1]); 

    hold on;

    [C, h] = m_contour(Lon, Lat, fall_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 

subplot(2,2,1) % winter
plot_SO_figure_subplot(NAN_randn_CMOD_WIND_RHO_MATRIX(:,:,2), ...
    Lon,...
    Lat,...
    redblue, ...
    [-1 1]); 
    hold on;

    [C, h] = m_contour(Lon, Lat, winter_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 

subplot(2,2,2) % spring
plot_SO_figure_subplot(NAN_randn_CMOD_WIND_RHO_MATRIX(:,:,3), ...
    Lon,...
    Lat,...
    redblue, ...
    [-1 1]); 

    hold on;
    [C, h] = m_contour(Lon, Lat, spring_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 

 subplot(2,2,3) % summer
plot_SO_figure_subplot(NAN_randn_CMOD_WIND_RHO_MATRIX(:,:,4), ...
    Lon,...
    Lat,...
    redblue, ...
    [-1 1]); 

    hold on;
    [C, h] = m_contour(Lon, Lat, summer_contour, 'linewi', 3, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 

 

    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;

 hy = ylabel(h, 'Random Number MOD & Wind Speed seasonal correlation coefficient', 'FontSize', 18);
% hy = ylabel(h, 'CMOD & Ice seasonal correlation coefficient', 'FontSize', 18);





%%
cd /Users/srishtidasarathy/Documents/Bowman/PhD_Phase_Two_SouthernOcean/Figures/Updated_figures_2007_2020_vars

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'RandomNumber_Seasonal_Spatial_Correlation_MOD_wind_1_1_masked_newcolormap.png','-dpng','-r300');       %  *// 300 dpi
%  

%%

% CMOD movie to check how clear data are and whether it needs interpolating 

% ntimes = length(CHLA_ONEDEGREERES_MONTHS(1,1,:));

% s = []; 

%       Now lets make movies to compare how these look

% contour_array = {winter_contour ; spring_contour ; summer_contour; fall_contour};

%%
figure(1), clf

vidfile = VideoWriter('ice_Season_TEST_1by1.mp4','MPEG-4');
vidfile.FrameRate = 1;
open(vidfile);

for i = 1 : 51
    disp(i)
    
    plot_SO_figure(SeaIce_SEASONS(:,:,i),...
        LON_aqua_step_edit,...
        LAT_aqua_step_edit,...
        'ice',...
        [0 1.2],...
        (['ice: ', datestr(times_seasons(i))]))
    
    
%     hold on;
%     [C, h] = m_contour(Lon, Lat, contour_array(i,:), 'linewi', 3, 'LineColor', [0 1 0.6]);
%     h.LevelList = 0.65  ;
%     %     clabel(C, h,'fontsize',13);
%     
%     
%     
%     hp4 = get(subplot(3,1, 3),'Position');
%     h = colorbar('Position', [hp4(1)+(hp4(3)-0.2)  0.33  0.02  0.4]);
%     h.FontWeight = 'bold';
%     h.FontSize = 15;
%     
    
    
    
    drawnow
    
    F(i) = getframe(gcf);
    writeVideo(vidfile,F(i));
    
    
    
end

close(vidfile)













