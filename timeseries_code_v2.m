
%%
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PhD_Phase_Two_SouthernOcean/Variables/2017_2020_vars

 load('Total_timetable_SO_MOD_NEW.mat')
%  load('Total_timetable_SO_Depol_Ratio.mat')

% I will instead be getting ice values by averaging the monthly grid that I
% made so that land values and water values are masked out :)
% This is in Monthly_gridding_MOD_Ice_v2.m

load('Southern_Ocean_Wind_Speed_Data_ONEDEGREERES.mat')
load('Aqua.mat')
% load('Monthly_grid_MOD_ice.mat');
load('Monthly_grid__NEW_MOD_ice_withopenocean.mat')
load('Monthly_grid__NEW_MOD_ice.mat')

% 
 Master_Wind_Speed(Master_Wind_Speed == -500) = nan;
% 
% save('Southern_Ocean_Wind_Speed_Data_ONEDEGREERES.mat', 'Master_Wind_Speed', 'LAT_wind', 'LON_wind', 'Master_filename', '-v7.3');

%%

t1 = datetime(2007,01,01);
t2 = datetime(2020,12,31);
times = t1:calmonths(1):t2; 

times = times';

%%
% Monthly averaged time series code below

Total_chl_a_monthly = zeros(1, 168); % from length of times array above
Total_wind_monthly = zeros(1, 168);
Total_ice_monthly = zeros(1, 168);
Total_ice_plusopenocean_monthly = zeros(1, 168);

for i = 1:168
    disp(i)
    Total_chl_a_monthly(i) = nanmean(Master_chl_a_full_res_SouthernOcean(:,:,i), [1 2]); % Find the mean of each page of X by specifying dimensions 1 and 2 as the operating dimensions.
    Total_wind_monthly(i)  = nanmean(Master_Wind_Speed(:,:,i), [1 2]);
    Total_ice_monthly(i)   = nanmean(SeaIce_MONTHS(:,:,i) , [1 2]);
    Total_ice_plusopenocean_monthly(i) = nanmean(SeaIce_openocean_MONTHS(:,:,i) , [1 2]);
end

% Mean absolute deviation calculation: 
chl_mad = zeros(1, 168);
wind_mad = zeros(1, 168);
ice_mad = zeros(1, 168);
ice_plusopenocean_mad = zeros(1, 168);


for i = 1:168
    disp(i)
    chl_mad(i) = mad(Master_chl_a_full_res_SouthernOcean(:,:,i), 0,[1 2]); % Find the mean of each page of X by specifying dimensions 1 and 2 as the operating dimensions.
    wind_mad(i)  = mad(Master_Wind_Speed(:,:,i), 0,[1 2]);
    ice_mad(i)   = mad(SeaIce_MONTHS(:,:,i) , 0,[1 2]);
    ice_plusopenocean_mad(i) = mad(SeaIce_openocean_MONTHS(:,:,i), 0,[1 2]);
end

%%
timetable_MOD_monthly_avg = retime(Total_timetable_SO_MOD_NEW, 'monthly', @nanmean);
MOD_Monthly_avg           = timetable_MOD_monthly_avg.CMOD_Surface;
MOD_Time_Months           = timetable_MOD_monthly_avg.Total_Profile_Time_New_Surface;
MOD_Lat_Months            = timetable_MOD_monthly_avg.Total_Latitude_Surface;
MOD_Lon_Months            = timetable_MOD_monthly_avg.Total_Longitude_Surface;


timetable_MOD_mad = retime(Total_timetable_SO_MOD_NEW, 'monthly', @mad);
MOD_mad           = timetable_MOD_monthly_avg.CMOD_Surface;
MOD_Time_Months   = timetable_MOD_monthly_avg.Total_Profile_Time_New_Surface;

chl_mad = chl_mad';
ice_plusopenocean_mad =  ice_plusopenocean_mad';
ice_mad = ice_mad';
wind_mad = wind_mad';

mad_timetable = timetable(MOD_Time_Months, MOD_mad, chl_mad, ice_plusopenocean_mad, ice_mad, wind_mad);

save('mad_timetable.mat' , 'mad_timetable', '-v7.3');

% 
% add_var = NaN;
% MOD_Monthly_avg = cat(1, MOD_Monthly_avg, add_var); % I'm adding a NaN because when i pulled data from calipso last month in 2020 was missing.

% Ice_only = Total_timetable_SO_Depol_Ratio.Total_Surface_532_Integrated_Depolarization_Ratio;
% Ice_only(Ice_only <= 0.15) = NaN; 
% 
% Total_timetable_SO_Depol_Ratio = addvars(Total_timetable_SO_Depol_Ratio, Ice_only); 
% 
% timetable_Depol_Ratio_monthly_avg = retime(Total_timetable_SO_Depol_Ratio, 'monthly', @nanmean);
% Depol_Ratio_Monthly_avg           = timetable_Depol_Ratio_monthly_avg.Total_Surface_532_Integrated_Depolarization_Ratio;
% Ice_only_Monthly_avg              = timetable_Depol_Ratio_monthly_avg.Ice_only;
% 
% Depol_Ratio_Monthly_avg = cat(1, Depol_Ratio_Monthly_avg, add_var); %
% Ice_only_Monthly_avg = cat(1, Ice_only_Monthly_avg, add_var); %
% 
% 
% Depol_Ratio_Time_Months           = timetable_Depol_Ratio_monthly_avg.Total_Profile_Time_New_Ice;
% Depol_Ratio_Lat_Months            = timetable_Depol_Ratio_monthly_avg.Total_Latitude_Ice;
% Depol_Ratio_Lon_Months            = timetable_Depol_Ratio_monthly_avg.Total_Longitude_Ice;

 


%%
% 
% save('monthly_avg_vars.mat', 'MOD_Monthly_avg', 'Total_chl_a_monthly', 'Total_ice_monthly', 'Total_wind_monthly', '-v7.3')
% writematrix(MOD_Monthly_avg, 'MOD.csv')
% writematrix(Total_chl_a_monthly, 'chl.csv')
% writematrix(Total_ice_monthly, 'seaice.csv')
% writematrix(Total_wind_monthly, 'windspeed.csv')
% writematrix(times, 'months.csv')


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
 
%  figure(1); clf; % open ocean masked for ice time series
 figure(2);clf;
 
 ax(1) = subplot(4,1,1);
 aa_splot(x, MOD_Monthly_avg, '-',...
     'linewidth', 1.5, ...
     'Color', black);
 
 xlim([x(1) x(end)]) 
 ylim([0 0.1])
 
 ylabel('MOD')
% l(1)= legend('MOD', 'Position',[0.76 0.92 0.1 0.0349]);

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
 
 % ICE, or switch out to variable that includes open ocean pixels
 ax(3) = subplot(4,1,3);
 aa_splot(x,  Total_ice_plusopenocean_monthly,'-',...
     'linewidth', 1.5,...
     'MarkerSize', 4,...
     'MarkerFaceColor', ice_blue,...
     'MarkerEdgeColor', ice_blue,...
     'Color', ice_blue);
  xlim([x(1)  x(end)]) 

%  ylim([0.4  0.8]) % this is for when open ocean is masked
 ylim([0.05 0.5]) % for when open ocean is included
 
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

print(gcf,'TimeSeries_openoceanincluded.png','-dpng','-r300');       %  *// 300 dpi


%%


cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication

load('Depol_Ratio_Monthly_avg_Vars.mat')
load('amsrmf_Monthly_avg_Vars.mat')
load('CMOD_Monthly_avg_Vars_Surface.mat') 


cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/Updated_all_Figures/
% 
% t1 = datetime(2006,06,01);
% t2 = datetime(2018,12,31);
% times = t1:calmonths(1):t2; 
% 


t1 = datetime(2007,01,01);
t2 = datetime(2007,12,31);
times_months = t1:calmonths(1):t2; 
times_months_num = month(times_months);

times = times'; 
Ice_Vector = double(Depol_Ratio_Monthly_avg(:,1));
% Ice_Vector(117,:) = [];
CMOD_Vector = CMOD_Monthly_avg_Surface(:,1);
% CMOD_Vector(117,:) = [];
Wind_Vector = amsrmf_Monthly_avg;
% times_without_feb = times;
% times_without_feb(117,:) = []; %
% times_for_wind = times;
% times_for_wind(117,:) =[];
% times_for_wind(64:73,:) = [];
% Wind_Vector(117,:) = [];
% Wind_Vector(64:73,:)=[];



[Ac_CMOD, tc_CMOD] = climatology(CMOD_Vector,times_without_feb, 'monthly'); 
[Ac_Ice, tc_Ice] = climatology(Ice_Vector, times_without_feb, 'monthly'); 
[Ac_Wind, tc_Wind] = climatology(Wind_Vector, times_for_wind, 'monthly'); 

Master_chl_a(isnan(Master_chl_a)) = 0; 
[Ac_Chl, tc_Chl] = climatology(Master_chl_a, times, 'monthly'); 

for i = 1:12
Ac_Chl_monthly(i) = mean2(Ac_Chl(:,:,i));
end

Ac_Chl_monthly(Ac_Chl_monthly <= 0) = NaN;


%%
% this is for standard deviation calculation, check cdt monthly to see what
% I did 

jan_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[1], @mad);
feb_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[2], @mad);
mar_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[3], @mad);
apr_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[4],  @mad);
may_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[5],  @mad);
jun_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[6],  @mad);
jul_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[7],  @mad);
aug_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[8],  @mad);
sep_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[9],  @mad);
oct_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[10],  @mad);
nov_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[11],  @mad);
dec_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[12],  @mad);

[Ac_Ice, tc_Ice] = climatology(Ice_Vector, times_without_feb, 'monthly'); 


jan_Ice_mad = monthly(Ice_Vector,times_without_feb,[1],  @mad);
feb_Ice_mad = monthly(Ice_Vector,times_without_feb,[2],  @mad);
mar_Ice_mad = monthly(Ice_Vector,times_without_feb,[3],  @mad);
apr_Ice_mad = monthly(Ice_Vector,times_without_feb,[4],  @mad);
may_Ice_mad = monthly(Ice_Vector,times_without_feb,[5],  @mad);
jun_Ice_mad = monthly(Ice_Vector,times_without_feb,[6],  @mad);
jul_Ice_mad = monthly(Ice_Vector,times_without_feb,[7],  @mad);
aug_Ice_mad = monthly(Ice_Vector,times_without_feb,[8],  @mad);
sep_Ice_mad = monthly(Ice_Vector,times_without_feb,[9],  @mad);
oct_Ice_mad = monthly(Ice_Vector,times_without_feb,[10],  @mad);
nov_Ice_mad = monthly(Ice_Vector,times_without_feb,[11],  @mad);
dec_Ice_mad = monthly(Ice_Vector,times_without_feb,[12],  @mad);


[Ac_Wind, tc_Wind] = climatology(Wind_Vector, times_for_wind, 'monthly'); 

jan_Wind_mad = monthly(Wind_Vector,times_for_wind,[1],  @mad);
feb_Wind_mad = monthly(Wind_Vector,times_for_wind,[2],  @mad);
mar_Wind_mad = monthly(Wind_Vector,times_for_wind,[3],  @mad);
apr_Wind_mad = monthly(Wind_Vector,times_for_wind,[4],  @mad);
may_Wind_mad = monthly(Wind_Vector,times_for_wind,[5],  @mad);
jun_Wind_mad = monthly(Wind_Vector,times_for_wind,[6],  @mad);
jul_Wind_mad = monthly(Wind_Vector,times_for_wind,[7],  @mad);
aug_Wind_mad = monthly(Wind_Vector,times_for_wind,[8],  @mad);
sep_Wind_mad = monthly(Wind_Vector,times_for_wind,[9],  @mad);
oct_Wind_mad = monthly(Wind_Vector,times_for_wind,[10],  @mad);
nov_Wind_mad = monthly(Wind_Vector,times_for_wind,[11],  @mad);
dec_Wind_mad = monthly(Wind_Vector,times_for_wind,[12],  @mad);

% chl and standard deviation calculation
Master_chl_a(isnan(Master_chl_a)) = 0; 
[Ac_Chl, tc_Chl] = climatology(Master_chl_a, times, 'monthly'); 

for i = 1:12
Ac_Chl_monthly(i) = mean2(Ac_Chl(:,:,i));
end

Ac_Chl_monthly(Ac_Chl_monthly <= 0) = NaN;

Ac_Chl_monthly = Ac_Chl_monthly';

for i = 1:12 
    disp(i)
    mad_chl_monthly_2d(:,:,i) = monthly(Master_chl_a, times, [i], 'dim', 3, @mad); 
end

    mad_chl_monthly = nanmean(mad_chl_monthly_2d, [1 2]);
    mad_chl_monthly = squeeze(mad_chl_monthly); 
    
% jan_chl_std = monthly(Master_chl_a,times,[1], 'omitnan', @std);
% jan_chl_std = 
% feb_chl_std = monthly(Master_chl_a,times,[2], 'omitnan', @std);
% mar_chl_std = monthly(Master_chl_a,times,[3], 'omitnan', @std);
% apr_chl_std = monthly(Master_chl_a,times,[4], 'omitnan', @std);
% may_chl_std = monthly(Master_chl_a,times,[5], 'omitnan', @std);
% jun_chl_std = monthly(Master_chl_a,times,[6], 'omitnan', @std);
% jul_chl_std = monthly(Master_chl_a,times,[7], 'omitnan', @std);
% aug_chl_std = monthly(Master_chl_a,times,[8], 'omitnan', @std);
% sep_chl_std = monthly(Master_chl_a,times,[9], 'omitnan', @std);
% oct_chl_std = monthly(Master_chl_a,times,[10], 'omitnan', @std);
% nov_chl_std = monthly(Master_chl_a,times,[11], 'omitnan', @std);
% dec_chl_std = monthly(Master_chl_a,times,[12], 'omitnan', @std);





%%
sienna = rgb('sienna'); 
gray = rgb('gray'); 
black = rgb('black'); 
grass_green = rgb('grass green');
ice_blue = rgb('sky blue'); 
wind_blue = rgb('scarlet');

% Create textbox
% pos = [x-start y-start x-width y-height]


%%

make_it_tight = true;
 subplot = @(m,n,p) subtightplot (m, n, p, [0.035 0.05], [0.1 0.05], [0.1 0.05]);
if ~make_it_tight,  clear subplot;  end


fig = figure(1), clf;

ax(1) = subplot(4,1,1);
grid on
aa_splot(1:12, Ac_CMOD, 'x-',...
    'linewidth', 2, ...
    'Color', black,...
    'MarkerSize', 9,...
    'MarkerFaceColor', black,...
    'MarkerEdgeColor', black);

ylabel('MOD')
%         xlabel('Months')
xlim([1,12])
ylim([0.01 0.09])
set(gca, 'ytick',0: 0.02: 0.08); 
AX=findall(0,'type','axes');
set(AX, 'FontSize', 18)
%                   xtickangle(20)


ax(2) = subplot(4,1,2);
grid on
aa_splot(1:12, Ac_Chl_monthly, 'v-',...
    'linewidth', 2, ...
    'Color', grass_green,...
    'MarkerSize', 7,...
    'MarkerFaceColor', grass_green,...
    'MarkerEdgeColor', grass_green);
% set(gca, 'xtick', 1:12,...
%     'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% title('Chl-{\ita} Climatology')
ylabel(sprintf('mg m^{-3}')),...
    %         xlabel('Months')
xlim([1,12])
ylim([0.0001 0.3])
set(gca, 'ytick',0: 0.1: 0.3); 


AX=findall(0,'type','axes');
set(AX, 'FontSize', 18)
% xtickangle(45)




ax(3) = subplot(4,1,3); 
grid on
aa_splot(1:12, Ac_Ice, 'o-',...
    'linewidth', 2, ...
    'Color', ice_blue,...
    'MarkerSize', 7,...
    'MarkerFaceColor', ice_blue,...
    'MarkerEdgeColor', ice_blue);
% set(gca, 'xtick', 1:12,...
%     'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% title('Ice Climatology')
ylabel('\delta',...
   'FontName','Helvetica Neue');%         xlabel('Months')
xlim([1,12])
ylim([0.4 0.7])
set(gca, 'ytick',0.4: 0.1: 0.7); 


AX=findall(0,'type','axes');
set(AX, 'FontSize', 18)
% xtickangle(45)




ax(4) = subplot(4,1,4);
grid on
aa_splot(1:12, Ac_Wind, 'd-',...
    'linewidth', 2, ...
    'Color', wind_blue,...
    'MarkerSize', 7,...
    'MarkerFaceColor', wind_blue,...
    'MarkerEdgeColor', wind_blue);
set(gca, 'xtick', 1:12,...
    'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% title('Wind Speed Climatology')
ylabel(sprintf('m s^{-1}'))
xlim([1,12])
ylim([6 13])

AX=findall(0,'type','axes');
set(AX, 'FontSize', 20)
 xtickangle(45)


[ax(1:3).XTickLabel] = deal([]);
ax(1).YTick(1) = []; 
ax(2).YTick(1) = []; 
ax(3).YTick(1) = []; 
ax(4).YTick(1) = [];
ax(4).FontSize = 23;
% ax(1) to see properties 

% from below: an(1).Position(2) = 0.87
% ax(1).Position(2) = 0.7673

% (0.7673 - should give the correct position for first subplot 

an(1) = annotation(gcf,'textbox',... % I drew this on in the figure 
  [0.11729695024077 0.915 0.16 0.0349],...
    'String','MOD',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(2) = annotation(gcf,'textbox',...
    [an(1).Position(1),...
    (ax(2).Position(2) + (an(1).Position(2) - ax(1).Position(2))),... % so this is the second axes, y position  + (y position of text box - y position of first axes)
    0.16 0.0349],...
    'String',{'Chl-{\ita}'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(3) = annotation(gcf, 'textbox', ...
    [an(1).Position(1),...
    (ax(3).Position(2) + (an(2).Position(2) - ax(2).Position(2))),...
    0.16 0.0349],...
    'String',{'Ice'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(4) = annotation(gcf, 'textbox', ...
    [an(1).Position(1),...
    (ax(4).Position(2) + (an(3).Position(2) - ax(3).Position(2))),...
    0.16 0.0349],...
    'String',{'Wind'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');


%%


% set(gcf,'PaperPositionMode','auto')
%  set(gcf,'PaperPosition','fillpage') 
orient(fig,'landscape')

print(fig,'Updated_Climatology_TimeSeries.png','-dpng','-r96');       %  *// 300 dpi
%  print(gcf, 'Ju


%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end
%% Upper and Lower Subplots with Titles
income = [3.2,4.1,5.0,5.6];
outgo = [2.5,4.0,3.35,4.9];
subplot(2,1,1); plot(income)
title('Income')
subplot(2,1,2); plot(outgo)
title('Outgo')

%% Subplots in Quadrants
figure
subplot(2,2,1)
text(.5,.5,{'subplot(2,2,1)';'or subplot 221'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,2)
text(.5,.5,{'subplot(2,2,2)';'or subplot 222'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,3)
text(.5,.5,{'subplot(2,2,3)';'or subplot 223'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,4)
text(.5,.5,{'subplot(2,2,4)';'or subplot 224'},...
    'FontSize',14,'HorizontalAlignment','center')













%%







figure(1), clf

vidfile = VideoWriter('algaetest.mp4','MPEG-4');
vidfile.FrameRate = 1;
open(vidfile);

for i = 1 : 168
    disp(i)
    
    plot_SO_figure(rot90(Master_chl_a_full_res_SouthernOcean(:,:,i)),...
        Longitude_Subset_SouthernOcean,...
        flip(Latitude_Subset_SouthernOcean),...
        'algae',...
        [0 1.2],...
        (['algae: ', datestr(times(i))]))
    
    
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





