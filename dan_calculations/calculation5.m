function [out_variable, out_variable_names] = calculation5(Data_in) ;
% we are calculating Mean daily temperature range here.
% can set two functions coppesding to 5%DB and 5%WB respectively

%% Define output variable names
out_variable_names = {'M_DB_Range', 'Monthly_DB_db', 'Monthly_DB_wb', ...
    'Monthly_WB_db', 'Monthly_WB_wb'} ; 

% You should not need to change anything below this point

% Deal out input variables
in_variable_names = {'yr', 'mo', 'hr', 'day', 'tmp', 'dew', 'wndmag', 'wnddir', 'tstar2'} ;
[yr, mo, day, hr] = datevec(Data_in.time-.25) ; 
for ivarn = 5:length(in_variable_names) ; 
    eval([in_variable_names{ivarn} ' = Data_in.' in_variable_names{ivarn} ';']) ; 
end 
temp = tmp ; 

tstar_f = tstar2 - 459.67 ; 
tstar_c = (tstar_f - 32).*(5/9) ;  % converting to celsius degree 

%% Mean daily range calculation corresponding to 5% DB/WB
% for month with 31 days 1,3,5,7,8,10,12
Daily_DB_31 = zeros(31*25,1) ; % 31days*25years
Daily_WB_31 = zeros(31*25,1) ;
Range_DB_31 = zeros(31*25,1) ;
Range_WB_31 = zeros(31*25,1) ;
max_t_31 = zeros(31*25,1) ;
max_tw_31 = zeros(31*25,1) ;
Monthly_DB_db_31 = zeros(7,1) ;
Monthly_DB__wb_31 = zeros(7,1) ;
Monthly_WB_db_31 = zeros(7,1) ;
Monthly_WB__wb_31 = zeros(7,1) ;
M_DB_Range_31 = zeros(7,1) ;
i = 1;
for imo = [1 3 5 7 8 10 12]
    j = 1;
    ind = find(mo==imo );
    t = temp(ind);  % hourly data for iyr imo
    tw = tstar_c(ind);
    
    % t for DB 
    ntim = sum(~isnan(t)) ;
    t_sort = sort(t(~isnan(t)), 'descend') ;
    threshold = floor(0.05*ntim) ;
    t_threshold = t_sort(threshold); % the value for 5% DB
    
    % tw for WB
    ntim = sum(~isnan(tw)) ;
    tw_sort = sort(tw(~isnan(tw)), 'descend') ;
    threshold = floor(0.05*ntim) ;
    tw_threshold = tw_sort(threshold); % the value for 5% DB
    
    for iyr = [1986:2010]
        for iday = 1:31  %31days fro July 
            indd = find(mo==imo & yr==iyr & day == iday) ;
            t = temp(indd) ; %hourly data for July
            tw = tstar_c(indd) ; 
            Range_DB_31(j) = max(t) - min(t);
            Range_WB_31(j) = max(tw) - min(tw);
            max_t_31(j) = max(t) ;
            max_tw_31(j) = max(tw) ;
            j = j+1;
        end 
        
    end
    M_DB_Range_31(i) = mean(Range_DB_31) ; 
    ind_D = find(max_t_31 >= t_threshold) ;  % 5% DB
    ind_W = find(max_tw_31 >= tw_threshold) ;  % 5% WB 
    
    Monthly_DB_db_31(i) = mean(Range_DB_31(ind_D));
    Monthly_DB_wb_31(i) = mean(Range_WB_31(ind_D));
    Monthly_WB_db_31(i) = mean(Range_DB_31(ind_W));
    Monthly_WB_wb_31(i) = mean(Range_WB_31(ind_W));
    
    i = i+1 ;
end

%% the months with 30 days: 4 6 9 11
Daily_DB_30 = zeros(30*25,1) ; 
Daily_WB_30 = zeros(30*25,1) ; 
Range_DB_30 = zeros(30*25,1) ; 
Range_WB_30 = zeros(30*25,1) ; 
max_t_30 = zeros(30*25,1) ;
max_tw_30 = zeros(30*25,1) ;
M_DB_Range_30 = zeros(4,1) ;
Monthly_DB__wb_30 = zeros(4,1) ;
Monthly_WB_db_30 = zeros(4,1) ;
Monthly_WB__wb_30 = zeros(4,1) ;
M_DB_Range_30 = zeros(4,1) ;

i = 1;
for imo = [4 6 9 11]
    j = 1;
    ind = find(mo==imo );
    t = temp(ind);  % hourly data for iyr imo
    tw = tstar_c(ind);
    
    % t for DB 
    ntim = sum(~isnan(t)) ;
    t_sort = sort(t(~isnan(t)), 'descend') ;
    threshold = floor(0.05*ntim) ;
    t_threshold = t_sort(threshold); % the value for 5% DB
    
    % tw for WB
    ntim = sum(~isnan(tw)) ;
    tw_sort = sort(tw(~isnan(tw)), 'descend') ;
    threshold = floor(0.05*ntim) ;
    tw_threshold = tw_sort(threshold); % the value for 5% DB
    
    for iyr = [1986:2010]
        for iday = 1:30  %31days fro July 
            indd = find(mo==imo & yr==iyr & day == iday) ;
            t = temp(indd) ; %hourly data for July
            tw = tstar_c(indd) ; 
            Range_DB_30(j) = max(t) - min(t);
            Range_WB_30(j) = max(tw) - min(tw);
            max_t_30(j) = max(t) ;
            max_tw_30(j) = max(tw) ;
            j = j+1;
        end 
        
    end
    M_DB_Range_30(i) = mean(Range_DB_30) ; 
    ind_D = find(max_t_30 >= t_threshold) ;  % 5% DB
    ind_W = find(max_tw_30 >= tw_threshold) ;  % 5% WB 
    
    Monthly_DB_db_30(i) = mean(Range_DB_30(ind_D));
    Monthly_DB_wb_30(i) = mean(Range_WB_30(ind_D));
    Monthly_WB_db_30(i) = mean(Range_DB_30(ind_W));
    Monthly_WB_wb_30(i) = mean(Range_WB_30(ind_W));
 
    i = i+1 ;
end

%% Feb  ndays = 28+(mod(year, 4) == 0)
days = zeros(25,1) ;  
i = 1 ;
for iyr = [1986:2010]
    days(i) = 28+(mod(iyr, 4) == 0);
    i = i+1 ; 
end
sum(days) ; % checking how many days in total 24 Feb


ind = find(mo==2 );
t = temp(ind);  % hourly data for iyr imo
tw = tstar_c(ind);
% t for DB 
ntim = sum(~isnan(t)) ;
t_sort = sort(t(~isnan(t)), 'descend') ;
threshold = floor(0.05*ntim) ;
t_threshold = t_sort(threshold); % the value for 5% DB   
% tw for WB
ntim = sum(~isnan(tw)) ;
tw_sort = sort(tw(~isnan(tw)), 'descend') ;
threshold = floor(0.05*ntim) ;
tw_threshold = tw_sort(threshold); % the value for 5% DB

Daily_DB = zeros(sum(days),1) ;
Daily_WB = zeros(sum(days),1) ;
Range_DB = zeros(sum(days),1) ;
Range_WB = zeros(sum(days),1) ; 
max_t = zeros(sum(days),1) ;
max_tw = zeros(sum(days),1) ;
j =  1 ; 
for iyr = [1986:2010]
    i = 1 ; 
    for iday  = 1:days(i)
        indd = find(mo==2 & yr==iyr & day == iday);
        t = temp(indd); 
        tw = tstar_c(indd);
        
        Range_DB(j) = max(t) - min(t);
        Range_WB(j) = max(tw) - min(tw);
        max_t(j) = max(t) ;
        max_tw(j) = max(tw) ;
        j = j+1;
    end 
    i = i+1 ; 
end

M_DB_Range_2 = mean(Range_DB);
ind_D = find(max_t >= t_threshold) ;  % 5% DB
ind_W = find(max_tw >= tw_threshold) ;  % 5% WB 
    
Monthly_DB_db_2 = mean(Range_DB(ind_D));
Monthly_DB_wb_2 = mean(Range_WB(ind_D));
Monthly_WB_db_2 = mean(Range_DB(ind_W));
Monthly_WB_wb_2 = mean(Range_WB(ind_W));
 
%%
M_DB_Range = nan(12,1) ; 
Monthly_DB_db = nan(12,1) ; 
Monthly_DB_wb = nan(12,1) ; 
Monthly_WB_db = nan(12,1) ; 
Monthly_WB_wb = nan(12,1) ; 

M_DB_Range([1 3 5 7 8 10 12]) = M_DB_Range_31 ; 
M_DB_Range([4 6 9 11]) = M_DB_Range_30 ; 
M_DB_Range(2) = M_DB_Range_2 ; 

Monthly_DB_db([1 3 5 7 8 10 12]) = Monthly_DB_db_31 ; 
Monthly_DB_db([4 6 9 11]) = Monthly_DB_db_30 ; 
Monthly_DB_db(2) = Monthly_DB_db_2 ; 

Monthly_DB_wb([1 3 5 7 8 10 12]) = Monthly_DB_wb_31 ; 
Monthly_DB_wb([4 6 9 11]) = Monthly_DB_wb_30 ; 
Monthly_DB_wb(2) = Monthly_DB_wb_2 ;

Monthly_WB_db([1 3 5 7 8 10 12]) = Monthly_WB_db_31 ; 
Monthly_WB_db([4 6 9 11]) = Monthly_WB_db_30 ; 
Monthly_WB_db(2) = Monthly_WB_db_2 ;

Monthly_WB_wb([1 3 5 7 8 10 12]) = Monthly_WB_wb_31 ; 
Monthly_WB_wb([4 6 9 11]) = Monthly_WB_wb_30 ; 
Monthly_WB_wb(2) = Monthly_WB_wb_2 ;

%% Deal output variables into out_variable

out_variable = cell(length(out_variable_names), 1) ; 
for ivarn = 1:length(out_variable_names) ; 
    eval(['out_variable{ivarn} = ' out_variable_names{ivarn} ';']) ; 
end
