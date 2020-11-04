function [out_variable, out_variable_names] = calculation(Data_in) ;
% we are calculating Parts of annual heating and humidification, 
% cooling and dehumidification design conditions, incluidng..
% hottest/coldest mont, hottest DB range,annual threahold exceedence...
% of WS, WD, DP and corresponding mean coincident DB.

%% Define output variable names
out_variable_names = {'m', 'db1', 'db2', 'ws1', 'mdb1', 'ws2', 'mdb2', 'mws1', 'mwd1', ...
    'x', 'mrange', 'mws3', 'mwd3'} ; 

% You should not need to change anything below this point

% Deal out input variables
in_variable_names = {'yr', 'mo', 'hr', 'day', 'tmp', 'dew', 'wndmag', 'wnddir'} ;
[yr, mo, day, hr] = datevec(Data_in.time-.25) ; 
for ivarn = 5:length(in_variable_names) ; 
    eval([in_variable_names{ivarn} ' = Data_in.' in_variable_names{ivarn} ';']) ; 
end 
temp = tmp ; 

%% Calculations : Coldes/hottest month 
mo_temp = zeros(12,1);
i = 1;
for imo = 1:12
    %for iyr = 1986:2010
    for iyr = 1986:2010
        ind = find(mo==imo & yr==iyr);
        t = temp(ind);  % hourly data for iyr imo
    end
    mo_temp(i) = mean(t);
    i = i+1;
end
[tem, x] = max(mo_temp(:)) ; % [tem, 7]:July
[tem, m] = min(mo_temp(:)) ; % [tem, 1]:January
%[x,y] = find(mo_temp==max(mo_temp(:))) ;% [7,1]:July
%[m,n] = find(mo_temp==min(mo_temp(:))) ;% [1,1]:January

%****************************************************************
% Hottest month DB range

if x==7&8
    DB_range = zeros(25*31,1);  %25years *31days
    j = 1;
    for iyr = [1986:2010]
        for iday = 1:31  %31days fro July 
            ind = find(mo==x & yr==iyr & day == iday);
            t = temp(ind);  %hourly data for July
            DB_range(j) = max(t) - min(t);
            j = j+1;
        end
    end
    mrange = mean(DB_range) ;
else
    DB_range = zeros(25*30,1);  %25years *31days
    j = 1;
    for iyr = [1986:2010]
        for iday = 1:30  %31days fro July 
            ind = find(mo==x & yr==iyr & day == iday);
            t = temp(ind);  %hourly data for July
            DB_range(j) = max(t) - min(t);
            j = j+1;
        end
    end
    mrange = mean(DB_range) ;
end
% coldest month temp/wind speed
ind_cold = find(mo == m) ;  
temp_cold = temp(ind_cold) ; 
wndmag_cold = wndmag(ind_cold) ; 

%% Calculations : Treshold_Exceedence
% 99.6% DB
[value] = temp ;
ntim = sum(~isnan(value)) ; %length(value) ;
value_sort = sort(value(~isnan(value)), 'descend') ;
threshold = floor(0.996*ntim) ;
value_threshold = value_sort(threshold);
db1 = value_threshold ;
% WS and WD cprresponding to 99.6% DB 
indd = find((value <= value_threshold + 5/18) & (value >= value_threshold - 5/18)) ;
mws1 = mean(wndmag(indd));
mwd1 = mean(wnddir(indd));


%%%%%%%%%%%%
for iyr = 1:nyr
             
    t = value(find(yr==years(iyr))); 
    ntim = sum(~isnan(value)) ; %length(value) ;
    value_sort = sort(value(~isnan(value)), 'descend') ;
    threshold = floor(0.996*ntim) ;
    value_threshold = value_sort(threshold);
    db1(iyr) = value_threshold ;
end

%%%%%%%%%%%


% 99% DB
[value] = temp ;
ntim = sum(~isnan(value)) ; %length(value) ;
value_sort = sort(value(~isnan(value)), 'descend') ;
threshold = floor(0.99*ntim) ;
value_threshold = value_sort(threshold);
db2 = value_threshold ;

% 0.4% DB
[value] = temp ;
ntim = sum(~isnan(value)) ; %length(value) ;
value_sort = sort(value(~isnan(value)), 'descend') ;
threshold = floor(0.004*ntim) ;
value_threshold = value_sort(threshold);
db3 = value_threshold ;
% WS and WD cprresponding to 0.4% DB 
indd = find((value <= value_threshold + 5/18) & (value >= value_threshold - 5/18)) ;
mws3 = mean(wndmag(indd))  ;
mwd3 = mean(wnddir(indd))  ;

% coldest_month WS/MCDB
[value] = wndmag_cold;
ntim = sum(~isnan(value)) ; %length(value) ;
value_sort = sort(value(~isnan(value)), 'descend') ;
threshold = floor(0.004*ntim) ;
value_threshold = value_sort(threshold);
ws1 = value_threshold ;
indd = find((value <= value_threshold + 5/18) & (value >= value_threshold - 5/18)) ;
mdb1 = mean(temp_cold(indd));

[value] = wndmag_cold;
ntim = sum(~isnan(value)) ; %length(value) ;
value_sort = sort(value(~isnan(value)), 'descend') ;
threshold = floor(0.01*ntim) ;
value_threshold = value_sort(threshold);
ws2 = value_threshold ;
indd = find((value <= value_threshold + 5/18) & (value >= value_threshold - 5/18)) ;
mdb2 = mean(temp_cold(indd));

%% Deal output variables into out_variable

out_variable = cell(length(out_variable_names), 1) ; 
for ivarn = 1:length(out_variable_names) ; 
    eval(['out_variable{ivarn} = ' out_variable_names{ivarn} ';']) ; 
end
