function [out_variable, out_variable_names] = calculation3(Data_in) ;
%% we are calculating extreme annual design conditions including ...
 % extreme annual wind speed, mean/std of max and min of extreme annual DB ...
 % n-year return period values of extreme DB in this code.
 % note: extreme max wet bulb temp will not be included in this code.

%% Define output variable names
out_variable_names = {'extreme_ws', 'mean_min', 'mean_max', 'sd_min', 'sd_max',...
                      'T_n_min', 'T_n_max'} ; 

% You should not need to change anything below this point
% Deal out input variables
in_variable_names = {'yr', 'mo', 'hr', 'day', 'tmp', 'dew', 'wndmag', 'wnddir'} ;
[yr, mo, day, hr] = datevec(Data_in.time-.25) ; 
for ivarn = 5:length(in_variable_names) ; 
    eval([in_variable_names{ivarn} ' = Data_in.' in_variable_names{ivarn} ';']) ; 
end 
temp = tmp ; 

%% Calculation: 
% Extreme annual wind speed 
extreme_ws = zeros(3,1) ; 
i = 1 ; 
for precent = [0.01 0.025 0.05]
    [value] = wndmag ;
    ntim = sum(~isnan(value)) ;
    value_sort = sort(value(~isnan(value)), 'descend') ;
    threshold = floor(precent*ntim) ;
    value_threshold = value_sort(threshold);
    extreme_ws(i) = value_threshold ;
    
    i = i+1 ; 
end
%extreme_ws

% Extrem Annual Design Conditions 
extreme = zeros(25,1) ; 
tmax = zeros(25,1) ;
tmin = zeros(25,1) ;
i = 1 ; 
for iyr = [1986:2010] 
        ind = find(yr==iyr);
        t = temp(ind) ;
        tmax(i) = max(t) ;
        tmin(i) = min(t) ;
        i = i+1 ; 
end 

mean_min = mean(tmin) ;
mean_max = mean(tmax) ;
sd_min = std(tmin) ;
sd_max = std(tmax) ;

% n-year return period values of extreme DB
T_n_min = zeros(4,1) ; 
i = 1 ; 
for n = [5 10 20 50]
    I = -1 ; % i for max and -1 for min
    M = mean(tmin) ;  % mean of max or min DB 
    s = std(tmin) ;  %std of max or min DB 
    ff = 0.5772 + log(log(n/(n-1))) ; 
    F = (-sqrt(6)/pi)*ff ;
    T_n_min(i) = M +I*F*s ;
    
    i = i+1 ;    
end 
%T_n_min

T_n_max = zeros(4,1) ; 
i = 1 ; 
for n = [5 10 20 50]
    I = 1 ; % i for max and -1 for min
    M = mean(tmax) ;  % mean of max or min DB 
    s = std(tmax) ;  %std of max or min DB 
    ff = 0.5772 + log(log(n/(n-1))) ; 
    F = (-sqrt(6)/pi)*ff ;
    T_n_max(i) = M +I*F*s ;
    
    i = i+1 ;    
end 
%T_n_max

%% Deal output variables into out_variable

out_variable = cell(length(out_variable_names), 1) ; 
for ivarn = 1:length(out_variable_names) ; 
    eval(['out_variable{ivarn} = ' out_variable_names{ivarn} ';']) ; 
end
