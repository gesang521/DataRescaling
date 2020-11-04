function [out_variable, out_variable_names] = calculation1(Data_in) ;
% we are calculating Cooling DB/MCWB, Evaporation WB/MCDB corresponding to ...
% 0.4% 1% 2% respectively, and Extreme Max WB in this code .

%% Define output variable names
out_variable_names = {'db', 'mcwb', 'wb', 'mcdb', 'extreme_max_wb'} ; 

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

extreme_max_wb = max(tstar_c) ;

%% Treshold exceedence claculation--Annually
% 0.4% 1% 2%  DB and MCWB
db = zeros(3,1) ; 
mcwb = zeros(3,1) ;
i = 1 ;
for percent = [0.004 0.01 0.02]
    ntim = sum(~isnan(temp)) ;
    temp_sort = sort(temp(~isnan(temp)), 'descend') ;
    threshold = floor(percent*ntim) ;
    temp_threshold = temp_sort(threshold) ;
    db(i) = temp_threshold ; 

    %mean coincident values calculation 
    indd = find((temp <= temp_threshold + 5/18) & (temp >= temp_threshold - 5/18)) ;
    mcwb(i) = mean(tstar_c(indd));
    
    i = i+1 ; 
end 

% 0.4% 1% 2%  WB and MCDB
wb = zeros(3,1) ;
mcdb = zeros(3,1) ;
i = 1; 
for percent = [0.004 0.01 0.02]
    ntim = sum(~isnan(tstar_c)) ;
    tstar_c_sort = sort(tstar_c(~isnan(tstar_c)), 'descend') ;
    threshold = floor(percent*ntim) ;
    tstar_c_threshold = tstar_c_sort(threshold) ;
    wb(i) = tstar_c_threshold ;

    %mean coincident values calculation 
    indd = find((tstar_c <= tstar_c_threshold + 5/18) & (tstar_c >= tstar_c_threshold - 5/18)) ;
    mcdb(i) = mean(temp(indd));
    
    i = i+1 ; 
end 

%% Deal output variables into out_variable

out_variable = cell(length(out_variable_names), 1) ; 
for ivarn = 1:length(out_variable_names) ; 
    eval(['out_variable{ivarn} = ' out_variable_names{ivarn} ';']) ; 
end

