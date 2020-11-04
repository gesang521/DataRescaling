function [out_variable, out_variable_names] = calculation7(Data_in) ;
% we are calculating calclate teh number of hours between 8am to 4pm with temp between 12.8 to 20.6 in this code.

%% Define output variable names
out_variable_names = {'hours8_4'} ;                                            

% You should not need to change anything below this point
% Deal out input variables
in_variable_names = {'yr', 'mo', 'hr', 'day', 'tmp', 'dew', 'wndmag', 'wnddir'} ;
[yr, mo, day, hr] = datevec(Data_in.time-.25) ; 
for ivarn = 5:length(in_variable_names) ; 
    eval([in_variable_names{ivarn} ' = Data_in.' in_variable_names{ivarn} ';']) ; 
end 
temp = tmp ; 

%% calclate teh number of hours between 8am to 4pm with temp between 12.8 to 20.6
% WI tme of 8am tp 4pm is UTC of 1pm to 9pm
ind_t = find(13 <= hr & hr <= 21) ;
ind_tmp = find(12.8<tmp(ind_t) & tmp(ind_t)<20.6) ;
hours8_4 = length(hr(ind_tmp))/24 ; % divded by 24 years, get the value for each year 
%% Deal output variables into out_variable

out_variable = cell(length(out_variable_names), 1) ;
for ivarn = 1:length(out_variable_names) ; 
    eval(['out_variable{ivarn} = ' out_variable_names{ivarn} ';']) ; 
end   