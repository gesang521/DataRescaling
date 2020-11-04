function [out_variable, out_variable_names] = calculation2(Data_in) ;
% we are calculating Humidification DP/MCDB and HR... 
% Dehumidification DP/MCDB and HR,corresponding to 99.6% 99%...
% 0.4% 1% 2% and Enthalpy/MCDB corresponding to 0.4% 1% 2%  in this code.

%% Define output variable names
out_variable_names = {'dp', 'hr', 'mcdb_dp', 'enthalpy', 'mcdb_enth'} ; 

% You should not need to change anything below this point
% Deal out input variables
in_variable_names = {'yr', 'mo', 'hr', 'day', 'tmp', 'dew', 'wndmag', 'wnddir', 'elv'} ;
[yr, mo, day, hr] = datevec(Data_in.time-.25) ; 
for ivarn = 5:length(in_variable_names) ; 
    eval([in_variable_names{ivarn} ' = Data_in.' in_variable_names{ivarn} ';']) ; 
end 
temp = tmp ; 
z = elv(1)
%% Parameter setting 
% z = 886.142 ;% unit:ft   264meters
p = 14.696*(1-6.8754*10^(-6)*z)^5.2559 ; %narometroc pressure psia
t_f = temp*9/5+32 ; %celsius to fahrenhelt 
dew_f = dew*9/5+32 ;

%% Saturation pressure calculation pw
% defining constant 
c1 = -1.0214165*10^4;
c2 = -4.8932428*10^0;
c3 = -5.3765794*10^(-3);
c4 = 1.9202377*10^(-7);
c5 = 3.5575832*10^(-10);
c6 = -9.0344688*10^(-14);
c7 = 4.1635019*10^(0);

c8 = -1.0440397*10^4;
c9 = -1.1294650*10^1;
c10 = -2.7022355*10^(-2);
c11 = 1.2890360*10^(-5);
c12 = -2.4780681*10^(-9);
c13 = 6.5459673*10^0;

%getting absolute temp by adding 456.67
t = t_f + 459.67 ;   % t: absolute temperature
td = dew_f + 459.67 ;  % t: absolute dew point
v = nan(length(temp),1) ; 

%calculating saturation pressure pws
ind_noNaN = find(~isnan(td) & ~isnan(t_f)) ;
ind = find(t_f(ind_noNaN) < 32) ; 
i = ind_noNaN(ind) ; 
v(i) = c1./td(i) + c2 + c3*td(i) + c4*td(i).^2 + c5*td(i).^3 + c6*td(i).^4 + c7*log(td(i));

ind = find(t_f(ind_noNaN) >= 32) ; 
i = ind_noNaN(ind) ; 
v(i) = c8./td(i) + c9 + c10*td(i) + c11*td(i).^2 + c12*td(i).^3 + c13*log(td(i));

pw = exp(v) ;
% Humidity ratio calculation  w
w = 0.621945*(pw./(p-pw)) ;

% Enthalpy for moisture air h
hd = 0.240.*(t_f-32);       % specific enthalpy for dry air
hg = 1061+0.444.*(t_f-32);  % specific enthalpy for saturated water vapor
h = hd + w.*hg ;
  
%% Threshold exceedence 
% Dehumidification DP/MCDB and HR
dp = zeros(5,1) ; 
hr = zeros(5,1) ; 
mcdb_dp = zeros(5,1); 
i = 1 ; 
for percent = [0.996 0.99 0.004 0.01 0.02]
    
    ntim = sum(~isnan(dew)) ;
    dew_sort = sort(dew(~isnan(dew)), 'descend') ;
    threshold = floor(percent*ntim) ;
    dew_threshold = dew_sort(threshold) ;
    dp(i) = dew_threshold ;
    % mean coincident values
    ind = find((dew <= dew_threshold + 5/18) & (dew >= dew_threshold - 5/18)) ;
    mcdb_dp(i) = mean(temp(ind)) ;
    
    [value] = w*1000 ;  % convert grains of moisture per lb of dry air to g/kg
    ntim = sum(~isnan(value)) ;
    value_sort = sort(value(~isnan(value)), 'descend') ;
    threshold = floor(percent*ntim) ;
    value_threshold = value_sort(threshold) ;
    hr(i) = value_threshold ;
    i = i+1 ; 
end
   
% Enthalpy/MCDB
enthalpy = zeros(3,1) ; 
mcdb_enth = zeros(3,1); 
i = 1 ; 
for percent = [0.004 0.01 0.02]
    [enth] = h.*1.05506./0.453592 ; %converting ttu/lb to kj/kg
    ntim = sum(~isnan(enth)) ;
    enth_sort = sort(enth(~isnan(enth)), 'descend') ;
    threshold = floor(percent*ntim) ;
    enth_threshold = enth_sort(threshold) ;
    enthalpy(i) = enth_threshold ;
    % mean coincident values
    ind = find((enth <= enth_threshold + 5/18) & (enth >= enth_threshold - 5/18)) ;
    mcdb_enth(i) = mean(temp(ind)) ;

    i = i+1 ; 
end

%% Deal output variables into out_variable

out_variable = cell(length(out_variable_names), 1) ; 
for ivarn = 1:length(out_variable_names) ; 
    eval(['out_variable{ivarn} = ' out_variable_names{ivarn} ';']) ; 
end
