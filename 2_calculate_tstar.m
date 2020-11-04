% this code is used to calculate wet bulb temp. Wet bulb temp will be
% calculated for each cities, each year under each scenario
% when this code is used for calculatio in other cities, need to change 
% city parameter at the begining and also the height z.
clear all

%% Set up parameters, directory paths, etc.

% choose rcp scenario 
scen_rcp = 'rcp45' ; 
scen_rcp_mod = 'rcp45' ;

% city 
% the order of cities cannot be changed since they are with the same index
% of latlon and elevation
CitiesLname = {'MADISON','CHICAGO','ATLANTA','BOSTON', 'DALLAS', 'HOUSTON', 'MIAMI', 'NASHVILLE',... 
    'OMAHA', 'STLOUIS', 'COLUMBUS', 'DENVER','MINNEAPOLIS',...
    'NEWYORK', 'RALEIGH', 'WASHINGTONDC',} ;
CitiesSname = {'MSN','MDW','ATL','BOS','DFW','IAH','MIA','BNA',...
    'OMA','STL','CMH','DNE','MSP',...
    'JFK','RDU','IAD'} ;

cityl = 'ATLANTA' ;
citys = 'ATL' ;

ind1 = find(strcmp(CitiesLname, cityl)) ;
ind2 = find(strcmp(CitiesSname, citys)) ;

% get elevation 
fgeo = '/data/shared/Projects/Gesang/Data/geoinfo.mat' ;
load(fgeo) ;
ind = find(strcmp(CitiesLname, cityl)) ;
elevation = Geoinfo(ind).elv ;
elevation = elevation*3.28084 ; % meter to ft

% scenario
rootdir = '/data/shared/Projects/Gesang' ; 
datadir = fullfile(rootdir, 'Data') ; 
routinedir = fullfile(rootdir, 'dan_calculations') ; 

scen85 = {'rcp85_2020','rcp85_2030','rcp85_2040','rcp85_2050','rcp85_2060','rcp85_2070','rcp85_2080','rcp85_2090'}; 
scen45 = {'rcp45_2020','rcp45_2030','rcp45_2040','rcp45_2050','rcp45_2060','rcp45_2070','rcp45_2080','rcp45_2090'}; 

% scen2 will include wet bulb
scen285 = {'rcp85b_2020','rcp85b_2030','rcp85b_2040','rcp85b_2050','rcp85b_2060','rcp85b_2070','rcp85b_2080','rcp85b_2090'}; 
scen245 = {'rcp45b_2020','rcp45b_2030','rcp45b_2040','rcp45b_2050','rcp45b_2060','rcp45b_2070','rcp45b_2080','rcp45b_2090'}; 

if strcmp(scen_rcp, 'rcp45')
    scen = scen45 ;
    scen2 = scen245 ;
else
    scen = scen85 ;
    scen2 = scen285 ;
end

% make data directory
for ideca = 1:length(scen)
    if exist(fullfile(datadir,cityl,scen_rcp,scen2{ideca})) == 0 ; 
        mkdir(fullfile(datadir,cityl,scen_rcp,scen2{ideca})) ; 
    end
end

% year
yrstr = {'2011_2030','2021_2040','2031_2050','2041_2060','2051_2070','2061_2080','2071_2090','2081_2100'} ;

%% Set up variables

% Load model names
load(fullfile(datadir, cityl, scen_rcp_mod, 'Model_Names.mat')) ; 
nmod = length(modnames) ; 

% Set up output data structure
Data.model = {} ; 
VarNames = {} ; 

%% Iterate through each routine, then each model
tic
for ideca = 1:2%length(yrstr)
    
    for imod = 1:nmod ; 
        Data(imod).model = modnames(imod).name ; 
        disp(Data(imod).model)
        filin = fullfile(datadir, cityl, scen_rcp, scen{ideca},...
            ['Hourly_' citys '.' Data(imod).model '.' yrstr{ideca} '.mat']) ; 
        filout = fullfile(datadir, cityl, scen_rcp, scen2{ideca}, ...
            ['Hourly_' citys '.' Data(imod).model '.' yrstr{ideca} '.mat']) ;
        load(filin) ; 

        % Deal  input input variables
        in_variable_names = {'yr', 'mo', 'hr', 'day', 'tmp', 'dew', 'wndmag', 'wnddir'} ;
        [yr, mo, day, hr] = datevec(Data_out.time-.25) ; 
        for ivarn = 5:length(in_variable_names) ; 
            eval([in_variable_names{ivarn} ' = Data_out.' in_variable_names{ivarn} ';']) ; 
        end 

        temp = tmp ; 

        z = elevation ;
        p = 14.696*(1-6.8754*10^(-6)*z)^5.2559 ;%narometroc pressure psia
        t_f = temp*9/5+32 ; %celsius to fahrenhelt 
        dew_f = dew*9/5+32 ;

        % Saturation pressure calculation pw
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

        % Humidity ratio calculation  wout
        wout = 0.621945*(pw./((p*1013/1013)-pw)) ;


        %% Calculation wet bulb temp
        tstar2 = nan(size(wout)) ; 

        % Define function for t_f >= 32 ; 
        pws1 = @(tstar) exp(c8/tstar + c9 + c10*tstar + c11*tstar.^2 + c12*tstar.^3 + c13*log(tstar)) ; 
        ws1 = @(tstar,p) 0.621945 * pws1(tstar) ./ (p - pws1(tstar)) ; 
        w1 = @(tstar,t_f,p) ((1093 - 0.556.*(tstar-459.67)).*ws1(tstar,p) - 0.240*(t_f - (tstar - 459.67))) ./ (1093 + .444*t_f - (tstar - 459.67)) ; 
        % Define function for t_f < 32 ; 
        pws2 = @(tstar) exp(c1/tstar + c2 + c3*tstar + c4*tstar.^2 + c5*tstar.^3 + c6*tstar.^4 + c7*log(tstar)) ; 
        ws2 = @(tstar,p) 0.621945 * pws2(tstar) ./ (p - pws2(tstar)) ; 
        w2 = @(tstar,t_f,p) ((1220 - 0.04.*(tstar-459.67)).*ws2(tstar,p) - 0.240*(t_f - (tstar - 459.67))) ./ (1220 + .444*t_f - 0.48*(tstar - 459.67)) ; 

        ind = find(~isnan(wout)) ; 
        for iind = 1:length(ind)
            i = ind(iind) ;
            if t_f(i) >= 32 
                tstar2(i) = fzero(@(tstar) w1(tstar, t_f(i), p) - wout(i), t_f(i)+459.67) ; 
            else
                tstar2(i) = fzero(@(tstar) w2(tstar, t_f(i), p) - wout(i), t_f(i)+459.67) ;    
            end
        end

        Data_out.tstar2 = tstar2 ; 
        save(filout, 'Data_out') ; 
    end
    
end
toc

%% ***********************************************
clc
clear all

%% Set up parameters, directory paths, etc.

% choose rcp scenario 
scen_rcp = 'rcp45' ; 
scen_rcp_mod = 'rcp45' ;

% city 
% the order of cities cannot be changed since they are with the same index
% of latlon and elevation
CitiesLname = {'MADISON','CHICAGO','ATLANTA','BOSTON', 'DALLAS', 'HOUSTON', 'MIAMI', 'NASHVILLE',... 
    'OMAHA', 'STLOUIS', 'COLUMBUS', 'DENVER','MINNEAPOLIS',...
    'NEWYORK', 'RALEIGH', 'WASHINGTONDC',} ;
CitiesSname = {'MSN','MDW','ATL','BOS','DFW','IAH','MIA','BNA',...
    'OMA','STL','CMH','DNE','MSP',...
    'JFK','RDU','IAD'} ;

cityl = 'ATLANTA' ;
citys = 'ATL' ;

ind1 = find(strcmp(CitiesLname, cityl)) ;
ind2 = find(strcmp(CitiesSname, citys)) ;

% get elevation 
fgeo = '/data/shared/Projects/Gesang/Data/geoinfo.mat' ;
load(fgeo) ;
ind = find(strcmp(CitiesLname, cityl)) ;
elevation = Geoinfo(ind).elv ;
elevation = elevation*3.28084 ; % meter to ft

% scenario
rootdir = '/data/shared/Projects/Gesang' ; 
datadir = fullfile(rootdir, 'Data') ; 
routinedir = fullfile(rootdir, 'dan_calculations') ; 

scen85 = {'rcp85_2020','rcp85_2030','rcp85_2040','rcp85_2050','rcp85_2060','rcp85_2070','rcp85_2080','rcp85_2090'}; 
scen45 = {'rcp45_2020','rcp45_2030','rcp45_2040','rcp45_2050','rcp45_2060','rcp45_2070','rcp45_2080','rcp45_2090'}; 

% scen2 will include wet bulb
scen285 = {'rcp85b_2020','rcp85b_2030','rcp85b_2040','rcp85b_2050','rcp85b_2060','rcp85b_2070','rcp85b_2080','rcp85b_2090'}; 
scen245 = {'rcp45b_2020','rcp45b_2030','rcp45b_2040','rcp45b_2050','rcp45b_2060','rcp45b_2070','rcp45b_2080','rcp45b_2090'}; 

if strcmp(scen_rcp, 'rcp45')
    scen = scen45 ;
    scen2 = scen245 ;
else
    scen = scen85 ;
    scen2 = scen285 ;
end
% year
yrstr = {'2011_2030','2021_2040','2031_2050','2041_2060','2051_2070','2061_2080','2071_2090','2081_2100'} ;

%% Set up variables

% Load model names
load(fullfile(datadir, cityl, scen_rcp_mod, 'Model_Names.mat')) ; 
nmod = length(modnames) ; 

% Set up output data structure
Data.model = {} ; 
VarNames = {} ; 

%% Iterate through each routine, then each model
tic
for ideca = 1:2%length(yrstr)
    
    for imod = 1:nmod ; 
        Data(imod).model = modnames(imod).name ; 
        disp(Data(imod).model)
        filin = fullfile(datadir, cityl, scen_rcp, scen{ideca},...
            ['Hourly_' citys '.' Data(imod).model '.' yrstr{ideca} '.mat']) ; 
        filout = fullfile(datadir, cityl, scen_rcp, scen2{ideca}, ...
            ['Hourly_' citys '.' Data(imod).model '.' yrstr{ideca} '.mat']) ;
        load(filin) ; 

        % Deal  input input variables
        in_variable_names = {'yr', 'mo', 'hr', 'day', 'tmp', 'dew', 'wndmag', 'wnddir'} ;
        [yr, mo, day, hr] = datevec(Data_out.time-.25) ; 
        for ivarn = 5:length(in_variable_names) ; 
            eval([in_variable_names{ivarn} ' = Data_out.' in_variable_names{ivarn} ';']) ; 
        end 

        temp = tmp ; 

        z = elevation ;
        p = 14.696*(1-6.8754*10^(-6)*z)^5.2559 ;%narometroc pressure psia
        t_f = temp*9/5+32 ; %celsius to fahrenhelt 
        dew_f = dew*9/5+32 ;

        % Saturation pressure calculation pw
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

        % Humidity ratio calculation  wout
        wout = 0.621945*(pw./((p*1013/1013)-pw)) ;


        %% Calculation wet bulb temp
        tstar2 = nan(size(wout)) ; 

        % Define function for t_f >= 32 ; 
        pws1 = @(tstar) exp(c8/tstar + c9 + c10*tstar + c11*tstar.^2 + c12*tstar.^3 + c13*log(tstar)) ; 
        ws1 = @(tstar,p) 0.621945 * pws1(tstar) ./ (p - pws1(tstar)) ; 
        w1 = @(tstar,t_f,p) ((1093 - 0.556.*(tstar-459.67)).*ws1(tstar,p) - 0.240*(t_f - (tstar - 459.67))) ./ (1093 + .444*t_f - (tstar - 459.67)) ; 
        % Define function for t_f < 32 ; 
        pws2 = @(tstar) exp(c1/tstar + c2 + c3*tstar + c4*tstar.^2 + c5*tstar.^3 + c6*tstar.^4 + c7*log(tstar)) ; 
        ws2 = @(tstar,p) 0.621945 * pws2(tstar) ./ (p - pws2(tstar)) ; 
        w2 = @(tstar,t_f,p) ((1220 - 0.04.*(tstar-459.67)).*ws2(tstar,p) - 0.240*(t_f - (tstar - 459.67))) ./ (1220 + .444*t_f - 0.48*(tstar - 459.67)) ; 

        ind = find(~isnan(wout)) ; 
        for iind = 1:length(ind)
            i = ind(iind) ;
            if t_f(i) >= 32 
                tstar2(i) = fzero(@(tstar) w1(tstar, t_f(i), p) - wout(i), t_f(i)+459.67) ; 
            else
                tstar2(i) = fzero(@(tstar) w2(tstar, t_f(i), p) - wout(i), t_f(i)+459.67) ;    
            end
        end

        Data_out.tstar2 = tstar2 ; 
        save(filout, 'Data_out') ; 
    end
    
end
toc
