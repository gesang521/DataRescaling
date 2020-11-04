%Rescale hourly data to future conditions
clc ;
clear all ;

%% Set parameters
% CDF files
rootdir = '/data1/djlorenz/downscaling3/final/cmip5' ; 
future_yr = {'2011_2030','2021_2040','2031_2050','2041_2060','2051_2070','2061_2080','2071_2090','2081_2100'} ; 

% Future scenarios
scen_rcp = 'rcp85' ;
scen_rcp_mod = 'rcp85' ;

scen85 = {'rcp85_2020','rcp85_2030','rcp85_2040','rcp85_2050','rcp85_2060','rcp85_2070','rcp85_2080','rcp85_2090'}; 
scen45 = {'rcp45_2020','rcp45_2030','rcp45_2040','rcp45_2050','rcp45_2060','rcp45_2070','rcp45_2080','rcp45_2090'}; 

if strcmp(scen_rcp, 'rcp45')
    scen = scen45 ;
else
    scen = scen85 ;
end

% cities 
CitiesLname = {'MADISON','CHICAGO','ATLANTA','BOSTON', 'DALLAS', 'HOUSTON', 'MIAMI', 'NASHVILLE',... 
    'OMAHA', 'STLOUIS', 'COLUMBUS', 'DENVER','MINNEAPOLIS',...
    'NEWYORK', 'RALEIGH', 'WASHINGTONDC'} ;
CitiesSname = {'MSN','MDW','ATL','BOS','DFW','IAH','MIA','BNA',...
    'OMA','STL','CMH','DNE','MSP',...
    'JFK','RDU','IAD'} ;
    
cityl = 'DENVER' ;
citys = 'DNE' ;

disp(cityl)

% get city latlon
fgeo = '/data/shared/Projects/Gesang/Data/geoinfo.mat' ; 
load(fgeo) ;
ind1 = find(strcmp(CitiesLname, cityl)) ;
ind2 = find(strcmp(CitiesSname, citys)) ;
Target_lat = Geoinfo(ind1).lat ; 
Target_lon = Geoinfo(ind1).lon ; 

%% Load data

% Historical data
datadir = '/data/shared/Projects/Gesang/Data/' ;
load(fullfile(datadir,cityl,['Hourly_' CitiesSname{ind2} '_Data.mat'])) ;
Data_obs = Data_out ;

% Make data directory
if exist(fullfile(datadir,cityl,scen_rcp)) == 0 ; 
    mkdir(fullfile(datadir,cityl,scen_rcp)) ; 
end

for ideca = 1:length(scen)
    if exist(fullfile(datadir,cityl,scen_rcp,scen{ideca})) == 0 ; 
        mkdir(fullfile(datadir,cityl,scen_rcp,scen{ideca})) ; 
    end
end

% We have one extra leap day, because there is no 1996

% Load model names - historical models ?1950-2005?
modname1 = dir(fullfile(rootdir, 'historical')) ; 
kp = [] ; 
for i = 1:length(modname1) ; 
    if ((modname1(i).isdir == 1) & ~any(strcmp(modname1(i).name, {'.', '..'}))) ;
        kp = [kp; i] ;
    end
end
modname1 = modname1(kp) ; 

% Load model names - future models ?2006-2100?
modname2 = dir(fullfile(rootdir, scen_rcp_mod)) ; 

% Get common models
kp = [] ; 
for i1 = 1:length(modname1) ; 
    for i2 = 1:length(modname2) ; 
        if strcmp(modname1(i1).name, modname2(i2).name)
            kp = [kp i1] ; 
        end
    end
end

% Retain models that are incommon
modnames = modname1(kp) ; 
nmod = length(modnames) ; 

%% Now, start calculations. 
% Start by getting tmax and tmin for each day

time = Data_obs.time ;
[yr, mon, day] = datevec(time) ; 

tmp = Data_obs.tmp ;
temp = reshape(tmp(7:(end-18)), [24 length(tmp)/24-1]) ; % 24hrs, total days in 25years 
%convert time form UTC to local time, and reshape hourly data to daily data 
tmax = max(temp) ; % daily max temp
tmin = min(temp) ; 

% Replicate the first and last day - to get end points
tmax = [tmax(1) tmax tmax(end)] ; 
tmin = [tmin(1) tmin tmin(end)] ; 
temp = [[nan(18,1); tmp(1:6)] temp [tmp(end+[-17:0]); nan(6, 1)]] ;

% get daily temporal scaling 
a1 = nan(size(temp)) ; 
for ihr = 1:24 ; 
    a1(ihr,:) = (temp(ihr,:) - tmin)./(tmax - tmin) ; 
end

% Get temporal scaling for CDFs. Note that:
% cdf_use = b1(i)*cdf(cdfmon(i)) + (1-b1(i))*cdfmon(b2(i)) ;
cdfyr = [yr(1)-1 yr(1:24:end)'] ; 
cdfday = [31 day(1:24:end)'] ;
cdfmon = [12 mon(1:24:end)'] ; 
b1 = nan(size(cdfyr)) ;
b2 = nan(size(cdfyr)) ; 
for i = 1:length(cdfyr) ; 
    dpm = [31 28+(~mod(cdfyr(i), 4)) 31 30 31 30 31 31 30 31 30 31] ; 
    midday = (dpm(cdfmon(i))+1)/2 ; 
    d2midday = cdfday(i) - midday ;
    b2(i) = mod( (cdfmon(i) - 1 + sign(d2midday)), 12) + 1 ; 
    dp = (dpm(b2(i)) + dpm(cdfmon(i)))/2 ; 
    b1(i) = 1 - abs(d2midday)/dp ;
end

%% Get CDF data

% Get lat / lon indices for loading
filin = fullfile(rootdir, 'historical', modnames(1).name, ...
    'r1i1p1', ['temp_cdf_1981_2010.nc']) ; 
lat = ncread(filin, 'lat') ;
lon = ncread(filin, 'lon') ; 
lev = ncread(filin, 'lev') ; 

[tem, ykp] = min(abs(lat - Target_lat)) ;
[tem, xkp] = min(abs(lon - Target_lon)) ; 
nlev = length(lev) ; 

start = [1 1 ykp xkp] ;%- 1 ; 
count = [12 nlev 1 1] ; 
start = fliplr(start) ; 
count = fliplr(count) ; 
%% Now, loop through the models
for ideca = 1:length(future_yr)
    for imod = 1:nmod

        disp(modnames(imod).name) ;

        % Initialize output data
        Data_out = Data_obs ; 

        % Get model names and data
        filin1 = fullfile(rootdir, 'historical', modnames(imod).name, ...
        'r1i1p1', ['temp_cdf_1981_2010.nc']) ; 
        filin2 = fullfile(rootdir, scen_rcp_mod, modnames(imod).name, ...
        'r1i1p1', ['temp_cdf_' future_yr{ideca} '.nc']) ; 

        tmax_cdf1 = ncread(filin1, 'cdftmax', start, count) ; 
        tmin_cdf1 = ncread(filin1, 'cdftmin', start, count) ; 
        tmax_cdf2 = ncread(filin2, 'cdftmax', start, count) ; 
        tmin_cdf2 = ncread(filin2, 'cdftmin', start, count) ; 
        tmax_cdf1 = squeeze(tmax_cdf1)' ; 
        tmax_cdf2 = squeeze(tmax_cdf2)' ; 
        tmin_cdf1 = squeeze(tmin_cdf1)' ; 
        tmin_cdf2 = squeeze(tmin_cdf2)' ; 


        tmin2 = nan(size(tmin)) ; 
        tmax2 = nan(size(tmax)) ; 

        % Scroll through dates
        for i = 1:length(tmax) ; 
            cdf_tmax1 = b1(i)*tmax_cdf1(cdfmon(i),:) + (1-b1(i))*tmax_cdf1(b2(i),:) ;
            cdf_tmin1 = b1(i)*tmin_cdf1(cdfmon(i),:) + (1-b1(i))*tmin_cdf1(b2(i),:) ;
            cdf_tmax2 = b1(i)*tmax_cdf2(cdfmon(i),:) + (1-b1(i))*tmax_cdf2(b2(i),:) ;
            cdf_tmin2 = b1(i)*tmin_cdf2(cdfmon(i),:) + (1-b1(i))*tmin_cdf2(b2(i),:) ;

            % Rescale tmax and tmin             
            ind = max(find(lev <= tmax(i))) ;
            pind = cdf_tmax1(ind) + (tmax(i)-lev(ind))*diff(cdf_tmax1(ind+[0 1]))/diff(lev(ind+[0 1])) ;
            ind = max(find(cdf_tmax2 <= pind)) ;
            tmax2(i) = lev(ind) + (pind - cdf_tmax2(ind))*diff(lev(ind+[0 1]))/diff(cdf_tmax2(ind+[0 1])) ;

            ind = max(find(lev <= tmin(i))) ; 
            pind = cdf_tmin1(ind) + (tmin(i)-lev(ind))*diff(cdf_tmin1(ind+[0 1]))/diff(lev(ind+[0 1])) ;
            ind = max(find(cdf_tmin2 <= pind)) ;
            tmin2(i) = lev(ind) + (pind - cdf_tmin2(ind))*diff(lev(ind+[0 1]))/diff(cdf_tmin2(ind+[0 1])) ;
        end

        % Scale temperature data
        temp2 = nan(size(temp)) ; 
        for ihr = 1:24 ;
            % Remember, t = (1 - a1)*tmin + a1*tmax ; 
            temp2(ihr,:) = (1 - a1(ihr,:)).*tmin2 + a1(ihr,:).*tmax2 ; 
        end
        tmp2 = reshape(temp2, length(tmp)+24, 1) ; 
        tmp2 = tmp2(19:(end - 6)) ; 

        % Dump into Data_out
        Data_out.tmp = tmp2 ; 

        % Dewpoint Temperature
        E0 = 0.611 ; 
        LdRv = 5423 ; 
        E = E0*exp(LdRv*(1/273 - 1./(Data_obs.dew+273.15))) ;
        Es = E0*exp(LdRv*(1/273 - 1./(Data_obs.tmp+273.15))) ; 
        RH = E./Es ; 

        Es2 = E0*exp(LdRv*(1/273 - 1./(Data_out.tmp + 273.15))) ; 
        E2 = RH.*Es2 ; 
        Data_out.dew = ((1/273.150) - log(E2/E0)/LdRv).^-1 - 273.15; 

        % Save results scen 

        filout = fullfile(datadir,cityl,scen_rcp,scen{ideca}, ...
            ['Hourly_' CitiesSname{ind2} '.' modnames(imod).name '.' future_yr{ideca} '.mat']) ; 
        save(filout, 'Data_out') ; 
    end  
end

filout2 = fullfile(datadir,cityl,scen_rcp,'Model_Names.mat') ; 
save(filout2, 'modnames') ; 

