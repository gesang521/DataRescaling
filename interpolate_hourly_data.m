%% Get hourly data into a format we can use
clc
clear all

% 'MADISON','CHICAGO','ATLANTA','MSN','MDW','ATL',
CitiesLname = {'MADISON','CHICAGO','ATLANTA','BOSTON', 'DALLAS', 'HOUSTON', 'MIAMI', 'NASHVILLE',... 
    'OMAHA', 'STLOUIS', 'COLUMBUS', 'DENVER', 'MINNEAPOLIS',...
    'NEWYORK', 'RALEIGH', 'WASHINGTONDC'}
CitiesSname = {'MSN','MDW','ATL','BOS','DFW','IAH','MIA','BNA',...
    'OMA','STL','CMH','DNE','MSP',...
    'JFK','RDU','IAD'}


%% Load data

dirin = '/data/shared/Projects/Gesang/Data' ;
%dirin = '/Users/dvimont/Documents/Research/Gesang/Data' ;

% This file contains full obs for the station, read in and saved as a .mat
% file.
for icity = 12%:length(CitiesLname)

    filin = fullfile( dirin,CitiesLname{icity},'Full_Data_1980_2010.mat' ) ;
    load(filin) ; 

    % We only want to process a subset of the data, so yr1 and yr2 are the
    % first year and last year that we will process (e.g. 1986-2010)
    
    if strcmp(CitiesLname{icity},'DENVER')  % DENVER'data start from 1994-07-18, we rescale form 1996 
      yr_original_data = 1996:2010 ; 
      yr1 = 1996; % Yr1 of data to keep for analysis
    else
      yr_original_data = 1980:2010 ; 
      yr1 = 1986 ; % Yr1 of data to keep for analysis
    end

    scale = [1 10 1 10 10 10] ;
%     yr_original_data = 1980:2010 ; 
%     yr1 = 1986 ;
    yr2 = 2010 ; % Yr2 of data to keep for analysis
    %yr_missing = [1994] ; % Missing years 
    yr_vec = yr1:yr2 ; % Define years that we will be using
    %yr_vec(ismember(yr_vec, yr_missing)) = [] ; % Get rid of missing years
    nyr = length(yr_vec) ;
    yr_ind = find(ismember(yr_original_data, yr_vec)) ; % Get indices from original data for years to keep

    %% Go through each year at a time

    % Define the total number of HOURS of data, so we can initialize variables
    nhrtot = (365*nyr + sum(~mod(yr_vec, 4)))*24 ; % Total number of hours - keep the 24hr at the end
    var_out = {'time', 'tmp', 'wnddir', 'wndmag', 'dew','slp'} ; 
    var_flag = {'', 'tmpflag', 'wnddirflag', 'wndmagflag', 'dewflag','slpflag'} ;
    missing_value = [0 9999 999 9999 9999 99999] ; 
    for i = 1:length(var_out) ; 
        eval(['Data_out.' var_out{i} ' = nan(nhrtot, 1) ; ']) ; 
    end

    yrind = 0 ; 
    for iyr = 1:nyr

        year = yr_vec(iyr) ; 
        disp(['Year = ' num2str(year)]) ; 

        % Define hourly time indices
        dpm = [31 28+(mod(year, 4) == 0) 31 30 31 30 31 31 30 31 30 31] ; 
        ntim_out = sum(dpm)*24 ; 
        yrind = yrind(end) + [1:ntim_out]' ; 
        timout = nan(ntim_out, 1) ; 
        ind = 0 ; 
        for imon = 1:12 ;
            nhr = dpm(imon)*24 ; 
            ind = ind(end) + [1:nhr] ;
            timout(ind) = datenum( ...
                repmat(yr_vec(iyr), [nhr, 1]), ... 
                repmat(imon, [nhr, 1]), ...
                reshape(repmat([1:dpm(imon)], [24, 1]), 24*dpm(imon), 1), ...
                repmat([0:23]', [dpm(imon), 1]), ...
                0, 0) ; 
        end

        % Get time indices from Data
        timdat = datenum(Data(yr_ind(iyr)).yr, Data(yr_ind(iyr)).mon, ...
            Data(yr_ind(iyr)).day, Data(yr_ind(iyr)).hr, Data(yr_ind(iyr)).min, 0) ; 

        % Get duplicate time indices
        [c, ia, ic] = unique(timdat) ; 
        duplicate_time = find(diff(ic) == 0) ; 

        % Scroll through variables 
        Data_out.time(yrind) = timout ; 
        for ivarn = 2:5 ; 
            eval(['tem = Data(yr_ind(iyr)).' var_out{ivarn} ';']) ;
            eval(['flag = Data(yr_ind(iyr)).' var_flag{ivarn} ';']) ; 
            tem(ismember(flag, double(['2367']'))) = NaN ; % Get rid of suspect or erroneous data
            tem(tem == missing_value(ivarn)) = NaN ; % Get rid of missing values

            % Replace duplicate times with average of duplicate times
            if ~isempty(duplicate_time) ; 
                for idup = 1:length(duplicate_time) ;
                    duplicate_ind = find(timdat == c(ic(duplicate_time(idup)))) ; 
                    tem(duplicate_ind(1)) = nanmean(tem(duplicate_ind)) ; 
                    tem(duplicate_ind(2:end)) = NaN ; 
                end ; 
            end 

            % Interpolate
            ind = find(~isnan(tem)) ; 
            tem2 = interp1(timdat(ind), tem(ind), timout) / scale(ivarn);

            % Dump data into Data_out
            eval(['Data_out.' var_out{ivarn} '(yrind) = tem2 ;']) ; 

        end

    end

    %% Save data
    Description = [ ...
        'Data_out is a structure that contains the following: '
        'Data_out.time = time (try datestr(Data_out.time(1)) )'
        'Data_out.tmp = temperature                           '
        'Data_out.wnddir = wind direction (between 10 and 360)'
        'Data_out.wndmag = wind magnitude                     '
        'Data_out.dew = dewpoint temp                         '
        'Data_out.slp = SLP                                   '] ;
    Data_out.Description = Description ; 

    filout = fullfile(dirin,CitiesLname{icity},['Hourly_' CitiesSname{icity} '_Data.mat']) ; 
    save(filout, 'Data_out') ; 

end
%% Example code for loading the data
% 
% if 0 ; 
% 
% clear
% 
 load 'Data/MADISON/Hourly_MSN_Data2.mat'
 var_out = {'time', 'tmp', 'wnddir', 'wndmag', 'dew', 'slp'} ; 
 for i = 1:length(var_out) ; 
     eval([var_out{i} ' = Data_out.' var_out{i} ';']) ; 
 end
 % Get year, month, day, etc
 [yr, mo, day] = datevec(time) ; 


 
% 
