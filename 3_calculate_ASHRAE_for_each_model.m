%% Run through all models to calculate ASHRAE variables
clc;
clear all; 

%% Set up cities and scenarios 
scen_rcp = 'rcp85' ;
scen_rcp_mod = 'rcp85' ;

% city info
CitiesLname = {'MADISON','CHICAGO','ATLANTA','BOSTON', 'DALLAS', 'HOUSTON', 'MIAMI', 'NASHVILLE',... 
    'OMAHA', 'STLOUIS', 'COLUMBUS', 'DENVER','MINNEAPOLIS',...
    'NEWYORK', 'RALEIGH', 'WASHINGTONDC',} ;
CitiesSname = {'MSN','MDW','ATL','BOS','DFW','IAH','MIA','BNA',...
    'OMA','STL','CMH','DNE','MSP',...
    'JFK','RDU','IAD'} ;

% scen
scen85 = {'rcp85_2020','rcp85_2030','rcp85_2040','rcp85_2050','rcp85_2060','rcp85_2070','rcp85_2080','rcp85_2090'}; 
scen45 = {'rcp45_2020','rcp45_2030','rcp45_2040','rcp45_2050','rcp45_2060','rcp45_2070','rcp45_2080','rcp45_2090'}; 

% scen2 will include wet bulb
scen285 = {'rcp85b_2020','rcp85b_2030','rcp85b_2040','rcp85b_2050','rcp85b_2060','rcp85b_2070','rcp85b_2080','rcp85b_2090'}; 
scen245 = {'rcp45b_2020','rcp45b_2030','rcp45b_2040','rcp45b_2050','rcp45b_2060','rcp45b_2070','rcp45b_2080','rcp45b_2090'}; 

yrstr = {'2011_2030','2021_2040','2031_2050','2041_2060','2051_2070','2061_2080','2071_2090','2081_2100',}; 
% yrstr = {'2031_2050'}
% run through all cities
for icity = 15 %length(CitiesLname)
    
    cityl = CitiesLname{icity}
    citys = CitiesSname{icity};

    if strcmp(scen_rcp, 'rcp45')
        scen = scen45 ;
        scen2 = scen245 ;
    else
        scen = scen85 ;
        scen2 = scen285 ;
    end

    %% Set up directory paths, etc.
    % directory paths
    rootdir = '/data/shared/Projects/Gesang' ; 
    datadir = fullfile(rootdir, 'Data') ; 
    routinedir = fullfile(rootdir, 'dan_calculations') ; 

    % make data directory to save ASHRAE variables
    if exist(fullfile(datadir,cityl,'ASHRAE_vars',scen_rcp,'ASHRAE_each_model')) == 0 ; 
        mkdir(fullfile(datadir,cityl,'ASHRAE_vars',scen_rcp,'ASHRAE_each_model')) ; 
    end

    % get calculation routines
    cd(routinedir) ; 
    routine_name = {'calculation', 'calculation1', 'calculation2', ...
         'calculation3', 'calculation4', 'calculation5', 'calculation6','calculation7'} ;
    %routine_name = {'calculation2','calculation4'} ; % routines with DP, MCD, and CDH
    nrout = length(routine_name) ; 

    %% Set up variables

    % Load model names
    load(fullfile(datadir, cityl, scen_rcp_mod, 'Model_Names.mat')) ; %'rcp85_2050' is the only file that include nodel names
    nmod = length(modnames) ; 

    % Set up output data structure
    Data.model = {} ; 
    VarNames = {} ; 

    %% Iterate through each routine, then each model 
    for ideca = 5:length(yrstr)
        totvar = 0 ; 
        for irout = 1:nrout  % loop for calculations

            for imod = 1:nmod ;  % loop for models
            tic

            Data(imod).model = modnames(imod).name ; 
            disp(['Routine: ' routine_name{irout} '  Model: ' Data(imod).model]) ;  
                filin = fullfile(datadir, cityl, scen_rcp, scen2{ideca}, ...
                    ['Hourly_' citys '.' Data(imod).model '.' yrstr{ideca} '.mat']) ;
            load(filin) ; 

            %% Insert routines here
            eval(['[data, varname] = ' routine_name{irout} '(Data_out) ;']) ; 
            nvar = length(varname) ;
                for ivarn = 1:nvar ;
                    eval(['Data(imod).' varname{ivarn} ' = data{ivarn} ;']) ; 
                end
                toc
            end

            for ivarn = 1:nvar ; 
                VarNames{totvar + ivarn} = varname{ivarn} ;
            end
            totvar = totvar + nvar ;    
        end

    %% Save results

        filout = fullfile(datadir, cityl, 'ASHRAE_vars', scen_rcp, 'ASHRAE_each_model', ['ASHRAE_Data.' scen_rcp '.' yrstr{ideca} '.mat']) ; 
        save(filout, 'Data', 'VarNames') ; 
    end
    cd('/home/gesang') ;
end
% exit

