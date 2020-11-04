%% Get average results from recalculations
clc ;
clear all;

%% Set up directory paths, etc.
cityl = 'RALEIGH' ;
citys = 'RDU' ;

yrstr = {'2011_2030','2021_2040','2031_2050','2041_2060','2051_2070','2061_2080','2071_2090','2081_2100',}; 
rootdir = '/data/shared/Projects/Gesang' ; 
datadir = fullfile(rootdir, 'Data') ; 

scens = {'rcp45','rcp85'} ;

for iscen  = 2%:length(scens)
    scen_rcp = scens{iscen} ;
    
    if exist(fullfile(datadir,cityl,'ASHRAE_vars',scen_rcp,'statistics')) == 0 ; 
        mkdir(fullfile(datadir,cityl,'ASHRAE_vars',scen_rcp,'statistics')) ; 
    end

    for ideca = 1:length(yrstr)
        load(fullfile(datadir,cityl,'ASHRAE_vars',scen_rcp, 'ASHRAE_each_model',['ASHRAE_Data.' scen_rcp '.' yrstr{ideca} '.mat'])) ;
        nmod = length(Data) ; 
        nvar = length(VarNames) ; 
        stattype = {'Mean', 'Mode', 'Raw'} ; nstat = length(stattype) ; 

        %% Get averages
        for istat = 1:nstat ;

          Data_stat(istat).stat = stattype{istat} ;  ; % Allow for adding new statistics later
          for ivar = 1:nvar  
            eval(['tem = Data(1).' VarNames{ivar} ';']) ; 
                temout = squeeze(nan([nmod size(tem)])) ; 

                for imod = 1:nmod  
                    eval(['temout(imod,:) = Data(imod).' VarNames{ivar} '(:);']) ; 
                end

                switch(lower(Data_stat(istat).stat))

                    case 'mean'

                        eval(['Data_stat(istat).' VarNames{ivar} ' = squeeze(mean(temout)) ; ']) ; 

                    case 'mode'

                        eval(['Data_stat(istat).' VarNames{ivar} ' = squeeze(mode(temout)) ; ']) ; 

                    case 'raw'

                        eval(['Data_stat(istat).' VarNames{ivar} ' = squeeze(temout) ; ']) ; 

                end % Switch 

            end % ivar loop

        end % istat loop


        %% Write output to .csv file
        varlength = zeros(length(VarNames),1) ; 
        for i = 1:length(VarNames) ; 
            eval(['varlength(i) = length(Data_stat(1).' VarNames{i} ');']) ; %Data_stat(1)-mean value
            %eval(['varlength(i) = length(Data_stat(3).' VarNames{i} ');']) ; %Data_stat(3)-raw data
        end

        csvfile_out = ['/data/shared/Projects/Gesang/Data/',cityl,'/ASHRAE_vars/',scen_rcp,'/statistics/ASHRAE.final.',yrstr{ideca},'.csv']; 

        fid = fopen(csvfile_out, 'w') ; 
        stattype = ones(length(VarNames)) ; 
        stattype([1 10]) = 2 ; 
        for i = 1:length(VarNames) ; 
            fprintf(fid, '%s,', VarNames{i}) ; 
            eval(['tem = Data_stat(' num2str(stattype(i)) ').' VarNames{i} ';']);
            fstring2 = repmat('%f, ', [1 length(tem)]) ;
            fstring2(end+[-1 0]) = '\n' ; 
            fprintf(fid, fstring2, tem) ; 
        end
        fclose(fid) ;
    end
end


%%
CitiesLname = {'MADISON','CHICAGO','ATLANTA','BOSTON', 'DALLAS', 'HOUSTON', 'MIAMI', 'NASHVILLE',... 
    'OMAHA', 'STLOUIS', 'COLUMBUS', 'DENVER','MINNEAPOLIS',...
    'NEWYORK', 'RALEIGH', 'WASHINGTONDC',} ;
get -r RALEIGH/ASHRAE_vars/rcp85/statistics /Users/gesang/Desktop/Dropbox/ClimateBuilding/cities/RALEIGH/rcp85

