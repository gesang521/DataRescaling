% this is code is used to get latlon and elevation info for cities
clear all

%% set parameters 
CitiesLname = {'MADISON','CHICAGO','ATLANTA','BOSTON', 'DALLAS', ... 
    'HOUSTON', 'MIAMI', 'NASHVILLE','OMAHA', 'STLOUIS', 'COLUMBUS', ...
    'DENVER', 'MINNEAPOLIS','NEWYORK', 'RALEIGH', 'WASHINGTONDC'} ;
 
CitiesSname = {'MSN','MDW','ATL','BOS','DFW',...
    'IAH','MIA','BNA','OMA','STL','CMH',...
    'DNE','MSP','JFK','RDU','IAD'} ;

ncities = length(CitiesLname) ;
varn = {'lat','lon','elv'}
 
%% load lat, lon and elev

dirin = '/data/shared/Projects/Gesang/Data' ;

for icity = 1:length(CitiesLname)
    
    disp(CitiesLname{icity}) ;
    % load lat and lon from processed mat file
    filin = fullfile(dirin,CitiesLname{icity},'Full_Data_1980_2010.mat' ) ;
    load(filin) ; 
    
    [latitude] = Data.lat ;
    [longitude] = Data.lon ;
    
    eval(['Geoinfo(icity).' varn{1} ' = latitude(1) ;']) ; 
    eval(['Geoinfo(icity).' varn{2} ' = longitude(1) ;']) ; 
    
    %% load elevation from the original csv file
    dirin2 = fullfile('/data','shared','ISD',CitiesLname{icity}) ;
    yr = 2000 ;
    
    elev_varlim =  {63:67, 65:69, 63:67, 63:65, 63:67, ...
        59:62, 63:65, 65:69, 63:67,63:67,62:66, ...
        64:69, 63:67, 63:65, 63:67, 65:68 } ;
    
    filin = fullfile(dirin2, [CitiesSname{icity} '_' num2str(yr) '_Data.csv']) ; 
    fid = fopen(filin, 'rt') ;
    g = textscan(fid, '%s', 'delimiter', '/n') ;
    fclose(fid);
    
    fid = fopen(filin,'rt') ;
    linetext = fgetl(fid) ;
    iline = 1 ; 
    linetext = fgetl(fid) ; 
    
    eval(['Geoinfo(icity).' varn{3} ' = ' linetext(elev_varlim{icity}) ';']) ;
    
    
end

%% Save results
% Description = [ ...
%         'Geoinfo is a structure that contains the following: '
%         'Geoinfo.lat = latitude for 16 cities                '
%         'Geoinfo.lon = longitude for 16 cities               '                           
%         'Geoinfo.elv = elevation for 16 cities               '] ;
%     
% Geoinfo.Description = Description ;   
dirout = fullfile('/data','shared','Projects','Gesang','Data') 
savePath = fullfile(dirout, '/geoinfo.mat');
save(savePath,'Geoinfo');  

