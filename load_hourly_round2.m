%% Load hourly data from .csv files
% Start by loading ALL data (including time stamp, QC, data, etc.)
% Save each year to a new .mat file


% CitiesLname = {'MADISON','CHICAGO','ATLANTA','BOSTON', 'DALLAS', 'HOUSTON', 'MIAMI', 'NASHVILLE',... 
%     'OMAHA', 'STLOUIS', 'COLUMBUS', 'DENVER', 'MINNEAPOLIS',...
%     'NEWYORK', 'RALEIGH', 'WASHINGTONDC'}
% 
% CitiesSname = {'MSN','MDW','ATL','BOS','DFW','IAH','MIA','BNA',...
%     'OMA','STL','CMH','DNE','MSP',...
%     'JFK','RDU','IAD'}

clear all
%% Section 1: Set Parameters

CitiesLname = {'DENVER'};
CitiesSname = {'DNE'};

for icity = 1:length(CitiesLname)
  
    dirin = fullfile('/data','shared','ISD',CitiesLname{icity}) ;
%     yr1 = 1980 ;
    yr1 = 1996 ; % FOR DENVER 
    yr2 = 2010 ; 
    yr_vec = yr1:yr2 ; 
    nyr = length(yr_vec) ; 

    % Initialize output
    varn = {'lat', 'lon', 'yr', 'mon', 'day', 'hr', 'min', ...
         'wnddir', 'wndmag','tmp','dew','slp'} ;
    varlim = {42:48, 52:60, 16:19, 21:22, 24:25, 27:28, 30:31, ...
          134:136, 142:145, 180:184, 190:194, 200:204} ;
   
    varflag = {'source', 'wnddirflag', 'wndobtype', 'wndmagflag', 'tmpflag', 'dewflag','slpflag'} ; 
    varflaglim = {38, 138, 140, 147, 186, 196, 206} ;

    for ivar = 1:length(varn) ; 
        eval(['Data.' varn{ivar} ' = [] ;']) ; 
    end
    for ivar = 1:length(varflag) ; 
        eval(['Data.' varflag{ivar} ' = [] ;']) ; 
    end

%% Determine how big a particular file is.
    
    for iyr = 1:nyr;
        tic
        disp(yr_vec(iyr)) ;

        % Get file name and open the file
        filin = fullfile(dirin, [CitiesSname{icity} '_' num2str(yr_vec(iyr)) '_Data.csv']) ; 
        fid = fopen(filin, 'rt') ; 
        g = textscan(fid, '%s', 'delimiter', '\n') ;
        fclose(fid) ;
        nline = length(g{1}) - 1 ; %g{1}gives all data, g{1}(1) gives the headline

        for ivar = 1:length(varn) ; 
            eval(['Data(iyr).' varn{ivar} ' = zeros(nline, 1) ;']) ; 
        end
        for ivar = 1:length(varflag) ; 
            eval(['Data(iyr).' varflag{ivar} ' = zeros(nline, 1) ;']) ; 
        end

        % Scroll through file to get data
        fid = fopen(filin, 'rt') ; 
        % Read the first line of text - this should be a header line, so we
        % will have to ignore it.
        linetext = fgetl(fid) ;
        iline = 1 ; 
        % Load next line - this is the first line of data
        linetext = fgetl(fid) ; 
        % Now, scroll through each line of the file
        while linetext > 0 ; 
            % Get data from string
            for ivar = 1:length(varn)                 
                eval(['Data(iyr).' varn{ivar} '(iline) = ' linetext(varlim{ivar}) ';']) ; 
            end
            for ivar = 1:length(varflag) ; 
                eval(['Data(iyr).' varflag{ivar} '(iline) = double(''' linetext(varflaglim{ivar}) ''');']) ;
            end
           
            % Load next line
            linetext = fgetl(fid) ;
            iline = iline + 1 ; 
        end
        toc
       
    end
    %% Save results
%     cd /data/shared/Projects/Gesang/Data/
%     save Full_Data_1980_2010.mat Data   
    dirout = fullfile('/data','shared','Projects','Gesang','Data') 
    savePath = fullfile(dirout,CitiesLname{icity}, '/Full_Data_1980_2010.mat');
    save(savePath,'Data');  
end

%% some side notes
%---------------this is for MINNEAPOLIS (data for 2010 has different varlim)
% while linetext > 0 ; 
%             % Get data from string
%             for ivar = 1:length(varn)  
%                 if iyr == nyr
%                     eval(['Data(iyr).' varn{ivar} '(iline) = ' linetext(varlim2{ivar}) ';']) ;   
%                 else
%                     eval(['Data(iyr).' varn{ivar} '(iline) = ' linetext(varlim{ivar}) ';']) ;   
%                 end
%             end
%             for ivar = 1:length(varflag) 
%                 if iyr == nyr
%                     eval(['Data(iyr).' varflag{ivar} '(iline) = double(''' linetext(varflaglim2{ivar}) ''');']) ;
%                 else
%                     eval(['Data(iyr).' varflag{ivar} '(iline) = double(''' linetext(varflaglim{ivar}) ''');']) ;
%                 end
%             end
% end
%---------------this is for WASHINGTONDC (data for 1986-1995 has different varlim)
%  while linetext > 0 ; 
%             % Get data from string
%             for ivar = 1:length(varn)  
%                 if (iyr >= 7) && (iyr < 17)
%                     eval(['Data(iyr).' varn{ivar} '(iline) = ' linetext(varlim2{ivar}) ';']) ;   
%                 
%                 else
%                     eval(['Data(iyr).' varn{ivar} '(iline) = ' linetext(varlim{ivar}) ';']) ;  
%                 end
%             end
%             for ivar = 1:length(varflag) 
%                 if (iyr >= 7) && (iyr < 17)    
%                     eval(['Data(iyr).' varflag{ivar} '(iline) = double(''' linetext(varflaglim2{ivar}) ''');']) ;
%                 else
%                     eval(['Data(iyr).' varflag{ivar} '(iline) = double(''' linetext(varflaglim{ivar}) ''');']) ;
%                 end   
%             end
% 
%             % Load next line
%             linetext = fgetl(fid) ;
%             iline = iline + 1 ; 
%         end