function [out_variable, out_variable_names] = calculation4(Data_in) ;
% Temperature, Degree-Days and Degree-Hours in this code. 

%% Define output variable names
out_variable_names = {'Annual_Tavg', 'Tavg', 'Sd', 'Annual_HDD_10', 'HDD_10', ...
     'Annual_HDD_183', 'HDD_183',  'Annual_CDD_10', 'CDD_10', 'Annual_CDD_183' ...
      'CDD_183', 'Annual_CDH_233','CDH_233', 'Annual_CDH_267','CDH_267'} ; 

% You should not need to change anything below this point
% Deal out input variables
in_variable_names = {'yr', 'mo', 'hr', 'day', 'tmp', 'dew', 'wndmag', 'wnddir'} ;
[yr, mo, day, hr] = datevec(Data_in.time-.25) ; 
for ivarn = 5:length(in_variable_names) ; 
    eval([in_variable_names{ivarn} ' = Data_in.' in_variable_names{ivarn} ';']) ; 
end 
temp = tmp ; 
%% Calculations : Temperature degree-days and degree hours 
% Heating degree days 
HDD_10 = nan(1,12) ; 
HDD_183= nan(1,12) ; 
for t_base = [10 18.3]
    % the months with 31 days: 1,3,5,7,8,10,12
    Daily_DB_31 = zeros(31*25,1) ; % 31days*25years
    Monthly_Mean_31 = zeros(7,1) ;
    Monthly_SD_31 = zeros(7,1) ; 
    HDD_31 = zeros(7,1) ; 
    CDD_31 = zeros(7,1) ; 
    i = 1;
    for imo = [1 3 5 7 8 10 12]
        j = 1;
        n = 0; 
        k = 1 ;
        for iyr = [1986:2010]
            for iday = 1:31  %31days fro July 
                ind = find(mo==imo & yr==iyr & day == iday);
                t = temp(ind); %hourly data for July
                Daily_DB_31(j) = (max(t)+min(t))./2; % mean daily temp for each month of 24 years
                a = t_base - Daily_DB_31(j) ; % HDD
                
                if a>0
                    n = n+a ;
                end 
                j = j+1;
            end 
        end
        HDD_31(i) = n/25 ;
        CDD_31(i) = n/25 ;
        Monthly_Mean_31(i) = mean(Daily_DB_31) ; 
        Monthly_SD_31(i) = std(Daily_DB_31);
        i = i+1 ;
    end

    % the months with 30 days: 4 6 9 11
    Daily_DB_30 = zeros(30*25,1) ; 
    Monthly_Mean_30 = zeros(4,1) ;
    Monthly_SD_30 = zeros(4,1) ; 
    HDD_30 = zeros(4,1) ;
    CDD_30 = zeros(4,1) ;
    i = 1;
    for imo = [4 6 9 11]
        j = 1;
        n = 0 ;
        for iyr = [1986:2010]
    
            for iday = 1:30  %31days fro July 
                ind = find(mo==imo & yr==iyr & day == iday);
                t = temp(ind); %hourly data for July
                Daily_DB_30(j) = (max(t)+min(t))./2;
                a = t_base - Daily_DB_30(j) ;
                
                if a>0
                   n = n+a ;
                end 
                j = j+1;
            end 
       
        end
        HDD_30(i) = n/25 ;
        CDD_30(i) = n/25 ; 
        Monthly_Mean_30(i) = mean(Daily_DB_30) ; 
        Monthly_SD_30(i) = std(Daily_DB_30);
        i = i+1 ;
    end
    
    % Feb  ndays = 28+(mod(year, 4) == 0)
    days = zeros(24,1) ;  
    i = 1 ;
    n = 0 ; 
    for iyr = [1986:1995,1997:2010]
        days(i) = 28+(mod(iyr, 4) == 0);
        i = i+1 ; 
    end
    sum(days) ; % checking how many days in total 25 Feb

    Daily_DB = zeros(677,1) ;
    j = 1 ; 
    n = 0 ;
    for iyr = [1986:2010]
        i = 1 ; 
        for iday  = 1:days(i)
            ind = find(mo==2 & yr==iyr & day == iday);
            t = temp(ind); %hourly data for July
            Daily_DB(j) = (max(t)+min(t))./2;
            a = t_base - Daily_DB(j) ;
            
            if a>0
               n = n+a ;
            end 
           j = j+1;
        end     
    end
    CDD_2 = n/25  ;
    HDD_2 = n/25 ;
    Monthly_Mean_2 = mean(Daily_DB) ;
    Monthly_SD_2 = std(Daily_DB) ;
    
    if t_base == 10 
        HDD_10([1 3 5 7 8 10 12]) = HDD_31 ; 
        HDD_10([4 6 9 11]) = HDD_30 ; 
        HDD_10(2) = HDD_2 ; 
    else
        HDD_183([1 3 5 7 8 10 12]) = HDD_31 ; 
        HDD_183([4 6 9 11]) = HDD_30 ; 
        HDD_183(2) = HDD_2 ; 
    end
end
Annual_HDD_10 = sum(HDD_10) ;
Annual_HDD_183 = sum(HDD_183) ;

%% Cooling degree days 
CDD_10 = nan(1,12) ; 
CDD_183= nan(1,12) ; 
for t_base = [10 18.3]
    % the months with 31 days: 1,3,5,7,8,10,12
    Daily_DB_31 = zeros(31*25,1) ; % 31days*24years
    Monthly_Mean_31 = zeros(7,1) ;
    Monthly_SD_31 = zeros(7,1) ; 
    HDD_31 = zeros(7,1) ; 
    CDD_31 = zeros(7,1) ; 
    i = 1;
    for imo = [1 3 5 7 8 10 12]
        j = 1;
        n = 0; 
        k = 1 ;
        for iyr = [1986:2010]
            for iday = 1:31  %31days fro July 
                ind = find(mo==imo & yr==iyr & day == iday);
                t = temp(ind); %hourly data for July
                Daily_DB_31(j) = (max(t)+min(t))./2; % mean daily temp for each month of 24 years
                a = t_base - Daily_DB_31(j) ; % HDD
                b = Daily_DB_31(j) - t_base ;% CDD
                if b>0
                    n = n+b ;
                end 
                j = j+1;
            end 
        end
        HDD_31(i) = n/25 ;
        CDD_31(i) = n/25 ;
        Monthly_Mean_31(i) = mean(Daily_DB_31) ; 
        Monthly_SD_31(i) = std(Daily_DB_31);
        i = i+1 ;
    end
    
    % the months with 30 days: 4 6 9 11
    Daily_DB_30 = zeros(30*25,1) ; 
    Monthly_Mean_30 = zeros(4,1) ;
    Monthly_SD_30 = zeros(4,1) ; 
    HDD_30 = zeros(4,1) ;
    CDD_30 = zeros(4,1) ;
    i = 1;
    for imo = [4 6 9 11]
        j = 1;
        n = 0 ;
        for iyr = [1986:2010]
    
            for iday = 1:30  %31days fro July 
                ind = find(mo==imo & yr==iyr & day == iday);
                t = temp(ind); %hourly data for July
                Daily_DB_30(j) = (max(t)+min(t))./2;
                a = t_base - Daily_DB_30(j) ;
                b =  Daily_DB_30(j) - t_base ;
                if b>0
                   n = n+b ;
                end 
                j = j+1;
            end 
       
        end
        HDD_30(i) = n/25 ;
        CDD_30(i) = n/25 ; 
        Monthly_Mean_30(i) = mean(Daily_DB_30) ; 
        Monthly_SD_30(i) = std(Daily_DB_30);
        i = i+1 ;
    end
    
    % Feb  ndays = 28+(mod(year, 4) == 0)
    days = zeros(25,1) ;  
    i = 1 ;
    n = 0 ; 
    for iyr = [1986:2010]
        days(i) = 28+(mod(iyr, 4) == 0);
        i = i+1 ; 
    end
    sum(days) ; % checking how many days in total 24 Feb

    Daily_DB = zeros(sum(days),1) ;
    j = 1 ; 
    n = 0 ;
    for iyr = [1986:1995,1997:2010]
        i = 1 ; 
        for iday  = 1:days(i)
            ind = find(mo==2 & yr==iyr & day == iday);
            t = temp(ind); %hourly data for July
            Daily_DB(j) = (max(t)+min(t))./2;
            a = t_base - Daily_DB(j) ;
            b =  Daily_DB(j) - t_base ;
            if b>0
               n = n+b ;
            end 
           j = j+1;
        end     
    end
    CDD_2 = n/25  ;
    HDD_2 = n/25 ;
    Monthly_Mean_2 = mean(Daily_DB) ;
    Monthly_SD_2 = std(Daily_DB) ;
    
    if t_base == 10 
        CDD_10([1 3 5 7 8 10 12]) = CDD_31 ; 
        CDD_10([4 6 9 11]) = CDD_30 ; 
        CDD_10(2) = CDD_2 ; 
    else
        CDD_183([1 3 5 7 8 10 12]) = CDD_31 ; 
        CDD_183([4 6 9 11]) = CDD_30 ; 
        CDD_183(2) = CDD_2 ; 
    end
end
Annual_CDD_10 = sum(CDD_10) ;
Annual_CDD_183 = sum(CDD_183) ;

%% Tavg, STD,  Annual DB
Tavg = nan(1, 12) ; 
Sd = nan(1, 12) ;
Tavg([1 3 5 7 8 10 12]) = Monthly_Mean_31 ; 
Tavg([4 6 9 11]) = Monthly_Mean_30 ; 
Tavg(2) = Monthly_Mean_2 ;
Annual_Tavg = mean(Tavg) ;

Sd([1 3 5 7 8 10 12]) = Monthly_SD_31 ; 
Sd([4 6 9 11]) = Monthly_SD_30 ; 
Sd(2) = Monthly_SD_2 ;

%% CDH  cooling degree hours 
CDH_233 = nan(1,12) ; 
CDH_267 = nan(1,12) ; 
for t_base_CDH = [23.3 26.7] ;
    CDH = zeros(12,1) ;
    i = 1;
    for imo = 1:12
        ind = find(mo==imo);  
        t = temp(ind); %hourly data for the month
        b =  t - t_base_CDH ;
        ind_CDH = find(b > 0) ;
        x = b(ind_CDH) ;
        CDH(i) = sum(x)/24 ;
    
        i = i+1 ;      
    end  
    if t_base_CDH == 23.3
        CDH_233 = CDH ; 
    else
        CDH_267 = CDH ; 
    end   
end
Annual_CDH_233 = sum(CDH_233) ;
Annual_CDH_267 = sum(CDH_267) ;

%% Deal output variables into out_variable

out_variable = cell(length(out_variable_names), 1) ; 
for ivarn = 1:length(out_variable_names) ; 
    eval(['out_variable{ivarn} = ' out_variable_names{ivarn} ';']) ; 
end
