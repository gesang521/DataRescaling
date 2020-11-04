function [out_variable, out_variable_names] = calculation6(Data_in) ;
% we are calculating monthly design BD/MCWB, WB/MCDB coppesding to...
% 0.4% 2% 5% 10% respectively

%% Define output variable names
out_variable_names = {'mo_db_04', 'mo_mcwb_04', 'mo_db_2', 'mo_mcwb_2', ...
                      'mo_db_5', 'mo_mcwb_5', 'mo_db_10', 'mo_mcwb_10',...
                      'mo_wb_04', 'mo_mcdb_04', 'mo_wb_2', 'mo_mcdb_2',...
                     'mo_wb_5', 'mo_mcdb_5','mo_wb_10', 'mo_mcdb_10'} ; 

% You should not need to change anything below this point

% Deal out input variables
in_variable_names = {'yr', 'mo', 'hr', 'day', 'tmp', 'dew', 'wndmag', 'wnddir', 'tstar2'} ;
[yr, mo, day, hr] = datevec(Data_in.time-.25) ; 
for ivarn = 5:length(in_variable_names) ; 
    eval([in_variable_names{ivarn} ' = Data_in.' in_variable_names{ivarn}]) ; 
end 
temp = tmp ; 

tstar_f = tstar2 - 459.67 ; 
tstar_c = (tstar_f - 32).*(5/9) ;  % converting to celsius degree 

%% Mean daily range calculation corresponding to 5% DB
mo_db_04 = nan(1,12) ; % (treshold, month)
mo_mcwb_04 = nan(1,12) ;
mo_wb_04 = nan(1,12) ; 
mo_mcdb_04 = nan(1,12) ;

mo_db_2 = nan(1,12) ;
mo_mcwb_2 = nan(1,12) ;
mo_wb_2 = nan(1,12) ;
mo_mcdb_2 = nan(1,12) ;

mo_db_5 = nan(1,12) ;
mo_mcwb_5 = nan(1,12) ;
mo_wb_5 = nan(1,12) ;
mo_mcdb_5 = nan(1,12) ;

mo_db_10 = nan(1,12) ;
mo_mcwb_10 = nan(1,12) ;
mo_wb_10 = nan(1,12) ;
mo_mcdb_10 = nan(1,12) ;

for percent = [0.004 0.02 0.05 0.1]
    j = 1;
    for imo = 1:12
        ind = find(mo==imo);
        t = temp(ind); 
        tw = tstar_c(ind);
        
        ntim = sum(~isnan(t)) ;
        t_sort = sort(t(~isnan(t)), 'descend') ;
        threshold = floor(percent*ntim) ;
        t_threshold = t_sort(threshold) ;
        mo_db(j) = t_threshold ;%;mean(x_31) ; % average of 24 years for each month 
        %mean coincident values calculation 
        indd = find((t <= t_threshold + 5/18) & (t >= t_threshold - 5/18)) ;
        mo_mcwb(j) = nanmean(tw(indd));
        
        ntim = sum(~isnan(tw)) ;
        tw_sort = sort(tw(~isnan(tw)), 'descend') ;
        threshold = floor(percent*ntim) ;
        tw_threshold = tw_sort(threshold) ;
        mo_wb(j) = tw_threshold ;%;mean(x_31) ; % average of 24 years for each month 
        %mean coincident values calculation 
        indd = find((tw <= tw_threshold + 5/18) & (tw >= tw_threshold - 5/18)) ;
        mo_mcdb(j) = nanmean(t(indd));
        
        j = j+1 ; 
    end
    
    if percent == 0.004
        mo_db_04 = mo_db  ; 
        mo_mcwb_04 = mo_mcwb  ;
        mo_wb_04 = mo_wb  ; 
        mo_mcdb_04 = mo_mcdb ;
    elseif percent == 0.02
        mo_db_2 = mo_db  ; 
        mo_mcwb_2 = mo_mcwb  ;
        mo_wb_2 = mo_wb  ; 
        mo_mcdb_2 = mo_mcdb ;
    elseif percent == 0.05
        mo_db_5 = mo_db  ; 
        mo_mcwb_5 = mo_mcwb  ;
        mo_wb_5 = mo_wb  ; 
        mo_mcdb_5 = mo_mcdb ;
    else
        mo_db_10 = mo_db  ; 
        mo_mcwb_10 = mo_mcwb ;
        mo_wb_10 = mo_wb  ; 
        mo_mcdb_10 = mo_mcdb ;
    end
        
end

%% Deal output variables into out_variable

out_variable = cell(length(out_variable_names), 1) ; 
for ivarn = 1:length(out_variable_names) ; 
    eval(['out_variable{ivarn} = ' out_variable_names{ivarn} ';']) ; 
end 
        
            
            
            
            
           

