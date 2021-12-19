%% weekly table function
%
% extract a weekly table from data
% 
% 
function table1 = week_tab(data)
    
    % read table
    dtable = readtable(data);
        
    % weekly cases
    dtable.day_cases = dtable.cases;
    dtable.day_cases(2:end) = dtable.cases(2:end)-dtable.cases(1:end-1);
    cases = dtable.day_cases;
    cases(end+1:7*ceil(length(cases)/7)) = 0;
    cases = reshape(cases,7,length(cases)/7);
    cases = sum(cases,1)';
    
    % weekly deaths
    dtable.day_deaths = dtable.deaths;
    dtable.day_deaths(2:end) = dtable.deaths(2:end)-dtable.deaths(1:end-1);
    deaths = dtable.day_deaths;
    deaths(end+1:7*ceil(length(deaths)/7)) = 0;
    deaths = reshape(deaths,7,length(deaths)/7);
    deaths = sum(deaths,1)';
    
    
    week = (1:length(cases))';
    table1 = table(week,cases,deaths);
end