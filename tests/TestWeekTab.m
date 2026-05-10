classdef TestWeekTab < matlab.unittest.TestCase
    methods (Test)
        function aggregatesWeeklyTotals(testCase)
            addpath(genpath('src'));
            weekly = week_tab('sample-data/sir_covid/al_covid.csv');
            testCase.verifyEqual(weekly.week(1), 1);
            testCase.verifyGreaterThan(height(weekly), 0);
            testCase.verifyEqual(width(weekly), 3);
        end
    end
end
