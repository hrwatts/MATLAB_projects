function results = run_smoke_checks()
%RUN_SMOKE_CHECKS Run a small automated smoke-test subset for the repo.

addpath(genpath('src'));
results = runtests({'tests/TestEulerMethod.m', ...
                    'tests/TestNewtonMethod.m', ...
                    'tests/TestWeekTab.m', ...
                    'tests/TestSirModel.m'});
disp(table(results))
end
