% Run GI motility simulation for 7 days with STC phenotype
% parameters

% MIT License
%
% Copyright (c) 2023 Raibatak Das and Applied BioMath, LLC
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

clear;

% Specify simulation inputs
parameter_file = 'parameter-table-stc.csv'; 
schedule_file = 'hapc-schedule.csv';
obs_interval = 10; % Interval between consecutive observations (min)
n_days = 7; % Number of days to run 

% Simulate time course
disp(sprintf('\nSimulating colonic motility model\n'))
disp(sprintf('\t Parameter file: %s', parameter_file))
disp(sprintf('\t Schedule file: %s', schedule_file))
disp(sprintf('\t Number of days: %d\n', n_days))

[simout, defsummary] = simulateGImotility(parameter_file, schedule_file, ...
                                          obs_interval, n_days);
% Plot simulation time course
plotGIsim(simout)
sgtitle('Simulated colonic transit for STC phenotype')

% Print defecation summary
disp(sprintf('Summary of defecation events:'))
disp('')
disp(defsummary)
