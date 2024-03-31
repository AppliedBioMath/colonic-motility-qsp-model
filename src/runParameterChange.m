% Demo script to show examples of simulating GI motility model with
% parameter updates

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

%% Load baseline parameters and HAPC schedules
parameter_file = 'parameter-table-healthy.csv'; 
basepars = loadParameters(parameter_file); % Creates parameter struct
disp(sprintf('Baseline parameters: \n'))
disp(basepars)

schedule_file = 'hapc-schedule.csv';
schedule = readtable(schedule_file);

%% Example 1: Change transport rate and simulate colon physiology module for 
%  healthy phenotype for 1 week 

% Update parameters - decrease transport rate by 50%
newpars = basepars;
newpars.ktrans = 0.5 * basepars.ktrans;
disp(sprintf('Updated parameters: \n'))
disp(newpars)

% Specify simulation parameters
obs_times = 0:10:1440; % minutes
n_days = 7;

% Simulate colon physiology module
[simout, defsummary] = simulateColonModule(newpars, schedule, obs_times, n_days);
                                      
% Plot simulation time course
plotGIsim(simout)
sgtitle('Simulated colonic transit for healthy phenotype with updated parameters')

% Print defecation summary
disp(sprintf('Summary of defecation events:'))
disp('')
disp(defsummary)

%% Example 2: Change binding rate of acetylcholine to muscarinic receptor and 
%  simulate dose response for healthy phenotype

% Update parameters - increase acetylcholine binding rate constant by 20%
newpars = basepars;
newpars.kbind_ach = 1.2 * basepars.kbind_ach;
disp(sprintf('Updated parameters: \n'))
disp(newpars)

% Simulate dose response with updated parameters
doses_mg = [0.1, 0.5, 1, 2, 5, 10, 20]; % Drug doses in mg

dose_response24h = simulateOralDoseResp24h(newpars, schedule, doses_mg);

disp('Predicted dose response:')
disp(dose_response24h)

% Plot dose response
plotOralDoseResp24h(dose_response24h)
sgtitle('Predicted single oral dose response (0-24h) for healthy phenotype with updated parameters')

