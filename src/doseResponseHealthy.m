% Simulate dose response in healthy phenotype 0-24h after a single ascending 
% oral dose of TAK-954

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

% Specity simulation inputs
parameter_file = 'parameter-table-healthy.csv'; 
parameters = loadParameters(parameter_file);

schedule_file = 'hapc-schedule.csv';
schedule = readtable(schedule_file);

doses_mg = [0.1, 0.5, 1, 2, 5, 10, 20]; % Drug doses in mg

% Generate dose response
disp(sprintf('\nSimulating dose response for healthy phenotype\n'))
disp(sprintf('\t 0-24 h after single oral dose of TAK-954'))
disp(sprintf('\t Parameter file: %s', parameter_file))
disp(sprintf('\t Schedule file: %s\n', schedule_file))
disp(sprintf('\t Drug dose (mg): %.2f\n', doses_mg))
disp('')

dose_response24h = simulateOralDoseResp24h(parameters, schedule, doses_mg);

disp('Predicted dose response:')
disp(dose_response24h)

% Plot dose response
plotOralDoseResp24h(dose_response24h)
sgtitle('Simulated single oral dose response (0-24h) for healthy phenotype')
