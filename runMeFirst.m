% Initialize GI motility simulations. RUN THIS AT THE START OF A SESSION!

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

clear; close all; clc

% Update paths
addpath('src', '-end')
addpath('parameters', '-end')

% Print license information
type('LICENSE')

% Print short description
disp('Initializing GI motility models')
disp('')
disp('------------------------------------------------')
disp('Available demo scripts (type script name to run)')
disp('------------------------------------------------')
disp('  runHealthy          - Simulate colonic motility for healthy phenotype')
disp('  runSTC              - Simulate colonic motility for STC phenotype')
disp('  doseResponseHealthy - Simulate dose response for healthy phenotype')
disp('  doseResponseSTC     - Simulate dose response for STC phenotype')
disp('  GImodelDemo         - Run all of the above')
disp(sprintf('\n'))
