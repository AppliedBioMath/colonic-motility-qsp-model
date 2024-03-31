function [simtable, defecations] = ...
            simulateGImotility(parameter_file, schedule_file, obs_interval_min, n_days)
    % Simulate GI motility model
    %
    % [simtable, defecations] = simulateGImotility(parameter_file, schedule_file, obs_interval, n_days)
    % 
    % Inputs:
    %   parameter_file: csv file with parameter table
    %   schedule_file: csv file with schedule of HAPC times and defecation
    %                  check times
    %   obs_interval_min: Time between consecutive observations (min)
    %   n_days: Number of days to run the simulation
    %
    % Outputs:
    %   simtable: Table of simulation time course
    %   defecations: Summary table of defecation events (if any)

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
    
    % Load simulation parameters
    pars = loadParameters(parameter_file);
    
    % Load HAPC schedule
    schedule = readtable(schedule_file);
    
    % Construct vector of observation times
    day2min = 24*60;
    t_obs = 0:obs_interval_min:day2min;
    
    % Simulate colon module
    [simtable, defecations] = simulateColonModule(pars, schedule, t_obs, n_days);
end
