function doseresp24h = simulateOralDoseResp24h(pars, schedule, doses_mg)
    % Simulate dose response over 0-24h after a single oral dose of TAK-954
    %
    % doseresp24h = simulateOralDoseResp24h(parameter_file, schedule_file, doses_mg)
    % 
    % Inputs:
    %   pars: Parameter struct
    %   schedule: Table of HAPC schedule
    %   doses_mg: TAK-954 dose in mg
    %
    % Output:
    %   doseresp24h: Table of predicted response

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
    
    %% Simulate 0-24h time course for signaling module for given doses
    n_hours = 24;
    [~, d] = simulateSignalingModule(pars, n_hours, doses_mg);
    
    %% Map receptor occupancy to colon physiology module input
    hill = @ (x, xmax, ec50, n) xmax * x.^n ./ (ec50^n + x.^n); % Hill function
    
    % Extract mapping parameters
    hmax = pars.hmax;     % Maximum # of HAPCs
    wmax = pars.wmax;     % Maximum liquid influx rate
    EC50 = pars.EC50;     % EC50
    n    = pars.HillCoef; % Hill coefficient
    
    % Compute HAPC frequency and water influx rate for each dose
    f = d.Mean_fbound_M;
    d.HAPCfrequency = round(hill(f, hmax, EC50, n));
    d.LiquidInfluxRate =  hill(f, wmax, EC50, n);
    
    %% Simulate colonic motility module for 24 h with mapped parameters
    nsims = numel(doses_mg);
    pnew = pars;
    defecations = zeros(nsims, 1); % Array to store defecation frequency
    liqfrac = zeros(nsims, 1); % Array to store liquid fraction
    for jj = 1:numel(doses_mg) 
        pnew.HAPCfrequency = d(jj, :).HAPCfrequency;
        pnew.kin_w = d(jj, :).LiquidInfluxRate;
        obs_times = [0, 1440];
        [~, defs] = simulateColonModule(pnew, schedule, obs_times, 1);
        ndefs = height(defs);
        defecations(jj) = ndefs;
        if ndefs > 0
            liqfrac(jj) = mean(defs.Liquid_fraction);
        else
            liqfrac(jj) = NaN;
        end
    end
    
    % Update dose reponse table
    d.nDefecations = defecations;
    d.LiquidFraction = liqfrac;
    d.ConsistencyScore = arrayfun(@lfrac2score, liqfrac);
    doseresp24h = d;
end
