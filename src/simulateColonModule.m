function [simtable, defecations] = ...
            simulateColonModule(parameters, schedule, obs_times_min, n_days)
    % Simulate time course for colon physiology module
    %
    % [simtable, defecations] = simulateColonModule(parameters, schedule, obs_times, n_days)
    % Inputs:
    %   parameters: struct with model parameters
    %   schedule: table of HAPC times and defecation check times
    %   obs_times_min: Observation times (min)
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
    
    % Create shortcuts
    p = parameters;
    kin_m         = p.kin_m; % Solid mass influx rate (g/min)
    kin_w         = p.kin_w; % Liquid mass influx rate (g/min)
    ktrans        = p.ktrans; % Colonic transport rate constant (min^-1)
    kloss         = p.kloss_m; % Solid mass loss rate constant (min^-1)
    kabs          = p.kabs_w; % Liquid mass absorption rate constant (min^-1)
    mthresh       = p.mthresh; % Defecation threshold mass (g)
    hapc_per_day  = p.HAPCfrequency; % HAPCs frequency (per day)
    hapc_duration = p.HAPCduration; % Duration of an HAPC per compartment (min)
    % Initial masses 
    y0 = [p.m1_0; p.w1_0; p.m2_0; p.w2_0; p.m3_0; p.w3_0; p.m4_0; p.w4_0];
    
    %% Construct HAPCs
    % Extract HAPC and check times
    [hapc_times, check_times] = getEventTimes(hapc_per_day, schedule);
    
    % Construct HAPC function 
    f = @(t) hapc(t, hapc_times, hapc_duration);
    
    %% Split day into chunks of times
    % Create vector of observation times
    day2min = 1440;
    tobs = unique([0, obs_times_min(:)', day2min, check_times, hapc_times, ...
                   hapc_times + hapc_duration, hapc_times + 2*hapc_duration, ...
                   hapc_times + 3*hapc_duration]);

    % Split into chunks using check times. Defecation check happens at the 
    % end of each chunk
    nchunks = numel(check_times);
    tchunks = cell(1, nchunks);
    breaks = cell(1, nchunks);
    obs_per_chunk = zeros(1, nchunks);
    nbreaks = zeros(1, nchunks);
    tsplit = unique([0, check_times]);
    % Compute observation times and break points within each chunk
    for jj = 1:nchunks
        idx = (tobs >= tsplit(jj)) & (tobs < tsplit(jj+1));
        tchunks{jj} = unique([tobs(idx), tsplit(jj+1)]);
        obs_per_chunk(jj) = numel(tchunks{jj});
        th = hapc_times((hapc_times >= tsplit(jj)) & ...
                        (hapc_times < tsplit(jj + 1)));
        breaks{jj} = unique([th, th + hapc_duration, th + 2*hapc_duration, ...
                             th + 3*hapc_duration, tsplit(jj + 1)]); 
        nbreaks(jj) = numel(breaks{jj});        
    end
    obs_per_day = sum(obs_per_chunk);
    indexer = cumsum([0, obs_per_chunk]);
    
    %% Simulate system
    simout = zeros(12, n_days * obs_per_day); % Preallocate output array
    for ii = 1:n_days
        t_offset = (ii-1) * day2min;
        idx_offset = (ii-1) * obs_per_day;
        % Initialize ODE solver
        sol = ode15s(@colonModuleODEs, [0, eps], y0);
        % Integrate over each chunk of time
        for jj = 1:nchunks
            tbreaks = breaks{jj};
            for kk = 1:nbreaks(jj)
                sol = odextend(sol, [], tbreaks(kk));
            end
            y = deval(sol, tchunks{jj});
            % Update simulation table
            % Specify array indices
            idx_start = idx_offset + (indexer(jj) + 1);
            idx_end = idx_offset + indexer(jj+1);
            idx = idx_start:idx_end;
            % Insert time points
            %t_rel = tchunks{jj};
            simout(1, idx) = t_offset + tchunks{jj};
            % Insert simulation output
            simout(2:9, idx) = y;
            % Perform defecation check and update initial value for next chunk
            y0 = sol.y(:, end);
            mass_in_comp4 = sum(y0(7:8));
            if mass_in_comp4 > mthresh
                % Insert defecation event
                simout(10:11, idx_end) = y0(7:8); % Defecation mass (solid + liquid)
                simout(12, idx_end) = y0(8)/mass_in_comp4; % liquid fraction
                y0(7:8) = [eps; 0]; % Update initial value for next chunk
                simout(8:9, idx_end) = [eps; 0]; % Empty last compartment
            end
            % Update ODE solver for next chunk
            t0 = sol.x(end);
            if jj < nchunks
                sol = ode15s(@colonModuleODEs, [t0, t0 + eps(t0)], y0);
            end
        end
    end
    
    %% Construct output table
    s = array2table(simout', 'VariableNames', ...
                              {'Time_min', 'm1_g', 'w1_g', 'm2_g', 'w2_g', ...
                               'm3_g', 'w3_g', 'm4_g', 'w4_g', ...
                               'mdef_g', 'wdef_g', 'ldef'});
    s = unique(s, 'rows', 'stable');
    s.Time_day = s.Time_min/day2min;
    
    % Compute liquid fraction in each compartment
    s.l1 = s.w1_g ./ (s.m1_g + s.w1_g);
    s.l2 = s.w2_g ./ (s.m2_g + s.w2_g);
    s.l3 = s.w3_g ./ (s.m3_g + s.w3_g);
    s.l4 = s.w4_g ./ (s.m4_g + s.w4_g);
    
    % Add columns for HAPCs
    f = @(t) hapc(mod(t, day2min), hapc_times, hapc_duration);
    s.f1 = arrayfun(f, s.Time_min);
    s.f2 = arrayfun(f, s.Time_min - hapc_duration);
    s.f3 = arrayfun(f, s.Time_min - 2*hapc_duration);
    
    simtable = s(:, {'Time_min', 'Time_day', 'm1_g', 'w1_g', 'l1', ...
                     'm2_g', 'w2_g', 'l2', 'm3_g', 'w3_g', 'l3', ...
                     'm4_g', 'w4_g', 'l4', 'mdef_g', 'wdef_g', 'ldef', ...
                     'f1', 'f2', 'f3'});
    
    % Summarize defecation events
    defs = s(s.mdef_g > 0, :);
    defecations = table();
    if isempty(defs)
        warning("No defecation observed");
    else
        defecations.Time_day = defs.Time_day;
        defecations.Defecation_mass_g = defs.mdef_g + defs.wdef_g;
        defecations.Liquid_fraction = defs.ldef;
    end
    
    %% Define ODEs
    function dydt = colonModuleODEs(t, y)
        % ODEs for colonic motility module. See manuscript for details.
        
        % State variables
        m1 = y(1); % Solid mass in compartment 1
        w1 = y(2); % Liquid mass in compartment 1
        m2 = y(3); % Solid mass in compartment 2
        w2 = y(4); % Liquid mass in compartment 2
        m3 = y(5); % Solid mass in compartment 3
        w3 = y(6); % Liquid mass in compartment 3
        m4 = y(7); % Solid mass in compartment 4
        w4 = y(8); % Liquid mass in compartment 4
        
        % Liquid fractions
        l1 = w1/(m1 + w1); % Liquid fraction in compartment 1
        l2 = w2/(m2 + w2); % Liquid fraction in compartment 2
        l3 = w3/(m3 + w3); % Liquid fraction in compartment 3
        
        % HAPCs
        f1 = f(t); % HAPC in compartment 1
        f2 = f(t - hapc_duration); % HAPC in compartment 2
        f3 = f(t - 2*hapc_duration); % HAPC in compartment 3
        
        % Rate laws
        % Compartment 1
        dm1dt = kin_m - kloss*m1 - f1*ktrans*l1*m1;
        dw1dt = kin_w - kabs*w1  - f1*ktrans*l1*w1;
        % Compartment 2
        dm2dt = f1*ktrans*l1*m1 - kloss*m2 - f2*ktrans*l2*m2;
        dw2dt = f1*ktrans*l1*w1 - kabs*w2  - f2*ktrans*l2*w2;
        % Compartment 3
        dm3dt = f2*ktrans*l2*m2 - kloss*m3 - f3*ktrans*l3*m3;
        dw3dt = f2*ktrans*l2*w2 - kabs*w3  - f3*ktrans*l3*w3;
        % Compartment 4
        dm4dt = f3*ktrans*l3*m3 - kloss*m4;
        dw4dt = f3*ktrans*l3*w3 - kabs*w4;
        
        dydt = [dm1dt; dw1dt; dm2dt; dw2dt; dm3dt; dw3dt; dm4dt; dw4dt];
    end
    
    %% Define HAPCs
    function h = hapc(t, tstart, delta_t)
        % Rectangular pulses at times tstart and of duration delta_t to 
        % represent HAPCs in colonic motility module
        n = numel(tstart);
        tvec = t * ones(1, n);
        % Check if t is within + delta_t of any start time 
        check = (tvec >= tstart) & (t < tstart + delta_t);
        if any(check)
            h = 1;
        else
            h = 0;
        end
    end 
end
