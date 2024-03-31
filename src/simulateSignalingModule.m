function [simtable, summary] = simulateSignalingModule(pars, n_hours, doses_mg)
    % Simulate 5-HT4 signaling module for supplied doses of TAK-954
    %
    % [simtable, summary] = simulateSignalingModule(parameter_file, n_hours, doses_mg)
    % 
    % Inputs:
    %   pars: Parameter struct
    %   n_hours: Number of hours to simulate
    %   doses_mg: TAK-954 dose in mg
    %
    % Outputs:
    %   simtable: Table of simulation time course
    %   summary: Summary table of average receptor occupancy of 5-HT4 and
    %            muscarinic receptors over simulation period
    
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
    % PK parameters
    drugMW = pars.drugMW; 
    Vc     = pars.Vcentral;      % Volume of central compartment (plasma)
    Vp     = pars.Vperipheral;   % Volume of peripheral compartment
    fb     = pars.fbound_pp;     % Plasma protein-bound fraction
    kabs   = pars.kabs_drug;          % Gut -> Central absorption rate constant
    kc2p   = pars.kc2p;          % Central -> Peripheral transport rate constant
    kp2c   = pars.kp2c;          % Peripheral -> Cental transport rate constant
    kel    = pars.kel;           % Elimination rate constant
    % PD parameters
    kbind_drug       = pars.kbind_drug; % Drug + 5-HT4 receptor binding rate constant
    kunbind_drug_ht4 = pars.kunbind_drug; % Drug:5-HT4 dissociation rate constant
    ksyn_5ht         = pars.ksyn_5ht; % 5-HT synthesis rate (represents constant influx)
    kbind_5ht        = pars.kbind_5ht; % 5-HT + 5-HT4 binding rate constant
    kunbind_5ht      = pars.kunbind_5ht; % 5-HT:5-HT4 dissociation rate constant
    kdeg_5ht         = pars.kdeg_5ht; % 5-HT degradation rate constant
    ksyn_5ht4r       = pars.ksyn_5ht4r; % 5-HT4 receptor synthesis rate
    kdeg_5ht4r       = pars.kdeg_5ht4r; % 5-HT4 receptor basal degradation rate constant
    kint_5ht4r       = pars.kint_5ht4r; % Bound 5-HT4 receptor internalization rate constant
    ksyn_achi        = pars.ksyn_achi; % Intracellular acetylcholine synthesis rate constant
    kdeg_achi        = pars.kdeg_achi; % Intracellular acetylcholine degradation rate constant
    krelease         = pars.krelease; % Acetylcholine release rate
    kbind_ach        = pars.kbind_ach; % ACh + M binding rate constant
    kunbind_ach      = pars.kunbind_ach; % ACh:M dissociation rate constant
    kclear_ach       = pars.kclear_ach; % Extracellular acetylcholine clearance rate constant
    ksyn_m           = pars.ksyn_m; % Muscarinic receptor synthesis rate
    kdeg_m           = pars.kdeg_m; % Muscarinic receptor degradation rate constant
    % Conversion factors
    day2sec = 24*60*60;
    hour2sec = 60*60;
    min2sec = 60;
    milli2nano = 1e6;
    
    % Compute protein binding/unbinding rates
    kbind_pp = pars.scaling_factor*fb;
    kunbind_pp = pars.scaling_factor*(1 - fb);
    
    % Convert doses from mg to nmol
    doses_nmol = doses_mg * milli2nano / drugMW;
    
    %% Simulate to drug-free steady state
    ode_opts = odeset('Events', @detectSteadyState);
    tmax = 7*day2sec;
    [~, ~,tss, yss, ie] = ode15s(@signalingModuleODEs, [0, tmax], ...
                                 zeros(12, 1), ode_opts);
    % disp(ie); % Termination flag
    % disp(tss/3600); % Hours to reach drug-free steady state
    % disp(yss'); % Drug-free steady state
    % Warn if model fails to reach drug-free steady state by tmax
    assert(~isempty(ie), ...
           "Failed to reach drug free steady state after 7 days! Please check model parameters.");
    
    % Use drug-free steady state as initial condition
    y0 = [0, 0, 0, yss(5:7), 0, yss(9:end)]';
   
    %% Simulate time course for each dose
    tf = n_hours*hour2sec;
    time_points = 0:min2sec:tf;
    nobs = numel(time_points);
    nsims = numel(doses_mg);
    sim = table();
    for jj = 1:nsims
        d0mg = doses_mg(jj);
        d0nmol = doses_nmol(jj);
        
        % Integrate using initial values from drug-free steady state
        sol = ode15s(@signalingModuleODEs, [0, tf], [d0nmol; y0]); 
        ypred = deval(sol, time_points);

        % Conver output to table
        s = array2table([d0mg*ones(nobs, 1), time_points(:), ypred'], ...
                        'VariableNames', {' Dose_mg', 'Time_sec', ...
                                          'dg', 'df', 'db', 'dp', ...
                                          'Serotonin_nM', 'Free_5HT4_nM', ...
                                          'Serotonin_5HT4_nM', 'Drug_5HT4_nM', ...
                                          'ACh_i_nM', 'ACh_nM', 'Free_Mrec_nM', ...
                                          'ACh_Mrec_nM'});
        sim = [sim; s];
    end
    
    %% Generate output table
    sim.Time_min = sim.Time_sec/min2sec;
    sim.fbound_5HT4 = (sim.Serotonin_5HT4_nM + sim.Drug_5HT4_nM) ./ ...
        (sim.Free_5HT4_nM + sim.Serotonin_5HT4_nM + sim.Drug_5HT4_nM);
    sim.fbound_Mrec = sim.ACh_Mrec_nM ./ (sim.Free_Mrec_nM + sim.ACh_Mrec_nM);
    sim.Total_drug_plasma_ngperL = (sim.df + sim.db) * drugMW / Vc ; % [Drug]_total in blood (ng/L)
    sim.Free_drug_plasma_ngperL = sim.df * drugMW / Vc; % [Drug]_free in blood (ng/L)
    sim.Peripheral_drug_ngperL = sim.dp * drugMW / Vp; % [Drug]_peripheral (ng/L)
    
    simtable = sim(:, {'Dose_mg', 'Time_sec', 'Time_min', ...
                       'Total_drug_plasma_ngperL', 'Free_drug_plasma_ngperL', ...
                       'Peripheral_drug_ngperL', 'Serotonin_nM', ...
                       'Free_5HT4_nM', 'Serotonin_5HT4_nM', 'Drug_5HT4_nM', 'fbound_5HT4', ...
                       'ACh_i_nM', 'ACh_nM', 'Free_Mrec_nM', 'ACh_Mrec_nM', ...
                       'fbound_Mrec'});
    
    %% Summarize mean receptor occupancy -vs- Dose (Integrated AUC / Interval length)
    s = sim(:, {'Dose_mg', 'fbound_5HT4', 'fbound_Mrec'});
    d = grpstats(s, 'Dose_mg', @ (x) trapz(x)/numel(x), ...
                            'VarNames', {'Dose_mg', 'Obs_count', 'Mean_fbound_5HT4', 'Mean_fbound_M'});
    
    summary = d(:, {'Dose_mg', 'Mean_fbound_5HT4', 'Mean_fbound_M'});

    %% ODEs for signaling module
    function dydt = signalingModuleODEs(t, y)
        % Species
        dg        = y(1);  % Drug dose (ng)
        df        = y(2);  % Free drug in plasma i.e central compartment (ng)
        db        = y(3);  % Plasma protein-bound drug (ng)
        dp        = y(4);  % Drug in peripheral compartment (ng)
        ht        = y(5);  % Serotonin (5-HT) (nM)
        ht4r      = y(6);  % 5-HT4 receptor (nM)
        ht_ht4r   = y(7);  % 5-HT bound 5-HT4 receptor (nM)
        drug_ht4r = y(8);  % Drug-bound 5-HT4 receptor (nM)
        ach_i     = y(9);  % Intracellular acetylcholine (nM)
        ach       = y(10); % Extracellular acetylcoline (nM)
        mr        = y(11); % Free muscarinic receptor (nM)
        ach_mr    = y(12); % Bound muscarinic receptor (nM)
        
        % Rate laws
        ddgdt        = -kabs*dg;
        ddfdt        =  kabs*dg - kel*df ...
                            - kbind_pp*df + kunbind_pp*db ...
                            - kc2p*df + kp2c*dp  ...
                            - kbind_drug*df*ht4r + kunbind_drug_ht4*drug_ht4r;
        ddbdt        = kbind_pp*df - kunbind_pp*db;
        ddpdt        = kc2p*df - kp2c*dp;
        dhtdt        = ksyn_5ht - kdeg_5ht*ht ...
                            - kbind_5ht*ht*ht4r + kunbind_5ht*ht_ht4r;
        dht4rdt      = ksyn_5ht4r - kdeg_5ht4r*ht4r ...
                            - kbind_5ht*ht*ht4r + kunbind_5ht*ht_ht4r ...
                            - kbind_drug*df*ht4r + kunbind_drug_ht4*drug_ht4r;
        dht_ht4rdt   = kbind_5ht*ht*ht4r - kunbind_5ht*ht_ht4r - kint_5ht4r*ht_ht4r;
        ddrug_ht4rdt = kbind_drug*df*ht4r - kunbind_drug_ht4*drug_ht4r - kint_5ht4r*drug_ht4r;
        dach_idt     = ksyn_achi - kdeg_achi*ach_i ...
                            - krelease*ach_i*(ht_ht4r + drug_ht4r);
        dachdt       = krelease*ach_i*(ht_ht4r + drug_ht4r) - kclear_ach*ach ...
                            - kbind_ach*ach*mr + kunbind_ach*ach_mr;
        dmrdt        = ksyn_m - kdeg_m*mr - kbind_ach*ach*mr + kunbind_ach*ach_mr;
        dach_mrdt    = kbind_ach*ach*mr - kunbind_ach*ach_mr - kdeg_m*ach_mr;
        
        % Derivative
        dydt = [ddgdt; ddfdt; ddbdt; ddpdt; dhtdt; dht4rdt; dht_ht4rdt; ...
                ddrug_ht4rdt; dach_idt; dachdt; dmrdt; dach_mrdt];
    end
    
    %% Event detection function to check for steady state
    function [value, isterminal, direction] = detectSteadyState(t, y)
        time_scale = 10;
        AbsTol = 1e-6; % |dydt * time_scale| <= AbsTol
        RelTol = 1e-9; % |dydt * time_scale / y| <= RelTol
        dydt = signalingModuleODEs(t, y);
        idx = [5, 6, 7, 9, 10, 11, 12]; % Variables to check for steady state
        deltas = abs(time_scale * dydt(idx));
        ys = abs(y(idx));
        ssAbs = (deltas <= AbsTol);
        ssRel = (deltas ./ ys <= RelTol);
        value = double(all(ssAbs | ssRel));
        isterminal = 1; % Stop integrating at steady state
        direction = [];
    end
    
end
