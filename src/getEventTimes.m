function [hapc_times, check_times] = getEventTimes(hapc_freq, schedule)
    % Extract HAPC times and defecation check times for given HAPC frequency
    %
    % [hapc_times, check_times] = getEventTimes(hapc_freq, schedule)
    %
    % Inputs
    %   hapc_freq: # of HAPCs per day
    %   schedule: Schedule table (load from csv schedule file)
    %
    % Outputs
    %   hapc_times: Vector of times for HAPC events
    %   check_times: Vector of times for defecation checks
    
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

    % Convert events to categorical variable
    schedule.Event = categorical(schedule.Event);

    % Extract correct event times from schedule
    t = schedule(schedule.Frequency == hapc_freq, :);
    t_hapc = t(t.Event == 'HAPC', :).Time';
    t_check = t(t.Event == 'Check', :).Time';

    % Return sorted lists after removing any duplicates
    hapc_times = unique(t_hapc);
    check_times = unique(t_check);
end
