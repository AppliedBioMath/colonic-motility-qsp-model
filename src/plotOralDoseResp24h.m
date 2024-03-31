function plotOralDoseResp24h(d)
    % Plot simulated dose response to TAK-954

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
    
    figure('Position', [100 100 800 600], 'Color', 'w')
    
    C1 = [ 44, 160,  44]/255;
    C2 = [140,  86,  75]/255;
    C3 = [214,  39,  40]/255;
    C4 = [148, 103, 189]/255;
    C5 = [ 31, 119, 180]/255;
    C6 = [255, 127,  14]/255;
    
    x = d.Dose_mg;
    
    % Muscarinic receptor bound fraction
    ax1 = subplot(321);
    y = d.Mean_fbound_M;
    plot(x, y, '-', 'Color', C1, 'LineWidth', 1.5)
    hold on
    plot(x, y, 'o', 'Color', 'w', 'MarkerSize', 10, ...
         'MarkerFaceColor', C1, 'LineWidth', 2)
    ylabel('Mean f_{bound} M')
    grid on
    
    % Number of defecations
    ax2 = subplot(322);
    y = d.nDefecations;
    plot(x, y, '-', 'Color', C2, 'LineWidth', 1.5)
    hold on
    plot(x, y, 'o', 'Color', 'w', 'MarkerSize', 10, ...
         'MarkerFaceColor', C2, 'LineWidth', 2)
    ylabel('# SBM (0-24 h)')
    grid on
    
    % HAPC count
    ax3 = subplot(323);
    y = d.HAPCfrequency;
    plot(x, y, '-', 'Color', C3, 'LineWidth', 1.5)
    hold on
    plot(x, y, 'o', 'Color', 'w', 'MarkerSize', 10, ...
         'MarkerFaceColor', C3, 'LineWidth', 2)
    ylabel('HAPCs (0-24 h)')
    grid on
    
    % Average liquid fraction
    ax4 = subplot(324);
    y = d.LiquidFraction;
    plot(x, y, '-', 'Color', C4, 'LineWidth', 1.5)
    hold on
    plot(x, y, 'o', 'Color', 'w', 'MarkerSize', 10, ...
         'MarkerFaceColor', C4, 'LineWidth', 2)
    ylabel('Liquid fraction')
    grid on
    
    % Liquid influx rate
    ax5 = subplot(325);
    y = d.LiquidInfluxRate;
    plot(x, y, '-', 'Color', C5, 'LineWidth', 1.5)
    hold on
    plot(x, y, 'o', 'Color', 'w', 'MarkerSize', 10, ...
         'MarkerFaceColor', C5, 'LineWidth', 2)
    xlabel('Oral dose (mg)')
    ylabel('Liquid influx rate (g/min)')
    grid on
    
    % Consistency score
    ax6 = subplot(326);
    y = d.ConsistencyScore;
    plot(x, y, '-', 'Color', C6, 'LineWidth', 1.5)
    hold on
    plot(x, y, 'o', 'Color', 'w', 'MarkerSize', 10, ...
         'MarkerFaceColor', C6, 'LineWidth', 2)
    xlabel('Oral dose (mg)')
    ylabel('Consistency score')
    grid on
    
    % Sync x-axes
    linkaxes([ax1, ax2, ax3, ax4, ax5, ax6], 'x');
end
