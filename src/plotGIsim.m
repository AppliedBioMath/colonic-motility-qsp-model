function plotGIsim(simout)
    % Plot colonic motility simulation time course

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
    
    % Define colors
    C1 = [ 31, 119, 180]/255;
    C4 = [214,  39,  40]/255;
    C5 = [140,  86,  75]/255;
    
    % Extract variables to plot
    s = simout;
    t = s.Time_day;
    t0 = min(t);
    tf = max(t);
    
    figure('Position', [100 100 800 600], 'Color', 'w')
    
    % HAPC episodes
    ax1 = subplot(6, 1, 1);
    x = s(s.f1 ~= 0, :).Time_day; % HAPC start times
    y = ones(size(x));
    scatter(x, y, 72, 'kd');
    ylabel('HAPC')
    yticks([])
    
    % Mass in compartment 1
    ax2 = subplot(6, 1, [2, 3]);
    m = s.m1_g;
    mt = s.m1_g + s.w1_g;
    x = [t; flip(t)];
    y = [mt; flip(m)];
    p1 = patch(x, y, C1);
    p1.EdgeColor = C1;
    p1.FaceAlpha = 0.6;
    hold on
    x = [t; tf; t0];
    y = [m; 0; 0];
    p2 = patch(x, y, C1);
    p2.FaceAlpha = 0.9;
    p2.EdgeColor = C1;
    ylabel('m_1, w_1 (g)')
    
    % Mass in compartment 4
    ax3 = subplot(6, 1, [4, 5]);
    m = s.m4_g;
    mt = s.m4_g + s.w4_g;
    x = [t; flip(t)];
    y = [mt; flip(m)];
    p1 = patch(x, y, C4);
    p1.EdgeColor = C4;
    p1.FaceAlpha = 0.6;
    hold on
    x = [t; tf; t0];
    y = [m; 0; 0];
    p2 = patch(x, y, C4);
    p2.FaceAlpha = 0.9;
    p2.EdgeColor = C4;
    ylabel('m_4, w_4 (g)')
    
    % Defecation events
    ax4 = subplot(6, 1, 6);
    d = s(s.mdef_g > 0, :);
    x = d.Time_day;
    y = ones(height(d), 1);
    m = d.mdef_g;
    mt = d.mdef_g + d.wdef_g;
    mul = 4;
    p1 = scatter(x, y, mt*mul, C5, 'MarkerFaceColor', C5);
    p1.MarkerEdgeColor = C5;
    p1.MarkerFaceAlpha = 0.6;
    hold on
    p2 = scatter(x, y, m*mul, C5, 'MarkerFaceColor', C5);
    p2.MarkerEdgeColor = 'w';
    p2.LineWidth = 1.5;
    xlabel('Time (day)');
    yticks([])
    ylabel('Defecation')    
    
    % Sync x-axes
    linkaxes([ax1, ax2, ax3, ax4], 'x');
    xlim([floor(t0), ceil(tf)]);
    
end
