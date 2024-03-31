function score = lfrac2score(l)
    % Convert liquid fraction to 7-point stool consistency score
    % (approximating the Bristol stool score)
    
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
    
    if (l > 0) & (l <= 0.5)
        score = 1;
    elseif (l > 0.5) & (l <= 0.6)
        score = 2;
    elseif (l > 0.6) & (l <= 0.7)
        score = 3;
    elseif (l > 0.7) & (l <= 0.8)
        score = 4;
    elseif (l > 0.8) & (l <= 0.85)
        score = 5;
    elseif (l > 0.85) & (l <= 0.9)
        score = 6;
    elseif (l > 0.9) & (l < 1)
        score = 7;
    else
        score = 0;
    end
end
