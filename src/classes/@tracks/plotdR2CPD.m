function varargout = plotdR2CPD(obj, ha, indices, delayIdx)
% plotdR2 Plots a 1-cummulative probability function of the squared 
% displacement per population
% 
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% plotSqDisp(obj, ha, indices, delayIdx, logbinFlag, normalization)
%   obj: allows >1 objs to be parsed
%   ha: handle to axes. If empty the current axis is taken
%   indices: indeces of the tracks to be displayed
%   delayIdx: the index/indices of the delay time(s) to be plotted
%
% h = plotdR2CPD(...) returns the handle to the line plotted.
%
% [h, ha] = plotCPD(...) also returns the handle of the axes in which 
% the curve was plotted.

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% 
% 
% -------------------------------------------------------------
% Copyright (C) 2021 J.J. Erik Maris
% Inorganic Chemistry & Catalysis, Utrecht University
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

if nargin < 2 || isempty(ha)
    ha = gca;
end

if nargin < 4 || isempty(delayIdx)
    delayIdx = 1;
end

nObj = numel(obj);

if nObj > 1 && numel(delayIdx) > 1
    error('tracks:plotSqDisp:multipleObjDelayIdx',...
        'Cannot plot multiple delay indices for multiple objects.')
end

all_indices = false;
if nargin < 3 || isempty(indices)
    all_indices = true;
end

if nObj > 1 % plot multiple populations; get cell with this
    dR2 = cell(nObj,1);
    for ii = 1:nObj
        if all_indices
            indices = 1:obj(ii).nTracks;
        end
        dR2{ii} = obj(ii).getdR2(delayIdx,indices);
        dR2{ii} = dR2{ii}{1}; % only one delayIdx
        [uF,uL] = obj(ii).getUnitFactor('pixelsize.^2');
        dR2{ii} = dR2{ii}.*uF;
    end
    nPlot = numel(obj);
else % plot multiple delay times; get cell with this
    if all_indices
        indices = 1:obj.nTracks;
    end
    [dR2, dt] = obj.getdR2(delayIdx,indices);
    [uF,uL] = obj.getUnitFactor('pixelsize.^2');
    dR2 = cellfun(@(x) x.*uF, dR2,'UniformOutput',false);
    nPlot = numel(delayIdx);
end

% data has already been sorted in getdR2 from low to high dR2
% compute the cumsum...
iCPD = cellfun(@(x) linspace(1,0,numel(x)), dR2,'UniformOutput',false);
% ...and normalize to get the cpd. Compute the inverse cpd for plotting
% iCPD = cellfun(@(x) 1-x./x(end), iCPD,'UniformOutput',false);


hps = gobjects(nPlot,1);
for ii = 1:nPlot
    
    if isempty(dR2{ii})
        warning('tracks:plotdRCPD:EmptySquaredDisplacement',...
        'Squared displacement of object %i has not been calculated. Skipping object.',ii)
        continue
    end

    if nObj > 1
        jj = ii;
        plotName = sprintf('Population %d', ii );
    else
        jj = 1;
        [uFtime, uLtime] = obj(jj).getUnitFactor('dt');
        plotName = sprintf('Delay time %f%s', dt(ii).*uFtime, uLtime);
    end
    
    colors = obj(jj).getColormap(nPlot);
    LW = obj(jj).lineWidth;
    
    hps(ii) = plot(ha,...
        dR2{ii},iCPD{ii}, ...
        '.', ...
        'Color', colors(ii,:), ...
        'LineWidth', LW, ...
        'DisplayName', plotName);

    hold(ha, 'on');
end

hold(ha, 'off');

set(ha,'YScale','log')

if nPlot > 1
    legend
end

xlabel(['Squared displacement (' uL ')'])
ylabel('CPD')


if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end
end
