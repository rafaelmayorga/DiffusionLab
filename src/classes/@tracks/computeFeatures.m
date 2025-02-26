function obj = computeFeatures(obj,preserveFlag,featureNames)
%COMPUTEFEATURES Computes the track features and stores them in a table
%       'firstFrame': first frame in which track appears
%       'firstFrameCoords': coordinates of the first frame of a track
%       'offFrames': frame numbers when no particle is detected during the
%           track. Is a result of a blinking > 0
%       'nTrackPointsOff': number of tracks points at which no particle is
%           detected
%       'nTrackPointsAll': number of tracks points between first and last
%           detection of a track.
%       'MinBoundCircleRadius'/vector: radius of minimum enclosing circle
%           of the track (XY)
%       'MinBoundCircleCenter'/matrix: center coordinates of minimum 
%           enclosing circle of the track (XY)
%       'CenterOfMass'/matrix: center of mass of the track
%       'MBCCminusCoM'/matrix: % get difference between CoM and center of 
%           minimum bounding circle in percent. Idea: if this difference is
%           small the track is a 'ball of wool', if the difference is large
%           it is an 'uncoiled ball of wool': in percent of 
%           'MinBoundCircleRadius' (i.e. between 0 to 1).
%
% computeEntropy Computes the Shannon entropy of a track
% 
% Track coordinates are set on a binary grid with the size of the minimum
% bounding cirle's diameter and divided in resolution^2 equally sized
% blocks. The entropy is Shannon's entropy of a binary image.

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
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

nObj = numel(obj);

if nargin < 2
    preserveFlag = false;
end

if nargin < 3
    featureNamesFlag = false;
else
    featureNamesFlag = true;
end

wb = waitbar(0,'Computing Features','Name','Computing all track features');

for ii = 1:nObj
    waitstr = ['Computing Features (' num2str(ii) '/' num2str(nObj) ')'];
    waitbar(0.2/nObj+(ii-1)/nObj,wb,waitstr)
    
    %% Dimensionality check
    
    nDim = obj(ii).nDim;
    if nDim == 0 % is empty
        continue
    end
    
    
    %% create new table if does not exist
    if isempty(obj(ii).features)
        T = table();
        T.trackID = (1:obj(ii).nTracks)'; % only create track IDs
        T.Properties.VariableUnits = {''};
        T.Properties.VariableDescriptions = {'Track ID'};
    else
        if preserveFlag
            % copy existing table
            T = obj(ii).features;
            T = sortrows(T,'trackID'); % code assumes sorted columns
            N = numel(T.Properties.VariableNames);
            if isempty(T.Properties.VariableUnits) % make sure it exists
                T.Properties.VariableUnits = cell(1,N);
            end
            if isempty(T.Properties.VariableDescriptions) % make sure it exists
                T.Properties.VariableDescriptions = cell(1,N);
            end
        else
            % make new
            T = table();
            T.trackID = (1:obj(ii).nTracks)'; % only create track IDs
            T.Properties.VariableUnits = {''};
            T.Properties.VariableDescriptions = {'Track ID'};
        end
    end
    
    if featureNamesFlag
        do = featureNames;
    else
        do = obj(ii).doFeatures;
    end
    
    % clear to-be-computed features
    bool = ismember(T.Properties.VariableNames,do);
    if any(bool)
        T = removevars(T,T.Properties.VariableNames{bool});
    end
    
    waitbar(0.3/nObj+(ii-1)/nObj,wb,waitstr)
    
    %% nTrackPoints
    
    thisFeatures = {'nTrackPoints'};
    
    if any(ismember(thisFeatures,do))
        if ismember('nTrackPoints',do); T.nTrackPoints = obj(ii).nTrackPoints;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Number of Points'; end
    end
    
    
    %% Blinking
    
    thisFeatures = {'nTrackPointsAll','firstFrame','onFraction','blinkingRate',...
        'OFFON','ONOFF'};
    
    if any(ismember(thisFeatures,do))
        if ~obj(ii).onTrace_valid
            obj(ii) = obj(ii).computeOnTrace;
        end
        
        % preallocation
        
        nTrackPointsAll = zeros(obj(ii).nTracks,1);
        firstFrame = zeros(obj(ii).nTracks,1);
        onFraction = zeros(obj(ii).nTracks,1);
        blinkingRate = zeros(obj(ii).nTracks,1);
        OFFON = cell(obj(ii).nTracks,1);
        ONOFF = cell(obj(ii).nTracks,1);
        nTrackPoints = obj(ii).nTrackPoints;
        
        % compute
        for jj = 1:obj(ii).nTracks
            nTrackPointsAll(jj)  =(obj(ii).time{jj}(end) - obj(ii).time{jj}(1) + 1);
            firstFrame(jj) = obj(ii).time{jj}(1); % take first frame of appearance from frames
            onFraction(jj) = nTrackPoints(jj)/nTrackPointsAll(jj);
            nShort = numel(obj(ii).onTrace{jj});
            % weirdly coded - check
            if obj(ii).onTrace{jj}(1) == 0
                nShortZ = nShort - 1; % compensate for zero to indicate ON
            else
                nShortZ = nShort;
            end
            % blinking is ON-OFF-ON
            if mod(nShort,2 == 0) % is even, ends on ON
                blinkingRate(jj) = floor(nShortZ-1/2)/nTrackPointsAll(jj);
            else % is odd, ends on OFF
                blinkingRate(jj) = floor(nShortZ/2)/nTrackPointsAll(jj);
            end
            
            % testing
            oT = obj(ii).onTrace{jj};
            coT = cumsum(oT);
            OFFON{jj} = coT(2:2:end); % the last time with OFF marks the transition
            ONOFF{jj} = coT(1:2:end);
            OFFON{jj} = OFFON{jj}(OFFON{jj} < obj(ii).nFrames);
            ONOFF{jj} = ONOFF{jj}(ONOFF{jj} > 0);
            
        end
        
        % write to table as table
        if ismember('nTrackPointsAll',do); T.nTrackPointsAll = nTrackPointsAll;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Number of all points'; end
        if ismember('firstFrame',do); T.firstFrame = firstFrame;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'First frame'; end
        if ismember('onFraction',do); T.onFraction = onFraction;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'ON-fraction'; end
        if ismember('blinkingRate',do); T.blinkingRate = blinkingRate;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Blinking rate'; end
        if ismember('OFFON',do); obj(ii).OFFON = OFFON; end
        if ismember('ONOFF',do); obj(ii).ONOFF = ONOFF; end
    end
    
    waitbar(0.4/nObj+(ii-1)/nObj,wb,waitstr)
    %% Center
    
    thisFeatures = {'MinBoundCircleRadius','MinBoundCircleCenter',...
        'CenterOfMass','MBCCminusCoM'};
    
    if any(ismember(thisFeatures,do))
        
        % preallocation
        
        MinBoundCircleRadius = zeros(obj(ii).nTracks,1);
        MinBoundCircleCenter = zeros(obj(ii).nTracks,2);
        CenterOfMass = zeros(obj(ii).nTracks,2);
        MBCCminusCoM = zeros(obj(ii).nTracks,1);
        
        % get the minimum bounding circle of the point cloud forming the
        % track (see: smallest-circle problem)
        
        if nDim < 2
            error('Number of dimensions is %i for population %i. MinBoundCircle, CenterOfMass, and MBCCminusCoM cannot be computed.',nDim,ii)
        elseif nDim > 2
            warning('Number of dimensions is %i for population %i. MinBoundCircle, CenterOfMass, and MBCCminusCoM are computed in the XY-plane.',nDim,ii)
        end
        
        for jj = 1:obj(ii).nTracks
            % --- minimum bounding circle 2D
            [C,R] = obj(ii).minboundcircle(obj(ii).coords{jj}(:,1),obj(ii).coords{jj}(:,2),1);
            MinBoundCircleRadius(jj) = R;
            MinBoundCircleCenter(jj,:) = C;
            % --- center of mass (CoM)
            CenterOfMass(jj,:) = mean(obj(ii).coords{jj}(:,1:2),1); % take mean over all rows
            % get difference between CoM and center of minimum bounding circle
            % in percent; idea: if this difference is small the track is a
            % 'ball of wool', if the difference is large it is an 'uncoiled ball of wool':
            % in percent of MinBoundCircleRadius (i.e. between 0 to 1)
            MBCCminusCoM(jj) = pdist2(MinBoundCircleCenter(jj,1:2),CenterOfMass(jj,1:2))/R; % default value for pdist2 is Eucledian
        end
        
        % write to table as table
        if ismember('MinBoundCircleRadius',do); T.MinBoundCircleRadius = MinBoundCircleRadius;...
                T.Properties.VariableUnits{end} = 'pixelsize'; T.Properties.VariableDescriptions{end} = 'MBC radius'; end
        if ismember('MinBoundCircleCenter',do); T.MinBoundCircleCenter = MinBoundCircleCenter;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'MCB center'; end
        if ismember('CenterOfMass',do); T.CenterOfMass = CenterOfMass;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Center of mass'; end
        if ismember('MBCCminusCoM',do); T.MBCCminusCoM = MBCCminusCoM;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'MBCC minus CoM'; end
        
    end
    
    waitbar(0.6/nObj+(ii-1)/nObj,wb,waitstr)
    %% Entropy
    
    thisFeatures = {'entropy'};
    
    if any(ismember(thisFeatures,do))
        
        if nDim < 2
            error('Number of dimensions is %i for population %i. Entropy cannot be computed.',nDim,ii)
        elseif nDim > 2
            warning('Number of dimensions is %i for population %i. Entropy is computed in the XY-plane.',nDim,ii)
        end
        
        if ~ismember('MinBoundCircleRadius',do) || ~ismember('MinBoundCircleCenter',do)
            objc = obj(ii).computeFeatures(false,{'MinBoundCircleRadius','MinBoundCircleCenter'});
            R = objc.features.MinBoundCircleRadius;
            CC = objc.features.MinBoundCircleCenter;
            clear objc
        else
            % use embounding circle radius as width and height
            R = T.MinBoundCircleRadius;
            CC = T.MinBoundCircleCenter;
        end
        
        % preallocation
        entropy_ = zeros(obj(ii).nTracks,1);
        
        for jj = 1:obj(ii).nTracks
            % idea: make it a binary image and get the entropy of the track
            imgSize = (round(R(jj)*obj(ii).resolution)*2)+2;
            I = zeros(imgSize,imgSize);
            % now place the points of the track in this image:
            
            relCoordX = uint8(round((obj(ii).coords{jj}(:,1) - (CC(jj,1)-R(jj)))*obj(ii).resolution)+1);
            relCoordY = uint8(round((obj(ii).coords{jj}(:,2) - (CC(jj,2)-R(jj)))*obj(ii).resolution)+1);
            I(sub2ind(size(I),relCoordY,relCoordX)) = 1; % avoid for loop
            entropy_(jj) = entropy(uint8(I));
        end
        
        % write to table as table
        if ismember('entropy',do); T.entropy = entropy_;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Entropy'; end
        
    end
    
    waitbar(0.8/nObj+(ii-1)/nObj,wb,waitstr)
    %% displacement
    
    thisFeatures = {'trackLength','tortuosity'};
    
    if any(ismember(thisFeatures,do))
        
        % preallocation
        
        trackLength = NaN(obj(ii).nTracks,1);
        tortuosity = NaN(obj(ii).nTracks,1);
        
        % compute
        for jj = 1:obj(ii).nTracks
            trackLength(jj) = sum(sqrt(sum((obj(ii).coords{jj}(1:end-1,:)-obj(ii).coords{jj}(2:end,:)).^2,2)));
            tortuosity(jj) = trackLength(jj) / sqrt(sum((obj(ii).coords{jj}(1,:)-obj(ii).coords{jj}(end,:)).^2));
        end
        
        % write to table as table
        if ismember('trackLength',do); T.trackLength = trackLength;...
                T.Properties.VariableUnits{end} = 'pixelsize'; T.Properties.VariableDescriptions{end} = 'Length'; end
        if ismember('tortuosity',do); T.tortuosity = tortuosity;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Tortuosity'; end
        
    end
    
    waitbar(0.9/nObj+(ii-1)/nObj,wb,waitstr)
    %% PCA
    
    thisFeatures = {'EVec1','EVec2','CVE1','CVE2','EV1angle'};
    
    if any(ismember(thisFeatures,do))
        
        if nDim < 2
            error('Number of dimensions is %i for population %i. Elongation (angle) cannot be computed.',nDim,ii)
        elseif nDim > 2
            warning('Number of dimensions is %i for population %i. Elongation (angle) is computed in the XY-plane.',nDim,ii)
        end
        
        % preallocation
        
        EVec1 = NaN(obj(ii).nTracks,2);
        EVec2 = NaN(obj(ii).nTracks,2);
        CVE1 = NaN(obj(ii).nTracks,1);
        CVE2 = NaN(obj(ii).nTracks,1);
        EV1angle = NaN(obj(ii).nTracks,1);
        
        % get the principal component of each track (direction of largest variance)
        for jj = 1:obj(ii).nTracks
            covM = cov(obj(ii).coords{jj}(:,1:2));
            [V,D]=eig(covM);
            
            % correct all EVecs to be in positive half space:
            % => if the y-part of EVec1 is negative, flip it by 180deg
            if V(2,1)<0
                V(2,1)=-V(2,1);
                V(1,1)=-V(1,1);
            end
            if V(2,2)<0
                V(2,2)=-V(2,2);
                V(1,2)=-V(1,2);
            end
            CVE(1,1) = D(1,1)/(D(1,1)+D(2,2));
            CVE(2,1) = D(2,2)/(D(1,1)+D(2,2));
            if D(1,1)>D(2,2)
                EVec1(jj,:) = V(:,1);
                EVec2(jj,:) = V(:,2);
                CVE1(jj) = CVE(1,1);
                CVE2(jj) = CVE(2,1);
            else
                EVec1(jj,:) = V(:,2);
                EVec2(jj,:) = V(:,1);
                CVE1(jj) = CVE(2,1);
                CVE2(jj) = CVE(1,1);
            end
            % get angle in degrees for each EVec1 (per definition:
            % angle between basis set vector [1 0] and EVec1)
            EV1angle(jj) = acosd(dot(EVec1(jj,:),[1 0])); % both vectors are of length 1
            
        end
        
        % compute normalized CVE1, which ranges from 0 to 1. In the 2D case,
        % CVE1 can be 1 (all points are on a line) or 0.5 (spherical cloud of
        % points). This scales 0.5-1 to 0-1.
        wEV1 = CVE1.*2 - 1;
        
        % write to table as table
        if ismember('EVec1',do); T.EVec1 = EVec1;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'EVec1'; end
        if ismember('EVec2',do); T.EVec2 = EVec2;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'EVec2'; end
        if ismember('CVE1',do); T.CVE1 = CVE1;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'CVE1'; end
        if ismember('CVE2',do); T.CVE2 = CVE2;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'CVE2'; end
        if ismember('wEV1',do); T.wEV1 = wEV1;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Elongation'; end
        if ismember('EV1angle',do); T.EV1angle = EV1angle;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Elongation angle'; end
        
    end
    
    waitbar(1/nObj+(ii-1)/nObj,wb,waitstr)
    %% Voronoi
    
    thisFeatures = {'voronoiSA'};
    
    if any(ismember(thisFeatures,do))
        
        if nDim < 2
            error('Number of dimensions is %i for population %i. Voronoi area cannot be computed.',nDim,ii)
        elseif nDim > 2
            warning('Number of dimensions is %i for population %i. Voronoi area is computed in the XY-plane.',nDim,ii)
        end
        
        if ~ismember('CenterOfMass',do)
            objc = obj(ii).computeFeatures(false,{'CenterOfMass'});
            COM = objc.features.CenterOfMass;
            clear objc
        else
            COM = T.CenterOfMass;
        end
        
        % preallocation
        voronoiSA = zeros(size(COM,1),1);
        
        % get the convex hull and extend it by 5 pixels; then add the hull points
        % to the data to cut of the voronoi diagram at the extended hull
        hPs = obj(ii).getHullPoints(COM,10);
        startOfhPs = size(COM,1)+1;
        COM = vertcat(COM,hPs);
        [v,c]=voronoin(COM);
        
        for jj = 1:length(c)
            if all(c{jj}~=1)  % If at least one of the indices is 1,
                % then it is an open region and we can't
                % patch that. These are introduced as hull points
                % to create a convex hull in the voronoi diagram
                if jj<startOfhPs % only patch real data points
                    voronoiSA(jj) = polyarea(v(c{jj},1),v(c{jj},2));
                end
            end
        end
        
        % write to table as table
        if ismember('voronoiSA',do); T.voronoiSA = voronoiSA;...
                T.Properties.VariableUnits{end} = 'pixelsize.^2'; T.Properties.VariableDescriptions{end} = 'Voronoi area'; end
        
    end
    
    %% Nearest Neighbor analysis
    
    thisFeatures = {'DensityOfPointsK1','DensityOfPointsR'};
    
    if any(ismember(thisFeatures,do))
        
        if nDim < 2
            error('Number of dimensions is %i for population %i. property cannot be computed.',nDim,ii)
        end
        
        %DensityOfPointsK1
        DensityOfPointsK1 = zeros(obj(ii).nTracks,1); %Initializing Variable
        if ~ismember('CenterOfMass',do)
            objc = obj(ii).computeFeatures(false,{'CenterOfMass'});
            COM = objc.features.CenterOfMass;
            clear objc
        else
            %COM = obj(ii).features.CenterOfMass;
            COM = T.CenterOfMass;
        end
        
        [~,r] = knnsearch(COM,COM,'K',2); %distance to the 2 closest track centers of mass (1st is itself)
        DensityOfPointsK1 = 1./(pi*r(:,2).^2); %density= 1/circle area where r is distance to the closest track center of mass
        DensityOfPointsK1=DensityOfPointsK1/obj(ii).nTracks; %divide point density by Number of Tracks
        %DensityOfPointsK1=DensityOfPointsK1/mean(DensityOfPointsK1); %divide point density by mean point density
        
        
        
        %DensityOfPointsR density of points within circle with radius R normalized by number of tracks
        R=obj(ii).DisplacementThreshold; % R is= Displacement Threshold
        AllPoints = COM;
        DensityOfPointsR = zeros(obj(ii).nTracks,1); %Initializing Variable
        for idx=1:length(AllPoints)
            Distances = sqrt( sum( (AllPoints-AllPoints(idx,:)).^2 ,2) );
            Ninside   = length( find(Distances<=R) );% number of points within each circle with radius R
            DensityOfPointsR(idx) = Ninside/(pi*R.^2);
        end
        DensityOfPointsR=DensityOfPointsR/obj(ii).nTracks; %Normalize density with number of Tracks
        %DensityOfPointsR=DensityOfPointsR/mean(DensityOfPointsR); %Normalize density with number of Tracks
        
        
        % T.CenterOfMass=COM;
        if ismember('DensityOfPointsK1',do); T.DensityOfPointsK1 = DensityOfPointsK1;...
                T.Properties.VariableUnits{end} = '1/pixelsize^2'; T.Properties.VariableDescriptions{end} = 'DensityOfPointsK1'; end
        if ismember('DensityOfPointsR',do); T.DensityOfPointsR = DensityOfPointsR;...
                T.Properties.VariableUnits{end} = '1/pixelsize^2'; T.Properties.VariableDescriptions{end} = 'DensityOfPointsR'; end
        
    end
    
    
    %% Immobile Seps fraction
    
    thisFeatures = {'FractionOfImmobileSteps','NImmobileSteps'};
    
    if any(ismember(thisFeatures,do))
        
        if nDim < 2
            error('Number of dimensions is %i for population %i. property cannot be computed.',nDim,ii)
        end
        
        
        % get squared displacements of individual tracks
        for jj=1:obj(ii).nTracks
            dummyCoords=obj(ii).coords{jj, 1};
            dR2{ii}{jj,1}=(diff(dummyCoords(:,1))).^2+(diff(dummyCoords(:,2))).^2;
        end
        
        for jj=1:obj(1,ii).nTracks
            dummy=dR2{ii}{jj,1};
            NImmobileSteps{ii}{jj,1}=numel(find(dummy<=obj(ii).DisplacementThreshold)); %number of steps below displacement trehesold
            FractionOfImmobileSteps{ii}{jj,1}=NImmobileSteps{ii}{jj,1}/numel(dummy)*100; %fraction (%) of steps below displacement threshold
        end
        
        
        if ismember('FractionOfImmobileSteps',do); T.FractionOfImmobileSteps = cell2mat(FractionOfImmobileSteps{1,ii});...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'FractionOfImmobileSteps'; end
        if ismember('NImmobileSteps',do); T.NImmobileSteps = cell2mat((NImmobileSteps{1,ii}));...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'NImmobileSteps'; end
    end
    
    %% Adsorption Time Analysis
    thisFeatures = {'MeanAdsorptionTime','MeanFlightTime','AdsorptionFligthTimeRatio'};
    if any(ismember(thisFeatures,do))
        if nDim < 2
            error('Number of dimensions is %i for population %i. property cannot be computed.',nDim,ii)
        end
        
        
        % get squared displacements of individual tracks
        for jj=1:obj(ii).nTracks
            dummyCoords=obj(ii).coords{jj, 1};
            dR2{ii}{jj,1}=(diff(dummyCoords(:,1))).^2+(diff(dummyCoords(:,2))).^2;
        end
        
        % Adsorption times
        for jj=1:obj(ii).nTracks
            dummy=dR2{ii}{jj,1};
            ImmobileStepsID{ii}{jj,1}=(find(dummy<=obj(ii).DisplacementThreshold)); %find displacements smaller than Disp Threshold
            ImmobileStepsTime{ii}{jj,1}=obj(1, ii).time{jj, 1}(ImmobileStepsID{1, ii}{jj, 1}); %times where particle is adsorbed
            dummy_sequence=ImmobileStepsID{ii}{jj,1};
            
            if isempty(dummy_sequence) % if no adsorption steps, no adsorption ID
                AdsorptionEventsID{ii}{jj,1}=[];
            else
                
                dummy_sequence(end+1) = -1; % Add a new isolating point
                D=diff(dummy_sequence); %calculate difference of immobile steps index (consecutive adsorption steps have a difference of 1)
                I_1 = find(D ~= 1); % Find indexes of isolating points
                [m,n] = size(I_1);
                Start_Idx = 1 ; % Set start index
                for tt = 1:m %for all adsorption events
                    End_Idx = I_1(tt); % Set end index
                    AdsorptionEventsID{ii}{jj,1}{tt} =  dummy_sequence(Start_Idx:End_Idx); % Find consecuative sequences
                    Start_Idx = End_Idx + 1;
                    % update start index for the next consecuitive sequence
                end
            end
            
        end
        
        FrameTime=obj(ii).getUnitFactor('dt'); %Frame time in in s
        for jj=1:obj(ii).nTracks
            % If no adsorption events Adsorption begining, end, Id and duration
            % are set empty
            if isempty(AdsorptionEventsID{ii}{jj,1})
                AdsoprtionDuration{ii}{jj,1}=[];
                AdsorptionEventsBeginingID{ii}{jj,1}=[];
                AdsorptionEventsEndID{ii}{jj,1}=[];
                AdsoprtionBeginingTime{ii}{jj,1}=[];
                AdsoprtionEndTime{ii}{jj,1}=[];
            else
                % for each adsorption event tt find begining id (minmum) end id
                % (max) and get times based on ids as well as duration of event
                % (difeerence begining and end times )
                for tt= 1:numel(AdsorptionEventsID{1, ii}{jj, 1})
                    AdsorptionEventsBeginingID{ii}{jj,1}{tt}= min(AdsorptionEventsID{1, ii}{jj, 1}{ tt});
                    AdsorptionEventsEndID{ii}{jj,1}{tt}= max(AdsorptionEventsID{1, ii}{jj, 1}{ tt});
                    AdsoprtionBeginingTime{ii}{jj,1}{tt}=obj(1, ii).time{jj, 1}(AdsorptionEventsBeginingID{ii}{jj,1}{tt});
                    AdsoprtionEndTime{ii}{jj,1}{tt}=obj(1, ii).time{jj, 1}(AdsorptionEventsEndID{ii}{jj,1}{tt});
                    if length(AdsorptionEventsID{ii}{jj,1}{tt})==1 %if adsorption event is only one frame long set duration=Framentime
                        AdsoprtionDuration{ii}{jj,1}{tt}=FrameTime;
                    else
                        AdsoprtionDuration{ii}{jj,1}{tt}=(AdsoprtionEndTime{ii}{jj,1}{tt}-AdsoprtionBeginingTime{ii}{jj,1}{tt})*FrameTime;
                    end
                end
            end
        end
        
        % get flight Times (time between two adsorption steps)
        for jj=1:obj(ii).nTracks
            if length(AdsoprtionDuration{ii}{jj,1})<2 % if there 1 adosrption event or less -> no flight time
                FlightTimes{ii}{jj,1}=[];
            else
                for tt=2:length(AdsoprtionDuration{1, ii}{jj, 1})
                    %time between adsorption end i and adsorption begining i+1
                    FlightTimes{ii}{jj,1}{tt-1}=(AdsoprtionBeginingTime{ii}{jj,1}{tt}-AdsoprtionEndTime{ii}{jj,1}{tt-1})*FrameTime;
                end
            end
        end
        
        %mean flight times of each track
        dummy_Flighttime=[];
        for jj=1:obj(ii).nTracks
            MeanFlightTimes{ii}(jj)=mean(cell2mat(FlightTimes{ii}{jj,1}));% mean flight times of each track
            for tt=1:length(FlightTimes{1, ii}{jj, 1})
                dummy_Flighttime=[dummy_Flighttime;FlightTimes{ii}{jj,1}{tt}];
                FlightTimesAllTracks{ii}=dummy_Flighttime;
            end
        end
        MeanFlightTimeAllTracks(ii)=mean(FlightTimesAllTracks{ii}); %mean flight time of all events and all tracks, not necessary here
        
        
        
        % Adsorption time correction - ignoring adsorption events at the begining/end of track
        %this is not necessary for the flight times as only flights between two
        %adsorption steps are considered
        AdsoprtionDurationUncorrected=AdsoprtionDuration;
        
        for jj=1:obj(ii).nTracks
            if not(isempty(AdsoprtionDuration{ii}{jj,1}))
                %if adsorption event beginning time is equal to track beginning time
                %-> adsorption duration not considered -> []
                if AdsoprtionBeginingTime{ii}{jj,1}{1,1}==obj(ii).time{jj, 1}(1)
                    AdsoprtionDuration{ii}{jj,1}{1,1}=[];
                end
                
                %if adsorption event end time is equal to track end time
                %-> adsorption duration not considered -> []
                if AdsoprtionEndTime{ii}{jj,1}{1,end}==obj(ii).time{jj, 1}(end)%AdsorptionTimesAllTracks
                    AdsoprtionDuration{ii}{jj,1}{1,end}=[];
                end
            end
        end
        
        %mean adsorption times of each track
        AdsorptionTimes=AdsoprtionDuration;
        dummy_adsorptiontime=[];
        for jj=1:obj(ii).nTracks
            MeanAdsorptionTimes{ii}(jj)=mean(cell2mat(AdsoprtionDuration{ii}{jj,1})); %mean adsorption times of each track
            for tt=1:length(AdsoprtionDuration{1, ii}{jj, 1})
                dummy_adsorptiontime=[dummy_adsorptiontime;AdsoprtionDuration{ii}{jj,1}{tt}];
                AdsorptionTimesAllTracks{ii}=dummy_adsorptiontime;
            end
            
        end
        MeanAdsorptionTimeAllTracks(ii)=mean(AdsorptionTimesAllTracks{ii}); %mean adsorption time of all events and all tracks, not necessary here
        
        
        for jj=1:obj(ii).nTracks
            %mean adsorption flight time ratio of each track
            AdsFlightRatio{ii}(jj)=MeanAdsorptionTimes{ii}(jj)./MeanFlightTimes{ii}(jj);
        end
        
        %convert proerties to mat formsat to fit in table
        MeanAdsorptionTime=((MeanAdsorptionTimes{1,ii}))';
        MeanFlightTime=((MeanFlightTimes{1,ii}))';
        AdsorptionFligthTimeRatio=((AdsFlightRatio{1,ii}))';
        
        
        
        if ismember('MeanAdsorptionTime',do); T.MeanAdsorptionTime = MeanAdsorptionTime;...
                T.Properties.VariableUnits{end} = 's'; T.Properties.VariableDescriptions{end} = 'MeanAdsorptionTime'; end
        if ismember('MeanFlightTime',do); T.MeanFlightTime = MeanFlightTime;...
                T.Properties.VariableUnits{end} = 's'; T.Properties.VariableDescriptions{end} = 'MeanFlightTime'; end
        if ismember('AdsorptionFligthTimeRatio',do); T.AdsorptionFligthTimeRatio = AdsorptionFligthTimeRatio;...
                T.Properties.VariableUnits{end} = 's'; T.Properties.VariableDescriptions{end} = 'AdsorptionFligthTimeRatio'; end
        
    end
    
    %% Distance to closet adsorbed particle
    
    thisFeatures = {'DistanceToClosestAdsorbedParticle'};
    if any(ismember(thisFeatures,do))
        if nDim < 2
            error('Number of dimensions is %i for population %i. property cannot be computed.',nDim,ii)
        end
        
        
        % How many videos in dataset?
        Datapaths=unique(obj(ii).TrackPath); %names of unique dataset paths
        DistanceToClosestAdsorbedParticle=NaN(obj(ii).nTracks,1); %initialize property
        %find tracks that belong to the same path
        
        if numel(Datapaths) ==1 %if there is only one video
            FullyImmobileIndex=obj(ii).features.FractionOfImmobileSteps==100; %find fully immobile tracks
            try % segment function doesn't work if there is nothing to segment
                objImmobile=obj(ii).segment(FullyImmobileIndex); %segment immobile tracks
                FullyAdsorbedCoords=objImmobile.features.CenterOfMass; %get centers of mass of immobile tracks
                COM_Tracks=obj(ii).features.CenterOfMass; % get centers of mass of other tracks
                [~,r] = knnsearch(FullyAdsorbedCoords,COM_Tracks,'K',1); % for each track find the closest immobile track
            catch % if there are no immobile tracks set r=NaN
                r=NaN(obj(ii).nTracks,1);
            end
            DistanceToClosestAdsorbedParticle=r;
        else % if the number of paths (videos) is >1 consider only adsorbed particles within the same video
            for v=1:numel(Datapaths)
                SamevideoIndex=strcmp(obj(ii).TrackPath,Datapaths(v)); % 1 if track belongs to same video, otherwise 0
                objSameVideo=obj(ii).segment(SamevideoIndex); %segment tracks that belong to the vth path
                FullyImmobileIndex=objSameVideo.features.FractionOfImmobileSteps==100; %find immobile tracks within video
                try % if there is nothing to segment there is an error
                    objSameVideoImmobile=objSameVideo.segment(FullyImmobileIndex); %segment tracks that belong to the same path
                    FullyAdsorbedCoords=objSameVideoImmobile.features.CenterOfMass; %coords of immonbile Track within the video
                    COM_Tracks=objSameVideo.features.CenterOfMass; %coords of tracks withing the video
                    [~,r] = knnsearch(FullyAdsorbedCoords,COM_Tracks,'K',1); %% for each track find the closest immobile track within the video
                catch %if there is nothing to segment (no immobile track) DistanceToClosestAdsorbedParticle=NaN
                    r=NaN(objSameVideo.nTracks,1);
                end
                DistanceToClosestAdsorbedParticle(SamevideoIndex)=r;
            end
        end
        
        if ismember('DistanceToClosestAdsorbedParticle',do); T.DistanceToClosestAdsorbedParticle = DistanceToClosestAdsorbedParticle;...
                T.Properties.VariableUnits{end} = 'Pixels'; T.Properties.VariableDescriptions{end} = 'MeanAdsorptionTime'; end
        
    end
    
    % end
    %% Save
    
    % store features
    obj(ii).features = T;
    % set validity
    obj(ii).features_valid = true;
    
end

close(wb)

end

%Template:
%
% thisFeatures = {};
%
% if any(ismember(thisFeatures,do))
%
%     % preallocation
%
%
%
%     % write to table as table
%     if ismember('clear',do); T.nTrackPointsAll = nTrackPointsAll; end
%
% end

