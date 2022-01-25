function [AdsorptionTimes,MeanAdsorptionTimes,FlightTimes,MeanFlightTimes,AdsFlightRatio,AdsorptionTimesAllTracks,FlightTimesAllTracks]=getAdsorptionTimes(obj);

for ii=1:numel(obj)    
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
        %MeanAdsorptionTime=((MeanAdsorptionTimes{1,ii}))';
        %MeanFlightTime=((MeanFlightTimes{1,ii}))';
        %AdsorptionFligthTimeRatio=((AdsFlightRatio{1,ii}))';
        
end

%% plotting ads times 
%{
if PlotFlag==1
    figure   %comparison fraction of immobile steps
    dummylegend=[];
    for ii = 1:numel(obj)
        %histogram(AdsorptionTimesAllTracks{ii},'Normalization','probability','BinWidth',FrameTime)
        histogram(AdsorptionTimesAllTracks{ii},'BinWidth',FrameTime,'Normalization','probability')
        hold on
        dummylegend=[dummylegend;ii];
    end
    xlabel ('Adsorption times (s)')
    ylabel ('probability (-)')
    %ylabel ('probability (-)')
    legend(num2str(dummylegend))
    set(gca,'YScale','log')
    
    
        figure   %comparison fraction of immobile steps
    dummylegend=[];
    for ii = 1:numel(obj)
        %histogram(MeanAdsorptionTimes{ii},'Normalization','probability','BinWidth',FrameTime)\
        histogram(MeanAdsorptionTimes{ii},'BinWidth',FrameTime,'Normalization','probability')
        hold on
        dummylegend=[dummylegend;ii];
    end
    xlabel ('mean Adsorption times (s)')
    ylabel ('probability (-)')
    %ylabel ('probability (-)')
    legend(num2str(dummylegend))
    set(gca,'YScale','log')
    
    figure
    bar(dummylegend,diag(MeanAdsorptionTimeAllTracks),'stacked')
    ylabel ('mean Adsorption times all tracks (s)')
    
    
    %% flight times
     figure   %comparison fraction of immobile steps
    dummylegend=[];
    for ii = 1:numel(obj)
       
        histogram(FlightTimesAllTracks{ii},'BinWidth',FrameTime,'Normalization','probability')
        hold on
        dummylegend=[dummylegend;ii];
    end
    xlabel ('Flight times (s)')
    ylabel ('probability (-)')
    %ylabel ('probability (-)')
    legend(num2str(dummylegend))
    set(gca,'YScale','log')
    
    
        figure   %comparison fraction of immobile steps
    dummylegend=[];
    for ii = 1:numel(obj)
        %histogram(MeanFlyTimes{ii},'Normalization','probability','BinWidth',FrameTime)\
        histogram(MeanFlightTimes{ii},'BinWidth',FrameTime,'Normalization','probability')
        hold on
        dummylegend=[dummylegend;ii];
    end
    xlabel ('mean flight times (s)')
    ylabel ('probability (-)')
    %ylabel ('probability (-)')
    legend(num2str(dummylegend))
    set(gca,'YScale','log')
    
    figure
    bar(dummylegend,diag(MeanFlightTimeAllTracks),'stacked')
    ylabel ('mean flight times all tracks (s)')
    
    figure
    bar(dummylegend,diag(MeanFlightTimeAllTracks),'stacked')
    ylabel ('mean flight times all tracks (s)')
    
     %% flight vs Ads times 
    
     figure
     for ii = 1:numel(obj)
       
         histogram(AdsFlightRatio{ii},'BinWidth',FrameTime*10,'Normalization','probability')
         hold on
         dummylegend=[dummylegend;ii];
     end
     xlabel ('Adsorption/Flight time ratio (-)')
     ylabel ('probability')
     legend(num2str(dummylegend))
     set(gca,'YScale','log')   
end
%}

end