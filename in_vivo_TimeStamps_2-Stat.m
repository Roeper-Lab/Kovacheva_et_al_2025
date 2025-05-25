%% Analysis of all periods for all neurons
% By Roeper's Lab

%% Import and align data

clear all                                           % clear the working space
ID = 'OHDA_late';                                     % assign the same ID as the GroupName from the 'TimeStamps_Import_Splited'
% 'ACSFmed_early'
cd(['ISI_' ID])                                     % goes to the previosly created for analysis folder 
addpath(pwd);                                       % adds the current path to matlab's search folders
DataFile = ['ISI_' ID '.mat'];                      % prepares the name of matlab's file with the imported timestamps

% class = 1;

load(DataFile);                                     % loads the data

lgLabels = {'Bursts', 'Pauses'};                    % variables needed for labeling of legends in graphs
scale = 's';                                        % addapts the scale to s (if needed change to ms)
% OutlierCriterion = 3;

Plot = 'yes1';                                       % for some graphs you can skip plotting them, so that the script runs faster

if strcmp(scale, 's') == 1;                         % adapting the bins to the time scale
    binsize = 0.01;
    xvalues = 0:binsize:60;
    xadd = 0.2;
    refpd = 0.005;
elseif strcmp(scale, 'ms') == 1;
    xvalues = 0:5:20^3;
    xadd = 200;
    refpd = 2.5;
end
   
    
% Set folders
mkdir(['Stat_' char(GroupName)]);                   % creates a new directory with for the analysis
oldFolder = cd(['Stat_' char(GroupName)]);          % Sets the return to old directory
addpath(pwd);                                       % adds the current path to matlab's search folders

%% Extract mean and median frequency, 25 and 75 percentiles, IQR and CV 

for c = 1:nCells
    Data(c).(GroupName).MeanFreq = mean(1./Data(c).(GroupName).ISI);             % Mean frequency
    Data(c).(GroupName).MedianFreq = median(1./Data(c).(GroupName).ISI);         % Median frequency
    
    Data(c).(GroupName).SDM = std(Data(c).(GroupName).ISI);                  % Standard deviation of mean
    Data(c).(GroupName).CV = Data(c).(GroupName).SDM/mean(Data(c).(GroupName).ISI);           % Coefficient of variation
    Data(c).(GroupName).CV2 = Data(c).(GroupName).CV^2;                                                    % CV squared
    
    Data(c).(GroupName).Q1 = prctile(Data(c).(GroupName).ISI,25);               % 25th percentile
    Data(c).(GroupName).Q3 = prctile(Data(c).(GroupName).ISI,75);               % 75th percentile
    Data(c).(GroupName).IQR = Data(c).(GroupName).Q3 - Data(c).(GroupName).Q1;    % 75th - 25th percentile (interquartile range)
    
    %% Determine outlier thresholds and extract upper and lower outlier ISIs
    % Control
    Data(c).(GroupName).UpStatThres = Data(c).(GroupName).Q3 + 3*(Data(c).(GroupName).IQR);    % Upper outilier threshold
    Data(c).(GroupName).LowStatThres = Data(c).(GroupName).Q1 - 3*(Data(c).(GroupName).IQR);    % Lower outilier threshold
    
    ISI = Data(c).(GroupName).ISI;
    
    UpOutIdx = find (ISI > Data(c).(GroupName).UpStatThres);         % Get upper outlier indices
    
    if isempty(UpOutIdx);
        Data(c).(GroupName).UpOutISI = 0;                    % Get upper outlier ISIs
    else
        Data(c).(GroupName).UpOutISI = ISI(UpOutIdx);
    end
    Data(c).(GroupName).UpOutIDx = UpOutIdx;
    
    LowOutIdx = find (ISI < Data(c).(GroupName).LowStatThres);       % Get lower outlier indices
    
    if isempty(LowOutIdx);
        Data(c).(GroupName).LowOutISI = 0;                  % Get lower outlier ISIs
    else
        Data(c).(GroupName).LowOutISI = ISI(LowOutIdx);
    end
    Data(c).(GroupName).LowOutIDx = LowOutIdx;
end



%% Find bursts
for c = 1:nCells
    % Set burst thresholds (Onset: 80 ms, Offset: 160 ms; Grace and Bunny 1984)
    if strcmp(scale,'ms');
        start_ISI= 80;
        end_ISI= 160;
    elseif strcmp(scale,'s');
        start_ISI= 0.08;
        end_ISI= 0.16;
    end
    
    ISI = Data(c).(GroupName).ISI;
    
    in_burst = false; % 'false' = not in burst; 'true' = in burst
    burst_count = 0; % counts the ISIs in a burst-event
    burstEventNidx = 0; % holds burst-as-a-whole-event count
    inISI = []; % vector with the ISIs from one burst-event for intra burst analysis
    burstN = []; % vector with number of ISIs per burst for eachs single burst-event
    
    intra_burst_freq_MEAN = 0; % vector with the mean freq for each burst-event
    intra_burst_freq_MAX = 0; % vector with the max freq for each burst-event
    burstDur = []; % vector with the burst durations of each and every burst-event
    burst_idx = []; % vector with the indeces of all ISIs in all bursts
    BurstISI=[]; % vector with the actual burst-ISIs from all bursts
    % burst_idx should have the same length as BurstISI
    preBrISI = []; postBrISI = []; % vectors for the pre- and post-burst-ISIs
    preBrISI_idx = []; postBrISI_idx = []; % vectors for the pre- and post-burst-ISIs
    twoSpBr = 0; threeSpBr = 0; mSpBr = 0; % counter for the burst-spikes
    spikesINburst = []; % vector with the number of spikes for each single burst-event
    
    for i=1:length(ISI);
        if  ~in_burst && ISI(i) < start_ISI; % the burst starts: the first ISI is less than 80ms (e.g.) and the in_burst is false
            burst_idx = [burst_idx i];  % saves the index of the burst ISI
            BurstISI = [BurstISI ISI(i)]; % Save actual burst ISI
            burst_count = burst_count+1; % frist 1 ISI in a burst
            inISI = [inISI ISI(i)]; % saves the ISIs in a burst for intra burst analysis
            in_burst = true; % declares that a burst has started
            if i > 1
                preBrISI = [preBrISI ISI(i-1)]; % in case that the spike train do not start with a burst, it saves the pre-burst ISI
                preBrISI_idx = [preBrISI_idx (i-1)];
            end
        elseif in_burst && ISI(i) < end_ISI % the burst continues: the next ISI is less than 160ms (e.g.) and the in_burst is true
            burst_idx = [burst_idx i]; % saves the next index
            BurstISI = [BurstISI ISI(i)]; % adds the next burst ISI
            burst_count = burst_count+1; % increases the burst with one more ISI
            inISI = [inISI  ISI(i)]; % adds the next ISI to the vector for the current burst
            in_burst = true; % declares that the burst is still ongoing
        elseif in_burst && ISI(i) >= end_ISI % the burst ends: the current ISI is bigger than 160ms (e.g.) and in_burst is true
               
                burstEventNidx = burstEventNidx+1; % icreases burst-as-a-whole-event with one 
                burstN = [burstN burst_count]; % Saves number of ISIs per burst for eachs single burst-event
                intra_burst_freq_MEAN(burstEventNidx) = mean(1./inISI); % Intra-burst frequency for each single burst-event - MEAN
                intra_burst_freq_MAX(burstEventNidx) = 1/min(inISI); % Intra-burst frequency for each single burst-event - MAX
                spikesINburst = [spikesINburst burst_count+1]; 
                burstDur = [burstDur sum(inISI)];
                
                % control check 
                if burst_count ~= length(inISI)
                    error (['Error in burst analysis:' int2str(c)]);
                end
                if length(burst_idx) ~= length(BurstISI)
                    error (['Error in burst analysis:' int2str(c)]);
                end
               
                % post-burst ISI data 
                if i < length(ISI) % in case that the whole spiketrain was not finished
                    postBrISI =  [postBrISI  ISI(i)]; % save the current ISI as a post-burst ISI
                    postBrISI_idx = [postBrISI_idx i];
                end
                
                % counting the different types of bursts
                if burst_count == 1 % in case of a 2-spiked burst
                    twoSpBr = twoSpBr + 1; % increase the counter for the 2-spiked burst
                elseif burst_count == 2 % in case of a 3-spiked burst
                    threeSpBr = threeSpBr + 1; % icrease the counter for the 3-spiked burst
                else mSpBr = mSpBr + 1; % in other case - increase the counter for the many-spiked burst
                end
                
                % zeroing the values
                in_burst = false;  % the burst has finished => in_burst is false/negative
                inISI = []; % empties the intra-burst ISIs for calculating the intra-burst frequency for each single burst-event
                burst_count=0; % restart the burst-ISI-counter 
                
        elseif ~in_burst && ISI(i) >= start_ISI % not a burst: the in_burst is negative/false and the ISI is bigger than 80ms
                
                in_burst=false; burst_count=0; % do not save this ISI, in_burst stases false and the burst-spike-counter is still 0
                
        end
    end  % all the ISIs for this period are analysed => save the analysis to the structure    
    
    if i == length(ISI) && burst_count >=1 && in_burst % in case that all the spikes are analysed and at the end there was an unfinished burst-event
        burstN = [burstN burst_count]; % adds the last unfinished-burst to the burst-event-spike-counter
        spikesINburst = [spikesINburst burst_count+1];
        intra_burst_freq_MEAN = [intra_burst_freq_MEAN mean(1./inISI)]; % Intra-burst frequency for each single burst-event - MEAN
        intra_burst_freq_MAX = [intra_burst_freq_MAX 1/min(inISI)]; % Intra-burst frequency for each single burst-event - MAX
        burstDur = [burstDur sum(inISI)];
        burstEventNidx = burstEventNidx + 1;
    end
 
    SFB = (sum(spikesINburst)*100)/(length(ISI)+1); % calculate the relative spikes fired in burst
    
    if isempty (BurstISI)   % in case that the cell/period had no bursts
        BurstISI = 0;
        SFB = 0; 
        burst_idx = 0;
    end
    
    % transering hte burst analysis to the struct
    
    Data(c).(GroupName).BurstISI = BurstISI; % vector with the  burst-ISIs from all burst-events
    Data(c).(GroupName).BurstStats.SFB = SFB;
    Data(c).(GroupName).BurstStats.indices = burst_idx;
    
    % if SFB>0
        Data(c).(GroupName).BurstStats.NIntra = burstN; % vector with number of ISIs per burst for eachs single burst-event
        Data(c).(GroupName).BurstStats.NspikesB = spikesINburst; % vector with number of spikes per burst
        Data(c).(GroupName).BurstStats.avgSpikesB = mean (Data(c).(GroupName).BurstStats.NspikesB); % average spikes per burst
        Data(c).(GroupName).BurstStats.N = burstEventNidx; % saves the number of burst-events
        Data(c).(GroupName).BurstStats.indices = burst_idx; % vector with all indeces of all ISI from all burst-events
        Data(c).(GroupName).BurstStats.IntraFreq = 1./BurstISI; % vetor with all frequencies from all burst-events

        Data(c).(GroupName).BurstStats.SingBurstEvDur = burstDur; % vector with all burst durations of each and every burst-events
        Data(c).(GroupName).BurstStats.MeanBurstDur = mean(burstDur); % mean of all burst-event-durations
        Data(c).(GroupName).BurstStats.preISI = preBrISI; % vector with all pre-burst-ISIs
        Data(c).(GroupName).BurstStats.indicesPRE = preBrISI_idx; % the indeces of all pre-burst-ISIs
        Data(c).(GroupName).BurstStats.postISI = postBrISI; % vector with all post-burst-ISIs
        Data(c).(GroupName).BurstStats.indicesPOST = postBrISI_idx; % the indeces of all post-burst-ISIs
        Data(c).(GroupName).BurstStats.MEANpreISI = mean(preBrISI); % calculating the mean of all pre-burst-ISIs
        Data(c).(GroupName).BurstStats.MEANpostISI = mean(postBrISI); % calculating the mean of all post-burst-ISIs
        Data(c).(GroupName).BurstStats.maxMaxIBfreq = 1/min(Data(c).(GroupName).BurstISI); % the single one maximal burst freq
        Data(c).(GroupName).BurstStats.maxIBfreq = intra_burst_freq_MAX; % vector with maximal frequencies for all burst-events
        Data(c).(GroupName).BurstStats.MEANmaxIBfreq = mean(intra_burst_freq_MAX); % mean of all max freq of all burst-events
        Data(c).(GroupName).BurstStats.meanIBfreq = intra_burst_freq_MEAN; % vector with the means of each burst-event freq
        Data(c).(GroupName).BurstStats.MEANmeanIBfreq = mean(intra_burst_freq_MEAN); % mean the means of each burst-event freq
        

        Data(c).(GroupName).BurstStats.twoSpikePerBurst = twoSpBr*100/length(burstN);
        Data(c).(GroupName).BurstStats.threeSpikePerBurst = threeSpBr*100/length(burstN);
        Data(c).(GroupName).BurstStats.manySpikePerBurst = mSpBr*100/length(burstN);
   % end
    
    % analysing the non-burst-ISIs data
    singleSpikeISI = []; % vector fro all non-(Pre-/post-)burst-ISIs,
    for p=1:length(ISI) % looking for non-burst, or non-pre-burst, or non-post-burst ISIs
       if ~ismember(p,Data(c).(GroupName).BurstStats.indices) && ~ismember(p,Data(c).(GroupName).BurstStats.indicesPRE) && ~ismember(p,Data(c).(GroupName).BurstStats.indicesPOST)
            singleSpikeISI = [ singleSpikeISI ISI(p)];
        end
    end
    
    Data(c).(GroupName).nonBurstISI = singleSpikeISI;
    Data(c).(GroupName).MeanFreqNONburstISI = mean(1./Data(c).(GroupName).nonBurstISI);             % Mean frequency
    Data(c).(GroupName).MedianFreqNONburstISI = median(1./Data(c).(GroupName).nonBurstISI);         % Median frequency
    
    Data(c).(GroupName).SDMnonBurstISI = std(Data(c).(GroupName).nonBurstISI);                  % Standard deviation of mean
    Data(c).(GroupName).CVnonBurstISI = Data(c).(GroupName).SDMnonBurstISI/mean(Data(c).(GroupName).nonBurstISI);           % Coefficient of variation
    
    singleSpikeISI = [];
    
end


%% Find pauses (Upper outliers that are not preceded by bursts)
for c = 1:nCells
    PauseSeedIdx = Data(c).(GroupName).UpOutIDx;
    
    Preced = PauseSeedIdx-1.; PauseN = 0; PauseIdx = 0; pbPauseN=0; pbPauseIdx=0;
    
    for i=1:length(Preced);
        if ismember(Preced(i), Data(c).(GroupName).BurstStats.indices) == 1     % Defines post-burst pauses
            pbPauseN =  pbPauseN+1;
            pbPauseIdx(pbPauseN) = PauseSeedIdx(i);
        else                                                                    % Defines putative hyperpolarization induced pauses
            PauseN = PauseN+1;
            PauseIdx(PauseN) = PauseSeedIdx(i);
        end
    end
    
    if PauseIdx > 0;
        PauseISI =  Data(c).(GroupName).ISI(PauseIdx');
        PIP = (length(PauseISI)*100)/length(ISI);    % Percentage of ISIs that are pauses
        
        Data(c).(GroupName).PauseISI = PauseISI;
        Data(c).(GroupName).PauseStats.indices = PauseIdx;
        Data(c).(GroupName).PauseStats.MeanPauseDur = mean (PauseISI);
        Data(c).(GroupName).PauseStats.MedianPauseDur = median (PauseISI);
        Data(c).(GroupName).PauseStats.N = PauseN;
        Data(c).(GroupName).PauseStats.PIP = PIP;
    else
        Data(c).(GroupName).PauseISI = 0;
        Data(c).(GroupName).PauseStats.indices = [];
        Data(c).(GroupName).PauseStats.MeanPauseDur = 0;
        Data(c).(GroupName).PauseStats.MedianPauseDur = 0;
        Data(c).(GroupName).PauseStats.N = 0;
        Data(c).(GroupName).PauseStats.PIP = 0;
    end
    
    if pbPauseIdx > 0;
        pbPauseISI =  Data(c).(GroupName).ISI(pbPauseIdx');
        pbPIP = (length(pbPauseISI)*100)/length(ISI);    % Percentage of ISIs that are pauses
        
        Data(c).(GroupName).pbPauseISI = pbPauseISI;
        Data(c).(GroupName).pbPauseStats.indices = pbPauseIdx;
        Data(c).(GroupName).pbPauseStats.MeanPauseDur = mean(pbPauseISI);
        Data(c).(GroupName).pbPauseStats.MedianPauseDur = median(pbPauseISI);
        Data(c).(GroupName).pbPauseStats.N = pbPauseN;
        Data(c).(GroupName).pbPauseStats.PIP = pbPIP;
    else
        Data(c).(GroupName).pbPauseISI = 0;
        Data(c).(GroupName).pbPauseStats.indices = [];
        Data(c).(GroupName).pbPauseStats.MeanPauseDur = 0;
        Data(c).(GroupName).pbPauseStats.MedianPauseDur = 0;
        Data(c).(GroupName).pbPauseStats.N = 0;
        Data(c).(GroupName).pbPauseStats.PIP = 0;
    end
        
    clear Remove PauseISI PauseIdx pbPauseISI pbPauseIdx
end
    
%% Plot histograms with color-coded bursts and pauses

if strcmp(Plot,'yes')
for c = 1:nCells
    ISI = Data(c).(GroupName).ISI;
    Pauses = Data(c).(GroupName).PauseISI;
    Bursts = Data(c).(GroupName).BurstISI; 
    
    xmax = max(ISI);
    [bincounts,binpositions] = hist(ISI,xvalues);
    binwidth = binpositions(2)-binpositions(1);
    histarea = binwidth*sum(bincounts);
    xplot = min(xvalues):binsize:max(xvalues);
    
    f = figure;
    set(f, 'Position', [50, 50, 800, 600]);
    hold on
    
    
    % Control Hist
    subplot('Position', [0.15, 0.1, 0.75, 0.75]);
    hold on
    [MainH,y] = hist(ISI,xvalues); bar(y,MainH,'EdgeColor','k', 'FaceColor','k');
    if Pauses ~= 0;
    [PauseH,y] = hist(Pauses,xvalues); lg2 = bar(y,PauseH,'EdgeColor',[0.2 0.4 1], 'FaceColor',[0.2 0.4 1]);
    else
        lg2 = bar(0,0,'EdgeColor',[0.2 0.4 1], 'FaceColor',[0.2 0.4 1]);
    end
    if Bursts ~= 0;
    [BurstH,y] = hist(Bursts,xvalues); lg1 = bar(y,BurstH,'EdgeColor',[1 0 0.2], 'FaceColor',[1 0 0.2]);
    else
        lg1 = bar(0,0,'EdgeColor',[1 0 0.2], 'FaceColor',[1 0 0.2]);
    end
    line(y,0,'Color', 'k', 'LineWidth', 1);
    
%     ht=title('Control');  % Defines title
%     set(ht,'interpreter','none', 'FontName','calibri','FontWeight','bold','FontSize',16);
    
    xlim([0 2]); % xlim([0 ceil(xmax+xadd)]);
    xlabel(['ISI (', char(scale),')'],'FontSize', 18,  'fontweight','bold','FontName','calibri');
    ylabel('Count (N)','FontSize', 18, 'fontweight','bold', 'FontName','calibri');
    set(gca,'Ticklength', [0 0.01],'FontSize',16,'LineWidth', 1, 'FontName','calibri');
    
    % General / colored ISI histo
    lp = subplot('Position', [0.375, 0.88, 0.3, 0.01]);
    p = get(lp,'position');
    lg = legend([lg1, lg2],lgLabels,'FontSize', 14,'FontName','calibri','Location','northoutside', 'Orientation', 'horizontal', 'box','off');
    set(lg,'position',p); axis(lp,'off');
    % v = get(lg,'title');
    % set(v,'string',texlabel([char(Data(c).(GroupName).CellName) ' - ' char(ID)],'literal'), 'FontSize', 18, 'fontweight','bold', 'FontName','calibri');
    t = title([Data(c).(GroupName).CellName ' - ' Data(c).(GroupName).Period ], 'FontSize', 18, 'fontweight','bold', 'FontName','calibri');
    set(t,'Interpreter','none'); 
    fname = [char(Data(c).(GroupName).CellName) '-' char(Data(c).(GroupName).Period) '_ISI-pauses-bursts'];
    saveas(f,fname,'jpg');
    savefig (fname);
    close all
    
    f = figure;
    set(f, 'Position', [50, 50, 800, 600]);
    hold on
    
    % black-white - ISI histo
    hold on
    [MainH,y] = hist(ISI,xvalues); bar(y,MainH,'EdgeColor','k', 'FaceColor','k');
    
    line(y,0,'Color', 'k', 'LineWidth', 1);
    xlim([0 2]); % xlim([0 ceil(xmax+xadd)]);
    xlabel(['ISI (', char(scale),')'],'FontSize', 18,  'fontweight','bold','FontName','calibri');
    ylabel('Count (N)','FontSize', 18, 'fontweight','bold', 'FontName','calibri');
    set(gca,'Ticklength', [0 0.01],'FontSize',16,'LineWidth', 1, 'FontName','calibri');
    
    t = title([Data(c).(GroupName).CellName ' - ' Data(c).(GroupName).Period ], 'FontSize', 18, 'fontweight','bold', 'FontName','calibri');
    set(t,'Interpreter','none'); 
    fname = [char(Data(c).(GroupName).CellName) '-' char(Data(c).(GroupName).Period)   '_ISI-histoBW'];
    
    saveas(f,fname,'jpg');
    savefig (fname);
    close all
    
    close all
end
else
end

%% Skewness & Kurtosis

for c = 1:nCells;
Data(c).(GroupName).Skewness = skewness(Data(c).(GroupName).ISI, 0);

Data(c).(GroupName).Kurtosis = kurtosis(Data(c).(GroupName).ISI, 0);

end


%% autocorrelation histogram 
if strcmp(Plot,'yes')
    for c = 1:nCells
        spk_t = cumsum(Data(c).(GroupName).ISI); % spike_times: [1 x nSpikes] double of spike times (in s)
        [x1, y1] = ACHplot (spk_t, 0.002, 2);
        [x2, y2] = ACHplot (spk_t, 0.01, 2);
        figure 
        hold on
        plot (y1, x1*5, 'Color', [0 0 0 0.2])
        plot (y2, x2, 'red', 'LineWidth', 1.2)
        set(gca,'FontWeight','bold')
        xlabel ('Delay [s]', 'FontSize', 12);
        ylabel ('Probability', 'FontSize', 12);


        t = title(['ACH: ' char(Data(c).(GroupName).CellName) '-' char(Data(c).(GroupName).Period)] , 'FontSize', 18, 'fontweight','bold', 'FontName','calibri');
        set(t,'Interpreter','none'); 

        %save
        filename = [char(Data(c).(GroupName).CellName) '-' char(Data(c).(GroupName).Period)  '_ACH'];

        savefig (filename);
        saveas(gca, filename,'jpg');
        savefig (filename);
        close all;
    end
end

%% fano factor

    for c = 1:nCells
        [Data(c).(GroupName).FanoFactor, Data(c).(GroupName).FanoFactorSHF] = FanoFactor(Data(c).(GroupName).ISI, Data(c).(GroupName).CellName, Data(c).(GroupName).Period, Plot);
            % , Fcat, FcatCVsq, Fcatshf, FcatshfCVsq, FcatFcatshf

    end


%% Normal Raster
if strcmp(Plot,'yes')
    for c = 1:nCells
        TS = Data(c).(GroupName).TimeStamps; 
        % Set parameters
        timeslice = 30;
        timeStartPoint = floor(min(TS));
        timeEndPoint = floor(max(TS)) + 1;
        
        timetotal = timeEndPoint - timeStartPoint;
        ls = (timetotal/timeslice);


        f = figure(1);
        set(f, 'Position', [50 50 800 ls*30]);
        hold on
        for l = 1:ls;
            t1 = timeslice*(l-1);
            t2 = timeslice*l;

            TSr = TS(TS > timetotal-t2+timeStartPoint & TS <= timetotal-t1+timeStartPoint);

            for i=1:length(TSr);
                RPlot = line([t2+TSr(i)-timetotal-timeStartPoint t2+TSr(i)-timetotal-timeStartPoint],[ls-l-0.42+timeStartPoint/timeslice ls-l+0.42+timeStartPoint/timeslice],'Color',[0 0 0], 'LineWidth',0.017); %
            end
        end
        set(gca,'FontName','arial','FontSize',12, 'YDir','reverse', 'TickDir','out', 'TickLength', [0.00 0]);
        ylim([timeStartPoint/timeslice-0.9 timeStartPoint/timeslice+ls-0.1]);
        
        xlim([-0.9 30.9]);
        xlabel('Time (s)'); % Time is in millisecond
        set(gca,'YTickLabel',[]);
        box on
        title(texlabel([char(Data(c).(GroupName).CellName) ': ' Data(c).(GroupName).Period],'literal')...
            ,'FontSize', 18, 'fontweight','bold', 'FontName','arial');

        fname = [char(Data(c).(GroupName).CellName) '-' char(Data(c).(GroupName).Period)  '-Raster.jpg'];
        saveas(f,fname,'jpg');
        savefig(fname);
        
        close all

        clear TS BurstIDX TSburst TSburstR TSr
    end
end



%% Plot histograms with color-coded bursts and pauses !!! NO x axis limit !!!
if strcmp(Plot,'yes')

    for c = 1:nCells
        ISI = Data(c).(GroupName).ISI;
        Pauses = Data(c).(GroupName).PauseISI;
        Bursts = Data(c).(GroupName).BurstISI; 

        xmax = max(ISI);
        [bincounts,binpositions] = hist(ISI,xvalues);
        binwidth = binpositions(2)-binpositions(1);
        histarea = binwidth*sum(bincounts);
        xplot = min(xvalues):binsize:max(xvalues);

        f = figure;
        set(f, 'Position', [50, 50, 800, 600]);
        hold on

        % Control Hist
        subplot('Position', [0.15, 0.1, 0.75, 0.75]);
        hold on
        [MainH,y] = hist(ISI,xvalues); bar(y,MainH,'EdgeColor','k', 'FaceColor','k');
        if Pauses ~= 0;
            [PauseH,y] = hist(Pauses,xvalues); lg2 = bar(y,PauseH,'EdgeColor','c', 'FaceColor','c');
        else
            lg2 = bar(0,0,'EdgeColor','c', 'FaceColor','c');
        end
        if Bursts ~= 0;
            [BurstH,y] = hist(Bursts,xvalues); lg1 = bar(y,BurstH,'EdgeColor','m', 'FaceColor','m');
        else
            lg1 = bar(0,0,'EdgeColor','m', 'FaceColor','m');
        end
        line(y,0,'Color', 'k', 'LineWidth', 1);

    %     ht=title('Control');  % Defines title
    %     set(ht,'interpreter','none', 'FontName','arial','FontWeight','bold','FontSize',16);

        xlim([0 ceil(xmax+xadd)]);
        xlabel(['ISI (', char(scale),')'],'FontSize', 18, 'FontName','arial');
        ylabel('Count (N)','FontSize', 18, 'FontName','arial');
        set(gca,'Ticklength', [0 0.01],'FontSize',16,'LineWidth', 1, 'FontName','arial');

        % General
        lp = subplot('Position', [0.375, 0.88, 0.3, 0.01]);
        p = get(lp,'position');
        lg = legend([lg1, lg2],lgLabels,'FontSize', 14,'FontName','arial','Location','northoutside', 'Orientation', 'horizontal', 'box','off');
        set(lg,'position',p); axis(lp,'off');

        t = title([Data(c).(GroupName).CellName ' - ' Data(c).(GroupName).Period ], 'FontSize', 18, 'fontweight','bold', 'FontName','calibri');
        set(t,'Interpreter','none'); 
        fname = [char(Data(c).(GroupName).CellName) '-' char(Data(c).(GroupName).Period)  '_ISIhisto-xUnlimited']; 
        saveas(f,fname,'jpg');
        savefig (fname);
        close all
    end
end

%% burst initiation rate
   
for  c = 1:nCells
    RecPeriodDelta = floor ( max (Data(c).(GroupName).TimeStamps) - min(Data(c).(GroupName).TimeStamps)) + 1;
    Data(c).(GroupName).BurstStats.BurstInitFreq = Data(c).(GroupName).BurstStats.N/RecPeriodDelta;
end


%% single burst analysis
if strcmp(Plot,'yes')
    for c = 1:nCells
        if Data(c).(GroupName).BurstStats.SFB > 0 

            tempSEM = std (Data(c).(GroupName).BurstStats.meanIBfreq) / sqrt(length(Data(c).(GroupName).BurstStats.meanIBfreq));
            Data(c).(GroupName).BurstStats.SEMmeanIBfreq = tempSEM;

            tempSEM = std (Data(c).(GroupName).BurstStats.maxIBfreq) / sqrt(length(Data(c).(GroupName).BurstStats.maxIBfreq));
            Data(c).(GroupName).BurstStats.SEMmaxIBfreq = tempSEM;

            % plot
            f = figure;

            edges = [0 0:1:150 150];
            histogram(Data(c).(GroupName).BurstStats.meanIBfreq, edges)
            ylim([0 120]);
            xlabel('Frequency (Hz)','FontSize', 18,  'fontweight','bold','FontName','calibri');
            ylabel('Count (N)','FontSize', 18, 'fontweight','bold', 'FontName','calibri');

            t = title( ['Single burst analysis - MEAN: ' Data(c).(GroupName).CellName ' - ' Data(c).(GroupName).Period ], 'FontSize', 18, 'fontweight','bold', 'FontName','calibri');
            set(t,'Interpreter','none'); 

            set(gca,'Ticklength', [0 0.01],'FontSize',16,'LineWidth', 1, 'FontName','calibri');
            fname = [char(Data(c).(GroupName).CellName) '-' char(Data(c).(GroupName).Period) '_mSingleIBFreq'];
            saveas(f,fname,'jpg');
            savefig (filename);
            close

            % plot
            f = figure;

            edges = [0 0:1:150 150];
            histogram(Data(c).(GroupName).BurstStats.maxIBfreq, edges)
            ylim([0 120]);
            xlabel('Frequency (Hz)','FontSize', 18,  'fontweight','bold','FontName','calibri');
            ylabel('Count (N)','FontSize', 18, 'fontweight','bold', 'FontName','calibri');

            t = title( ['Single burst analysis - MAX: ' Data(c).(GroupName).CellName ' - ' Data(c).(GroupName).Period ], 'FontSize', 18, 'fontweight','bold', 'FontName','calibri');
            set(t,'Interpreter','none'); 

            set(gca,'Ticklength', [0 0.01],'FontSize',16,'LineWidth', 1, 'FontName','calibri');
            fname = [char(Data(c).(GroupName).CellName) '-' char(Data(c).(GroupName).Period)  '_maxSingleIBFreq'];
            saveas(f,fname,'jpg');
            savefig (fname);
            close
        end
    end
end

%% histo all mean single IB freq
if strcmp(Plot,'yes')
    f = figure;
    hold on
    col=0;
    for c = 1:nCells
            if Data(c).(GroupName).BurstStats.SFB > 0

                % plot
                edges = [0 0:1:150 150];
                cc = col/255;
                histogram(Data(c).(GroupName).BurstStats.meanIBfreq, edges, 'FaceColor', [0 0 cc], 'FaceAlpha',0.2)
                ylim([0 150]);
                xlabel('Frequency (Hz)','FontSize', 18,  'fontweight','bold','FontName','calibri');
                ylabel('Count (N)','FontSize', 18, 'fontweight','bold', 'FontName','calibri');

                set(gca,'Ticklength', [0 0.01],'FontSize',16,'LineWidth', 1, 'FontName','calibri');
                col = col + 1;
            end
    end

    titleBL  = 'Mean intra-burst Frequencies for all burst-events';
    title(titleBL, 'FontSize', 18, 'fontweight','bold', 'FontName','calibri');

    fname = 'mSingleIBFreq_allCells';
    saveas(f,fname,'jpg');
    savefig (fname);
    close
    hold off
end
%% histo group mean single IB freq
if strcmp(Plot,'yes')
    temp = [];
    for c = 1:nCells
        if Data(c).(GroupName).BurstStats.SFB > 0
            temp  = [temp Data(c).(GroupName).BurstStats.meanIBfreq];
        end
    end

    f = figure;
    % plot
    edges = [0 0:1:200 200];
    histogram(temp, edges)
    ylim([0 600]);
    xlabel('Frequency (Hz)','FontSize', 18,  'fontweight','bold','FontName','calibri');
    ylabel('Count (N)','FontSize', 18, 'fontweight','bold', 'FontName','calibri');
    titleBL  = 'Group mean IB Freq';
    title(titleBL, 'FontSize', 18, 'fontweight','bold', 'FontName','calibri');
    set(gca,'Ticklength', [0 0.01],'FontSize',16,'LineWidth', 1, 'FontName','calibri');

    fname = 'mSingleIBFreq_group';
    saveas(f,fname,'jpg');
    savefig (fname);
    close
end
%% histo all max single IB freq
if strcmp(Plot,'yes')
    f = figure;
    hold on
    col=0;
    for c = 1:nCells
    if Data(c).(GroupName).BurstStats.SFB > 0
        % plot
        edges = [0 0:1:150 150];
        cc = col/255;
        histogram(Data(c).(GroupName).BurstStats.maxIBfreq, edges, 'FaceColor', [0 0 cc], 'FaceAlpha',0.2)
        ylim([0 150]);
        xlabel('Frequency (Hz)','FontSize', 18,  'fontweight','bold','FontName','calibri');
        ylabel('Count (N)','FontSize', 18, 'fontweight','bold', 'FontName','calibri');
        set(gca,'Ticklength', [0 0.01],'FontSize',16,'LineWidth', 1, 'FontName','calibri');
        col = col + 1;
    end
    end
    titleBL  = 'max IB Freq';
    title(titleBL, 'FontSize', 18, 'fontweight','bold', 'FontName','calibri');
    fname = 'maxSingleIBFreq_allCells';
    saveas(f,fname,'jpg');
    savefig (fname);
    close
    hold off
end
%% histo group max single IB freq
if strcmp(Plot,'yes')
    temp = [];
    for c = 1:nCells
        if Data(c).(GroupName).BurstStats.SFB > 0
            temp  = [temp Data(c).(GroupName).BurstStats.maxIBfreq];
        end
    end

    f = figure;
    % plot
    edges = [0 0:1:150 150];
    histogram(temp, edges);
    xlabel('Frequency (Hz)','FontSize', 18,  'fontweight','bold','FontName','calibri');
    ylabel('Count (N)','FontSize', 18, 'fontweight','bold', 'FontName','calibri');
    titleBL  = 'Group max IB Freq';
    title(titleBL, 'FontSize', 18, 'fontweight','bold', 'FontName','calibri');
    set(gca,'Ticklength', [0 0.01],'FontSize',16,'LineWidth', 1, 'FontName','calibri');

    fname = 'maxSingleIBFreq_group';
    saveas(f,fname,'jpg');
    savefig (fname);
    close
end


%% export data in excel file

tabl = zeros(nCells,15);

for c = 1:nCells
    % tabl (c, 1) = c;
    tabl (c+1, 3) = Data(c).(GroupName).MeanFreq;
    tabl (c+1, 4) = Data(c).(GroupName).MedianFreq;
    tabl (c+1, 5) = Data(c).(GroupName).SDM;  % Standard deviation of mean
    tabl (c+1, 6) = Data(c).(GroupName).CV*100;  % Coefficient of variation
    tabl (c+1, 7) = Data(c).(GroupName).BurstStats.SFB;
    tabl (c+1, 8) = Data(c).(GroupName).Skewness;
    tabl (c+1, 9) = Data(c).(GroupName).Kurtosis;
    
    if Data(c).(GroupName).BurstStats.SFB > 0 
        tabl (c+1, 10) = Data(c).(GroupName).BurstStats.MEANmeanIBfreq;
        tabl (c+1, 11) = Data(c).(GroupName).BurstStats.MEANmaxIBfreq;
        tabl (c+1, 12) = Data(c).(GroupName).BurstStats.maxMaxIBfreq; 
        tabl (c+1, 13) = Data(c).(GroupName).BurstStats.MeanBurstDur;
        tabl (c+1, 14) = Data(c).(GroupName).BurstStats.MEANpreISI;
        tabl (c+1, 15) = Data(c).(GroupName).BurstStats.MEANpostISI;
        tabl (c+1, 16) = Data(c).(GroupName).BurstStats.twoSpikePerBurst;
        tabl (c+1, 17) = Data(c).(GroupName).BurstStats.threeSpikePerBurst;
        tabl (c+1, 18) = Data(c).(GroupName).BurstStats.manySpikePerBurst;
        tabl (c+1, 19) = Data(c).(GroupName).BurstStats.avgSpikesB;
        tabl (c+1, 20) = Data(c).(GroupName).BurstStats.BurstInitFreq;
    end
    
    tabl (c+1, 21) = Data(c).(GroupName).FanoFactor;
    tabl (c+1, 22) = Data(c).(GroupName).FanoFactorSHF;
    tabl (c+1, 23) = Data(c).(GroupName).PauseStats.MeanPauseDur;
    tabl (c+1, 24) = Data(c).(GroupName).PauseStats.MedianPauseDur;
    tabl (c+1, 25) = Data(c).(GroupName).PauseStats.N;
    tabl (c+1, 26) = Data(c).(GroupName).PauseStats.PIP;
    
    tabl (c+1, 27) = Data(c).(GroupName).MeanFreqNONburstISI;  % Mean frequency without bursts isi
    tabl (c+1, 28) = Data(c).(GroupName).MedianFreqNONburstISI;         % Median frequency without bursts isi
    tabl (c+1, 29) = Data(c).(GroupName).CVnonBurstISI;
    
end
    
filename = [GroupName '.xlsx'];
xlswrite (filename, tabl);

tabl2 = zeros(5000,nCells);
for c = 1:nCells
    if Data(c).(GroupName).BurstStats.SFB > 0
        tabl2(1:length(Data(c).(GroupName).BurstStats.meanIBfreq),c) = Data(c).(GroupName).BurstStats.meanIBfreq;
    end
end

filename = [GroupName '_mSingleIBF.xlsx'];
xlswrite (filename, tabl2);

tabl3 = zeros(5000,nCells);
for c = 1:nCells
    if Data(c).(GroupName).BurstStats.SFB > 0
        tabl3(1:length(Data(c).(GroupName).BurstStats.maxIBfreq),c) = Data(c).(GroupName).BurstStats.maxIBfreq;
    end
end

filename = [GroupName '_maxSingleIBF.xlsx'];
xlswrite (filename, tabl3);

tabl4CellNames = repmat (blanks(5),nCells,1);
for c = 1:nCells
    NameLength = length (char(Data(c).(GroupName).CellName));
    tabl4CellNames (c, 1:NameLength) =char(Data(c).(GroupName).CellName);
end


fileID = fopen('CellName.txt','wt');
for r=1:size(tabl4CellNames,1)
    fprintf(fileID,'%s\n',tabl4CellNames(r,:));
end
fclose(fileID);


% tabl5CellPeriod = repmat (blanks(10),nCells,1);
% for c = 1:nCells
%     PeriodLength = length (Data(c).(GroupName).Period);
%     tabl5CellPeriod (c, 1:PeriodLength) = Data(c).(GroupName).Period;
% end
% 
% fileID = fopen('CellPeriods.txt','wt');
% for r=1:size(tabl5CellPeriod,1)
%     fprintf(fileID,'%s\n',tabl5CellPeriod(r,:));
% end
% fclose(fileID);

%% End function

close all
clearvars -except Data DrugData nCells GroupName GroupName DrugGroupName oldFolder OutlierCriterion;
save (['ISI_stats_Single_' char(GroupName)]);
cd(oldFolder)
clearvars


%% functions

%% ACH 
function [ac_plot,xbin_plot] = ACHplot(spk_t,binsize,max_t)
 
    xbin_centers = - max_t - binsize:binsize:max_t + binsize; % first and last bins are to be deleted later; symmetric - from minus to plus
    ac = zeros(size(xbin_centers)); % results

    for iSpk = 1:length(spk_t)

       relative_spk_t = spk_t - spk_t(iSpk); % calculate the differences between the spike train with each of it components

       ac = ac + hist(relative_spk_t,xbin_centers); % note that hist() puts all spikes outside the bin centers in the first and last bins! delete later.

    end

    xbin = xbin_centers(2:end-1); % remove unwanted bins
    ac = ac(2:end-1);

    ac = ac./max(ac)*100; % normalize

    % to show only the positive values (right side)
    xbin_plot = xbin(xbin > 0);
    ac_plot = ac (length(xbin_plot)+2:end);

end

%% Fano factor

function  [F , Fshf] = FanoFactor (ISI, colheaders, period, pplot)
   
    RecTime = 600;                              % Maximum recording time 
    data = cumsum(ISI);
    Tinc = 0.1;                                 % Interval increment
    T = [Tinc:Tinc:RecTime/2] ;                 % Intervals for counting number of spikes (s)
    F = zeros(1, length(T));                       % Opens a matrix for the compilation of all FF functions
    
    % CV and CV^2 calculation

    display('Computing CV...')
    ISIsdm = nanstd(ISI);           % Calculates standard deviation of ISIs
    ISImean = nanmean(ISI);         % Calculates mean of ISIs
    CV = ISIsdm./ISImean;           % Calculates CV
    CVsq = CV.^2;                   % Calculates CV squared
    CV = CV'; CVsq = CVsq'; 
    
    
    % Fano factor calculation

    %%% Time delimited time count
    
         
    
   for T = Tinc:Tinc:(RecTime/2);     % Intervals for counting number of spikes (s)
        index = round(T/Tinc);              % Defines index for stacking data in Fcat        
        interval = [0:T:600];               % Sets interval limits
        Ni = length(interval)-1;            % Number of intervals
        countstack = zeros(Ni,1);           % Creates matrix that will receive spike counts

        for i = 1:Ni
            datawindow = (data >= interval(i) & data <= interval(i)+T); % Logical operation for defining spike counts
            count = sum(datawindow==1);                                 % Counts the number of spikes in interval 
            countstack(i,:) = count(1,:);                               % Sums all counts in one matrix
            clear datawindow;
        end

        clear Ni; clear interval; 

        %%% Calculating FF 

        F = var(countstack(:,1))/mean(countstack(:,1));    % Calculates the Fano Factor (F) for each neuron
         
             Fcat(:,index) = F;                                    % Compiles all individual FF functions in one matrix
% 
        if      index == 1;
                display('Computing Fano Factor...');

        elseif  index == ((RecTime/2)/100)*25;
                display('25% complete');

        elseif  index == ((RecTime/2)/100)*50;
                display('50% complete');

        elseif  index == ((RecTime/2)/100)*75;
                display('75% complete');

        elseif  index == RecTime/2;
                display('Plotting FF(T)functions...');

        end
    end

    %%% Ploting F(T)

    T = [Tinc:Tinc:RecTime/2] ;                 % Defines T variable
    % mkdir(['FF.T_' char(GroupName)]);           % Creats directory for saving figures
    % oldFolder = cd(['FF.T_' char(GroupName)]);  % Sets the return to old directory
   
    figH=zeros(3); % Defines plot handles 

    if strcmp(pplot,'yes')                   
        x = figure;
        figH = plot(T, Fcat, 'b');                                % Plots the desired FF function for a given figure handle
        set(figH,'LineWidth',1.5);                                     % Defines line width
        set(gca,'FontName','arial','FontWeight','bold','FontSize',16);      % Defines font properties
        set(gca,'XScale','log', 'YScale', 'log');                           % Defines log scale
        axis([0.1 50 0.01 60])                                              % Defines axes ranges
        xlabel('T (s)');ylabel('FF (N spikes)');                            % Defines axes labels
        ht=title(colheaders);                                          % Defines title
        set(ht,'interpreter','none')                                        % Allows underscores in the title
        fname=['FF.T_' char(colheaders) '-' period '.jpg'];                        % Defines figure title (neuron code from colheaders)
        saveas(x,fname,'jpg');
        savefig (fname);               % Saves figure
                                                     
        close all
    end

%%% Plotting FF/CV^2(T)

    display('Plotting FF/CV^2(T)functions...')
   
    FcatCVsq(1,:) = Fcat(1,:)./CVsq(1);
    if strcmp(pplot,'yes')
        x = figure;
        figH= plot(T, FcatCVsq(1,:), 'b');                                % Plots the desired FF function for a given figure handle
        set(figH(:,1),'LineWidth',1.5);                                     % Defines line width
        set(gca,'FontName','arial','FontWeight','bold','FontSize',16);      % Defines font properties
        set(gca,'XScale','log', 'YScale', 'log');                           % Defines log scale
        axis([0.1 50 0.01 60])                                              % Defines axes ranges
        xlabel('T (s)');ylabel('FF/CV^2');                                  % Defines axes labels
        ht=title(colheaders);                                          % Defines title
        set(ht,'interpreter','none')                                        % Allows underscores in the title
        fname=['FF_CV2_' char(colheaders(1,:)) '-' period '.jpg'];                               % Defines figure title (neuron code from colheaders)
        saveas(x,fname,'jpg');
        savefig (fname);                                             % Saves figure
        close all
    end
    
    
%%% Shuffled Fano Factor calculation

    display('Shuffling ISIs...')            % ISI matrix to be shuffled
    ix = randperm(length(ISI));             % randomly indexing the ISI array
    ISIshf = ISI(ix);                       % applying the random indexing of the ISI 
    datashf = cumsum(ISIshf);        % datashf = nancumsum(ISIshf,1,2); 
    
    datashf(:,1) = sortrows(datashf(:,1));  % Orders the spike time matrix


    for T = [Tinc:Tinc:RecTime/2] ;     % Intervals for counting number of spikes (s)
        index = round(T/Tinc);              % Defines index for stacking data in Fcat        
        interval = [0:T:600];               % Sets interval limits
        Ni = length(interval)-1;            % Number of intervals
        countstack = zeros(Ni,1);           % Creates matrix that will receive spike counts

        for i = 1:Ni
            datawindow = (datashf >= interval(i) & datashf <= interval(i)+T);   % Logical operation for defining spike counts
            count = sum(datawindow==1);                                         % Counts the number of spikes in interval 
            countstack(i,:) = count(1,:);                                       % Sums all counts in one matrix
            clear datawindow;
        end

        clear Ni; clear interval; 

    %%% Calculating shuffled FF 

    
        Fshf(1) = var(countstack(:,1))/mean(countstack(:,1));    % Calculates the shuffled Fano Factor (F) for each neuron
     
        Fcatshf(:,index) = Fshf;                                 % Compiles all individual shuffled FF functions in one matrix

        if      index == 1;
                display('Computing shuffled Fano Factor...');

        elseif  index == ((RecTime/2)/100)*25;
                display('25% complete');

        elseif  index == ((RecTime/2)/100)*50;
                display('50% complete');

        elseif  index == ((RecTime/2)/100)*75;
                display('75% complete');

        elseif  index == RecTime/2;
                display('Plotting FFshf/CV^2(T)functions...');       
        end
    end
   
   
    %%% Ploting FFshf(T)

    T = [Tinc:Tinc:RecTime/2] ;                     % Defines T variable
%     mkdir(['FFshf.T_' char(GroupName)]);            % Creats directory for saving figures
%     oldFolder = cd(['FFshf.T_' char(GroupName)]);   % Sets the return to old directory
    figH=zeros(3,1);                                % Defines plot handles 

    if strcmp(pplot,'yes')                    
        x = figure;
        figH (:,1) = plot(T, Fcatshf(1,:), 'b');                                % Plots the desired FF function for a given figure handle
        set(figH(:,1),'LineWidth',1.5);                                     % Defines line width
        set(gca,'FontName','arial','FontWeight','bold','FontSize',16);      % Defines font properties
        set(gca,'XScale','log', 'YScale', 'log');                           % Defines log scale
        axis([0.1 50 0.01 60])                                              % Defines axes ranges
        xlabel('T (s)');ylabel('FF_{shuffled} (N spikes)');                   % Defines axes labels
        ht=title(colheaders);                                          % Defines title
        set(ht,'interpreter','none')                                        % Allows underscores in the title
        fname=['FFshf.T_' char(colheaders) '-' period '.jpg'];                    % Defines figure title (neuron code from colheaders)
        saveas(x,fname,'jpg');
        % saveas(x,fname,'fig');                                              % Saves figure
        close all
    end
   

    %%% Plotting FFshf/CV^2(T)

    
      
        FcatshfCVsq(1,:) = Fcatshf(1,:)./CVsq(1);
        if strcmp(pplot,'yes')
            x = figure;
            figH (:,1) = plot(T, FcatshfCVsq(1,:), 'b');                                % Plots the desired FF function for a given figure handle
            set(figH(:,1),'LineWidth',1.5);                                     % Defines line width
            set(gca,'FontName','arial','FontWeight','bold','FontSize',16);      % Defines font properties
            set(gca,'XScale','log', 'YScale', 'log');                           % Defines log scale
            axis([0.1 50 0.01 60])                                              % Defines axes ranges
            xlabel('T (s)');ylabel('FF_{shuffled}/CV^2');                               % Defines axes labels
            ht=title(colheaders);                                          % Defines title
            set(ht,'interpreter','none')                                        % Allows underscores in the title
            fname=['FFshf.CV.T_' char(colheaders) '-' period '.jpg'];                 % Defines figure title (neuron code from colheaders)
            saveas(x,fname,'jpg');
            % saveas(x,fname,'fig');                                              % Saves figure
            close all
        end
    
    

    %%% Plotting FF/FFshf(T)

    display('Plotting FF/FFshf(T)functions...');

%     mkdir(['FF.FFshf.T_' char(GroupName)]);             % Creats directory for saving figures
%     oldFolder = cd(['FF.FFshf.T_' char(GroupName)]);    % Sets the return to old directory
    
    FcatFcatshf = Fcat./Fcatshf;                        % Elementwise division of FF and shuffled FF matrices

    if strcmp(pplot,'yes')
        x = figure;
        figH (:,1)= plot(T, FcatFcatshf(1,:), 'b');                                % Plots the desired FF function for a given figure handle
        set(figH(:,1),'LineWidth',1.5);                                     % Defines line width
        set(gca,'FontName','arial','FontWeight','bold','FontSize',16);      % Defines font properties
        set(gca,'XScale','log', 'YScale', 'log');                           % Defines log scale
        axis([0.1 50 0.01 60])                                              % Defines axes ranges
        xlabel('T (s)');ylabel('FF/FF_{shuffled}');                           % Defines axes labels
        ht=title(colheaders);                                          % Defines title
        set(ht,'interpreter','none')                                        % Allows underscores in the title
        fname=['FF.FFshf.T_' char(colheaders) '-' period '.jpg'];                 % Defines figure title (neuron code from colheaders)
        saveas(x,fname,'jpg');
        % saveas(x,fname,'fig');                                              % Saves figure
        close all
    end
    
end




