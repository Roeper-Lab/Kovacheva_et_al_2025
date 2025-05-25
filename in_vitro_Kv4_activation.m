clear all
close all

%% load data



files = dir ('in_vitro_Kv4_activation_example.mat');
GroupName = 'Kv4_activation';  % aktueller Name der Gruppe
[nCells,~] =size(files);

for i=1:nCells
    temp = load(files(i).name);
    name = files(i).name;
    tempTracNames = fieldnames(temp);
    
    % checks if the traces are exported double as _1 and _2 + removes _2
    for j=1:length(tempTracNames)
        a = char(tempTracNames(j));
        a(end);
        if a(end) == '2'            
            temp = rmfield(temp, a);
        end
    end
       
    ds(i).cell = temp;
    ds(i).cell.nTr = length(fieldnames(ds(i).cell)); % number of traces
    ds(i).cell.name = name(1:end-4);
    tt = fieldnames(ds(i).cell);
    tt2 = cell2mat(tt(1));
    tfind = strfind(tt2, '_');
    ds(i).cell.TrName =  tt2(1:tfind(3));  % initial part of the trace name, example 'Trace_1_2_'
    % ds(i).cell.stim = ds(i).cell.trace (:,3);
    
    
end

clear files i temp name tt tt2


% Set folders
mkdir(['Stat_' char(GroupName)]);                   % creates a new directory with for the analysis
oldFolder = cd(['Stat_' char(GroupName)]);          % Sets the return to old directory
addpath(pwd);  

%% plot raw data
p=1;
for k=1:nCells
    % 
    n_tr = ds(k).cell.nTr; % number of traces

    f = figure('Position',[100 50 700 1000],'Renderer', 'painters');
    
    hold on
    title(ds(k).cell.name);
    for i=1:n_tr-1
        subplot (n_tr,2,p);
        
        TS = ds(k).cell.(char([ds(k).cell.TrName int2str(i) '_1']))(:,1); % time stamps of the current trace
        V = ds(k).cell.(char([ds(k).cell.TrName int2str(i) '_1']))(:,2); % voltage of the current trace
        plot(TS, V);
        ylim([-0.8e-9 15.0e-9]);
        title([ds(k).cell.name ' - ' int2str(i)]);
        
        p=p+1;
        
        subplot (n_tr,2,p);
        TS = ds(k).cell.(char([ds(k).cell.TrName int2str(i) '_1']))(:,1); % time stamps of the current trace
        V = ds(k).cell.(char([ds(k).cell.TrName int2str(i) '_1']))(:,3); % voltage of the current trace
        plot(TS, V);
        ylim([-0.12 0.0]);
        
        p=p+1;
    end
%     fname = [GroupName '_' ds(k).cell.name '_raw1'];
%     saveas (f, fname, 'emf')
%     saveas (f, fname, 'jpg')
    
    p=1;
end

%%
for k=1:nCells
    
    n_tr = ds(k).cell.nTr; % number of traces

    f = figure('Position',[100 200 1400 800], 'Renderer', 'painters');
    hold on
    for i=1:n_tr
        subplot (2,1,1)
        hold on
            plot(ds(k).cell.(char([ds(k).cell.TrName int2str(i) '_1']))(:,1), ds(k).cell.(char([ds(k).cell.TrName int2str(i) '_1']))(:,2), 'k')
            xlim([1 6]);
            ylim([-5.0e-9 15.0e-9]);
            title(ds(k).cell.name);
        subplot (2,1,2)
        hold on
            plot(ds(k).cell.(char([ds(k).cell.TrName int2str(i) '_1']))(:,1), ds(k).cell.(char([ds(k).cell.TrName int2str(i) '_1']))(:,3), 'r')
            xlim([1 6]);
    end

    hold off
    
%     fname = [GroupName '_' ds(k).cell.name '_raw2'];
%     saveas (f, fname, 'emf')
%     saveas (f, fname, 'jpg')
end


%% find peaks

ProtV = [-60:5:-15]; % pulse protocol voltage: -15 mV, -20mV ....


for k=1:nCells
    % 
    n_tr = ds(k).cell.nTr; % number of traces

    f = figure('Position',[100 200 350 200],'Renderer', 'painters');
    hold on
    for i=1:n_tr
        TS = ds(k).cell.(char([ds(k).cell.TrName int2str(i) '_1']))(:,1); % time stamps of the current trace
        V = ds(k).cell.(char([ds(k).cell.TrName int2str(i) '_1']))(:,2); % voltage of the current trace
        Fs = length(TS)./TS(end); % sample rate
        
        idx = 3.002*Fs; % index of 3rd sec 
        idx2 = 3.5*Fs; % index of 6th sec 
        [pks,locs] = max( V(idx:idx2) );
        ds(k).cell.(char(['PeakV_' int2str(i)])) = pks;
        ds(k).cell.(char(['PeakTS_' int2str(i)])) = locs+round(idx);
        ds(k).cell.(char(['PeakDelay_ms_' int2str(i)])) =  (TS(locs+round(idx))-3)*1000;
        
        if ds(k).cell.name(end-3:end) == 'op11' | ds(k).cell.name(end-3:end) == 'op12' | ds(k).cell.name(end-3:end) == 'op13' | ds(k).cell.name(end-3:end) == 'op14' | ds(k).cell.name(end-3:end) == 'op15'
            idx3 = 3.898*Fs; % index for begin of last 100ms
            idx4 = 3.998*Fs; % index for end of 100ms
        else
            idx3 = 4.898*Fs; % index for begin of last 100ms
            idx4 = 4.998*Fs; % index for end of 100ms
        end


        
        I_steady = mean(V(idx3:idx4));
        ds(k).cell.(char(['Isteady_' int2str(i)])) = I_steady;
                
        delta = pks-I_steady;
        ds(k).cell.(char(['delta_' int2str(i)])) = delta;
        
        
        
        
        Itau_lim = pks-abs(delta * 1/exp(1)); % where current is 1/e  // before 2/exp(1)
         
        t = locs+round(idx);
        while V(t)>Itau_lim
            t=t+1;
        end
        
        ds(k).cell.(char(['Itau_A' int2str(i)])) = V(t);
        ds(k).cell.(char(['ItauTS' int2str(i)])) = t;
        ds(k).cell.(char(['tau_s' int2str(i)])) = TS(round(t))-TS(locs+round(idx));
        ds(k).cell.(char(['conduct_' int2str(i)])) = delta/(100+ProtV(i));
        
        
        
        
        % plot raw data with peaks
        plot(TS, V, 'k')
        scatter(TS(locs+round(idx)), pks, 'o');
        
        scatter(TS(t), V(t), 'o');
        
        line([TS(round(idx3)) TS(round(idx4))], [I_steady I_steady], 'Color' , 'r');
        xlim([2.99 5.04]);
        ylim([-0.8e-9 10.0e-9]); %  15.0e-9]);
        
    end
    hold off
    fname = [GroupName '_' ds(k).cell.name '_pks'];
    title(ds(k).cell.name)
    saveas (f, fname, 'emf')
    saveas (f, fname, 'jpg')
    
    
end

f = figure('Renderer', 'painters');
hold on
t = 1/nCells;
faceCol = [0 1 .8];
for  k=1:nCells
    n_tr = ds(k).cell.nTr; % number of traces
    faceCol(1) = faceCol(1)+t;
    faceCol(2) = faceCol(2)-t;
    for i=1:n_tr
        scatter(ProtV(i), ds(k).cell.(char(['PeakV_' int2str(i)])),  'MarkerFaceColor','none','MarkerEdgeAlpha',.5, 'MarkerEdgeColor',abs(faceCol))
    end
end

hold off
fname = [GroupName '_all_pksAmplitude'];
title(fname)
saveas (f, fname, 'emf')
saveas (f, fname, 'jpg')

%% export data


%export current peak [pA]
tabl = zeros(nCells,10);

for k = 1:nCells
    n_tr = ds(k).cell.nTr; % number of traces
    for i=1:n_tr
        tabl (k+1, i+1) = ds(k).cell.(char(['PeakV_' int2str(i)]))*10^12;
    end
    
end

filename = [GroupName '_peaks_groupAnalysis.xlsx'];
xlswrite (filename, tabl);

% export current steady [pA]

tabl = zeros(nCells,10);

for k = 1:nCells
    n_tr = ds(k).cell.nTr; % number of traces
    for i=1:n_tr
        tabl (k+1, i+1) = ds(k).cell.(char(['Isteady_' int2str(i)]))*10^12;
    end
    
end

filename = [GroupName '_Isteady_groupAnalysis.xlsx'];
xlswrite (filename, tabl);

% export delta amplitude [pA]
tabl = zeros(nCells,10);

for k = 1:nCells
    n_tr = ds(k).cell.nTr; % number of traces
    for i=1:n_tr
        tabl (k+1, i+1) = ds(k).cell.(char(['delta_' int2str(i)]))*10^12;
    end
    
end

filename = [GroupName '_delta_groupAnalysis.xlsx'];
xlswrite (filename, tabl);





% export conductance[S]
tabl = zeros(nCells,10);

for k = 1:nCells
    n_tr = ds(k).cell.nTr; % number of traces
    for i=1:n_tr
        tabl (k+1, i+1) = ds(k).cell.(char(['conduct_' int2str(i)]))*10^12;
    end
    
end

filename = [GroupName '_conduct_groupAnalysis.xlsx'];
xlswrite (filename, tabl);






% export tau in ms

tabl = zeros(nCells,10);

for k = 1:nCells
    n_tr = ds(k).cell.nTr; % number of traces
    for i=1:n_tr
        tabl (k+1, i+1) =ds(k).cell.(char(['tau_s' int2str(i)]))*1000;
    end
    
end

filename = [GroupName '_tau_groupAnalysis.xlsx'];
xlswrite (filename, tabl);

%ds(k).cell.(char(['tau_s' int2str(i)]));

% names expport

tabl4CellNames = repmat (blanks(12),nCells,1);
for i = 1:nCells
    NameLength = length (ds(i).cell.name);
    tabl4CellNames (i, 1:NameLength) = ds(i).cell.name;
end

fileID = fopen('CellName.txt','wt');
for r=1:size(tabl4CellNames,1)
    fprintf(fileID,'%s\n',tabl4CellNames(r,:));
end
fclose(fileID);

save ([GroupName '_groupAnalysis.mat'])
close all