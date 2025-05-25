% v3 - different symbol size in dependence to rot degree
close all
clear all

load (['behavior_rotSequence.mat'])

%%
for k = 1:nSessions
    session = char(SessionNames(k));
    f = figure ('Position', [0 0 1400 1000]);
    
    for i=1:nAnimals
        subplot(5, 5, i);
        hold on 
            xlim([0 600]);
            line ([0 600], [0 0],'Color', [.7 .7 .7])
            % ylim([-2 2]);
            r_counter = 1; l_counter = - 1; 
            rot_PlotValue = [];
            for j = 1:RotData(i).animal.(session).nRot
                curr_degree = RotData(i).animal.(session).rotDegree(j);
                if RotData(i).animal.(session).rotDir(j) == 'Left'
                    scatter(RotData(i).animal.(session).timeStamps(j), l_counter, (2+30*(curr_degree/500)),...
                        'MarkerFaceColor', [1 .6 .4], 'MarkerEdgeColor', 'none');
                    line ([RotData(i).animal.(session).timeStamps(j) RotData(i).animal.(session).timeStamps(j)], [0 l_counter],...
                        'Color', [.8 .8 .8],'LineStyle',':');
                    rot_PlotValue = [rot_PlotValue l_counter];
                    r_counter = 1;
                    l_counter = l_counter - 1;
                else 
                    scatter(RotData(i).animal.(session).timeStamps(j), r_counter,(2+30*(curr_degree/500)), ...
                        'MarkerFaceColor', [0 .8 .6], 'MarkerEdgeColor', 'none');
                    line ([RotData(i).animal.(session).timeStamps(j) RotData(i).animal.(session).timeStamps(j)], [0 r_counter],...
                        'Color', [.8 .8 .8],'LineStyle',':')
                    rot_PlotValue = [rot_PlotValue r_counter];
                    r_counter = r_counter + 1;
                    l_counter = - 1;
                end
            end
            plot(RotData(i).animal.(session).timeStamps, rot_PlotValue, 'Color', [.8 .8 .8]);

            area(RotData(i).animal.(session).timeStamps, rot_PlotValue, 'LineStyle', 'none', 'FaceColor', 'k', 'FaceAlpha', .07);
            
            title([RotData(i).animal.Aname ' - ' session], 'interpreter', 'none');
        hold off
        
        RotData(i).animal.(session).rot_yPlotValue = rot_PlotValue;
        RotData(i).animal.(session).AUC = trapz(RotData(i).animal.(session).timeStamps, rot_PlotValue);

        transitions = 0;
        for j = 1:RotData(i).animal.(session).nRot-1
            if (RotData(i).animal.(session).rot_yPlotValue(j)>0 & RotData(i).animal.(session).rot_yPlotValue(j+1)<0) || ...
                    (RotData(i).animal.(session).rot_yPlotValue(j)<0 & RotData(i).animal.(session).rot_yPlotValue(j+1)>0)
                transitions = transitions +1;
            end
        end
        
        RotData(i).animal.(session).transitions = transitions;
    end
    
    subplot(5,5, i+1)
    hold on
        xlim([0 600]);
        line ([0 600], [0 0],'Color', [.7 .7 .7])
        ylim([-6 6]);
        scatter (500, -5.8, (2+30*(450/500)), 'MarkerFaceColor', [0 .8 .6],  'MarkerEdgeColor', 'none') % example of 450 ? , 'MarkerEdgeColor', 'k'
        text (510, -5.8, '450')
        scatter (500, -5, (2+30*(360/500)), 'MarkerFaceColor',[0 .8 .6],  'MarkerEdgeColor', 'none') % example of 360 ?
        text (510, -5, '360')
        scatter (500, -4.5, (2+30*(270/500)), 'MarkerFaceColor', [0 .8 .6],  'MarkerEdgeColor', 'none') % example of 270 ?
        text (510, -4.5, '270')
        scatter (500, -4, (2+30*(180/500)), 'MarkerFaceColor', [0 .8 .6],  'MarkerEdgeColor', 'none') % example of 180?
        text (510, -4, '180')
        scatter (500, -3.5,(2+30*(90/500)), 'MarkerFaceColor',[0 .8 .6],  'MarkerEdgeColor', 'none') % exampl of 90?
        text (510, -3.5, '90')
    hold off 
    fname = [session '_rotSequence'];
    % saveas (f, fname, 'jpg'); % to save as jpg
    if strcmp(session ,'BL2') ||  strcmp(session ,'BL3') ||  strcmp(session ,'POP4') ||  strcmp(session ,'POP68') || ...
            strcmp(session ,'POP64') || strcmp(session ,'POP28')|| strcmp(session ,'POP20') 
        % saveas (f, fname, 'fig'); % to save as matlab figure
    end
    % close(f) 
    
end

% save (['_rotSequence_analysis.mat'])


%% export data 
% uncomment to export data as excel files
% 
% tabl = zeros(nSessions+1, nAnimals+1);
% tabl2 = zeros(nSessions+1, nAnimals+1);
% for i = 1:nAnimals 
%     for j = 1:nSessions
%         session = char(SessionNames(j));
%         tabl (j+1, i+1) = RotData(i).animal.(session).transitions;
%         tabl2 (j+1, i+1) = RotData(i).animal.(session).AUC;
%     end
% end
% 
% filename = ['_rotSeq_transitions.xlsx'];
% xlswrite (filename, tabl);
% filename = ['_rotSeq_AUC.xlsx'];
% xlswrite (filename, tabl2);
% 
% % Animal names export
% tabl4animalNames = repmat (blanks(12),nAnimals,1);
% row = 1;
% for i = 1:nAnimals
%     NameLength = length (RotData(i).animal.Aname);
%     tabl4animalNames (i, 1:NameLength) = RotData(i).animal.Aname; % 
% end
%      
% fileID = fopen('_Animal_names.txt','wt');
% for r=1:size(tabl4animalNames,1)
%      fprintf(fileID,'%s\n',tabl4animalNames(r,:));
% end
% fclose(fileID);
% 
% % session names export
% tabl4sessionNames = repmat (blanks(12),nSessions,1);
% row = 1;
% for i = 1:nSessions
%     NameLength = length (char(SessionNames(i)));
%     tabl4sessionNames (i, 1:NameLength) = char(SessionNames(i)); % 
% end
%      
% fileID = fopen('_Session_names.txt','wt');
% for r=1:size(tabl4sessionNames,1)
%      fprintf(fileID,'%s\n',tabl4sessionNames(r,:));
% end
% fclose(fileID);