%% Analysis of soundwalk
clear all 
clc

%load all data for questions 1-4
load('sw_2022_group_1A');
load('sw_2022_group_1B');
load('sw_2022_group_2A');
load('sw_2022_group_2B');

%redefine names for simplicity
groups.group_1A = cat(3, sw_2022_group_1A(:).Questions_1234);
groups.group_1B = cat(3, sw_2022_group_1B(:).Questions_1234);
groups.group_2A = cat(3, sw_2022_group_2A(:).Questions_1234);
groups.group_2B = cat(3, sw_2022_group_2B(:).Questions_1234);
groups.group_1A_2A = cat(3, groups.group_1A, groups.group_2A);
groups.group_1B_2B = cat(3, groups.group_1B, groups.group_2B);
groups.group_1A_1B_2A_2B = cat(3, groups.group_1A, groups.group_1B, groups.group_2A, groups.group_2B);

%% PROBLEM 1

%Calculate aritmetic mean for all group combinations (see function below)
part_1_arithmetic_mean_Q_1 = a_mean_all_groups(groups,1).*4./12.4;
save('part_1_arithmetic_mean_Q_1.mat' , 'part_1_arithmetic_mean_Q_1');

part_1_arithmetic_mean_Q_2 = a_mean_all_groups(groups,2).*4./12.4;
save('part_1_arithmetic_mean_Q_2.mat' , 'part_1_arithmetic_mean_Q_2');

part_1_arithmetic_mean_Q_3 = a_mean_all_groups(groups,3).*4./12.4;
save('part_1_arithmetic_mean_Q_3.mat' , 'part_1_arithmetic_mean_Q_3');


%Calculate 95% confidense intervals for the first three qs for all groups
%included
part_1_confidence_intervals_Q_1 = conf_int(groups.group_1A_1B_2A_2B, 1, 0.05);
save('part_1_confidence_intervals_Q_1.mat', 'part_1_confidence_intervals_Q_1');

part_1_confidence_intervals_Q_2 = conf_int(groups.group_1A_1B_2A_2B, 2, 0.05);
save('part_1_confidence_intervals_Q_2.mat', 'part_1_confidence_intervals_Q_2');

part_1_confidence_intervals_Q_3 = conf_int(groups.group_1A_1B_2A_2B, 3, 0.05);
save('part_1_confidence_intervals_Q_3.mat', 'part_1_confidence_intervals_Q_3');

%% PROBLEM 2
%Show the results in a clear graphical representation (arithmetic means for all groups and all 
%combinations and the 95% confidence intervals only for the combination of groups 1A & 1B & 2A & 2B).
%Hint: look at the Matlab errorbar function.

groupnames = ["Group 1A", "Group 1B", "Group 2A", "Group 2B", "Group 1A & 2A", "Group 1B & 2B", "All groups"];

figure(1);
sgtitle('How loud is it here? [0-4](Perceived Loudness)');

for i = 1:7;
    subplot(4,2,i);
    bar([1:1:8],part_1_arithmetic_mean_Q_1(i,:))
    ylim([0 4]);
    subtitle(groupnames(i));
    xlabel('Location');
    grid on
    if groupnames(i) == "All groups";
        hold on
        nerr = abs(part_1_arithmetic_mean_Q_1(i,:) - part_1_confidence_intervals_Q_1(1,:));
        perr = abs(part_1_arithmetic_mean_Q_1(i,:) - part_1_confidence_intervals_Q_1(2,:));
        er = errorbar([1:1:8],part_1_arithmetic_mean_Q_1(i,:), nerr, perr);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none'; 
        hold off
    end
end

figure(2);
sgtitle('How appropriate is the sound to the surrounding? [0-4]');

for i = 1:7;
    subplot(4,2,i);
    bar([1:1:8],part_1_arithmetic_mean_Q_2(i,:))
    ylim([0 4]);
    subtitle(groupnames(i));
    xlabel('Location');
    grid on
    if groupnames(i) == "All groups";
        hold on
        nerr = abs(part_1_arithmetic_mean_Q_2(i,:) - part_1_confidence_intervals_Q_2(1,:));
        perr = abs(part_1_arithmetic_mean_Q_2(i,:) - part_1_confidence_intervals_Q_2(2,:));
        er = errorbar([1:1:8],part_1_arithmetic_mean_Q_2(i,:), nerr, perr);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none'; 
        hold off
    end
end

figure(3);
sgtitle('How often would you like to visit this place again? [0-4]');

for i = 1:7;
    subplot(4,2,i);
    bar([1:1:8],part_1_arithmetic_mean_Q_3(i,:))
    ylim([0 4]);
    subtitle(groupnames(i));
    xlabel('Location');
    grid on
    if groupnames(i) == "All groups";
        hold on
        nerr = abs(part_1_arithmetic_mean_Q_3(i,:) - part_1_confidence_intervals_Q_3(1,:));
        perr = abs(part_1_arithmetic_mean_Q_3(i,:) - part_1_confidence_intervals_Q_3(2,:));
        er = errorbar([1:1:8],part_1_arithmetic_mean_Q_3(i,:), nerr, perr);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none'; 
        hold off
    end
end

%% PROBLEM 3/4
letters = ['a', 'b','c', 'd', 'e', 'f', 'g', 'h'];

%aritmetic means
for i = 1:length(letters);
    data = a_mean_all_groups(groups,i+3);
    filename = strcat('part_1_arithmetic_mean_Q_4', letters(i), '.mat');
    save(filename, 'data');
end


%confidense intervalls
for i = 1:length(letters);
    data = conf_int(groups.group_1A_1B_2A_2B, i+3, 0.05);    
    filename = strcat('part_1_confidence_intervals_Q_4', letters(i), '.mat');
    save(filename, 'data'); 
end

pause(1);

%load data
for i = 1:length(letters);
    load(strcat('part_1_arithmetic_mean_Q_4', letters(i)));
    load(strcat('part_1_confidence_intervals_Q_4', letters(i)));
end







%% FUNCTION DECLARATIONS



function plsntns = pleasantness(p, a, ca, ch, v, m);

end

function CI = conf_int(group, qs, p);
    %returns the confidence intervall for all locations
    
    %returns the mean for a group for all locations
    number_of_locations = length(group(:,qs,1));
    
    for i = 1:number_of_locations;
        x = squeeze(group(i,qs,:)).*4./12.4;
        %standard deviation
        STD = std(x);

        % Standard error
        SEM = STD/sqrt(length(x)); 
        % T-score
        ts = tinv(1-p/2,length(x)-1);
       
        % Confidence interval
        CI(:,i) = [mean(x) - ts*SEM, mean(x) + ts*SEM];
    end
end


function mean_all = a_mean_all_groups(groups, qs);
    %returns an (n_groups x n_locations) array with aritmetic means 
    fn = fieldnames(groups);
    for i=1:length(fn)
        mean_all(i,:) = aritmetic_mean(groups.(fn{i}),qs);
    end
end 


function a_mean = aritmetic_mean(group, qs)
    %returns the mean for a group for all locations
    number_of_locations = length(group(:,qs,1));
    for i = 1:number_of_locations;
        a_mean(i) = mean(squeeze(group(i,qs,:)));
    end
end



