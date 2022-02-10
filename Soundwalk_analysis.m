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
    my_field = strcat('part_1_arithmetic_mean_Q_4', letters(i));
    aritm_means_Q_4.(my_field) = a_mean_all_groups(groups,i+3);
    data = aritm_means_Q_4.(my_field);
    filename = strcat('part_1_arithmetic_mean_Q_4', letters(i), '.mat');
    save(filename, 'data');
end


%confidense intervalls
for i = 1:length(letters);
    my_field = strcat('part_1_confidence_intervals_Q_4', letters(i));
    conf_ints_Q_4.(my_field) = conf_int(groups.group_1A_1B_2A_2B, i+3, 0.05);
    data = conf_ints_Q_4.(my_field);
    filename = strcat('part_1_confidence_intervals_Q_4', letters(i), '.mat');
    save(filename, 'data');
end


%% PROBLEM 5/6
p = aritm_means_Q_4.part_1_arithmetic_mean_Q_4a;
ch = aritm_means_Q_4.part_1_arithmetic_mean_Q_4b;
v = aritm_means_Q_4.part_1_arithmetic_mean_Q_4c;
u = aritm_means_Q_4.part_1_arithmetic_mean_Q_4d;
ca = aritm_means_Q_4.part_1_arithmetic_mean_Q_4e;
a = aritm_means_Q_4.part_1_arithmetic_mean_Q_4f;
e = aritm_means_Q_4.part_1_arithmetic_mean_Q_4g;
m = aritm_means_Q_4.part_1_arithmetic_mean_Q_4h;

part_1_pleasantness = pleasantness(p, a, ca, ch, v, m);
save('part_1_pleasantness.mat', 'part_1_pleasantness');

part_1_eventfulness = eventfulness(e, u, ca, ch, v, m);
save('part_1_eventfulness.mat', 'part_1_eventfulness');


%% PROBLEM 7
figure(4);
scatter(part_1_pleasantness, part_1_eventfulness);
legend('Location 1', 'Location 2', 'Location 3', 'Location 4', 'Location 5', 'Location 6', 'Location 7', 'Location 8');
grid on;
xlabel('Pleasantness');
ylabel('Eventfulness');
xlim([-1 1]);
ylim([-1 1]);



%% COMPARISON TO PREVIOUS YEARS

load('hrsv22PreviousSoundwalkData.mat');
clc

years17_22 = cat(3, previousYears(:).Q1234AM);
CI_years17_22 = cat(3, previousYears(:).Q1234CI);



markers = ['o','d', 'p', '>', '^', '<', 's'];
xlocations = {'Central Station', 'Trädgårdsföreningen', 'Stora Teatern', 'Stora Nygatan', 'Kungsportplatsen', 'Domkyrkan', 'Brunnsparken', 'Nordstan'};

%Question 1
figure(5);
for i = 1:7;
    if i == 6;
        err = abs(part_1_arithmetic_mean_Q_1(7,:) - part_1_confidence_intervals_Q_1(1,:));
        er = errorbar([0.6+i*0.1:1:7.6+i*0.1], part_1_arithmetic_mean_Q_1(7,:), err, err);
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
        hold on;
        h(i) = scatter([0.6+i*0.1:1:7.6+i*0.1], part_1_arithmetic_mean_Q_1(7,:), 75, markers(6), 'filled');
    elseif i == 7;
        h(i) = scatter([0.6+i*0.1:1:7.6+i*0.1], years17_22(:,1,6), 75, markers(i), 'filled');
        hold on;
    else;
        err = abs(years17_22(:,1,i) - CI_years17_22(:,1,2*i));
        er = errorbar([0.6+i*0.1:1:7.6+i*0.1], years17_22(:,1,i), err, err);
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
        hold on;
        h(i) = scatter([0.6+i*0.1:1:7.6+i*0.1], years17_22(:,1,i), 75, markers(i), 'filled');
    end;
    hold on;
end
hold off;
xlim([0 9]);
ylim([0 4]);
xticks(1:1:8);
xticklabels(xlocations);
grid on;
ylabel('Perceived loudness [0-4]');
legend(h, '2017', '2018', '2019', '2020', '2021', '2022', 'All years');
title('Answers to question 1: How loud is it here? - year 2017-2022');
ax = gca; 
ax.FontSize = 16;


%Question 2
figure(6);
for i = 1:7;
    if i == 6;
        err = abs(part_1_arithmetic_mean_Q_2(7,:) - part_1_confidence_intervals_Q_2(1,:));
        er = errorbar([0.6+i*0.1:1:7.6+i*0.1], part_1_arithmetic_mean_Q_2(7,:), err, err);
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
        hold on;
        h(i) = scatter([0.6+i*0.1:1:7.6+i*0.1], part_1_arithmetic_mean_Q_2(7,:), 75, markers(6), 'filled');
    elseif i == 7;
        h(i) = scatter([0.6+i*0.1:1:7.6+i*0.1], years17_22(:,2,6), 75, markers(i), 'filled');
        hold on;
    else;
        err = abs(years17_22(:,2,i) - CI_years17_22(:,2,2*i));
        er = errorbar([0.6+i*0.1:1:7.6+i*0.1], years17_22(:,2,i), err, err);
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
        hold on;
        h(i) = scatter([0.6+i*0.1:1:7.6+i*0.1], years17_22(:,2,i), 75, markers(i), 'filled');
    end;
    hold on;
end
hold off;
xlim([0 9]);
ylim([0 4]);
xticks(1:1:8);
xticklabels(xlocations);
grid on;
ylabel('How appropriate is the sound to the surrounding [0-4]');
legend(h, '2017', '2018', '2019', '2020', '2021', '2022', 'All years');
title('Answers to question 2: How appropriate is the sound to the surrounding? - year 2017-2022');
ax = gca; 
ax.FontSize = 16;


%Question 3
figure(6);
for i = 1:7;
    if i == 6;
        err = abs(part_1_arithmetic_mean_Q_3(7,:) - part_1_confidence_intervals_Q_3(1,:));
        er = errorbar([0.6+i*0.1:1:7.6+i*0.1], part_1_arithmetic_mean_Q_3(7,:), err, err);
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
        hold on;
        h(i) = scatter([0.6+i*0.1:1:7.6+i*0.1], part_1_arithmetic_mean_Q_3(7,:), 75, markers(6), 'filled');
    elseif i == 7;
        h(i) = scatter([0.6+i*0.1:1:7.6+i*0.1], years17_22(:,3,6), 75, markers(i), 'filled');
        hold on;
    else;
        err = abs(years17_22(:,3,i) - CI_years17_22(:,3,2*i));
        er = errorbar([0.6+i*0.1:1:7.6+i*0.1], years17_22(:,3,i), err, err);
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
        hold on;
        h(i) = scatter([0.6+i*0.1:1:7.6+i*0.1], years17_22(:,3,i), 75, markers(i), 'filled');
    end;
    hold on;
end
hold off;
xlim([0 9]);
ylim([0 4]);
xticks(1:1:8);
xticklabels(xlocations);
grid on;
ylabel('How often would you like to visit this place again [0-4]');
legend(h, '2017', '2018', '2019', '2020', '2021', '2022', 'All years');
title('Answers to question 3: How often would you like to visit this place again? - year 2017-2022');
ax = gca; 
ax.FontSize = 16;





%% FUNCTION DECLARATIONS



function P = pleasantness(p, a, ca, ch, v, m);
    P = ((sqrt(2)*(p-a))+(ca-ch)+(v-m))/(4+sqrt(8));
end

function E = eventfulness(e, u, ca, ch, v, m);
    E = ((sqrt(2)*(e-u))-(ca-ch)+(v-m))/(4+sqrt(8));
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


