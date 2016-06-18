%This file plots the trajectories, the mean trajectory and the full
%variance after convergence
%
%In the plot red indicates high variance, black low variance.
%

%INIT
clear;
clf;


load('exampledata.mat');

nrnodes = 40;
nriter = 20;

map = somdtw(stroke{1}, nrnodes);
map.wl = 50;
map.ewmal = 0.998;


%Rename variables
for i = 2:10
        %Plot

        plot(stroke{i}(1,:), stroke{i}(2,:)), hold on


        
        %scatter(map.nodes(1,:), map.nodes(2,:), 'b'), hold on;
        
        

        
        
        
        
        
        map = map.add(stroke{i}); %Add stroke to SOM data set
        pause(0.005);
        %Update map a few iterations

        
end


[map fullmean fullvar] = map.fullmean;


map.lr = 0.05;



for j = 1:10 %Can be optimized muchbetter
    map = map.adaptDTW(25);
    map = map.adapt(300);
end

%Code below is for plot

axis([0, 0.5, -1.2, 0.8]);
pause(0.05);

for i = 1:(map.nrnodes-1)
    varmax = max(fullvar);
    varmin = min(fullvar);
    vardelta = varmax-varmin;
    %colour =  [max(min( ((map.segvar(i)-varmin)/vardelta)^(1/2),1),0), -max(min(((map.segvar(i)-varmin)/vardelta)^(1/2), 1),0) + 1,0];
    colour =  [max(min( ((fullvar(i)-varmin)/vardelta)^(1/2),1),0), 0,0];
    thisplot = plot([map.nodes(1,i),map.nodes(1,i+1)], [map.nodes(2,i),map.nodes(2,i+1)]);
    hold on;
    set(thisplot,'Color', colour, 'LineWidth',3);
end





