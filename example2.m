%This file compares continuous EWMA variance and trajectory estimate with
%convergence from full data
%
%In the plot red indicates high variance, black low variance.
%

%INIT
clear;

load('exampledata.mat');

nrnodes = 40;
nriter = 5;

map = somdtw(stroke{1}, nrnodes);
map.wl = 90;
map.ewmal = 0.99;
lr = 0.1;


%Online calculation:

for i = 30:60 %Set i for determining the used training set, i = 89 max
        %Plot
        clf;
        plot(stroke{i}(1,:), stroke{i}(2,:)), hold on


        
        %scatter(map.nodes(1,:), map.nodes(2,:), 'b'), hold on;
        
        
        
        for j = 1:(map.nrnodes-1)
            varmax = max(map.segvar);
            varmin = min(map.segvar);
            vardelta = varmax-varmin;
            %colour =  [max(min( ((map.segvar(i)-varmin)/vardelta)^(1/2),1),0), -max(min(((map.segvar(i)-varmin)/vardelta)^(1/2), 1),0) + 1,0];
            colour =  [max(min( ((map.segvar(j)-varmin)/vardelta)^(1/2),1),0), 0,0];
            thisplot = plot([map.nodes(1,j),map.nodes(1,j+1)], [map.nodes(2,j),map.nodes(2,j+1)]);
            set(thisplot,'Color', colour, 'LineWidth',3);
        end
        axis([0, 0.5, -1.2, 0.8]);
        pause(0.05);
        
        
        
        
        
        map = map.add(stroke{i}); %Add stroke to SOM data set
        
        %Update map a few iterations
        for j = 1:nriter
            map = map.adaptDTW(1);
            map = map.adapt(20);
        end
        
end


%Full convergence:
[map fullmean fullvar] = map.fullmean;


lr = 0.01; %Better precision for full variance



for j = 1:100
    map = map.adaptDTW(5);
    map = map.adapt(20);
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





