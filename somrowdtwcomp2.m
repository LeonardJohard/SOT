

%INIT
clear;
clf;
%load('\\10.50.128.199\store\projects\SKILLS\ROW\ML\Energy analysis\Sessions\session_2_3.mat');
%load('\\10.50.128.199\store\projects\SKILLS\ROW\ML\Energy analysis\Sessions\energies_2_3.mat');


load('d:\users\leonard\Desktop\SOM\1_TOM\FreeRace.mat');

%load('d:\users\leonard\Desktop\SOM\session_1_2.mat');
%load('d:\users\leonard\Desktop\SOM\energies_1_2.mat');
rawdata{1} = ENERGY2_data([2,4],:);
rawdata{1} = addphase(rawdata{1});

load('d:\users\leonard\Desktop\SOM\1_TOM\Guidedrace.mat');

rawdata{2} = ENERGY2_data([2,4],:);
rawdata{2} = addphase(rawdata{2});

load('d:\users\leonard\Desktop\SOM\3_DAN\FreeRace.mat');

rawdata{3} = ENERGY2_data([2,4],:);
rawdata{3} = addphase(rawdata{3});

load('d:\users\leonard\Desktop\SOM\3_DAN\Guidedrace.mat');

rawdata{4} = ENERGY2_data([2,4],:);
rawdata{4} = addphase(rawdata{4});

tindex = [5000 69000 5000 60000 5000 60000 5000 60000];
for dataset = 1:4
subplot(2,2,dataset)
%Rename variables for compatibility with code
L_alpha = rawdata{dataset}(1,:);
L_phi = rawdata{dataset}(2,:);
energy.phases = rawdata{dataset}(3,:); 

is_map = 0;
recordstart = 0;
thisstroke = [];


%alphai = 
%phii =

varmax = 0.01;
nrnodes = 40;

%Run data with classification



for tpos = tindex(2*(dataset-1)+1):tindex(2*(dataset-1)+2)
    if((energy.phases(tpos) == 4) && (energy.phases(tpos+1) == 1))
        
        if length(thisstroke) >= nrnodes

        if recordstart == 0  
            recordstart = 1;
            
        else
            if is_map == 0;
            %Init som
            
            somlength = length(thisstroke);
            
            map = somdtw(thisstroke, nrnodes);
            map.wl = 50;
            map.ewmal = 0.998;
            

            plot(thisstroke(1,:), thisstroke(2,:))
            axis([0, 0.5, -1.2, 0.8]);
            pause(0.05);
            end
            %Save stroke and train
            map = map.add(thisstroke);
            for ii = 1:30
            map = map.adaptDTW(25);
            map = map.adapt(300);
            end
            is_map = 1;
                plot(thisstroke(1,:), thisstroke(2,:)), hold on
                axis([0, 0.5, -1.2, 0.8]);
                pause(0.05);

          
            
        end
        end
        %Start new stroke
        thisstroke = []; 
    end
    
    thisstroke = [thisstroke,[L_alpha(tpos),L_phi(tpos)]'];
end

map.segmean = zeros(1,nrnodes-1);
map.segvar = 0.0001*ones(1,nrnodes-1);
map.ewmal = 0.98;

savedmap{dataset} = map;
end

allvar = [];
map = {};
for plotnr = 1:4
    [map{plotnr} fullmean{plotnr} fullvar{plotnr}] = savedmap{plotnr}.fullmean;
    allvar = [allvar fullvar{plotnr}];
end

                        varmax = max(allvar);
                    varmin = min(allvar);
                    vardelta = varmax-varmin;
for plotnr = 1:4
                    
                subplot(2,2,plotnr)
 
                axis([0, 0.5, -1.2, 0.8]);
                pause(0.05);
                
                %scatter(map.nodes(1,:), map.nodes(2,:), 'b'), hold on;
                

                
                for i = 1:(map{plotnr}.nrnodes-1)

                    %colour =  [max(min( ((map.segvar(i)-varmin)/vardelta)^(1/2),1),0), -max(min(((map.segvar(i)-varmin)/vardelta)^(1/2), 1),0) + 1,0];
                    colour =  [max(min( ((fullvar{plotnr}(i)-varmin)/vardelta)^(1/2),1),0), 0,0];
                    thisplot = plot([map{plotnr}.nodes(1,i),map{plotnr}.nodes(1,i+1)], [map{plotnr}.nodes(2,i),map{plotnr}.nodes(2,i+1)]);
                    hold on;
                    set(thisplot,'Color', colour, 'LineWidth',3);
                end
end