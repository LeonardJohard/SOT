

%INIT
clear;
%load('\\10.50.128.199\store\projects\SKILLS\ROW\ML\Energy analysis\Sessions\session_2_3.mat');
%load('\\10.50.128.199\store\projects\SKILLS\ROW\ML\Energy analysis\Sessions\energies_2_3.mat');


load('1_TOM\FreeRace.mat');

%load('d:\users\leonard\Desktop\SOM\session_1_2.mat');
%load('d:\users\leonard\Desktop\SOM\energies_1_2.mat');

    
rawdata = ENERGY2_data([2,4],:);
rawdata = addphase(rawdata);

%Rename variables for compatibility with code
L_alpha = rawdata(1,:);
L_phi = rawdata(2,:);
energy.phases = rawdata(3,:); 

exists_map = 0;
recordstart = 0;
thisstroke = [];


%alphai = 
%phii =

varmax = 0.01;
nrnodes = 40;

%Run data with classification

for tpos = 40000:42000
    if((energy.phases(tpos) == 4) && (energy.phases(tpos+1) == 1))
        
        if length(thisstroke) >= nrnodes

        if recordstart == 0  
            recordstart = 1;
            
        else
            if exists_map == 0;
                %Init som

                somlength = length(thisstroke);

                map = somdtw(thisstroke, nrnodes);
                map.wl = 50;
                map.ewmal = 0.998;

                clf;
                plot(thisstroke(1,:), thisstroke(2,:))
                axis([0, 0.5, -1.2, 0.8]);
                pause(0.05);
            end
            
            
            
            
            %Save stroke and train
            map = map.add(thisstroke);
            for ii = 1:1
            map = map.adaptDTW(25);
            map = map.adapt(300);
            end
            exists_map = 1;
                plot(thisstroke(1,:), thisstroke(2,:)), hold on
                axis([0, 0.5, -1.2, 0.8]);
                pause(0.05);

          
            
        end
        end
        
        %Start new stroke
        thisstroke = []; 
    end
    
    thisstroke = [thisstroke,[L_alpha(tpos),L_phi(tpos)]']; %Add point to current stroke
end


%map.segmean = zeros(1,nrnodes-1);
%map.segvar = 0.0001*ones(1,nrnodes-1);
%map.ewmal = 0.98;
                
[map fullmean fullvar] = map.fullmean;
                axis([0, 0.5, -1.2, 0.8]);
                pause(0.05);
                
                %scatter(map.nodes(1,:), map.nodes(2,:), 'b'), hold on;
                

                
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