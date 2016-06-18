

%INIT
clear;
load('session_1_2.mat');
load('energies_1_2.mat');



is_map = 0;
recordstart = 0;
thisstroke = [];

sprint = session.sprint;
labellist = genvarname(labels(sprint));
    
for k = 1:32
   temp = data(sprint);
   eval([labellist{k} '= temp(k,:);']);
end                   




%Run data with classification

for i = 2000:(length(energy.phases))
    if((energy.phases(i-1) == 4) && (energy.phases(i) == 1))
        if recordstart == 0  
            recordstart = 1;
            
        else
            if is_map == 0;
            %Init som
            nrnodes = 20;
            somlength = length(thisstroke);
            
            map = som(thisstroke, nrnodes);
            
            
            
            scatter(map.nodes(1,:), map.nodes(2,:)), hold on
            plot(thisstroke(1,:), thisstroke(2,:))
            axis([0, 0.3, -1.2, 0.8]);
            pause(0.05);
            end
            %Save stroke and train
            map = map.add(thisstroke);
            map = map.adapt(100);
            is_map = 1;
                clf;
                scatter(map.nodes(1,:), map.nodes(2,:)), hold on
                plot([map.nodes(1,:),map.nodes(1,1) ], [map.nodes(2,:),map.nodes(2,1)], 'g'), hold on
                plot(thisstroke(1,:), thisstroke(2,:))
                axis([0, 0.3, -1.2, 0.8]);
                pause(0.05);

          
            
        end
        %Start new stroke
        thisstroke = []; 
    end
    thisstroke = [thisstroke,[L_alpha(i),L_phi(i)]'];
end