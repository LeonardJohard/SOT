

%INIT
clear;
load('\\10.50.128.199\store\projects\SKILLS\ROW\ML\Energy analysis\Sessions\session_2_1.mat');
load('\\10.50.128.199\store\projects\SKILLS\ROW\ML\Energy analysis\Sessions\energies_2_1.mat');



is_map = 0;
recordstart = 0;
thisstroke = [];

sprint = session.sprint;
labellist = genvarname(labels(sprint));
    
for k = 1:32
   temp = data(sprint);
   eval([labellist{k} '= temp(k,:);']);
end                   


nrnodes = 60;

%Run data with classification

for i = 60000:80000
    if((energy.phases(i-1) == 4) && (energy.phases(i) == 1))
        
        if length(thisstroke) >= nrnodes

        if recordstart == 0  
            recordstart = 1;
            
        else
            if is_map == 0;
            %Init som
            
            somlength = length(thisstroke);
            
            map = somdtw(thisstroke, nrnodes);
            map.wl = 100;
            
            clf;
            scatter(map.nodes(1,:), map.nodes(2,:)), hold on
            plot(thisstroke(1,:), thisstroke(2,:))
            axis([0, 0.5, -1.2, 0.8]);
            pause(0.05);
            clf;
            end
            %Save stroke and train
            map = map.add(thisstroke);
            for ii = 1:1
            map = map.adaptDTW(2);
            map = map.adapt(1000);
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
    thisstroke = [thisstroke,[L_alpha(i),L_phi(i)]'];
end

                map = map.adaptDTW(5);
                map = map.adapt(1000);
                scatter(map.nodes(1,:), map.nodes(2,:), 'r'), hold on
                plot([map.nodes(1,:),map.nodes(1,1) ], [map.nodes(2,:),map.nodes(2,1)], 'r', 'LineWidth',3), hold on
