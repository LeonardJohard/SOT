

%INIT
clear;
%load('\\10.50.128.199\store\projects\SKILLS\ROW\ML\Energy analysis\Sessions\session_2_2.mat');
%load('\\10.50.128.199\store\projects\SKILLS\ROW\ML\Energy analysis\Sessions\energies_2_2.mat');
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

varmax = 0.01;
nrnodes = 50;

%Run data with classification

for tpos = 2000:(length(energy.phases))
    if((energy.phases(tpos-1) == 4) && (energy.phases(tpos) == 1))
        
        if length(thisstroke) >= nrnodes

        if recordstart == 0  
            recordstart = 1;
            
        else
            if is_map == 0;
            %Init som
            
            somlength = length(thisstroke);
            
            map = somdtw(thisstroke, nrnodes);
            map.wl = 30;
            map.ewmal = 0.96;
            map.lr = 0.03;
            
            clf;
            scatter(map.nodes(1,:), map.nodes(2,:)), hold on;
            plot(thisstroke(1,:), thisstroke(2,:));
            axis([0, 0.5, -1.2, 0.8]);
            pause(0.01);
            end
            %Save stroke and train
            map = map.add(thisstroke);
            for ii = 1:1
            map = map.adaptDTW(2);
            map = map.adapt(10);
            end
            is_map = 1;
                clf;
                plot(thisstroke(1,:), thisstroke(2,:)), hold on;
                scatter(map.nodes(1,:), map.nodes(2,:),'b');
                length(thisstroke(1,:))
                for i = 1:(map.nrnodes-1)
                    
                   varmax = max(map.segvar);
                   varmin = min(map.segvar);
                   vardelta = varmax-varmin;
                   colour =  [max(min( ((map.segvar(i)-varmin)/vardelta)^(1/2),1),0), 0,0];
                   thisplot = plot([map.nodes(1,i),map.nodes(1,i+1)], [map.nodes(2,i),map.nodes(2,i+1)]);
                   set(thisplot,'Color', colour, 'LineWidth',3);
                end
                axis([0, 0.5, -1.2, 0.8]);
                pause(0.001);

          
            
        end
        end
        %Start new stroke
        thisstroke = []; 
    end
    
    thisstroke = [thisstroke,[L_alpha(tpos),L_phi(tpos)]'];
end
for i = 1:1
                clf;
                map.lr = 0.001;
                map = map.adaptDTW(75);

                
                
                plot(thisstroke(1,:), thisstroke(2,:));
                scatter(map.nodes(1,:), map.nodes(2,:), 'b'), hold on;
                for i = 1:(map.nrnodes-1)
                    varmax = max(map.segvar);
                    varmin = min(map.segvar);
                    vardelta = varmax-varmin;
                    colour =  [max(min( ((map.segvar(i)-varmin)/vardelta)^(1/2),1),0), 0,0];
                    thisplot = plot([map.nodes(1,i),map.nodes(1,i+1)], [map.nodes(2,i),map.nodes(2,i+1)]), hold on;
                    set(thisplot,'Color', colour, 'LineWidth',3);
                end
                axis([0, 0.5, -1.2, 0.8]);
                pause(0.01);
end

