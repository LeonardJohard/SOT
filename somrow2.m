

clear;
load('\\10.50.128.199\store\projects\SKILLS\ROW\ML\Energy analysis\Sessions\session_1_1.mat');
load('\\10.50.128.199\store\projects\SKILLS\ROW\ML\Energy analysis\Sessions\energies_1_1.mat');



is_map = 0;
recordstart = 0;
thisstroke = [];

sprint = session.sprint;
labellist = genvarname(labels(sprint));
    
for k = 1:32
   temp = data(sprint);
   eval([labellist{k} '= temp(k,:);']);
end                   




for i = 2:(length(energy.phases))
    if((energy.phases(i-1) == 4) && (energy.phases(i) == 1))
        if recordstart == 0
            recordstart = 1;
            
        else
            if is_map == 0;
                %Init som
                strlength = length(thisstroke);
                somsize = 20;
                somline = newsom([-0.3, 0.5; - 1.5, 1], [somsize, 1], 'gridtop', 'linkdist', 2, 3);
                somline.IW{1,1} = thisstroke(:,floor(1:((strlength-1)/(somsize-1)):strlength))';
                somline.adaptFcn = 'trains';
                %somline.layers{1,1}.distances = min(somline.layers{1,1}.distances, (somsize - somline.layers{1,1}.distances));
                is_map = 1;
            else
                [somline]= adapt(somline,thisstroke);
            
            
            end
        
        
        end
        thisstroke = [];
    end
    thisstroke = [thisstroke,[L_alpha(i),L_phi(i)]'];
end