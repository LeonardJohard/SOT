classdef somdtw
    
    %SOMDTW - 1-D Self-oganizing maps with dynamic time warping
    %
    %Provides methods for variance analysis of trajectories
    %
    %     Usage: Add() all trajectories, then run adapt() and adaptDTW() until
    % convergence. Get variance from segvar estimate or run fullmean for
    % complete analysis.
    %Functions:
    %somdtw(examplepos, nrnodes): Create a map from example
    %
    %add(pos): Add a trajectory to the stored trajectory set
    %
    %adaptDTW(obj, steps): Updates the assignment of trajectory points to
    %nodes. Draws a random stored trajectory and iterates for a number of 'steps'.
    %
    %adapt(obj, steps): Updates node positions to assigned trajecory
    %points.  Draws a random stored trajectory and iterates for a number of 'steps'.
    %
    %ewma(obj, pos, q): Update the moving average and variance estimates
    %(run adaptDTW to trigger this function externally)
    %
    %fullmean(obj): Calculates mean and variance of all stored trajectories
    %(run adapt and adaptDTW until convergence of node positions before)
    %
    %
    %
    %To do: elliptic variance, debugging of > 2 indim, built-in plotting
    %functions, reset C and find other solution for surjective mapping,
    %improved standard parameters
    %
    %By Leonard Johard (leonardjohard@yahoo.se)
    
   properties
      nodes = [];
      nrnodes = 0;
      %The raw data series
      rawdata = {};
      indim = 2; %Input dimensions
      nrseries = 0;
      
      segmean = []; %EWMA mean position estimate
      segvar = []; %EWMA variance estimate
      
      %The DTW-transformed data points
      traindata = [];
      trnodes = []; %The closest node, per latest DTW update
      
      trindex = {}; %Index of the transformed point (for updating the transform)
      nrtrdata = 0;
      SampleNumber = 0;
      
      ewmal = 0.99; %lambda exponential weighting for ewma
      wl = 50; %Window length: Maximum number of strokes saved in the map. The stored strokes are used for calculating adapt, adaptdtw and fullmean.
      
        
      

      dw = [];
      lr = 0.005;
      spreadmatrix = [];
      
   end
   methods 
      function  obj = somdtw(examplepos, nrnodes)
          obj.nrnodes = nrnodes;
          obj.indim = size(examplepos,1);
          strlength = size(examplepos, 2);
          indexvec = floor(1:((strlength-1)/(obj.nrnodes-1)):strlength);
          obj.nodes = examplepos(:,indexvec);
          obj.segmean = zeros(obj.indim,nrnodes-1);
          obj.segvar = 0.001*ones(1,nrnodes-1);
      end
       
      function obj = add(obj, pos)
        
        obj.nrseries = obj.nrseries + 1;
        obj.rawdata{obj.nrseries}  = pos;
       
        obj.trindex{obj.nrseries} = [size(obj.traindata, 2) + 1, size(obj.traindata, 2) + obj.nrnodes];
        
        C = [1 1 1.0;1 0 1.0]; %Prevents same datapoint being assigned to more than 1 node. (Other solution?)
        
        [~, M] = somdistance(pos, obj.nodes);
        [p,q,D,sc] = dpfast(M, C);
        

        traindatatemp = zeros(obj.indim,obj.nrnodes);
        nodebins = zeros(1,obj.nrnodes);
                
        for i = 1:length(q)
           traindatatemp(:,q(i)) = traindatatemp(:,q(i)) + pos(:,i); 
           nodebins(1,q(i)) = nodebins(1,q(i)) + 1;
        end
        
        for i = 1:obj.indim
            traindatatemp(i, :) = traindatatemp(i, :) ./ nodebins;
        end
        
        
        obj.trnodes = [obj.trnodes, 1:obj.nrnodes];
        obj.traindata = [obj.traindata, traindatatemp];
        obj.nrtrdata = size(obj.traindata, 2);
        
        %Reducing tainingdata to time window
        while obj.nrseries > obj.wl
            obj.rawdata = obj.rawdata(1 + 1:obj.wl + 1);
            obj.trindex = obj.trindex(1 + 1:obj.wl + 1);

            obj.trnodes = obj.trnodes(:,obj.nrnodes + 1:obj.nrnodes*(obj.wl + 1));
            obj.traindata = obj.traindata(:,obj.nrnodes + 1:obj.nrnodes*(obj.wl + 1));
            for i = length(obj.trindex)
                obj.trindex{i} = obj.trindex{i} - obj.nrnodes;
            end
            
            obj.nrtrdata = size(obj.traindata, 2);

            obj.nrseries = obj.nrseries - 1;
        end
        
        
      end
      
      function obj = adaptDTW(obj, steps) %Updates SOM DTW minimal point-to-node assignments, + mean + var. Draws trajectories randomly from the stored set.
          steps = min(steps, obj.nrseries);
          seriesselect = randi(obj.nrseries, 1, steps);
          for i = 1:steps
            pos = obj.rawdata{seriesselect(i)};
            [~, M] = somdistance(pos, obj.nodes);
            
            C = [1 1 1.0;1 0 1.0]; %Prevents same datapoint being assigned to more than 1 node. (To do: Other solution)
            
            [p,q,D,sc] = dpfast(M, C); %C can be removed for default matrix
            index = obj.trindex{seriesselect(i)};
            
        traindatatemp = zeros(obj.indim,obj.nrnodes);
        nodebins = zeros(1,obj.nrnodes);
        for j = 1:length(q)
           traindatatemp(:,q(j)) = traindatatemp(:,q(j)) + pos(:,j); 
           nodebins(1,q(j)) = nodebins(1,q(j)) + 1;
        end
        
        for j = 1:obj.indim
            traindatatemp(j, :) = traindatatemp(j, :) ./ nodebins;
        end
        
        obj = obj.ewma(pos,q); %updates mean + var (is computationally cheap here with q)

        obj.traindata(:,index(1):index(2)) = traindatatemp;
        obj.nrtrdata = size(obj.traindata, 2);
            
          end
      end
      
      
      
      
      function obj = adapt(obj, steps) %Draws trajectories randomly from the stored set.
        totalsteps = steps;
        stepcount = 0;
        if obj.nrtrdata == 0;
            error('Error: No trainingdata in adapt');
        end
        
        while stepcount < totalsteps
            steps = min(steps, obj.nrtrdata);
            stepcount = stepcount + steps;
        
        obj.dw = zeros(obj.indim,obj.nrnodes);
        
        dataselect = randi(obj.nrtrdata, 1, steps);
        pos = obj.traindata(:,dataselect);
        nodenr = obj.trnodes(:,dataselect);

        %Adapting positions
        for i=1:steps
            obj.dw(:,nodenr(i)) = obj.dw(:,nodenr(i)) + (pos(:,i) - obj.nodes(:, nodenr(i)));
            mini = mod(nodenr(i)-2,obj.nrnodes)+1;
            plusi = mod(nodenr(i),obj.nrnodes)+1;
            obj.dw(:,plusi) = obj.dw(:,plusi) + 0.3*(pos(:,i) - obj.nodes(:,plusi));
            obj.dw(:,mini) = obj.dw(:,mini) +  0.3*(pos(:,i) - obj.nodes(:,mini));
        

        
        end

        obj.nodes = obj.nodes + obj.lr*obj.dw;
        
        end
      end
     
      function obj = ewma(obj, pos, q) %exponentially weighted moving average
          neighbornodes = [];
          neighbornodes(1,:) = q-1;
          neighbornodes(2,:) = q+1;
          nodebins = zeros(obj.indim,obj.nrnodes-1); %Number of datapoints per segment
          nodesums = zeros(obj.indim,obj.nrnodes-1); %Sum per segment, for calculating mean/segment
          
          for i = 1:length(q)
            currentnodes = [];
            currentseg = [];
            
            %Collecting three adjacent nodes, if exist
            if obj.nrnodes >= neighbornodes(1,i) && neighbornodes(1,i) > 0
                currentnodes = [currentnodes, obj.nodes(:,neighbornodes(1,i))];
                currentseg = [currentseg, neighbornodes(1,i)];
            end
                        
            currentnodes = [currentnodes, obj.nodes(:,q(i))];
            
            if obj.nrnodes >= neighbornodes(2,i) && neighbornodes(2,i) > 0
                currentnodes = [currentnodes, obj.nodes(:,neighbornodes(2,i))];
                currentseg = [currentseg, q(i)];
            end
    
           [winnerseg, distance, r] = linesegdist(pos(:,i), currentnodes);

           if 1 >= r && r >= 0 %Point within line segment
             nodebins(:,currentseg(winnerseg)) = nodebins(:,currentseg(winnerseg)) + ones(obj.indim,1);
             nodesums(:,currentseg(winnerseg)) = nodesums(:,currentseg(winnerseg)) + distance;
           end
           
          end
          
          %-Calculating mean distance per segment of curve finished
          
          nodemeanpos = [];
          sq_dist = [];
          for i = 1:(obj.nrnodes-1) %Updating stored avg + var
              if nodebins(1,i) > 0
                nodemeanpos(:,i) = nodesums(:,i) ./ nodebins(:,i);
                sq_dist(i) = sum((nodemeanpos(:,i) - obj.segmean(:,i)).^2);
                obj.segmean(:,i) = obj.ewmal*obj.segmean(:,i) + (1-obj.ewmal)*nodemeanpos(:,i);
                obj.segvar(i) = obj.ewmal*obj.segvar(i) + (1-obj.ewmal)*sq_dist(i);
              end
          end
          
      end
      
      
      
      function [obj totalnodeavg totalnodevar] = fullmean(obj)
          
          %Init
          totalnodebins = zeros(obj.indim,obj.nrnodes-1); %Number of curves
          totalnodesums = zeros(obj.indim,obj.nrnodes-1); %Aaverage of the means of the curves
          totalnodevarbins = zeros(1,obj.nrnodes-1); %Number of curves
          totalnodevarsums = zeros(1,obj.nrnodes-1); %Average of the means of the curves
          
          %Main loop
          for i = 1:obj.nrseries
              pos = obj.rawdata{i};
              [~, M] = somdistance(pos, obj.nodes);
              
              C = [1 1 1.0;1 0 1.0]; %Prevents same datapoint being assigned to more than 1 node. (Other solution?)
              
              %DTW
              [p,q,D,sc] = dpfast(M, C); %C can be removed for default matrix
              
              obj = obj.ewma(pos,q); %updates ewma mean + var
              
              
              neighbornodes = [];
              neighbornodes(1,:) = (q-1);
              neighbornodes(2,:) = (q+1);
              nodebins = zeros(obj.indim,obj.nrnodes-1); %Number of datapoints per segment
              
              nodesums = zeros(obj.indim,obj.nrnodes-1); %Sum per segment, for calculating mean/segment
              
              
              for ii = 1:length(q)
                  currentnodes = [];
                  currentseg = [];
                  
                  
                  %Collecting three adjacent nodes to optimum, if exist
                  if (obj.nrnodes > neighbornodes(1,ii)) && (neighbornodes(1,ii) > 0)
                      currentnodes = [currentnodes, obj.nodes(:,neighbornodes(1,ii))];
                      currentseg = [currentseg, neighbornodes(1,ii)];
                  end
                  
                  currentnodes = [currentnodes, obj.nodes(:,q(ii))];
                  
                  if (obj.nrnodes > neighbornodes(2,ii)) && (neighbornodes(2,ii) > 0)
                      currentnodes = [currentnodes, obj.nodes(:,neighbornodes(2,ii))];
                      currentseg = [currentseg, q(ii)];
                  end
                  
                  %Decide closest line segment.
                  [winnerseg, distance, r] = linesegdist(pos(:,ii), currentnodes);

                  
                  if 1 >= r && r >= 0 %Point within line segment
                      nodebins(:,currentseg(winnerseg)) = nodebins(:,currentseg(winnerseg)) + ones(obj.indim,1);
                      nodesums(:,currentseg(winnerseg)) = nodesums(:,currentseg(winnerseg)) + distance;
                  end
                  
              end
              
              nodemean{i}(:,:) = nan(obj.indim,obj.nrnodes-1);
              
              %Mean position of curve
              for ii = 1:(obj.nrnodes-1) 
                  
                  if nodebins(1,ii) > 0
                      nodemean{i}(:,ii) = nodesums(:,ii) ./ nodebins(:,ii);
                      if ii == 1
                          1;
                      end
                  else
                      nodemean{i}(:,ii) = NaN;
                  end
                  
              end
              
              %Mean over all curve averages
              for ii = 1:(obj.nrnodes-1) %Mean per curve
                  if (0 < nodebins(1,ii))
                      totalnodebins(:,ii) = totalnodebins(:,ii) + 1;
                      totalnodesums(:,ii) =  totalnodesums(:,ii) + (nodemean{i}(:,ii));
                  end
              end
              
          end
              
          
            for ii = 1:(obj.nrnodes-1)
            if totalnodebins(1,ii) > 0
                totalnodeavg(:,ii) = totalnodesums(:,ii) ./ totalnodebins(:,ii);
            end
          
          
          
            
            end
            
            %Variance over all curve averages
            
            for i = 1:obj.nrseries
                for ii = 1:obj.nrnodes-1
                    if ~isnan(nodemean{i}(1,ii))
                        totalnodevarbins(ii) = totalnodevarbins(ii) + 1;
                        totalnodevarsums(ii) = totalnodevarsums(ii) + sum((nodemean{i}(:,ii) - totalnodeavg(:,ii)).^2);

                    end
                end
            end
            
            
            
          for ii = 1:(obj.nrnodes-1)
            if totalnodevarbins(ii) > 1
                totalnodevar(ii) = totalnodevarsums(ii) ./ (totalnodevarbins(ii)-1);
            end
          end
          
          

 
 
      end
      

          
   end
end