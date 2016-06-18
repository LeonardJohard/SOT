classdef somdtw
   properties
      nodes = [];
      nrnodes = 0;
      %The raw data series
      rawdata = {};
      indim = 2;
      nrseries = 0;
      segmean = [];
      segvar = [];
      %The DTW-transformed data points
      traindata = [];
      trnodes = []; %The closest node, per latest DTW update
      
      trindex = {}; %Index of the transformed point (for updating the transform)
      nrtrdata = 0;
      SampleNumber = 0;
      
      ewmal = 0.99; %lambda exponential weighting for ewma
      wl = 50; %Window length: Number of strokes used to calculate the mean
      
        
      

      dw = [];
      lr = 0.01;
      spreadmatrix = [];
      
   end
   methods 
      function  obj = somdtw(examplepos, nrnodes)
          obj.nrnodes = nrnodes;
          strlength = size(examplepos, 2);
          indexvec = floor(1:((strlength-1)/(obj.nrnodes-1)):strlength);
          obj.nodes = examplepos(:,indexvec);
          obj.segmean = zeros(1,nrnodes-1);
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
            obj.traindata = obj.traindata(:,obj.nrnodes + 1:obj.nrnodes*(obj.wl + 1))
            for i = length(obj.trindex)
                obj.trindex{i} = obj.trindex{i} - obj.nrnodes;
            end
            
            obj.nrtrdata = size(obj.traindata, 2);

            obj.nrseries = obj.nrseries - 1;
        end
        
        
      end
      
      function obj = adaptDTW(obj, steps) %Updates SOM assignments, + mean + var
          steps = min(steps, obj.nrseries);
          seriesselect = randi(obj.nrseries, 1, steps);
          for i = 1:steps
            pos = obj.rawdata{seriesselect(i)};
            [~, M] = somdistance(pos, obj.nodes);
            
            C = [1 1 1.0;1 0 1.0]; %Prevents same datapoint being assigned to more than 1 node. (Other solution?)
            
            [p,q,D,sc] = dpfast(M, C); %C can be removed for default matrix
            index = obj.trindex{seriesselect(i)};
            
        traindatatemp = zeros(obj.indim,obj.nrnodes);
        nodebins = zeros(1,obj.nrnodes);
        for i = 1:length(q)
           traindatatemp(:,q(i)) = traindatatemp(:,q(i)) + pos(:,i); 
           nodebins(1,q(i)) = nodebins(1,q(i)) + 1;
        end
        
        for i = 1:obj.indim
            traindatatemp(i, :) = traindatatemp(i, :) ./ nodebins;
        end
        
        obj = obj.ewma(pos,q); %updates mean + var

        obj.traindata(:,index(1):index(2)) = traindatatemp;
        obj.nrtrdata = size(obj.traindata, 2);
            
          end
      end
      
      
      
      
      function obj = adapt(obj, steps)
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
        
        %nodenr = obj.trnodes somdistance(pos, obj.nodes); %Get min distance node

        %Training
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
     
      function obj = ewma(obj, pos, q) %exponentiall weighted moving average
          dsegmean = zeros(obj.indim,obj.nrnodes-1);
          dsegvar = zeros(obj.indim,obj.nrnodes-1);
          neighbornodes = [];
          neighbornodes(1,:) = q-1;
          neighbornodes(2,:) = q+1;
          nodebins = zeros(1,obj.nrnodes-1); %Number of datapoints per segment
          nodesums = zeros(1,obj.nrnodes-1); %Sum per segment, for calculating mean/segment
          
          for i = 1:length(q)
            currentnodes = [];
            currentseg = [];
            
            %Collecting three adjacent nodes, if exist
            if obj.nrnodes > neighbornodes(1,i) && neighbornodes(1,i) > 0
                currentnodes = [currentnodes, obj.nodes(:,neighbornodes(1,i))];
                currentseg = [currentseg, neighbornodes(1,i)];
            end
                        
            currentnodes = [currentnodes, pos(:,q(i))];
            
            if obj.nrnodes > neighbornodes(2,i) && neighbornodes(2,i) > 0
                currentnodes = [currentnodes, obj.nodes(:,neighbornodes(2,i))];
                currentseg = [currentseg, q(i)];
            end
    
           [winnernode, distance, r] = linesegdist(pos(:,i), currentnodes);

           if 1 >= r && r >= 0 %Point within line segment
             nodebins(1,currentseg(winnernode)) = nodebins(1,currentseg(winnernode)) + 1;
             nodesums(1,currentseg(winnernode)) = nodesums(1,currentseg(winnernode)) + distance;
           end
           
          end
          
          %-Calculating mean distance per segment of curve finished
          
          nodemeanpos = [];
          sq_dist = [];
          for i = 1:(obj.nrnodes-1) %Updating stored avg + var
              if nodebins(1,i) > 0
                nodemeanpos(i) = nodesums(1,i) / nodebins(1,i);
                sq_dist(i) = (nodemeanpos(i) - obj.segmean(i))^2;
                obj.segmean(i) = obj.ewmal*obj.segmean(i) + (1-obj.ewmal)*nodemeanpos(i);
                obj.segvar(i) = obj.ewmal*obj.segvar(i) + (1-obj.ewmal)*sq_dist(i);
                'mean, var:'
                obj.segmean(i)
                obj.segvar(i)
              end
          end
          
      end
      
      
      
      function [obj totalnodeavg totalnodevar] = fullmean(obj)
          totalnodebins = zeros(1,obj.nrnodes-1); %Number of curves
          totalnodesums = zeros(1,obj.nrnodes-1); %As above, but for average of the means of the curves
          totalnodevarbins = zeros(1,obj.nrnodes-1); %Number of curves
          totalnodevarsums = zeros(1,obj.nrnodes-1); %As above, but for average of the means of the curves
          
          for i = 1:obj.nrseries
            pos = obj.rawdata{i};
            [dummy, M] = somdistance(pos, obj.nodes);
            
            C = [1 1 1.0;1 0 1.0]; %Prevents same datapoint being assigned to more than 1 node. (Other solution?)
            
            [p,q,D,sc] = dpfast(M, C); %C can be removed for default matrix
            
             obj = obj.ewma(pos,q); %updates mean + var

            qdebug{i} = q;
            
            dsegmean = zeros(obj.indim,obj.nrnodes-1);
          dsegvar = zeros(obj.indim,obj.nrnodes-1);
          neighbornodes = [];
          neighbornodes(1,:) = (q-1);
          neighbornodes(2,:) = (q+1);
          nodebins = zeros(1,obj.nrnodes-1); %Number of datapoints per segment

          nodesums = zeros(1,obj.nrnodes-1); %Sum per segment, for calculating mean/segment
          
          
          for ii = 1:length(q)
            currentnodes = [];
            currentseg = [];
            
            
            
            %Cop paste & edit from ewma code
            %Collecting three adjacent nodes, if exist
            if obj.nrnodes > neighbornodes(1,ii) && neighbornodes(1,ii) > 0
                currentnodes = [currentnodes, obj.nodes(:,neighbornodes(1,ii))];
                currentseg = [currentseg, neighbornodes(1,ii)];
            end
                        
            currentnodes = [currentnodes, pos(:,q(ii))];
            
            if obj.nrnodes > neighbornodes(2,ii) && neighbornodes(2,ii) > 0
                currentnodes = [currentnodes, obj.nodes(:,neighbornodes(2,ii))];
                currentseg = [currentseg, q(ii)];
            end
    
           [winnernode, distance, r] = linesegdist(pos(:,ii), currentnodes);
            rdebug{i}(ii) = r;
           if 1 >= r && r >= 0 %Point within line segment
             nodebins(1,currentseg(winnernode)) = nodebins(1,currentseg(winnernode)) + 1;
             nodesums(1,currentseg(winnernode)) = nodesums(1,currentseg(winnernode)) + distance;
           end
           
          end
          
          nodemean(i,:) = nan(1,obj.nrnodes-1);
            for ii = 1:(obj.nrnodes-1) %Mean per curve
               
              if nodebins(1,ii) > 0 
                nodemean(i,ii) = nodesums(ii) ./ nodebins(ii);
              else
                  nodemean(i,ii) = NaN;
              end
          
            end
          
             %Preparing mean over all curves
             for ii = 1:(obj.nrnodes-1) %Mean per curve
                 if (0 < nodebins(ii))
                    totalnodebins(ii) = totalnodebins(ii) + 1;
                    totalnodesums(ii) =  totalnodesums(ii) + (nodemean(i,ii));
                 end
             end
    
          end
              
          
            for ii = 1:(obj.nrnodes-1)
            if totalnodebins(ii) > 0
                totalnodeavg(ii) = totalnodesums(ii) ./ totalnodebins(ii);
            end
          
          
          
            end
            
            
            %Variance - calculated afterwards, when avg is determined
            
            for i = 1:obj.nrseries
                for ii = 1:obj.nrnodes-1
                    if ~isnan(nodemean(i,ii))
                        totalnodevarbins(ii) = totalnodevarbins(ii) + 1;
                        totalnodevarsums(ii) = totalnodevarsums(ii) + (nodemean(i,ii) - totalnodeavg(ii))^2;

                    end
                end
            end
            
            
            
          for ii = 1:(obj.nrnodes-1)
            if totalnodevarbins(ii) > 1
                totalnodevar(ii) = totalnodevarsums(ii) ./ (totalnodevarbins(ii)-1);
            end
          end
          
          

 
 
      end
      

          
   end% methods
end