classdef som
   properties
      nodes = [];
      nrnodes = 0;
      traindata = [];
      nrtrdata = 0;
      SampleNumber = 0;
      Stress
      Strain
      Modulus = 0;
      dw = [];
      lr = 0.01;
      spreadmatrix = [];
   end
   methods 
      function  obj = som(examplepos, nrnodes)
          obj.nrnodes = nrnodes;
          strlength = size(examplepos, 2);
          indexvec = floor(1:((strlength-1)/(obj.nrnodes-1)):strlength);
          obj.nodes = examplepos(:,indexvec);
      end
       
      function obj = add(obj, pos)
          
        obj.traindata = [obj.traindata, pos];
        obj.nrtrdata = size(obj.traindata, 2);
      end
      
      function obj = adapt(obj, steps)
        totalsteps = steps;
        stepcount = 0;
        if steps == 0;
            cout << "error, no trainingdata";
        end
        
        while stepcount < totalsteps
            steps = min(steps, obj.nrtrdata);
            stepcount = stepcount + steps;
        
        obj.dw = zeros(2,obj.nrnodes);
        
        dataselect = randi(obj.nrtrdata, 1, steps);
        pos = obj.traindata(:,dataselect);
        
        
        [nodenr,dummy] = somdistance(pos, obj.nodes); %Get min distance node
        
        %Training
        for i=1:steps
            obj.dw(:,nodenr(i)) = obj.dw(:,nodenr(i)) + (pos(:,i) - obj.nodes(:, nodenr(i)));
            mini = mod(nodenr(i)-2,obj.nrnodes)+1;
            plusi = mod(nodenr(i),obj.nrnodes)+1;
            obj.dw(:,plusi) = obj.dw(:,plusi) + 0.5*(pos(:,i) - obj.nodes(:,plusi));
            obj.dw(:,mini) = obj.dw(:,mini) +  0.5*(pos(:,i) - obj.nodes(:,mini));
        end
        
        obj.nodes = obj.nodes + obj.lr*obj.dw;
        
        end
     end

   end% methods
end