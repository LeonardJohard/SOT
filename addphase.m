function  addeddata = addphase(rawdata)
        load phasenet.mat;
splinesize_back2 = 3; %Data point just before recognition point (RP)
splinesize_back = 30; %Second data point before RP, at longer distance
splinesize_forward =15; %Last data point after RP used
splinesize_forward2 =3; %Data point just after RP

    raw_size = length(rawdata);
    indata = zeros(10, raw_size - (splinesize_back + splinesize_forward));
    
   
    for i = 1:raw_size - (splinesize_back + splinesize_forward)

        %Collecting data from time window

        coord = rawdata(:,i:(i+splinesize_back + splinesize_forward));

        
        %Feature extraction
        xcoef = [coord(1, 1)-coord(1,splinesize_back), ...
            coord(1, splinesize_back-splinesize_back2)-coord(1,splinesize_back), ...
            mean(coord(1,splinesize_back-splinesize_back2:end))-coord(1,splinesize_back), ...
            coord(1,end)-coord(1,splinesize_back), ...
            coord(1,splinesize_back+splinesize_forward2)-coord(1,splinesize_back)];
        ycoef = [coord(2, 1)-coord(2,splinesize_back), coord(2, ...
            splinesize_back-splinesize_back2)-coord(2,splinesize_back), ...
            mean(coord(2,splinesize_back-splinesize_back2:end))-coord(2,splinesize_back), ...
            coord(2,end)-coord(2,splinesize_back), ...
            coord(2,splinesize_back+splinesize_forward2)-coord(2,splinesize_back)];
        coefs = [xcoef, ycoef]';
        indata(:,i) = coefs;
    end
    
    classtemp = sim(net, indata); %Resulting unfiltered classification
    [~, class] = max(classtemp);
    
%State machine filter
     newclass = zeros(1,length(class));
    i = 2;
    while i < length(class)-1;
         if class(i)  ~= class(i + 1) && mod(class(i + 1) - newclass(i-1), 4 ) <= 2

             if class(i + 2) == class(i + 1)
             newclass(i:i+19) = class(i + 1);
             i = i + 20;
                else
             newclass(i) = newclass(i-1);
             i = i + 1;
             end
 
         else
             newclass(i) = class(i+1);
             i = i + 1;
             end

     end
    
    class = [class(1) newclass(1:length(class)-1)];
    
    origdata_temp = [rawdata;zeros(1,splinesize_back), class,zeros(1,splinesize_forward)];
    addeddata = origdata_temp;