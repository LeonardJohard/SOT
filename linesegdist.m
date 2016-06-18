function  [node, distance, rout] = linesegdist(p, nodevec)
%AB = line, C = point,
%http://www.codeguru.com/forum/showthread.php?t=194400

vlength = size(nodevec, 2) - 1;
r=-ones(1,vlength);
L = sqrt(sum((nodevec(:,2:end) - nodevec(:,1:end-1)).^2, 1));
distvec = [];
for i = 1:vlength
    AB(:,i) = nodevec(:,i+1) - nodevec(:,i);
    if sum(AB(:,i)) ~= 0;
    AC(:,i) = p - nodevec(:,i);
    r(i) = ( AC(:,i)' * AB(:,i) ) / L(i).^2;
    if (0 <= r(i))  && (r(i) <= 1)
        distvec(:,i) = AC(:,i) - (r(i)*AB(:,i)); %Position relative to line
  
    else
        if r(i) < 0
            distvec(:,i) = AC(:,i);
        else
            distvec(:,i) = p - nodevec(:,i+1); %BC
        end
    end
    end
end

if ~ isempty(distvec)

[~, segindex] = min(sum(distvec.^2,1));

node = segindex;
distance = distvec(:,segindex);
rout = r(segindex); %if r < 0 or > 1, then should not be used for variance calculations



else
    distance = -1;
    node = -1;
    rout = -1;
end