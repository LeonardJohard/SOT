% Class for managing structured data expressed in matrix form and
% with labels.
%
% Construction: structy(value N x M,labels,column ranges,optional units of
% measure)
% w = structy([1,2,3;4 5 7;10 12 13],{'a','c'},{1:2,3},{'m','s'});
%   Two fields 'a', and 'c', of which first is two dimensional
% w.a 
%   Access of field
% size(w)
%   Size of data [3 x 3]
% size(w,'a')
%   Size of field [3 x 2]
% labels(w)
%   List of labels
% units(w)
%   List of units
% Operations:
% [a;b;...]
%   Concatenates data vertically (like along time)
% [a,b,...]
%   Joins data horizontally. Requires non conflicting names (use addprefix)
% project(a,field or field list)
%   Extract fields as new structure. Field can be negative for removing fields
% project(a,field or field list,samedata)
%   Allows to keep the original data
% asstruct(a)
%   Return struct with labels and data. It can be used for compatibility storage on MAT files
% w = structy(storedstruct); 
%    Another construction using structure from asstruct(a)
% sort(a,fieldname)
%   Sort by field
% resample(a,timefield,step_or_timerrange,mode)
%    Uses timefield as time information and performs data interpolation
% get(a,fieldname)
% set(a,fieldname,value)
% a.name
% a(range)
% a(range).name
% a(range).name(columns)
% a.name(range)
% a.name(range,columns)
% a.name = data
% a(range) = data
% a(:) = [];  % removes all data
% rename(a,oldfield,newfield)
% isfield(a,name)
% addprefix(a,prefix)
% double(a) => compacted data
%
% Emanuele Ruffaldi 2009/12/15 (updated 2009/01/19)
classdef structy
    properties(GetAccess = private)
        labs = {} % list of fields names
        c = struct()  % field with indices 
        u = {}  % units
        d = []  % data
        totsize = 0  % total size for optimization
    end    
    properties 
        meta = struct()
    end    
    methods
        function obj = structy(c,labels,indices,units)
            if nargin > 0
                % copy constructor
                if isa(c,'structy')
                    obj.labs = c.labs;
                    obj.d = c.d;
                    obj.c = c.c;
                    obj.u = c.u;
                    obj.totsize = c.totsize;
                elseif isstruct(c)
                    obj.labs = c.labels;
                    obj.u = c.units;
                    obj.c = struct();
                    fn = c.labels;
                    totsize = 0;         
                    rows = 0;
                    for I=1:length(fn)
                        q = getfield(c,fn{I});
                        if I == 1
                            rows = size(q,1);
                        elseif rows ~= size(q,1)
                            error('not same rows');
                        end
                        obj.c = setfield(obj.c,fn{I},(totsize+1):(totsize+size(q,2)));
                        totsize = totsize + size(q,2);
                    end                        
                    d = zeros(rows,totsize);
                    tgt = 1;
                    for I=1:length(fn)
                        q = getfield(c,fn{I});
                        d(:,tgt:(tgt+size(q,2)-1)) = q;
                        tgt = tgt + size(q,2);
                    end
                    obj.d = d;                        
                    obj.totsize = totsize;                
                else
                    % empty data
                    if nargin == 1
                        obj.d = c;
                    else
                        % labels specified
                        if length(labels) ~= size(labels,1)
                            labels = labels';
                        end
                        % analyze data and labels
                        data = c;
                        if nargin < 3 || isempty(indices)
                            if size(labels,2) == 1
                                indices = num2cell(1:length(labels));    
                            else
                                indices = labels(:,2);
                            end
                        elseif length(indices) ~= length(labels)
                            error('Number of labels is different than number of indices %d vs %d',length(indices),length(labels));
                        end
                        if max(max(cell2mat(indices))) > size(data,2)
                            error('Data Columns are less than label columns');
                        end
                        F = struct;
                        k = 0;
                        for I=1:length(labels)
                            F.(labels{I}) = indices{I};
                            k = k + length(indices{I});
                        end
                        obj.c = F;
                        obj.d = data;
                        obj.labs = labels;               
                        if nargin >= 4
                            obj.u = units;
                        else
                            obj.u = cell(1,length(labels));
                            for I=1:length(labels)
                                obj.u{1,I} = '';
                            end
                        end
                        obj.totsize = k;
                    end
                end
            end
        end

        % list of indices of low level data that are valid
        function idx = indices(r)
            fn = fieldnames(r.c); 
            idx = []; % up to r.totsize
            for I=1:length(fn)
                idx = [idx r.c.(fn{I})]; 
            end
        end 
        
        % sort(obj,fieldname,mode)
        % sort(obj,fieldname)
        function obj = sort(obj,field,mode)
            if nargin < 3
                mode = 'ascend';
            end             
            [Y,I] = sort(obj.get(field),1,mode);
            obj.d = obj.d(I,:);            
        end
        
        % concatenates two lists with same fields
        function r = vertcat(a, varargin)
            r = a;
            for J=1:length(varargin)                
                b = varargin{J};
                if isempty(b)
                    continue;
                end
                if ~isequal(a.c,b.c)
                    error('not same labels');
                end
                if size(a.d,2) ~= size(b.d,2)
                    a = compact(a);
                    b = compact(b);
                    if size(a.d,2) ~= size(b.d,2)
                        error('not same data size also after compactification');
                    end
                end
                r = structy(a);
                r.d = [a.d;b.d];
            end
        end
        
        function r = iscompact(a)
            r = size(a.d,2) == a.totsize;
        end
        
        function r = length(a)
            r = size(a.d,1);
        end
        
        function r = isempty(a)
            r = isempty(a.d);
        end
        
        function r = max(a)
            if iscompact(a)
                r = max(a.d);
            else
                r = max(a.d(:,indices(a)));
            end
        end
        
        function r = min(a)
            if iscompact(a)
                r = min(a.d);
            else
                r = min(a.d(:,indices(a)));
            end
        end
        
        function r = median(a)
            if iscompact(a)
                r = median(a.d);
            else
                r = median(a.d(:,indices(a)));
            end
        end
        
        function r = std(a)
            if iscompact(a)
                r = std(a.d);
            else
                r = std(a.d(:,indices(a)));
            end
        end
        
        function r = mean(a)
            if iscompact(a)
                r = mean(a.d);
            else
                r = mean(a.d(:,indices(a)));
            end
        end
        % size(obj, fieldname)
        % size(obj, idx) as in matrices (rows or columns)
        function r = size(a,b)
            if nargin == 2
                if ischar(b)
                        r = [size(a.d,1),length(a.c.(b))];
                else
                    r = [size(a.d,1),a.totsize]; % size should be here
                    r = r(b);
                end
            else
                r = [size(a.d,1),a.totsize]; % size should be here
            end            
        end
        
        % adds a prefix to all fields
        function obj = addprefix(obj,prefix)
            for I=1:length(obj.labs)
                oldname = obj.labs{I};
                newname = [prefix obj.labs{I}];
                idx = getfield(obj.c,oldname);
                obj.c = setfield(rmfield(obj.c,oldname),newname,idx);
                obj.labs{I} = newname;
            end            
        end
        
        % renames a field
        function obj = rename(obj,oldname,newname)
            if ~strcmp(oldname,newname)
                mo = strmatch(oldname,obj.labs);                
                mn = strmatch(newname,obj.labs);               
                if isempty(mo)
                    error('name not found');
                elseif ~isempty(mn)
                    error('new name exists');
                else
                    obj.labs{mo} = newname;
                    idx = getfield(obj.c,oldname);
                    obj.c = setfield(rmfield(obj.c,oldname),newname,idx);
                end
            end
        end
        
        % returns true if field exists
        function r = isfield(obj,name)
            r = ~isempty(strmatch(name,obj.labs));
        end
        
        
        % joints horizontally lists if they same number of rows
        % names should be not conflicting
        function r = horzcat(a,varargin)
            r = a;
            for J=1:length(varargin)
                b = varargin{J};
                if isempty(b)
                    continue
                end
                if size(a.d,1) ~= size(b.d,1)
                    error('not same number of elements');
                end
                if nargin  == 2
                    prefix = '';
                end
                d= [a.d,b.d];
                c = a.c;
                labels = a.labs;
                J = size(a.d,2);
                for I=1:length(b.labs)
                    bname = b.labs{I};
                    name = bname;
                    if isfield(c,name) 
                        error('duplicate name in structure, use suffix for second'); 
                    end
                    labels{end+1} = name;                
                    c = setfield(c,name,b.c.(bname)+J); 
                end
                r = structy();
                r.labs = labels;
                r.d = d;
                r.c = c;
                r.totsize = a.totsize+b.totsize;
                r.u = [a.u,b.u]; 
            end
        end 
       
       % returns pure data
       function c = data(obj)
            c = obj.d;
       end 
        
       % returns units
       function c = units(obj)
            c = obj.u;
       end 
       
       % returns indices
       function c = indexes(obj)
            c = obj.c;
       end 
        
       % return labels
       function c = labels(obj)
            c = obj.labs;
       end 
        
        % return packed data
        function c = double(obj)
            c = obj.d(:,obj.indices());
        end % double
        
        function b = subsref(a,s)
        % SUBSREF Implement obj([1 ...])
        % A.name
        % A(range)
        % A(range).name
        % A(range).name(columns)
        % A.name(range)           
        % A.name(range,columns)
            ind = [];
            subi = [];
            if length(s) == 3
                if strcmp(s(1).type,'()') && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
                    if length(s(1).subs) ~= 1 || length(s(3).subs) ~= 1
                        error('only supported A(range).name(columns)');
                    end
                    subi = s(3).subs{1};
                    s = s(1:2);                    
                end
            end
            if length(s) == 2
                if strcmp(s(1).type,'()') && strcmp(s(2).type,'.')
                    if length(s(1).subs) > 1
                        error('Supported only single dimensional (time) ranges');
                    else
                        ind = s(1).subs{1};
                        s = s(2);
                    end
                elseif strcmp(s(1).type,'.') && strcmp(s(2).type,'()')
                    switch length(s(2).subs)
                    case 1
                        ind = s(2).subs{1};
                        s = s(1);
                    case 2
                        ind = s(2).subs{1};
                        subi = s(2).subs{2};
                        s = s(1);
                    otherwise
                        error('Supported only single dimensional (time) ranges');
                    end
                else
                    error('Supported only var(range).name or var.name(range)');
                end
                switch s.type
                case '()'
                    error('Double range subscript is not allowed');
                case '.'
                   if isfield(a.c,s.subs)
                       if isempty(ind)
                           b = a.d(:,a.c.(s.subs)); 
                       else
                           if isempty(subi)
                               b = a.d(ind,a.c.(s.subs)); 
                           else
                               q = a.c.(s.subs);
                               b = a.d(ind,q(subi)); 
                           end
                       end
                   else
                       error(['''' s.subs '''' ' is not a property'])
                   end
                otherwise
                    error('Specify value for x as obj(x)')
                end
            elseif length(s) == 1
                switch s.type
                case '()'
                    ind = s.subs{:};
                    b = structy(a);
                    b.d = b.d(ind,:);
                case '.'
                   if isfield(a.c,s.subs)
                       b = a.d(:,a.c.(s.subs)); 
                   else
                       a.labels
                       error(['''' s.subs '''' ' is not a property'])
                   end
                otherwise
                    error('Specify value for x as obj(x)')
                end
            else
                 error('Supported only var(range) or var.name or var(range).name');     
            end
        end
        function a = subsasgn(a,s,b)
        % a(range) = value
        % a.name = value
            if length(s) > 1
                error('only var(range) = values or var.name = values are allowed');
            end
            switch s.type
            case '()'
                % compactify if needed
                if size(a.d,2) ~= a.totsize
                    a = compact(a);
                end
                ind = s.subs{:};
                a.d(ind,:) = b;
            case '.'
               switch(s.subs)
                   otherwise
                   if isfield(a.c,s.subs)
                       a.d(:,a.c.(s.subs)) = b; 
                   else
                       error(['''' s.subs '''' ' is not a property'])
                   end
               end
            otherwise
                error('Specify value for x as obj(x)')
            end
        end % subsref

        % sets the field
        function a = set(a,name,b)
           if isfield(a.c,name)
               a.d(:,a.c.(name)) = b; 
           else
               error(['''' name '''' ' is not a property'])
           end
        end % subsref
        
        % gets the field explicitly
        function r = get(a,name)
           if isfield(a.c,name)
               r = a.d(:,a.c.(name)); 
           else
               error(['''' name '''' ' is not a property'])
           end
        end % subsref
        
        % compact data structure
        function r = compact(obj)
            if obj.totsize ~= size(obj.d,2)
                fn = fieldnames(obj.c);
                idx = [];
                r = structy();
                for I=1:length(fn)
                    q = obj.c.(fn{I});
                    idx = [idx q]; 
                    r.c = obj.c.(fn{I})((length(idx)-length(q)+1):length(idx)); 
                end
                r.d = obj.d(idx,:);            
                r.totsize = length(idx);
            else
                r  = obj;
            end
        end 
        
        % reduces the data by selected labels. 
        %
        % xlabels are the labels, used as prefix. If xlabels starts with ~ it is removed
        % samedata (default 0) preserves original data
        function r = project(obj,xlabels,samedata)
            if nargin == 2
                samedata = 0;
            end
            if ischar(xlabels)
                xlabels = {xlabels};
            end
            labels = {};
            for I=1:length(xlabels)
                x = xlabels{I};
                if iscell(x)
                    labels = [labels x];
                else
                    labels{end+1} = x;
                end
            end
            r = structy();

            ql = obj.labs';
            r.labs = {};
            r.c = struct;
            
            % first compute resulting label indices
            posindices = [];
            negindices = 1:length(ql);
            for I=1:length(labels)
                l = labels{I}
                if l(1) == '~'
                    idx = strmatch(l(2:end),ql);
                    if isempty(idx)
                        error(['cannot find label ',l]);
                    end
                    negindices = setdiff(negindices,idx);
                else
                    idx = strmatch(l,ql);
                    if isempty(idx)
                        error(['cannot find label ',l]);
                    end
                    posindices = [posindices,idx];
                end
            end
            if length(negindices) ~= length(ql)
                if ~isempty(posindices)
                    error('cannot mix positive with negative selectors in project');
                end
                selected = negindices;
            else
                selected = posindices;
            end
            if samedata == 0
                OI = [];
                for J=1:length(selected)
                    ll = ql{selected(J)}; %label
                    id = getfield(obj.c,ll); 
                    uu = obj.u{selected(J)};
                    r.labs{end+1} = ll;
                    r.u{end+1} = uu;
                    OI = [OI id];
                    r.c.(ll) = length(OI); % same pos
                end
                r.d = obj.d(:,OI);
                r.totsize = length(OI);
            else
                k = 0;
                for J=1:length(selected)
                    uu = obj.u{selected(J)};
                    ll = ql{selected(J)}; %label
                    id = getfield(obj.c,ll); 
                    k  = k + length(id);
                    r.labs{end+1} = ll;
                    r.u{end+1} = uu;
                    r.c.(ll) = id; 
                end
                r.d = obj.d;
                r.totsize = k;
            end
        end
        
        % uses interp1 for interpolating data in the time range using given timefield
        % resample(a,timefield,step_or_timerange,mode)
        % step_or_timerange can be new timestep and it is computed in range min:max of the time field
        %
        % all the other fields are being interpolated using the given mode
        function a = resample(a,timefield,step_or_timerange,mode)
            if nargin < 4
                mode = [];
            end
            % compact for not performing excessive works
            a = compact(a);
            if ischar(timefield)
                idx = getfield(a.c,timefield);
            else
                idx = timefield;
            end
            if length(idx) > 1
                error('Time is vector');
            end
            timedata = get(a,timefield);
            if length(step_or_timerange) == 1
                mt = min(timedata);
                Mt = max(timedata);
                outtimedata = (mt:step_or_timerange:Mt)';
            else
                outtimedata = step_or_timerange;
                if size(outtimedata,2) ~= 2
                    outtimedata = outtimedata';
                end
            end
            if idx > 1
                pre = interp1(timedata,a.d(:,1:(idx-1)),outtimedata,mode);
            else
                pre = [];
            end
            if idx < size(a.d,1)
                post = interp1(timedata,a.d(:,(idx+1):end),outtimedata,mode);
                if idx > 1
                    a.d = [pre,outtimedata,post];
                else
                    a.d = [outtimedata,post];
                end
            else            
                if idx > 1
                    a.d = [pre,outtimedata];
                else
                    a.d = outtimedata;
                end
            end
        end
        
        function out = asstruct(r)
            out = struct();
            out.labels = r.labs;
            out.units = r.u;
            fn = fieldnames(r.c); 
            for I=1:length(fn)
                out = setfield(out,fn{I},r.d(:,getfield(r.c,fn{I})));
            end
        end
            
        function disp(obj)
            disp(sprintf('structy sized %d/%d rows %d',obj.totsize,size(obj.d,2),size(obj.d,1)));
            disp(obj.labels);
	    end % disp
    end % methods
        
end % classdef
