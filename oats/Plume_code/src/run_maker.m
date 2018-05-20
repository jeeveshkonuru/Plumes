function datout = run_maker(datin)
%Andrew Rzeznik 8/10/216
%Input: a cell array of structs, where each field in the struct is a cell array. 
%Output: a cell array of structs, which has been lengthed by replacing each
%single struct from the input with a set of structs corresponding to its
%first non-unit length cell array entry.



numstructs=length(datin); % Current number of structs in the cell
unpacknum=ones(1,numstructs); %number of runs each struct unpacks to
newstructs{numstructs}={};

%% Pre-Sanititon of inputs to be worked with
for i=1:numstructs
    names=fieldnames(datin{i});
    numnames=length(names);
    for j=1:numnames
        if iscell(datin{i}.(names{j})) 
        elseif ischar(datin{i}.(names{j}))
            datin{i}.(names{j})={datin{i}.(names{j})};
        elseif isfloat (datin{i}.(names{j}))
            datin{i}.(names{j}) = num2cell(datin{i}.(names{j}));
        end
    end
end

%% Main rollout loop
for i=1:numstructs
    curstruct=datin{i};
    names=fieldnames(curstruct);
    numnames=length(names);
    for j=1:numnames
        exnum=length(curstruct.(names{j})); %length of current struct field
        if exnum > 1 
            tempstruct=curstruct;
            tempholdstructs{exnum}={};
            for k=1:exnum
                tempstruct.(names{j}) = curstruct.(names{j})(k);
                tempholdstructs{k}=run_maker({tempstruct});
            end
            newstructs{i}=tempholdstructs;
            unpacknum(i)=exnum;
            break
        end
        newstructs{i}={{curstruct}};
    end
end


%% Recombination loop
index=1;

for i=1:numstructs
    for j=1:unpacknum(i)
        for k=1:length(newstructs{i}{j})
            datout{index}=newstructs{i}{j}{k};
            index=index+1;
        end
    end
end

%% De-sanitization of all struct fields to basic types
% Code to de-cell all of the variables in all the structs
for i=1:(index-1)
    names=fieldnames(datout{i});
    numnames=length(names);
    for j=1:numnames
        if iscell(datout{i}.(names{j})) 
            datout{i}.(names{j})=datout{i}.(names{j}){1};
        end
    end
end
