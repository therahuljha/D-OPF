%% This is the code to read the data from txt file and do the power flow

clc
clear all
filename = 'ieee123.txt';
% function [bus_data, branch_data] = read_file(filename)

%% read the load data for each node

fid = fopen(filename);
tline = fgetl(fid);
flag = 0;
while flag == 0
    if (strcmp(tline(1:16),'BUS DATA FOLLOWS')) == 1
        flag = 1;
    end
    tline = fgetl(fid);
end
flag = 0;
index_counter = 1;

while flag == 0
    if (strcmp(tline(1:4),'-999')) == 1
        flag = 1;
    end
    if flag == 0
     eachelement =   strsplit(tline);
     bus_name1 = eachelement(3);
     bus_name(index_counter,:) = bus_name1;
   string_bus_dataP  = eachelement(end-1);
   bus_dataP(index_counter,:) = string_bus_dataP;
   string_bus_dataQ  = eachelement(end);
   bus_dataQ(index_counter,:) = string_bus_dataQ;
    index_counter = index_counter +1;
    tline = fgetl(fid);
    end
end
bus_data = [bus_name (bus_dataP) (bus_dataQ)];
% for i = 1:size(bus_dataP)
%     bus_number(i,1) = i;
% end
% bus_data = [bus_number str2double(bus_dataP) str2double(bus_dataQ)];
%% read the branch data for each edge

tline = fgetl(fid);
tline = fgetl(fid);
flag = 0;
index_counter = 1;
while flag == 0
    if (strcmp(tline(1:4),'-999')) == 1
        flag = 1;
    end
    if flag == 0
        eachelement =   strsplit(tline);
        frombus(index_counter,:) =  eachelement(end-3);
        tobus(index_counter,:) =  eachelement(end-2);
        branch_dataR(index_counter,:) = eachelement(end-1);
        branch_dataX(index_counter,:) = eachelement(end);
        index_counter = index_counter + 1;
        tline = fgetl(fid);
    end
end
fclose(fid);
branch_data =[(frombus) (tobus) (branch_dataR) (branch_dataX)];


%% formation of tree 
G = graph(frombus,tobus);
nnodes = G.Nodes;
nedges = G.Edges;
T=dfsearch(G,'150','edgetonew');  
frbustree = T(:,1);
tobustree = T(:,2);

bkVA = 1000;                    % base kVA
bKV = 4.16/sqrt(3);                      % base kV
bZ = ((bKV)^2)*1000/bkVA;                % base Z
Vs = 1.0;

%% defining variables and solving the power flow
prob = optimproblem;

bigM = 100;
Pij = optimvar('Pij',size(nedges,1),1,'LowerBound',-bigM, 'UpperBound',bigM);
Qij = optimvar('Qij',size(nedges,1),1,'LowerBound',-bigM, 'UpperBound',bigM);
v = optimvar('v',size(nedges,1),1,'LowerBound',0.025, 'UpperBound',1.1);

prob.Objective = 0;

Pbal = optimconstr(size(nedges,1),1);
Qbal = optimconstr(size(nedges,1),1);
vequa = optimconstr(size(nedges,1),1);

for i = 1:size(nedges,1)

     lineindx_fr = find(strcmp(branch_data(:,1), frbustree(i)));
     lineindx_to = find(strcmp(branch_data(:,2), tobustree(i)));
      
     resistance = (str2double(branch_data(lineindx_to,3)))/bZ;
     reactance = (str2double(branch_data(lineindx_to,4)))/bZ;
     
   
     pa = find(strcmp(tobustree(:,1), tobustree(i)));
     ch = find(strcmp(frbustree(:,1), tobustree(i)));
     
     self = find(strcmp(T(:,2), tobustree(i)));
     parent = find(strcmp(T(pa,1),T(:,2) ));
     
     indx = find(strcmp(bus_data(:,1), tobustree(i)));
     Pkw = str2double(bus_data(indx,2))/bkVA;
     Qkw = str2double(bus_data(indx,3))/bkVA;
     
     %% P12- summation(P23) = PL2 and Q12- summation(Q23) = QL2-QC2 
     
     if isempty(ch) 
         Pbalance(i,1)  = Pij(pa)  -Pkw ;
         Qbalance(i,1)  = Qij(pa)  -Qkw ;
     else
         Pbalance(i,1) = Pij(pa) -Pij(ch(1)) -Pkw ;
         Qbalance(i,1) = Qij(pa) -Qij(ch(1)) -Qkw ;
         for j = 1:(size(ch,1)-1)
            Pbalance(i,1)  = Pbalance(i,1) - Pij(ch(j+1));
            Qbalance(i,1)  = Qbalance(i,1) - Qij(ch(j+1));
         end
     end
    Pbal(i,1) = Pbalance(i,1) == 0;
    Qbal(i,1) = Qbalance(i,1) == 0;
    
    %% volatge equation V2 -V1 + 2rP+2xP =0 
    
    if i == 1
        voltagebalance(i,1) =   v(self)- Vs + 2*resistance*Pij(pa) + 2*reactance*Qij(pa);
    else
        voltagebalance(i,1) =   v(self) - v(parent) + 2*resistance*Pij(pa) + 2*reactance*Qij(pa);
    end
     vequa (i,1) = voltagebalance(i,1) ==0;
     
end


prob.Constraints.Pbal = Pbal;
prob.Constraints.Qbal = Qbal;
prob.Constraints.vequa = vequa;

% showproblem(prob);

sol = solve(prob);
    
