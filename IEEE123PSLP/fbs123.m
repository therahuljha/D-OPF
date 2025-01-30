% clc;
% clear all;
function [I_angle] = fbs123
bKVA = 1000;
bKV = 4.16/sqrt(3);
bZ = ((bKV)^2)*1000/bKVA;     %%5.768
load  linedata124C.txt;
fb = linedata124C(:,1);
tb = linedata124C(:,2);
G = graph(fb,tb);
T=dfsearch(G,1,'edgetonew');

%%% Formation of Z matrix for the three phase load
for k = 1:length(linedata124C(:,1))

     m = 3;
    for i = 1:3
        for j = 1:3
              if i <= j
                    Z{k}(i,j) = complex(linedata124C(k,m),linedata124C(k,m+6))/bZ;
                    m = m+1;
                else 
                    Z{k}(i,j) = Z{k}(j,i);
                end
           end
    end
end

load powerdata1.txt;
powerdata = powerdata1;
% load powerdata2.txt;
% powerdata = powerdata2;
% load powerdata3.txt;
% powerdata = powerdata3;
nb = size(powerdata,1);
%%% formation of complex power in matrix form (3,1)
for i = 1:nb
    m=2;
    for j = 1:3
        S{i}(j,1) =  complex(powerdata(i,m),powerdata(i,m+1))/bKVA;
        m = m+2;
    end
end

Vln = [1+0j; -0.5-0.866j; -0.5+0.866j ];   %%% defining the intial voltage

%%% converting load into constant impedance load
for i = 1:nb
    
    Z_load{i} = ((abs(Vln)).^2)./conj(S{i});
   
end

A = [1 0 0; 0 1 0; 0 0 1];      %% defing A for V2 =AV1-BI2   forward sweep

for i = 1:nb        %% defining intial voltages to all nodes
    V_old {i} = Vln;
end
iter = 1;
while iter <= 2
    
    
    for k = 1:nb
        I_l{k} = V_old{k}./Z_load{k} ;        %% calculating current for all the nodes
    end
    
nl = size(linedata124C,1);            %% size of line


%%% backward sweep current calculation I = V/Z

 for k = 1:nl
I_f{k} = zeros(3,1);            %%% intilizing all the flow current
 end
for k = nb:-1:2             %%% finding currents for all the end nodes

    child = find(k== T(:,1));          %% finding child for all the nodes
    Parent = find(k == T(:,2));
    if isempty(child)
        I_f{Parent} = I_l{k} ;
    end
end

for k = nb:-1:2                     %%finding current for each lineflow iter 1
    child = find(k== T(:,1));          %% finding child for all the nodes
    Parent = find(k == T(:,2));
    nchild(k)=size(find((k)==T(:,1)),1) ;
    if (nchild(k)==1)
        I_f{Parent} = I_l{k} ;
            
        for i=1:length(child)
            I_f{Parent} = I_f{Parent} + I_f{child(i)} ;
            
        end
   
    
    elseif (nchild(k)==2) | (nchild(k)==3)
        I_f{Parent} = I_l{k} ;
            
        for i=1:length(child)
            I_f{Parent} = I_f{Parent} + I_f{child(i)} ;
            
        end
   
    end
end

for k = nb:-1:2                          %%finding current for each lineflow iter 2
    child = find(k== T(:,1));          %% finding child for all the nodes
    Parent = find(k == T(:,2));
    nchild(k)=size(find((k)==T(:,1)),1) ;
    if (nchild(k)==1)
        I_f{Parent} = I_l{k} ;
            
        for i=1:length(child)
            I_f{Parent} = I_f{Parent} + I_f{child(i)} ;
            
        end
   
    
    elseif (nchild(k)==2) | (nchild(k)==3)
        I_f{Parent} = I_l{k} ;
            
        for i=1:length(child)
            I_f{Parent} = I_f{Parent} + I_f{child(i)} ;
            
        end
   
    end
end
for k = nb:-1:2                          %%finding current for each lineflow iter 3
    child = find(k== T(:,1));          %% finding child for all the nodes
    Parent = find(k == T(:,2));
    nchild(k)=size(find((k)==T(:,1)),1) ;
    if (nchild(k)==1)
        I_f{Parent} = I_l{k} ;
            
        for i=1:length(child)
            I_f{Parent} = I_f{Parent} + I_f{child(i)} ;
            
        end
   
    
    elseif (nchild(k)==2) | (nchild(k)==3)
        I_f{Parent} = I_l{k} ;
            
        for i=1:length(child)
            I_f{Parent} = I_f{Parent} + I_f{child(i)} ;
            
        end
   
    end
end

for k = nb:-1:2                          %%finding current for each lineflow iter 4
    child = find(k== T(:,1));          %% finding child for all the nodes
    Parent = find(k == T(:,2));
    nchild(k)=size(find((k)==T(:,1)),1) ;
    if (nchild(k)==1)
        I_f{Parent} = I_l{k} ;
            
        for i=1:length(child)
            I_f{Parent} = I_f{Parent} + I_f{child(i)} ;
            
        end
   
    
    elseif (nchild(k)==2) | (nchild(k)==3)
        I_f{Parent} = I_l{k} ;
            
        for i=1:length(child)
            I_f{Parent} = I_f{Parent} + I_f{child(i)} ;
            
        end
   
    end
end


%%% forward loop for voltages V2 = AV1-I12*Z12

for k = 2:nb
     Parent = find(k == T(:,2));
     child = find(k== T(:,1));
    V_old{1} = Vln;
    V_old{k} = A*V_old{T(Parent)} - Z{k-1}*I_f{(k-1)};
end

iter = iter+1;
end
iter;

for k = 1:nl
    I_angle{k} = radtodeg(angle(I_f{k}));
end
