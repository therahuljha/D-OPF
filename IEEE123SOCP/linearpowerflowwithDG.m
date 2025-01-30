 function [VA, VB, VC, x_cvr,fval_cvr,TableA,TableB,TableC,VolttableA,VolttableB,VolttableC] = linearpowerflowwithDG(t)
% clc;
% clear all;
% t=97;

load linedata124C.txt               %%load linedata
load branchphA.txt;                 %%load connection among the line for phase A,B,C
load branchphB.txt;
load branchphCC.txt;
lineA= branchphA;
lineB= branchphB;
lineC= branchphCC;
branch=linedata124C;

% load powerdata1.txt;
% powerdata = powerdata1;

% load powerdata2.txt;
% powerdata = powerdata2;

load powerdata3.txt;
powerdata = powerdata3;

load loadshape15.txt;
mult = loadshape15(:,2);
mult1 = loadshape15(:,3);


fb = linedata124C(:,1);
tb = linedata124C(:,2);
fbA = lineA(:,1);
tbA = lineA(:,2);
G = graph(fb,tb);
nbusA =  length(lineA) +1;
GA = graph(fbA,tbA);
VA = dfsearch(GA,1);

fbB = lineB(:,1);
tbB = lineB(:,2);

nbusB =  length(lineB) +1;
GB = graph(fbB,tbB);
VB = dfsearch(GB,1);


fbC = lineC(:,1);
tbC = lineC(:,2);
nbusC =  length(lineC) +1;
GC = graph(fbC,tbC);
VC = dfsearch(GC,1);

tnbA = length(fbA);
tnbB = length(fbB);
tnbC = length(fbC);

nbus = [nbusA , nbusB , nbusC];
nb = size((powerdata),1);
% nb = 128;
bKVA = 1000;
bKV = 4.16/sqrt(3);
bZ = ((bKV)^2)*1000/bKVA;
rAA = ((branch(:,3))/bZ);
rAB = ((branch(:,4))/bZ);
rAC = ((branch(:,5))/bZ);
rBB = ((branch(:,6))/bZ);
rBC = ((branch(:,7))/bZ);
rCC = ((branch(:,8))/bZ);
xAA = ((branch(:,9))/bZ);
xAB = ((branch(:,10))/bZ);
xAC = ((branch(:,11))/bZ);
xBB = ((branch(:,12))/bZ);
xBC = ((branch(:,13))/bZ);
xCC = ((branch(:,14))/bZ);
PLA = (powerdata(:,2).*mult(t))/bKVA;
PLB = (powerdata(:,4).*mult(t))/bKVA;
PLC = (powerdata(:,6).*mult(t))/bKVA;
QLA = (powerdata(:,3).*mult(t))/bKVA;
QLB = (powerdata(:,5).*mult(t))/bKVA;
QLC = (powerdata(:,7).*mult(t))/bKVA;
QCA = powerdata(:,8)/bKVA;
QCB = powerdata(:,9)/bKVA;
QCC = powerdata(:,10)/bKVA;
PGA = (powerdata(:,11).*mult1(t))/bKVA;
PGB = (powerdata(:,12).*mult1(t))/bKVA;
PGC = (powerdata(:,13).*mult1(t))/bKVA;
T=dfsearch(G,1,'edgetonew');
TA=dfsearch(GA,1,'edgetonew');
TB=dfsearch(GB,1,'edgetonew');
TC=dfsearch(GC,1,'edgetonew');

%%% formation of impedance matrix
RAA = zeros(nb);
XAA  = zeros(nb);
RAB = zeros(nb);
XAB  = zeros(nb);
RAC = zeros(nb);
XAC  = zeros(nb);
RBB = zeros(nb);
XBB  = zeros(nb);
RBC = zeros(nb);
XBC  = zeros(nb);
RCC = zeros(nb);
XCC  = zeros(nb);
for i = 1:(nb-1)
RAA(fb(i), tb(i)) = rAA(i);
RAA(tb(i) ,fb(i)) = RAA(fb(i), tb(i)) ;
XAA(fb(i), tb(i))= xAA(i);
XAA(tb(i) ,fb(i)) = XAA(fb(i), tb(i)) ;
RAB(fb(i), tb(i)) = rAB(i);
RAB(tb(i) ,fb(i)) = RAB(fb(i), tb(i)) ;
XAB(fb(i), tb(i))=  xAB(i);
XAB(tb(i) ,fb(i)) = XAB(fb(i), tb(i)) ;
RAC(fb(i), tb(i)) = rAC(i);
RAC(tb(i) ,fb(i)) = RAC(fb(i), tb(i)) ;
XAC(fb(i), tb(i))=  xAC(i);
XAC(tb(i) ,fb(i)) = XAC(fb(i), tb(i)) ;
RBB(fb(i), tb(i)) = rBB(i);
RBB(tb(i) ,fb(i)) = RBB(fb(i), tb(i)) ;
XBB(fb(i), tb(i))= xBB(i);
XBB(tb(i) ,fb(i)) = XBB(fb(i), tb(i)) ;
RBC(fb(i), tb(i)) = rBC(i);
RBC(tb(i) ,fb(i)) = RBC(fb(i), tb(i)) ;
XBC(fb(i), tb(i))= xBC(i);
XBC(tb(i) ,fb(i)) = XBC(fb(i), tb(i)) ;
RCC(fb(i), tb(i)) = rCC(i);
RCC(tb(i) ,fb(i)) = RCC(fb(i), tb(i)) ;
XCC(fb(i), tb(i))= xCC(i);
XCC(tb(i) ,fb(i)) = XCC(fb(i), tb(i)) ;
end

 Ap = 1:nbusA-1;                          % defining the unknowns for phaseA
 Aq = nbusA:2*(nbusA-1);
 Av = 2*(nbusA):3*(nbusA-1)+1;

 Bp = (3*(nbusA-1)+1)+1:(3*(nbusA-1)+1)+(nbusB-1);                    % defining the unknowns for phaseB
 Bq = (3*(nbusA-1)+1)+nbusB :(3*(nbusA-1)+1) + 2*(nbusB-1);
 Bv = (3*(nbusA-1)+1)+2*(nbusB):(3*(nbusA-1)+1) + 3*(nbusB-1) +1 ;


Cp = ((3*(nbusA-1)+1) + 3*(nbusB-1) +1)+1:((3*(nbusA-1)+1) + 3*(nbusB-1) +1) + (nbusC-1);             % defining the unknowns for phaseC
Cq = ((3*(nbusA-1)+1) + 3*(nbusB-1) +1) + nbusC :((3*(nbusA-1)+1) + 3*(nbusB-1) +1) + 2*(nbusC-1);
Cv = ((3*(nbusA-1)+1) + 3*(nbusB-1) +1)+2*(nbusC):((3*(nbusA-1)+1) + 3*(nbusB-1) +1) + 3*(nbusC-1) +1 ;


TableA = [TA(:,1) TA(:,2) Ap'  Aq'   Av'];
TableB = [TB(:,1) TB(:,2) Bp'  Bq'   Bv'];
TableC = [TC(:,1) TC(:,2) Cp'  Cq'   Cv'];

Da = 2*(nbusA)-1:3*(nbusA-1)+1;
Db = (3*(nbusA-1)+1)+2*(nbusB)-1:(3*(nbusA-1)+1) + 3*(nbusB-1) +1;
Dc = ((3*(nbusA-1)+1) + 3*(nbusB-1) +1)+2*(nbusC)-1:((3*(nbusA-1)+1) + 3*(nbusB-1) +1) + 3*(nbusC-1) +1 ;

VolttableA = Da';
VolttableB = Db';
VolttableC = Dc';


for i = 1:(nb)                       % number of child of a node for Phase A
zA(i)=size(find((i)==TA(:,1)),1) ;
row = find(i == TA(:,1));
end
zA;

for i = 1:(nb)                        % number of child of a node for Phase B
zB(i)=size(find((i)==TB(:,1)),1) ;
row = find(i == TB(:,1));
end
zB;

for i = 1:(nb)                            % number of child of a node for Phase C
zC(i)=size(find((i)==TC(:,1)),1) ;
row = find(i == TC(:,1));
end
zC;

%%% formation of equality matrix. Here, the column of the Aeq will
%%% be maximum number of unkonowns and rows will be number of equality
%%% equations. 
%%
totalvarpqv = VolttableC(end);
Vs = 1.1025;
% CVR_P = 0.6;                %%% CVR factor for P = 0.6
% CVR_Q = 3;                  %%% CVR factor for Q = 3

CVR_P = 2.0;                %%% CVR factor for P = 0.6
CVR_Q = 2.0; 

Aeq = zeros(totalvarpqv-3,totalvarpqv);
beq = zeros(totalvarpqv-3,1);
 for i =2:nb
    k = zA(i);
    row = find(i == TA(:,1));
 
  if isempty(row)
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
 if isempty(ParentA)  
Aeq = Aeq;

 elseif ((isempty(ParentB)) && (isempty(ParentC)) ) 
           Poc = find(TA(ParentA,1) == TA(:,1));
   Aeq(ParentA,TableA(ParentA,3))= 1;                                   %% if variable is present it will be 1
   Aeq(ParentA,TableA(ParentA,5))= -(CVR_P/2)*PLA(TableA(ParentA,2));         %% As the load is voltage control hence  it will be of this form where the CVR value is taken 
   Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,5)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(Poc(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2));
   
 elseif isempty(ParentB)   
Poc = find(TA(ParentA,1) == TA(:,1));
   Aeq(ParentA,TableA(ParentA,3))= 1;
   Aeq(ParentA,TableA(ParentA,5))= -(CVR_P/2)*PLA(TableA(ParentA,2));
   Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,5)))= -(CVR_Q/2)*QLA(TableA(ParentA,2));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(Poc(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2));
   
elseif isempty(ParentC) 
 Poc = find(TA(ParentA,1) == TA(:,1));
   Aeq(ParentA,TableA(ParentA,3))= 1;
   Aeq(ParentA,TableA(ParentA,5))= -(CVR_P/2)*PLA(TableA(ParentA,2));
   Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,5)))= -(CVR_Q/2)*QLA(TableA(ParentA,2));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(Poc(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2));
     
else
    
  Poc = find(TA(ParentA,1) == TA(:,1));
   Aeq(ParentA,TableA(ParentA,3))= 1;
   Aeq(ParentA,TableA(ParentA,5))= -(CVR_P/2)*PLA(TableA(ParentA,2));
   Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,5)))= -(CVR_Q/2)*QLA(TableA(ParentA,2));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(Poc(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2));
   
 end
 
    else  
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    
 if isempty(ParentA)  
Aeq = Aeq;

 elseif ((isempty(ParentB)) && (isempty(ParentC)) ) 
      Poc = find(TA(ParentA,1) == TA(:,1));
   Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
       end
   Aeq(ParentA,TableA(ParentA,5))= -(CVR_P/2)*PLA(TableA(ParentA,2));
   Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
      Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -  1;
    end
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,5)))= -(CVR_Q/2)*QLA(TableA(ParentA,2));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(Poc(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2));
   
 elseif isempty(ParentB)   
Poc = find(TA(ParentA,1) == TA(:,1));
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
       end
    Aeq(ParentA,TableA(ParentA,5))= -(CVR_P/2)*PLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -  1;
    end
    
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,5)))= -(CVR_Q/2)*QLA(TableA(ParentA,2));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(Poc(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2)) ;
   
elseif isempty(ParentC) 
 Poc = find(TA(ParentA,1) == TA(:,1));
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
       end
    Aeq(ParentA,TableA(ParentA,5))= -(CVR_P/2)*PLA(TableA(ParentA,2)); 
    
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -  1;
    end
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,5)))= -(CVR_Q/2)*QLA(TableA(ParentA,2));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(Poc(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
   beq(ParentA+(nbusA-1)) =  (1-(CVR_Q/2))*QLA(TableA(ParentA,2));
     
else
    
  Poc = find(TA(ParentA,1) == TA(:,1));
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
       end
    Aeq(ParentA,TableA(ParentA,5))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -  1;
    end
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,5)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(Poc(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2));
end
    end
 Aeq(3*(nbusA-1)+1,VolttableA(1)) = 1;
 beq(3*(nbusA-1)+1) = Vs;
 
 end
Aeq;

for i =2:nb
    k = zB(i) ;
    row = find(i == TB(:,1));
 
  if isempty(row)
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
 if isempty(ParentB)  
Aeq = Aeq;

 elseif ((isempty(ParentA)) && (isempty(ParentC)) )  
     Poc = find(TB(ParentB,1) == TB(:,1));
   Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
   Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(CVR_P/2)*PLB(TableB(ParentB,2));
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,5)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2));
   
 elseif isempty(ParentA)   
 Poc = find(TB(ParentB,1) == TB(:,1));
   Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
   Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(CVR_P/2)*PLB(TableB(ParentB,2));
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,5)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   beq(3*(nbusA-1)+1+ParentB) = (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2));
   
elseif isempty(ParentC) 
 Poc = find(TB(ParentB,1) == TB(:,1));
   Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
   Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(CVR_P/2)*PLB(TableB(ParentB,2));
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,5)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2));   
else
    Poc = find(TB(ParentB,1) == TB(:,1));
   Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
   Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(CVR_P/2)*PLB(TableB(ParentB,2));
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,5)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2));

 end
 
 
    else  
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
 if isempty(ParentB)  
    Aeq = Aeq;

 elseif ((isempty(ParentA)) && (isempty(ParentC)) ) 
     Poc = find(TB(ParentB,1) == TB(:,1));
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
       end
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(CVR_P/2)*PLB(TableB(ParentB,2));    
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB++(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -  1;
    end
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,5)))= -(CVR_Q/2)*QLB(TableB(ParentB,2));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) =  (1-(CVR_Q/2))*QLB(TableB(ParentB,2));
    
 elseif isempty(ParentA)   
     Poc = find(TB(ParentB,1) == TB(:,1));
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
       end
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(CVR_P/2)*PLB(TableB(ParentB,2));   
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -  1;
    end
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,5)))= -(CVR_Q/2)*QLB(TableB(ParentB,2));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   beq(3*(nbusA-1)+1+ParentB)=(1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2));
  
elseif isempty(ParentC) 
Poc = find(TB(ParentB,1) == TB(:,1));
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
       end
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(CVR_P/2)*PLB(TableB(ParentB,2));    
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -  1;
    end
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,5)))= -(CVR_Q/2)*QLB(TableB(ParentB,2));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   beq(3*(nbusA-1)+1+ParentB)=(1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2));
  
     
else
    
  Poc = find(TB(ParentB,1) == TB(:,1));
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
       end
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(CVR_P/2)*PLB(TableB(ParentB,2));    
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -  1;
    end
   Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,5)))= -(CVR_Q/2)*QLB(TableB(ParentB,2));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2));
   
end
  end
   
 Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1,VolttableB(1)) = 1;
 beq(3*(nbusA-1)+1+3*(nbusB-1)+1) = Vs;
 end
Aeq;
for i =2:nb
    k = zC(i) ;
    row = find(i == TC(:,1));
 
  if isempty(row)
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
 if isempty(ParentC)  
    Aeq = Aeq;

 elseif ((isempty(ParentA)) && (isempty(ParentB)) ) 
   Poc = find(TC(ParentC,1) == TC(:,1));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(CVR_P/2)*PLC(TableC(ParentC,2));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,5)))= -(CVR_Q/2)*QLC(TableC(ParentC,2));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2));
    
 elseif isempty(ParentA)   
 Poc = find(TC(ParentC,1) == TC(:,1));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(CVR_P/2)*PLC(TableC(ParentC,2)); 
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,5)))= -(CVR_Q/2)*QLC(TableC(ParentC,2));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2));
   
elseif isempty(ParentB) 
     Poc = find(TC(ParentC,1) == TC(:,1));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(CVR_P/2)*PLC(TableC(ParentC,2));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,5)))= -(CVR_Q/2)*QLC(TableC(ParentC,2));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2)); 
   
else
    
 Poc = find(TC(ParentC,1) == TC(:,1));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(CVR_P/2)*PLC(TableC(ParentC,2));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,5)))= -(CVR_Q/2)*QLC(TableC(ParentC,2));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2));
 end
 

    else  
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
 if isempty(ParentC)  
    Aeq = Aeq;

 elseif ((isempty(ParentA)) && (isempty(ParentB)) ) 
     Poc = find(TC(ParentC,1) == TC(:,1));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
       end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(CVR_P/2)*PLC(TableC(ParentC,2));     
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -1;
    end
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,5)))= -(CVR_Q/2)*QLC(TableC(ParentC,2));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2));
    
 elseif isempty(ParentA)   
Poc = find(TC(ParentC,1) == TC(:,1));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
       end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(CVR_P/2)*PLC(TableC(ParentC,2));   
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -  1;
    end
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,5)))= -(CVR_Q/2)*QLC(TableC(ParentC,2));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2));
   
elseif isempty(ParentB) 
Poc = find(TC(ParentC,1) == TC(:,1));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
       end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -  1;
    end
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,5)))= -(CVR_Q/2)*QLC(TableC(ParentC,2));
   
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) =  (1-(CVR_Q/2))*QLC(TableC(ParentC,2));
     
else
    
 Poc = find(TC(ParentC,1) == TC(:,1));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
       end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(CVR_P/2)*PLC(TableC(ParentC,2));     
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -1;
    end
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,5)))= -(CVR_Q/2)*QLC(TableC(ParentC,2));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(Poc(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) =  (1-(CVR_Q/2))*QLC(TableC(ParentC,2));
end
  end
  
  Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+3*(nbusC-1)+1,VolttableC(1)) = 1;
  beq(3*(nbusA-1)+1+3*(nbusB-1)+1+3*(nbusC-1)+1) = Vs;
end
 Aeq;
 beq;

%% addition of DG Q variable

totalpfvar = size(Aeq,2);
p =1;

for k = 1:size(powerdata,1)
    if powerdata(k,11) ~= 0 && powerdata(k,12) ~= 0 && powerdata(k,13) ~= 0
         ParentA = find(k == TA(:,2));
         Aeq(TableA(ParentA ,4),totalpfvar+p) = 1;
         p = p+1;
         ParentB = find(k == TB(:,2));
         Aeq(TableB(ParentB ,4),totalpfvar+p) = 1;
         p = p+1;
         ParentC = find(k == TC(:,2));
         Aeq(TableC(ParentC ,4),totalpfvar+p) = 1;
         p = p+1;
    elseif powerdata(k,11) ~= 0
        ParentA = find(k == TA(:,2));
        Aeq(TableA(ParentA ,4),totalpfvar+p) = 1;
        p = p+1;
       
    elseif powerdata(k,12) ~= 0
        ParentB = find(k == TB(:,2));
        Aeq(TableB(ParentB ,4),totalpfvar+p) = 1;
        p = p+1;
        
    elseif powerdata(k,13) ~= 0
        ParentC = find(k == TC(:,2));
        Aeq(TableC(ParentC ,4),totalpfvar+p) = 1;
        p = p+1;
    end
        
end

totalvar = size(Aeq,2);

Dgvar = totalvar- totalpfvar;
%%%formation of objective function

f = zeros(totalvar,1);
f(TableA(1,3)) = 1;
f(TableB(1,3)) = 1;
f(TableC(1,3)) = 1;
f;                           %% objective is min of P from substation

%%
for i=1:totalvar             %% total # of P,Q V
ctype(i)='C';
end


pvmult = 1;

S_1 = 1.15*pvmult*0.04;

Q_limit_1 = sqrt(S_1^2 - PGA(11)^2); 


lb1(1:182,1)= -inf;
lb1(183:274,1)= (0.955^2*ones(92,1));
lb1(275:432,1)= -inf;
lb1(433:512,1)= (0.955^2*ones(80,1));
lb1(513:688,1)= -inf;
lb1(689:777,1)= (0.955^2*ones(89,1));
lb3 = (-Q_limit_1*ones(Dgvar,1));
lb =  [lb1;lb3];

ub1(1:182,1)= inf;
ub1(183:274,1)= (1.1025*ones(92,1));
ub1(275:432,1)= inf;
ub1(433:512,1)= (1.1025*ones(80,1));
ub1(513:688,1)= inf;
ub1(689:777,1)= (1.1025*ones(89,1));
ub3 = (Q_limit_1*ones(Dgvar,1));
ub =  [ub1;ub3];
%%
% optimoptions(@intlinprog,'Display','iter');
% [x_cvr,fval_cvr,exitflag,output] = intlinprog(f,intcon,[],[],Aeq,beq,lb,ub);

 options = cplexoptimset;
    options.Display = 'on';
  [x_cvr,fval_cvr] = cplexmilp(f,[],[],Aeq,beq,[],[],[],lb,ub,ctype);




VA = sqrt(x_cvr(183:274));
VB = sqrt(x_cvr(433:512));
VC = sqrt(x_cvr(689:777));
VA = [sqrt(Vs);VA];
VB = [sqrt(Vs);VB];
VC = [sqrt(Vs);VC];