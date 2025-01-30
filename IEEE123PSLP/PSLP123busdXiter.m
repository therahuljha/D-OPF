clc;
clear all;
%%% here the the deltaX formulation is obtained by using the linear data
%%% and iteration is done on it
function [xnonlin,xlin,fvallin, fvalnonlin, iteration] = PSLP123busdXiter(t)
 

tic
load linedata124C.txt
load branchphA.txt;
load branchphB.txt;
load branchphCC.txt;
lineA= branchphA;
lineB= branchphB;
lineC= branchphCC;
branch=linedata124C;
load powerdata1.txt;
powerdata = powerdata1;
% load powerdata2.txt;
% powerdata = powerdata2;
% load powerdata3.txt;
% powerdata = powerdata3;
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
PLA = powerdata(:,2)/bKVA;
PLB = powerdata(:,4)/bKVA;
PLC = powerdata(:,6)/bKVA;
QLA = powerdata(:,3)/bKVA;
QLB = powerdata(:,5)/bKVA;
QLC = powerdata(:,7)/bKVA;
QCA = powerdata(:,8)/bKVA;
QCB = powerdata(:,9)/bKVA;
QCC = powerdata(:,10)/bKVA;
PGA = powerdata(:,11)/bKVA;
PGB = powerdata(:,12)/bKVA;
PGC = powerdata(:,13)/bKVA;

PLA = PLA.*mult(t);
QLA = QLA.*mult(t);
PLB = PLB.*mult(t);
QLB = QLB.*mult(t);
PLC = PLC.*mult(t);
QLC = QLC.*mult(t);


PGA = PGA.*mult1(t);
PGB = PGB.*mult1(t);
PGC = PGC.*mult1(t);


T=dfsearch(G,1,'edgetonew');
TA=dfsearch(GA,1,'edgetonew');
TB=dfsearch(GB,1,'edgetonew');
TC=dfsearch(GC,1,'edgetonew');

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

 Ap=1:nbusA-1;                          % defining the unknowns for phaseA
 Aq=nbusA:2*(nbusA-1);
 Ai=2*(nbusA)-1:3*(nbusA-1);
 Av=3*(nbusA)-1:4*(nbusA-1)+1;

 Bp=(4*(nbusA-1)+1)+1:(4*(nbusA-1)+1)+(nbusB-1);                    % defining the unknowns for phaseB
 Bq=(4*(nbusA-1)+1)+nbusB :(4*(nbusA-1)+1) + 2*(nbusB-1);
 Bi=(4*(nbusA-1)+1)+ 2*(nbusB)-1:(4*(nbusA-1)+1) + 3*(nbusB-1) ;
 Bv=(4*(nbusA-1)+1)+3*(nbusB)-1:(4*(nbusA-1)+1) + 4*(nbusB-1) +1 ;


Cp=((4*(nbusA-1)+1) + 4*(nbusB-1) +1)+1:((4*(nbusA-1)+1) + 4*(nbusB-1) +1) + (nbusC-1);             % defining the unknowns for phaseC
Cq=((4*(nbusA-1)+1) + 4*(nbusB-1) +1) + nbusC :((4*(nbusA-1)+1) + 4*(nbusB-1) +1) + 2*(nbusC-1);
Ci=((4*(nbusA-1)+1) + 4*(nbusB-1) +1)+ 2*(nbusC)-1:((4*(nbusA-1)+1) + 4*(nbusB-1) +1) + 3*(nbusC-1) ;
Cv=((4*(nbusA-1)+1) + 4*(nbusB-1) +1)+3*(nbusC)-1:((4*(nbusA-1)+1) + 4*(nbusB-1) +1) + 4*(nbusC-1) +1 ;

TableA = [TA(:,1) TA(:,2) Ap'  Aq'  Ai' Av'];
TableB = [TB(:,1) TB(:,2) Bp'  Bq'  Bi' Bv'];
TableC = [TC(:,1) TC(:,2) Cp'  Cq'  Ci' Cv'];

Da = 3*(nbusA)-2:4*(nbusA-1)+1;
Db = (4*(nbusA-1)+1)+3*(nbusB)-2:(4*(nbusA-1)+1) + 4*(nbusB-1) +1;
Dc = ((4*(nbusA-1)+1) + 4*(nbusB-1) +1)+3*(nbusC)-2:((4*(nbusA-1)+1) + 4*(nbusB-1) +1) + 4*(nbusC-1) +1 ;
VolttableA = Da';
VolttableB = Db';
VolttableC = Dc';

nii= zeros(nb,3) ;
p = ((4*(nbusA-1)+1) + 4*(nbusB-1) +1) + 4*(nbusC-1) +1;
for i = 2:nb
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    if (~isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC))
        for j = 1:3
            row = find(i == T(:,2));
           nii (row,j) = p +j;
        end      
        p = nii (row,j);
    elseif (isempty(ParentB))&& (~isempty(ParentA)) && (~isempty(ParentC))        
        row = find(i == T(:,2));
        nii(row,2)=p+1;
        p = nii(row,2);
    elseif (isempty(ParentC))&& (~isempty(ParentA)) && (~isempty(ParentB))
        row = find(i == T(:,2));
        nii(row,1)=p+1;
        p = nii(row,1);
    elseif (isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC))
         row = find(i == T(:,2));
        nii(row,3)=p+1;
        p = nii(row,3);
    end
end
nii;
nii(nb,:) =[]  ;    %% number of variables due to IaIb , IaIc, IbIc
Tableii = [T nii];

 %% calling angle 
 [I_angle] = fbs123;
 I_angle = cell2mat(I_angle);

 %% calling the linear data
[VAlin, VBlin, VClin, xlin, fvallin,TableAlin,TableBlin,TableClin, VolttableAlin,VolttableBlin,VolttableClin] = linearpowerflowwithDG(t);

vs = 1.1025;

%% intitilizing the P,Q,V
xintial(1:182) = xlin(1:182);
xintial(274:365) = xlin(183:274);
xintial(366:523) = xlin(275:432);
xintial(603:682) = xlin(433:512);
xintial(683:858) = xlin(513:688);
xintial(947:1035) = xlin(689:777);

%% current initial value

for i =2:nb
ParentA = find(i == TA(:,2));
 if (~isempty(ParentA))
PocA = find(TA(ParentA,1) == TA(:,1));
xintial(TableA(ParentA,5)) =  (xlin(TableAlin(ParentA,3))^2 +  xlin(TableAlin(ParentA,4))^2)/xlin(VolttableAlin(PocA(1))) ;
 end

ParentB = find(i == TB(:,2));
if (~isempty(ParentB))
PocB = find(TB(ParentB,1) == TB(:,1));
xintial(TableB(ParentB,5)) =  (xlin(TableBlin(ParentB,3))^2 + xlin(TableBlin(ParentB,4))^2)/xlin(VolttableBlin(PocB(1))) ;
end
ParentC = find(i == TC(:,2));
if (~isempty(ParentC))
PocC = find(TC(ParentC,1) == TC(:,1));
xintial(TableC(ParentC,5)) =  (xlin(TableClin(ParentC,3))^2 +  xlin(TableClin(ParentC,4))^2)/xlin(VolttableClin(PocC(1))) ;
end
end

%% IaIb, IaIb intial value

for i =2:nb
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    parent = find(i==T(:,2));
     if (~isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC))      

       xintial(Tableii(parent,3)) = xintial(TableA(ParentA,5))*xintial(TableB(ParentB,5));     
       xintial(Tableii(parent,4)) = xintial(TableA(ParentA,5))*xintial(TableC(ParentC,5));  
       xintial(Tableii(parent,5)) = xintial(TableB(ParentB,5))*xintial(TableC(ParentC,5));
    
      elseif (isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC))
        xintial(Tableii(parent,5)) = xintial(TableB(ParentB,5))*xintial(TableC(ParentC,5));
        
      elseif (isempty(ParentB))&& (~isempty(ParentA)) && (~isempty(ParentC))
        xintial(Tableii(parent,4)) = xintial(TableA(ParentA,5))*xintial(TableC(ParentC,5));    
  
     elseif (isempty(ParentC))&& (~isempty(ParentA)) && (~isempty(ParentB))
        xintial(Tableii(parent,3)) = xintial(TableA(ParentA,5))*xintial(TableB(ParentB,5));     
     
     end
end

% dgvar = 10;   
% xintial(3705+1:3705+dgvar) = xlin(778:end);   %%with 30%DG

% dgvar = 30;   
% xintial(3705+1:3705+dgvar) = xlin(778:end);   %%with 30%DG

xintial(3706:3750) = xlin(778:end);   %%with 50%DG

apred = (xintial(TableA(1,3))+ xintial(TableB(1,3))+xintial(TableC(1,3)));
%% iteration will start from here

iteration=0;
tolerance = 0.001;
error = 1;
delta = 0.04;

while error > tolerance
    
%% intitilization of IaIb constant 
ii0 = [T,zeros(nb-1,3)];
for i =2:nb
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    parent = find(i==T(:,2));
     if (~isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC))
          PocA = find(TA(ParentA,1) == TA(:,1));
          PocB = find(TB(ParentB,1) == TB(:,1));
          PocC = find(TC(ParentC,1) == TC(:,1));
          
         ii0(parent,3)= (xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) + ...
            (xintial(TableA(ParentA,3))^2*xintial(TableBlin(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))+...
            (xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) + ...
            (xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))));
         ii0(parent,4) = (xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) + ...
            (xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+...
            (xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) + ...
            (xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))));
         ii0(parent,5) = (xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,3))^2)/(xintial(VolttableB(PocB(1)))*xintial(VolttableC(PocC(1)))) + ...
            (xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2)/(xintial(VolttableB(PocB(1)))*xintial(VolttableC(PocC(1))))+...
            (xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2)/(xintial(VolttableB(PocB(1)))*xintial(VolttableC(PocC(1)))) + ...
            (xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,4))^2)/(xintial(VolttableB(PocB(1)))*xintial(VolttableC(PocC(1))));        
            
      elseif (isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC))
          PocB = find(TB(ParentB,1) == TB(:,1));
          PocC = find(TC(ParentC,1) == TC(:,1));
         ii0(parent,5) = (xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,3))^2)/(xintial(VolttableB(PocB(1)))*xintial(VolttableC(PocC(1)))) + ...
            (xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2)/(xintial(VolttableB(PocB(1)))*xintial(VolttableC(PocC(1))))+...
            (xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2)/(xintial(VolttableB(PocB(1)))*xintial(VolttableC(PocC(1)))) + ...
            (xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,4))^2)/(xintial(VolttableB(PocB(1)))*xintial(VolttableC(PocC(1))));        
         
      elseif (isempty(ParentB))&& (~isempty(ParentA)) && (~isempty(ParentC))
           PocA = find(TA(ParentA,1) == TA(:,1));
           PocC = find(TC(ParentC,1) == TC(:,1));
        ii0 (parent,4) = (xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) + ...
            (xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+...
            (xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) + ...
            (xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))));       
  
     elseif (isempty(ParentC))&& (~isempty(ParentA)) && (~isempty(ParentB))
           PocA = find(TA(ParentA,1) == TA(:,1));
           PocB = find(TB(ParentB,1) == TB(:,1));
        ii0(parent,3) = (xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) + ...
            (xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))+...
            (xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) + ...
            (xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))));     
     end
end

Tableii_0 =  ii0; 


%% Aeq and beq formation

% CVR_P = 0.6;                %%% CVR factor for P = 0.6
% CVR_Q = 3;                  %%% CVR factor for Q = 3

CVR_P = 2.0;                %%% CVR factor for P = 0.6
CVR_Q = 2.0; 

Aeq = zeros(1235,1235);
beq = zeros(1235,1);
 for i =2:nb
    k = zA(i) ;
    row = find(i == TA(:,1));
 
  if isempty(row)
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    Parent = find(i == T(:,2));
 if isempty(ParentA)  
    Aeq = Aeq;

 elseif ((isempty(ParentB)) && (isempty(ParentC)) ) 
    PocA = find(TA(ParentA,1) == TA(:,1));
    
   %%% for P
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));     
           
    %%% for Q
   Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
   Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2));  
   
     %%% for I
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= -2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= -2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= ((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2);   
   
   %%% for V   
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+3*(nbusA-1),VolttableA(PocA(1)))= -1 ;
   Aeq(ParentA+3*(nbusA-1),(TableA(ParentA,3)))= 2*RAA(TA((ParentA),1),TA((ParentA),2));
   Aeq(ParentA+3*(nbusA-1),(TableA(ParentA,4)))= 2*XAA(TA((ParentA),1),TA((ParentA),2));
   Aeq(ParentA+3*(nbusA-1),(TableA(ParentA,5)))= -(RAA(TA((ParentA),1),TA((ParentA),2))^2+ XAA(TA((ParentA),1),TA((ParentA),2))^2);
       
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2))-(xintial(TableA(ParentA,3))+(-(RAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,5))...
       +(-(CVR_P/2)*PLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6)));   
   beq(ParentA+(nbusA-1)) =  (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2))-(xintial(TableA(ParentA,4))+(-(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,5))...
       +(-(CVR_Q/2)*QLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))) ;   
   beq(ParentA+2*(nbusA-1)) = -(xintial(TableA(ParentA,5))+(-2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,3))...
       +(-2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,4))...
       +(((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2))*xintial(VolttableA(PocA(1))));
   beq(ParentA+3*(nbusA-1)) = -(xintial(TableA(ParentA,6))+(-xintial(VolttableA(PocA(1))))+(2*RAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,3))...
       +(2*XAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,4))+(-(RAA(TA((ParentA),1),TA((ParentA),2))^2+ XAA(TA((ParentA),1),TA((ParentA),2))^2))*xintial(TableA(ParentA,5)));
   
   
 elseif isempty(ParentB)   
    PocA = find(TA(ParentA,1) == TA(:,1));
    PocC = find(TC(ParentC,1) == TC(:,1));    

 %%% for P
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(ParentA,5))= -RAA(TA((ParentA),1),TA((ParentA),2));
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));     
    Aeq(ParentA,Tableii(Parent,4))= -(RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));

  %%% for Q
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,4))= - (XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));  
    
    %%% for I
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= -2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= -2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= ((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2);    
    
  %%% for V   
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+3*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2);
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,5))= -((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2);
   Aeq(ParentA+3*(nbusA-1),Tableii(Parent,4))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
    
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2))-(xintial(TableA(ParentA,3))+(-RAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,5))...
       +(-(CVR_P/2)*PLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
       +(-(RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,4)));  
   beq(ParentA+(nbusA-1)) =  (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2))-(xintial(TableA(ParentA,4))+(-(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,5))...
       +(-(CVR_Q/2)*QLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
       +(-(XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,4)));
    beq(ParentA+2*(nbusA-1)) = -(xintial(TableA(ParentA,5))+(-2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,3))...
       +(-2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,4))...
       +(((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2))*xintial(VolttableA(PocA(1))));
   beq(ParentA+3*(nbusA-1)) = -(xintial(TableA(ParentA,6))+(-xintial(VolttableA(PocA(1))))+(2*RAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,3))...
       +(2*XAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,4))+(-(RAA(TA((ParentA),1),TA((ParentA),2))^2+ XAA(TA((ParentA),1),TA((ParentA),2))^2))*xintial(TableA(ParentA,5))...
       +(-RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2))))*xintial(TableC(ParentC,3))+(-XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2))))*xintial(TableC(ParentC,4))...
       +(-((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableC(ParentC,5))...
       +(- 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,4)));
   
elseif isempty(ParentC) 
    
  PocA = find(TA(ParentA,1) == TA(:,1));
  PocB = find(TB(ParentB,1) == TB(:,1));
   
 %%% for P
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(ParentA,5))= -RAA(TA((ParentA),1),TA((ParentA),2));
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));     
    Aeq(ParentA,Tableii(Parent,3))= -(RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));

  %%% for Q
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,3))= - (XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
    
    %%% for I
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= -2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= -2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= ((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2);    
    
  %%% for V   
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+3*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2);
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,5))= -((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2);
   Aeq(ParentA+3*(nbusA-1),Tableii(Parent,3))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));

   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2))-(xintial(TableA(ParentA,3))+(-RAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,5))...
       +(-(CVR_P/2)*PLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
       +(-(RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,3)));    
   beq(ParentA+(nbusA-1)) =  (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2))-(xintial(TableA(ParentA,4))+(-(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,5))...
       +(-(CVR_Q/2)*QLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
       +(-(XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,3)));
   beq(ParentA+2*(nbusA-1)) = -(xintial(TableA(ParentA,5))+(-2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,3))...
       +(-2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,4))...
       +(((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2))*xintial(VolttableA(PocA(1))));
   beq(ParentA+3*(nbusA-1)) = -(xintial(TableA(ParentA,6))+(-xintial(VolttableA(PocA(1))))+(2*RAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,3))...
       +(2*XAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,4))+(-(RAA(TA((ParentA),1),TA((ParentA),2))^2+ XAA(TA((ParentA),1),TA((ParentA),2))^2))*xintial(TableA(ParentA,5))...
       +(-RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2))))*xintial(TableB(ParentB,3))+(-XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2))))*xintial(TableB(ParentB,4))...
       +(-((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableB(ParentB,5))...
       +(- 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent))))))*xintial(Tableii(Parent,3)));
   
  else
    
  PocA = find(TA(ParentA,1) == TA(:,1));
  PocB = find(TB(ParentB,1) == TB(:,1));
  PocC = find(TC(ParentC,1) == TC(:,1));

  %%% for P
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(ParentA,5))= -RAA(TA((ParentA),1),TA((ParentA),2));
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));     
    Aeq(ParentA,Tableii(Parent,3))= -(RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
    Aeq(ParentA,Tableii(Parent,4))= -(RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));

  %%% for Q
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),TableA(ParentA,6))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,3))= - (XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,4))= - (XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));  
    
    %%% for I
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= -2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= -2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= ((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2);    
    
  %%% for V   
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+3*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2);
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,5))= -((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2);
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,5))= -((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2);
   Aeq(ParentA+3*(nbusA-1),Tableii(Parent,3))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(ParentA+3*(nbusA-1),Tableii(Parent,4))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(ParentA+3*(nbusA-1),Tableii(Parent,5))= - 2*((((RAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));  
    
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2))-(xintial(TableA(ParentA,3))+(-RAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,5))...
       +(-(CVR_P/2)*PLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
       +(-(RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,3))...
       +(-(RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,4)));     
   beq(ParentA+(nbusA-1)) =  (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2))-(xintial(TableA(ParentA,4))+(-(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,5))...
       +(-(CVR_Q/2)*QLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
       +(-(XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,3))...
       +(- (XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,4)));
   beq(ParentA+2*(nbusA-1)) = -(xintial(TableA(ParentA,5))+(-2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,3))...
       +(-2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,4))...
       +(((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2))*xintial(VolttableA(PocA(1))));
   beq(ParentA+3*(nbusA-1)) = -(xintial(TableA(ParentA,6))+(-xintial(VolttableA(PocA(1))))+(2*RAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,3))...
       +(2*XAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,4))+(-(RAA(TA((ParentA),1),TA((ParentA),2))^2+ XAA(TA((ParentA),1),TA((ParentA),2))^2))*xintial(TableA(ParentA,5))...
       +(-RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2))))*xintial(TableB(ParentB,3))+(-XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2))))*xintial(TableB(ParentB,4))...
       +(-((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableB(ParentB,5))...
       +(-RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2))))*xintial(TableC(ParentC,3))+ (-XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2))))*xintial(TableC(ParentC,4)) ...
       +( -((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableC(ParentC,5)) ...
       +(- 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent))))))*xintial(Tableii(Parent,3))...
       +(- 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,4))...
       +(- 2*((((RAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,5))); 
  
 end
 
    else  
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    Parent = find(i == T(:,2));
    
 if isempty(ParentA)  
        Aeq = Aeq;

 elseif ((isempty(ParentB)) && (isempty(ParentC)) ) 
       PocA = find(TA(ParentA,1) == TA(:,1));
     
    %%% for P
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
        end
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
   
    
    %%% for Q
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -1;
    end
    Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
  
           %%% for I
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= -2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= -2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= ((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2);    
   
   %%% for V   
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+3*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
  
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2))-(xintial(TableA(ParentA,3))+(-RAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,5))+ (-(CVR_P/2)*PLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
       +(-xintial(TableA(row(1),3))));           
        for j = 1:length(row)-1
            beq(ParentA) = beq(ParentA)+ xintial(TableA(row(j+1),3));
        end
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2))-(xintial(TableA(ParentA,4))+(-(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,5))+(-(CVR_Q/2)*QLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
       +(- xintial(TableA(row(1),4))));      
         for j = 1:length(row)-1
             beq(ParentA+(nbusA-1)) =  beq(ParentA+(nbusA-1))+xintial(TableA(row(j+1),4));
         end
   beq(ParentA+2*(nbusA-1)) = -(xintial(TableA(ParentA,5))+(-2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,3))+(((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2))*xintial(VolttableA(PocA(1)))...
       +(-2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,4)));
   beq(ParentA+3*(nbusA-1)) = -(xintial(TableA(ParentA,6))-xintial(VolttableA(PocA(1))) + (2*(RAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,3)) + (2*(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,4))...
       +(-((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableA(ParentA,5))) ;
                        
 elseif isempty(ParentB)   
    PocA = find(TA(ParentA,1) == TA(:,1));
    PocC = find(TC(ParentC,1) == TC(:,1));
    
    %%% for P
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
        end
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA,Tableii(Parent,4))= - (RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));       
    
    %%% for Q
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -1;
    end
    Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
     Aeq(ParentA+(nbusA-1),Tableii(Parent,4))= - (XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));
  
           %%% for I
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= -2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= -2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= ((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2);    
   
   %%% for V
   
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+3*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));  
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,5))= -((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+3*(nbusA-1),Tableii(Parent,4))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));

  
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2))-(xintial(TableA(ParentA,3))+(-RAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,5))+ (-(CVR_P/2)*PLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
      +(- (RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,4)) +(-xintial(TableA(row(1),3))));  
         for j = 1:length(row)-1
             beq(ParentA) =  beq(ParentA)+xintial(TableA(row(j+1),3));
         end
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2))-(xintial(TableA(ParentA,4))+(-(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,5))+(-(CVR_Q/2)*QLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
       +(- (XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,4))+(- xintial(TableA(row(1),4))));         
        for j = 1:length(row)-1
             beq(ParentA+(nbusA-1)) =  beq(ParentA+(nbusA-1))+xintial(TableA(row(j+1),4));
        end
   beq(ParentA+2*(nbusA-1)) = -(xintial(TableA(ParentA,5))+(-2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,3))+(((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2))*xintial(VolttableA(PocA(1)))...
       +(-2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,4)));
   beq(ParentA+3*(nbusA-1)) = -(xintial(TableA(ParentA,6))-xintial(VolttableA(PocA(1))) + (2*(RAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,3)) + (2*(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,4))...
       +(-((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableA(ParentA,5))+ (-RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2))))*xintial(TableC(ParentC,3))...
       +(-XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2))))*xintial(TableC(ParentC,4))+(-((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableC(ParentC,5))...
       +(- 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,4))) ;     
        
        
elseif isempty(ParentC) 
  PocA = find(TA(ParentA,1) == TA(:,1));
  PocB = find(TB(ParentB,1) == TB(:,1));
  
     %%% for P
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
        end
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA,Tableii(Parent,3))= - (RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
    
    
    %%% for Q
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -1;
    end
    Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,3))= - (XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
   
  
           %%% for I
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= -2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= -2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= ((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2);    
   
   %%% for V   
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+3*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,5))= -((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+3*(nbusA-1),Tableii(Parent,3))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));

  
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2))-(xintial(TableA(ParentA,3))+(-RAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,5))+ (-(CVR_P/2)*PLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
      +(- (RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,3)) +(-xintial(TableA(row(1),3))));  
         for j = 1:length(row)-1
             beq(ParentA) =  beq(ParentA)+xintial(TableA(row(j+1),3));
         end
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2))-(xintial(TableA(ParentA,4))+(-(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,5))+(-(CVR_Q/2)*QLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
       +(-(XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,3))+(- xintial(TableA(row(1),4))));         
        for j = 1:length(row)-1
             beq(ParentA+(nbusA-1)) =  beq(ParentA+(nbusA-1))+xintial(TableA(row(j+1),4));
        end
   beq(ParentA+2*(nbusA-1)) = -(xintial(TableA(ParentA,5))+(-2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,3))+(((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2))*xintial(VolttableA(PocA(1)))...
       +(-2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,4)));
   beq(ParentA+3*(nbusA-1)) = -(xintial(TableA(ParentA,6))-xintial(VolttableA(PocA(1))) + (2*(RAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,3)) + (2*(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,4))...
       +(-((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableA(ParentA,5))+ ( -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2))))*xintial(TableB(ParentB,3))...
       +(-XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2))))*xintial(TableB(ParentB,4))+(-((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableB(ParentB,5))...
       +( - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent))))))*xintial(Tableii(Parent,3))) ;         
     
else
    
  PocA = find(TA(ParentA,1) == TA(:,1));
  PocB = find(TB(ParentB,1) == TB(:,1));
  PocC = find(TC(ParentC,1) == TC(:,1));
  
    %%% for P
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
        end
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA,Tableii(Parent,3))= - (RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
    Aeq(ParentA,Tableii(Parent,4))= - (RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));       
    
    %%% for Q
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -1;
    end
    Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,3))= - (XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,4))= - (XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));
  
           %%% for I
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= -2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= -2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1 ;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= ((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2);    
   
   %%% for V
   
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+3*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+3*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+3*(nbusA-1),TableB(ParentB,5))= -((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+3*(nbusA-1),TableC(ParentC,5))= -((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+3*(nbusA-1),Tableii(Parent,3))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(ParentA+3*(nbusA-1),Tableii(Parent,4))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(ParentA+3*(nbusA-1),Tableii(Parent,5))= - 2*((((RAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2))-(xintial(TableA(ParentA,3))+(-RAA(TA((ParentA),1),TA((ParentA),2)))*xintial(TableA(ParentA,5))+ (-(CVR_P/2)*PLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
      +(- (RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,3)) ...
      +(- (RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,4))+(-xintial(TableA(row(1),3))));  
         for j = 1:length(row)-1
             beq(ParentA) =  beq(ParentA)+xintial(TableA(row(j+1),3));
         end
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2))-(xintial(TableA(ParentA,4))+(-(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,5))+(-(CVR_Q/2)*QLA(TableA(ParentA,2)))*xintial(TableA(ParentA,6))...
       +(-(XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,3))...
       +(- (XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent)))))*xintial(Tableii(Parent,4))+(- xintial(TableA(row(1),4))));         
        for j = 1:length(row)-1
             beq(ParentA+(nbusA-1)) =  beq(ParentA+(nbusA-1))+xintial(TableA(row(j+1),4));
        end
   beq(ParentA+2*(nbusA-1)) = -(xintial(TableA(ParentA,5))+(-2*(xintial(TableA(ParentA,3)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,3))+(((xintial(TableA(ParentA,3)))^2*xintial(VolttableA(PocA(1)))^-2 + (xintial(TableA(ParentA,4)))^2*xintial(VolttableA(PocA(1)))^-2))*xintial(VolttableA(PocA(1)))...
       +(-2*(xintial(TableA(ParentA,4)))*xintial(VolttableA(PocA(1)))^-1)*xintial(TableA(ParentA,4)));
   beq(ParentA+3*(nbusA-1)) = -(xintial(TableA(ParentA,6))-xintial(VolttableA(PocA(1))) + (2*(RAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,3)) + (2*(XAA(TA((ParentA),1),TA((ParentA),2))))*xintial(TableA(ParentA,4))...
       +(-((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableA(ParentA,5))+ ( -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2))))*xintial(TableB(ParentB,3))...
       +(-XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2))))*xintial(TableB(ParentB,4))+(-((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableB(ParentB,5))...
       +(-RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2))))*xintial(TableC(ParentC,3)) + (-XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2))))*xintial(TableC(ParentC,4))...
       +(-((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2))*xintial(TableC(ParentC,4)) ...
       +( - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent))))))*xintial(Tableii(Parent,3))...
       +( 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,4))...
       +(- 2*((((RAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,5))) ;            
%    end
        end
    end
 Aeq(4*(nbusA-1)+1,VolttableA(1)) = 1;
 beq(4*(nbusA-1)+1) = vs-xintial(VolttableA(1));
 
end
Aeq;

%% Phase B
 for i =2:nb
    k = zB(i) ;
    row = find(i == TB(:,1));
 
  if isempty(row)
      Parent = find(i == T(:,2));
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
 if isempty(ParentB)  
    Aeq = Aeq;
 elseif ((isempty(ParentA)) && (isempty(ParentC)) )  
     
     %%%P     
   PocB = find(TB(ParentB,1) == TB(:,1));
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,3)) = 1;
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));   
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
   
     %%%Q  
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    
      %%% for I
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= -2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= -2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= ((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2);   
    
   %%%V  
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   
   
   beq(4*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2))-(xintial(TableB(ParentB,3))+(-(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_P/2)*PLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6)));
   beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2))-(xintial(TableB(ParentB,4))+(-(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_Q/2)*QLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6)));
   beq(4*(nbusA-1)+1+ParentB+2*(nbusB-1)) = -(xintial(TableB(ParentB,5))+(-2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,3))+(((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1)))...
       +(-2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,4)));
   beq(4*(nbusA-1)+1+ParentB+3*(nbusB-1)) =  -(xintial(TableB(ParentB,6))-xintial(VolttableB(PocB(1))) + (2*(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,3)) + (2*(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,4))...
       +(-((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) )*xintial(TableB(ParentB,5)));
   
   
 elseif isempty(ParentA)   
    PocB = find(TB(ParentB,1) == TB(:,1));
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,3)) = 1;
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));   
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB,Tableii(Parent,5))= -(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,5))= - (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
      %%% for I
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= -2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= -2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= ((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2);   
  
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,5))= -((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),Tableii(Parent,5))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
   beq(4*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2))-(xintial(TableB(ParentB,3))+(-(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_P/2)*PLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +(-(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,5)));
   beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2))-(xintial(TableB(ParentB,4))+(-(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_Q/2)*QLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +(- (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,5)));
   beq(4*(nbusA-1)+1+ParentB+2*(nbusB-1)) = -(xintial(TableB(ParentB,5))+(-2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,3))+(((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1)))...
       +(-2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,4)));
   beq(4*(nbusA-1)+1+ParentB+3*(nbusB-1)) =  -(xintial(TableB(ParentB,6))-xintial(VolttableB(PocB(1))) + (2*(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,3)) + (2*(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,4))...
       +(-((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) )*xintial(TableB(ParentB,5))+ (-RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2))))*xintial(TableC(ParentC,3))...
       +(-XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2))))*xintial(TableC(ParentC,4))+(-((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2))*xintial(TableC(ParentC,5))...
       +(- 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,5)));
   
elseif isempty(ParentC) 
   PocB = find(TB(ParentB,1) == TB(:,1));
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,3)) = 1;
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));   
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB,Tableii(Parent,3))= - (RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
   
    
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,3))= - (XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
   
     %%% for I
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= -2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= -2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= ((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2);   
    
    
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,5))= -((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),Tableii(Parent,3))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   
   beq(4*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2))-(xintial(TableB(ParentB,3))+(-(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_P/2)*PLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +(- (RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,3)));
   beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2))-(xintial(TableB(ParentB,4))+(-(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_Q/2)*QLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +(- (XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,3)));
   beq(4*(nbusA-1)+1+ParentB+2*(nbusB-1)) = -(xintial(TableB(ParentB,5))+(-2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,3))+(((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1)))...
       +(-2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,4)));
   beq(4*(nbusA-1)+1+ParentB+3*(nbusB-1)) =  -(xintial(TableB(ParentB,6))-xintial(VolttableB(PocB(1))) + (2*(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,3)) + (2*(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,4))...
       +(-((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) )*xintial(TableB(ParentB,5))+ (-RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableA(ParentA,3))...
       +(-XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableA(ParentA,4))+(-((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2))*xintial(TableA(ParentA,5))...
       +(- 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent))))))*xintial(Tableii(Parent,3)));
        
  
 else  
     
   PocB = find(TB(ParentB,1) == TB(:,1));
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,3)) = 1;
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));   
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB,Tableii(Parent,3))= - (RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
    Aeq(4*(nbusA-1)+1+ParentB,Tableii(Parent,5))= -(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,3))= - (XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,5))= - (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
      %%% for I
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= -2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= -2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= ((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2);   
  
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,5))= -((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,5))= -((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),Tableii(Parent,3))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),Tableii(Parent,4))= - 2*((((RAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),Tableii(Parent,5))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
  beq(4*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2))-(xintial(TableB(ParentB,3))+(-(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_P/2)*PLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +(- (RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,3))...
       +(-(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,5)));
   beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2))-(xintial(TableB(ParentB,4))+(-(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_Q/2)*QLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +(- (XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,3))...
       +(- (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,5)));
   beq(4*(nbusA-1)+1+ParentB+2*(nbusB-1)) = -(xintial(TableB(ParentB,5))+(-2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,3))+(((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1)))...
       +(-2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,4)));
   beq(4*(nbusA-1)+1+ParentB+3*(nbusB-1)) =  -(xintial(TableB(ParentB,6))-xintial(VolttableB(PocB(1))) + (2*(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,3)) + (2*(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,4))...
       +(-((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) )*xintial(TableB(ParentB,5))+ (-RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableA(ParentA,3))...
       +(-XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableA(ParentA,4))+(-((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2))*xintial(TableA(ParentA,5))...
       +(-RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2))))*xintial(TableC(ParentC,3))+(-XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2))))*xintial(TableC(ParentC,4)) ...
       +(-((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2))*xintial(TableC(ParentC,5))...
       +(- 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent))))))*xintial(Tableii(Parent,3))...
       +(- 2*((((RAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,4))...
  +(- 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,5)) );
   
 end
 
  else  
    Parent = find(i == T(:,2));
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    
 if isempty(ParentB)  
    Aeq = Aeq;

 elseif ((isempty(ParentA)) && (isempty(ParentC)))  
    PocB = find(TB(ParentB,1) == TB(:,1));
%%% P equations  
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(4*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
        end
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));    
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    
    
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -1;
    end
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    
      %%% for I
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= -2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= -2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= ((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2);   
   
  %%% V equations  
    
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;

     
   beq(4*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2))-(xintial(TableB(ParentB,3))+(-(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_P/2)*PLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))+(-xintial(TableB(row(1),3))));
           for j = 1:length(row)-1
             beq(4*(nbusA-1)+1+ParentB) =  beq(4*(nbusA-1)+1+ParentB)+xintial(TableB(row(j+1),3));
           end
         
   beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2))-(xintial(TableB(ParentB,4))+(-(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_Q/2)*QLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))+((-xintial((nbusB-1)+TableB(row(1),3)))));
            for j = 1:length(row)-1
               beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) =  beq(4*(nbusA-1)+1+ParentB+(nbusB-1))+xintial((nbusB-1)+TableB(row(j+1),3));
            end
    
   beq(4*(nbusA-1)+1+ParentB+2*(nbusB-1)) = -(xintial(TableB(ParentB,5))+(-2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,3))+(((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1)))...
       +(-2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,4)));
   
    beq(4*(nbusA-1)+1+ParentB+3*(nbusB-1)) = -(xintial(TableB(ParentB,6))-xintial(VolttableB(PocB(1))) + (2*(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,3)) + (2*(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,4))...
       +(-((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) )*xintial(TableB(ParentB,5)));
   
 elseif isempty(ParentA)   

    PocB = find(TB(ParentB,1) == TB(:,1));
%%% P equations  
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(4*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
        end
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));    
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB,Tableii(Parent,5))= -(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -1;
    end
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,5))= - (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
      %%% for I
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= -2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= -2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= ((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2);   
   
  %%% V equations  
    
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,5))= -((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),Tableii(Parent,5))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
   beq(4*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2))-(xintial(TableB(ParentB,3))+(-(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_P/2)*PLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +(-(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,5))+ (-xintial(TableB(row(1),3))));
           for j = 1:length(row)-1
                     beq(4*(nbusA-1)+1+ParentB) =  beq(4*(nbusA-1)+1+ParentB)+xintial(TableB(row(j+1),3));
           end
         
   beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2))-(xintial(TableB(ParentB,4))+(-(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_Q/2)*QLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +(- (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,5))+((-xintial((nbusB-1)+TableB(row(1),3)))));
    for j = 1:length(row)-1
       beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) =  beq(4*(nbusA-1)+1+ParentB+(nbusB-1))+xintial((nbusB-1)+TableB(row(j+1),3));
    end
    
   beq(4*(nbusA-1)+1+ParentB+2*(nbusB-1)) = -(xintial(TableB(ParentB,5))+(-2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,3))+(((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1)))...
       +(-2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,4)));
   
    beq(4*(nbusA-1)+1+ParentB+3*(nbusB-1)) = -(xintial(TableB(ParentB,6))-xintial(VolttableB(PocB(1))) + (2*(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,3)) + (2*(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,4))...
       +(-((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) )*xintial(TableB(ParentB,5))...
       +(-RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2))))*xintial(TableC(ParentC,3))+(-XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2))))*xintial(TableC(ParentC,4)) ...
       +(-((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2))*xintial(TableC(ParentC,5))...
       +(- 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,5)));
     
   
elseif isempty(ParentC) 
    
PocB = find(TB(ParentB,1) == TB(:,1));
%%% P equations  
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(4*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
        end
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));    
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB,Tableii(Parent,3))= -(RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
   
    
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -1;
    end
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,3))= -(XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
    
      %%% for I
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= -2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= -2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= ((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2);   
    
  %%% V equations  
  
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,5))= -((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),Tableii(Parent,3))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
  
   beq(4*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2))-(xintial(TableB(ParentB,3))+(-(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_P/2)*PLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +( -(RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,3))+ (-xintial(TableB(row(1),3))));
           for j = 1:length(row)-1
                     beq(4*(nbusA-1)+1+ParentB) =  beq(4*(nbusA-1)+1+ParentB)+xintial(TableB(row(j+1),3));
           end
         
   beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2))-(xintial(TableB(ParentB,4))+(-(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_Q/2)*QLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +(-(XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,3))+((-xintial((nbusB-1)+TableB(row(1),3)))));
    for j = 1:length(row)-1
       beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) =  beq(4*(nbusA-1)+1+ParentB+(nbusB-1))+xintial((nbusB-1)+TableB(row(j+1),3));
    end
    
   beq(4*(nbusA-1)+1+ParentB+2*(nbusB-1)) = -(xintial(TableB(ParentB,5))+(-2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,3))+(((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1)))...
       +(-2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,4)));
   
    beq(4*(nbusA-1)+1+ParentB+3*(nbusB-1)) = -(xintial(TableB(ParentB,6))-xintial(VolttableB(PocB(1))) + (2*(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,3)) + (2*(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,4))...
       +(-((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) )*xintial(TableB(ParentB,5))...
       +(-RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableA(ParentA,3))+(-XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableA(ParentA,4)) ...
       +(-((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2))*xintial(TableA(ParentA,5))...
       +(- 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent))))))*xintial(Tableii(Parent,3)));
   
 
   
else

PocB = find(TB(ParentB,1) == TB(:,1));

%%% P equations  
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(4*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
        end
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));    
    Aeq(4*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB,Tableii(Parent,3))= -(RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
    Aeq(4*(nbusA-1)+1+ParentB,Tableii(Parent,5))= -(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -1;
    end
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,3))= -(XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
    Aeq(4*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,5))= - (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
      %%% for I
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= -2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= -2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= ((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2);   
   
  %%% V equations     
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableA(ParentA,5))= -((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),TableC(ParentC,5))= -((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),Tableii(Parent,3))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),Tableii(Parent,4))= - 2*((((RAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(4*(nbusA-1)+1+ParentB+3*(nbusB-1),Tableii(Parent,5))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
 beq(4*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2))-(xintial(TableB(ParentB,3))+(-(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_P/2)*PLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +( -(RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,3))...
       +(-(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,5))+ (-xintial(TableB(row(1),3))));
           for j = 1:length(row)-1
                     beq(4*(nbusA-1)+1+ParentB) =  beq(4*(nbusA-1)+1+ParentB)+xintial(TableB(row(j+1),3));
           end
         
   beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2))-(xintial(TableB(ParentB,4))+(-(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,5))+(-(CVR_Q/2)*QLB(TableB(ParentB,2)))*xintial(TableB(ParentB,6))...
       +(-(XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,3))...
       +(- (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent)))))*xintial(Tableii(Parent,5))+((-xintial((nbusB-1)+TableB(row(1),3)))));
    for j = 1:length(row)-1
       beq(4*(nbusA-1)+1+ParentB+(nbusB-1)) =  beq(4*(nbusA-1)+1+ParentB+(nbusB-1))+xintial((nbusB-1)+TableB(row(j+1),3));
    end
    
   beq(4*(nbusA-1)+1+ParentB+2*(nbusB-1)) = -(xintial(TableB(ParentB,5))+(-2*(xintial(TableB(ParentB,3)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,3))+(((xintial(TableB(ParentB,3)))^2*xintial(VolttableB(PocB(1)))^-2 + (xintial(TableB(ParentB,4)))^2*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1)))...
       +(-2*(xintial(TableB(ParentB,4)))*xintial(VolttableB(PocB(1)))^-1)*xintial(TableB(ParentB,4)));
   
    beq(4*(nbusA-1)+1+ParentB+3*(nbusB-1)) = -(xintial(TableB(ParentB,6))-xintial(VolttableB(PocB(1))) + (2*(RBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,3)) + (2*(XBB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableB(ParentB,4))...
       +(-((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) )*xintial(TableB(ParentB,5))...
       +(-RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableA(ParentA,3))+(-XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2))))*xintial(TableA(ParentA,4)) ...
       +(-((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2))*xintial(TableA(ParentA,5))...
       +( -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2))))*xintial(TableC(ParentC,3))+ (-XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2))))*xintial(TableC(ParentC,4))...
       +(-((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2))*xintial(TableC(ParentC,5)) ...
       +(- 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent))))))*xintial(Tableii(Parent,3))...
       +(- 2*((((RAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,4))...
       +(- 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,5)));


end
  end
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1,VolttableB(1)) = 1;
    beq(4*(nbusA-1)+1+4*(nbusB-1)+1) = vs-xintial(VolttableB(1));

 
 end
Aeq;

%% Phase C
 for i =2:nb
    k = zC(i) ;
    row = find(i == TC(:,1));
 
  if isempty(row)
    Parent = find(i == T(:,2));
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
 if isempty(ParentC)  
    Aeq = Aeq;
 elseif ((isempty(ParentB)) && (isempty(ParentA)) )  
    PocC = find(TC(ParentC,1) == TC(:,1));
%%% P equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
   
    
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    
          %%% for I
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= -2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= -2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= ((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2);   
   
    
  %%% V equations     
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   

   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2))-(xintial(TableC(ParentC,3))+(-(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_P/2)*PLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2))-(xintial(TableC(ParentC,4))+(-(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_Q/2)*QLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1)) = -(xintial(TableC(ParentC,5))+(-2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,3))+(((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1)))...
       +(-2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,4)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1)) =  -(xintial(TableC(ParentC,6))-xintial(VolttableC(PocC(1))) + (2*(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,3)) + (2*(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,4))...
       +(-((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableC(ParentC,5)));
   
 elseif isempty(ParentA)   
  PocC = find(TC(ParentC,1) == TC(:,1));
     
    %%% P equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,Tableii(Parent,5))= -(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,5))= -(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
              %%% for I
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= -2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= -2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= ((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2);   
    
  %%% V equations     
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,5))= -((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),Tableii(Parent,5))= - 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  

 beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2))-(xintial(TableC(ParentC,3))+(-(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_P/2)*PLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
        +(-(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,5)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2))-(xintial(TableC(ParentC,4))+(-(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_Q/2)*QLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
       +(-(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,5)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1)) = -(xintial(TableC(ParentC,5))+(-2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,3))+(((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1)))...
       +(-2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,4)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1)) =  -(xintial(TableC(ParentC,6))-xintial(VolttableC(PocC(1))) + (2*(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,3)) + (2*(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,4))...
       +(-((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableC(ParentC,5))  ...
       +(-RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableB(ParentB,3))+(-XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableB(ParentB,4)) ...
       +(-((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableB(ParentB,5))...
       +(- 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,5)) );
   
elseif isempty(ParentB) 
   PocC = find(TC(ParentC,1) == TC(:,1));
%%% P equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,Tableii(Parent,4))= -(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
        
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,4))= -(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))) + RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    
              %%% for I
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= -2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= -2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= ((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2);   
   
  %%% V equations     
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,5))= -((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),Tableii(Parent,4))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
  
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2))-(xintial(TableC(ParentC,3))+(-(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_P/2)*PLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
       +(-(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,4)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2))-(xintial(TableC(ParentC,4))+(-(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_Q/2)*QLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
       +(-(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))) + RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,4)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1)) = -(xintial(TableC(ParentC,5))+(-2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,3))+(((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1)))...
       +(-2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,4)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1)) =  -(xintial(TableC(ParentC,6))-xintial(VolttableC(PocC(1))) + (2*(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,3)) + (2*(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,4))...
       +(-((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableC(ParentC,5))+ (-RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableA(ParentA,3))...
       +(-XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableA(ParentA,4))+(-((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableA(ParentA,5))...
       +(- 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,4))); 
    
  
 else
          
  PocC = find(TC(ParentC,1) == TC(:,1));
    
  %%% P equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,Tableii(Parent,4))= -(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,Tableii(Parent,5))= -(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentB,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,4))= -(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))) + RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,5))= -(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
   %%% for I
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= -2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= -2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= ((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2);   
   
   %%% V equations     
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,5))= -((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,5))= -((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),Tableii(Parent,3))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))-(RAB(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),Tableii(Parent,4))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),Tableii(Parent,5))= - 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
 beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2))-(xintial(TableC(ParentC,3))+(-(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_P/2)*PLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
       +(-(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,4))...
       +(-(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,5)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2))-(xintial(TableC(ParentC,4))+(-(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_Q/2)*QLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
       +(-(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))) + RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,4))...
       +(-(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,5)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1)) = -(xintial(TableC(ParentC,5))+(-2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,3))+(((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1)))...
       +(-2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,4)));
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1)) =  -(xintial(TableC(ParentC,6))-xintial(VolttableC(PocC(1))) + (2*(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,3)) + (2*(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,4))...
       +(-((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableC(ParentC,5))+ (-RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableA(ParentA,3))...
       +(-XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableA(ParentA,4))+(-((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableA(ParentA,5))...
       +(-RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableB(ParentB,3))+(-XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableB(ParentB,4)) ...
       +(-((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableB(ParentB,5))...
       +(- 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))-(RAB(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent))))))*xintial(Tableii(Parent,3))...
       +(- 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,4))...
  +(- 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,5)) );

 end
 
  else 
    Parent = find(i == T(:,2));
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    
 if isempty(ParentC)  
    Aeq = Aeq;

 elseif ((isempty(ParentA)) && (isempty(ParentB)))  
    PocC = find(TC(ParentC,1) == TC(:,1));
    %%% P equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
        end
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));    
    
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;    
    for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -1;
    end
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    
   %%% for I
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= -2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= -2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= ((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2);   
   
  %%% V equations     
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   

  beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2))-(xintial(TableC(ParentC,3))+(-(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_P/2)*PLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
    + ((-xintial(TableC(row(1),3)))));
           for j = 1:length(row)-1
               beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC) = beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)+xintial(TableC(row(j+1),3));
           end
           
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2))-(xintial(TableC(ParentC,4))+(-(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_Q/2)*QLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
     + ((-xintial(TableC(row(1),4)))));
           for j = 1:length(row)-1
                beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) + xintial(TableC(row(j+1),4));
           end
   
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1)) = -(xintial(TableC(ParentC,5))+(-2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,3))+(((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1)))...
       +(-2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,4)));
   
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1)) =  -(xintial(TableC(ParentC,6))-xintial(VolttableC(PocC(1))) + (2*(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,3)) + (2*(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,4))...
       +(-((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableC(ParentC,5)));
   
 elseif isempty(ParentA)   

PocC = find(TC(ParentC,1) == TC(:,1));
%%% P equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
        end
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,Tableii(Parent,5))= -(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -1;
    end
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,5))= -(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
              %%% for I
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= -2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= -2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= ((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2);   
   
  %%% V equations     
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,5))= -((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),Tableii(Parent,5))= - 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
  beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2))-(xintial(TableC(ParentC,3))+(-(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_P/2)*PLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
      +(-(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,5)) + ((-xintial(TableC(row(1),3)))));
           for j = 1:length(row)-1
               beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC) = beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)+xintial(TableC(row(j+1),3));
           end
           
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2))-(xintial(TableC(ParentC,4))+(-(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_Q/2)*QLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
      +(-(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,5))+ ((-xintial(TableC(row(1),4)))));
           for j = 1:length(row)-1
                beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) + xintial(TableC(row(j+1),4));
           end
   
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1)) = -(xintial(TableC(ParentC,5))+(-2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,3))+(((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1)))...
       +(-2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,4)));
   
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1)) =  -(xintial(TableC(ParentC,6))-xintial(VolttableC(PocC(1))) + (2*(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,3)) + (2*(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,4))...
       +(-((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableC(ParentC,5)) ...
       +(-RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableB(ParentB,3))+(-XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableB(ParentB,4)) ...
       +(-((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableB(ParentB,5))...
       +(- 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,5)) );
   
elseif isempty(ParentB) 
PocC = find(TC(ParentC,1) == TC(:,1));
%%% P equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
        end
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,Tableii(Parent,4))= -(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    
    
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;    
    for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -1;
    end
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,4))= -(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))+ RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    
              %%% for I
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= -2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= -2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= ((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2);   
    
  %%% V equations     
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,5))= -((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),Tableii(Parent,4))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
  

  beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2))-(xintial(TableC(ParentC,3))+(-(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_P/2)*PLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
       +(-(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,4))+ ((-xintial(TableC(row(1),3)))));
           for j = 1:length(row)-1
               beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC) = beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)+xintial(TableC(row(j+1),3));
           end
           
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2))-(xintial(TableC(ParentC,4))+(-(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_Q/2)*QLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
       +(-(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))) + RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,4))+ ((-xintial(TableC(row(1),4)))));
           for j = 1:length(row)-1
                beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) + xintial(TableC(row(j+1),4));
           end
   
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1)) = -(xintial(TableC(ParentC,5))+(-2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,3))+(((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1)))...
       +(-2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,4)));
   
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1)) =  -(xintial(TableC(ParentC,6))-xintial(VolttableC(PocC(1))) + (2*(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,3)) + (2*(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,4))...
       +(-((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableC(ParentC,5))+ (-RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableA(ParentA,3))...
       +(-XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableA(ParentA,4))+(-((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableA(ParentA,5))...
       +(- 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,4)));   
   
   
else

PocC = find(TC(ParentC,1) == TC(:,1));
    %%% P equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
        end
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,Tableii(Parent,4))= -(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC,Tableii(Parent,5))= -(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
    %%% Q equations  
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1), TableC(row(1),4))= -1;    
    for j = 1:length(row)-1
        Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(row(j+1),4)) =  -1;
    end
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,4))= -(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))) + RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,5))= -(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
   %%% for I
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= -2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= -2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1 ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= ((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2);   
   
   %%% V equations     
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableA(ParentA,5))= -((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),TableB(ParentB,5))= -((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),Tableii(Parent,3))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))-(RAB(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),Tableii(Parent,4))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1),Tableii(Parent,5))= - 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  

  beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2))-(xintial(TableC(ParentC,3))+(-(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_P/2)*PLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
       +(-(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,4))...
       +(-(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,5)) + ((-xintial(TableC(row(1),3)))));
           for j = 1:length(row)-1
               beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC) = beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC)+xintial(TableC(row(j+1),3));
           end
           
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2))-(xintial(TableC(ParentC,4))+(-(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,5))+(-(CVR_Q/2)*QLC(TableC(ParentC,2)))*xintial(TableC(ParentC,6))...
       +(-(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))) + RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,4))...
       +(-(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))))*xintial(Tableii(Parent,5))+ ((-xintial(TableC(row(1),4)))));
           for j = 1:length(row)-1
                beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) = beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+(nbusC-1)) + xintial(TableC(row(j+1),4));
           end
   
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+2*(nbusC-1)) = -(xintial(TableC(ParentC,5))+(-2*(xintial(TableC(ParentC,3)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,3))+(((xintial(TableC(ParentC,3)))^2*xintial(VolttableC(PocC(1)))^-2 + (xintial(TableC(ParentC,4)))^2*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1)))...
       +(-2*(xintial(TableC(ParentC,4)))*xintial(VolttableC(PocC(1)))^-1)*xintial(TableC(ParentC,4)));
   
   beq(4*(nbusA-1)+1+4*(nbusB-1)+1+ParentC+3*(nbusC-1)) =  -(xintial(TableC(ParentC,6))-xintial(VolttableC(PocC(1))) + (2*(RCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,3)) + (2*(XCC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableC(ParentC,4))...
       +(-((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableC(ParentC,5))+ (-RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableA(ParentA,3))...
       +(-XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableA(ParentA,4))+(-((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableA(ParentA,5))...
       +(-RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableB(ParentB,3))+(-XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2))))*xintial(TableB(ParentB,4)) ...
       +(-((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2))*xintial(TableB(ParentB,5))...
       +(- 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))-(RAB(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent))))))*xintial(Tableii(Parent,3))...
       +(- 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,4))...
  +(- 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent))))))*xintial(Tableii(Parent,5)) );

end
  end
   
  Aeq(4*(nbusA-1)+1+4*(nbusB-1)+1+4*(nbusC-1)+1,VolttableC(1)) = 1;
  beq(4*(nbusA-1)+1+4*(nbusB-1)+1+4*(nbusC-1)+1) = vs-xintial(VolttableC(1));
 
 end
Aeq;

%% IaIb, IaIc, IbIc equations

for i = 2:nb
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    Parent = find(i == T(:,2));
    
    if ((isempty(ParentA)) && (isempty(ParentB))) 
        Aeq= Aeq;     
        
    elseif ((isempty(ParentB)) && (isempty(ParentC))) 
        Aeq= Aeq;    
        
     elseif ((isempty(ParentA))&&(isempty(ParentB)) && (isempty(ParentC)))      
         Aeq= Aeq;     
         
    elseif ((isempty(ParentA)) && (isempty(ParentC))) 
        Aeq= Aeq;    
        
    elseif isempty(ParentA)  
         PocB = find(TB(ParentB,1) == TB(:,1));
         PocC = find(TC(ParentC,1) == TC(:,1));
        %%%%IbIc
   
   Aeq(Tableii(Parent,5), Tableii(Parent,5)) = 1;
   Aeq(Tableii(Parent,5), TableB(ParentB,3)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), TableB(ParentB,4)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), VolttableB(PocB(1))) = (0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2));
  
   Aeq(Tableii(Parent,5), TableC(ParentC,3)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
     
   Aeq(Tableii(Parent,5), TableC(ParentC,4)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), VolttableC(PocC(1))) = (0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1));
  
  beq(Tableii(Parent,5)) = -(xintial(Tableii(Parent,5))+(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,3))+ ( -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,4)) ...
      + ((0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2)))*xintial(VolttableB(PocB(1)))...
      +(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))+ (2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableC(ParentC,3))...
      +(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))+ (2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableC(ParentC,4))...
      +((0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1)))*xintial(VolttableC(PocC(1))));
   
    elseif isempty(ParentB)
        
          PocA = find(TA(ParentA,1) == TA(:,1));
          PocC = find(TC(ParentC,1) == TC(:,1));
                
                %%%%IaIc
       Aeq(Tableii(Parent,4), Tableii(Parent,4)) = 1;
       Aeq(Tableii(Parent,4), TableA(ParentA,3)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
  
       Aeq(Tableii(Parent,4), TableA(ParentA,4)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
  
       Aeq(Tableii(Parent,4), VolttableA(PocA(1))) = (0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1);
  
      Aeq(Tableii(Parent,4), TableC(ParentC,3)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
     
      Aeq(Tableii(Parent,4), TableC(ParentC,4)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
  
      Aeq(Tableii(Parent,4), VolttableC(PocC(1))) = (0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2);
  
      beq(Tableii(Parent,4)) = -(xintial(Tableii(Parent,4))+(-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) )))*xintial(TableA(ParentA,3)) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableA(ParentA,4)) ...
      + ((0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1))*xintial(VolttableA(PocA(1))) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableC(ParentC,3)) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableC(ParentC,4)) ...
      + ((0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1))));
  
    elseif isempty(ParentC)
        
          PocA = find(TA(ParentA,1) == TA(:,1));
         PocB = find(TB(ParentB,1) == TB(:,1));
              
         %%%%IaIb
        Aeq(Tableii(Parent,3), Tableii(Parent,3)) = 1;
        Aeq(Tableii(Parent,3), TableA(ParentA,3)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
       Aeq(Tableii(Parent,3), TableA(ParentA,4)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
       Aeq(Tableii(Parent,3), VolttableA(PocA(1))) = (0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1);
  
      Aeq(Tableii(Parent,3), TableB(ParentB,3)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
     
      Aeq(Tableii(Parent,3), TableB(ParentB,4)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
      Aeq(Tableii(Parent,3), VolttableB(PocB(1))) = (0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2);
  
      beq(Tableii(Parent,3)) = -(xintial(Tableii(Parent,3))+ (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableA(ParentA,3))...
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableA(ParentA,4))...
         + ((0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
         + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
         + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1))*xintial(VolttableA(PocA(1)))... 
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,3))...
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,4))...
         + ((0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
         + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
         + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1))));
  
    else
         PocA = find(TA(ParentA,1) == TA(:,1));
         PocB = find(TB(ParentB,1) == TB(:,1));
         PocC = find(TC(ParentC,1) == TC(:,1));
          if  ((Tableii_0(Parent,3) == 0) && (Tableii_0(Parent,4) == 0) && (Tableii_0(Parent,5) == 0))               
                 %%%%IaIb
        Aeq(Tableii(Parent,3), Tableii(Parent,3)) = 1;
        Aeq(Tableii(Parent,3), TableA(ParentA,3)) = 0;  
        Aeq(Tableii(Parent,3), TableA(ParentA,4)) = 0;  
        Aeq(Tableii(Parent,3), VolttableA(PocA(1))) = 0;  
        Aeq(Tableii(Parent,3), TableB(ParentB,3)) = 0;     
        Aeq(Tableii(Parent,3), TableB(ParentB,4)) = 0;  
        Aeq(Tableii(Parent,3), VolttableB(PocB(1))) = 0;
        
        beq(Tableii(Parent,3)) = 0;
  
    %%%%IaIc
       Aeq(Tableii(Parent,4), Tableii(Parent,4)) = 1;
       Aeq(Tableii(Parent,4), TableA(ParentA,3)) = 0;  
       Aeq(Tableii(Parent,4), TableA(ParentA,4)) = 0;  
       Aeq(Tableii(Parent,4), VolttableA(PocA(1))) = 0;  
       Aeq(Tableii(Parent,4), TableC(ParentC,3)) = 0;     
       Aeq(Tableii(Parent,4), TableC(ParentC,4)) = 0;  
       Aeq(Tableii(Parent,4), VolttableC(PocC(1))) = 0;
  
        beq(Tableii(Parent,4)) = 0;
        
   %%%%IbIc
   
       Aeq(Tableii(Parent,5), Tableii(Parent,5)) = 1;
       Aeq(Tableii(Parent,5), TableB(ParentB,3)) = 0;  
       Aeq(Tableii(Parent,5), TableB(ParentB,4)) = 0;  
       Aeq(Tableii(Parent,5), VolttableB(PocB(1))) = 0;  
       Aeq(Tableii(Parent,5), TableC(ParentC,3)) = 0;     
       Aeq(Tableii(Parent,5), TableC(ParentC,4)) = 0;  
       Aeq(Tableii(Parent,5), VolttableC(PocC(1))) = 0; 
        beq(Tableii(Parent,5)) = 0;
        
        
          elseif  ((Tableii_0(Parent,4) == 0) && (Tableii_0(Parent,5) == 0))  
              
               %%%%IaIb
        Aeq(Tableii(Parent,3), Tableii(Parent,3)) = 1;
        Aeq(Tableii(Parent,3), TableA(ParentA,3)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
       Aeq(Tableii(Parent,3), TableA(ParentA,4)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
       Aeq(Tableii(Parent,3), VolttableA(PocA(1))) = (0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1);
  
      Aeq(Tableii(Parent,3), TableB(ParentB,3)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
     
      Aeq(Tableii(Parent,3), TableB(ParentB,4)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
      Aeq(Tableii(Parent,3), VolttableB(PocB(1))) = (0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2);
  
       beq(Tableii(Parent,3)) = -(xintial(Tableii(Parent,3))+ (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableA(ParentA,3))...
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableA(ParentA,4))...
         + ((0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
         + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
         + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1))*xintial(VolttableA(PocA(1)))... 
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,3))...
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,4))...
         + ((0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
         + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
         + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1))));
  
    %%%%IaIc
       Aeq(Tableii(Parent,4), Tableii(Parent,4)) = 1;
       Aeq(Tableii(Parent,4), TableA(ParentA,3)) = 0;  
       Aeq(Tableii(Parent,4), TableA(ParentA,4)) = 0;  
       Aeq(Tableii(Parent,4), VolttableA(PocA(1))) = 0;  
      Aeq(Tableii(Parent,4), TableC(ParentC,3)) = 0;     
      Aeq(Tableii(Parent,4), TableC(ParentC,4)) = 0;  
      Aeq(Tableii(Parent,4), VolttableC(PocC(1))) = 0;  
      beq(Tableii(Parent,4)) = 0;
  
   %%%%IbIc
   
   Aeq(Tableii(Parent,5), Tableii(Parent,5)) = 1;
   Aeq(Tableii(Parent,5), TableB(ParentB,3)) = 0;  
   Aeq(Tableii(Parent,5), TableB(ParentB,4)) = 0;  
   Aeq(Tableii(Parent,5), VolttableB(PocB(1))) = 0;  
   Aeq(Tableii(Parent,5), TableC(ParentC,3)) = 0;     
   Aeq(Tableii(Parent,5), TableC(ParentC,4)) = 0;  
   Aeq(Tableii(Parent,5), VolttableC(PocC(1))) = 0;
   beq(Tableii(Parent,5)) = 0;
        
          elseif ((Tableii_0(Parent,3) == 0)  && (Tableii_0(Parent,5) == 0))   
              
              %%%%IaIb
              Aeq(Tableii(Parent,3), Tableii(Parent,3)) = 1;
              Aeq(Tableii(Parent,3), TableA(ParentA,3)) = 0;
              Aeq(Tableii(Parent,3), TableA(ParentA,4)) = 0;
              Aeq(Tableii(Parent,3), VolttableA(PocA(1))) = 0;
              Aeq(Tableii(Parent,3), TableB(ParentB,3)) = 0;
              Aeq(Tableii(Parent,3), TableB(ParentB,4)) = 0;
              Aeq(Tableii(Parent,3), VolttableB(PocB(1))) = 0;
              beq(Tableii(Parent,3)) = 0;
  
    %%%%IaIc
       Aeq(Tableii(Parent,4), Tableii(Parent,4)) = 1;
       Aeq(Tableii(Parent,4), TableA(ParentA,3)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ));
  
       Aeq(Tableii(Parent,4), TableA(ParentA,4)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
  
       Aeq(Tableii(Parent,4), VolttableA(PocA(1))) = (0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1);
  
      Aeq(Tableii(Parent,4), TableC(ParentC,3)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
     
      Aeq(Tableii(Parent,4), TableC(ParentC,4)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
  
      Aeq(Tableii(Parent,4), VolttableC(PocC(1))) = (0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2);
  
      beq(Tableii(Parent,4)) = -(xintial(Tableii(Parent,4))+(-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) )))*xintial(TableA(ParentA,3)) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableA(ParentA,4)) ...
      + ((0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1))*xintial(VolttableA(PocA(1))) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableC(ParentC,3)) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableC(ParentC,4)) ...
      + ((0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1))));
  
              %%%%IbIc

              Aeq(Tableii(Parent,5), Tableii(Parent,5)) = 1;
              Aeq(Tableii(Parent,5), TableB(ParentB,3)) = 0;
              Aeq(Tableii(Parent,5), TableB(ParentB,4)) = 0;
              Aeq(Tableii(Parent,5), VolttableB(PocB(1))) = 0;
              Aeq(Tableii(Parent,5), TableC(ParentC,3)) = 0;
              Aeq(Tableii(Parent,5), TableC(ParentC,4)) = 0;
              Aeq(Tableii(Parent,5), VolttableC(PocC(1))) = 0;
              beq(Tableii(Parent,5)) = 0;
   
          elseif ((Tableii_0(Parent,3) == 0)  && (Tableii_0(Parent,4) == 0))   
              
               %%%%IaIb
        Aeq(Tableii(Parent,3), Tableii(Parent,3)) = 1;
        Aeq(Tableii(Parent,3), TableA(ParentA,3)) = 0;  
        Aeq(Tableii(Parent,3), TableA(ParentA,4)) = 0;  
        Aeq(Tableii(Parent,3), VolttableA(PocA(1))) = 0; 
        Aeq(Tableii(Parent,3), TableB(ParentB,3)) = 0;     
        Aeq(Tableii(Parent,3), TableB(ParentB,4)) = 0;  
        Aeq(Tableii(Parent,3), VolttableB(PocB(1))) = 0;  
        beq(Tableii(Parent,3)) = 0;
  
    %%%%IaIc
       Aeq(Tableii(Parent,4), Tableii(Parent,4)) = 1;
       Aeq(Tableii(Parent,4), TableA(ParentA,3)) = 0;  
       Aeq(Tableii(Parent,4), TableA(ParentA,4)) = 0;  
       Aeq(Tableii(Parent,4), VolttableA(PocA(1))) = 0;  
       Aeq(Tableii(Parent,4), TableC(ParentC,3)) = 0;     
       Aeq(Tableii(Parent,4), TableC(ParentC,4)) = 0;  
       Aeq(Tableii(Parent,4), VolttableC(PocC(1))) = 0;  
       beq(Tableii(Parent,4)) = 0;
  
   %%%%IbIc
   
   Aeq(Tableii(Parent,5), Tableii(Parent,5)) = 1;
   Aeq(Tableii(Parent,5), TableB(ParentB,3)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), TableB(ParentB,4)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), VolttableB(PocB(1))) = (0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2));
  
   Aeq(Tableii(Parent,5), TableC(ParentC,3)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
     
   Aeq(Tableii(Parent,5), TableC(ParentC,4)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), VolttableC(PocC(1))) = (0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1));
  
  
   beq(Tableii(Parent,5)) = -(xintial(Tableii(Parent,5))+(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,3))+ ( -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,4)) ...
      + ((0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2)))*xintial(VolttableB(PocB(1)))...
      +(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))+ (2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableC(ParentC,3))...
      +(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))+ (2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableC(ParentC,4))...
      +((0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1)))*xintial(VolttableC(PocC(1))));
  
  
   
          elseif ((Tableii_0(Parent,5) == 0))   
              
              %%%%IaIb
        Aeq(Tableii(Parent,3), Tableii(Parent,3)) = 1;
        Aeq(Tableii(Parent,3), TableA(ParentA,3)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
       Aeq(Tableii(Parent,3), TableA(ParentA,4)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
       Aeq(Tableii(Parent,3), VolttableA(PocA(1))) = (0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1);
  
      Aeq(Tableii(Parent,3), TableB(ParentB,3)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
     
      Aeq(Tableii(Parent,3), TableB(ParentB,4)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
      Aeq(Tableii(Parent,3), VolttableB(PocB(1))) = (0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2);
  
       beq(Tableii(Parent,3)) = -(xintial(Tableii(Parent,3))+ (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableA(ParentA,3))...
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableA(ParentA,4))...
         + ((0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
         + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
         + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1))*xintial(VolttableA(PocA(1)))... 
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,3))...
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,4))...
         + ((0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
         + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
         + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1))));
  
    %%%%IaIc
       Aeq(Tableii(Parent,4), Tableii(Parent,4)) = 1;
       Aeq(Tableii(Parent,4), TableA(ParentA,3)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ));
  
       Aeq(Tableii(Parent,4), TableA(ParentA,4)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
  
       Aeq(Tableii(Parent,4), VolttableA(PocA(1))) = (0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1);
  
      Aeq(Tableii(Parent,4), TableC(ParentC,3)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
     
      Aeq(Tableii(Parent,4), TableC(ParentC,4)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
  
      Aeq(Tableii(Parent,4), VolttableC(PocC(1))) = (0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2);
  
      beq(Tableii(Parent,4)) = -(xintial(Tableii(Parent,4))+(-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) )))*xintial(TableA(ParentA,3)) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableA(ParentA,4)) ...
      + ((0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1))*xintial(VolttableA(PocA(1))) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableC(ParentC,3)) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableC(ParentC,4)) ...
      + ((0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1))));
  
   %%%%IbIc
   
   Aeq(Tableii(Parent,5), Tableii(Parent,5)) = 1;
   Aeq(Tableii(Parent,5), TableB(ParentB,3)) = 0;  
   Aeq(Tableii(Parent,5), TableB(ParentB,4)) = 0;  
   Aeq(Tableii(Parent,5), VolttableB(PocB(1))) = 0;  
   Aeq(Tableii(Parent,5), TableC(ParentC,3)) = 0;     
   Aeq(Tableii(Parent,5), TableC(ParentC,4)) = 0;  
   Aeq(Tableii(Parent,5), VolttableC(PocC(1))) = 0;   
   beq(Tableii(Parent,5)) = 0;
   
  
          elseif ((Tableii_0(Parent,4) == 0))   
              
              %%%%IaIb
        Aeq(Tableii(Parent,3), Tableii(Parent,3)) = 1;
        Aeq(Tableii(Parent,3), TableA(ParentA,3)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
       Aeq(Tableii(Parent,3), TableA(ParentA,4)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
       Aeq(Tableii(Parent,3), VolttableA(PocA(1))) = (0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1);
  
      Aeq(Tableii(Parent,3), TableB(ParentB,3)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
     
      Aeq(Tableii(Parent,3), TableB(ParentB,4)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
      Aeq(Tableii(Parent,3), VolttableB(PocB(1))) = (0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2);
  
       beq(Tableii(Parent,3)) = -(xintial(Tableii(Parent,3))+ (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableA(ParentA,3))...
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableA(ParentA,4))...
         + ((0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
         + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
         + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1))*xintial(VolttableA(PocA(1)))... 
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,3))...
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,4))...
         + ((0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
         + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
         + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1))));
  
    %%%%IaIc
       Aeq(Tableii(Parent,4), Tableii(Parent,4)) = 1;
       Aeq(Tableii(Parent,4), TableA(ParentA,3)) = 0;  
       Aeq(Tableii(Parent,4), TableA(ParentA,4)) = 0;  
       Aeq(Tableii(Parent,4), VolttableA(PocA(1))) = 0;  
       Aeq(Tableii(Parent,4), TableC(ParentC,3)) = 0;     
       Aeq(Tableii(Parent,4), TableC(ParentC,4)) =0;  
       Aeq(Tableii(Parent,4), VolttableC(PocC(1))) = 0;  
       beq(Tableii(Parent,4)) = 0;
  
   %%%%IbIc
   
   Aeq(Tableii(Parent,5), Tableii(Parent,5)) = 1;
   Aeq(Tableii(Parent,5), TableB(ParentB,3)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), TableB(ParentB,4)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), VolttableB(PocB(1))) = (0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2));
  
   Aeq(Tableii(Parent,5), TableC(ParentC,3)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
     
   Aeq(Tableii(Parent,5), TableC(ParentC,4)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), VolttableC(PocC(1))) = (0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1));
  
  
   beq(Tableii(Parent,5)) = -(xintial(Tableii(Parent,5))+(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,3))+ ( -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,4)) ...
      + ((0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2)))*xintial(VolttableB(PocB(1)))...
      +(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))+ (2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableC(ParentC,3))...
      +(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))+ (2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableC(ParentC,4))...
      +((0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1)))*xintial(VolttableC(PocC(1))));
  
   
          elseif ((Tableii_0(Parent,3) == 0))   
  
  %%%%IaIb
        Aeq(Tableii(Parent,3), Tableii(Parent,3)) = 1;
        Aeq(Tableii(Parent,3), TableA(ParentA,3)) = 0;  
        Aeq(Tableii(Parent,3), TableA(ParentA,4)) = 0;  
        Aeq(Tableii(Parent,3), VolttableA(PocA(1))) = 0;  
        Aeq(Tableii(Parent,3), TableB(ParentB,3)) = 0;     
        Aeq(Tableii(Parent,3), TableB(ParentB,4)) = 0;  
        Aeq(Tableii(Parent,3), VolttableB(PocB(1))) = 0;  
       beq(Tableii(Parent,3)) = 0;
  
    %%%%IaIc
       Aeq(Tableii(Parent,4), Tableii(Parent,4)) = 1;
       Aeq(Tableii(Parent,4), TableA(ParentA,3)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ));
  
       Aeq(Tableii(Parent,4), TableA(ParentA,4)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
  
       Aeq(Tableii(Parent,4), VolttableA(PocA(1))) = (0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1);
  
      Aeq(Tableii(Parent,4), TableC(ParentC,3)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
     
      Aeq(Tableii(Parent,4), TableC(ParentC,4)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
  
      Aeq(Tableii(Parent,4), VolttableC(PocC(1))) = (0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2);
  
      beq(Tableii(Parent,4)) = -(xintial(Tableii(Parent,4))+(-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) )))*xintial(TableA(ParentA,3)) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableA(ParentA,4)) ...
      + ((0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1))*xintial(VolttableA(PocA(1))) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableC(ParentC,3)) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableC(ParentC,4)) ...
      + ((0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1))));
  
   %%%%IbIc
   
   Aeq(Tableii(Parent,5), Tableii(Parent,5)) = 1;
   Aeq(Tableii(Parent,5), TableB(ParentB,3)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), TableB(ParentB,4)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), VolttableB(PocB(1))) = (0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2));
  
   Aeq(Tableii(Parent,5), TableC(ParentC,3)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
     
   Aeq(Tableii(Parent,5), TableC(ParentC,4)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), VolttableC(PocC(1))) = (0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1));
  
  
   beq(Tableii(Parent,5)) = -(xintial(Tableii(Parent,5))+(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,3))+ ( -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,4)) ...
      + ((0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2)))*xintial(VolttableB(PocB(1)))...
      +(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))+ (2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableC(ParentC,3))...
      +(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))+ (2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableC(ParentC,4))...
      +((0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1)))*xintial(VolttableC(PocC(1))));
        
          else
              
  %%%%IaIb
        Aeq(Tableii(Parent,3), Tableii(Parent,3)) = 1;
        Aeq(Tableii(Parent,3), TableA(ParentA,3)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
       Aeq(Tableii(Parent,3), TableA(ParentA,4)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
       Aeq(Tableii(Parent,3), VolttableA(PocA(1))) = (0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1);
  
      Aeq(Tableii(Parent,3), TableB(ParentB,3)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
     
      Aeq(Tableii(Parent,3), TableB(ParentB,4)) = -((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1))))));
  
      Aeq(Tableii(Parent,3), VolttableB(PocB(1))) = (0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2);
  
       beq(Tableii(Parent,3)) = -(xintial(Tableii(Parent,3))+ (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableA(ParentA,3))...
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableA(ParentA,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableA(ParentA,4))...
         + ((0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
         + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
         + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableB(PocB(1)))^-1))*xintial(VolttableA(PocA(1)))... 
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableB(ParentB,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,3))...
         + (-((0.5*Tableii_0(Parent,3)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))) ...
         + (2*xintial(TableB(ParentB,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,4))...
         + ((0.5*Tableii_0(Parent,3)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
         + xintial(TableA(ParentA,3))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
         + xintial(TableA(ParentA,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableB(PocB(1)))^-2))*xintial(VolttableB(PocB(1))));
  
    %%%%IaIc
       Aeq(Tableii(Parent,4), Tableii(Parent,4)) = 1;
       Aeq(Tableii(Parent,4), TableA(ParentA,3)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ));
  
       Aeq(Tableii(Parent,4), TableA(ParentA,4)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
  
       Aeq(Tableii(Parent,4), VolttableA(PocA(1))) = (0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1);
  
      Aeq(Tableii(Parent,4), TableC(ParentC,3)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
     
      Aeq(Tableii(Parent,4), TableC(ParentC,4)) = -((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))));
  
      Aeq(Tableii(Parent,4), VolttableC(PocC(1))) = (0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2);
  
      beq(Tableii(Parent,4)) = -(xintial(Tableii(Parent,4))+(-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) ...
      + (2*xintial(TableA(ParentA,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))) )))*xintial(TableA(ParentA,3)) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableA(ParentA,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableA(ParentA,4)) ...
      + ((0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-2*xintial(VolttableC(PocC(1)))^-1))*xintial(VolttableA(PocA(1))) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableC(ParentC,3))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableC(ParentC,3)) ...
      + (-((0.5*Tableii_0(Parent,4)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,3))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1))))+ (2*xintial(TableC(ParentC,4))*xintial(TableA(ParentA,4))^2)/(xintial(VolttableA(PocA(1)))*xintial(VolttableC(PocC(1)))))))*xintial(TableC(ParentC,4)) ...
      + ((0.5*Tableii_0(Parent,4)^(-0.5))*(xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2 ...
      + xintial(TableA(ParentA,4))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableA(PocA(1)))^-1*xintial(VolttableC(PocC(1)))^-2))*xintial(VolttableC(PocC(1))));
  
   %%%%IbIc
   
   Aeq(Tableii(Parent,5), Tableii(Parent,5)) = 1;
   Aeq(Tableii(Parent,5), TableB(ParentB,3)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), TableB(ParentB,4)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), VolttableB(PocB(1))) = (0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2));
  
   Aeq(Tableii(Parent,5), TableC(ParentC,3)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
     
   Aeq(Tableii(Parent,5), TableC(ParentC,4)) = -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))));
  
   Aeq(Tableii(Parent,5), VolttableC(PocC(1))) = (0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1));
  
  
   beq(Tableii(Parent,5)) = -(xintial(Tableii(Parent,5))+(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,3))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,3))+ ( -((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))) ...
      + (2*xintial(TableB(ParentB,4))*xintial(TableC(ParentC,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableB(ParentB,4)) ...
      + ((0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-1*xintial(VolttableB(PocB(1)))^-2)))*xintial(VolttableB(PocB(1)))...
      +(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))+ (2*xintial(TableC(ParentC,3))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableC(ParentC,3))...
      +(-((0.5*Tableii_0(Parent,5)^(-0.5))*((2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,3))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1))))+ (2*xintial(TableC(ParentC,4))*xintial(TableB(ParentB,4))^2)/(xintial(VolttableC(PocC(1)))*xintial(VolttableB(PocB(1)))))))*xintial(TableC(ParentC,4))...
      +((0.5*Tableii_0(Parent,5)^(-0.5))*((xintial(TableC(ParentC,3))^2*xintial(TableB(ParentB,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableB(ParentB,3))^2*xintial(TableC(ParentC,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 + xintial(TableB(ParentB,4))^2*xintial(TableC(ParentC,3))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1 ...
      + xintial(TableC(ParentC,4))^2*xintial(TableB(ParentB,4))^2*xintial(VolttableC(PocC(1)))^-2*xintial(VolttableB(PocB(1)))^-1)))*xintial(VolttableC(PocC(1))));
   
          end
   
    end
 
end

 Aeq;
 beq;
 
indexnewvar = size(Aeq,2);
%%  

newvar =   zeros(size(Aeq,1),2*size(Aeq,1));         %% addition of columns due to addition of control variable p and n
Aeq = horzcat(Aeq,newvar);
for i = 1:size(Aeq,1)
    Aeq(i, indexnewvar+i) = -1;
    Aeq(i, indexnewvar+(size(newvar,2))/2+i) = 1;
end

Aeq(274,1509) = 0;
Aeq(274,2286) = 0;
Aeq(512,1747) = 0;
Aeq(512,2524) = 0;
Aeq(777,2012) = 0;
Aeq(777,2789) = 0;

dgindx = size(Aeq,2);

%% addition of DG Q variable 3706:3735

p =1;
for k = 1:size(powerdata,1)
    if powerdata(k,11) ~= 0 && powerdata(k,12) ~= 0 && powerdata(k,13) ~= 0
        ParentA = find(k == TA(:,2));
         Aeq(TableA(ParentA,4),dgindx+p) = 1;
         beq(TableA(ParentA,4)) = beq(TableA(ParentA,4))- xintial(dgindx+p);
         p = p+1;
         ParentB = find(k == TB(:,2));
         Aeq(TableB(ParentB,4),dgindx+p) = 1;
         beq(TableB(ParentB,4)) = beq(TableB(ParentB,4))- xintial(dgindx+p);
         p = p+1;
         ParentC = find(k == TC(:,2));
         Aeq(TableC(ParentC ,4),dgindx+p) = 1;
         beq(TableC(ParentC,4)) = beq(TableA(ParentC,4))- xintial(dgindx+p);
         p = p+1;
    elseif powerdata(k,11) ~= 0
        ParentA = find(k == TA(:,2));
        Aeq(TableA(ParentA ,4),dgindx+p) = 1;
        beq(TableA(ParentA,4)) = beq(TableA(ParentA,4))- xintial(dgindx+p);
        p = p+1;
       
    elseif powerdata(k,12) ~= 0
        ParentB = find(k == TB(:,2));
        Aeq(TableB(ParentB ,4),dgindx+p) = 1;
        beq(TableB(ParentB,4)) = beq(TableB(ParentB,4))- xintial(dgindx+p);
        p = p+1;
        
    elseif powerdata(k,13) ~= 0
        ParentC = find(k == TC(:,2));
        Aeq(TableC(ParentC ,4),dgindx+p) = 1;
        beq(TableC(ParentC,4)) = beq(TableC(ParentC,4))- xintial(dgindx+p);
        p = p+1;
    end
        
end

totalvar = size(Aeq,2);
tcontvar = totalvar - dgindx;
%% bound on Q due to original bound
pvmult = 1;

S_1 = 1.3*pvmult*0.04;

Q_limit = sqrt(S_1^2 - PGA(11)^2); 
% 

% A= zeros(3026,totalvar);
% b= zeros(3026,1);
% A= zeros(3106,totalvar);
% b= zeros(3106,1);
A= zeros(3166,totalvar);
b= zeros(3166,1);

p=1;
indexA = 0;
for i = 1:tcontvar
    
A(indexA+p,dgindx+i) = 1;
A(indexA+p+1,dgindx+i) = -1;
b(indexA+p) = Q_limit;
b(indexA+p+1) = Q_limit; 
p = p+2;
    
end

%% bound on Q using delta
p=1;
% indexA = 20;
% indexA = 60;
indexA = 90;
for i = 1:tcontvar
    
A(indexA+p,dgindx+i) = 1;
A(indexA+p+1,dgindx+i) = -1;
b(indexA+p) = delta;
b(indexA+p+1) = delta; 
p = p+2;    
end
% indexA= 40;
% indexA= 120;
indexA= 180;


%% bound on V

for i =2:nb
   ParentA = find(i == TA(:,2));                                                         
   A(indexA+ParentA,TableA(ParentA,6)) = -1;
   b(indexA+ParentA)= xintial(TableA(ParentA,6)) -0.9025;
end
% indexA= 131;
% indexA= 211;
indexA= 271;

for i =2:nb
   ParentA = find(i == TA(:,2));
%    Poc = find(T(Parent,1) == T(:,1));
   A(ParentA+indexA,TableA(ParentA,6)) = 1;
   b(ParentA+indexA)=1.1025- xintial(TableA(ParentA,6));
end
% indexA= 222;
% indexA= 302;
indexA= 362;

for i =2:nb
   ParentB = find(i == TB(:,2));                                                         
   A(indexA+ParentB,TableB(ParentB,6)) = -1;
   b(indexA+ParentB)= xintial(TableB(ParentB,6)) -0.9025;
end
% indexA= 301;
% indexA= 381;
indexA= 441;

for i =2:nb
   ParentB = find(i == TB(:,2));
   A(ParentB+indexA,TableB(ParentB,6)) = 1;
   b(ParentB+indexA)=1.1025- xintial(TableB(ParentB,6));
end
% indexA= 380;
% indexA= 460;
indexA=520;
for i =2:nb
   ParentC = find(i == TC(:,2));                                                         
   A(indexA+ParentC,TableC(ParentC,6)) = -1;
   b(indexA+ParentC)= xintial(TableC(ParentC,6)) -0.9025;
end
% indexA= 468;
% indexA= 548;
indexA= 608;

for i =2:nb
   ParentC = find(i == TC(:,2));
   A(ParentC+indexA,TableC(ParentC,6)) = 1;
   b(ParentC+indexA)=1.1025- xintial(TableC(ParentC,6));
end
% indexA= 556;
% indexA= 636;
indexA= 696;

for i = 1:size( newvar,2)
    A(indexA+i,indexnewvar+i ) = -1;
    b(indexA+i) = 0;
end
    
% totalvar = size(Aeq,2);
f = zeros(totalvar,1);
f(TableA(1,3)) = 1;
f(TableB(1,3)) = 1;
f(TableC(1,3)) = 1;

for i = 1:size( newvar,2)
    f(indexnewvar+i)=10000;
end
f;

for i=1:totalvar             %% total # of P,Q V
ctype(i)='C';
end
% 
options = cplexoptimset;
    options.Display = 'on';
  [delx,dfvalnlin] = cplexmilp(f,A,b',Aeq,beq,[],[],[],[],[],ctype);
  

 
%% iteration 0 is needed as linear to nonlinear

if iteration == 0
 apred = ((xintial(TableA(1,3))+delx(TableA(1,3))+ xintial(TableB(1,3))+delx(TableB(1,3))+xintial(TableC(1,3)))+delx(TableC(1,3)))- ((xintial(TableA(1,3))+ xintial(TableB(1,3))+xintial(TableC(1,3))));
else
apred = (((xintial(TableA(1,3))+ xintial(TableB(1,3))+xintial(TableC(1,3))))-((xintial(TableA(1,3))+delx(TableA(1,3))+ xintial(TableB(1,3))+delx(TableB(1,3))+xintial(TableC(1,3)))+delx(TableC(1,3))));
end

iteration= iteration+1;
alpha = 1;
xnew = xintial(1:1235)+alpha*delx(1:1235)';
xnew(3705+1:3705+tcontvar) = xintial(3705+1:3705+tcontvar)+alpha*delx(3705+1:3705+tcontvar)';

if iteration == 0
    error =  abs(((xintial(TableA(1,3))+delx(TableA(1,3))+ xintial(TableB(1,3))+delx(TableB(1,3))+xintial(TableC(1,3)))+delx(TableC(1,3)))- ((xintial(TableA(1,3))+ xintial(TableB(1,3))+xintial(TableC(1,3)))));
else
    error = abs((((xintial(TableA(1,3))+ xintial(TableB(1,3))+xintial(TableC(1,3))))-((xintial(TableA(1,3))+delx(TableA(1,3))+ xintial(TableB(1,3))+delx(TableB(1,3))+xintial(TableC(1,3)))+delx(TableC(1,3)))));
end

xintial = xnew;
if apred > 0
     delta = delta/2;
elseif apred < 0
    delta = 2*delta;
else
    delta = delta;
end

fvalnonlin = xnew(TableA(1,3))+xnew(TableB(1,3))+xnew(TableC(1,3));  %% obtaining fvalnonlinear
xnonlin = xnew;

end

 toc
 
