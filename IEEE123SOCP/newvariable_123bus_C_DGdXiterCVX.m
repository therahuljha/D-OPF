clc;
clear all;
t=97;


tic
load linedata124C.txt
load branchphA.txt;
load branchphB.txt;
load branchphCC.txt;
lineA= branchphA;
lineB= branchphB;
lineC= branchphCC;
branch=linedata124C;

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
% plot(GA)

fbB = lineB(:,1);
tbB = lineB(:,2);
nbusB =  length(lineB) +1;
GB = graph(fbB,tbB);
VB = dfsearch(GB,1);
% figure
% plot(GB)

fbC = lineC(:,1);
tbC = lineC(:,2);
nbusC =  length(lineC) +1;

GC = graph(fbC,tbC);
VC = dfsearch(GC,1);
% figure
% plot(GC)
tnbA = length(fbA);
tnbB = length(fbB);
tnbC = length(fbC);
nbusC =  length(lineC) +1;
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

%%
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

Aii= zeros(nb-1,3) ;
nii= zeros(nb,3) ;
p = ((4*(nbusA-1)+1) + 4*(nbusB-1) +1) + 4*(nbusC-1) +1;
Aii_1 =0;
for i = 2:nb
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    if (~isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC))
        for j = 1:3
            row = find(i == T(:,2));
           nii (row,j) = p +j;
           Aii(row,j) = Aii_1+j;
        end      
        p = nii (row,j);
        Aii_1 = Aii(row,j);
    elseif (isempty(ParentB))&& (~isempty(ParentA)) && (~isempty(ParentC))        
        row = find(i == T(:,2));
        nii(row,2)=p+1;
        Aii(row,2) =  Aii_1+1;
        p = nii(row,2);
        Aii_1 = Aii(row,2);
    elseif (isempty(ParentC))&& (~isempty(ParentA)) && (~isempty(ParentB))
        row = find(i == T(:,2));
        nii(row,1)=p+1;
        Aii(row,1) =  Aii_1+1;
        p = nii(row,1);
        Aii_1 = Aii(row,1);
    elseif (isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC))
         row = find(i == T(:,2));
        nii(row,3)=p+1;
        Aii(row,3) =  Aii_1+1;
        p = nii(row,3);
        Aii_1 = Aii(row,3);
    end
end
nii;
nii(nb,:) =[]  ;    %% number of variables due to IaIb , IaIc, IbIc
Tableii = [T nii];
Tableii_i = [T Aii ];

%%
[xnl,VAlin, VBlin, VClin, xlin, fvallin,TableAlin,TableBlin,TableClin, VolttableAlin,VolttableBlin,VolttableClin,error1,fvalue1] = firstSOCPforintialpointusingCVX(t);

 %% calling angle 
 [I_angle] = fbs123;
 I_angle = cell2mat(I_angle);

vs = 1.1025;

%% intitilizing the P,Q,V
xintial = xnl;

%%

error = 1;
iteration=0;
tolerance = 0.001;

while error > tolerance

% CVR_P = 0.6;                %%% CVR factor for P = 0.6
% CVR_Q = 3;                  %%% CVR factor for Q = 3

CVR_P = 2.0;                %%% CVR factor for P = 0.6
CVR_Q = 2.0; 

Aeq = zeros(777,1235);
beq = zeros(777,1);

 for i =2:nb
    k = zA(i) ;
    row = find(i == TA(:,1));
 
  if isempty(row)
    Parent = find(i == T(:,2));
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
 if isempty(ParentA)  
    Aeq = Aeq;
 elseif ((isempty(ParentB)) && (isempty(ParentC)) )  
    PocA = find(TA(ParentA,1) == TA(:,1));
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
       
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
      
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
   
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
   beq(ParentA+(nbusA-1))= (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2));
   beq(ParentA+2*(nbusA-1))= 0;
   
 elseif isempty(ParentB)   
     PocA = find(TA(ParentA,1) == TA(:,1));
    PocC = find(TC(ParentC,1) == TC(:,1)); 
    
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA,Tableii(Parent,4))= -(RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));
    
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,4))= - (XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));
  
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,5))= -((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),Tableii(Parent,4))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   
  beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));  
   beq(ParentA+(nbusA-1)) =  (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2));
   beq(ParentA+2*(nbusA-1)) = 0;
   
elseif isempty(ParentC) 
   PocA = find(TA(ParentA,1) == TA(:,1));
   PocB = find(TB(ParentB,1) == TB(:,1));
  
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA,Tableii(Parent,3))= -(RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));    
    
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,3))= - (XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));   
  
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,5))= -((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),Tableii(Parent,3))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));    
   beq(ParentA+(nbusA-1)) =  (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2));
   beq(ParentA+2*(nbusA-1)) = 0;
     
  
else
  PocA = find(TA(ParentA,1) == TA(:,1));
  PocB = find(TB(ParentB,1) == TB(:,1));
  PocC = find(TC(ParentC,1) == TC(:,1));
  
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA,Tableii(Parent,3))= -(RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
    Aeq(ParentA,Tableii(Parent,4))= -(RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));
    
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,3))= - (XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,4))= - (XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));
  
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,5))= -((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,5))= -((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),Tableii(Parent,3))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(ParentA+2*(nbusA-1),Tableii(Parent,4))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(ParentA+2*(nbusA-1),Tableii(Parent,5))= - 2*((((RAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));     
   beq(ParentA+(nbusA-1)) =  (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2));
   beq(ParentA+2*(nbusA-1)) = 0; 
     
 end
 
  else  
    Parent = find(i == T(:,2));
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    
 if isempty(ParentA)  
    Aeq = Aeq;
    
 elseif ((isempty(ParentB)) && (isempty(ParentC)))  
    PocA = find(TA(ParentA,1) == TA(:,1));
    
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
        end
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
    
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -1;
    end
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
   Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
     
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
     
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2));
         
   beq(ParentA+2*(nbusA-1)) = 0;
   
 elseif isempty(ParentB)   
    PocA = find(TA(ParentA,1) == TA(:,1));
    PocC = find(TC(ParentC,1) == TC(:,1));
    
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
        end
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA,Tableii(Parent,4))= - (RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));
    
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -1;
    end
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,4))= - (XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));
  
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,5))= -((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),Tableii(Parent,4))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
         
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2));
   
   beq(ParentA+2*(nbusA-1)) = 0;     
    
   
elseif isempty(ParentC) 
  PocA = find(TA(ParentA,1) == TA(:,1));
  PocB = find(TB(ParentB,1) == TB(:,1)); 
  
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
        end
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA,Tableii(Parent,3))= - (RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
   
    
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -1;
    end
   Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,3))= - (XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
    
  
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,5))= -((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),Tableii(Parent,3))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
       
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2));

   beq(ParentA+2*(nbusA-1)) = 0 ;         
     
   
 else
 
    PocA = find(TA(ParentA,1) == TA(:,1));
    PocB = find(TB(ParentB,1) == TB(:,1));
    PocC = find(TC(ParentC,1) == TC(:,1));
  
    Aeq(ParentA,TableA(ParentA,3))= 1;
    Aeq(ParentA,TableA(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
        end
    Aeq(ParentA,TableA(ParentA,6))= -(CVR_P/2)*PLA(TableA(ParentA,2));    
    Aeq(ParentA,TableA(ParentA,5))= -(RAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA,Tableii(Parent,3))= -(RAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) - XAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
    Aeq(ParentA,Tableii(Parent,4))= -(RAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) - XAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));
    
    Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
    Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -1;
    end
    Aeq(ParentA+(nbusA-1),(TableA(ParentA,6)))= -(CVR_Q/2)*QLA(TableA(ParentA,2)); 
    Aeq(ParentA+(nbusA-1),TableA(ParentA,5))= -(XAA(TA((ParentA),1),TA((ParentA),2)));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,3))= - (XAB(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(2,Parent))- (I_angle(1,Parent))) + RAB(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(2,Parent))- (I_angle(1,Parent))));
    Aeq(ParentA+(nbusA-1),Tableii(Parent,4))= - (XAC(TA((ParentA),1),TA((ParentA),2))*cosd((I_angle(3,Parent))- (I_angle(1,Parent))) + RAC(TA((ParentA),1),TA((ParentA),2))*sind((I_angle(3,Parent))- (I_angle(1,Parent))));
  
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,6))= 1;
   Aeq(ParentA+2*(nbusA-1),VolttableA(PocA(1)))= -1;
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
   Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= -((RAA(TA((ParentA),1),TA((ParentA),2)))^2 + (XAA(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),TableB(ParentB,5))= -((RAB(TA((ParentA),1),TA((ParentA),2)))^2 + (XAB(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),TableC(ParentC,5))= -((RAC(TA((ParentA),1),TA((ParentA),2)))^2 + (XAC(TA((ParentA),1),TA((ParentA),2)))^2) ;
   Aeq(ParentA+2*(nbusA-1),Tableii(Parent,3))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAB(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAB(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(ParentA+2*(nbusA-1),Tableii(Parent,4))= - 2*((((RAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAA(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAA(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(ParentA+2*(nbusA-1),Tableii(Parent,5))= - 2*((((RAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))+(XAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TA((ParentA),1),TA((ParentA),2)))*(RAC(TA((ParentA),1),TA((ParentA),2)))-(RAB(TA((ParentA),1),TA((ParentA),2)))*(XAC(TA((ParentA),1),TA((ParentA),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
   beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
      
   beq(ParentA+(nbusA-1)) = (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2));
    
   beq(ParentA+2*(nbusA-1)) = 0;            


end
    end
 Aeq(3*(nbusA-1)+1,VolttableA(1)) = 1;
 beq(3*(nbusA-1)+1) = vs;
 
 end
Aeq;

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
   
     PocB = find(TB(ParentB,1) == TB(:,1));
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3)) = 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));   
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
   
    
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    
  
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+2*(nbusB-1)) = 0;
   
 elseif isempty(ParentA)   
    PocB = find(TB(ParentB,1) == TB(:,1));
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3)) = 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));   
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB,Tableii(Parent,5))= -(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,5))= - (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
  
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,5))= -((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),Tableii(Parent,5))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2)); 
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+2*(nbusB-1)) = 0;
   
   
elseif isempty(ParentC) 
   PocB = find(TB(ParentB,1) == TB(:,1));
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3)) = 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));   
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB,Tableii(Parent,3))= - (RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
   
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,3))= - (XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
     
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,5))= -((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),Tableii(Parent,3))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+2*(nbusB-1)) = 0;
  
  
 else
     
   PocB = find(TB(ParentB,1) == TB(:,1));
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3)) = 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));   
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB,Tableii(Parent,3))= - (RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
    Aeq(3*(nbusA-1)+1+ParentB,Tableii(Parent,5))= -(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,3))= - (XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,5))= - (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
  
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,5))= -((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,5))= -((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),Tableii(Parent,3))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),Tableii(Parent,4))= - 2*((((RAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),Tableii(Parent,5))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2)); 
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2));
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) =  0;
   
   
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
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
            Aeq(3*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =  -1;
        end
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));    
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    
    
    %%% Q equations  
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
        for j = 1:length(row)-1
            Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -1;
        end
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
   
  %%% V equations  
    
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,5))= -((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2) ;
     
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
          
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2));
          
   beq(3*(nbusA-1)+1+ParentB+2*(nbusB-1)) = 0;
   
 elseif isempty(ParentA)   

    PocB = find(TB(ParentB,1) == TB(:,1));
%%% P equations  
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
        end
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));    
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB,Tableii(Parent,5))= -(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
    %%% Q equations  
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -1;
    end
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,5))= - (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
  
    %%% V equations  
    
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,5))= -((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),Tableii(Parent,5))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
   beq(3*(nbusA-1)+1+ParentB) = (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
      
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2));

   beq(3*(nbusA-1)+1+ParentB+2*(nbusB-1)) = 0;
     
elseif isempty(ParentC) 
    
PocB = find(TB(ParentB,1) == TB(:,1));
%%% P equations  
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
           Aeq(3*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
        end
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));    
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RAA(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB,Tableii(Parent,3))= -(RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
   
    
    %%% Q equations  
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -1;
    end
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,3))= -(XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
    
  %%% V equations  
  
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,5))= -((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),Tableii(Parent,3))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
  
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
    
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2));
    
   beq(3*(nbusA-1)+1+ParentB+2*(nbusB-1)) =  0;
   
   
else

PocB = find(TB(ParentB,1) == TB(:,1));
%%% P equations  
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;
    Aeq(3*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
        end
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,6))= -(CVR_P/2)*PLB(TableB(ParentB,2));    
    Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(RBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB,Tableii(Parent,3))= -(RAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) - XAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
    Aeq(3*(nbusA-1)+1+ParentB,Tableii(Parent,5))= -(RBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) - XBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
    
    %%% Q equations  
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -1;
    end
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,6)))= -(CVR_Q/2)*QLB(TableB(ParentB,2)); 
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,5))= -(XBB(TB((ParentB),1),TB((ParentB),2)));
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,3))= -(XAB(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))) + RAB(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(1,Parent))- (I_angle(2,Parent))));
    Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),Tableii(Parent,5))= - (XBC(TB((ParentB),1),TB((ParentB),2))*cosd((I_angle(3,Parent))- (I_angle(2,Parent))) + RBC(TB((ParentB),1),TB((ParentB),2))*sind((I_angle(3,Parent))- (I_angle(2,Parent))));
  %%% V equations  
  
  
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,6))= 1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(PocB(1)))= -1;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= -((RBB(TB((ParentB),1),TB((ParentB),2)))^2 + (XBB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,5))= -((RAB(TB((ParentB),1),TB((ParentB),2)))^2 + (XAB(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,5))= -((RBC(TB((ParentB),1),TB((ParentB),2)))^2 + (XBC(TB((ParentB),1),TB((ParentB),2)))^2) ;
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),Tableii(Parent,3))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RAB(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XAB(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBB(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBB(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),Tableii(Parent,4))= - 2*((((RAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RAB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),Tableii(Parent,5))= - 2*((((RBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))+(XBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBB(TB((ParentB),1),TB((ParentB),2)))*(RBC(TB((ParentB),1),TB((ParentB),2)))-(RBB(TB((ParentB),1),TB((ParentB),2)))*(XBC(TB((ParentB),1),TB((ParentB),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  
   beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
         
   beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2));
    
   beq(3*(nbusA-1)+1+ParentB+2*(nbusB-1)) = 0;

end
  end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1,VolttableB(1)) = 1;
    beq(3*(nbusA-1)+1+3*(nbusB-1)+1) = vs;

 
 end
Aeq;

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
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
   
    
    %%% Q equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
   
    
  %%% V equations     
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
 
   

   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+ 2*(nbusC-1)) = 0;
   
 elseif isempty(ParentA)   
  PocC = find(TC(ParentC,1) == TC(:,1));
%%% P equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,Tableii(Parent,5))= -(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
    %%% Q equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,5))= -(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
  %%% V equations     
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,5))= -((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),Tableii(Parent,5))= - 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  

   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2)); 
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1)) = 0;
   
   
elseif isempty(ParentB) 
   PocC = find(TC(ParentC,1) == TC(:,1));
%%% P equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,Tableii(Parent,4))= -(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
        
    %%% Q equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,4))= -(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))) + RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    
    
  %%% V equations     
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,5))= -((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),Tableii(Parent,4))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
  
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1)) = 0; 
           
  
 else
     
  PocC = find(TC(ParentC,1) == TC(:,1));
%%% P equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,Tableii(Parent,4))= -(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,Tableii(Parent,5))= -(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
    %%% Q equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentB,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,4))= -(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))) + RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,5))= -(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
  %%% V equations     
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,5))= -((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,5))= -((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),Tableii(Parent,3))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))-(RAB(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),Tableii(Parent,4))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),Tableii(Parent,5))= - 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  

   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) =  (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1)) = 0;


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
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
        end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));    
    
    %%% Q equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -1;
    end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));    
    
  %%% V equations     
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)=(1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
          
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2));
          
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1)) = 0;
   
 elseif isempty(ParentA)   

PocC = find(TC(ParentC,1) == TC(:,1));
%%% P equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
        end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,Tableii(Parent,5))= -(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
    %%% Q equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -1;
    end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,5))= -(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
  %%% V equations     
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,5))= -((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),Tableii(Parent,5))= - 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  

   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
      
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2));
      
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1)) = 0;
  
   
elseif isempty(ParentB) 
    PocC = find(TC(ParentC,1) == TC(:,1));
%%% P equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
        end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,Tableii(Parent,4))= -(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    
     %%% Q equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -1;
    end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,4))= -(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))+ RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    
    
  %%% V equations     
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,5))= -((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),Tableii(Parent,4))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
  
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2)); 
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2));        
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1)) = 0;   
   
else

PocC = find(TC(ParentC,1) == TC(:,1));
%%% P equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
        for j = 1:length(row)-1
            Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   - 1;
        end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,6))= -(CVR_P/2)*PLC(TableC(ParentC,2));    
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(RCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,Tableii(Parent,4))= -(RAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent)))-XAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,Tableii(Parent,5))= -(RBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) - XBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
    %%% Q equations  
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;    
    for j = 1:length(row)-1
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -1;
    end
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,6)))= -(CVR_Q/2)*QLC(TableC(ParentC,2)); 
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,5))= -(XCC(TC((ParentC),1),TC((ParentC),2)));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,4))= -(XAC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))) + RAC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(1,Parent))- (I_angle(3,Parent))));
    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),Tableii(Parent,5))= -(XBC(TC((ParentC),1),TC((ParentC),2))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))) + RBC(TC((ParentC),1),TC((ParentC),2))*sind((I_angle(2,Parent))- (I_angle(3,Parent))));
    
  %%% V equations     
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,6))= 1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(PocC(1)))= -1;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= -((RCC(TC((ParentC),1),TC((ParentC),2)))^2 + (XCC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,5))= -((RAC(TC((ParentC),1),TC((ParentC),2)))^2 + (XAC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,5))= -((RBC(TC((ParentC),1),TC((ParentC),2)))^2 + (XBC(TC((ParentC),1),TC((ParentC),2)))^2) ;
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),Tableii(Parent,3))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(2,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RBC(TC((ParentC),1),TC((ParentC),2)))-(RAB(TC((ParentC),1),TC((ParentC),2)))*(XBC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(2,Parent)))));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),Tableii(Parent,4))= - 2*((((RAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(1,Parent))- (I_angle(3,Parent))))...
       - (((XAC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RAC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(1,Parent))- (I_angle(3,Parent)))));
   Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),Tableii(Parent,5))= - 2*((((RBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))+(XBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*cosd((I_angle(2,Parent))- (I_angle(3,Parent))))...
       - (((XBC(TC((ParentC),1),TC((ParentC),2)))*(RCC(TC((ParentC),1),TC((ParentC),2)))-(RBC(TC((ParentC),1),TC((ParentC),2)))*(XCC(TC((ParentC),1),TC((ParentC),2))))*sind((I_angle(2,Parent))- (I_angle(3,Parent)))));
  

   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
 
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2));
   
   beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1)) = 0;

    
end
  end
   
  Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+3*(nbusC-1)+1,VolttableC(1)) = 1;
  beq(3*(nbusA-1)+1+3*(nbusB-1)+1+3*(nbusC-1)+1) = vs;
 
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
%% addition of DG Q variable

p =1;
for k = 1:size(powerdata,1)
    if powerdata(k,11) ~= 0 && powerdata(k,12) ~= 0 && powerdata(k,13) ~= 0
        ParentA = find(k == TA(:,2));
         Aeq(TableAlin(ParentA ,4),dgindx+p) = 1;
         p = p+1;
         ParentB = find(k == TB(:,2));
         Aeq(TableBlin(ParentB ,4),dgindx+p) = 1;
         p = p+1;
         ParentC = find(k == TC(:,2));
         Aeq(TableClin(ParentC ,4),dgindx+p) = 1;
         p = p+1;
    elseif powerdata(k,11) ~= 0
        ParentA = find(k == TA(:,2));
        Aeq(TableAlin(ParentA ,4),dgindx+p) = 1;
        p = p+1;
       
    elseif powerdata(k,12) ~= 0
        ParentB = find(k == TB(:,2));
        Aeq(TableBlin(ParentB ,4),dgindx+p) = 1;
        p = p+1;
        
    elseif powerdata(k,13) ~= 0
        ParentC = find(k == TC(:,2));
        Aeq(TableClin(ParentC ,4),dgindx+p) = 1;
        p = p+1;
    end
        
end

totalvar = size(Aeq,2);
tcontvar = totalvar - dgindx;


%%

A = zeros(458,totalvar);
b = zeros(458,1);

gamma = 0.9;
% gamma = 0.85;
for i =2:nb    
    ParentA = find(i == TA(:,2));
 
    if (~isempty(ParentA))  
       PocA = find(TA(ParentA,1) == TA(:,1));       
        A(ParentA,TableA(ParentA,3)) = -2*xintial(TableA(ParentA,3));
        A(ParentA,TableA(ParentA,4)) = -2*xintial(TableA(ParentA,4));
        A(ParentA,TableA(ParentA,5)) = xintial(VolttableA(PocA(1)));
        A(ParentA,VolttableA(PocA(1))) = xintial(TableA(ParentA,5));
        b(ParentA)= -(gamma + 1)*(xintial(TableA(ParentA,3))^2+xintial(TableA(ParentA,4))^2 -xintial(TableA(ParentA,5))*xintial(VolttableA(PocA(1))));

      end

    ParentB = find(i == TB(:,2));      
      if (~isempty(ParentB))
          PocB = find(TB(ParentB,1) == TB(:,1));
          A((nbusA-1)+ParentB,TableB(ParentB,3)) = -2*xintial(TableB(ParentB,3));
          A((nbusA-1)+ParentB,TableB(ParentB,4)) = -2*xintial(TableB(ParentB,4));
          A((nbusA-1)+ParentB,TableB(ParentB,5)) = xintial(VolttableB(PocB(1)));
          A((nbusA-1)+ParentB,VolttableB(PocB(1))) =  xintial(TableB(ParentB,5));
          b((nbusA-1)+ParentB)= -(gamma + 1)*(xintial(TableB(ParentB,3))^2+xintial(TableB(ParentB,4))^2 -xintial(TableB(ParentB,5))*xintial(VolttableB(PocB(1))));
      end

    ParentC = find(i == TC(:,2));      
      if (~isempty(ParentC))
          PocC = find(TC(ParentC,1) == TC(:,1));
          A((nbusA-1)+(nbusB-1)+ParentC,TableC(ParentC,3)) = -2*xintial(TableC(ParentC,3));
          A((nbusA-1)+(nbusB-1)+ParentC,TableC(ParentC,4)) = -2*xintial(TableC(ParentC,4));
          A((nbusA-1)+(nbusB-1)+ParentC,TableC(ParentC,5)) = xintial(VolttableC(PocC(1)));
          A((nbusA-1)+(nbusB-1)+ParentC,VolttableC(PocC(1))) =  xintial(TableC(ParentC,5));
          b((nbusA-1)+(nbusB-1)+ParentC)= -(gamma + 1)*(xintial(TableC(ParentC,3))^2+xintial(TableC(ParentC,4))^2 -xintial(TableC(ParentC,5))*xintial(VolttableC(PocC(1))));
      end

    Parent = find(i == T(:,2));
    if ((~isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC)))
         A(Tableii_i(Parent,3)+(nbusA-1)+(nbusB-1)+(nbusC-1),Tableii(Parent,3))  = -2*xintial(Tableii(Parent,3)) ;
         A(Tableii_i(Parent,3)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableA(ParentA,5))  = xintial(TableB(ParentB,5));
         A(Tableii_i(Parent,3)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableB(ParentB,5))  = xintial(TableA(ParentA,5));
         b(Tableii_i(Parent,3)+(nbusA-1)+(nbusB-1)+(nbusC-1))= -(gamma + 1)*(xintial(Tableii(Parent,3))^2 -xintial(TableA(ParentA,5))*xintial(TableB(ParentB,5)));
         
         A(Tableii_i(Parent,4)+(nbusA-1)+(nbusB-1)+(nbusC-1),Tableii(Parent,4))  = -2*xintial(Tableii(Parent,4)) ;
         A(Tableii_i(Parent,4)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableA(ParentA,5))  = xintial(TableC(ParentC,5));
         A(Tableii_i(Parent,4)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableC(ParentC,5))  = xintial(TableA(ParentA,5));
         b(Tableii_i(Parent,4)+(nbusA-1)+(nbusB-1)+(nbusC-1))= -(gamma + 1)*(xintial(Tableii(Parent,4))^2 -xintial(TableA(ParentA,5))*xintial(TableC(ParentC,5)));
         
         A(Tableii_i(Parent,5)+(nbusA-1)+(nbusB-1)+(nbusC-1),Tableii(Parent,5))  = -2*xintial(Tableii(Parent,5)) ;
         A(Tableii_i(Parent,5)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableB(ParentB,5))  = xintial(TableC(ParentC,5));
         A(Tableii_i(Parent,5)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableC(ParentC,5))  = xintial(TableB(ParentB,5));
         b(Tableii_i(Parent,5)+(nbusA-1)+(nbusB-1)+(nbusC-1))= -(gamma + 1)*(xintial(Tableii(Parent,5))^2 -xintial(TableB(ParentB,5))*xintial(TableC(ParentC,5)));         
         
    elseif ((~isempty(ParentA)) && (~isempty(ParentB)))
         A(Tableii_i(Parent,3)+(nbusA-1)+(nbusB-1)+(nbusC-1),Tableii(Parent,3))  = -2*xintial(Tableii(Parent,3)) ;
         A(Tableii_i(Parent,3)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableA(ParentA,5))  = xintial(TableB(ParentB,5));
         A(Tableii_i(Parent,3)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableB(ParentB,5))  = xintial(TableA(ParentA,5));
         b(Tableii_i(Parent,3)+(nbusA-1)+(nbusB-1)+(nbusC-1))= -(gamma + 1)*(xintial(Tableii(Parent,3))^2 -xintial(TableA(ParentA,5))*xintial(TableB(ParentB,5)));
         
    elseif ((~isempty(ParentA)) && (~isempty(ParentC)))
         A(Tableii_i(Parent,4)+(nbusA-1)+(nbusB-1)+(nbusC-1),Tableii(Parent,4))  = -2*xintial(Tableii(Parent,4)) ;
         A(Tableii_i(Parent,4)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableA(ParentA,5))  = xintial(TableC(ParentC,5));
         A(Tableii_i(Parent,4)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableC(ParentC,5))  = xintial(TableA(ParentA,5));
         b(Tableii_i(Parent,4)+(nbusA-1)+(nbusB-1)+(nbusC-1))= -(gamma + 1)*(xintial(Tableii(Parent,4))^2 -xintial(TableA(ParentA,5))*xintial(TableC(ParentC,5)));
         
    elseif ((~isempty(ParentB)) && (~isempty(ParentC)))
         A(Tableii_i(Parent,5)+(nbusA-1)+(nbusB-1)+(nbusC-1),Tableii(Parent,5))  = -2*xintial(Tableii(Parent,5)) ;
         A(Tableii_i(Parent,5)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableB(ParentB,5))  = xintial(TableC(ParentC,5));
         A(Tableii_i(Parent,5)+(nbusA-1)+(nbusB-1)+(nbusC-1),TableC(ParentC,5))  = xintial(TableB(ParentB,5));
         b(Tableii_i(Parent,5)+(nbusA-1)+(nbusB-1)+(nbusC-1))= -(gamma + 1)*(xintial(Tableii(Parent,5))^2 -xintial(TableB(ParentB,5))*xintial(TableC(ParentC,5))); 
    end
    
    A1 = zeros(4038,totalvar);
b1 = zeros(4038,1);
 
end

%% with 10%DG
A1 = zeros(4038,totalvar);
b1 = zeros(4038,1);

%% with 30%DG
% A1 = zeros(4078,totalvar);
% b1 = zeros(4078,1);

%% with 50%DG
% A1 = zeros(4108,totalvar);
% b1 = zeros(4108,1);

delP = 0.2;
% delP = 0.1;
indexA=0;
for i =2:nb
   ParentA = find(i == TA(:,2));                                                         
   A1(indexA+ParentA,TableA(ParentA,3)) = -1;
   b1(indexA+ParentA)= -xintial(TableA(ParentA,3)) + delP;
end

indexA=91;
for i =2:nb
   ParentA = find(i == TA(:,2));
   A1(ParentA+indexA,TableA(ParentA,3)) = 1;
   b1(ParentA+indexA)= delP + xintial(TableA(ParentA,3));
end

indexA=182;
for i =2:nb
   ParentA = find(i == TA(:,2));                                                         
   A1(indexA+ParentA,TableA(ParentA,4)) = -1;
   b1(indexA+ParentA)= -xintial(TableA(ParentA,4)) + delP;
end

 indexA=273;
for i =2:nb
   ParentA = find(i == TA(:,2));
   A1(ParentA+indexA,TableA(ParentA,4)) = 1;
   b1(ParentA+indexA)= delP + xintial(TableA(ParentA,4));
end

indexA=364;
for i =2:nb
   ParentA = find(i == TA(:,2));                                                         
   A1(indexA+ParentA,TableA(ParentA,5)) = -1;
    b1(indexA+ParentA)= -xintial(TableA(ParentA,5)) + delP;
end

 indexA=455;
for i =2:nb
   ParentA = find(i == TA(:,2));
   A1(ParentA+indexA,TableA(ParentA,5)) = 1;
   b1(ParentA+indexA)= delP + xintial(TableA(ParentA,5));
end

indexA=546;
for i =2:nb
   ParentB = find(i == TB(:,2));                                                         
   A1(indexA+ParentB,TableB(ParentB,3)) = -1;
   b1(indexA+ParentB)= -xintial(TableB(ParentB,3)) + delP;
end

indexA=625;
for i =2:nb
   ParentB = find(i == TB(:,2));
   A1(ParentB+indexA,TableB(ParentB,3)) = 1;
   b1(ParentB+indexA)= delP + xintial(TableB(ParentB,3));
end

indexA=704;
for i =2:nb
   ParentB = find(i == TB(:,2));                                                         
   A1(indexA+ParentB,TableB(ParentB,4)) = -1;
   b1(indexA+ParentB)= -xintial(TableB(ParentB,4)) + delP;
end

 indexA=783;
for i =2:nb
   ParentB = find(i == TB(:,2));
   A1(ParentB+indexA,TableB(ParentB,4)) = 1;
   b1(ParentB+indexA)= delP + xintial(TableB(ParentB,4));
end

indexA=862;
for i =2:nb
   ParentB = find(i == TB(:,2));                                                         
   A1(indexA+ParentB,TableB(ParentB,5)) = -1;
    b1(indexA+ParentB)= -xintial(TableB(ParentB,5)) + delP;

end

 indexA=941;
for i =2:nb
   ParentB = find(i == TB(:,2));
   A1(ParentB+indexA,TableB(ParentB,5)) = 1;
   b1(ParentB+indexA)= delP + xintial(TableB(ParentB,5));
end

indexA=1020;
for i =2:nb
   ParentC = find(i == TC(:,2));                                                         
   A1(indexA+ParentC,TableC(ParentC,3)) = -1;
   b1(indexA+ParentC)= -xintial(TableC(ParentC,3)) + delP;
end

indexA=1108;
for i =2:nb
   ParentC = find(i == TC(:,2));
   A1(ParentC+indexA,TableC(ParentC,3)) = 1;
   b1(ParentC+indexA)= delP + xintial(TableC(ParentC,3));
end

 
indexA=1196;
for i =2:nb
   ParentC = find(i == TC(:,2));                                                         
   A1(indexA+ParentC,TableC(ParentC,4)) = -1;
   b1(indexA+ParentC)= -xintial(TableC(ParentC,4)) + delP;
end

indexA=1284;
for i =2:nb
   ParentC = find(i == TC(:,2));
   A1(ParentC+indexA,TableC(ParentC,4)) = 1;
   b1(ParentC+indexA)= delP + xintial(TableC(ParentC,4));
end

deli= 0.50;
% deli= 0.30;
indexA=1372;
for i =2:nb
   ParentC = find(i == TC(:,2));                                                         
   A1(indexA+ParentC,TableC(ParentC,5)) = -1;
   b1(indexA+ParentC)= -xintial(TableC(ParentC,5)) + deli;
end

indexA=1460;
for i =2:nb
   ParentC = find(i == TC(:,2));
   A1(ParentC+indexA,TableC(ParentC,5)) = 1;
   b1(ParentC+indexA)= deli + xintial(TableC(ParentC,5));
end

indexA= 1548;
for i =2:nb
   ParentA = find(i == TA(:,2));                                                         
   A1(indexA+ParentA,TableA(ParentA,6)) = -1;
 b1(indexA+ParentA)= -0.9025;
end

indexA= 1639;

for i =2:nb
   ParentA = find(i == TA(:,2));
   A1(ParentA+indexA,TableA(ParentA,6)) = 1;
   b1(ParentA+indexA)= 1.1025;
end

indexA= 1730;

for i =2:nb
   ParentB = find(i == TB(:,2));                                                         
   A1(indexA+ParentB,TableB(ParentB,6)) = -1;
b1(indexA+ParentB)= -0.9025;
end
  
 indexA= 1809;

for i =2:nb
   ParentB = find(i == TB(:,2));
   A1(ParentB+indexA,TableB(ParentB,6)) = 1;
   b1(ParentB+indexA)=1.1025;
end

 indexA= 1888;

for i =2:nb
   ParentC = find(i == TC(:,2));                                                         
   A1(indexA+ParentC,TableC(ParentC,6)) = -1;
    b1(indexA+ParentC)= -0.9025;
end

 indexA= 1976;

for i =2:nb
   ParentC = find(i == TC(:,2));
   A1(ParentC+indexA,TableC(ParentC,6)) = 1;
   b1(ParentC+indexA)=1.1025;
end

iidev = 0.50;
% iidev = 0.30;
indexA=2064;
for i =2:nb
   Parent = find(i == T(:,2));
   ParentA = find(i == TA(:,2));
   ParentB = find(i == TB(:,2));
   ParentC = find(i == TC(:,2));
   
 if ((~isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC)))
   A1(Tableii_i(Parent,3)+indexA,Tableii(Parent,3)) = -1;
   b1(Tableii_i(Parent,3)+indexA)=iidev- xintial(Tableii(Parent,3));
   A1(Tableii_i(Parent,4)+indexA,Tableii(Parent,4)) = -1;
   b1(Tableii_i(Parent,4)+indexA)=iidev- xintial(Tableii(Parent,4));
   A1(Tableii_i(Parent,5)+indexA,Tableii(Parent,5)) = -1;
   b1(Tableii_i(Parent,5)+indexA)=iidev- xintial(Tableii(Parent,5));
   
    elseif ((~isempty(ParentA)) && (~isempty(ParentB)))
       A1(Tableii_i(Parent,3)+indexA,Tableii(Parent,3)) = -1;
       b1(Tableii_i(Parent,3)+indexA)=iidev- xintial(Tableii(Parent,3));

    elseif ((~isempty(ParentA))  && (~isempty(ParentC)))
        A1(Tableii_i(Parent,4)+indexA,Tableii(Parent,4)) = -1;
        b1(Tableii_i(Parent,4)+indexA)=iidev- xintial(Tableii(Parent,4));
        
    elseif ((~isempty(ParentB)) && (~isempty(ParentC)))
        A1(Tableii_i(Parent,5)+indexA,Tableii(Parent,5)) = -1;
        b1(Tableii_i(Parent,5)+indexA)=iidev- xintial(Tableii(Parent,5));
        
    end
end

indexA=2264;
for i =2:nb
   Parent = find(i == T(:,2));
   ParentA = find(i == TA(:,2));
   ParentB = find(i == TB(:,2));
   ParentC = find(i == TC(:,2));
   
    if ((~isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC)))
       A1(Tableii_i(Parent,3)+indexA,Tableii(Parent,3)) = 1;
       b1(Tableii_i(Parent,3)+indexA)=iidev+ xintial(Tableii(Parent,3));
       A1(Tableii_i(Parent,4)+indexA,Tableii(Parent,4)) = 1;
       b1(Tableii_i(Parent,4)+indexA)=iidev+ xintial(Tableii(Parent,4));
       A1(Tableii_i(Parent,5)+indexA,Tableii(Parent,5)) = 1;
       b1(Tableii_i(Parent,5)+indexA)=iidev+ xintial(Tableii(Parent,5));
   
    elseif ((~isempty(ParentA)) && (~isempty(ParentB)))
       A1(Tableii_i(Parent,3)+indexA,Tableii(Parent,3)) = 1;
       b1(Tableii_i(Parent,3)+indexA)=iidev+xintial(Tableii(Parent,3));

    elseif ((~isempty(ParentA))  && (~isempty(ParentC)))
        A1(Tableii_i(Parent,4)+indexA,Tableii(Parent,4)) = 1;
        b1(Tableii_i(Parent,4)+indexA)=iidev+ xintial(Tableii(Parent,4));
        
    elseif ((~isempty(ParentB)) && (~isempty(ParentC)))
        A1(Tableii_i(Parent,5)+indexA,Tableii(Parent,5)) = 1;
        b1(Tableii_i(Parent,5)+indexA)=iidev+ xintial(Tableii(Parent,5));
        
    end
end

indexA=2464;
pvmult = 1.0;
S = 1.3*pvmult*0.04;

Q_limit = sqrt(S^2 - PGA(11)^2); 

p=1;

for i = 1:tcontvar
    
A1(indexA+p,dgindx+i) = 1;
A1(indexA+p+1,dgindx+i) = -1;
b1(indexA+p) = Q_limit;
b1(indexA+p+1) = Q_limit; 
p = p+2;
    
end


 
indexA=2484;   %% with 10%DG
% indexA=2524;   %% with 30%DG

% indexA=2554;    %% with 50%DG

for i = 1:1554
    A1(indexA+i, indexnewvar+i) = -1;
    b1(indexA+i) = 0;
end


c1 = zeros(1,totalvar);
c1(1,1)=1;
c1(1,366)=1;
c1(1,683)=1;

for i = indexnewvar+1:indexnewvar+1554
    c1(1,i) = 1000000;
end


tic
%         cvx_solver MOSEK
     cvx_solver Gurobi
cvx_begin  

   variable x(totalvar);
    minimize(c1*(x));
    
subject to
for i =2:nb
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    Parent = find(i == T(:,2));
        if (~isempty(ParentA))  
           PocA = find(TA(ParentA,1) == TA(:,1));
           {[x(TableA(ParentA,3));  x(TableA(ParentA,4))], x(VolttableA(PocA(1))), x(TableA(ParentA,5))} == rotated_lorentz(2)
        end
        if (~isempty(ParentB))  
           PocB = find(TB(ParentB,1) == TB(:,1));
           {[x(TableB(ParentB,3));  x(TableB(ParentB,4))], x(VolttableB(PocB(1))), x(TableB(ParentB,5))} == rotated_lorentz(2)
        end
        if (~isempty(ParentC))  
           PocC = find(TC(ParentC,1) == TC(:,1));
           {[x(TableC(ParentC,3));  x(TableC(ParentC,4))], x(VolttableC(PocC(1))), x(TableC(ParentC,5))} == rotated_lorentz(2)
        end
        if ((~isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC)))
            {x(Tableii(Parent,3)),x(TableA(ParentA,5)),x(TableB(ParentB,5))}  == rotated_lorentz(1)
            {x(Tableii(Parent,4)),x(TableA(ParentA,5)),x(TableC(ParentC,5))}  == rotated_lorentz(1)
            {x(Tableii(Parent,5)),x(TableB(ParentB,5)),x(TableC(ParentC,5))}  == rotated_lorentz(1)
        end
        if ((~isempty(ParentA)) && (~isempty(ParentB)))
            {x(Tableii(Parent,3)),x(TableA(ParentA,5)),x(TableB(ParentB,5))}  == rotated_lorentz(1)
        end
        if ((~isempty(ParentA)) && (~isempty(ParentC)))
            {x(Tableii(Parent,4)),x(TableA(ParentA,5)),x(TableC(ParentC,5))}  == rotated_lorentz(1)
        end
        if ((~isempty(ParentB)) && (~isempty(ParentC)))
            {x(Tableii(Parent,5)),x(TableB(ParentB,5)),x(TableC(ParentC,5))}  == rotated_lorentz(1)
        end
end
    Aeq*x == beq;
    A*x <= b;
      A1*x <= b1;
%     x >= lb;
%     x <= ub;
    cvx_solver 
    
cvx_end
toc


telapse(:,1) = toc;

xintial(1:indexnewvar) = x(1:indexnewvar);

xintial(dgindx+1:dgindx+tcontvar) = x(dgindx+1:dgindx+tcontvar);


[c, ceq] = ineqcons1(x);

error = max(abs(c))

 iteration= iteration+1;
 

fvalue=xintial(1)+xintial(366)+xintial(683)
Qvalue=xintial(92)+xintial(445)+xintial(771);
error1 = [error1 error];
end
VAnon = sqrt(x(274:365));
VBnon = sqrt(x(603:682));
VCnon = sqrt(x(947:1035));