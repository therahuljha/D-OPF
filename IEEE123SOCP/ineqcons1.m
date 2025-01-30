 
function [c, ceq] = ineqcons1(x)
load linedata124C.txt
load branchphA.txt;
load branchphB.txt;
load branchphCC.txt;
lineA= branchphA;
lineB= branchphB;
lineC= branchphCC;
branch=linedata124C;

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
nb = 125;

T=dfsearch(G,1,'edgetonew');
TA=dfsearch(GA,1,'edgetonew');
TB=dfsearch(GB,1,'edgetonew');
TC=dfsearch(GC,1,'edgetonew');

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

Aii= zeros(nb,3) ;
nii = zeros(nb,3) ;
p = ((4*(nbusA-1)+1) + 4*(nbusB-1) +1) + 4*(nbusC-1) +1;
nii_1 =0;
for i = 2:nb
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
  if (~isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC))
     for j = 1:3
         row = find(i == T(:,2));
      Aii (row,j) = p +j;
      nii(row,j) = nii_1+j;
     end
        p = Aii (row,j);
        nii_1 = nii(row,j);
  elseif (isempty(ParentB))&& (~isempty(ParentA)) && (~isempty(ParentC))
    row = find(i == T(:,2));
      Aii(row,2)=p+1;
      p = Aii(row,2);
      nii(row,2) =  nii_1+1;
      nii_1 = nii(row,2);
  elseif (isempty(ParentC))&& (~isempty(ParentA)) && (~isempty(ParentB))
    row = find(i == T(:,2));
      Aii(row,1)=p+1;
      p = Aii(row,1);
      nii(row,1) =  nii_1+1;
      nii_1 = nii(row,1);
  elseif (isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC))
      Aii(row,3)=p+1;
      p = Aii(row,3);
      nii(row,3) =  nii_1+1;
      nii_1 = nii(row,3);
  end
end

Aii;
Aii(nb,:) =[]  ;    %% number of variables due to IaIb , IaIc, IbIc
Tableii = [T Aii];
 
nii(nb,:)=[];
Tableii_i = [T nii];


 ceq = [];
for i =2:nb
    
    ParentA = find(i == TA(:,2));
    if isempty(ParentA)  
       c = c;
    else
       PocA = find(TA(ParentA,1) == TA(:,1));
        c(ParentA) =(x(TableA(ParentA,3))^2 +  x(TableA(ParentA,4))^2) - x(TableA(ParentA,5))*x(VolttableA(PocA(1)));
      end

    ParentB = find(i == TB(:,2));      
      if isempty(ParentB)
          c = c;
      else
          PocB = find(TB(ParentB,1) == TB(:,1));
          c((nbusA-1)+ParentB) = (x(TableB(ParentB,3))^2 +  x(TableB(ParentB,4))^2) -x(TableB(ParentB,5))*x(VolttableB(PocB(1)));
      end

    ParentC = find(i == TC(:,2));      
      if isempty(ParentC)
          c = c;
      else
          PocC = find(TC(ParentC,1) == TC(:,1));
          c((nbusA-1)+(nbusB-1)+ParentC) = (x(TableC(ParentC,3))^2 +  x(TableC(ParentC,4))^2) - x(TableC(ParentC,5))*x(VolttableC(PocC(1)));
      end

    Parent = find(i == T(:,2));
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    if ((~isempty(ParentA)) && (~isempty(ParentB)) && (~isempty(ParentC)))
         c(Tableii_i(Parent,3)+(nbusA-1)+(nbusB-1)+(nbusC-1))  = (x(Tableii(Parent,3))^2 - x(TableA(ParentA,5))*x(TableB(ParentB,5))) ;
         c(Tableii_i(Parent,4)+(nbusA-1)+(nbusB-1)+(nbusC-1))  = (x(Tableii(Parent,4))^2 - x(TableA(ParentA,5))*x(TableC(ParentC,5))) ;
         c(Tableii_i(Parent,5)+(nbusA-1)+(nbusB-1)+(nbusC-1))  = (x(Tableii(Parent,5))^2 - x(TableC(ParentC,5))*x(TableB(ParentB,5))) ;
    elseif ((~isempty(ParentA)) && (~isempty(ParentB)))
         c(Tableii_i(Parent,3)+(nbusA-1)+(nbusB-1)+(nbusC-1))  = (x(Tableii(Parent,3))^2 - x(TableA(ParentA,5))*x(TableB(ParentB,5))) ;
    elseif ((~isempty(ParentA)) && (~isempty(ParentC)))
         c(Tableii_i(Parent,4)+(nbusA-1)+(nbusB-1)+(nbusC-1))  = (x(Tableii(Parent,4))^2 - x(TableA(ParentA,5))*x(TableC(ParentC,5))) ;
    elseif ((~isempty(ParentB)) && (~isempty(ParentC)))
         c(Tableii_i(Parent,5)+(nbusA-1)+(nbusB-1)+(nbusC-1))  = (x(Tableii(Parent,5))^2 - x(TableC(ParentC,5))*x(TableB(ParentB,5)));
    end

 
end
  
  
  
  
  
 
 
 