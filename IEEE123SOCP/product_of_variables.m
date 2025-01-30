function[A, b]=product_of_variables(A,b, vec1P,vec2P,from ,to,s)
% For zp1>=0;
    zp1=zeros(1,s);
    p=1;
    for i=from:to
        zp1(p,i)=-1;
        p=p+1;
    end
    A=[A;zp1];
    b=[b; zeros(32,1)]; 
    
    %% z>= A-(1-x)M
     p=1;
     M =100000;
    for i=from:to
        zp1(p,i)=-1;        
        zp1(p,vec2P(1))=1;
        zp1(p,vec1P(p))=M;
        p=p+1;
    end
    A=[A;zp1];
    b=[b; M*ones(32,1)];
    %% 
    % For zp1<=M*Binary;
    zp1=zeros(1,s);             %% phase A 
    p=1;
    for i=from:to
        zp1(p,i)=1;
        zp1(p,vec1P(p))=-100000;
         p=p+1;
    end
    A=[A;zp1];
    b=[b; zeros(32,1)];
    %% 
    %% 
    % For zp1<=A;
    zp1=zeros(1,s);             %% phase A 
    p=1;
    for i=from:to
        zp1(p,i)=1;
        zp1(p,vec2P(1))=-1;
         p=p+1;
    end
    A=[A;zp1];
    b=[b; zeros(32,1)];  
    
    
end








