function X=getptsonsphere(mu,P,N)
n=length(mu);
if N>n
    error('N has to less than n')
end
    
A=sqrtm(eye(n));
%***************  PA   ******************************
if N==1
    A=eye(n);
    X=zeros(2*n,n);
    for i=1:1:n
            
            X(i,:)=A(:,i);
            X(i+n,:)=-A(:,i);
            
   
    end
    
    
end

if N==n
    
    %*******generating the CA direction***********
    index=GenerateIndex(n,n*ones(1,n));
    [roww,coll]=size(index);
    dr=[];
    for i=1:1:roww
        if length(find(index(i,:)>2))==0
            dr=vertcat(dr,index(i,:));
        end
    end
    
    %**************** CA - Space diagnols**************
    X=zeros(2^n,n);
    mo=-1*ones(1,n);
    for i=1:1:2^n
        rr=mo.^dr(i,:);
        sig=0;
        for j=1:1:n
            sig=sig+rr(j)*A(:,j);
        end
        X(i,:)=sig;
    end
    
elseif N~=1 % plane diagonal
    
    
    index=GenerateIndex(n,n*ones(1,n));
    
    dr=[];
    for i=N+1:1:n
        index(find(index==i))=0;
    end
    
    [roww,coll]=size(index);
    for i=1:1:roww
        if length(find(index(i,:)==0))==n-N
            dr=vertcat(dr,index(i,:));
        end
    end
    [rowwdr,coll]=size(dr);
    drr=dr(1,:);
    for i=1:1:rowwdr
        [rdr,coll]=size(drr);
        dd=0;
        for j=1:1:rdr
            dd(j)=sum(abs(drr(j,:)-dr(i,:)));
        end
        
        if length(find(dd==0))==0
            drr=[drr;dr(i,:)];
        end
    end
    drr(find(drr==2))=-1;
    X=zeros(size(drr,1),n);
    for i=1:1:size(drr,1)
        sig=0;
        for j=1:1:n
            sig=sig+drr(i,j)*A(:,j);
        end
        
        X(i,:)=sig;
    end
    
    
    
end

X=X./repmat(sqrt(sum(X.^2,2)),1,n);

A=sqrtm(P);
for i=1:size(X,1)
   X(i,:)= A*X(i,:)'+mu;
end