% function
% divide users into G groups and calculate WW

function [X,WW]=usergrouping(rho,A,G)
[N,M]=size(A);

% initialize WW
%---------------------------------------------------------------------------
% Firstly, we randomly selected a user k1 and set WW_1 = A_k1 . 
% Then we calculate the metrics of the correlation between w? 1 and all the other users. 
% The user with the smallest metric is chosen to be k2 and set w? 2 = ak2 . 
% Again, the metrics about w? 1 and w? 2 are measured and the least correlated user k3 is determined. 
% Repeat these steps until G users are found.
ind=zeros(1,G);
ind(1)=randi(M);
for g=2:G
    met=zeros(g-1,M);
    for k=1:M
        for i=1:g-1
            met(i,k)=rho(k)^2*norm(A(:,k)'*A(:,ind(i)))^2-rho(k)^2*N;%n = norm(v) returns the 2-norm or Euclidean norm of vector v.
        end
    end
    if size(met,1)>1
        met=max(met);
    end
    ind(g)=find(met==min(met));%便于理解：请看Reference [7-22]    
end

WW=A(:,ind);
%---------------------------------------------------------------------------

eta=zeros(G,M);
flag=1;
while flag
    X=zeros(G,M);
    for k=1:M
        for g=1:G
            eta(g,k)=rho(k)^2*norm(A(:,k)'*WW(:,g))^2-rho(k)^2*N;
        end
        [val,lab]=max(eta(:,k));
        X(lab,k)=1;
    end
    
    % update WW
    WW0=centerupdate(WW,rho,A,X);
    
    % converge?
    if norm(WW-WW0,'fro')<1e-1
        flag=0;
    end
    
    WW=WW0;
end







