% 平流层项目 sum-rate
% zhangqian

clear all;
clc

% setting
N=100;      % # of antennas
M=200;     % # of users
P=10^(46/10);   % power in mW
sigma=sqrt(10^(-169/10)*10^7);  %  noise covariance  -169dBm/Hz  10MHz
height=20;  % antenna height in km
radius=20;  % user location radius in km
freq=2.4*10^9;  % frequency 2.4GHz
lambda=3*10^8/freq;  % wave length in m
antd=4*lambda;  % antenna spacing in m

% G_seq=2:N;
G_seq=10:10:N; 

tau_seq=1:-0.2:0.8;

ave_pro_sr=zeros(length(tau_seq),length(G_seq));
ave_ran_sr=zeros(1,length(G_seq));
ave_dist_sr=zeros(length(tau_seq),length(G_seq));

ave_pro_num=zeros(length(tau_seq),length(G_seq));
ave_dist_num=zeros(length(tau_seq),length(G_seq));

ave_pro_cput=zeros(1,length(G_seq));
ave_dist_cput=zeros(1,length(G_seq));


count=10;
runs_out=10;
runs_in=200;
runs=10;

for ite=1:count
        
    type=2;
    [loca,lab]=generateloca(M,radius,type);  % generate user location 2-by-M
    dist=sqrt(loca(1,:).^2+loca(2,:).^2+height^2);
    ele=atan(height./sqrt(loca(1,:).^2+loca(2,:).^2));
    theta=2*pi*antd/lambda.*(loca(2,:)./sqrt(loca(1,:).^2+loca(2,:).^2+height^2));
    K=0.8888*exp(0.03272*(ele./pi.*180));  % K-factor  equation obtained by the curve fitting of the results in the paper
    
    H=zeros(N,M);
    alpha=1./(4*pi*dist*10^3/lambda).^2;  % large scale fading
    rho=sqrt(K./(1+K));
    A=exp(sqrt(-1).*kron(theta,[1:N]'));  % N-by-M
    B=randn(N,M);
    for k=1:M
        H(:,k)=sqrt(alpha(k))*(rho(k).*A(:,k)+sqrt(1-rho(k)^2).*B(:,k));
    end
    pro_sr_G=zeros(length(tau_seq),length(G_seq));
    ran_sr_G=zeros(1,length(G_seq));
    dist_sr_G=zeros(length(tau_seq),length(G_seq));
    
    for i=1:length(G_seq)
        
        G=G_seq(i);
        
        % proposed method
        %---------------------------------------------------------------------------
        ave_sr=zeros(1,length(tau_seq));
        
        for r=1:runs_out
            
            pro_times=zeros(M,length(G_seq),length(tau_seq));
            
            t1=clock;
            [X,WW]=usergrouping(rho,A,G);  % divide users into G groups   X: G-by-M  WW: N-by-G
            t2=clock;
            ave_pro_cput(i)=ave_pro_cput(i)+etime(t2,t1);
            
            for tau_ite=1:length(tau_seq)
                
                tau=tau_seq(tau_ite);
            
                % user selection
                sumrate=zeros(1,runs_in);%runs_in是多次取均值
                for j=1:runs_in
                    W=zeros(N,G);
                    for g=1:G
                        w=tau.*WW(:,g)+(1-tau^2).*randn(N,1);
                        W(:,g)=sqrt(P/G)/norm(w).*w;
                    end
                    sinr=zeros(1,M);
                    for k=1:M
                        g=find(X(:,k));
                        sinr(k)=norm(H(:,k)'*W(:,g))^2/(norm(H(:,k)'*W)^2-norm(H(:,k)'*W(:,g))^2+sigma^2);
                    end
                    select=zeros(1,G);
                    rate=zeros(1,G);
                    for g=1:G
                        [val,ind]=max(sinr.*X(g,:));
                        select(g)=ind;
                        pro_times(ind,i,tau_ite)=pro_times(ind,i,tau_ite)+1;
                        rate(g)=log2(1+val);
                    end
                    sumrate(j)=sum(rate);
                end
                ave_sr(tau_ite)=ave_sr(tau_ite)+mean(sumrate);
            
            end
            
            ave_pro_times=sum((pro_times>0),1);

            ave_pro_num(:,i)=ave_pro_num(:,i)+reshape(ave_pro_times(:,i,:)./G_seq(i),length(tau_seq),1);
            
        end
        ave_sr=ave_sr./runs_out;
        
        pro_sr_G(:,i)=ave_sr';
        %---------------------------------------------------------------------------
        % random user selection
        %---------------------------------------------------------------------------
        ran_sr=0;
        for j=1:runs
            ran_ind=randi(M,1,G);
            ran_W=sqrt(P/G)/sqrt(N).*A(:,ran_ind);
            ran_sinr=zeros(1,G);
            for g=1:G
                ran_sinr(g)=norm(H(:,ran_ind(g))'*ran_W(:,g))^2/(norm(H(:,ran_ind(g))'*ran_W)^2-norm(H(:,ran_ind(g))'*ran_W(:,g))^2+sigma^2);
            end
            ran_rate=log2(1+ran_sinr);
            ran_sr=ran_sr+sum(ran_rate);
        end
        ran_sr=ran_sr/runs;

        ran_sr_G(i)=ran_sr;
        %---------------------------------------------------------------------------
        % distance based method
        %---------------------------------------------------------------------------
        ave_sr=zeros(1,length(tau_seq));
        
        for r=1:runs_out
            
            dist_times=zeros(M,length(G_seq),length(tau_seq));
            
            t1=clock;
            [X,WW]=usergrouping_dist(alpha.*10^12,rho,A,G);  % divide users into G groups   X: G-by-M  WW: N-by-G
            t2=clock;
            ave_dist_cput(i)=ave_dist_cput(i)+etime(t2,t1);
            
            for tau_ite=1:length(tau_seq)
                
                tau=tau_seq(tau_ite);
            
                % user selection
                sumrate=zeros(1,runs_in);
                for j=1:runs_in
                    W=zeros(N,G);
                    for g=1:G
                        w=tau.*WW(:,g)+(1-tau^2).*randn(N,1);
                        W(:,g)=sqrt(P/G)/norm(w).*w;
                    end
                    sinr=zeros(1,M);
                    for k=1:M
                        g=find(X(:,k));
                        sinr(k)=norm(H(:,k)'*W(:,g))^2/(norm(H(:,k)'*W)^2-norm(H(:,k)'*W(:,g))^2+sigma^2);
                    end
                    select=zeros(1,G);
                    rate=zeros(1,G);
                    for g=1:G
                        [val,ind]=max(sinr.*X(g,:));
                        select(g)=ind;
                        dist_times(ind,i,tau_ite)=dist_times(ind,i,tau_ite)+1;
                        rate(g)=log2(1+val);
                    end
                    sumrate(j)=sum(rate);
                end
                ave_sr(tau_ite)=ave_sr(tau_ite)+mean(sumrate);
            
            end
            
            ave_dist_times=sum((dist_times>0),1);

            ave_dist_num(:,i)=ave_dist_num(:,i)+reshape(ave_dist_times(:,i,:)./G_seq(i),length(tau_seq),1);
            
        end
        ave_sr=ave_sr./runs_out;
        
        dist_sr_G(:,i)=ave_sr';
        %---------------------------------------------------------------------------
        
    end
    
    ave_pro_sr=ave_pro_sr+pro_sr_G;
    ave_ran_sr=ave_ran_sr+ran_sr_G;
    ave_dist_sr=ave_dist_sr+dist_sr_G;

end

ave_pro_sr=ave_pro_sr./count;
ave_ran_sr=ave_ran_sr./count;
ave_dist_sr=ave_dist_sr./count;
ave_pro_num=ave_pro_num./(count*runs_out);
ave_dist_num=ave_dist_num./(count*runs_out);
ave_pro_cput=ave_pro_cput./(count*runs_out);
ave_dist_cput=ave_dist_cput./(count*runs_out);







