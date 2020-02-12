
function [w,eta]=solve_w(a,b)

[N,Mg]=size(a);

[r,res]=mosekopt('symbcon echo(0)');
% [r,res]=mosekopt('symbcon');

% specify the non-conic part of the problem
prob=[];
num_var=2*N+1+1;
num_lcons=Mg;
la=zeros(num_lcons,num_var);
lc=zeros(num_lcons,1);
uc=zeros(num_lcons,1);

prob.c=[zeros(1,num_var-1) -1]';

num_c=0;
% real(a'w+b)-eta.*ones(Mg,1)>=0
A=[real(a') -imag(a')];
la(num_c+1:num_c+Mg,:)=[A zeros(Mg,1) -ones(Mg,1)];
lc(num_c+1:num_c+Mg)=-b';
uc(num_c+1:num_c+Mg)=inf.*ones(Mg,1);
num_c=num_c+Mg;

% % imag(a'w)
% la(num_c+1:num_c+Mg,:)=[imag(a') real(a') zeros(Mg,2)];
% num_c=num_c+Mg;

% variable cons
lx=-inf*ones(num_var,1);
ux=inf*ones(num_var,1);

% q=sqrt(N)
lx(num_var-1)=sqrt(N);
ux(num_var-1)=sqrt(N);

% % eta>=0
% lx(num_var)=0;

prob.a=sparse(la);
prob.blc=lc';
prob.buc=uc';
prob.blx=lx';
prob.bux=ux';


% set up the cone information
tp1=res.symbcon.MSK_CT_QUAD;
% tp2=res.symbcon.MSK_CT_RQUAD;
prob.cones.type=tp1;

sub1=zeros(1,num_var-1);
sub2=1;
% norm(w)<=q
sub1=[num_var-1 1:2*N];

prob.cones.sub=sub1;
prob.cones.subptr=sub2;

[r,res]=mosekopt('minimize echo(0)',prob);
% [r,res]=mosekopt('minimize',prob);

if strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
    % obtain the solution
    xsol=res.sol.itr.xx;
    w=xsol(1:N)+sqrt(-1).*xsol(N+1:2*N);
    eta=xsol(num_var);
else
    cvx_begin quiet
    variable w(N) complex
    variable eta
    maximize eta
    subject to
    norm(w) <= sqrt(N);
    for i=1:Mg
        real(a(:,i)'*w)+b(i) >= eta;
    end
    cvx_end
end


