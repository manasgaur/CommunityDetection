function [cutvalue,evalue]=ncuts(A)

r=size(A,1);
D=zeros(r,r);

for i=1:r
    A(i,i)=0;
end

for i=1:r
    sumA=sum(A(i,:));
    D(i,i)=sumA;
end

Q0=0;Q1=1;
Y1=ones(r,1);
itr=1;
qarr=zeros(1000,1);
while norm((Q1-Q0),'fro') > 0.005
    Q0=Q1;
    Q1=Y1'*(D-A)*Y1;
    tmp=(D-A)*Y1;
    Y1=tmp/norm(tmp,'fro');
    qarr(itr)=Q0;
    itr=itr+1;
    %disp(Q0);
end
figure, hold on
loglog(1:itr,qarr(1:itr,1),'-.r*');
hold off
evalue=Y1'*(D-A)*Y1;
cutvalue=Y1(r);
end