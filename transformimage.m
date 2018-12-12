function J=transformimage(I, H, s)

Jh=s(1); Jw=s(2);

[X2,Y2]=meshgrid(1:Jw,1:Jh) ;

invH = inv(H);

X=zeros(Jh,Jw);

Y=zeros(Jh,Jw);

for n=1:numel(X2)

P2=[X2(n); Y2(n); 1];

P = invH*P2;

X(n)=P(1)/P(3); Y(n)=P(2)/P(3);

end

J=zeros(Jh,Jw,size(I,3));

for k=1:size(I,3)

J(:,:,k) = interp2(I(:,:,k), X,Y, 'linear',0);

end