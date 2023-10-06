function [z,lambda]=psplineBspline_mod(x,y,m)
%
% Definizione dei nodi delle B-splines
% m (number of breakpoints) <<n (number of knots)
% Il primo nodo coincidente con il primo dato
% l' ultimo   nodo coincidente con l'ultimo dato
%
% Per creare basi considero nodi equispaziati (Bspline stessa altezza)
knots=linspace(x(1),x(end),m-6); % perchè escludo i sei nodi esterni

%VETTORE DEI NODI AMPLIATO (i miei xi vedi presentazione)
k=m; %aggiungo 3 nodi a sx e 3 a dx per completare le basi considerando di
      %generare cubic B-splines (ordine 4 e passano per 5 nodi)

%distanza tra i primi due nodi per costruire h ampiezza
meshsizes=knots(2)-knots(1);

x0=min(knots)-meshsizes;
xNp1=max(knots)+meshsizes;
%
% Costruzione dei nodi esterni
% i 3 a sx equispaziati
z(1)=x0-2*meshsizes;
z(2)=x0-meshsizes;
z(3)=x0;
% nodi interni coincidono con quelli di partenza
z(4:k-3)=knots;
% i 3 a dx
z(k-2)=xNp1;
z(k-1)=xNp1+meshsizes;
z(k)=xNp1+2*meshsizes;
%questo vettore contiene tutti i nodi necessari per costruire le B-splines
break_point=z; %  

x=x(:);
y=y(:);
n=length(x) ;

%Matrice pesi coincide con matrice identica
W=eye(n,n);

% Costruzione matrice (B^T B+lambda D^T D)

%1)	costruzione della matrice dei valori delle B-spline Bj 
% costruite sui nodi breakpoints e valutate nei dati x

for i=1:n
    for j=1:k-4 %perchè l'ultima b-s verrà costruita nel penultimo brpoint
        if  x(i)>break_point(j)&&  x(i)< break_point(j+4)
            %Cubic B-spline
            tmpb=bspline(break_point(j:j+4));
            B(i,j)= fnval(tmpb, x(i));
        end
    end

end

%2)	costruzione della matrice D della discretizzazione di L2

D=zeros(m-6,m-4);
for i=1:m-6
    for j=1:m-4
        if i==j
            D(i,j)=1;
        else if j==i+1
                D(i,j)=-2;
        else if j==i+2
                D(i,j)=1;
        end
        end
        end
    end
end


%% calcolo di lambda

% Blambda=B;
% Dlambda=D;
% [U,sm] =cgsvd(Blambda,Dlambda);
% lambda = l_curve(U,sm,y);
lambda=1;

A=(B'*W*B+lambda^2*D'*D);

%risoluzione di  A*c=B'*y
y1=W*y;
y2=B'*y1;
c=A\y2;

x_ampl=break_point;
mc=length(break_point) ;
n=mc-6;
xx=x;
for k=1:length(xx)
    yy(k)=0;
    for i=1:mc-4
        if xx(k)>x_ampl(i)&& xx(k)< x_ampl(i+4)
            tmpb=bspline(x_ampl(i:i+4));
            yy(k)=yy(k)+c(i)*fnval(tmpb,xx(k));
        end
    end
end

z=yy;
end

