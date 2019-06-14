

close all; figure
fun =@(x)exp(3*x)*sin(200*x^2)/(1+20*x^2);
computefromL(fun,7)
function computefromL(fun, n)
L=zeros(n-1);
for i = 1:n-1
    for j = 1:n-1
        if i==j && i==1
            L(i,j)=2;
        end
        if i==j+1 && i~=1
            L(i,j)=1/L(j,j);
        end
        if i==j && i>1
            L(i,j)=sqrt(4-L(j,i-1)^2);
        end
    end
end
T=L*L';
x=linspace(0,1,n+1)';
h=x(2)-x(1);
y=arrayfun(fun,x);
der = matlabFunction(diff(sym(fun)));
D=arrayfun(der,x);
theta = zeros(n-1,1);
for i = 2:n

    d1= (y(i+1)-y(i))/(x(i+1)-x(i));
    d2= (y(i)-y(i-1))/(x(i)-x(i-1));
    theta(i-1)=(3/h)*(d1-d2);
    
end
%find c using product matrix, then same thing using lower triangular and
%its transpose to check they're roughly equivalent
c=(L'*L)\theta;
c=L'\(L\theta)
c=vertcat(0,vertcat(c,0))

b=zeros(n+1,1);c=zeros(n+1,1);d=zeros(n,1);
for i = 1:n
   a=y(i);
%    b=D(i);
%    c = 3*(y(i+1)-y(i)) - 2*D(i)-D(i+1) 
%    d = 2*(y(i)-y(i+1))+D(i)+D(i+1)
   d=(c(i+1)-c(i))/(3*h);
   b=(y(i+1)-y(i))/(x(i+1)-x(i))-(h/3)*(2*c(i)+c(i+1));
   cube =@(X) a + b*(X-x(i))+c(i)*(X-x(i))^2+d*(X-x(i))^3;
   cutx=linspace((i-1)*h,h*i,100);
   cuty=arrayfun(cube, cutx);
   plot(cutx,cuty,'linewidth',2)
   hold on
   
end
hold on

realxpts=linspace(0,1,1000);
realypts=arrayfun(fun,realxpts);
plot(realxpts,realypts)


end

