
%Script to recursively construct and plot the first 6 "Cardinal B Splines"
%and the corresponding wavelets, and the normalized refineable functions
%For the definition of Cardinal B Splines see:  https://en.wikipedia.org/wiki/B-spline#Cardinal_B-spline
% For an explanation of the normalized refineable function, see https://www.math.cuhk.edu.hk/~djfeng/fengpapers/ACHA/refinable.pdf
close all

%set y-values by calling function explicitly defining the 3rd Cardinal
%B-Spline (same as defined on Wikipedia)
x=0:0.01:4;
y=arrayfun(@B3,x);
plot(x,y,'linewidth',2)
title('B3(x)')

%Again use explicitly defined formula for 4th Cardinal B-Spline
figure
b4=@(x)NthRecursiveSpline(4,x);
y4=arrayfun(b4,x);
plot(x,y4,'linewidth',2)
title('B4(x)')

%Instead we switch to our recursively constructed version and plot the
%different-degree splines on the same plot, and then add the corresponding "wavelet"
figure
hold on
b3=@(x)NthRecursiveSpline(3,x);
y3=arrayfun(b3,x);
plot(x,y3)
b4=@(x)NthRecursiveSpline(4,x);
y4=arrayfun(b4,x);
plot(x,y4)
hold on
%We also add the "refinable" function
plotRefinableFunction(1)
plotRefinableFunction(2)
plotRefinableFunction(3)
plotRefinableFunction(4)
legend('Spline3','Spline4','f1','f2','f3','f4')
figure;hold on;
plotWavelet(2)
plotWavelet(3)
plotWavelet(4)
plotWavelet(5)
legend('s2','s3','s4','s5')

%function to find the values of mth wavelet at each point and plot it 
function plotWavelet(m)
x=0:0.01:9;
pt=find(x==((2*m-1)/2));
y=zeros(1,length(x));
for i = 1:length(x)
    y(i)=wavelet(m,x(i));
end
plot(x,y)
plot(x(pt),y(pt),'*','linewidth',3,'HandleVisibility','off')
end
%recursive construction of mth wavelet, for each x value
function a = wavelet(m,x)
 a=0;
for k = 1:3*m-2+1

    sum = 0;
    for l = 1:m+1
        sum=sum+Bk(2*m,k-l+1);
    end
    q=(((-1)^(k-1))/(2^(m-1)))*sum;
    a=a+q*NthRecursiveSpline(m,2*x-k);
    
end
end
%Recursive defenition of Cardinal B-Splines with two base cases
function a = NthRecursiveSpline(m,x)
if m==1
    a=B1(x);
end
if m==2
   a=B2(x); 
end
if m>2
    a=(x/(m-1))*NthRecursiveSpline(m-1,x)+((m-x)/(m-1))*NthRecursiveSpline(m-1,x-1);
end

end
%Helper function for calculating the wavelet, this is the "B2" term in the
%formula
function a = Bk(m,k)

if m==2
    a=B2k(k);
else
    a=(k/(m-1))*Bk(m-1,k)+((m-k)/(m-1))*Bk(m-1,k-1);
    
end
end
%base case for calculating the "b2" term in the wavelet formula
function a = B2k(k)
a=0;
if k==1
a=1;
end
end
function plotRefinableFunction(m)
x=0:0.01:4;
y=zeros(1,length(x));
for i = 1:length(x)
    y(i)=0;
    for j=1:m+1
        cj=(1/(2^m-1))*choose(m,j);
        y(i)=y(i)+cj*NthRecursiveSpline(m,2*x(i)-j);
    end
end
plot(x,y,'linewidth',2)
end
function a = B1(x)
a=0;
if x>=0 && x<1
  a=1;
end
end

%explicit definition of 2nd CArdinal B-spline
function a = B2(x)
a=0;
if x>=0 && x<1
  a=x;
end
if x>=1 && x<2
    a=2-x;
end
end
%explicit definition of 3rd Cardinal B-spline
function a = B3(x)
a=0;
if x>=0 && x<1
  a=(x^2)/2;
end
if x>=1 && x<2
    a=(1/2)*(-2*(x^2)+6*x-3);
end
if x>=2 && x<3
    a=(1/2)*(x-3)^2;
end
end
%explicit definition of 4th Cardinal B-spline
function a = B4(x)
a=0;
if x>=0 && x<1
  a=(x^3)/6;
end
if x>=1 && x<2
    a=-(x^4)/18+(x^2)/2-(5*x)/18;
end
if x>=2 && x<3
    a=(1/6)*(3*x^3-24*x^2+60*x-44);
end
if x>=3 && x<4
    a=-((x^3)/6) + 2*(x^2)-(8*x)+(32/3);
end
end
% easy helper function for the refineable function
function a = choose(n,k)
if k<=n
     a=nchoosek(n,k);
else
    a=0;
end
end

