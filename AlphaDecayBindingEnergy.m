%Alpha decay:
%Plotting the difference in binding energy from losing one nucleon
m=300
x=(1:m)
y=zeros(1,m);
yd=zeros(1,m);
q=15.8;
w=17.8;
e=23.7;
r=.711;


for  i=1:1:m
    N=x(i);
    f= 4*e*N/(8*e*N+2*r*(N^(5/3)))
    BE = q*N-w*(N^(2/3))-(r*f^2*N^2)/(N^(1/3))-(e*(2*f*N-N)^2)/N;
    y(i)=BE;
end
for  i=4:1:m
    N=x(i)-4;
    f= 4*e*N/(8*e*N+2*r*(N^(5/3)))
    BED = q*N-w*(N^(2/3))-(r*f^2*N^2)/(N^(1/3))-(e*(2*f*N-N)^2)/N;
    yd(i)=BED;
end
plot(x,y-yd-28.3)
xlabel('Number of Nucleons (A)')
ylabel('Difference in Binding Energy Minus Single a-particle binding energy')