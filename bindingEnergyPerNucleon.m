%script using the Semi-Empirical Mass formula to plot predicted binding energy vs
%number of nucleons, showing that the nuclei are less bound per nucleon for
%the intermediate values of A(number of nucleons)
%This can be explained by noting that the larger nucleii will have a high
%BE/A value making them likely to undergo fission, while the smaller ones
%are similarly more likely to undergo fusion

%In reality  56Fe should have the maximum binding energy, so this gives a
%rough approximation but not an exact prediction

x=(1:240);
y=zeros(1,240);
q=15.8;
w=17.8;
e=23.7;
r=.711;


for  i=1:1:240
    N=x(i);
    f= 4*e*N/(8*e*N+2*r*(N^(5/3)));
    BE = q*N-w*(N^(2/3))-(r*f^2*N^2)/(N^(1/3))-(e*(2*f*N-N)^2)/N;
    y(i)=BE/N;
end
plot(x,y)
xlabel('Number of Nucleons (A)')
ylabel('Binding Energy per Nucleon in MeV')