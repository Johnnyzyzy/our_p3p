function y= xrand(m,n,r)
%function description:
%generate random numbers ranging from r(1) to r(2);
y= r(1)+rand(m,n)*(r(2)-r(1));

return