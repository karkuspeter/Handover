function q = mpi2pi(q)
% transforms the staning qin values to -pi to pi

while any(any(q > pi))
    i=find(q>pi); 
    q(i)=q(i)-2*pi; 
end    

while any(any(q < -pi))
    i=find(q < -pi); 
    q(i)=q(i)+2*pi; 
end