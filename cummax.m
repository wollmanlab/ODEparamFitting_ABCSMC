function cmx = cummax(A)

cmx=nan(size(A)); 
for i=1:size(A,1)
    cmx(i,:)=max(A(1:i,:)); 
end