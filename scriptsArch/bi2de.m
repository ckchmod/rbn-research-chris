function b = bi2de(vec) 

b= 0;
for i = 1:length(vec)
   b = b + vec(i)*2^(length(vec)-i);
end