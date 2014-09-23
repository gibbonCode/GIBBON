function opt=fmin_fvalfind_ppform(pp,x,y)

yf=ppval(pp,x);
opt=abs(yf-y)^2;

end