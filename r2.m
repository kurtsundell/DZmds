function [r2] = r2(x, y)

x = x;
y = y;

lengthx = length(x);

xmean = mean(x);
ymean = mean(y);

xcov = zeros(1,length(x));
ycov = zeros(1,length(y));

for i = 1:lengthx;
    xcov(i) = x(i) - xmean;
end

for i = 1:lengthx;
    ycov(i) = y(i) - ymean;
end

xcov = transpose (xcov);
ycov = transpose (ycov);

mult = xcov.*ycov;
numerator = sum(mult);

xcov2 = xcov.*xcov;
sumxcov2 = sum(xcov2);
ycov2 = ycov.*ycov;
sumycov2 = sum(ycov2);
mult2 = sumxcov2*sumycov2;
denominator = sqrt(mult2);

r = numerator/denominator;

r2 = r*r;