function [s,chk] = memd(a1,a2,b1,b2,begin,ndeg,res)

dr = 0.0174533;

% for i = begin : 1 : begin+ndeg-1
%    if i > 360
%      dir = i - 360;
%    else
%      dir = i;
%    end
%    s(dir)=0;
% end
s = zeros(1,360);

% switch to Lygre & Krogstad notation
d1 = a1;
d2 = b1;
d3 = a2;
d4 = b2;

c1 = d1 + 1i*d2;
c2 = d3 + 1i*d4;

p1 = (c1 - c2*conj(c1)) / (1-abs(c1).^2);
p2 =  c2 - c1*p1;
x = 1 - p1*conj(c1) - p2*conj(c2);

% sum over 'ndeg' in steps, get distribution with 'res' degree resolution
tot = 0;
offset = 0.5 * (1.0 - res);
for rn = begin-offset : res : begin+ndeg-1+offset
    a = rn*dr;
    e1 = cos(1*a) - 1i*sin(1*a);
    e2 = cos(2*a) - 1i*sin(2*a);
    y = abs(1 - p1*e1 - p2*e2)^2;

    % put in proper 1 deg directional band
    ndir = round(rn);
    if ndir > 360 
        ndir = ndir-360;
    end
%     ndir = mod(ndir,360);

    % normalize by 360/(step size) if not running full 360 degrees
    if ndeg ~= 360
      s(ndir) = s(ndir) + abs(x/y)/(360./res);
    else
      s(ndir)=s(ndir)+abs(x/y);
    end
    tot = tot+abs(x/y);
end

% normalize spectrum for full 360 degree run
if ndeg ~= 360
    for i = 1 : 1 : 360
        s(i) = s(i)/tot;
    end 
    chk = 1;
else
% tot should = 360.  If directional peak is extremely narrow then
% 1 deg resolution may be insufficient and tot ~= 360
    chk = tot/(360./res);
end
