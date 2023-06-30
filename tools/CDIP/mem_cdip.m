function [s,chk] = mem(a1,a2,b1,b2)
    begin = 1;
    ndeg = 360;
    res = 1;
    [s,chk] = memd(a1,a2,b1,b2,begin,ndeg,res);
    
