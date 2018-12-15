function num = digits(in)
    num=0;
    if max(size(in)) == 1
        num = double(str2num(in)>5);
    else
        num = double(str2num(in(1))>5) + digits(in(2:max(size(in))));
end
end