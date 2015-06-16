function [ value ] = GetTimeValue( n,size,time, array )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
for i = 1:size
    if time(i) >= n
        value = array(i);
        break;
    end
end;
end

