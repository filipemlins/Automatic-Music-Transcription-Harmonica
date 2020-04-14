function [ output_frequency ] = correct_frequency( input_frequency )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
   n1 = floor(12*log(input_frequency/442)/log(2)); 
  output_frequency = 442*2^(n1/12);
end

