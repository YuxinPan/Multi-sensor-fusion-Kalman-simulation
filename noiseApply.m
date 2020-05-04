% as main function in the file, locate at top of the file
function r = noiseApply(coordinate,scaleUniform,scaleNormal)

    r(1) = coordinate(1)+noiseGen(scaleUniform,scaleNormal);
    r(2) = coordinate(2)+noiseGen(scaleUniform,scaleNormal);
    r(3) = coordinate(3)+noiseGen(scaleUniform,scaleNormal);
    
end

% add both uniform and normal noise
function r = noiseGen(scaleUniform,scaleNormal)

    uniformNoise = rand()-0.5; % uniform distribution -0.5 to 0.5
    normalNoise = randn(); % normal distribution with mean 0 and standard deviation 1 
    r = scaleUniform*uniformNoise + scaleNormal*normalNoise;
    
end
