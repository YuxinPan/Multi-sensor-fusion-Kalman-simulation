% generate impulse/disturbance % distinguish disturbance and noise
function impulse = impulseGen(scanRate)
  duration =  rand()+2/scanRate; % in unit scans
  t = (0:1/scanRate:duration);
   
  unitstep = t>=0;

  a = randn();
  b = randn();
  c = randn();

  poly = a* (t).^4.*unitstep-b*(t).^2.*unitstep-c*(t).*unitstep;

  polyN = a* (-t+duration).^4.*unitstep-b*(-t+duration).^2.*unitstep-c*(-t+duration).*unitstep;

  impulse = [poly(1:floor(duration*scanRate/2)) polyN(floor(duration*scanRate/2):floor(duration*scanRate))];
  impulse = impulse';
 end