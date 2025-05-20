function [devia] = fitPhase(driftFit,logKSpace,drift)
% better use rem of imag+pi
devia=logKSpace-driftFit*drift;
devia=real(devia)+1i*(rem(imag(devia)+pi,2*pi)-pi);
devia=sum(sum(abs(devia).^2));
end