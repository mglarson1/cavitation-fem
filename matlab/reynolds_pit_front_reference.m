function row=reynolds_pit_front_reference(right)
% Independent 1D quadrature: H^3 P'=H-H(xc), P(xc)=P'(xc)=P(right)=0.
H=@(x)1+exp(-(x-12).^2);
F=@(xc)integral(@(x)(H(x)-H(xc))./H(x).^3,xc,right,'AbsTol',1e-13,'RelTol',1e-12);
xc=fzero(F,[8,12],optimset('TolX',1e-12));xpeak=24-xc;
ppeak=6*integral(@(x)(H(x)-H(xc))./H(x).^3,xc,xpeak,'AbsTol',1e-13,'RelTol',1e-12);
row=table(right,xc,xc-12,xpeak,ppeak,F(xc),'VariableNames', ...
 {'right','front_x','front_relative_to_pit','peak_x','peak_pressure','closure_residual'});
end
