load spec;
ws = spec(:,7);
vxc = -10.455963;
elda = -5.956967;
#vxc = -11.249807;
#elda= 6.009187;

ReSig = spec(:,8)-vxc;
ImSig = spec(:,9);
Sigma = ReSig + I*ImSig;

A = 1/pi * abs(ImSig) ./ ( (ws-elda-ReSig).^2 + ImSig.^2);

Gamma = abs( ImSig )/pi;
dw = ws(2)-ws(1);

if( min(ws) < elda && elda < max(ws))
  IndElda = (elda - min(ws))/dw + 1;
  distA   = IndElda - floor(IndElda);
  distB   = ceil(IndElda) - IndElda;
  GammaElda = distB*Gamma(floor(IndElda)) + distA*Gamma(ceil(IndElda));
  DGammaA = Gamma( floor(IndElda)+1) - Gamma( floor(IndElda) -1);
  DGammaA /= 2*dw;
  DGammaB = Gamma( ceil(IndElda)+1) - Gamma( ceil(IndElda) -1);
  DGammaB /= 2*dw;
  DGammaElda = distB*DGammaA + distA*DGammaB;

  deltaE  = distB*ReSig(floor(IndElda)) + distA*ReSig(ceil(IndElda));
  etak = pi*GammaElda;

  DSigmaA = Sigma( floor(IndElda)+1) - Sigma( floor(IndElda) -1);
  DSigmaA /= 2*dw;
  DSigmaB = Sigma( ceil(IndElda)+1) - Sigma( ceil(IndElda) -1);
  DSigmaB /= 2*dw;
  DSigmaElda = distB*DSigmaA + distA*DSigmaB;

  alphak  =  imag(DSigmaElda);
  gammak  = -real(DSigmaElda);
else
  printf("error - bad freq. range!");
  break;
endif;

ts   = [-10 : 0.01 : 10]';
#ts   = [-60: 0.01 : 60]';
dt   = ts(2)-ts(1);
wt   = ws*ts';

%# calculate qp spectral function following Hedin/Almbladh
Cqp   = -I*ts*deltaE - etak*abs(ts) + I*alphak*sign(ts) - gammak;
expCqp= exp(-I*elda*ts) .* exp(Cqp);
exptw = exp(I*ts*ws');

%# compute full spectral function following Gunnarsson/Aryasetiawan
%# have to make choice for chemical potential - middle of gap?
eF = 6.009187 + 0.6;
CS = Gamma.*( (ws-elda) < eF )-GammaElda-(ws-elda)*DGammaElda;
%# get rid of divergence at ws=elda
minInd = floor(IndElda) - 10;
maxInd = floor(IndElda) + 10;
P= polyfit((ws-elda)(minInd:maxInd),CS(minInd:maxInd),2);
PN  = [P(1) 0 0];
CSZ = polyval(PN,ws-elda);
#wsZ = (ws-elda)(minInd:maxInd);
#CSF = CS(minInd:maxInd);
#a   = sum( CSF.*wsZ.^2 ) / sum(wsZ.^4);
#CSZ2= a*(ws-elda).^2;

wCut= 1.5;
fCut= exp(-(ws-elda).^2/wCut^2);
CSN = CS  .* (1-fCut) + CSZ.*fCut;
CS2 = CSN ./ (ws-elda).^2;

%# Fourier transform to real time
CSt = dw*sum( dmult(CS2, exp(-I*wt)), 1);
CSt = transpose(CSt);
CSt = exp(I*elda*ts) .* CSt;
expC= exp(CSt) .* expCqp;

%# Fourier transform back to freq. space
wsN = [-80:0.1:80]';
exptwN = exp(I*ts*wsN');

%# quasiparticle spectral function
Aqp  = sum( dmult(expCqp, exptwN ),1);
Aqp  = transpose(Aqp);
Aqp *= dt/2/pi;

%# full spectral function
Agun= sum( dmult(expC, exptwN ),1);
Agun= transpose(Agun);
Agun *= dt/2/pi;

%# first order expansion
dAS = dt* sum( dmult( expCqp .* CSt, exptwN), 1);
dAS = transpose(dAS);
dAS /= 2*pi;
