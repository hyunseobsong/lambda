function tblOut = lambda_v0p31_crowdsourced(tbl,phspan,wrt,dataDescrp)

nCpd=size(tbl,1); % # of compounds

%% add 'Z' column (electron charge) if not provided from the input data
% -- make a reduced-size table

MolForm=tbl.MolForm;
C=tbl.C;
H=tbl.H;
N=tbl.N;
O=tbl.O;
P=tbl.P;
S=tbl.S;

%-- check if Z column exists and otherwise add
idx=find(strcmp('Z',tbl.Properties.VariableNames));
if ~isempty(idx) % 'Z' column exists
    Z=tbl.Z;
else % 'Z' column does not exist
    Z=zeros(nCpd,1);
end

tblCpd=table(MolForm,C,H,N,O,P,S,Z);

%% initialize output variables 
%-- thermodynamic functions 
delGcox=zeros(nCpd,1); % electron donor half rxn, kJ/C-mol
delGd=zeros(nCpd,1); % electron donor half rxn, kJ/mol
delGcat=zeros(nCpd,1); % catabolic rxn, kJ/mol
delGan=zeros(nCpd,1); % anabolic rxn, kJ/mol
delGdis=zeros(nCpd,1); % dissipation energy, kJ/mol
lambda=zeros(nCpd,1); % dimensionless 
%-- CUE and stoichiometric coeffs in rxns 
CUE=zeros(nCpd,1); % carbon use efficiency (=C-mole biomass yield) at pH=7
NUE=zeros(nCpd,1); % nitrogen use efficiency (=C-mole biomass yield) at pH=7
TER=zeros(nCpd,1); % threhold element ratio (=C-mole biomass yield) at pH=7
stoichD=zeros(nCpd,10); % electron donor half rxn
stoichA=zeros(nCpd,10); % electron acceptor half rxn
stoichCat=zeros(nCpd,10); % catabolic rxn
stoichAn=zeros(nCpd,10); % anabolic rxn
stoichMet=zeros(nCpd,10); % metabolic rxn

%% get outpout table
% phspan=0:1:14;
% phspan=7;
nph=length(phspan);
CUEph=zeros(nCpd,nph);
NUEph=zeros(nCpd,nph);
TERph=zeros(nCpd,nph);
collate = cell(nph,2);
for iph=1:nph
    ph=phspan(iph);
    for iCpd=1:nCpd
        abcdefz=[tblCpd.C(iCpd),tblCpd.H(iCpd),tblCpd.N(iCpd),tblCpd.O(iCpd),tblCpd.P(iCpd),tblCpd.S(iCpd),tblCpd.Z(iCpd)];
        [delGcox_,delGd_,delGcat_,delGan_,delGdis_,lambda_,...
            CUE_,NUE_,TER_,stoichD_,stoichA_,stoichCat_,stoichAn_,stoichMet_] = ...
            getThermoStoich(abcdefz,ph);

        % store results for individual compounds
        delGcox(iCpd)=delGcox_;
        delGd(iCpd)=delGd_;
        delGcat(iCpd)=delGcat_;
        delGan(iCpd)=delGan_;
        delGdis(iCpd)=delGdis_;
        lambda(iCpd)=lambda_;
        CUE(iCpd)=CUE_;
        NUE(iCpd)=NUE_;
        TER(iCpd)=TER_;
        stoichD(iCpd,:)=stoichD_;
        stoichA(iCpd,:)=stoichA_;
        stoichCat(iCpd,:)=stoichCat_;
        stoichAn(iCpd,:)=stoichAn_;
        stoichMet(iCpd,:)=stoichMet_;
    end    
    
    CUEph(:,iph)=CUE;
    NUEph(:,iph)=NUE;
    TERph(:,iph)=TER;

    % store the results as tables
    tblThermo=table(delGcox,delGd,delGcat,delGan,delGdis,lambda);
    tblStoich=table(CUE,NUE,TER,stoichD,stoichA,stoichCat,stoichAn,stoichMet); %[e-donor,h2o,hco3,nh4,hpo4,hs,h,e,e-acceptor,biom]
    tblOut=[tblCpd tblThermo tblStoich];

    % Fix negative lambda
    idxZero = find(tblOut.lambda<0);
    tblOut.lambda(idxZero) = 0;

    % save the results
    if strcmpi(wrt,'y')
        save(sprintf("%s_out_pH%s.mat",dataDescrp,num2str(ph)))
        write(tblOut,sprintf("%s_out_pH%s.csv",dataDescrp,num2str(ph)))
    end

    % collate the results
    collate{iph,1} = ph;
    collate{iph,2} = tblOut;

end

% function output
tblOut = cell2table(collate,"VariableNames",["pH","tblOut"]);

% figure
% x=phspan;
% y=CUEph';
% 
% plot(x,y(:,1:100),'-o')
% set(gca,'fontsize',14,'linewidth',1.5)
% xlabel('pH')
% ylabel('CUE')


function [delGcox,delGd,delGcat,delGan,delGdis,lambda,...
        CUE,NUE,TER,stoichD,stoichA,stoichCat,stoichAn,stoichMet] = ...
        getThermoStoich(abcdefz,pH)

a=abcdefz(1);
b=abcdefz(2);
c=abcdefz(3);
d=abcdefz(4);
e=abcdefz(5);
f=abcdefz(6);
z=abcdefz(7);

% Step 1a) stoichD: stoichiometries for an electron donor
ySource=-1;
yH2o=-(3*a+4*e-d);
yHco3=a;
yNh4=c;
yHpo4=e;
yHs=f;
yH=5*a+b-4*c-2*d+7*e-f;
yE=-z+4*a+b-3*c-2*d+5*e-2*f;
stoichD=[ySource,yH2o,yHco3,yNh4,yHpo4,yHs,yH,yE];
stoichD([9,10])=[0,0]; % add additional components: e-acceptor and biomass

% Step 1b) stoichA: stoichiometries for an electron acceptor (i.e., oxygen)
stoichA=zeros(1,10);
stoichA(9)=-1; % oxygen
stoichA(7)=-4; % h+
stoichA(8)=-4; % e-
stoichA(2)=2; % h2o

% Step 1c) stoichCat: stoichiometries for catabolic reaciton 
yEd=stoichD(8);
yEa=stoichA(8);
stoichCat=stoichD-(yEd/yEa)*stoichA;

% Step 2a) stoichAnStar: stoichiometries for anabolic reaciton 
%          (N source = NH4+)

abcdefzBiom=[1 1.8 0.2 0.5 0 0 0]; % C H_1.8 N_0.2 O_0.5
aB=abcdefzBiom(1);
bB=abcdefzBiom(2);
cB=abcdefzBiom(3);
dB=abcdefzBiom(4);
eB=abcdefzBiom(5);
fB=abcdefzBiom(6);
zB=abcdefzBiom(7);

ySource=-1;
yH2o=-(3*aB+4*eB-dB);
yHco3=aB;
yNh4=cB;
yHpo4=eB;
yHs=fB;
yH=5*aB+bB-4*cB-2*dB+7*eB-fB;
yE=-zB+4*aB+bB-3*cB-2*dB+5*eB-2*fB;
stoichAnStarB=[ySource,yH2o,yHco3,yNh4,yHpo4,yHs,yH,yE];
stoichAnStarB([9,10])=[0,0]; % add additional components: e-acceptor and biomass
stoichAnStarB=-stoichAnStarB;
stoichAnStarB(10)=stoichAnStarB(1);
stoichAnStarB(1)=0;

% Step 2b) "overall" anabolic reaction
stoichAnStar=stoichAnStarB+(1/a)*stoichD;
yEana=stoichAnStar(8);
if yEana>0
    stoichAn=stoichAnStar-yEana/yEa*stoichA;
elseif yEana<0
    stoichAn=stoichAnStar-yEana/yEd*stoichD;
else
    stoichAn=stoichAnStar; 
end

% Step 3: get lambda

% - estimate delGd0 using LaRowe and Van Cappellen (2011)
ne=-z+4*a+b-3*c-2*d+5*e-2*f; % number of electrons transferred in D 
nosc=-ne/a+4; % nominal oxidataion state of carbon 
delGcox0=60.3-28.5*nosc; % kJ/C-mol
delGd0=delGcox0*a*abs(stoichD(1)); % kJ/rxn

% - estimate delGf0 for electron donor
delGf0_D_zero=0;
delGf0_zero=[delGf0_D_zero -237.2 -586.9 -79.5 -1089.1 12.0 0 0 16.5 -67];
delGcox0_zero=delGf0_zero*stoichD';
delGf0_D_est=(delGd0-delGcox0_zero)/stoichD(1);
% - finally, delGf0
delGf0=delGf0_zero;
delGf0(1)=delGf0_D_est;

% - standard delG at pH=0
delGcat0=delGf0*stoichCat';
delGan0=delGf0*stoichAn';

% - stadard delG at pH=7
R=0.008314; % kJ/(K.mol)
T=298.15; % K
iProton=7; % [eD,h2o,hco3-,nh4+,hpo4^2-,hs-,h+,e-,eA,biom]

if pH==0
    delGd=delGd0;
    delGcox=delGcox0;
    delGcat=delGcat0;
    delGan=delGan0;
elseif pH>0
    actHplus=10^(-pH);
    delGd=delGd0+R*T*stoichD(iProton)*log(actHplus);
    delGcox=delGd/a;
    delGcat=delGcat0+R*T*stoichCat(iProton)*log(actHplus);
    delGan=delGan0+R*T*stoichAn(iProton)*log(actHplus);
else
    error('pH cannot be negative!')
end

% The Thermodynamic Electron Equivalents Model (TEEM)
% --------
eta=0.43; 
delGsyn=200; % kJ/(mol.X) (BOX 9 in Kleerebezem and Van Loosdrecht, 2010)

if delGan<0
    m=1;
else
    m=-1;
end
% In calculating lambda across pH values, we assume delGsyn0=delGsyn 
% because chemical composition of biomass builidng block(X)
% is assumeed to be the same as biomass (i.e., CH1.8O0.5N0.2 for both)
% in this case stoich coeff. of h+ is zero, i.e., no pH effect)
% (see pp. 30-31 in Kleerebezem and Van Loosdrecht (2010)

lambda=(delGan*eta^m+delGsyn)/(-delGcat*eta);

if lambda>0
    stoichMet=lambda*stoichCat+stoichAn;
else
    stoichMet=stoichAn;
end
%**************************************************************
% delGdis0=delGf0*stoichMet'; % <---- it should be stoichMet0'
delGdis = lambda*(-delGcat) - delGan; % from the dissipation method 
%**************************************************************
% % delGdis=delGf0*stoichMet'+R*T*stoichMet(iProton)*log(actHplus);
% delGdis=delGdis0+R*T*stoichMet(iProton)*log(actHplus);
% delGdis0=-delGdis0;
% delGdis=-delGdis;

% delGdis=200+18*(6-a)^1.8 + exp(((-0.2+nosc)^2)^0.16*(3.6+0.4*a));
% delGsyn to be equal to the minimum value for delGdis suggested by 
% the dissipation method (200 kJ mol X−1). 
% (p. 35 in Kleerebezem and Van Loosdrecht, 2010)

% updates in v0p3 *******************************************************
%-- add CUE0 and CUE (Carbon Use Efficiency):
% (the order of chemicals in stoich: [e-donor,h2o,hco3,nh4,hpo4,hs,h,e,e-acceptor,biom]
CUE = stoichMet(10)*1/(abs(stoichMet(1))*a);
NUE = stoichMet(10)*0.2/(abs(stoichMet(1))*c+abs(stoichMet(4))*1);
TER = (NUE/CUE)*(1/0.2);
%==========================================================================
