#include "src/common.h"
#include "src/invariants.h"

int main()
{

//-----------------------------------------------------------------SM1Loop----------------------------------------------------------------------------
/*
Param g1("g1", 1);
Param g2("g2", 1);
Param g3("g3", 1);
Param Yu("Yu", 3); Yu.botRight();
Param Yud = Dagger(Yu);
Param Yd("Yd", 3); Yd.botRight();
Param Ydd = Dagger(Yd);
Param Ye("Ye", 3); Ye.botRight();
Param Yed = Dagger(Ye);
Param m2("m2", 1);
Param lab("lab",1);

Function Y2(1);
Y2 = Tr(3*Yu*Yud + 3*Yd*Ydd + Ye*Yed);
Function H(1);
H = Tr(3*Yu*Yud*Yu*Yud + 3*Yd*Ydd*Yd*Ydd + Ye*Yed*Ye*Yed);

vector<BetaFunc> SM1;
vector<BetaFunc> SM2;

BetaFunc bg1(g1);
bg1 = ir(41,10)*g1*g1*g1;
SM1.push_back(bg1);

BetaFunc bg2(g2);
bg2 = ir(-19,6)*g2*g2*g2;
SM1.push_back(bg2);

BetaFunc bg3(g3);
bg3 = -7*g3*g3*g3;
SM1.push_back(bg3);

BetaFunc bYu(Yu);
bYu = ir(3,2)*Yu*Yud*Yu - ir(3,2)*Yu*Ydd*Yd + Y2*Yu- ir(17,20)*Yu*g1*g1 - ir(9,4)*Yu*g2*g2 - 8*Yu*g3*g3;
SM1.push_back(bYu);

BetaFunc bYuc = Conjugate(bYu);
SM1.push_back(bYuc);

BetaFunc bYd(Yd);
bYd = ir(3,2)*Yd*Ydd*Yd - ir(3,2)*Yd*Yud*Yu + Y2*Yd - ir(1,4)*Yd*g1*g1 - ir(9,4)*Yd*g2*g2 - 8*Yd*g3*g3;
SM1.push_back(bYd);

BetaFunc bYdc = Conjugate(bYd);
SM1.push_back(bYdc);

BetaFunc bYe(Ye);
bYe = ir(3,2)*Ye*Yed*Ye + Y2*Ye - ir(9,4)*Ye*g1*g1 - ir(9,4)*Ye*g2*g2;
SM1.push_back(bYe);

BetaFunc bYec = Conjugate(bYe);
SM1.push_back(bYec);

BetaFunc bm2(m2);
bm2 = 6*m2*lab + 2*Y2*m2 - ir(9,10)*m2*g1*g1 - ir(9,2)*m2*g2*g2;
SM1.push_back(bm2);

BetaFunc blab(lab);
blab = 12*lab*lab - ir(9,5)*lab*g1*g1 - 9*lab*g2*g2 + ir(27,100)*g1*g1*g1*g1 + ir(9,10)*g1*g1*g2*g2 + ir(9,4)*g2*g2 + 4*Y2*lab - 4*H;
SM1.push_back(blab);

findInvariants(SM1, SM2);

*/

//-------------------------------------------------------------------MSSMParameters------------------------------------------------------------------------

Param g1("g1", 1);
Param g2("g2", 1);
Param g3("g3", 1);
Param M1("M1", 1);
Param M2("M2", 1);
Param M3("M3", 1);
Param mu("mu", 1);

Param Yt("Yu", 3); Yt.botRight();
Param Ytd = Dagger(Yt);
Param Yb("Yd", 3); Yb.botRight();
Param Ybd = Dagger(Yb);
Param Ye("Ye", 3); Ye.botRight();
Param Yed = Dagger(Ye);

Param ht("hu", 3); ht.botRight();
Param htd = Dagger(ht);
Param hb("hd", 3); hb.botRight();
Param hbd = Dagger(hb);
Param he("he", 3); he.botRight();
Param hed = Dagger(he);

Param B("B", 1);
Param Bd = Dagger(B);
Param mHu("mHu", 1);
Param mHd("mHd", 1);
Param mQ("mQ", 3); mQ.botRightDiag();
Param mL("mL", 3); mL.botRightDiag();
Param mt("mt", 3); mt.botRightDiag();
Param mb("mb", 3); mb.botRightDiag();
Param me("me", 3); me.botRightDiag();

Function S = Tr(mQ) - Tr(mL) - 2*Tr(mt) + Tr(mb) + Tr(me) + Tr(mHu) - Tr(mHd);

//-------------------------------------------------------------------MSSM Beta-Funcs------------------------------------------------------------------------
vector<BetaFunc> MSSM1;
vector<BetaFunc> MSSM2;

BetaFunc bg1_1(g1);
BetaFunc bg1_2(g1);
bg1_1 = g1*g1*g1*ir(33,5);
bg1_2 = g1*g1*g1*g1*g1*ir(199,25) + g1*g1*g1*g2*g2*ir(27,5) + g1*g1*g1*g3*g3*ir(88,5) + g1*g1*g1*(ir(-26,5)*Tr(Ytd*Yt)-ir(14,5)*Tr(Ybd*Yb)-ir(18,5)*Tr(Yed*Ye));
MSSM1.push_back(bg1_1);
MSSM2.push_back(bg1_2);

BetaFunc bg2_1(g2);
BetaFunc bg2_2(g2);
bg2_1 = g2*g2*g2*(1);
bg2_2 = g2*g1*g1*g2*g2*ir(9,5) + g2*g2*g2*g2*g2*(25) + g2*g2*g2*g3*g3*(24) + g2*g2*g2*(-6*Tr(Ytd*Yt)-6*Tr(Ybd*Yb)-2*Tr(Yed*Ye));
MSSM1.push_back(bg2_1);
MSSM2.push_back(bg2_2);

BetaFunc bg3_1(g3);
BetaFunc bg3_2(g3);
bg3_1 = g3*g3*g3*(-3);
bg3_2 = g3*g1*g1*g3*g3*ir(11,5) + g3*g2*g2*g3*g3*(9) + g3*g3*g3*g3*g3*(14) + g3*g3*g3*(-4*Tr(Ytd*Yt)-4*Tr(Ybd*Yb));
MSSM1.push_back(bg3_1);
MSSM2.push_back(bg3_2);

BetaFunc bM1_1(M1);
BetaFunc bM1_2(M1);
bM1_1 = g1*g1*M1*ir(66,5);
bM1_2 = g1*g1*g1*g1*M1*ir(796,25) + g1*g1*g2*g2*M1*ir(54,5) + g1*g1*g2*g2*M2*ir(54,5) + g1*g1*g3*g3*M1*ir(176,5) + g1*g1*g3*g3*M3*ir(176,5) +
	g1*g1*M1*(ir(-52,5)*Tr(Ytd*Yt)-ir(28,5)*Tr(Ybd*Yb)-ir(36,5)*Tr(Yed*Ye)) + g1*g1*(ir(52,5)*Tr(ht*Ytd) + ir(28,5)*Tr(hb*Ybd) + ir(36,5)*Tr(he*Yed));
MSSM1.push_back(bM1_1);
MSSM2.push_back(bM1_2);

BetaFunc bM2_1(M2);
BetaFunc bM2_2(M2);
bM2_1 = g2*g2*M2*(2);
bM2_2 = g1*g1*g2*g2*M1*ir(18,5) + g1*g1*g2*g2*M2*ir(18,5) + g2*g2*g2*g2*M2*(100) + g2*g2*g3*g3*M2*(48) + g2*g2*g3*g3*M3*(48) + g2*g2*M2*(-12*Tr(Ytd*Yt)-12*Tr(Ybd*Yb)-4*Tr(Yed*Ye)) +
	g2*g2*(12*Tr(ht*Ytd) + 12*Tr(hb*Ybd) + 4*Tr(he*Yed));
MSSM1.push_back(bM2_1);
MSSM2.push_back(bM2_2);

BetaFunc bM3_1(M3);
BetaFunc bM3_2(M3);
bM3_1 = g3*g3*M3*(-6);
bM3_2 = g1*g1*g3*g3*M1*ir(22,5) + g1*g1*g3*g3*M3*ir(22,5) + g2*g2*g3*g3*M2*(18) + g2*g2*g3*g3*M3*(18) + g3*g3*g3*g3*M3*(56) + g3*g3*M3*(-8*Tr(Ytd*Yt)-8*Tr(Ybd*Yb)) + g3*g3*(8*Tr(ht*Ytd) +
	8*Tr(hb*Ybd));
MSSM1.push_back(bM3_1);
MSSM2.push_back(bM3_2);

BetaFunc bYt_1(Yt);
BetaFunc bYt_2(Yt);
bYt_1 = Yt*g1*g1*ir(-13,15) + Yt*g2*g2*(-3) + Yt*g3*g3*ir(-16,3) + Yt*Ytd*Yt*(3) + Yt*Ybd*Yb*(1) + Yt*(3*Tr(Ytd*Yt));
bYt_2 = Yt*g1*g1*g1*g1*ir(2743,450) + Yt*g1*g1*g2*g2*(1) + Yt*g1*g1*g3*g3*ir(136,45) + Yt*g1*g1*(ir(4,5)*Tr(Ytd*Yt)) + Yt*g2*g2*g2*g2*ir(15,2) + Yt*g2*g2*g3*g3*(8) + Yt*g3*g3*g3*g3*ir(-16,9) +
	Yt*g3*g3*(16*Tr(Ytd*Yt)) + Yt*Ytd*Yt*g1*g1*ir(2,5) + Yt*Ytd*Yt*g2*g2*(6) + Yt*Ytd*Yt*Ytd*Yt*(-4) + Yt*Ytd*Yt*(-9*Tr(Ytd*Yt)) + Yt*Ybd*Yb*g1*g1*ir(2,5) +
	Yt*Ybd*Yb*Ytd*Yt*(-2) + Yt*Ybd*Yb*Ybd*Yb*(-2) + Yt*Ybd*Yb*(-3*Tr(Ybd*Yb)-Tr(Yed*Ye)) + Yt*(-9*Tr(Ytd*Yt*Ytd*Yt)-3*Tr(Ytd*Yt*Ybd*Yb));
MSSM1.push_back(bYt_1);
MSSM2.push_back(bYt_2);

BetaFunc bYtc_1 = Conjugate(bYt_1);
BetaFunc bYtc_2 = Conjugate(bYt_2);
MSSM1.push_back(bYtc_1);
MSSM2.push_back(bYtc_2);

BetaFunc bYb_1(Yb);
BetaFunc bYb_2(Yb);
bYb_1 = Yb*g1*g1*ir(-7,15) + Yb*g2*g2*(-3) + Yb*g3*g3*ir(-16,3) + Yb*Ytd*Yt*(1) + Yb*Ybd*Yb*(3) + Yb*(3*Tr(Ybd*Yb) + Tr(Yed*Ye));
bYb_2 = Yb*g1*g1*g1*g1*ir(287,90) + Yb*g1*g1*g2*g2*(1) + Yb*g1*g1*g3*g3*ir(8,9) + Yb*g1*g1*(ir(-2,5)*Tr(Ybd*Yb) + ir(6,5)*Tr(Yed*Ye)) + Yb*g2*g2*g2*g2*ir(15,2) + Yb*g2*g2*g3*g3*(8) +
	Yb*g3*g3*g3*g3*ir(-16,9) + Yb*g3*g3*(16*Tr(Ybd*Yb)) + Yb*Ytd*Yt*g1*g1*ir(4,5) + Yb*Ytd*Yt*Ytd*Yt*(-2) + Yb*Ytd*Yt*Ybd*Yb*(-2) + Yb*Ytd*Yt*(-3*Tr(Ytd*Yt)) + Yb*Ybd*Yb*g1*g1*ir(4,5) +
	Yb*Ybd*Yb*g2*g2*(6) + Yb*Ybd*Yb*Ybd*Yb*(-4) + Yb*Ybd*Yb*(-9*Tr(Ybd*Yb)-3*Tr(Yed*Ye)) + Yb*(-3*Tr(Ytd*Yt*Ybd*Yb)-9*Tr(Ybd*Yb*Ybd*Yb)-3*Tr(Yed*Ye*Yed*Ye));
MSSM1.push_back(bYb_1);
MSSM2.push_back(bYb_2);

BetaFunc bYbc_1 = Conjugate(bYb_1);
BetaFunc bYbc_2 = Conjugate(bYb_2);
MSSM1.push_back(bYbc_1);
MSSM2.push_back(bYbc_2);

BetaFunc bYe_1(Ye);
BetaFunc bYe_2(Ye);
bYe_1 = Ye*g1*g1*ir(-9,5) + Ye*g2*g2*(-3) + Ye*Yed*Ye*(3) + Ye*(3*Tr(Ybd*Yb) + Tr(Yed*Ye));
bYe_2 = Ye*g1*g1*g1*g1*ir(27,2) + Ye*g1*g1*g2*g2*ir(9,5) + Ye*g1*g1*(ir(-2,5)*Tr(Ybd*Yb) + ir(6,5)*Tr(Yed*Ye)) + Ye*g2*g2*g2*g2*ir(15,2) + Ye*g3*g3*(16*Tr(Ybd*Yb)) + Ye*Yed*Ye*g2*g2*(6) +
	Ye*Yed*Ye*Yed*Ye*(-4) + Ye*Yed*Ye*(-9*Tr(Ybd*Yb)-3*Tr(Yed*Ye)) + Ye*(-3*Tr(Ytd*Yt*Ybd*Yb)-9*Tr(Ybd*Yb*Ybd*Yb)-3*Tr(Yed*Ye*Yed*Ye));
MSSM1.push_back(bYe_1);
MSSM2.push_back(bYe_2);

BetaFunc bYec_1 = Conjugate(bYe_1);
BetaFunc bYec_2 = Conjugate(bYe_2);
MSSM1.push_back(bYec_1);
MSSM2.push_back(bYec_2);

BetaFunc bht_1(ht);
BetaFunc bht_2(ht);
bht_1 = ht*g1*g1*ir(-13,15) + ht*g2*g2*(-3) + ht*g3*g3*ir(-16,3) + ht*Ytd*Yt*(5) + ht*Ybd*Yb*(1) + ht*(3*Tr(Ytd*Yt)) +
	Yt*g1*g1*M1*ir(26,15) + Yt*g2*g2*M2*(6) + Yt*g3*g3*M3*ir(32,3) + Yt*Ytd*ht*(4) + Yt*Ybd*hb*(2) + Yt*(6*Tr(ht*Ytd));
bht_2 = ht*g1*g1*g1*g1*ir(2743,450) + ht*g1*g1*g2*g2*(1) + ht*g1*g1*g3*g3*ir(136,45) + ht*g1*g1*(ir(4,5)*Tr(Ytd*Yt)) + ht*g2*g2*g2*g2*ir(15,2) + ht*g2*g2*g3*g3*(8) +
	ht*g3*g3*g3*g3*ir(-16,9) + ht*g3*g3*(16*Tr(Ytd*Yt)) + ht*Ytd*Yt*g2*g2*(12) + ht*Ytd*Yt*Ytd*Yt*(-6) + ht*Ytd*Yt*(-15*Tr(Ytd*Yt)) + ht*Ybd*Yb*g1*g1*ir(2,5) +
	ht*Ybd*Yb*Ytd*Yt*(-4) + ht*Ybd*Yb*Ybd*Yb*(-2) + ht*Ybd*Yb*(-3*Tr(Ybd*Yb)-Tr(Yed*Ye)) + ht*(-9*Tr(Ytd*Yt*Ytd*Yt)-3*Tr(Ytd*Yt*Ybd*Yb)) + Yt*g1*g1*g1*g1*M1*ir(-5486,225) +
	Yt*g1*g1*g2*g2*M1*(-2) + Yt*g1*g1*g2*g2*M2*(-2) + Yt*g1*g1*g3*g3*M1*ir(-272,45) + Yt*g1*g1*g3*g3*M3*ir(-272,45) + Yt*g1*g1*M1*(ir(-8,5)*Tr(Ytd*Yt)) + Yt*g1*g1*(ir(8,5)*Tr(ht*Ytd)) +
	Yt*g2*g2*g2*g2*M2*(-30) + Yt*g2*g2*g3*g3*M2*(-16) + Yt*g2*g2*g3*g3*M3*(-16) + Yt*g3*g3*g3*g3*M3*ir(64,9) + Yt*g3*g3*M3*(-32*Tr(Ytd*Yt)) + Yt*g3*g3*(32*Tr(ht*Ytd)) +
	Yt*Ytd*ht*g1*g1*ir(6,5) + Yt*Ytd*ht*g2*g2*(6) + Yt*Ytd*ht*Ytd*Yt*(-8) + Yt*Ytd*ht*(-12*Tr(Ytd*Yt)) + Yt*Ytd*Yt*g1*g1*M1*ir(-4,5) + Yt*Ytd*Yt*g2*g2*M2*(-12) + Yt*Ytd*Yt*Ytd*ht*(-6) +
	Yt*Ytd*Yt*(-18*Tr(ht*Ytd)) + Yt*Ybd*hb*g1*g1*ir(4,5) + Yt*Ybd*hb*Ytd*Yt*(-4) + Yt*Ybd*hb*Ybd*Yb*(-4) + Yt*Ybd*hb*(-6*Tr(Ybd*Yb)-2*Tr(Yed*Ye)) + Yt*Ybd*Yb*g1*g1*M1*ir(-4,5) +
	Yt*Ybd*Yb*Ytd*ht*(-2) + Yt*Ybd*Yb*Ybd*hb*(-4) + Yt*Ybd*Yb*(-6*Tr(hb*Ybd)-2*Tr(he*Yed)) + Yt*(-6*Tr(hb*Ytd*Yt*Ybd)-36*Tr(Ytd*ht*Ytd*Yt)-6*Tr(Ytd*ht*Ybd*Yb));
MSSM1.push_back(bht_1);
MSSM2.push_back(bht_2);

BetaFunc bhtd_1 = Conjugate(bht_1);
BetaFunc bhtd_2 = Conjugate(bht_2);
MSSM1.push_back(bhtd_1);
MSSM2.push_back(bhtd_2);

BetaFunc bhb_1(hb);
BetaFunc bhb_2(hb);
bhb_1 =	hb*g1*g1*ir(-7,15) + hb*g2*g2*(-3) + hb*g3*g3*ir(-16,3) + hb*Ytd*Yt*(1) + hb*Ybd*Yb*(5) + hb*(3*Tr(Ybd*Yb) + Tr(Yed*Ye)) + Yb*g1*g1*M1*ir(14,15) +
	Yb*g2*g2*M2*(6) + Yb*g3*g3*M3*ir(32,3) + Yb*Ytd*ht*(2) + Yb*Ybd*hb*(4) + Yb*(6*Tr(hb*Ybd) + 2*Tr(he*Yed));
bhb_2 =	hb*g1*g1*g1*g1*ir(287,90) + hb*g1*g1*g2*g2*(1) + hb*g1*g1*g3*g3*ir(8,9) + hb*g1*g1*(ir(-2,5)*Tr(Ybd*Yb) + ir(6,5)*Tr(Yed*Ye)) + hb*g2*g2*g2*g2*ir(15,2) + hb*g2*g2*g3*g3*(8) +
	hb*g3*g3*g3*g3*ir(-16,9) + hb*g3*g3*(16*Tr(Ybd*Yb)) + hb*Ytd*Yt*g1*g1*ir(4,5) + hb*Ytd*Yt*Ytd*Yt*(-2) + hb*Ytd*Yt*Ybd*Yb*(-4) + hb*Ytd*Yt*(-3*Tr(Ytd*Yt)) + hb*Ybd*Yb*g1*g1*ir(6,5) +
 	hb*Ybd*Yb*g2*g2*(12) + hb*Ybd*Yb*Ybd*Yb*(-6) + hb*Ybd*Yb*(-15*Tr(Ybd*Yb)-5*Tr(Yed*Ye)) + hb*(-3*Tr(Ytd*Yt*Ybd*Yb)-9*Tr(Ybd*Yb*Ybd*Yb)-3*Tr(Yed*Ye*Yed*Ye)) + Yb*g1*g1*g1*g1*M1*ir(-574,45) +
	Yb*g1*g1*g2*g2*M1*(-2) + Yb*g1*g1*g2*g2*M2*(-2) + Yb*g1*g1*g3*g3*M1*ir(-16,9) + Yb*g1*g1*g3*g3*M3*ir(-16,9) + Yb*g1*g1*M1*(ir(4,5)*Tr(Ybd*Yb)-ir(12,5)*Tr(Yed*Ye)) +
	Yb*g1*g1*(ir(-4,5)*Tr(hb*Ybd) + ir(12,5)*Tr(he*Yed)) + Yb*g2*g2*g2*g2*M2*(-30) + Yb*g2*g2*g3*g3*M2*(-16) + Yb*g2*g2*g3*g3*M3*(-16) + Yb*g3*g3*g3*g3*M3*ir(64,9) + Yb*g3*g3*M3*(-32*Tr(Ybd*Yb)) +
	Yb*g3*g3*(32*Tr(hb*Ybd)) + Yb*Ytd*ht*g1*g1*ir(8,5) + Yb*Ytd*ht*Ytd*Yt*(-4) + Yb*Ytd*ht*Ybd*Yb*(-4) + Yb*Ytd*ht*(-6*Tr(Ytd*Yt)) + Yb*Ytd*Yt*g1*g1*M1*ir(-8,5) + Yb*Ytd*Yt*Ytd*ht*(-4) +
	Yb*Ytd*Yt*Ybd*hb*(-2) + Yb*Ytd*Yt*(-6*Tr(ht*Ytd)) + Yb*Ybd*hb*g1*g1*ir(6,5) + Yb*Ybd*hb*g2*g2*(6) + Yb*Ybd*hb*Ybd*Yb*(-8) + Yb*Ybd*hb*(-12*Tr(Ybd*Yb)-4*Tr(Yed*Ye)) +
	Yb*Ybd*Yb*g1*g1*M1*ir(-8,5) + Yb*Ybd*Yb*g2*g2*M2*(-12) + Yb*Ybd*Yb*Ybd*hb*(-6) + Yb*Ybd*Yb*(-18*Tr(hb*Ybd)-6*Tr(he*Yed)) +
	Yb*(-6*Tr(hb*Ytd*Yt*Ybd)-6*Tr(Ytd*ht*Ybd*Yb)-36*Tr(Ybd*hb*Ybd*Yb)-12*Tr(Yed*he*Yed*Ye));
MSSM1.push_back(bhb_1);
MSSM2.push_back(bhb_2);

BetaFunc bhbd_1 = Conjugate(bhb_1);
BetaFunc bhbd_2 = Conjugate(bhb_2);
MSSM1.push_back(bhbd_1);
MSSM2.push_back(bhbd_2);

BetaFunc bhe_1(he);
BetaFunc bhe_2(he);
bhe_1 = he*g1*g1*ir(-9,5) + he*g2*g2*(-3) + he*Yed*Ye*(5) + he*(3*Tr(Ybd*Yb) + Tr(Yed*Ye)) + Ye*g1*g1*M1*ir(18,5) +
	Ye*g2*g2*M2*(6) + Ye*Yed*he*(4) + Ye*(6*Tr(hb*Ybd) + 2*Tr(he*Yed));
bhe_2 = he*g1*g1*g1*g1*ir(27,2) + he*g1*g1*g2*g2*ir(9,5) + he*g1*g1*(ir(-2,5)*Tr(Ybd*Yb) + ir(6,5)*Tr(Yed*Ye)) + he*g2*g2*g2*g2*ir(15,2) + he*g3*g3*(16*Tr(Ybd*Yb)) + he*Yed*Ye*g1*g1*ir(-6,5) +
	he*Yed*Ye*g2*g2*(12) + he*Yed*Ye*Yed*Ye*(-6) + he*Yed*Ye*(-15*Tr(Ybd*Yb)-5*Tr(Yed*Ye)) + he*(-3*Tr(Ytd*Yt*Ybd*Yb)-9*Tr(Ybd*Yb*Ybd*Yb)-3*Tr(Yed*Ye*Yed*Ye)) + Ye*g1*g1*g1*g1*M1*(-54) +
	Ye*g1*g1*g2*g2*M1*ir(-18,5) + Ye*g1*g1*g2*g2*M2*ir(-18,5) + Ye*g1*g1*M1*(ir(4,5)*Tr(Ybd*Yb)-ir(12,5)*Tr(Yed*Ye)) + Ye*g1*g1*(ir(-4,5)*Tr(hb*Ybd) + ir(12,5)*Tr(he*Yed)) +
	Ye*g2*g2*g2*g2*M2*(-30) + Ye*g3*g3*M3*(-32*Tr(Ybd*Yb)) + Ye*g3*g3*(32*Tr(hb*Ybd)) + Ye*Yed*he*g1*g1*ir(6,5) + Ye*Yed*he*g2*g2*(6) + Ye*Yed*he*Yed*Ye*(-8) +
	Ye*Yed*he*(-12*Tr(Ybd*Yb)-4*Tr(Yed*Ye)) + Ye*Yed*Ye*g2*g2*M2*(-12) + Ye*Yed*Ye*Yed*he*(-6) + Ye*Yed*Ye*(-18*Tr(hb*Ybd)-6*Tr(he*Yed)) +
	Ye*(-6*Tr(hb*Ytd*Yt*Ybd)-6*Tr(Ytd*ht*Ybd*Yb)-36*Tr(Ybd*hb*Ybd*Yb)-12*Tr(Yed*he*Yed*Ye));
MSSM1.push_back(bhe_1);
MSSM2.push_back(bhe_2);

BetaFunc bhed_1 = Conjugate(bhe_1);
BetaFunc bhed_2 = Conjugate(bhe_2);
MSSM1.push_back(bhed_1);
MSSM2.push_back(bhed_2);

BetaFunc bmHu_1(mHu);
BetaFunc bmHu_2(mHu);
bmHu_1 = g1*g1*M1*M1*ir(-6,5) + g2*g2*M2*M2*(-6) + 6*Tr(ht*htd) + 6*Tr(Ytd*Yt)*mHu + 6*Tr(Ytd*Yt*mQ) + 6*Tr(Yt*Ytd*mt) + ir(3,5)*g1*g1*S;
bmHu_2 =g1*g1*g1*g1*M1*M1*ir(621,25) + g1*g1*g1*g1*(ir(24,25)*Tr(mt) + ir(6,25)*Tr(mb) + ir(3,25)*Tr(mQ) + ir(18,25)*Tr(me) + ir(9,25)*Tr(mL) + ir(9,25)*mHd + ir(9,25)*mHu) +
	g1*g1*g2*g2*M1*M1*ir(18,5) + g1*g1*g2*g2*M1*M2*ir(18,5) + g1*g1*g2*g2*M2*M2*ir(18,5) + g1*g1*M1*M1*(ir(16,5)*Tr(Ytd*Yt)) + g1*g1*M1*(ir(-8,5)*Tr(ht*Ytd)-ir(8,5)*Tr(htd*Yt)) +
	g1*g1*(ir(8,5)*Tr(ht*htd) + ir(8,5)*Tr(Ytd*Yt)*mHu + ir(8,5)*Tr(Ytd*Yt*mQ) + ir(8,5)*Tr(Yt*Ytd*mt)) + g2*g2*g2*g2*M2*M2*(33) + g2*g2*g2*g2*(9*Tr(mQ) + 3*Tr(mL) + 3*mHd + 3*mHu) +
	g3*g3*M3*M3*(64*Tr(Ytd*Yt)) + g3*g3*M3*(-32*Tr(ht*Ytd)-32*Tr(htd*Yt)) + g3*g3*(32*Tr(ht*htd) + 32*Tr(Ytd*Yt)*mHu + 32*Tr(Ytd*Yt*mQ) + 32*Tr(Yt*Ytd*mt)) -
	6*Tr(hb*htd*Yt*Ybd)-36*Tr(htd*ht*Ytd*Yt)-6*Tr(htd*ht*Ybd*Yb)-6*Tr(hbd*hb*Ytd*Yt)-36*Tr(Ytd*ht*htd*Yt)-6*Tr(Ytd*ht*hbd*Yb)-36*Tr(Ytd*Yt*Ytd*Yt)*mHu-36*Tr(Ytd*Yt*Ytd*Yt*mQ)-
	6*Tr(Ytd*Yt*Ybd*Yb)*mHd-6*Tr(Ytd*Yt*Ybd*Yb)*mHu-6*Tr(Ytd*Yt*Ybd*Yb*mQ)-36*Tr(Yt*Ytd*Yt*Ytd*mt)-6*Tr(Yt*Ybd*Yb*Ytd*mt)-6*Tr(Ybd*Yb*Ytd*Yt*mQ)-6*Tr(Yb*Ytd*Yt*Ybd*mb);
MSSM1.push_back(bmHu_1);
MSSM2.push_back(bmHu_2);

BetaFunc bmHd_1(mHd);
BetaFunc bmHd_2(mHd);
bmHd_1 = g1*g1*M1*M1*ir(-6,5) + g2*g2*M2*M2*(-6) + 6*Tr(hb*hbd) + 2*Tr(he*hed) + 6*Tr(Ybd*Yb)*mHd + 6*Tr(Ybd*Yb*mQ) +
	6*Tr(Yb*Ybd*mb) + 2*Tr(Yed*Ye)*mHd + 2*Tr(Yed*Ye*mL) + 2*Tr(Ye*Yed*me) + ir(-3,5)*g1*g1*S;
bmHd_2 =g1*g1*g1*g1*M1*M1*ir(621,25) + g1*g1*g1*g1*(ir(24,25)*Tr(mt) + ir(6,25)*Tr(mb) + ir(3,25)*Tr(mQ) + ir(18,25)*Tr(me) + ir(9,25)*Tr(mL) + ir(9,25)*mHd + ir(9,25)*mHu) +
	g1*g1*g2*g2*M1*M1*ir(18,5) + g1*g1*g2*g2*M1*M2*ir(18,5) + g1*g1*g2*g2*M2*M2*ir(18,5) + g1*g1*M1*M1*(ir(-8,5)*Tr(Ybd*Yb) + ir(24,5)*Tr(Yed*Ye)) +
	g1*g1*M1*(ir(4,5)*Tr(hb*Ybd)-ir(12,5)*Tr(he*Yed) + ir(4,5)*Tr(hbd*Yb)-ir(12,5)*Tr(hed*Ye)) + g1*g1*(ir(-4,5)*Tr(hb*hbd) + ir(12,5)*Tr(he*hed)-ir(4,5)*Tr(Ybd*Yb)*mHd-
	ir(4,5)*Tr(Ybd*Yb*mQ)-ir(4,5)*Tr(Yb*Ybd*mb) + ir(12,5)*Tr(Yed*Ye)*mHd + ir(12,5)*Tr(Yed*Ye*mL) + ir(12,5)*Tr(Ye*Yed*me)) + g2*g2*g2*g2*M2*M2*(33) + g2*g2*g2*g2*(9*Tr(mQ) +
	3*Tr(mL) + 3*mHd + 3*mHu) + g3*g3*M3*M3*(64*Tr(Ybd*Yb)) + g3*g3*M3*(-32*Tr(hb*Ybd)-32*Tr(hbd*Yb)) + g3*g3*(32*Tr(hb*hbd) + 32*Tr(Ybd*Yb)*mHd + 32*Tr(Ybd*Yb*mQ) + 32*Tr(Yb*Ybd*mb)) -
	6*Tr(hb*htd*Yt*Ybd)-6*Tr(htd*ht*Ybd*Yb)-6*Tr(hbd*hb*Ytd*Yt)-36*Tr(hbd*hb*Ybd*Yb)-12*Tr(hed*he*Yed*Ye)-6*Tr(Ytd*ht*hbd*Yb)-6*Tr(Ytd*Yt*Ybd*Yb)*mHd-6*Tr(Ytd*Yt*Ybd*Yb)*mHu-
	6*Tr(Ytd*Yt*Ybd*Yb*mQ)-6*Tr(Yt*Ybd*Yb*Ytd*mt)-36*Tr(Ybd*hb*hbd*Yb)-6*Tr(Ybd*Yb*Ytd*Yt*mQ)-36*Tr(Ybd*Yb*Ybd*Yb)*mHd-36*Tr(Ybd*Yb*Ybd*Yb*mQ)-6*Tr(Yb*Ytd*Yt*Ybd*mb)-
	36*Tr(Yb*Ybd*Yb*Ybd*mb)-12*Tr(Yed*he*hed*Ye)-12*Tr(Yed*Ye*Yed*Ye)*mHd-12*Tr(Yed*Ye*Yed*Ye*mL)-12*Tr(Ye*Yed*Ye*Yed*me);
MSSM1.push_back(bmHd_1);
MSSM2.push_back(bmHd_2);

BetaFunc bmt_1(mt);
BetaFunc bmt_2(mt);
bmt_1 = g1*g1*M1*M1*ir(-32,15) + g3*g3*M3*M3*ir(-32,3) + ht*(4*htd) + Yt*Ytd*mt*(2) + Yt*Ytd*(4*mHu) + Yt*mQ*Ytd*(4) + mt*Yt*Ytd*(2) + ir(-4,5)*g1*g1*S;
bmt_2 = g1*g1*g1*g1*M1*M1*ir(3424,75) + g1*g1*g1*g1*(ir(128,75)*Tr(mt) + ir(32,75)*Tr(mb) + ir(16,75)*Tr(mQ) + ir(32,25)*Tr(me) + ir(16,25)*Tr(mL) + ir(16,25)*mHd + ir(16,25)*mHu) +
	g1*g1*g3*g3*M1*M1*(ir(512,45)) + g1*g1*g3*g3*M1*M3*(ir(512,45)) + g1*g1*g3*g3*M3*M3*ir(512,45) + g3*g3*g3*g3*M3*M3*ir(-128,3) + g3*g3*g3*g3*(ir(16,3)*Tr(mt) + ir(16,3)*Tr(mb) +
	ir(32,3)*Tr(mQ)) + ht*g1*g1*(ir(-4,5)*htd) + ht*g2*g2*(12*htd) + ht*htd*Yt*Ytd*(-4) + ht*hbd*Yb*Ytd*(-4) + ht*Ytd*g1*g1*M1*ir(4,5) + ht*Ytd*g2*g2*M2*(-12) + ht*Ytd*Yt*(-4*htd) +
	ht*Ytd*(-12*Tr(htd*Yt)) + ht*Ybd*Yb*(-4*htd) + ht*(-12*htd*Tr(Ytd*Yt)) + Yt*g1*g1*M1*(ir(4,5)*htd) + Yt*g2*g2*M2*(-12*htd) + Yt*htd*ht*Ytd*(-4) + Yt*hbd*hb*Ytd*(-4) +
	Yt*Ytd*g1*g1*M1*M1*ir(-8,5) + Yt*Ytd*g1*g1*(ir(-4,5)*mHu) + Yt*Ytd*g2*g2*M2*M2*(24) + Yt*Ytd*g2*g2*(12*mHu) + Yt*Ytd*ht*(- 4*htd) + Yt*Ytd*Yt*Ytd*mt*(- 2) +
	Yt*Ytd*Yt*Ytd*(-8*mHu) + Yt*Ytd*Yt*mQ*Ytd*(-4) + Yt*Ytd*mt*g1*g1*(ir(-2,5)) + Yt*Ytd*mt*g2*g2*(6) + Yt*Ytd*mt*Yt*Ytd*(-4) + Yt*Ytd*mt*(-6*Tr(Ytd*Yt)) + Yt*Ytd*(-12*Tr(ht*htd) -
	24*Tr(Ytd*Yt)*mHu - 12*Tr(Ytd*Yt*mQ) - 12*Tr(Yt*Ytd*mt)) + Yt*Ybd*hb*(-4*htd) + Yt*Ybd*Yb*Ytd*mt*(-2) + Yt*Ybd*Yb*Ytd*(-4*mHd - 4*mHu) + Yt*Ybd*Yb*mQ*Ytd*(-4) +
	Yt*Ybd*mb*Yb*Ytd*(-4) + Yt*mQ*Ytd*g1*g1*ir(-4,5) + Yt*mQ*Ytd*g2*g2*(12) + Yt*mQ*Ytd*Yt*Ytd*(-4) + Yt*mQ*Ytd*(-12*Tr(Ytd*Yt)) + Yt*mQ*Ybd*Yb*Ytd*(-4) +
	Yt*(-12*htd*Tr(ht*Ytd)) + mt*Yt*Ytd*g1*g1*ir(-2,5) + mt*Yt*Ytd*g2*g2*(6) + mt*Yt*Ytd*Yt*Ytd*(-2) + mt*Yt*Ytd*(-6*Tr(Ytd*Yt)) + mt*Yt*Ybd*Yb*Ytd*(-2);
MSSM1.push_back(bmt_1);
MSSM2.push_back(bmt_2);

BetaFunc bmb_1(mb);
BetaFunc bmb_2(mb);
bmb_1 = g1*g1*M1*M1*ir(-8,15) + g3*g3*M3*M3*ir(-32,3) + hb*(4*hbd) + Yb*Ybd*mb*(2) + Yb*Ybd*(4*mHd) + Yb*mQ*Ybd*(4) + mb*Yb*Ybd*(2) + ir(2,5)*g1*g1*S;
bmb_2 = g1*g1*g1*g1*M1*M1*ir(808,75) + g1*g1*g1*g1*(ir(32,75)*Tr(mt) + ir(8,75)*Tr(mb) + ir(4,75)*Tr(mQ) + ir(8,25)*Tr(me) + ir(4,25)*Tr(mL) + ir(4,25)*mHd + ir(4,25)*mHu) +
	g1*g1*g3*g3*M1*M1*ir(128,45) + g1*g1*g3*g3*M1*M3*ir(128,45) + g1*g1*g3*g3*M3*M3*ir(128,45) + g3*g3*g3*g3*M3*M3*ir(-128,3) + g3*g3*g3*g3*(ir(16,3)*Tr(mt) + ir(16,3)*Tr(mb) + ir(32,3)*Tr(mQ)) +
 	hb*g1*g1*(ir(4,5)*hbd) + hb*g2*g2*(12*hbd) + hb*htd*Yt*Ybd*(-4) + hb*hbd*Yb*Ybd*(-4) + hb*Ytd*Yt*(-4*hbd) + hb*Ybd*g1*g1*M1*ir(-4,5) + hb*Ybd*g2*g2*M2*(-12) + hb*Ybd*Yb*(-4*hbd) +
	hb*Ybd*(-12*Tr(hbd*Yb) - 4*Tr(hed*Ye)) + hb*(-12*hbd*Tr(Ybd*Yb) - 4*hbd*Tr(Yed*Ye)) + Yb*g1*g1*M1*(ir(-4,5)*hbd) + Yb*g2*g2*M2*(-12*hbd) + Yb*htd*ht*Ybd*(-4) + Yb*hbd*hb*Ybd*(-4) +
	Yb*Ytd*ht*(-4*hbd) + Yb*Ytd*Yt*Ybd*mb*(-2) + Yb*Ytd*Yt*Ybd*(-4*mHd - 4*mHu) + Yb*Ytd*Yt*mQ*Ybd*(-4) + Yb*Ytd*mt*Yt*Ybd*(-4) + Yb*Ybd*g1*g1*M1*M1*ir(8,5) + Yb*Ybd*g1*g1*(ir(4,5)*mHd) +
	Yb*Ybd*g2*g2*M2*M2*(24) + Yb*Ybd*g2*g2*(12*mHd) + Yb*Ybd*hb*(-4*hbd) + Yb*Ybd*Yb*Ybd*mb*(-2) + Yb*Ybd*Yb*Ybd*(-8*mHd) + Yb*Ybd*Yb*mQ*Ybd*(-4) + Yb*Ybd*mb*g1*g1*ir(2,5) + Yb*Ybd*mb*g2*g2*(6) +
	Yb*Ybd*mb*Yb*Ybd*(-4) + Yb*Ybd*mb*(-6*Tr(Ybd*Yb) - 2*Tr(Yed*Ye)) + Yb*Ybd*(-12*Tr(hb*hbd) - 4*Tr(he*hed) - 24*Tr(Ybd*Yb)*mHd - 12*Tr(Ybd*Yb*mQ) - 12*Tr(Yb*Ybd*mb) - 8*Tr(Yed*Ye)*mHd -
	4*Tr(Yed*Ye*mL)- 4*Tr(Ye*Yed*me)) + Yb*mQ*Ytd*Yt*Ybd*(-4) + Yb*mQ*Ybd*g1*g1*ir(4,5) + Yb*mQ*Ybd*g2*g2*(12) + Yb*mQ*Ybd*Yb*Ybd*(-4) + Yb*mQ*Ybd*(-12*Tr(Ybd*Yb) - 4*Tr(Yed*Ye)) +
	Yb*(-12*hbd*Tr(hb*Ybd) - 4*hbd*Tr(he*Yed)) + mb*Yb*Ytd*Yt*Ybd*(-2) + mb*Yb*Ybd*g1*g1*ir(2,5) + mb*Yb*Ybd*g2*g2*(6) + mb*Yb*Ybd*Yb*Ybd*(-2) + mb*Yb*Ybd*(-6*Tr(Ybd*Yb) - 2*Tr(Yed*Ye));
MSSM1.push_back(bmb_1);
MSSM2.push_back(bmb_2);

BetaFunc bmQ_1(mQ);
BetaFunc bmQ_2(mQ);
bmQ_1 =	g1*g1*M1*M1*ir(-2,15) + g2*g2*M2*M2*(-6) + g3*g3*M3*M3*ir(-32,3) + htd*ht*(2) + hbd*hb*(2) + Ytd*Yt*mQ*(1) + Ytd*Yt*(2*mHu) + Ytd*mt*Yt*(2) +
	Ybd*Yb*mQ*(1) + Ybd*Yb*(2*mHd) + Ybd*mb*Yb*(2) + mQ*Ytd*Yt*(1) + mQ*Ybd*Yb*(1) + ir(1,5)*g1*g1*S;
bmQ_2 =	g1*g1*g1*g1*M1*M1*ir(199,75) + g1*g1*g1*g1*(ir(8,75)*Tr(mt) + ir(2,75)*Tr(mb) + ir(1,75)*Tr(mQ) + ir(2,25)*Tr(me) + ir(1,25)*Tr(mL) + ir(1,25)*mHd + ir(1,25)*mHu) + g1*g1*g2*g2*M1*M1*ir(2,5) +
	g1*g1*g2*g2*M1*M2*ir(2,5) + g1*g1*g2*g2*M2*M2*ir(2,5) + g1*g1*g3*g3*M1*M1*ir(32,45) + g1*g1*g3*g3*M1*M3*ir(32,45) + g1*g1*g3*g3*M3*M3*ir(32,45) + g2*g2*g2*g2*M2*M2*(33) +
	g2*g2*g2*g2*(9*Tr(mQ) + 3*Tr(mL) + 3*mHd + 3*mHu) + g2*g2*g3*g3*M2*M2*(32) + g2*g2*g3*g3*M2*M3*(32) + g2*g2*g3*g3*M3*M3*(32) + g3*g3*g3*g3*M3*M3*ir(-128,3) + g3*g3*g3*g3*(ir(16,3)*Tr(mt) +
	ir(16,3)*Tr(mb) + ir(32,3)*Tr(mQ)) + htd*ht*g1*g1*ir(8,5) + htd*ht*Ytd*Yt*(-4) + htd*ht*(-6*Tr(Ytd*Yt)) + htd*Yt*g1*g1*M1*ir(-8,5) + htd*Yt*Ytd*ht*(-4) + htd*Yt*(-6*Tr(ht*Ytd)) +
	hbd*hb*g1*g1*ir(4,5) + hbd*hb*Ybd*Yb*(-4) + hbd*hb*(-6*Tr(Ybd*Yb)-2*Tr(Yed*Ye)) + hbd*Yb*g1*g1*M1*ir(-4,5) + hbd*Yb*Ybd*hb*(-4) + hbd*Yb*(-6*Tr(hb*Ybd)-2*Tr(he*Yed)) +
	Ytd*ht*g1*g1*M1*ir(-8,5) + Ytd*ht*htd*Yt*(-4) + Ytd*ht*(-6*Tr(htd*Yt)) + Ytd*Yt*g1*g1*M1*M1*ir(16,5) + Ytd*Yt*g1*g1*(ir(8,5)*mHu) + Ytd*Yt*htd*ht*(-4) + Ytd*Yt*Ytd*Yt*mQ*(-2) +
	Ytd*Yt*Ytd*Yt*(-8*mHu) + Ytd*Yt*Ytd*mt*Yt*(-4) + Ytd*Yt*mQ*g1*g1*ir(4,5) + Ytd*Yt*mQ*Ytd*Yt*(-4) + Ytd*Yt*mQ*(-3*Tr(Ytd*Yt)) +
	Ytd*Yt*(-6*Tr(ht*htd)-12*Tr(Ytd*Yt)*mHu-6*Tr(Ytd*Yt*mQ)-6*Tr(Yt*Ytd*mt)) + Ytd*mt*Yt*g1*g1*ir(8,5) + Ytd*mt*Yt*Ytd*Yt*(-4) + Ytd*mt*Yt*(-6*Tr(Ytd*Yt)) + Ybd*hb*g1*g1*M1*ir(-4,5) +
	Ybd*hb*hbd*Yb*(-4) + Ybd*hb*(-6*Tr(hbd*Yb)-2*Tr(hed*Ye)) + Ybd*Yb*g1*g1*M1*M1*ir(8,5) + Ybd*Yb*g1*g1*(ir(4,5)*mHd) + Ybd*Yb*hbd*hb*(-4) + Ybd*Yb*Ybd*Yb*mQ*(-2) + Ybd*Yb*Ybd*Yb*(-8*mHd) +
	Ybd*Yb*Ybd*mb*Yb*(-4) + Ybd*Yb*mQ*g1*g1*ir(2,5) + Ybd*Yb*mQ*Ybd*Yb*(-4) + Ybd*Yb*mQ*(-3*Tr(Ybd*Yb)-Tr(Yed*Ye)) +
	Ybd*Yb*(-6*Tr(hb*hbd)-2*Tr(he*hed)-12*Tr(Ybd*Yb)*mHd-6*Tr(Ybd*Yb*mQ)-6*Tr(Yb*Ybd*mb)-4*Tr(Yed*Ye)*mHd-2*Tr(Yed*Ye*mL)-6*Tr(Ye*Yed*me)) +
	Ybd*mb*Yb*g1*g1*ir(4,5) + Ybd*mb*Yb*Ybd*Yb*(-4) + Ybd*mb*Yb*(-6*Tr(Ybd*Yb)-2*Tr(Yed*Ye)) + mQ*Ytd*Yt*g1*g1*ir(4,5) + mQ*Ytd*Yt*Ytd*Yt*(-2) + mQ*Ytd*Yt*(-3*Tr(Ytd*Yt)) +
	mQ*Ybd*Yb*g1*g1*ir(2,5) + mQ*Ybd*Yb*Ybd*Yb*(-2) + mQ*Ybd*Yb*(-3*Tr(Ybd*Yb)-Tr(Yed*Ye));
MSSM1.push_back(bmQ_1);
MSSM2.push_back(bmQ_2);

BetaFunc bme_1(me);
BetaFunc bme_2(me);
bme_1 =	g1*g1*M1*M1*ir(-24,5) + he*(4*hed) + Ye*Yed*me*(2) + Ye*Yed*(4*mHd) + Ye*mL*Yed*(4) + me*Ye*Yed*(2) + ir(6,5)*g1*g1*S;
bme_2 =	g1*g1*g1*g1*M1*M1*ir(2808,25) + g1*g1*g1*g1*(ir(96,25)*Tr(mt) + ir(24,25)*Tr(mb) + ir(12,25)*Tr(mQ) + ir(72,25)*Tr(me) + ir(36,25)*Tr(mL) + ir(36,25)*mHd + ir(36,25)*mHu) +
	he*g1*g1*(ir(-12,5)*hed) + he*g2*g2*(12*hed) + he*hed*Ye*Yed*(-4) + he*Yed*g1*g1*M1*ir(12,5) + he*Yed*g2*g2*M2*(-12) + he*Yed*Ye*(-4*hed) + he*Yed*(-12*Tr(hbd*Yb)-4*Tr(hed*Ye)) +
	he*(-12*hed*Tr(Ybd*Yb)-4*hed*Tr(Yed*Ye)) + Ye*g1*g1*M1*(ir(12,5)*hed) + Ye*g2*g2*M2*(-12*hed) + Ye*hed*he*Yed*(-4) + Ye*Yed*g1*g1*M1*M1*ir(-24,5) + Ye*Yed*g1*g1*(ir(-12,5)*mHd) +
	Ye*Yed*g2*g2*M2*M2*(24) + Ye*Yed*g2*g2*(12*mHd) + Ye*Yed*he*(-4*hed) + Ye*Yed*Ye*Yed*me*(-2) + Ye*Yed*Ye*Yed*(-8*mHd) + Ye*Yed*Ye*mL*Yed*(-4) + Ye*Yed*me*g1*g1*ir(-6,5) +
	Ye*Yed*me*g2*g2*(6) + Ye*Yed*me*Ye*Yed*(-4) + Ye*Yed*me*(-6*Tr(Ybd*Yb)-2*Tr(Yed*Ye)) +
	Ye*Yed*(-12*Tr(hb*hbd)-4*Tr(he*hed)-24*Tr(Ybd*Yb)*mHd-12*Tr(Ybd*Yb*mQ)-12*Tr(Yb*Ybd*mb)-8*Tr(Yed*Ye)*mHd-4*Tr(Yed*Ye*mL)-4*Tr(Ye*Yed*me)) +
 	Ye*mL*Yed*g1*g1*ir(-12,5) + Ye*mL*Yed*g2*g2*(12) + Ye*mL*Yed*Ye*Yed*(-4) + Ye*mL*Yed*(-12*Tr(Ybd*Yb)-4*Tr(Yed*Ye)) + Ye*(-12*hed*Tr(hb*Ybd)-4*hed*Tr(he*Yed)) +
	me*Ye*Yed*g1*g1*ir(-6,5) + me*Ye*Yed*g2*g2*(6) + me*Ye*Yed*Ye*Yed*(-2) + me*Ye*Yed*(-6*Tr(Ybd*Yb)-2*Tr(Yed*Ye));
MSSM1.push_back(bme_1);
MSSM2.push_back(bme_2);

BetaFunc bmL_1(mL);
BetaFunc bmL_2(mL);
bmL_1 = g1*g1*M1*M1*ir(-6,5) + g2*g2*M2*M2*(-6) + hed*he*(2) + Yed*Ye*mL*(1) + Yed*Ye*(2*mHd) + Yed*me*Ye*(2) + mL*Yed*Ye*(1) + ir(-3,5)*g1*g1*S;
bmL_2 = 	g1*g1*g1*g1*M1*M1*ir(621,25) + g1*g1*g1*g1*(ir(24,25)*Tr(mt) + ir(6,25)*Tr(mb) + ir(3,25)*Tr(mQ) + ir(18,25)*Tr(me) + ir(9,25)*Tr(mL) + ir(9,25)*mHd + ir(9,25)*mHu) +
	g1*g1*g2*g2*M1*M1*ir(18,5) + g1*g1*g2*g2*M1*M2*ir(18,5) + g1*g1*g2*g2*M2*M2*ir(18,5) + g2*g2*g2*g2*M2*M2*(33) + g2*g2*g2*g2*(9*Tr(mQ) + 3*Tr(mL) + 3*mHd + 3*mHu) +
	hed*he*g1*g1*ir(12,5) + hed*he*Yed*Ye*(-4) + hed*he*(-6*Tr(Ybd*Yb)-2*Tr(Yed*Ye)) + hed*Ye*g1*g1*M1*ir(-12,5) + hed*Ye*Yed*he*(-4) + hed*Ye*(-6*Tr(hb*Ybd)-2*Tr(he*Yed)) +
	Yed*he*g1*g1*M1*ir(-12,5) + Yed*he*hed*Ye*(-4) + Yed*he*(-6*Tr(hbd*Yb)-2*Tr(hed*Ye)) + Yed*Ye*g1*g1*M1*M1*ir(24,5) + Yed*Ye*g1*g1*(ir(12,5)*mHd) + Yed*Ye*hed*he*(-4) +
	Yed*Ye*Yed*Ye*mL*(-2) + Yed*Ye*Yed*Ye*(-8*mHd) + Yed*Ye*Yed*me*Ye*(-4) + Yed*Ye*mL*g1*g1*ir(6,5) + Yed*Ye*mL*Yed*Ye*(-4) + Yed*Ye*mL*(-3*Tr(Ybd*Yb)-Tr(Yed*Ye)) +
	Yed*Ye*(-6*Tr(hb*hbd)-2*Tr(he*hed)-12*Tr(Ybd*Yb)*mHd-6*Tr(Ybd*Yb*mQ)-6*Tr(Yb*Ybd*mb)-4*Tr(Yed*Ye)*mHd-2*Tr(Yed*Ye*mL)-2*Tr(Ye*Yed*me)) + Yed*me*Ye*g1*g1*ir(12,5) +
	Yed*me*Ye*Yed*Ye*(-4) + Yed*me*Ye*(-6*Tr(Ybd*Yb)-2*Tr(Yed*Ye)) + mL*Yed*Ye*g1*g1*ir(6,5) + mL*Yed*Ye*Yed*Ye*(-2) + mL*Yed*Ye*(-3*Tr(Ybd*Yb)-Tr(Yed*Ye));
MSSM1.push_back(bmL_1);
MSSM2.push_back(bmL_2);

BetaFunc bmu_1(mu);
BetaFunc bmu_2(mu);
bmu_1 =	mu*Tr(3*Yt*Ytd + 3*Yb*Ybd + Ye*Yed) - 3*mu*g2*g2 - ir(3,5)*mu*g1*g1;
bmu_2 =	mu*Tr(-3*Yt*Ytd*Yt*Ytd - 3*Yb*Ybd*Yb*Ybd - 2*Yt*Ytd*Yb*Ybd - Ye*Yed*Ye*Yed) + mu*Tr(Yt*Ytd)*(16*g3*g3 + ir(4,5)*g1*g1) + mu*Tr(Yb*Ybd)*(16*g3*g3 - ir(2,5)*g1*g1) +
	mu*ir(6,5)*g1*g1*Tr(Ye*Yed) + mu*ir(15,2)*g2*g2*g2*g2 + mu*ir(9,5)*g1*g1*g2*g2 + mu*ir(207,50)*g1*g1*g1*g1;
MSSM1.push_back(bmu_1);
MSSM2.push_back(bmu_2);

BetaFunc bB_1(B);
BetaFunc bB_2(B);
bB_1 = B*Tr(3*Yt*Ytd + 3*Yb*Ybd + Ye*Yed) - 3*B*g2*g2 - ir(3,5)*B*g1*g1 +
	mu*Tr(6*ht*Ytd + 6*hb*Ybd + 2*he*Yed) + 6*mu*M2*g2*g2 + ir(6,5)*mu*M1*g1*g1;
bB_2 =  B*Tr(-3*Yt*Ytd*Yt*Ytd - 3*Yb*Ybd*Yb*Ybd - 2*Yt*Ytd*Yb*Ybd - Ye*Yed*Ye*Yed) + B*Tr(Yt*Ytd)*(16*g3*g3 + ir(4,5)*g1*g1) + B*Tr(Yb*Ybd)*(16*g3*g3 - ir(2,5)*g1*g1) +
	B*ir(6,5)*g1*g1*Tr(Ye*Yed) + B*ir(15,2)*g2*g2*g2*g2 + B*ir(9,5)*g1*g1*g2*g2 + B*ir(207,50)*g1*g1*g1*g1 +
	12*mu*Tr(-3*ht*Ytd*Yt*Ytd - 3*hb*Ybd*Yb*Ybd - ht*Ybd*Yb*Ytd - hb*Ytd*Yt*Ybd - he*Yed*Ye*Yed) + mu*Tr(ht*Ytd)*(32*g3*g3 + ir(8,5)*g1*g1) + mu*Tr(hb*Ybd)*(32*g3*g3 - ir(4,5)*g1*g1) +
	ir(12,5)*mu*g1*g1*Tr(he*Yed) - mu*Tr(Yt*Ytd)*(32*g3*g3*M3 + ir(8,5)*g1*g1*M1) - mu*Tr(Yb*Ybd)*(32*g3*g3*M3 - ir(4,5)*g1*g1*M1) - ir(12,5)*mu*g1*g1*M1*Tr(Ye*Yed) - 30*mu*g2*g2*g2*g2*M2 -
	ir(18,5)*mu*g1*g1*g2*g2*M1 + ir(18,5)*mu*g1*g1*g2*g2*M2 - ir(414,25)*g1*g1*g1*g1*M1;
MSSM1.push_back(bB_1);
MSSM2.push_back(bB_2);

BetaFunc bBc_1 = Conjugate(bB_1);
BetaFunc bBc_2 = Conjugate(bB_2);
MSSM1.push_back(bBc_1);
MSSM2.push_back(bB_2);


findInvariants(MSSM1, MSSM2);

//--------------------------------------------------------------------------Test System----------------------------------------------------------
/*
Param x("x", 1, 0, 1);
Param y("y", 1, 0, 2);

vector<BetaFunc> test_1;
vector<BetaFunc> test_2;

BetaFunc bx_1(x);
BetaFunc bx_2(x);
bx_1 = 3*x*x + 3*y;
bx_2 = -1*x*x*x+5*x*y;
test_1.push_back(bx_1);
test_2.push_back(bx_2);

BetaFunc by_1(y);
BetaFunc by_2(y);
by_1 = -2*x*x*x - 2*x*y;
by_2 = 2*x*x*x*x - 4*x*x*y-2*y*y;
test_1.push_back(by_1);
test_2.push_back(by_2);

findInvariants(test_1, test_2);
*/

//--------------------------------------------------------------------------Test System 2----------------------------------------------------------
/*
Param x("x", 1, 0, 1);
Param y("y", 1, 0, 1);
Param z("z", 1, 0, 2);

vector<BetaFunc> test_1;
vector<BetaFunc> test_2;

BetaFunc bx_1(x);
BetaFunc bx_2(x);
bx_1 = -2*x*x + 2*y*y;
bx_2 = -1*x*x*x - 1*x*y*y;
test_1.push_back(bx_1);
test_2.push_back(bx_2);

BetaFunc by_1(y);
BetaFunc by_2(y);
by_1 = 1*x*x - 1*y*y;
by_2 = 1*x*y*y;
test_1.push_back(by_1);
test_2.push_back(by_2);

BetaFunc bz_1(z);
BetaFunc bz_2(z);
bz_1 = -5*x*x*x + 2*x*x*y + 5*x*y*y - 2*y*y*y;
bz_2 = 2*x*x*x*x + 13*y*y*y*y;
test_1.push_back(bz_1);
test_2.push_back(bz_2);

findInvariants(test_1, test_2);
*/
}
