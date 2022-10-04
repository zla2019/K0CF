#include <iostream>
#include "utils.h"

double CFLL(double* x, double* par);
std::complex<double> amplitude(double x);
std::complex<double> f0(double x);
std::complex<double> f1(double x);
double F1(double z);
double F2(double x);
double f0refMom(double x);
double a0refMom(double x);

const double hbarc = 0.197327;
const double MassK0s = 0.497611;
const double MassPi = 0.13957039;
const double MassEta = 0.547862;

////Antonelli
//const double MassF0 = 0.973;
//const double GammaF0KK = 2.763;
//const double GammaF0PiPi = 0.5283;
//const double MassA0 = 0.985;
//const double GammaA0KK = 0.4038;
//const double GammaA0PiPi = 0.3711;
////Achasov2001
//const double MassF0 = 0.996;
//const double GammaF0KK = 1.305;
//const double GammaF0PiPi = 0.2684;
//const double MassA0 = 0.992;
//const double GammaA0KK = 0.5555;
//const double GammaA0PiPi = 0.4401;
////Achasov2003
//const double MassF0 = 0.996;
//const double GammaF0KK = 1.305;
//const double GammaF0PiPi = 0.2684;
//const double MassA0 = 1.003;
//const double GammaA0KK = 0.8365;
//const double GammaA0PiPi = 0.4580;
//Martin
const double MassF0 = 0.978;
const double GammaF0KK = 0.792;
const double GammaF0PiPi = 0.1990;
const double MassA0 = 0.974;
const double GammaA0KK = 0.3330;
const double GammaA0PiPi = 0.2220;

void CFLLFunc()
{
	TF1* fCFMidLL;

	fCFMidLL = new TF1(Form("fCFMidLL"), CFLL, 0, 0.4, 2);
	fCFMidLL->SetParameter(0, 0.771);	//LL par.
	fCFMidLL->SetParameter(1, 4.6);
	//fCFMidLL->SetParameter(0, 0.802);	//gaussian par.
	//fCFMidLL->SetParameter(1, 5.84);
	std::cout << fCFMidLL->Eval(0.05) << std::endl;

	TCanvas* caCFMid = new TCanvas("caCFMid", "caCFMid", 640, 820);
	fCFMidLL->SetLineColor(kBlue);
	fCFMidLL->SetMaximum(2.5);
	fCFMidLL->SetMinimum(0.8);
	fCFMidLL->Draw();
}

double CFLLSI(double* x, double* par)
{
	double tKstar = 0.5 * x[0] / hbarc;
	double fsiTerm1 = 0.5 * par[0] * pow((std::abs(amplitude(tKstar * hbarc)) / par[1]), 2);
	double fsiTerm2 = 0.5 * par[0] * 4 / sqrt(TMath::Pi()) / par[1] * amplitude(tKstar * hbarc).real() * F1(2. * tKstar * par[1]);
	double fsiTerm3 = - 0.5 * par[0] * 2 / par[1] * amplitude(tKstar * hbarc).imag() * F2(2. * tKstar * par[1]);
	return fsiTerm1 + fsiTerm2 + fsiTerm3 + 1;
}

double CFLL(double* x, double* par)
{
	double qsTerm = par[0] * (TMath::Exp(-1 * par[1]*par[1] * x[0]*x[0] / (hbarc*hbarc)));
	double tKstar = 0.5 * x[0] / hbarc;
	double fsiTerm1 = 0.5 * par[0] * pow(std::abs(amplitude(tKstar * hbarc) / par[1]), 2);
	double fsiTerm2 = 0.5 * par[0] * 4 / sqrt(TMath::Pi()) / par[1] * amplitude(tKstar * hbarc).real() * F1(2. * tKstar * par[1]);
	double fsiTerm3 = - 0.5 * par[0] * 2 / par[1] * amplitude(tKstar * hbarc).imag() * F2(2. * tKstar * par[1]);
	return qsTerm + fsiTerm1 + fsiTerm2 + fsiTerm3 + 1;
}

std::complex<double> amplitude(double tKstar)
{
	std::complex<double> ImI = std::complex<double>(0, 1.);
	//return pow((1 / 17. + 0.5 * 2.7 * tKstar*tKstar - ImI*tKstar), -1);
	return 0.5 * (f0(tKstar) + f1(tKstar)) * hbarc;
}

std::complex<double> f0(double tKstar)	//f0(980)
{
	double s = 4 * (MassK0s*MassK0s + tKstar*tKstar);
	std::complex<double> a = std::complex<double>(0, GammaF0KK * tKstar);
	//std::complex<double> b = std::complex<double>(0, 0.05 * 0.973);
	std::complex<double> b = std::complex<double>(0, GammaF0PiPi * f0refMom(tKstar));
	//std::complex<double> b = std::complex<double>(0, 1 * 0.973);
	//std::complex<double> b = std::complex<double>(0, 5 * 0.5283);
	std::complex<double> result = GammaF0KK / (MassF0*MassF0 - s - a - b);
	return result;
}

std::complex<double> f1(double tKstar)	//a0(980)
{
	double s = 4 * (MassK0s*MassK0s + tKstar*tKstar);
	std::complex<double> a = std::complex<double>(0, GammaA0KK * tKstar);
	//std::complex<double> b = std::complex<double>(0, 0.06 * 0.985);
	std::complex<double> b = std::complex<double>(0, GammaA0PiPi * a0refMom(tKstar));
	//std::complex<double> b = std::complex<double>(0, 5 * 0.4038);
	std::complex<double> result = GammaA0KK / (MassA0*MassA0 - s - a - b);
	return result;
}

double F1(double z)
{
	TF1* fTmp = new TF1("fTmp", "exp(x*x - [0]*[0]) / [0]", 0, z);
	fTmp->SetParameter(0, z);
	double result = fTmp->Integral(0, z);
	delete fTmp;
	return result;
}

double F2(double x)
{
	return (1 - TMath::Exp(-1 * x*x)) / x;
}

double f0refMom(double tKstar)
{
	return sqrt(tKstar*tKstar + MassK0s*MassK0s - MassPi*MassPi);
}

double a0refMom(double tKstar)
{
	double s = 4 * (MassK0s*MassK0s + tKstar*tKstar);
	return sqrt(pow(MassPi, 4) + pow(MassEta, 4) + s*s - 2 * (MassPi*MassPi * MassEta*MassEta + MassPi*MassPi * s + MassEta*MassEta * s)) / (2 * sqrt(s));
}
