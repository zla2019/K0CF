#include <iostream>
#include <string.h>
#ifndef _UTILS_H_
#define _UTILS_H_

void			shiftGr(TGraphErrors* gr, float shift = 0.2);
double			getYield(TGraphErrors *gr, TF1 *f, double& err);	//get yield from the dN/dpt/dy spectrum, and unmeaured part calculate by TF1 function
double			getYield(TGraphErrors *gr, TH1F *h, double& err);
TGraphErrors*		timePt2(TGraphErrors* gr);
void 			timePt(TGraphErrors* gr);
void 			timePt(TH1F* h);
void 			dividePt(TGraphErrors* gr);
TH1F* 			gr2h(TGraphErrors* gr);
TGraphErrors* 		h2gr(TH1F* h);
TH1F* 			f2h(TF1* f, float binWidth, float lower, float upper);
TH1F* 			cutHist(float lowerEdge = 0, float upperEdge = 1., TH1F* h = NULL);
TH1F* 			cutHist2(float lowerEdge = 0, float upperEdge = 1., TH1F* h = NULL);
double 			combineIntegral(TH1F* h, TF1* f, float& err, float lower = 0, float upper = 1);
TObject* 		getCopy(TFile* input, const char* name);
TH1F* 			projectionZ(TH3F* _3h, float lowerX, float upperX, float lowerY, float upperY, const char* name = "");
TH1F* 			projectionX(TH3F* _3h, float lowerY, float upperY, float lowerZ, float upperZ, const char* name = "");
TH1F* 			projectionY(TH3F* _3h, float lowerX, float upperX, float lowerZ, float upperZ, const char* name = "");
TH1F* 			projectionX(TH2F* _2h, float lowerY, float upperY, const char* name = "");
TH1F* 			projectionY(TH2F* _2h, float lowerX, float upperX, const char* name = "");
TH1F*			extractSig(TH1F* tot, TH1F* bg, float normalLower, float normalUpper, const char* opt = "my");
			//extract signal by tot - bg inv mass plot
TH1F*			extractResidual(TH1F* sig, TF1* fitFunc, TF1* residualBg, int n, float lower, float upper);	//subtract by tf1 value  in each bin center
TH1F*			extractResidual2(TH1F* sig, TF1* fitFunc, TF1* residualBg, int n, float lower, float upper);	//subtract by integral tf1 in each bins range
void			setStyle(TH1F* h, int markerStyle, float markerSize, Color_t markerColor, int lineWidth = 1, Color_t lineColor = kBlack);
void			setStyle(TGraphErrors* gr, int markerStyle, float markerSize, Color_t markerColor, int lineWidth = 1, Color_t lineColor = kBlack);
void			drawText(float x, float y, float size, const char* txt1, const char* txt2 = NULL, const char* txt3 = NULL, const char* txt4 = NULL, const char* txt5 = NULL);
void			drawText(float x, float y, float size, Color_t color, const char* txt1, const char* txt2 = NULL, const char* txt3 = NULL, const char* txt4 = NULL, const char* txt5 = NULL);
void			drawXBaseLine(double xpoint, TPad* pad, Color_t color = kBlack, int lineStyle = 1, short lineWidth = 1);	//developing
void			drawYBaseLine(double ypoint, TPad* pad, Color_t color = kBlack, int lineStyle = 1, short lineWidth = 1);	//developing
void                    drawXRegion(double xpoint1, double xpoint2, TPad* pad, Color_t color = kRed, int lineStyle = 2, short lineWidth = 2);
void                    drawYRegion(double ypoint1, double ypoint2, TPad* pad, Color_t color = kRed, int lineStyle = 2, short lineWidth = 2);
float			getIntegral(TH1* h, float lower, float upper, const char* opt = "");
float			getIntegral(TH1* h, float lower, float upper, double& err, const char* opt = "");
void			setXYTitle(TH1* h, const char* xtitle, const char* ytitle);
void			removeNegtive(TH1F* h);
void			drawBox(float x1, float y1, float x2, float y2, int lineWidth = 2, int lineStyle = 2, Color_t lineColor = kRed);
TH1F*			foldHist(TH1F* h, float centralLine = 0);
void			addHist(TH1F*& h1, TH1F* h2, const char* name, float scale = 1, float baseline = 0);
void			addHist(TH1F*& h1, TH1F* h2, const char* name, TH1F* hScale, float baseline = 0);
void			addHist(TH3F*& _3h1, TH3F* _3h2, const char* name);
void			setMarker(TH1F* h1, int markerStyle = 20, float markerSize = 1, int markerColor = 1);
void			setMarker(TGraphErrors* gr, int markerStyle = 20, float markerSize = 1, int markerColor = 1);
void			setLine(TH1F* h1, int lineStyle = 1, int lineWidth = 1, int lineColor = 1);
float			getMax(TH1F* h1, TH1F* h2);
void			scale(TH1F* h1, float scale, float baseline = 0);
void			scale(TH1F* h1, TH1F* hScale, float baseline = 0);
void			divide(TH1F* h1, TH1F* h2, float baseline = 0, float baseline2 = 0);
void			scaleX(TH1F*& h1, float scale = 1);

void scaleX(TH1F*& h1, float scale)
{
	float xLower = h1->GetBinLowEdge(1);
	float xUpper = h1->GetBinLowEdge(h1->GetNbinsX()) + h1->GetBinWidth(h1->GetNbinsX());
	TH1F* h2 = new TH1F(h1->GetName(), h1->GetTitle(), h1->GetNbinsX(), xLower * scale, xUpper * scale);
	for(int ibin = 0; ibin < h1->GetNbinsX(); ++ibin) {
		h2->SetBinContent(ibin + 1, h1->GetBinContent(ibin + 1));
		h2->SetBinError(ibin + 1, h1->GetBinError(ibin + 1));
	}
	delete h1;
	h1 = NULL;
	h1 = h2;
}

void divide(TH1F* h1, TH1F* h2, float baseline, float baseline2)
{
	for(int ibin = 0; ibin < h1->GetNbinsX(); ++ibin) {
		float binContent = h1->GetBinContent(ibin + 1) - baseline;
		float binError = h1->GetBinError(ibin + 1);
		float scale = h2->GetBinContent(ibin + 1) - baseline2;
		float scaleError = h2->GetBinError(ibin + 1);
		binError = binContent / scale * sqrt(pow(binError / binContent, 2) + pow(scaleError / scale, 2));
		binContent /= scale;
		h1->SetBinContent(ibin + 1, binContent + baseline);
		h1->SetBinError(ibin + 1, binError);
	}
}

void scale(TH1F* h1, TH1F* hScale, float baseline)
{
	for(int ibin = 0; ibin < h1->GetNbinsX(); ++ibin) {
		float binContent = h1->GetBinContent(ibin + 1) - baseline;
		float binError = h1->GetBinError(ibin + 1);
		float scale = hScale->GetBinContent(ibin + 1);
		float scaleError = hScale->GetBinError(ibin + 1);
		binError = binContent * scale * sqrt(pow(binError / binContent, 2) + pow(scaleError / scale, 2));
		binContent *= scale;
		h1->SetBinContent(ibin + 1, binContent + baseline);
		h1->SetBinError(ibin + 1, binError);
	}
}

void scale(TH1F* h1, float scale, float baseline)
{
	for(int ibin = 0; ibin < h1->GetNbinsX(); ++ibin) {
		float binContent = h1->GetBinContent(ibin + 1);
		float binError = h1->GetBinError(ibin + 1);
		binContent = (binContent - baseline) * scale + baseline;
		h1->SetBinContent(ibin + 1, binContent);
		h1->SetBinError(ibin + 1, binError);
	}
}

void addHist(TH3F*& _3h1, TH3F* _3h2, const char* name)
{
	if(_3h1) {
		_3h1->Add(_3h2);
	} else if(!_3h1) {
		_3h1 = (TH3F*)_3h2->Clone();
		_3h1->SetName(name);
	}
}

void addHist(TH1F*& h1, TH1F* h2, const char* name, TH1F* hScale, float baseline)
{
	if(h1) {
		for(int ibin = 0; ibin < h1->GetNbinsX(); ++ibin) {
			float binContent = h1->GetBinContent(ibin + 1) - baseline;
			float binError = h1->GetBinError(ibin + 1);
			float binContent2 = h2->GetBinContent(ibin + 1) - baseline;
			float binError2 = h2->GetBinError(ibin + 1);
			float weight = hScale->GetBinContent(ibin + 1);
			float weightError = hScale->GetBinError(ibin + 1);
			float totBinContent = binContent + binContent2 * weight;
			float totBinError = binContent2 * weight * sqrt(pow(binError2 / binContent2, 2) + pow(weightError / weight, 2));
			totBinError = sqrt(binError*binError + totBinError*totBinError);
			h1->SetBinContent(ibin + 1, totBinContent + baseline);
			h1->SetBinError(ibin + 1, totBinError);
		}
	} else if(!h1) {
		h1 = (TH1F*)h2->Clone();
		h1->SetName(name);
		for(int ibin = 0; ibin < h1->GetNbinsX(); ++ibin) {
			float binContent = h2->GetBinContent(ibin + 1) - baseline;
			float binError = h2->GetBinError(ibin + 1);
			float weight = hScale->GetBinContent(ibin + 1);
			float weightError = hScale->GetBinError(ibin + 1);
			float totBinContent = binContent * weight;
			float totBinError = totBinContent * sqrt(pow(binError / binContent, 2) + pow(weightError / weight, 2));
			h1->SetBinContent(ibin + 1, totBinContent + baseline);
			h1->SetBinError(ibin + 1, totBinError);
		}
	}
}

float getMax(TH1F* h1, TH1F* h2)
{
	return h1->GetMaximum() > h2->GetMaximum() ? h1->GetMaximum() : h2->GetMaximum();
}

void setLine(TH1F* h1, int lineStyle, int lineWidth, int lineColor)
{
	h1->SetLineStyle(lineStyle);
	h1->SetLineWidth(lineWidth);
	h1->SetLineColor(lineColor);
}

void setMarker(TH1F* h1, int markerStyle, float markerSize, int markerColor)
{
	h1->SetMarkerStyle(markerStyle);
	h1->SetMarkerSize(markerSize);
	h1->SetMarkerColor(markerColor);
}

void setMarker(TGraphErrors* gr, int markerStyle, float markerSize, int markerColor)
{
	gr->SetMarkerStyle(markerStyle);
	gr->SetMarkerSize(markerSize);
	gr->SetMarkerColor(markerColor);
}

void addHist(TH1F*& h1, TH1F* h2, const char* name, float scale, float baseline)
{
	if(h1) {
		//h1->Add(h2, scale);
		for(int ibin = 0; ibin < h1->GetNbinsX(); ++ibin) {
			float binContent = h1->GetBinContent(ibin + 1) - baseline;
			float binError = h1->GetBinError(ibin + 1);
			float binContent2 = h2->GetBinContent(ibin + 1) - baseline;
			float binError2 = h2->GetBinError(ibin + 1);
			float totBinContent = binContent + binContent2 * scale + baseline;
			float totBinError = sqrt(binError*binError + binError2*binError2);
			h1->SetBinContent(ibin + 1, totBinContent);
			h1->SetBinError(ibin + 1, totBinError);
		}
	} else {
		h1 = (TH1F*)h2->Clone();
		h1->SetName(name);
		for(int ibin = 0; ibin < h1->GetNbinsX(); ++ibin) {
			float binContent = h2->GetBinContent(ibin + 1) - baseline;
			float binError = h2->GetBinError(ibin + 1);
			binContent *= scale;
			h1->SetBinContent(ibin + 1, binContent + baseline);
			h1->SetBinError(ibin + 1, binError);
		}
	}
}

TH1F* foldHist(TH1F* h, float centralLine) {
	int centralBin = h->FindBin(centralLine + 1e-4);
	int nBins = h->GetNbinsX() - centralBin + 1;
	float lowerEdge = centralLine, upperEdge = h->GetBinLowEdge(h->GetNbinsX()) + h->GetBinWidth(h->GetNbinsX());
	TH1F* h2 = new TH1F(h->GetName(), h->GetTitle(), nBins, lowerEdge, upperEdge);
	for(int ibin = 0; ibin < nBins; ++ibin) {
		float binContent = h->GetBinContent(centralBin + ibin);
		binContent += h->GetBinContent(centralBin - 1 - ibin);
		float binError = h->GetBinError(centralBin + ibin);
		binError = sqrt(binError*binError + h->GetBinError(centralBin - 1 - (2 * ibin)));
		h2->SetBinContent(ibin + 1, binContent);
		h2->SetBinError(ibin + 1, binError);
	}
	delete h;
	return h2;
}

void drawBox(float x1, float y1, float x2, float y2, int lineWidth, int lineStyle, Color_t lineColor)
{
	TLine* line[4];
	line[0] = new TLine(x1, y1, x1, y2);
	line[1] = new TLine(x1, y2, x2, y2);
	line[2] = new TLine(x2, y2, x2, y1);
	line[3] = new TLine(x2, y1, x1, y1);
	for(int iedge = 0; iedge < 4; ++iedge) {
		line[iedge]->SetLineWidth(lineWidth);
		line[iedge]->SetLineStyle(lineStyle);
		line[iedge]->SetLineColor(lineColor);
		line[iedge]->Draw("same");
	}
}

void removeNegtive(TH1F* h)
{
	int nBins = h->GetNbinsX();
	for(int ibin = 0; ibin < nBins; ++ibin) {
		if(h->GetBinContent(ibin + 1) < 0) {
			h->SetBinContent(ibin + 1, 0);
		}
	}
}

void setXYTitle(TH1* h, const char* xtitle, const char* ytitle)
{
	h->GetXaxis()->SetTitle(xtitle);
	h->GetYaxis()->SetTitle(ytitle);
}

float getIntegral(TH1* h, float lower, float upper, double& err, const char* opt)
{
	int lowerBin = h->GetXaxis()->FindBin(lower + 0.000001);
	int upperBin = h->GetXaxis()->FindBin(upper - 0.000001);
	return h->IntegralAndError(lowerBin, upperBin, err, opt);
}

float getIntegral(TH1* h, float lower, float upper, const char* opt)
{
	int lowerBin = h->GetXaxis()->FindBin(lower + 0.000001);
	int upperBin = h->GetXaxis()->FindBin(upper - 0.000001);
	return h->Integral(lowerBin, upperBin, opt);
}

void drawYBaseLine(double ypoint, TPad* pad, Color_t color, int lineStyle, short lineWidth) //developing
{
        if(!pad) {
                return NULL;
        }
        pad->Update();
        pad->cd();
        float xLower = 0, xUpper = 0;
        if(pad->GetLogx() == 1) {
                xLower = pow(10, pad->GetUxmin());
                xUpper = pow(10, pad->GetUxmax());
        } else if(pad->GetLogx() == 0) {
                xLower = pad->GetUxmin();
                xUpper = pad->GetUxmax();
        }
        TLine* line = new TLine(xLower, ypoint, xUpper, ypoint);
        line->SetLineStyle(lineStyle);
        line->SetLineWidth(lineWidth);
        line->SetLineColor(color);
        line->Draw("same");
}	//developing

void drawXBaseLine(double xpoint, TPad* pad, Color_t color, int lineStyle, short lineWidth) //developing
{
        if(!pad) {
                return NULL;
        }
        pad->Update();
        pad->cd();
        float yLower = 0, yUpper = 0;
        if(pad->GetLogy() == 1) {
                yLower = pow(10, pad->GetUymin());
                yUpper = pow(10, pad->GetUymax());
        } else if(pad->GetLogy() == 0) {
                yLower = pad->GetUymin();
                yUpper = pad->GetUymax();
        }
        TLine* line = new TLine(xpoint, yLower, xpoint, yUpper);
        line->SetLineStyle(lineStyle);
        line->SetLineWidth(lineWidth);
        line->SetLineColor(color);
        line->Draw("same");
}	//developing

void drawXRegion(double xpoint1, double xpoint2, TPad* pad, Color_t color, int lineStyle, short lineWidth)
{
        drawXBaseLine(xpoint1, pad, color, lineStyle, lineWidth);
        drawXBaseLine(xpoint2, pad, color, lineStyle, lineWidth);
}

void drawYRegion(double ypoint1, double ypoint2, TPad* pad, Color_t color, int lineStyle, short lineWidth)
{
        drawYBaseLine(ypoint1, pad, color, lineStyle, lineWidth);
        drawYBaseLine(ypoint2, pad, color, lineStyle, lineWidth);
}

void drawText(float x, float y, float size, const char* txt1, const char* txt2, const char* txt3, const char* txt4, const char* txt5)
{
	TLatex* ltx = new TLatex(x, y, txt1);
	ltx->SetTextSize(size);
	ltx->SetNDC();
	ltx->DrawLatex(x, y, txt1);
	if(txt2) {
		ltx->DrawLatex(x, y - size * 1.5, txt2);
	}
	if(txt3) {
		ltx->DrawLatex(x, y - size * 1.5 * 2, txt3);
	}
	if(txt4) {
		ltx->DrawLatex(x, y - size * 1.5 * 3, txt4);
	}
	if(txt5) {
		ltx->DrawLatex(x, y - size * 1.5 * 4, txt5);
	}
}

void drawText(float x, float y, float size, Color_t color, const char* txt1, const char* txt2, const char* txt3, const char* txt4, const char* txt5)
{
	TLatex* ltx = new TLatex(x, y, txt1);
	ltx->SetTextSize(size);
	ltx->SetTextColor(color);
	ltx->SetNDC();
	ltx->DrawLatex(x, y, txt1);
	if(txt2) {
		ltx->DrawLatex(x, y - size * 1.5, txt2);
	}
	if(txt3) {
		ltx->DrawLatex(x, y - size * 1.5 * 2, txt3);
	}
	if(txt4) {
		ltx->DrawLatex(x, y - size * 1.5 * 3, txt4);
	}
	if(txt5) {
		ltx->DrawLatex(x, y - size * 1.5 * 4, txt5);
	}
}

void setStyle(TGraphErrors* gr, int markerStyle, float markerSize, Color_t markerColor, int lineWidth, Color_t lineColor)
{
	gr->SetMarkerStyle(markerStyle);
	gr->SetMarkerSize(markerSize);
	gr->SetMarkerColor(markerColor);
	gr->SetLineWidth(lineWidth);
	gr->SetLineColor(lineColor);
}

void setStyle(TH1F* h, int markerStyle, float markerSize, Color_t markerColor, int lineWidth, Color_t lineColor)
{
	h->SetMarkerStyle(markerStyle);
	h->SetMarkerSize(markerSize);
	h->SetMarkerColor(markerColor);
	h->SetLineWidth(lineWidth);
	h->SetLineColor(lineColor);
}

TH1F* extractResidual(TH1F* sig, TF1* fitFunc, TF1* residualBg, int n, float lower, float upper)
{
	sig->Fit(fitFunc, "MO0", "", lower, upper);
	TH1F* finalPlot = (TH1F*)sig->Clone();
	residualBg->SetParameters(&(fitFunc->GetParameters()[n]));
	finalPlot->Add(residualBg, -1);
	return finalPlot;
}

TH1F* extractResidual2(TH1F* sig, TF1* fitFunc, TF1* residualBg, int n, float lower, float upper)
{
	sig->Fit(fitFunc, "MO0", "", lower, upper);
	TH1F* finalPlot = (TH1F*)sig->Clone();
	residualBg->SetParameters(&(fitFunc->GetParameters()[n]));
	int nBins = finalPlot->GetNbinsX();
	for(int ibin = 0; ibin < nBins; ++ibin) {
		float histContent = sig->GetBinContent(ibin + 1);
		float funcContent = residualBg->Integral(sig->GetBinLowEdge(ibin + 1), sig->GetBinLowEdge(ibin + 1) + sig->GetBinWidth(ibin + 1));
		float binContent = histContent - funcContent;

		float histError = sig->GetBinError(ibin + 1);
		float funcError = sqrt(abs(funcContent / sig->GetBinWidth(ibin + 1)));
		float binError = sqrt(histError*histError + funcError*funcError);
		if(binContent < 0) {
			binContent = 0;
		}
		finalPlot->SetBinContent(ibin + 1, binContent);
		finalPlot->SetBinError(ibin + 1, binError);
	}
	return finalPlot;
}

TH1F* extractSig(TH1F* tot, TH1F* bg, float normalLower, float normalUpper, const char* opt)
{
	TH1F* sig = new TH1F(*tot);
	TH1F* bgTmp = new TH1F(*bg);
	double scaler = 0, sigCount = 0, bgCount = 0;
	int lowerNormalBin = sig->FindBin(normalLower + 0.000001);
	int upperNormalBin = sig->FindBin(normalUpper - 0.000001);
	sigCount = sig->Integral(lowerNormalBin, upperNormalBin);
	bgCount = bgTmp->Integral(lowerNormalBin, upperNormalBin);
	scaler = sigCount / bgCount;
	bgTmp->Scale(-scaler);
	if(!strcmp(opt, "my")) {
		int nBins = sig->GetNbinsX();
		for(int ibin = 0; ibin < nBins; ++ibin) {
			sig->AddBinContent(ibin + 1, bgTmp->GetBinContent(ibin + 1));
		}
	} else {
		sig->Add(bgTmp);
	}

	delete bgTmp;
	return sig;
}

TH1F* projectionZ(TH3F* _3h, float lowerX, float upperX, float lowerY, float upperY, const char* name)
{
	int lowerXbin = _3h->GetXaxis()->FindBin(lowerX + 0.000001);
	int upperXbin = _3h->GetXaxis()->FindBin(upperX - 0.000001);
	int lowerYbin = _3h->GetYaxis()->FindBin(lowerY + 0.000001);
	int upperYbin = _3h->GetYaxis()->FindBin(upperY - 0.000001);

	TH1F* h = (TH1F*)_3h->ProjectionZ(name, lowerXbin, upperXbin, lowerYbin, upperYbin);
	h->Sumw2();
	return h;
}

TH1F* projectionX(TH3F* _3h, float lowerY, float upperY, float lowerZ, float upperZ, const char* name)
{
	int lowerYbin = _3h->GetYaxis()->FindBin(lowerY + 0.000001);
	int upperYbin = _3h->GetYaxis()->FindBin(upperY - 0.000001);
	int lowerZbin = _3h->GetZaxis()->FindBin(lowerZ + 0.000001);
	int upperZbin = _3h->GetZaxis()->FindBin(upperZ - 0.000001);

	TH1F* h = (TH1F*)_3h->ProjectionX(name, lowerYbin, upperYbin, lowerZbin, upperZbin);
	h->Sumw2();
	return h;
}

TH1F* projectionY(TH3F* _3h, float lowerX, float upperX, float lowerZ, float upperZ, const char* name)
{
	int lowerXbin = _3h->GetXaxis()->FindBin(lowerX + 0.000001);
	int upperXbin = _3h->GetXaxis()->FindBin(upperX - 0.000001);
	int lowerZbin = _3h->GetZaxis()->FindBin(lowerZ + 0.000001);
	int upperZbin = _3h->GetZaxis()->FindBin(upperZ - 0.000001);

	TH1F* h = (TH1F*)_3h->ProjectionY(name, lowerXbin, upperXbin, lowerZbin, upperZbin);
	h->Sumw2();
	return h;
}

TH1F* projectionX(TH2F* _2h, float lowerY, float upperY, const char* name)
{
	int lowerYbin = _2h->GetYaxis()->FindBin(lowerY + 0.000001);
	int upperYbin = _2h->GetYaxis()->FindBin(upperY - 0.000001);

	TH1F* h = (TH1F*)_2h->ProjectionX(name, lowerYbin, upperYbin);
	h->Sumw2();
	return h;
}

TH1F* projectionY(TH2F* _2h, float lowerX, float upperX, const char* name)
{
	int lowerXbin = _2h->GetXaxis()->FindBin(lowerX + 0.000001);
	int upperXbin = _2h->GetXaxis()->FindBin(upperX - 0.000001);

	TH1F* h = (TH1F*)_2h->ProjectionY(name, lowerXbin, upperXbin);
	h->Sumw2();
	return h;
}

TObject* getCopy(TFile* input, const char* name)
{
	TObject* obj = (TObject*)input->Get(name);
	if(!obj) {
		std::cout << name << " not found" << std::endl;
	}
	TObject* cloneObj = obj->Clone();
	return cloneObj;
}

double getYield(TGraphErrors *gr, TF1 *f, double& err)
{
	int nPoint = gr->GetN();
	double histBoundary[2] = { gr->GetPointX(0) - gr->GetErrorX(0),  gr->GetPointX(nPoint - 1) + gr->GetErrorX(nPoint - 1) };
	double integralPart = f->Integral(0, histBoundary[0]) + f->Integral(histBoundary[1], 3.);
	double integralPartErr = sqrt(pow(f->IntegralError(0, histBoundary[0]), 2) + pow(f->IntegralError(histBoundary[1], 3.), 2));
	double histPart = 0;
	double histPartErr = 0;
	for(int ipoint = 0; ipoint < nPoint; ++ipoint) {
		histPart += (gr->GetPointY(ipoint) * 2 * gr->GetErrorX(ipoint));
		histPartErr = sqrt(histPartErr*histPartErr + pow(gr->GetErrorY(ipoint), 2));
	}
	err = sqrt(histPartErr*histPartErr + integralPartErr*integralPartErr);
	return histPart + integralPart;
}

double getYield(TGraphErrors *gr, TH1F *h, double& err)
{
	int nPoint = gr->GetN();
	double histBoundary[2] = { gr->GetPointX(0) - gr->GetErrorX(0),  gr->GetPointX(nPoint - 1) + gr->GetErrorX(nPoint - 1) };

	double grPart = 0;
	double grPartErr = 0;
	for(int ipoint = 0; ipoint < nPoint; ++ipoint) {
		grPart += (gr->GetPointY(ipoint) * 2 * gr->GetErrorX(ipoint));
		grPartErr = sqrt(grPartErr*grPartErr + pow(gr->GetErrorY(ipoint), 2));
	}

	double histPart = 0;
	double histPartErr = 0;
	int leftBin = h->FindBin(histBoundary[0] - 0.000001);
	int rightBin = h->FindBin(histBoundary[1] + 0.000001);
	histPart = h->Integral(0, leftBin, "width") + h->Integral(rightBin, h->GetNbinsX(), "width");

	err = grPartErr;

	return grPart + histPart;
}

TGraphErrors* timePt2(TGraphErrors* gr)
{
	int nPoint = gr->GetN();
	float pointX[nPoint];
	float pointXErr[nPoint];
	float pointY[nPoint];
	float pointYErr[nPoint];
	for(int ipoint = 0; ipoint < nPoint; ++ipoint) {
		pointX[ipoint] = gr->GetPointX(ipoint);
		pointXErr[ipoint] = gr->GetErrorX(ipoint);
		pointY[ipoint] = gr->GetPointY(ipoint) * pointX[ipoint];
		pointYErr[ipoint] = gr->GetErrorY(ipoint) * pointX[ipoint];
	}
	TGraphErrors* gr2 = new TGraphErrors(nPoint, pointX, pointY, pointXErr, pointYErr);
	return gr2;
}

void timePt(TGraphErrors* gr)
{
	int nPoint = gr->GetN();
	float pointX[nPoint];
	float pointXErr[nPoint];
	float pointY[nPoint];
	float pointYErr[nPoint];
	for(int ipoint = 0; ipoint < nPoint; ++ipoint) {
		pointX[ipoint] = gr->GetPointX(ipoint);
		pointXErr[ipoint] = gr->GetErrorX(ipoint);
		pointY[ipoint] = gr->GetPointY(ipoint) * pointX[ipoint];
		pointYErr[ipoint] = gr->GetErrorY(ipoint) * pointX[ipoint];

		gr->SetPoint(ipoint, pointX[ipoint], pointY[ipoint]);
		gr->SetPointError(ipoint, pointXErr[ipoint], pointYErr[ipoint]);
	}
}

void timePt(TH1F* h)
{
	int nBins = h->GetNbinsX();
	for(int ibin = 0; ibin < nBins; ++ibin) {
		float newBinContent = h->GetBinContent(ibin + 1) * h->GetBinCenter(ibin + 1);
		float newBinError = h->GetBinError(ibin + 1) * h->GetBinCenter(ibin + 1);
		h->SetBinContent(ibin + 1, newBinContent);
		h->SetBinError(ibin + 1, newBinError);
	}
}

void dividePt(TGraphErrors* gr)
{
	int nPoint = gr->GetN();
	float pointX[nPoint];
	float pointXErr[nPoint];
	float pointY[nPoint];
	float pointYErr[nPoint];
	for(int ipoint = 0; ipoint < nPoint; ++ipoint) {
		pointX[ipoint] = gr->GetPointX(ipoint);
		pointXErr[ipoint] = gr->GetErrorX(ipoint);
		pointY[ipoint] = gr->GetPointY(ipoint) / pointX[ipoint];
		pointYErr[ipoint] = gr->GetErrorY(ipoint) / pointX[ipoint];

		gr->SetPoint(ipoint, pointX[ipoint], pointY[ipoint]);
		gr->SetPointError(ipoint, pointXErr[ipoint], pointYErr[ipoint]);
	}
}

TH1F* gr2h(TGraphErrors* gr)
{
	int nBins = gr->GetN();
	float binBoundary[nBins + 1];
	binBoundary[0] = gr->GetPointX(0) - gr->GetErrorXlow(0);
	for(int ibin = 0; ibin < nBins; ++ibin) {
		binBoundary[ibin + 1] = gr->GetPointX(ibin) + gr->GetErrorXhigh(ibin);
	}
	TH1F* h = new TH1F(gr->GetName(), gr->GetTitle(), nBins, binBoundary);
	for(int ibin = 0; ibin < nBins; ++ibin) {
		h->SetBinContent(ibin + 1, gr->GetPointY(ibin));
		h->SetBinError(ibin + 1, gr->GetErrorY(ibin));
	}
	return h;
}

TGraphErrors* h2gr(TH1F* h)
{
	int nPoint = h->GetNbinsX();
	float pointX[nPoint];
	float pointXErr[nPoint];
	float pointY[nPoint];
	float pointYErr[nPoint];
	for(int ipoint = 0; ipoint < nPoint; ++ipoint) {
		pointX[ipoint] = h->GetBinCenter(ipoint + 1);
		pointXErr[ipoint] = h->GetBinWidth(ipoint + 1) / 2;
		pointY[ipoint] = h->GetBinContent(ipoint + 1);
		pointYErr[ipoint] = h->GetBinError(ipoint + 1);
	}
	TGraphErrors* gr = new TGraphErrors(nPoint, pointX, pointY, pointXErr, pointYErr);
	return gr;
}

TH1F* f2h(TF1* f, float binWidth, float lower, float upper)
{
	int nBins = (upper - lower) / binWidth;
	TH1F* h = new TH1F(Form("%s_hist", f->GetName()), Form("%s_hist", f->GetName()), nBins, lower, upper);
	for(int ibin = 0; ibin < nBins; ++ibin) {
		h->SetBinContent(ibin + 1, f->Integral(lower + ibin * binWidth, lower + ibin * binWidth + binWidth) / binWidth);
		h->SetBinError(ibin + 1, f->IntegralError(lower + ibin * binWidth, lower + ibin * binWidth + binWidth) / binWidth);
	}
	return h;
}

TH1F* cutHist(float lowerEdge, float upperEdge, TH1F* h)
{
	if(!h) {
		std::cout << "ZLALOG >> missing target object" << std::endl;
		return NULL;
	}
	int lowerBin = h->GetXaxis()->FindBin(lowerEdge + 0.000001);
	int upperBin = h->GetXaxis()->FindBin(upperEdge - 0.000001);
	int nBins = upperBin - lowerBin + 1;
	float binBoundary[nBins + 1];
	for(int ibin = lowerBin; ibin <= upperBin + 1; ++ibin) {
		binBoundary[ibin - lowerBin] = h->GetBinLowEdge(ibin);
	}
	TH1F* h2 = new TH1F(h->GetName(), h->GetTitle(), nBins, binBoundary);
	for(int ibin = 0; ibin < nBins; ++ibin) {
		h2->SetBinContent(ibin + 1, h->GetBinContent(lowerBin + ibin));
		h2->SetBinError(ibin + 1, h->GetBinError(lowerBin + ibin));
	}
	delete h;
	h = NULL;
	return h2;
}

TH1F* cutHist2(float lowerEdge, float upperEdge, TH1F* h)
{
	if(!h) {
		std::cout << "ZLALOG >> missing target object" << std::endl;
		return NULL;
	}
	int lowerBin = h->GetXaxis()->FindBin(lowerEdge + 0.000001);
	int upperBin = h->GetXaxis()->FindBin(upperEdge - 0.000001);
	int nBins = upperBin - lowerBin + 1;
	TH1F* h2 = new TH1F(h->GetName(), h->GetTitle(), nBins, lowerEdge, upperEdge);
	for(int ibin = 0; ibin < nBins; ++ibin) {
		h2->SetBinContent(ibin + 1, h->GetBinContent(lowerBin + ibin));
		h2->SetBinError(ibin + 1, h->GetBinError(lowerBin + ibin));
	}
	delete h;
	h = NULL;
	return h2;
}

TH1F* reflectHist(TH1F* h, float centLine)
{
	int centBin = h->FindBin(centLine - 0.000001);
	int nBins = centBin * 2;
	float lowerEdge = h->GetBinLowEdge(1), upperEdge = 2 * centLine - lowerEdge;
	TH1F* hReflect = new TH1F("", "", nBins, lowerEdge, upperEdge);
	for(int ibin = 0; ibin < centBin; ++ibin) {
		hReflect->SetBinContent(ibin + 1, h->GetBinContent(ibin + 1));
		hReflect->SetBinError(ibin + 1, h->GetBinError(ibin + 1));
		hReflect->SetBinContent(2 * centBin - ibin, h->GetBinContent(ibin + 1));
		hReflect->SetBinError(2 * centBin - ibin, h->GetBinError(ibin + 1));
	}
	return hReflect;
}

double combineIntegral(TH1F* h, TF1* f, float& err, float lower, float upper)
{
	int histStartBin = 1, histEndBin = h->GetNbinsX();
	for(int ibin = 0; ibin < histEndBin; ++ibin) {
		if(h->GetBinContent(ibin + 1) == 0) {
			histStartBin = ibin + 1;
		} else {
			break;
		}
	}

	for(int ibin = h->GetNbinsX(); ibin > 0; --ibin) {
		if(h->GetBinContent(ibin) == 0) {
			histEndBin = ibin;
		} else {
			histEndBin--;
			break;
		}
	}

	float histStart = h->GetBinLowEdge(histStartBin), histEnd = h->GetBinLowEdge(histEndBin) + h->GetBinWidth(histEndBin);
	double errHistPart = 0, errFuncPart = 0;
	double integral = 0, integralHistPart = 0, integralFuncPart = 0;
	integralHistPart = h->IntegralAndError(histStartBin, histEndBin, errHistPart, "width");
	integralFuncPart = f->Integral(lower, histStart) + f->Integral(histEnd, upper);
	errFuncPart = sqrt(pow(f->IntegralError(lower, histStart), 2) + pow(f->IntegralError(histEnd, upper), 2));
	integral = integralHistPart + integralFuncPart;
	err = sqrt(errHistPart*errHistPart + errFuncPart*errFuncPart);
	return integral;
}

void shiftGr(TGraphErrors* gr, float shift)
{
	int nPoint = gr->GetN();
	for(int ipoint = 0; ipoint < nPoint; ++ipoint) {
		gr->SetPointX(ipoint, gr->GetPointX(ipoint) + shift);
	}
}
#endif
