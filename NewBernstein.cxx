
#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>
#include "TMath.h"
#include "NewBernstein.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"

using namespace std;

ClassImp(NewBernstein)
;


//_____________________________________________________________________________
NewBernstein::NewBernstein()
{
}


//_____________________________________________________________________________
NewBernstein::NewBernstein(const char* name, const char* title,
                           RooAbsReal& x, const RooArgList& coefList):
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefficients","List of coefficients",this)
{
  // Constructor
  TIterator* coefIter = coefList.createIterator() ;
  RooAbsArg* coef ;
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      cout << "NewBernstein::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName()
	   << " is not of type RooAbsReal" << endl ;
      assert(0) ;
    }
    _coefList.add(*coef) ;
  }
  delete coefIter ;
  _xmin = _x.min();
  _xmax = _x.max();
}



//_____________________________________________________________________________
NewBernstein::NewBernstein(const NewBernstein& other, const char* name) :
  RooAbsPdf(other, name),
  _x("x", this, other._x),
  _coefList("coefList",this,other._coefList)
{
  _xmin = other._xmin;
  _xmax = other._xmax;
}


//_____________________________________________________________________________
Double_t NewBernstein::evaluate() const
{

  //Double_t xmin = _x.min();
  Double_t xmin = _xmin;
  Double_t xmax = _xmax;
  //Double_t x = (_x - xmin) / (_x.max() - xmin); // rescale to [0,1]
  Double_t x = (_x - xmin) / (xmax - xmin);
  Int_t degree = _coefList.getSize() - 1; // n+1 polys of degree n
  RooFIter iter = _coefList.fwdIterator();

  if(degree == 0) {

    return ((RooAbsReal *)iter.next())->getVal();

  } else if(degree == 1) {

    Double_t a0 = ((RooAbsReal *)iter.next())->getVal(); // c0
    Double_t a1 = ((RooAbsReal *)iter.next())->getVal() - a0; // c1 - c0
    return a1 * x + a0;

  } else if(degree == 2) {

    Double_t a0 = ((RooAbsReal *)iter.next())->getVal(); // c0
    Double_t a1 = 2 * (((RooAbsReal *)iter.next())->getVal() - a0); // 2 * (c1 - c0)
    Double_t a2 = ((RooAbsReal *)iter.next())->getVal() - a1 - a0; // c0 - 2 * c1 + c2
    return (a2 * x + a1) * x + a0;

  } else if(degree > 2) {

    Double_t t = x;
    Double_t s = 1 - x;

    Double_t result = ((RooAbsReal *)iter.next())->getVal() * s;
    for(Int_t i = 1; i < degree; i++) {
      result = (result + t * TMath::Binomial(degree, i) * ((RooAbsReal *)iter.next())->getVal()) * s;
      t *= x;
    }
    result += t * ((RooAbsReal *)iter.next())->getVal();

    return result;
  }

  // in case list of arguments passed is empty
  return TMath::SignalingNaN();
}


//_____________________________________________________________________________
Int_t NewBernstein::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  // No analytical calculation available (yet) of integrals over subranges
  // if (rangeName && strlen(rangeName)) {
  //   return 0 ;
  // }

  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}


//_____________________________________________________________________________
Double_t NewBernstein::analyticalIntegral(Int_t code, const char* rangeName) const
{
  (void)code;
  assert(code==1) ;
  const Double_t xminfull(_xmin), xmaxfull(_xmax);
  Double_t xmin = _x.min(rangeName); Double_t xmax = _x.max(rangeName);
  
  Double_t norm(0) ;

  Double_t fullRange = xmaxfull - xminfull;

  Double_t minScaled = (xmin - xminfull) / fullRange;
  Double_t maxScaled = 1.- (xmaxfull - xmax) / fullRange;

  double val =  fullRange * (evalAnaInt(maxScaled) - evalAnaInt(minScaled));
  return val;
}

Double_t NewBernstein::evalAnaInt(const Double_t x) const {
  Int_t degree= _coefList.getSize()-1; // n+1 polys of degree n
  RooFIter iter = _coefList.fwdIterator() ;
  Double_t temp=0, val=0;
  for (int i=0; i<=degree; ++i){
    // for each of the i Bernstein basis polynomials
    // represent it in the 'power basis' (the naive polynomial basis)
    // where the integral is straight forward.
    temp = 0;
    for (int j=i; j<=degree; ++j){ // power basis≈ß
      temp += pow(-1.,j-i) * TMath::Binomial(degree, j) * TMath::Binomial(j,i) * pow(x,j+1) / (j+1);
    }
    temp *= ((RooAbsReal*)iter.next())->getVal(); // include coeff
    val += temp; // add this basis's contribution to total
  }
  return val;
}