#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;

class NewBernstein : public RooAbsPdf {
public:

  NewBernstein() ;
  NewBernstein(const char *name, const char *title,
               RooAbsReal& _x, const RooArgList& _coefList) ;

  NewBernstein(const NewBernstein& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new NewBernstein(*this, newname); }
  inline virtual ~NewBernstein() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

private:

  RooRealProxy _x;
  RooListProxy _coefList ;

  Double_t _xmin, _xmax;

  Double_t evaluate() const;
  Double_t evalAnaInt(const Double_t x) const;
  ClassDef(NewBernstein,1) // Bernstein polynomial PDF
};

#endif
