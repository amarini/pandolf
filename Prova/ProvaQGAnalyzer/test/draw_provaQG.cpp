#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLegend.h"
#include "DrawBase.h"


void drawOne( DrawBase* db, const std::string& suffix );



int main() {

  DrawBase* db = new DrawBase("provaQG");
  db->set_outputdir("prova");

  db->set_rebin(2);

  TFile* file = TFile::Open("provaQG.root");
  db->add_mcFile(file, "prova", "prova");

  drawOne( db, "quark_pt3050" );
  drawOne( db, "gluon_pt3050" );

  drawOne( db, "quark_pt80120" );
  drawOne( db, "gluon_pt80120" );

  drawOne( db, "quark_pt200300" );
  drawOne( db, "gluon_pt200300" );

  return 0;

}


void drawOne( DrawBase* db, const std::string& suffix ) {

  std::vector< HistoAndName > hn;

  HistoAndName hn_old;
  hn_old.histoName = "qglOLD_" + suffix;
  hn_old.legendName = "Histogram based LD";
  hn.push_back(hn_old);

  HistoAndName hn_new;
  hn_new.histoName = "qglNEW_" + suffix;
  hn_new.legendName = "Fit based LD";
  hn_new.markerStyle = 20;
  hn.push_back(hn_new);

  db->compareDifferentHistos( hn, suffix, "QG LD" );

}
