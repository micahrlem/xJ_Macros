#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

#include "xj_functions.h"
#include "read_binning.h"

#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TClass.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>


void Draw_Unfolded_Hists_Reweight_Test()

{

  /* read_binning rb(configfile.c_str());
  Int_t read_nbins = rb.get_nbins();
  const int nbins = read_nbins;
  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
    }

   float truth_leading_cut = rb.get_truth_leading_cut();
  float truth_subleading_cut = rb.get_truth_subleading_cut();

  float reco_leading_cut = rb.get_reco_leading_cut();
  float reco_subleading_cut = rb.get_reco_subleading_cut();

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();

  int truth_leading_bin = rb.get_truth_leading_bin();
  int truth_subleading_bin = rb.get_truth_subleading_bin();

  int reco_leading_bin = rb.get_reco_leading_bin();
  int reco_subleading_bin = rb.get_reco_subleading_bin();

  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();

   std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;*/

// List of all input ROOT files
std::vector<std::string> filenames = {
    // Pythia
    "UnfoldedxJHists_calofit_pythia8_R2.root",
    "UnfoldedxJHists_calofit_pythia8_R4.root",
    "UnfoldedxJHists_calofit_pythia8_negJER_R2.root",
    "UnfoldedxJHists_calofit_pythia8_negJER_R4.root",
    "UnfoldedxJHists_calofit_pythia8_posJER_R2.root",
    "UnfoldedxJHists_calofit_pythia8_posJER_R4.root",
    "UnfoldedxJHists_calofit_pythia8_negJES_R2.root",
    "UnfoldedxJHists_calofit_pythia8_negJES_R4.root",
    "UnfoldedxJHists_calofit_pythia8_posJES_R2.root",
    "UnfoldedxJHists_calofit_pythia8_posJES_R4.root",

    // Herwig Reweight (with test prefix)
    "UnfoldedxJHists_calofit_Herwig_Reweight_R2.root",
    "UnfoldedxJHists_calofit_Herwig_Reweight_R4.root",
    "UnfoldedxJHists_calofit_Herwig_Reweight_negJER_R2.root",
    "UnfoldedxJHists_calofit_Herwig_Reweight_negJER_R4.root",
    "UnfoldedxJHists_calofit_Herwig_Reweight_posJER_R2.root",
    "UnfoldedxJHists_calofit_Herwig_Reweight_posJER_R4.root",
    "UnfoldedxJHists_calofit_Herwig_Reweight_negJES_R2.root",
    "UnfoldedxJHists_calofit_Herwig_Reweight_negJES_R4.root",
    "UnfoldedxJHists_calofit_Herwig_Reweight_posJES_R2.root",
    "UnfoldedxJHists_calofit_Herwig_Reweight_posJES_R4.root"
};


// Output: map of histogram name → histogram pointer
std::map<std::string, TH1D*> named_histograms;

// Loop over all ROOT files
for (const auto& fname : filenames) {
    TFile* file = TFile::Open(fname.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open " << fname << std::endl;
        continue;
    }

    // Create a tag like "pythia8_R2" from the filename
    std::string tag = fname;
    tag.erase(0, std::string("UnfoldedxJHists_").length());
    size_t pos = tag.rfind(".root");
    if (pos != std::string::npos) tag.erase(pos);

    TIter next(file->GetListOfKeys());
    TKey* key;

    // Loop over all objects in the file
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (obj->IsA()->InheritsFrom("TH1D")) {
            std::string hname = obj->GetName();
            std::string fullkey = hname + "_" + tag;
            named_histograms[fullkey] = (TH1D*)obj;
        }
    }

   
 }
    
SetsPhenixStyle();
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

//Canvas making
 // Declare the canvas map
std::map<std::string, TCanvas*> canvases;

// Configuration
 
std::vector<std::string> closures = {"half", "full"};
std::vector<std::string> pt_ranges_closure = {"20_60"};
std::vector<std::string> pt_ranges_variation = {"20_30", "30_40", "40_60"};
std::vector<std::string> Rvals = {"R2", "R4"};
std::vector<std::string> gens = {"calofit_pythia8", "calofit_Herwig_Reweight"};

// 1. Half and Full Closure Tests (only for 20_60)
for (const auto& gen : gens) {
    for (const auto& R : Rvals) {
        std::string base = gen + "_" + R + "_";
        canvases["canvas_halfclosure_" + base + "20_60"] = new TCanvas(("canvas_halfclosure_" + base + "20_60").c_str(), "Half Closure", 700, 800);
        canvases["canvas_fullclosure_" + base + "20_60"] = new TCanvas(("canvas_fullclosure_" + base + "20_60").c_str(), "Full Closure", 700, 800);
    }
}

// 2. Base Unfolding, JER, JES Variations (for 20_30, 30_40, 40_60)
for (const auto& range : pt_ranges_variation) {
    for (const auto& gen : gens) {
        for (const auto& R : Rvals) {
            std::string tag = gen + "_" + R + "_" + range;
            canvases["canvas_base_" + tag] = new TCanvas(("canvas_base_" + tag).c_str(), "Base Unfolding", 700, 800);
            canvases["canvas_JER_" + tag]  = new TCanvas(("canvas_JER_" + tag).c_str(),  "JER Variation", 700, 800);
            canvases["canvas_JES_" + tag]  = new TCanvas(("canvas_JES_" + tag).c_str(),  "JES Variation", 700, 800);
        }
    }
}

 // Add to the same map from earlier

std::vector<std::string> pt_ranges_compare = {"20_30", "30_40", "40_60", "20_60"};

// Only one canvas per pt range, used for both R and generator comparisons
for (const auto& range : pt_ranges_compare) {
    std::string name = "canvas_compare_" + range;
    canvases[name] = new TCanvas(name.c_str(), ("Comparison canvas for " + range).c_str(), 700, 800);
}



 std::cout << "=== Canvas Keys ===" << std::endl;
for (const auto& [key, canvas] : canvases) {
    std::cout << key;
    if (!canvas) std::cout << "  [nullptr]";
    std::cout << std::endl;
}

std::cout << "=== Histogram Keys ===" << std::endl;
for (const auto& [key, hist] : named_histograms) {
    std::cout << key;
    if (!hist) std::cout << "  [nullptr]";
    std::cout << std::endl;
}

  std::map<std::string, std::string> pt_labels = {
    {"20_30", "20.9 < p_{T,1} < 31.2 GeV"},
    {"30_40", "31.2 < p_{T,1} < 40.7 GeV"},
    {"40_60", "40.7 < p_{T,1} < 60.8 GeV"}
};


// Comparison Plots


for (const auto& pt_range : pt_ranges_compare) {
    std::string canvas_key = "canvas_compare_" + pt_range;

    std::string tag_r4_pythia = pt_range + "_calofit_pythia8_R4";
    std::string tag_r4_herwig = pt_range + "_calofit_Herwig_Reweight_R4";
    std::string tag_r2_pythia = pt_range + "_calofit_pythia8_R2";
    std::string tag_r2_herwig = pt_range + "_calofit_Herwig_Reweight_R2";

    std::string h_r4_pythia = "h_xj_unfold_data_" + tag_r4_pythia;
    std::string h_r4_herwig = "h_xj_unfold_data_" + tag_r4_herwig;
    std::string h_r2_pythia = "h_xj_unfold_data_" + tag_r2_pythia;
    std::string h_r2_herwig = "h_xj_unfold_data_" + tag_r2_herwig;

    if (!canvases.count(canvas_key) || !canvases[canvas_key]) {
        std::cerr << "Missing canvas: " << canvas_key << std::endl;
        continue;
    }

    bool missing = false;
    if (!named_histograms.count(h_r4_pythia) || !named_histograms[h_r4_pythia]) {
        std::cerr << "Missing histogram: " << h_r4_pythia << std::endl;
        missing = true;
    }
    if (!named_histograms.count(h_r4_herwig) || !named_histograms[h_r4_herwig]) {
        std::cerr << "Missing histogram: " << h_r4_herwig << std::endl;
        missing = true;
    }
    if (!named_histograms.count(h_r2_pythia) || !named_histograms[h_r2_pythia]) {
        std::cerr << "Missing histogram: " << h_r2_pythia << std::endl;
        missing = true;
    }
    if (!named_histograms.count(h_r2_herwig) || !named_histograms[h_r2_herwig]) {
        std::cerr << "Missing histogram: " << h_r2_herwig << std::endl;
        missing = true;
    }
    if (missing) continue;

    xj_functions::plotxj_and_compare_new(
        canvases[canvas_key],
        named_histograms[h_r4_pythia],
        named_histograms[h_r4_herwig],
        named_histograms[h_r2_pythia],
        named_histograms[h_r2_herwig],
        "x_{J}", "PYTHIA/HERWIG",
        0, 1.0,
        1, 20, 2, 53, 6, 20, 9, 53);

    TLegend* leg_main = nullptr;
    std::string legend_key = "leg_compare_" + pt_range;
    xj_functions::draw_xj_legend(leg_main, legend_key.c_str(), 0.45, 0.75, 0.71, 0.90,
        named_histograms[h_r4_pythia],
        named_histograms[h_r4_herwig],
        named_histograms[h_r2_pythia],
        named_histograms[h_r2_herwig],
        "R=0.4 Data Unfolded w/ PYTHIA8", "R=0.4 Data Unfolded w/ HERWIG 7", "R=0.2 Data Unfolded w/ PYTHIA8", "R=0.2 Data Unfolded w/ HERWIG7", true);

    TLegend* leg_sphenix = nullptr;
    std::string leg_sphenix_key = "leg_sphenix_compare_" + pt_range;
    std::string pt_range_str = (pt_range == "20_30") ? "20.9 < p_{T,1} < 31.2 GeV" :
                               (pt_range == "30_40") ? "31.2 < p_{T,1} < 40.7 GeV" :
                               (pt_range == "40_60") ? "40.7 < p_{T,1} < 60.8 GeV" :
                                                        "20.9 < p_{T,1} < 60.8 GeV";

    xj_functions::DrawsPHENIXLegendNew(leg_sphenix, leg_sphenix_key.c_str(),
        0.13, 0.52, 0.43, 0.90, "p+p #sqrt{s}=200 GeV",
        "anti-#it{k}_{#it{t}} #it{R} = 0.2, 0.4", pt_range_str.c_str(), "Matched Truth Dijets",
        true, true, false, true, false);

    std::string outpath = "xj_plot_calofit_pdfs/" + canvas_key + ".pdf";
    canvases[canvas_key]->SaveAs(outpath.c_str());
    std::cout << "Saved: " << outpath << std::endl;
}



 

// Individual Plots


//Closure Plots
for (const auto& closure : closures) {
    for (const auto& pt_range : pt_ranges_closure) {
        for (const auto& gen : gens) {
            for (const auto& R : Rvals) {

                std::string tag = closure + "_" + pt_range + "_" + gen + "_" + R;
                std::string canvas_key = "canvas_" + closure + "closure_" + gen + "_" + R + "_" + pt_range;

                std::string h_reco   = "h_xj_reco_"   + tag;
                std::string h_unfold = "h_xj_unfold_" + tag;
                std::string h_truth  = "h_xj_truth_"  + tag;

                if (!canvases.count(canvas_key) || !canvases[canvas_key]) {
                    std::cerr << "Missing canvas: " << canvas_key << std::endl;
                    continue;
                }
                if (!named_histograms.count(h_reco) || !named_histograms[h_reco] ||
                    !named_histograms.count(h_unfold) || !named_histograms[h_unfold] ||
                    !named_histograms.count(h_truth) || !named_histograms[h_truth]) {
                    std::cerr << "Missing one or more histograms for: " << tag << std::endl;
                    continue;
                }

                xj_functions::plotxj_and_ratio(canvases[canvas_key],
                    named_histograms[h_reco],
                    named_histograms[h_unfold],
                    named_histograms[h_truth],
                    nullptr, false, "x_{J}", "Unfold/Truth",
                    0, 1.0, 4, 20, 4, 53, kGreen+1, 20, 8, 53);

                TLegend* leg_main = nullptr;
                std::string legend_key = "leg_" + canvas_key;
                std::string simtype     = (gen == "calofit_pythia8") ? "PYTHIA 8" : "HERWIG 7";
                std::string simtypeReco = simtype + " Reco";

                xj_functions::draw_xj_legend(leg_main, legend_key.c_str(), 0.61, 0.80, 0.71, 0.90,
                    named_histograms[h_truth],
                    named_histograms[h_reco],
                    named_histograms[h_unfold],
                    named_histograms[h_unfold],
                    simtype.c_str(), simtypeReco.c_str(), "Reco Unfolded", " ");

                TLegend* leg_sphenix = nullptr;
                std::string leg_sphenix_key = "leg_sphenix_" + canvas_key;
               // Map R code to actual radius string

		std::string R_value_str = (R == "R2") ? "anti-#it{k}_{#it{t}} #it{R} = 0.2" : "anti-#it{k}_{#it{t}} #it{R} = 0.4";
xj_functions::DrawsPHENIXLegendNew(leg_sphenix, leg_sphenix_key.c_str(), 0.13, 0.52, 0.43, 0.90, "p+p #sqrt{s}=200 GeV", R_value_str.c_str(), "20.9 < p_{T,1} < 60.8 GeV", "Matched Truth Dijets", true, true, false, true, false);
leg_sphenix->AddEntry("", (closure == "half" ? "Half Closure Test" : "Full Closure Test"), "");




                std::string outpath = "xj_plot_calofit_pdfs/" + canvas_key + ".pdf";
                canvases[canvas_key]->SaveAs(outpath.c_str());

                std::cout << "Saved: " << outpath << std::endl;
            }
        }
    }
}



//Base


for (const auto& pt_range : pt_ranges_variation) {
    for (const auto& gen : gens) {
        for (const auto& R : Rvals) {

            std::string tag = pt_range + "_" + gen + "_" + R;
            std::string canvas_key = "canvas_base_" + gen + "_" + R + "_" + pt_range;

            std::string h_reco_half     = "h_xj_reco_half_"     + tag;
            std::string h_unfold_data   = "h_xj_unfold_data_"   + tag;
            std::string h_unfold_half   = "h_xj_unfold_half_"   + tag;
            std::string h_data          = "h_xj_data_"          + tag;
            std::string h_reco_full     = "h_xj_reco_full_"     + tag;

            if (!canvases.count(canvas_key) || !canvases[canvas_key]) {
                std::cerr << "Missing canvas: " << canvas_key << std::endl;
                continue;
            }
            if (!named_histograms.count(h_reco_half)   || !named_histograms[h_reco_half] ||
                !named_histograms.count(h_unfold_data) || !named_histograms[h_unfold_data] ||
                !named_histograms.count(h_unfold_half) || !named_histograms[h_unfold_half] ||
                !named_histograms.count(h_data)        || !named_histograms[h_data] ||
                !named_histograms.count(h_reco_full)   || !named_histograms[h_reco_full]) {
                std::cerr << "Missing one or more histograms for: " << tag << std::endl;
                continue;
            }

            // Main Plot
            xj_functions::plotxj_and_ratio(canvases[canvas_key],
                named_histograms[h_reco_half],
                named_histograms[h_unfold_data],
                named_histograms[h_unfold_half],
                named_histograms[h_data],
                true, "x_{J, Data}", "Data/Reco",
                0, 1.0, 4, 20, 1, 53, 4, 53, 1, 20);

            // Main Legend
            TLegend* leg_main = nullptr;
            std::string legend_key = "leg_" + canvas_key;
            std::string simtype            = (gen == "calofit_pythia8") ? "PYTHIA 8" : "HERWIG 7";
            std::string simtypeReco        = simtype + " Reco";
            std::string simtypeRecoUnfold  = simtype + " Reco Unfold";

            xj_functions::draw_xj_legend(leg_main, legend_key.c_str(), 0.61, 0.75, 0.71, 0.90,
                named_histograms[h_reco_half],
                named_histograms[h_unfold_half],
                named_histograms[h_data],
                named_histograms[h_unfold_data],
                simtypeReco.c_str(), simtypeRecoUnfold.c_str(), "Data", "Data Unfolded", true);

            // sPHENIX Legend
            TLegend* leg_sphenix = nullptr;
            std::string leg_sphenix_key = "leg_sphenix_" + canvas_key;
            std::string R_value_str = (R == "R2") ? "anti-#it{k}_{#it{t}} #it{R} = 0.2" : "anti-#it{k}_{#it{t}} #it{R} = 0.4";

            xj_functions::DrawsPHENIXLegendNew(leg_sphenix, leg_sphenix_key.c_str(),
                0.13, 0.52, 0.43, 0.90,
                "p+p #sqrt{s}=200 GeV", R_value_str.c_str(), pt_labels[pt_range].c_str(), "Matched Truth Dijets", true, true, false, true, false);

            // Save Canvas
            std::string outpath = "xj_plot_calofit_pdfs/" + canvas_key + ".pdf";
            canvases[canvas_key]->SaveAs(outpath.c_str());
            std::cout << "Saved: " << outpath << std::endl;
        }
    }
}



//JER

 
for (const auto& pt_range : pt_ranges_variation) {
    for (const auto& gen : gens) {
        for (const auto& R : Rvals) {

            std::string base_tag = pt_range + "_" + gen + "_" + R;
            std::string jer_tag_pos = pt_range + "_" + gen + "_posJER_" + R;
            std::string jer_tag_neg = pt_range + "_" + gen + "_negJER_" + R;

            std::string canvas_jer_key  = "canvas_JER_" + gen + "_" + R + "_" + pt_range;
            std::string h_jer_pos       = "h_xj_unfold_data_" + jer_tag_pos;
            std::string h_jer_neg       = "h_xj_unfold_data_" + jer_tag_neg;
            std::string h_unfold        = "h_xj_unfold_data_"       + base_tag;
            std::string h_truth         = "h_xj_truth_full_"   + base_tag;

            if (!canvases.count(canvas_jer_key) || !canvases[canvas_jer_key]) {
                std::cerr << "Missing canvas: " << canvas_jer_key << std::endl;
                continue;
            }

            bool missing = false;
            if (!named_histograms.count(h_truth) || !named_histograms[h_truth]) {
                std::cerr << "Missing histogram: " << h_truth << std::endl;
                missing = true;
            }
            if (!named_histograms.count(h_unfold) || !named_histograms[h_unfold]) {
                std::cerr << "Missing histogram: " << h_unfold << std::endl;
                missing = true;
            }
            if (!named_histograms.count(h_jer_pos) || !named_histograms[h_jer_pos]) {
                std::cerr << "Missing histogram: " << h_jer_pos << std::endl;
                missing = true;
            }
            if (!named_histograms.count(h_jer_neg) || !named_histograms[h_jer_neg]) {
                std::cerr << "Missing histogram: " << h_jer_neg << std::endl;
                missing = true;
            }
            if (missing) continue;

            // Plot JER uncertainty
            xj_functions::plotxj_and_tworatios(canvases[canvas_jer_key],
                named_histograms[h_truth],      // h1
                named_histograms[h_unfold],     // horig
                named_histograms[h_jer_pos],    // hdiv (Unfold + JER)
                named_histograms[h_jer_neg],    // hextra (Unfold - JER)
                true, "x_{J}", "Unfold / (Unfold #pm JER)",
                0, 1.0, kGreen + 3, 20, 1, 20, kMagenta - 2, 20, 2, 20);

            std::string simtype = (gen == "calofit_pythia8") ? "PYTHIA 8" : "HERWIG 7";
            std::string legend_key = "leg_" + canvas_jer_key;

            TLegend* leg_main = nullptr;
            xj_functions::draw_xj_legend(leg_main, legend_key.c_str(), 0.61, 0.75, 0.71, 0.90,
                named_histograms[h_truth],
                named_histograms[h_unfold],
                named_histograms[h_jer_pos],
                named_histograms[h_jer_neg],
					 simtype.c_str(), "Unfold", "Unfold + 5% JER", "Unfold - 5% JER",true);

            TLegend* leg_sphenix = nullptr;
            std::string leg_sphenix_key = "leg_sphenix_" + canvas_jer_key;
            std::string R_value_str = (R == "R2") ? "anti-#it{k}_{#it{t}} #it{R} = 0.2" : "anti-#it{k}_{#it{t}} #it{R} = 0.4";
            std::string pt_range_str = (pt_range == "20_30") ? "20.9 < p_{T,1} < 31.2 GeV" :
                                       (pt_range == "30_40") ? "31.2 < p_{T,1} < 40.7 GeV" :
                                                               "40.7 < p_{T,1} < 60.8 GeV";

            xj_functions::DrawsPHENIXLegendNew(leg_sphenix, leg_sphenix_key.c_str(),
                0.13, 0.52, 0.43, 0.90, "p+p #sqrt{s}=200 GeV",
                R_value_str.c_str(), pt_range_str.c_str(), "Matched Truth Dijets",
                true, true, false, true, false);

            std::string outpath_jer = "xj_plot_calofit_pdfs/" + canvas_jer_key + ".pdf";
            canvases[canvas_jer_key]->SaveAs(outpath_jer.c_str());
            std::cout << "Saved: " << outpath_jer << std::endl;
        }
    }
}


//JES
for (const auto& pt_range : pt_ranges_variation) {
    for (const auto& gen : gens) {
        for (const auto& R : Rvals) {

            std::string base_tag    = pt_range + "_" + gen + "_" + R;
            std::string jes_tag_pos = pt_range + "_" + gen + "_posJES_" + R;
            std::string jes_tag_neg = pt_range + "_" + gen + "_negJES_" + R;

            std::string canvas_jes_key  = "canvas_JES_" + gen + "_" + R + "_" + pt_range;
            std::string h_jes_pos       = "h_xj_unfold_data_" + jes_tag_pos;
            std::string h_jes_neg       = "h_xj_unfold_data_" + jes_tag_neg;
            std::string h_unfold        = "h_xj_unfold_data_" + base_tag;
            std::string h_truth         = "h_xj_truth_full_"  + base_tag;

            if (!canvases.count(canvas_jes_key) || !canvases[canvas_jes_key]) {
                std::cerr << "Missing canvas: " << canvas_jes_key << std::endl;
                continue;
            }

            bool missing = false;
            if (!named_histograms.count(h_truth) || !named_histograms[h_truth]) {
                std::cerr << "Missing histogram: " << h_truth << std::endl;
                missing = true;
            }
            if (!named_histograms.count(h_unfold) || !named_histograms[h_unfold]) {
                std::cerr << "Missing histogram: " << h_unfold << std::endl;
                missing = true;
            }
            if (!named_histograms.count(h_jes_pos) || !named_histograms[h_jes_pos]) {
                std::cerr << "Missing histogram: " << h_jes_pos << std::endl;
                missing = true;
            }
            if (!named_histograms.count(h_jes_neg) || !named_histograms[h_jes_neg]) {
                std::cerr << "Missing histogram: " << h_jes_neg << std::endl;
                missing = true;
            }
            if (missing) continue;

            // Plot JES uncertainty
            xj_functions::plotxj_and_tworatios(canvases[canvas_jes_key],
                named_histograms[h_truth],      // h1
                named_histograms[h_unfold],     // horig
                named_histograms[h_jes_pos],    // hdiv (Unfold + JES)
                named_histograms[h_jes_neg],    // hextra (Unfold - JES)
                true, "x_{J}", "Unfold / (Unfold #pm JES)",
                0, 1.0, kGreen + 3, 20, 1, 20, kMagenta - 2, 20, 2, 20);

            std::string simtype = (gen == "calofit_pythia8") ? "PYTHIA 8" : "HERWIG 7";
            std::string legend_key = "leg_" + canvas_jes_key;

            TLegend* leg_main = nullptr;
            xj_functions::draw_xj_legend(leg_main, legend_key.c_str(), 0.61, 0.75, 0.71, 0.90,
                named_histograms[h_truth],
                named_histograms[h_unfold],
                named_histograms[h_jes_pos],
                named_histograms[h_jes_neg],
                simtype.c_str(), "Unfold", "Unfold + 5% JES", "Unfold - 5% JES", true);

            TLegend* leg_sphenix = nullptr;
            std::string leg_sphenix_key = "leg_sphenix_" + canvas_jes_key;
            std::string R_value_str = (R == "R2") ? "anti-#it{k}_{#it{t}} #it{R} = 0.2" : "anti-#it{k}_{#it{t}} #it{R} = 0.4";
            std::string pt_range_str = (pt_range == "20_30") ? "20.9 < p_{T,1} < 31.2 GeV" :
                                       (pt_range == "30_40") ? "31.2 < p_{T,1} < 40.7 GeV" :
                                                               "40.7 < p_{T,1} < 60.8 GeV";

            xj_functions::DrawsPHENIXLegendNew(leg_sphenix, leg_sphenix_key.c_str(),
                0.13, 0.52, 0.43, 0.90, "p+p #sqrt{s}=200 GeV",
                R_value_str.c_str(), pt_range_str.c_str(), "Matched Truth Dijets",
                true, true, false, true, false);

            std::string outpath_jes = "xj_plot_calofit_pdfs/" + canvas_jes_key + ".pdf";
            canvases[canvas_jes_key]->SaveAs(outpath_jes.c_str());
            std::cout << "Saved: " << outpath_jes << std::endl;
        }
    }
}

 





 
}//close macro
