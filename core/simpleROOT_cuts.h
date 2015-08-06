#ifndef simpleROOT_cuts_h
#define simpleROOT_cuts_h

TCut sel_ee("lepID[0]*lepID[1]==-11*11");
TCut sel_mumu("lepID[0]*lepID[1]==-13*13");
TCut sel_emu("lepID[0]*lepID[1]==-11*13");
TCut sel_basic("goodVtx && lepPt[0]>25 && lepPt[1]>20 && l1l2DR>0.3 && !(abs(lepEta[0]) > 1.4 && abs(lepEta[0]) < 1.6) && !(abs(lepEta[1]) > 1.4 && abs(lepEta[1]) < 1.6)");
TCut sel_CE("abs(lepEta[0])<1.4 && abs(lepEta[1])<1.4");
TCut sel_FW = sel_basic && !sel_CE;
TCut sel_SF = sel_ee || sel_mumu;
TCut sel_OF = sel_emu;
TCut sel_M2070("l1l2M>20 && l1l2M<70");
TCut sel_basic_sync("lepPt[0]>25 && lepPt[1]>20 && l1l2DR>0.3 && !(abs(lepEta[0]) > 1.4 && abs(lepEta[0]) < 1.6) && !(abs(lepEta[1]) > 1.4 && abs(lepEta[1]) < 1.6)");
TCut sel_lep_iso("lepIso[0]<0.15 && lepIso[1]<0.15");
TCut sel_M70110("l1l2M>70 && l1l2M<110");
TCut sel_M81101("l1l2M>81 && l1l2M<101");
TCut sel_M0150("l1l2M<150");
TCut sel_M50("l1l2M>50");
TCut sel_jetsMET8TeV("(njets == 2 && t1met > 150) || (njets >=3 && t1met > 100)");
TCut sel_ej0("njets==0");
TCut sel_ej1("njets==1");
TCut sel_ej2("njets==2");
TCut sel_ej3("njets==3");
TCut sel_ej4("njets==4");
TCut sel_ij0("njets>=0");
TCut sel_ij1("njets>=1");
TCut sel_ij2("njets>=2");
TCut sel_ij3("njets>=3");
TCut sel_ij4("njets>=4");


TCut sel_FSSub("((lepID[0]*lepID[1] == -11*13) ? -1 : 0) + ((lepID[0]*lepID[1] == -11*11) ? 1 : 0) + ((lepID[0]*lepID[1] == -13*13) ? 1 : 0)");

TCut sel_ee_trig_fire("HLT_e1e2 == 1");
TCut sel_mumu_trig_fire("HLT_mu1mu2 == 1");
TCut sel_emu_trig_fire("HLT_mu1e2 == 1 || HLT_e1mu2 == 1");
TCut sel_prescale_check("!(HLT_e1e2 > 1 || HLT_mu1mu2 > 1 || HLT_mu1e2 > 1 || HLT_e1mu2 > 1)"); // should be false if at least one of the diepton paths is zero

TCut sel_ee_trig = sel_ee && sel_ee_trig_fire && sel_prescale_check;
TCut sel_mumu_trig = sel_mumu && sel_mumu_trig_fire && sel_prescale_check;
TCut sel_emu_trig = sel_emu && sel_emu_trig_fire && sel_prescale_check;
TCut sel_SF_trig = sel_ee_trig || sel_mumu_trig;
TCut sel_OF_trig = sel_emu_trig;

TCut sel_dileptontrigger_nleps2 = (sel_ee_trig || sel_mumu_trig || sel_emu_trig);
TCut sel_dileptontrigger_nleps2_ij1 = (sel_ee_trig || sel_mumu_trig || sel_emu_trig) && sel_ij1;
TCut sel_dileptontrigger_nleps2_ij2 = (sel_ee_trig || sel_mumu_trig || sel_emu_trig) && sel_ij2;
TCut sel_dileptontrigger_nleps2_ij3 = (sel_ee_trig || sel_mumu_trig || sel_emu_trig) && sel_ij3;
// events->SetAlias("jetHT","Sum$(jetPt)")


TCut sel_gen_ee("genlepID[0]*genlepID[1]==-11*11");
TCut sel_gen_mumu("genlepID[0]*genlepID[1]==-13*13");
TCut sel_gen_emu("genlepID[0]*genlepID[1]==-11*13");
TCut sel_gen_basic("genlepPt[0]>25 && genlepPt[1]>20 && genl1l2DR>0.3");
TCut sel_gen_CE("abs(genlepEta[0])<1.4 && abs(genlepEta[1])<1.4");
TCut sel_gen_SF = sel_gen_ee || sel_gen_mumu;
TCut sel_gen_OF = sel_gen_emu;
TCut sel_gen_M2070("genl1l2M>20 && genl1l2M<70");

TCut sel_ZZ_correct("(genlepID[0]*genlepID[1]==-11*11 && abs(genlepID[2])==13) || (genlepID[0]*genlepID[1]==-13*13 && abs(genlepID[2])==11) ");

TCut sel_recoMatchGen("lepGenMatchIndex[0]!=-1 && lepGenMatchIndex[1]!=-1"); //both leptons have been matched
TCut sel_hierarchyMatch("genlepID[lepGenMatchIndex[0]] == lepID[0] && genlepID[lepGenMatchIndex[1]] == lepID[1]"); //both leptons have been matched with correct hierarchy


// --- signle leptons cuts

TCut sel_mu1("abs(lepID[0]) == 13");
TCut sel_gen_mu1("abs(lepID[0]) == 13");
TCut sel_recoMatchGen_mu1("lepGenMatchIndex[0]!=-1 && genlepID[lepGenMatchIndex[0]] == lepID[0]");

TCut sel_mu2("abs(lepID[1]) == 13");
TCut sel_gen_mu2("abs(lepID[1]) == 13");
TCut sel_recoMatchGen_mu2("lepGenMatchIndex[1]!=-1 && genlepID[lepGenMatchIndex[1]] == lepID[1]");

TCut sel_e1("abs(lepID[0]) == 11");
TCut sel_gen_e1("abs(lepID[0]) == 11");
TCut sel_recoMatchGen_e1("lepGenMatchIndex[0]!=-1 && genlepID[lepGenMatchIndex[0]] == lepID[0]");

TCut sel_e2("abs(lepID[1]) == 11");
TCut sel_gen_e2("abs(lepID[1]) == 11");
TCut sel_recoMatchGen_e2("lepGenMatchIndex[1]!=-1 && genlepID[lepGenMatchIndex[1]] == lepID[1]");
#endif
