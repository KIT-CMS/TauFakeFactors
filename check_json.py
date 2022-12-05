import os
import correctionlib

work_folder = "new_tests_wo_emb_v6"
ch = "mt"
sf_dir = "workdir/{}/{}/".format(work_folder, ch)
json = correctionlib.CorrectionSet.from_file(os.path.join(sf_dir, "fake_factors.json"))


tau_pt = 56.
njets = 5.
nbtags = 1.
lep_mt = 23.
dR = 1.4

QCD_ff = json["QCD_fake_factors"].evaluate(tau_pt, njets)
print("QCD FF part:", QCD_ff)

Wjets_ff = json["Wjets_fake_factors"].evaluate(tau_pt, njets, dR)
print("Wjets FF part:", Wjets_ff)

ttbar_ff = json["ttbar_fake_factors"].evaluate(tau_pt, njets)
print("ttbar FF part:", ttbar_ff)

QCD_frac = json["process_fractions"].evaluate("QCD", lep_mt, nbtags)
print("QCD fraction:", QCD_frac)

Wjets_frac = json["process_fractions"].evaluate("Wjets", lep_mt, nbtags)
print("Wjets fraction:", Wjets_frac)

ttbar_frac = json["process_fractions"].evaluate("ttbar", lep_mt, nbtags)
print("ttbar fraction:", ttbar_frac)

fake_factor = QCD_ff*QCD_frac + Wjets_ff*Wjets_frac + ttbar_ff*ttbar_frac
print("full FF:", fake_factor)
