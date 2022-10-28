
def QCD_SS_estimate(hists):
    qcd = hists["data"].Clone()

    for sample in hists:
        if sample not in ["data", "data_substracted", "QCD"] and "_T" not in sample:
            qcd.Add(hists[sample], -1)
    # check for negative bins
    for i in range(qcd.GetNbinsX()):
        if qcd.GetBinContent(i) < 0.:
            qcd.SetBinContent(i, 0.)
    return qcd