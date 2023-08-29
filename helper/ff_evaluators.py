import ROOT
import correctionlib


class FakeFactorEvaluator:
    def __init__(self, config, process, for_DRtoSR=False):
        if not for_DRtoSR:
            self.for_DRtoSR = ""
            self.ff_path = (
                "workdir/"
                + config["workdir_name"]
                + "/"
                + config["era"]
                + "/fake_factors_{}.json".format(config["channel"])
            )
        else:
            self.for_DRtoSR = "for_DRtoSR"
            self.ff_path = (
                "workdir/"
                + config["workdir_name"]
                + "/"
                + config["era"]
                + "/corrections/{ch}/fake_factors_{ch}_for_corrections.json".format(
                    ch=config["channel"]
                )
            )

        self.process = process
        correctionlib.register_pyroot_binding()
        print(f"Loading fake factor file {self.ff_path} for process {self.process}")
        if self.process not in ["QCD", "QCD_subleading", "Wjets", "ttbar"]:
            raise ValueError(f"Unknown process {self.process}")
        ROOT.gInterpreter.Declare(
            f'auto {self.process}_{self.for_DRtoSR} = correction::CorrectionSet::from_file("{self.ff_path}")->at("{self.process}_fake_factors");'
        )

    def evaluate_subleading_lep_pt_njets(self, rdf):
        rdf = rdf.Define(
            f"{self.process}_fake_factor",
            f'{self.process}_{self.for_DRtoSR}->evaluate({{pt_2, (float)njets, "nominal"}})',
        )
        return rdf
    
    def evaluate_subleading_boosted_lep_pt_njets(self, rdf):
        rdf = rdf.Define(
            f"{self.process}_fake_factor",
            f'{self.process}_{self.for_DRtoSR}->evaluate({{boosted_pt_2, (float)njets, "nominal"}})',
        )
        return rdf

    def evaluate_leading_lep_pt_njets(self, rdf):
        rdf = rdf.Define(
            f"{self.process}_fake_factor",
            f'{self.process}_{self.for_DRtoSR}->evaluate({{pt_1, (float)njets, "nominal"}})',
        )
        return rdf
    
    def evaluate_leading_boosted_lep_pt_njets(self, rdf):
        rdf = rdf.Define(
            f"{self.process}_fake_factor",
            f'{self.process}_{self.for_DRtoSR}->evaluate({{boosted_pt_1, (float)njets, "nominal"}})',
        )
        return rdf

    def evaluate_subleading_lep_pt_njets_deltaR(self, rdf):
        rdf = rdf.Define(
            f"{self.process}_fake_factor",
            f'{self.process}_{self.for_DRtoSR}->evaluate({{pt_2, (float)njets, deltaR_ditaupair, "nominal"}})',
        )
        return rdf
    
    def evaluate_subleading_boosted_lep_pt_njets_deltaR(self, rdf):
        rdf = rdf.Define(
            f"{self.process}_fake_factor",
            f'{self.process}_{self.for_DRtoSR}->evaluate({{pt_2, (float)njets, deltaR_ditaupair, "nominal"}})',
        )
        return rdf


class FakeFactorCorrectionEvaluator:
    def __init__(self, config, process, for_DRtoSR=True):
        if not for_DRtoSR:
            self.for_DRtoSR = ""
            self.corr_path = (
                "workdir/"
                + config["workdir_name"]
                + "/"
                + config["era"]
                + "/FF_corrections_{}.json".format(config["channel"])
            )
        else:
            self.for_DRtoSR = "for_DRtoSR"
            self.corr_path = (
                "workdir/"
                + config["workdir_name"]
                + "/"
                + config["era"]
                + "/corrections/{ch}/FF_corrections_{ch}_for_DRtoSR.json".format(
                    ch=config["channel"]
                )
            )

        self.process = process
        correctionlib.register_pyroot_binding()
        print(
            f"Loading fake factor correction file {self.corr_path} for process {self.process}"
        )
        if process not in ["QCD", "QCD_subleading", "Wjets", "ttbar"]:
            raise ValueError(f"Unknown process {self.process}")
        if process == "QCD" and config["channel"] == "tt":
            ROOT.gInterpreter.Declare(
                f'auto {self.process}_corr_{self.for_DRtoSR} = correction::CorrectionSet::from_file("{self.corr_path}")->at("{self.process}_non_closure_subleading_lep_pt_correction");'
            )
        else:
            ROOT.gInterpreter.Declare(
                f'auto {self.process}_corr_{self.for_DRtoSR} = correction::CorrectionSet::from_file("{self.corr_path}")->at("{self.process}_non_closure_leading_lep_pt_correction");'
            )

    def evaluate_leading_lep_pt(self, rdf):
        rdf = rdf.Define(
            f"{self.process}_ff_corr",
            f'{self.process}_corr_{self.for_DRtoSR}->evaluate({{pt_1, "nominal"}})',
        )
        return rdf
    
    def evaluate_leading_boosted_lep_pt(self, rdf):
        rdf = rdf.Define(
            f"{self.process}_ff_corr",
            f'{self.process}_corr_{self.for_DRtoSR}->evaluate({{boosted_pt_1, "nominal"}})',
        )
        return rdf

    def evaluate_subleading_lep_pt(self, rdf):
        rdf = rdf.Define(
            f"{self.process}_ff_corr",
            f'{self.process}_corr_{self.for_DRtoSR}->evaluate({{pt_2, "nominal"}})',
        )
        return rdf
    
    def evaluate_subleading_boosted_lep_pt(self, rdf):
        rdf = rdf.Define(
            f"{self.process}_ff_corr",
            f'{self.process}_corr_{self.for_DRtoSR}->evaluate({{boosted_pt_2, "nominal"}})',
        )
        return rdf