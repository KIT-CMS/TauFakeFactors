templates:
  split_schemes:
    3j: &3j_split
      split_categories:
        njets: ["==0", "==1", ">=2"]
      split_categories_binedges:
        njets: [-0.5, 0.5, 1.5, 22.5]
      correction_option:
        "==0": "smoothed"
        "==1": "smoothed"
        ">=2": "smoothed"
      bandwidth:
        "==0": 1.0
        "==1": 1.0
        ">=2": 1.0
    2j: &2j_split
      split_categories:
        njets: ["<=1", ">=2"]
      split_categories_binedges:
        njets: [-0.5, 1.5, 22.5]
      correction_option:
        "<=1": "smoothed"
        ">=2": "smoothed"
      bandwidth:
        "<=1": 1
        ">=2": 1
  var_dependence_n_bins:
    QCD:
      equipopulated_binning_options: &QCD_var_dependence_n_bins
        var_dependence_n_bins:
          "==0": 11
          "==1": 9
          ">=2": 7
    Wjets:
      equipopulated_binning_options: &Wjets_var_dependence_n_bins
        var_dependence_n_bins:
          "==0": 15
          "==1": 11
          ">=2": 9
    ttbar:
      equipopulated_binning_options: &ttbar_var_dependence_n_bins
        var_dependence_n_bins:
          "<=1": 11
          ">=2": 15
  correction_variations__with_mc_subtraction_shift: &correction_variations__with_mc_subtraction_shift
    correction_variations:
      - StatShift
      - SystMCShift
      - SystBandAsym
  correction_variations__without_mc_subtraction_shift: &correction_variations__without_mc_subtraction_shift
    correction_variations:
      - StatShift
      - SystBandAsym
  variables:
    eta_1:
      var_dependence: eta_1
      equipopulated_binning_options:
        variable_config:
          eta_1:
            min: -2.1
            max: +2.1
            rounding: 2
      correction_option: "smoothed"
      bandwidth: 1.0
    eta_2:
      var_dependence: eta_2
      equipopulated_binning_options:
        variable_config:
          eta_2:
            min: -2.5
            max: +2.5
            rounding: 2
      correction_option: "smoothed"
      bandwidth: 1.5
    deltaEta_ditaupair: &deltaEta_ditaupair
      var_dependence: deltaEta_ditaupair
      equipopulated_binning_options: &deltaEta_ditaupair__eq_bin_opt
        variable_config:
          deltaEta_ditaupair:
            min: -4.9
            max: +4.9
            rounding: 2
      correction_option: "smoothed"
      bandwidth: 1.75
    deltaR_ditaupair: &deltaR_ditaupair
      var_dependence: deltaR_ditaupair
      equipopulated_binning_options: &deltaR_ditaupair__eq_bin_opt
        variable_config:
          deltaR_ditaupair:
            min: 0.5
            max: 5.0
            rounding: 4
      correction_option: "smoothed"
    jeta_1:
      var_dependence: jeta_1
      equipopulated_binning_options:
        variable_config:
          jeta_1:
            min: -5.0
            max: +5.0
            rounding: 2
      correction_option: "smoothed"
      bandwidth: 2.0
    jeta_2:
      var_dependence: jeta_2
      equipopulated_binning_options:
        variable_config:
          jeta_2:
            min: -5.0
            max: +5.0
            rounding: 2
      correction_option: "smoothed"
      bandwidth: 2.0
    jpt_1:
      var_dependence: jpt_1
      equipopulated_binning_options:
        variable_config:
          jpt_1:
            min: 30.0
            max: 150.0
            rounding: 2
      correction_option: "smoothed"
    jpt_2:
      var_dependence: jpt_2
      equipopulated_binning_options:
        variable_config:
          jpt_2:
            min: 30.0
            max: 150.0
            rounding: 2
      correction_option: "smoothed"
    met:
      var_dependence: met
      equipopulated_binning_options:
        variable_config:
          met:
            min: 0.0
            max: 150.0
            rounding: 2
      correction_option: "smoothed"
      bandwidth: 25
    mt_1:
      var_dependence: mt_1
      equipopulated_binning_options:
        variable_config:
          mt_1:
            min: 0.0
            max: 70.0
            rounding: 2
      correction_option: "smoothed"
      bandwidth: 20
    pt_tt:
      var_dependence: pt_tt
      equipopulated_binning_options:
        variable_config:
          pt_tt:
            min: 0.0
            max: 150.0
            rounding: 2
      correction_option: "smoothed"
      bandwidth: 25
    pt_ttjj: &pt_ttjj
      var_dependence: pt_ttjj
      equipopulated_binning_options: &pt_ttjj__eq_bin_opt
        variable_config:
          pt_ttjj:
            min: 0.0
            max: 150.0
            rounding: 2
      correction_option: "smoothed"
      bandwidth: 25
    mt_tot:
      var_dependence: mt_tot
      equipopulated_binning_options:
        variable_config:
          mt_tot:
            min: 0.0
            max: 250.0
            rounding: 2
      correction_option: "smoothed"
      bandwidth: 40
    mass_2: &mass_2
      var_dependence: mass_2
      equipopulated_binning_options: &mass_2__eq_bin_opt
        variable_config:
          mass_2:
            min: 0.2
            max: 2.0
            rounding: 4
        add_left: [0.0]
      correction_option: "binwise#[0]+smoothed"
      bandwidth: 0.3
    tau_decaymode_2: &tau_decaymode_2
      var_dependence: tau_decaymode_2
      correction_option: "binwise"
      bandwidth: 1.0
      var_bins: [-0.5, 0.5, 9.5, 10.5, 11.5]
    nbtag:
      var_dependence: nbtag
      correction_option: "binwise"
      bandwidth: 1.0
      var_bins: [-0.5, 0.5, 1.5, 2.5, 22.5]
    iso_1: &iso_1
      var_dependence: iso_1
      equipopulated_binning_options: &iso_1__eq_bin_opt
        variable_config:
          iso_1:
            min: 0.00005
            max: 0.15
            rounding: 6
        add_left: [0.0]
      correction_option: "binwise#[0]+smoothed"
      bandwidth: 0.02
    deltaR_1j1: &deltaR_1j1
      var_dependence: deltaR_1j1
      equipopulated_binning_options: &deltaR_1j1__eq_bin_opt
        variable_config:
          deltaR_1j1:
            min: 0.5
            max: 7.0
            rounding: 3
        correction_option: "smoothed"
        bandwidth: 0.9
    deltaR_12j1:
      var_dependence: deltaR_12j1
      equipopulated_binning_options:
        variable_config:
          deltaR_12j1:
            min: 0.0
            max: 10.0
        correction_option: "smoothed"
        bandwidth: 1.25
    m_vis: &m_vis
      var_dependence: m_vis
      equipopulated_binning_options: &m_vis__eq_bin_opt
        variable_config:
          m_vis:
            min: 40.0
            max: 250.0
            rounding: 2
      correction_option: "smoothed"
      bandwidth: 30

channel: mt
target_processes:
  QCD:
    chain_DR_SR_to_non_closure: false

    non_closure:
      tau_decaymode_2:
        <<: [*tau_decaymode_2, *3j_split]
        correction_variations: ["StatShift", "SystMCShift"]
      deltaEta_ditaupair:
        <<: [*deltaEta_ditaupair, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*deltaEta_ditaupair__eq_bin_opt, *QCD_var_dependence_n_bins]
        bandwidth: [1.25, 1.5, 2.0]
        var_bins:
          "==0": [-4.9, -2.13, -1.42, -0.93, -0.53, -0.18, 0.16, 0.56, 0.96, 1.43, 2.08, 4.9]
          "==1": [-4.9, -2.0, -1.23, -0.65, -0.23, 0.25, 0.72, 1.28, 2.04, 4.9]
          ">=2": [-4.9, -1.9, -1.03, -0.39, 0.26, 0.87, 1.76, 4.9]
      deltaR_ditaupair:
        <<: [*deltaR_ditaupair, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*deltaR_ditaupair__eq_bin_opt, *QCD_var_dependence_n_bins]
        bandwidth: [0.8, 1.1, 1.2]
        var_bins:
          "==0": [0.5, 2.6814, 2.869, 2.9787, 3.0526, 3.1104, 3.1696, 3.2609, 3.3865, 3.5747, 3.9186, 5.0]
          "==1": [0.5, 1.653, 2.227, 2.5243, 2.7304, 2.9125, 3.0715, 3.2696, 3.6723, 5.0]
          ">=2": [0.5, 1.3291, 1.9185, 2.3753, 2.7661, 3.0644, 3.4088, 5.0]
      deltaR_1j1:
        <<: [*deltaR_1j1, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*deltaR_1j1__eq_bin_opt, *QCD_var_dependence_n_bins]
        correction_option: ["skip", "smoothed", "smoothed"]
        var_bins:
          "==0": [0.5, 7.0]
          "==1": [0.5, 1.67, 2.205, 2.51, 2.774, 3.019, 3.235, 3.586, 4.235, 7.0]
          ">=2": [0.5, 1.825, 2.382, 2.761, 3.009, 3.272, 3.759, 7.0]
      pt_ttjj:
        <<: [*pt_ttjj, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*pt_ttjj__eq_bin_opt, *QCD_var_dependence_n_bins]
        correction_option: ["skip", "skip", "smoothed"]
        bandwidth: [1.0, 1.0, 35.0]
        var_bins:
          "==0": [0.0, 150.0]
          "==1": [0.0, 150.0]
          ">=2": [0.0, 15.52, 22.28, 30.06, 38.0, 49.31, 69.77, 150.0]
      mass_2:
        <<: [*mass_2, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*mass_2__eq_bin_opt, *QCD_var_dependence_n_bins]
        bandwidth: [0.25, 0.35, 0.4]
        var_bins:
          "==0": [0.0, 0.2, 0.5107, 0.6216, 0.7388, 0.835, 0.9121, 0.9941, 1.0732, 1.1543, 1.25, 1.3623, 2.0]
          "==1": [0.0, 0.2, 0.5322, 0.665, 0.79, 0.896, 1.0088, 1.0957, 1.1982, 1.3271, 2.0]
          ">=2": [0.0, 0.2, 0.564, 0.7314, 0.8735, 0.9932, 1.1123, 1.2715, 2.0]
      iso_1:
        <<: [*iso_1, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: *QCD_var_dependence_n_bins
          variable_config:
            iso_1:
              min: 0.05
              max: 0.15
              rounding: 6
        bandwidth: [0.04, 0.04, 0.04]
        correction_option: "smoothed"
        var_bins:
          "==0": [0.05, 0.056882, 0.063543, 0.070468, 0.078181, 0.086647, 0.095377, 0.105744, 0.115638, 0.126327, 0.137486, 0.15]
          "==1": [0.05, 0.057096, 0.065105, 0.073987, 0.083517, 0.096737, 0.108642, 0.122276, 0.13489, 0.15]
          ">=2": [0.05, 0.059351, 0.068657, 0.081325, 0.096132, 0.110882, 0.128247, 0.15]

  Wjets:
    chain_DR_SR_to_non_closure: false

    non_closure:
      tau_decaymode_2:
        <<: [*tau_decaymode_2, *3j_split]
        correction_variations: ["StatShift", "SystMCShift"]
      deltaEta_ditaupair:
        <<: [*deltaEta_ditaupair, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*deltaEta_ditaupair__eq_bin_opt, *Wjets_var_dependence_n_bins]
        bandwidth: [0.8, 1.1, 1.25]
        var_bins:
          "==0": [-4.9, -1.9, -1.36, -1.02, -0.75, -0.52, -0.31, -0.11, 0.09, 0.3, 0.51, 0.74, 1.0, 1.34, 1.87, 4.9]
          "==1": [-4.9, -1.75, -1.16, -0.75, -0.43, -0.14, 0.14, 0.43, 0.74, 1.14, 1.74, 4.9]
          ">=2": [-4.9, -1.56, -0.96, -0.53, -0.16, 0.19, 0.54, 0.95, 1.54, 4.9]
      deltaR_ditaupair:
        <<: [*deltaR_ditaupair, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*deltaR_ditaupair__eq_bin_opt, *Wjets_var_dependence_n_bins]
        bandwidth: [0.5, 0.6, 0.7]
        var_bins:
          "==0": [0.5, 1.291, 1.6664, 1.9482, 2.1706, 2.3528, 2.5196, 2.6606, 2.783, 2.8902, 2.9821, 3.0643, 3.1321, 3.2287, 3.4481, 5.0]
          "==1": [0.5, 1.1946, 1.5966, 1.9154, 2.1746, 2.4058, 2.6139, 2.801, 2.9638, 3.1049, 3.3191, 5.0]
          ">=2": [0.5, 1.1567, 1.5816, 1.9338, 2.2206, 2.4933, 2.7366, 2.9631, 3.1951, 5.0]
      deltaR_1j1:
        <<: [*deltaR_1j1, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*deltaR_1j1__eq_bin_opt, *Wjets_var_dependence_n_bins]
        correction_option: ["skip", "smoothed", "smoothed"]
        bandwidth: [1.0, 1.0, 1.2]
        var_bins:
          "==0": [0.5, 7.0]
          "==1": [0.5, 1.257, 1.688, 2.026, 2.312, 2.572, 2.811, 3.037, 3.282, 3.69, 4.342, 7.0]
          ">=2": [0.5, 1.242, 1.733, 2.12, 2.446, 2.72, 2.978, 3.261, 3.842, 7.0]
      pt_ttjj:
        <<: [*pt_ttjj, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*pt_ttjj__eq_bin_opt, *Wjets_var_dependence_n_bins]
        correction_option: ["skip", "skip", "smoothed"]
        bandwidth: [1.0, 1.0, 30]
        var_bins:
          "==0": [0.0, 150.0]
          "==1": [0.0, 150.0]
          ">=2": [0.0, 13.13, 19.85, 25.75, 31.71, 38.56, 46.7, 57.79, 75.16, 150.0]
      mass_2:
        <<: [*mass_2, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*mass_2__eq_bin_opt, *Wjets_var_dependence_n_bins]
        bandwidth: [0.225, 0.225, 0.25]
        var_bins:
          "==0": [0.0, 0.2, 0.4827, 0.5635, 0.6396, 0.707, 0.7681, 0.8281, 0.8857, 0.9453, 1.0078, 1.0713, 1.1377, 1.2109, 1.2949, 1.4092, 2.0]
          "==1": [0.0, 0.2, 0.5088, 0.6177, 0.7119, 0.7959, 0.8696, 0.9551, 1.0381, 1.1279, 1.2285, 1.3545, 2.0]
          ">=2": [0.0, 0.2, 0.5308, 0.6577, 0.7612, 0.8506, 0.957, 1.0625, 1.1758, 1.3076, 2.0]
      iso_1:
        <<: [*iso_1, *3j_split, *correction_variations__with_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*iso_1__eq_bin_opt, *Wjets_var_dependence_n_bins]
        bandwidth: [0.03, 0.04, 0.05]
        var_bins:
          "==0": [0.0, 0.00005, 0.004932, 0.007516, 0.010087, 0.012785, 0.015814, 0.019493, 0.023479, 0.028324, 0.034053, 0.040962, 0.049714, 0.060662, 0.075963, 0.100649, 0.15]
          "==1": [0.0, 0.00005, 0.005621, 0.009072, 0.012884, 0.017279, 0.022713, 0.029417, 0.037605, 0.049066, 0.064637, 0.090462, 0.15]
          ">=2": [0.0, 0.00005, 0.005602, 0.009668, 0.014723, 0.020748, 0.028099, 0.039254, 0.055474, 0.082078, 0.15]

  ttbar:
    non_closure:
      tau_decaymode_2:
        <<: [*tau_decaymode_2, *2j_split]
        correction_variations: ["StatShift"]
      deltaEta_ditaupair:
        <<: [*deltaEta_ditaupair, *2j_split, *correction_variations__without_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*deltaEta_ditaupair__eq_bin_opt, *ttbar_var_dependence_n_bins]
        bandwidth: [1.5, 1.75]
        var_bins:
          "<=1": [-4.9, -1.57, -1.06, -0.71, -0.43, -0.17, 0.11, 0.39, 0.69, 1.05, 1.56, 4.9]
          ">=2": [-4.9, -2.0, -1.47, -1.1, -0.82, -0.58, -0.35, -0.12, 0.1, 0.33, 0.55, 0.79, 1.07, 1.43, 1.99, 4.9]
      deltaR_ditaupair:
        <<: [*deltaR_ditaupair, *2j_split, *correction_variations__without_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*deltaR_ditaupair__eq_bin_opt, *ttbar_var_dependence_n_bins]
        bandwidth: [0.75, 1.0]
        var_bins:
          "<=1": [0.5, 1.2589, 1.6931, 2.0017, 2.272, 2.4783, 2.6662, 2.8297, 2.9717, 3.103, 3.283, 5.0]
          ">=2": [0.5, 0.9472, 1.2341, 1.4726, 1.6983, 1.9055, 2.0972, 2.2769, 2.4468, 2.6064, 2.7588, 2.9015, 3.0423, 3.1826, 3.465, 5.0]
      deltaR_1j1:
        <<: [*deltaR_1j1, *2j_split, *correction_variations__without_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*deltaR_1j1__eq_bin_opt, *ttbar_var_dependence_n_bins]
        var_bins:
          "<=1": [0.5, 1.343, 1.785, 2.09, 2.344, 2.552, 2.723, 2.882, 3.037, 3.185, 3.478, 7.0]
          ">=2": [0.5, 1.158, 1.535, 1.836, 2.086, 2.291, 2.467, 2.618, 2.754, 2.872, 2.984, 3.087, 3.201, 3.395, 3.767, 7.0]
      pt_ttjj:
        <<: [*pt_ttjj, *2j_split, *correction_variations__without_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*pt_ttjj__eq_bin_opt, *ttbar_var_dependence_n_bins]
        correction_option: ["skip", "smoothed"]
        var_bins:
          "<=1": [0.0, 150.0]
          ">=2": [0.0, 15.23, 22.09, 27.98, 33.27, 38.42, 43.61, 48.82, 54.23, 60.19, 66.83, 74.34, 83.25, 95.43, 112.85, 150.0]
      mass_2:
        <<: [*mass_2, *2j_split, *correction_variations__without_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*mass_2__eq_bin_opt, *ttbar_var_dependence_n_bins]
        bandwidth: [0.225, 0.175]
        var_bins:
          "<=1": [0.0, 0.2, 0.5542, 0.7139, 0.8491, 0.9507, 1.0391, 1.1172, 1.1963, 1.2773, 1.374, 1.4785, 2.0]
          ">=2": [0.0, 0.2, 0.4985, 0.6196, 0.7344, 0.8311, 0.9072, 0.9736, 1.0352, 1.0957, 1.1553, 1.2168, 1.2793, 1.3447, 1.416, 1.5039, 2.0]
      iso_1:
        <<: [*iso_1, *2j_split, *correction_variations__without_mc_subtraction_shift]
        equipopulated_binning_options:
          <<: [*iso_1__eq_bin_opt, *ttbar_var_dependence_n_bins]
        bandwidth: [0.05, 0.05]
        var_bins:
          "<=1": [0.0, 0.00005, 0.006133, 0.009499, 0.01269, 0.017176, 0.022144, 0.027969, 0.0369, 0.048008, 0.063692, 0.090369, 0.15]
          ">=2": [0.0, 0.00005, 0.004676, 0.007108, 0.009587, 0.012164, 0.015111, 0.018535, 0.022638, 0.027259, 0.032977, 0.039821, 0.048317, 0.060037, 0.075678, 0.101079, 0.15]
