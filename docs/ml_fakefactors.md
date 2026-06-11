# ML-Based Fake Factor Evaluation

The framework supports the evaluation of fake factors based on machine learning models in the ONNX format instead of the classic approach of measuring them (see [Fake Factor Calculation](fakefactors.md)). ONNXRuntime is utilized to to evaluate the models within ROOT RDataFrames.

## ML Training

Not include yet. Files/models have to be provided externally.

## Configuration

To use ML-based fake factors, the configuration must specify the ONNX model details for each process. These details are defined in a YAML configuration file, which should be located in the `workdir/TAG` directory. The configuration includes:

- **Model Path**: The path to the ONNX model file.
- **Model Inputs**: A list of input variable names required by the model.
- **Define Columns**: Optional expressions to define missing columns dynamically.

Example YAML configuration:
```yaml
target_processes:
  QCD:
    model_path: "models/qcd_fake_factors.onnx"
    model_input:
      - "event_parity"
      - "pt_1"
      - "njets"
      ...
    define_columns:
      event_parity: event % 2
      ...
  Wjets:
    model_path: "models/wjets_fake_factors.onnx"
    model_input:
      - "pt_2"
      - "nbtag"
      ...
  ...
```

## Notes

- Ensure that the ONNX model is compatible with ONNXRuntime and that the input variable names match those in the configuration.