import argparse
import os
from collections import OrderedDict
from contextlib import suppress
from typing import Any

import yaml
import yaml.constructor

args = argparse.ArgumentParser(description='Config Adjustment')
args.add_argument('--directory', type=str, default='configs/smhtt_ul', help='config directroy path')
args.add_argument("--eras", type=str, nargs='+', default=['2016', '2017', '2018'], help='eras')
args.add_argument("--channels", type=str, nargs='+', default=['et', 'mt', 'tt'])
args.add_argument("--ntuple-path", type=str, help='ntuple path used for preselection')
args.add_argument("--output-file-path", type=str, help='output path after preselection')
args.add_argument("--workdir-name", type=str, help='Storage directory of calculation results')


def yaml_setup() -> None:
    def _keep_dict_order(self, data: Any) -> Any:
        return self.represent_mapping('tag:yaml.org,2002:omap', data.items())
    yaml.add_representer(OrderedDict, _keep_dict_order)    


if __name__ == '__main__':
    args = args.parse_args()
    yaml_setup()
    for era in args.eras:
        for file in os.listdir(f'{args.directory}/{era}'):
            if any(it in file for it in args.channels):
                with open(f'{args.directory}/{era}/{file}', 'r') as f:
                    config = OrderedDict(yaml.load(f, Loader=yaml.FullLoader))
                    initial_order = list(config.keys())
                    
                for key, replacement in [
                    ("ntuple_path", args.ntuple_path),
                    ("output_path", args.output_file_path),
                    ("file_paths", args.output_file_path),
                    ("workdir_name", args.workdir_name),
                ]:

                    with suppress(KeyError):
                        config[key] = replacement

                with open(f'{args.directory}/{era}/{file}', 'w') as yaml_file:
                    yaml.dump(
                        OrderedDict({k: config[k] for k in initial_order}),
                        yaml_file,
                        sort_keys=False,
                        indent=4,
                        default_flow_style=False,
                        width=float("inf"),
                    )

                with open(f'{args.directory}/{era}/{file}', 'r') as yaml_file:
                    data = yaml_file.read().splitlines(True)
                with open(f'{args.directory}/{era}/{file}', 'w') as yaml_file:
                    yaml_file.writelines(data[1:])
