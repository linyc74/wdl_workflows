import json
import pandas as pd
from typing import Dict, Any


class Main:

    TEMPLATE_JSON = f'{__file__[:-3]}.json'
    SAMPLES_CSV = 'samples.csv'
    R1_SUFFIX = '_R1.fastq.gz'
    R2_SUFFIX = '_R2.fastq.gz'

    inputs: Dict[str, Any]
    df: pd.DataFrame

    def main(self):
        self.read_files()
        self.fill_in_samples()
        self.write_json()

    def read_files(self):
        self.inputs = json.load(open(self.TEMPLATE_JSON))
        self.df = pd.read_csv(self.SAMPLES_CSV)

    def fill_in_samples(self):
        normal_names = self.df['Normal Sample Name'].tolist()
        normal_fastqs = [
            [f'{x}{self.R1_SUFFIX}', f'{x}{self.R2_SUFFIX}'] for x in normal_names
        ]

        tumor_names = self.df['Tumor Sample Name'].tolist()
        tumor_fastqs = [
            [f'{x}{self.R1_SUFFIX}', f'{x}{self.R2_SUFFIX}'] for x in tumor_names
        ]

        self.inputs['SomaticTNPaired.normalSampleNames'] = normal_names
        self.inputs['SomaticTNPaired.inFileNormalFastqPairs'] = normal_fastqs
        self.inputs['SomaticTNPaired.tumorSampleNames'] = tumor_names
        self.inputs['SomaticTNPaired.inFileTumorFastqPairs'] = tumor_fastqs

        self.inputs['SomaticTNPaired.refIntervalBeds'] = self.df['BED File'].tolist()

    def write_json(self):
        with open('inputs.json', 'w') as f:
            json.dump(self.inputs, f, indent=4)


if __name__ == '__main__':
    Main().main()
