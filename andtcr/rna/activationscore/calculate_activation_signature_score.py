from optparse import OptionParser, make_option
import os

import numpy as np
import pandas as pd


class ActivationOptionParser(OptionParser):
    options = [
        make_option('-f', '--file_name', action='store',
                    type='string', dest='file_name',
                    help='Path to expression file for processing.'),
    ]

    id_column = 'gene symbol'

    def __init__(self, option_list=None, **kwargs):
        OptionParser.__init__(self,
                              option_list=option_list or self.options,
                              **kwargs)

    def get_signature(self):
        dir = os.path.dirname(__file__)
        data_file = os.path.join(dir, 'activation_signature.txt')
        data = pd.io.parsers.read_csv(data_file,
                                      sep='\t', header=0)
        data = data.fillna(0)
        return data

    def get_input_data(self, filename):
        data = pd.io.parsers.read_csv(filename,
                                      sep='\t', header=0)

        # Error check.
        colnames = [s.lower() for s in data.columns]
        data.columns = colnames

        try:
            data.index = data['gene symbol'].str.lower()
            data = data.drop('gene symbol', axis=1)
            self.id_column = 'gene symbol'
        except KeyError:
            try:
                data.index = data['refseq'].str.lower()
                data = data.drop('refseq', axis=1)
                self.id_column = 'refseq'
            except KeyError:
                raise Exception('Columns must contain either "gene symbol" ' +
                                'or "refseq." Make sure one of these two ' +
                                'identifier columns is present and that the ' +
                                'header is correctly labeled.')

        # Replace strings
        def replace_strings(val):
            try:
                val < 0
            except TypeError:
                try:
                    val = float(val)
                except ValueError:
                    val = 0
            return val

        for val_col in data.columns:
            data[val_col] = data[val_col].map(replace_strings)

        data = data.fillna(0)
        return data

    def normalize_input(self, input_data):
        '''
        Mean-center input data.
        '''
        input_data = input_data.sub(input_data.mean(axis=1), axis=0)

        return input_data

    def merge_data(self, sig_data, input_data, filename):
        '''
        Given our two data sets, we want to extract the activation 
        signature genes from the input data and return the rows of interest
        for further processing. We join on either Refseq or Gene symbol.
        '''
        sig_data.index = sig_data[self.id_column].str.lower()
        merged = pd.merge(sig_data[['pc1']], input_data,
                          how='inner', left_index=True,
                          right_index=True)

        if len(merged) < len(sig_data):
            print('Warning! {} of {} signature genes found.'.format(
                len(merged), len(sig_data)))

        return merged

    def calculate_scores(self, merged):
        '''
        Use the PC1 value to calculate the total score for the input data.
        '''
        data_cols = merged.columns.tolist()
        data_cols.remove('pc1')
        scores = np.dot(merged['pc1'], merged[data_cols])
        return pd.Series(scores, index=data_cols)

    def scores_from_file(self, filename):
        sig_data = self.get_signature()
        input_data = self.get_input_data(filename)

        input_data = self.normalize_input(input_data)
        merged = self.merge_data(sig_data, input_data, filename)

        return self.calculate_scores(merged)

if __name__ == '__main__':
    parser = ActivationOptionParser()
    options, args = parser.parse_args()

    score = parser.scores_from_file(options.file_name)

    print('Activation scores:\n{}'.format(score))
