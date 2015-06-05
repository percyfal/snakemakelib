# Copyright (C) 2015 by Per Unneberg
import os
import re
import pandas as pd
from snakemakelib.exceptions import SamplesException, OutputFilesException
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

class Results(dict):
    _inputfiles = []
    _samples = []
    _keys = ()

    def __init__(self, inputs=(), samples=(), keys=(), *args, **kw):
        """Init Results.

        Args:
          inputs (tuple): tuple consisting of (inputfilename, sample)
            tuple pairs or input filenames; if the latter, require that
            samples be present
          samples (tuple): tuple of sample names
          keys (tuple): tuple of keys

        """
        super(Results, self).__init__()
        for x in inputs:
            if isinstance(x, tuple):
                (inputfile, sample) = x
                self._inputfiles.append(inputfile)
                self._samples.append(sample)
            else:
                self._inputfiles = inputs
                self._samples = samples
                break
        if (len(self._samples) != len(self._inputfiles)):
            raise SamplesException("sample list must be as long as input file list")
        self._keys = keys
        for k in self._keys:
            self[k] = pd.DataFrame()
        self._collect_results()

    def _collect_results(self):
        raise NotImplementedError(__name__ + ": abstract base class Result: implement in subclass")

    def save(self, outputfiles, **kw):
        """Save data frames to outputfiles.

        Results keys will be sorted by name so that implicitly the
        outputfiles correspond to the key names in sorted orde

        Args:
          outputfiles (list): list of output file names. Length must correspond to number of keys.
          kw (dict): keyword arguments to pass to pandas DataFrame.to_csv function
        
        """
        if (len(self.keys()) != len(outputfiles)):
            raise OutputFilesException("wrong number of outputfiles; must be equal to number of keys: {klen}".format(klen=len(self)))
        for (key, ofile) in zip(sorted(list(self.keys())), outputfiles):
            smllogger.info("saving data frame {df}".format(df=key))
            df = self.get(key, None)
            if df is None:
                smllogger.warn("No data for {df}; skipping".format(df=key))
                continue
            df.to_csv(ofile, **kw)

    def load(self, infile, sep="\t", newline="\n", comment_char="#"):
        """Load a data file.

        Args:
          infile (str): input file name
          sep (str): field separator
          newline (str): newline character
          comment_char (str): comment character

        """
        with open(infile) as fh:
            data = [x.strip(newline).split(sep) for x in fh.readlines() if not x.strip() == ""]
        return data

    def load_data_frame(self, infile, default_loader="read_csv", **kwargs):
        """Load a data frame.

        Use pandas interface directly to read input file.

        Args:
          infile (str): input file name
          default_loader (str): explicitly tell pandas what loader function to use
          kwargs: keyword arguments to pass to pandas

        """
        try:
            ext = os.path.splitext(infile)
            loader = getattr(pd, "read_" + ext.rstrip("."))
        except:
            smllogger.warn("Couldn't determine file format; falling back on default loader '{}'".format(default_loader))
            loader = getattr(pd, default_loader)
        return loader(infile, **kwargs)

    def parse_data(self, data, comment_char="#", rs=(None, None), skip=0):
        """Parse a data set.

        Args:
          data (list): list of data lines split into columns
          comment_char (str): comment character
          rs (tuple): record separator, a 2-tuple (start_string,
            end_string) inbetween which data lines will be parsed
          skip (int): skip number of lines after first record separator match before reading

        Returns:
          DataFrame: pandas DataFrame object with parsed data

        """
        indices = (0 if rs[0] is None else next((i for i in range(len(data)) if rs[0] in data[i]), -1),
                   len(data) if rs[1] is None else next((i for i in range(len(data)) if rs[0] in data[i]), -1))
        return pd.DataFrame (data[indices[0] + skip:indices[1]])
