#!/usr/bin/env python
"""Pack ASCII injection files into a psd.xml file."""

import os
from argparse import SUPPRESS,  FileType

import glue.ligolw.utils
import lal
import lal.series
import numpy as np
from ligo.skymap.tool import ArgumentParser, register_to_xmldoc


def read_o5_txt_file(filepath, config_column):
    """
    Read O5 format .txt file and extract one configuration.

    Parameters:
        filepath: Path to the .txt file
        config_column: Which column to read (1=O5a, 2=O5b, 3=O5c)

    Returns:
        Tuple: (frequency, asd_data)
    """
    data = np.loadtxt(filepath, skiprows=1)
    frequency = data[:, 0]
    asd_data = data[:, config_column]
    return frequency, asd_data


#  Config columns mapping (LIGO O5 configurations)
config_columns = {
    "O5a": 1,  # Column 1 = O5aStrain
    "O5b": 2,  # Column 2 = O5bStrain
    "O5c": 3,  # Column 3 = O5cStrain
}



# Command line interface
detector_names = [d.frDetector.prefix for d in lal.CachedDetectors]
detector_long_names = [d.frDetector.name for d in lal.CachedDetectors]
parser = ArgumentParser()
parser.add_argument(
    "-o",
    "--output",
    metavar="OUT.xml[.gz]",
    type=FileType("wb"),
    default="-",
    help="Name of output file [default: stdout]",
)
for name, long_name in zip(detector_names, detector_long_names):
    parser.add_argument(
        "--" + name,
        metavar="PSD.txt",
        type=FileType("r"),
        default=SUPPRESS,
        help="PSD function for {0} detector".format(long_name),
    )

parser.add_argument(
    "--config",
    metavar="CONFIG",
    type=str,
    default=None,
    help="Configuration name (O5a/O5b/O5c for O5 mode)",
)

args = parser.parse_args()

psds = {}

# Process each detector
for name in detector_names:
    psd_file = getattr(args, name, None)
    if psd_file is None:
        continue

    # Check if O5 mode: config in config_columns and detector is H1 or L1
    if args.config in config_columns and name in ["H1", "L1"]:
        # Read specific column
        column = config_columns[args.config]
        f, asd = read_o5_txt_file(psd_file, column)
    else:
        # Standard mode: read normally
        f, asd = np.loadtxt(psd_file).T

    # Process PSD data
    psd = np.square(asd)

    f0 = 10.0
    fmax = 4096.0
    df = 1.0

    fgrid = np.arange(f0, fmax, df)
    series = lal.CreateREAL8FrequencySeries(
        (psd_file.name).split(".")[0], 0, f0, df, lal.SecondUnit, len(fgrid)
    )
    series.data.data = np.exp(np.interp(np.log(fgrid), np.log(f), np.log(psd)))

    psds[name] = series

xmldoc = lal.series.make_psd_xmldoc(psds)
register_to_xmldoc(xmldoc, parser, args)

with glue.ligolw.utils.SignalsTrap():
    glue.ligolw.utils.write_fileobj(
        xmldoc, args.output, gz=(os.path.splitext(args.output.name)[-1] == ".gz")
    )
