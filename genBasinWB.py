import argparse
import sys
import os
import parallelMod
import wbMod
import datetime
import pandas as pd

def main():
    """
    Main calling program to execute water balance calculations from an NWM simulation. Necessary arguments are:
    1.) A list of gage ID's, with associated NHD COM IDs. If no gage ID, number appriately. Must be in CSV format...
    2.) A path to the forcing directory.
    3.) A path to the output directory where NWM output resides.
    4.) An output directory to place the final output file, which will be in NetCDF format.
    :return:
    """
    # Parse out the configuration file.
    parser = argparse.ArgumentParser(description='Main calling program to run water budget analysis on NWM output')
    parser.add_argument('gage_list', metavar='gage_list', type=str, nargs='+',
                        help='CSV file with list of gage IDs along with NHD comID values to run analysis on')
    parser.add_argument('begDate', metavar='begDate', type=int, nargs='+',
                        help='Beginning date of analysis in YYYYMMDDHH format')
    parser.add_argument('endDate', metavar='endDate', type=int, nargs='+',
                        help='Ending date of analysis in YYYYMMDDHH format')
    parser.add_argument('force_dir', metavar='force_dir', type=str, nargs='+',
                        help='Directory containing forcing files. This is necessary to compute the inputs into '
                             'the basin')
    parser.add_argument('model_dir', metavar='model_dir', type=str, nargs='+',
                        help='Directory containing necessary NWM output files to be read in')
    parser.add_argument('out_dir', metavar='output directory where water balance values will be outputted to '
                                           'a NetCDF file')

    # Process the input arguments into the program.
    args = parser.parse_args()

    # Make sure arguments (files/directories) exist for sanity checking.
    if not os.path.isfile(args.gage_list[0]):
        print("Expected input list of gages to process: " + args.gage_list[0] + " not found.")
        sys.exit(-1)

    if not os.path.isdir(args.force_dir[0]):
        print("Expected forcing directory for the NWM simulation: " + args.force_dir[0] + " not found.")
        sys.exit(-1)

    if not os.path.isdir(args.out_dir[0]):
        print("Expected output directory to hold water balance output files: " + args.out_dir[0] + " not found.")
        sys.exit(-1)

    # Initialize the MPI objects necessary to parallize the processing.
    mpiMeta = parallelMod.MpiConfig()
    try:
        mpiMeta.initialize_comm()
    except:
        print("Unable to launch MPI objects.")
        sys.exit(1)

    # Initialize the water balance object.
    wb_data = wbMod.wbObj()

    # Initialize datetime objects
    try:
        wb_data.bDateGlobal = datetime.datetime.strptime(str(args.begDate[0]), '%Y%m%d%H')
    except:
        print("Unable to initialize the global beginning date on rank: " + str(mpiMeta.rank))
        mpiMeta.comm.Abort()
        sys.exit(1)

    try:
        wb_data.eDateGlobal = datetime.datetime.strptime(str(args.endDate[0]), '%Y%m%d%H')
    except:
        print("Unable to initialize the global ending date on rank: " + str(mpiMeta.rank))
        mpiMeta.comm.Abort()
        sys.exit(1)

    # Read in the list of gages to process.
    try:
        gagesTmp = pd.read_csv(args.gage_list[0])
    except:
        print("Unable to read in CSV file: " + str(args.gage_list[0]) + " on rank: " + str(mpiMeta.rank))
        mpiMeta.comm.Abort()
        sys.exit(1)

    if "gage_id" not in gagesTmp.keys():
        print("Expected column header 'gage_id' not found in: " + str(args.gage_list[0]))
        mpiMeta.comm.Abort()
        sys.exit(1)

    if "link" not in gagesTmp.keys():
        print("Expected column header 'link' not found in: " + str(args.gage_list[0]))
        mpiMeta.comm.Abort()
        sys.exit(1)

    wb_data.basinsGlobal = gagesTmp.gage_id
    wb_data.linksGlobal = gagesTmp.link
    wb_data.nGlobalBasins = len(gagesTmp.gage_id)
    dtTmp = wb_data.eDateGlobal - wb_data.bDateGlobal
    wb_data.nGlobalSteps = int(dtTmp.days * 24) + int(dtTmp.seconds/3600.0)

    # Reset variables for memory purposes.
    gagesTmp = None
    dtTmp = None

    # Initialize local and global arrays. Global arrays are only calculated on
    # rank 0 as this is sent to the output files.



if __name__ == "__main__":
    main()