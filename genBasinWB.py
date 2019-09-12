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
    parser.add_argument('out_dir', metavar='out_dir', type=str, nargs='+',
                        help='Directory containing all of the output files.')
    parser.add_argument('geoGrid', metavar='geoGrid', type=str, nargs='+',
                        help='Geogrid file defining the land surface modeling grid.')
    parser.add_argument('fullDom', metavar='fullDom', type=str, nargs='+',
                        help='Fulldom file defining the 2D routing grid.')
    parser.add_argument('rtLink', metavar='rtLink', type=str, nargs='+',
                        help='Route link path.')
    parser.add_argument('spWtFile', metavar='spWtFile', type=str, nargs='+',
                        help='Spatial weight file used in generation of basin masks.')


    # Process the input arguments into the program.
    args = parser.parse_args()

    # Make sure arguments (files/directories) exist for sanity checking.
    if not os.path.isfile(args.gage_list[0]):
        print("Expected input list of gages to process: " + args.gage_list[0] + " not found.")
        sys.exit(-1)

    if not os.path.isdir(args.force_dir[0]):
        print("Expected forcing directory for the NWM simulation: " + args.force_dir[0] + " not found.")
        sys.exit(-1)

    if not os.path.isdir(args.model_dir[0]):
        print("Expected model output directory: " + args.model_dir[0] + " not found.")
        sys.exit(-1)

    if not os.path.isdir(args.out_dir[0]):
        print("Expected output directory to hold water balance output files: " + args.out_dir[0] + " not found.")
        sys.exit(-1)

    if not os.path.isfile(args.geoGrid[0]):
        print("Expected geogrid file: " + args.geoGrid[0] + " not found.")
        sys.exit(1)

    if not os.path.isfile(args.fullDom[0]):
        print("Expected fullDom file: " + args.fullDom[0] + " not found.")
        sys.exit(1)

    if not os.path.isfile(args.rtLink[0]):
        print("Expected route link file: " + args.rtLink[0] + " not found.")
        sys.exit(1)

    if not os.path.isfile(args.spWtFile[0]):
        print("Expected spatial weight file: " + args.spWtFile[0] + " not found.")
        sys.exit(1)

    # Initialize the MPI objects necessary to parallize the processing.
    mpiMeta = parallelMod.MpiConfig()
    try:
        mpiMeta.initialize_comm()
    except:
        print("Unable to launch MPI objects.")
        sys.exit(1)

    # Initialize the water balance object.
    wb_data = wbMod.wbObj()

    # Assign file paths.
    wb_data.fullDomPath = args.fullDom[0]
    wb_data.geoPath = args.geoGrid[0]
    wb_data.rtLinkPath = args.rtLink[0]
    wb_data.spWtPath = args.spWtFile[0]
    wb_data.modelDir = args.model_dir[0]
    wb_data.outDir = args.out_dir[0]

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

    # Initialize local information on which boundaries and basins to process. Local
    # arrays holding water budget variables will also be initialized.
    try:
        mpiMeta.calc_boundaries(wb_data)
    except:
        print("Unable to calculate local basins and dates")
        mpiMeta.comm.Abort()
        sys.exit(1)

    # Perform upstream tracing and calculation of masks on the hydro and land grids.
    wb_data.calcGeoParams(mpiMeta)

    # Begin looping through each time step, read in the land/hydro/gw variables, aggregate 2D fields using the
    # basin masks, and store data in the local arrays.
    for step in range(mpiMeta.bInd.shape[0]):
        dCurrent = wb_data.bDateGlobal + datetime.timedelta(seconds=3600*int(mpiMeta.dInd[step]))
        basinCurrent = mpiMeta.bInd[step]

        # Read in LDASOUT variables and aggregate to the basins.

        # Place output into a final NetCDF file.
        wb_data.readLdasOut(step, dCurrent, basinCurrent, mpiMeta)

    # Create output file containing final data.
    wb_data.outputWb(mpiMeta)


if __name__ == "__main__":
    main()