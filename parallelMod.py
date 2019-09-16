from mpi4py import MPI
import numpy as np
import math

class MpiConfig:
    """
    Abstract class for defining the MPI parameters,
    along with initialization of the MPI communication
    handle from mpi4py.
    """
    def __init__(self):
        """
        Initialize the MPI abstract class that will contain basic
        information and communication handles.
        """
        self.comm = None
        self.rank = None
        self.size = None
        self.dInd = None
        self.bInd = None

    def initialize_comm(self):
        """
        Initial function to initialize MPI.
        :return:
        """
        try:
            self.comm = MPI.COMM_WORLD
        except:
            print("Unable to initialize the MPI Communicator object")
            raise Exception()

        try:
            self.size = self.comm.Get_size()
        except:
            print("Unable to retrieve the MPI size.")
            raise  Exception()

        try:
            self.rank = self.comm.Get_rank()
        except:
            print("Unable to retrieve the MPI processor rank")
            raise Exception()

    def calc_boundaries(self, wbObj):
        """
        Function to break up the processing amongst the different processors.
        Processing will be broken up depending on the time period, and number of basins.
        :param wbObj:
        :return:
        """
        # First, multiply the total number of steps by the number of basins. This gives us
        # the total number of steps to process.
        all_steps = wbObj.nGlobalSteps * wbObj.nGlobalBasins

        # Next, divide by the number of processors to get a local number of steps to process.
        local_steps_constant = math.floor(all_steps / self.size)

        # Calculate the remainder steps, these will be processed by the last rank.
        remainder = all_steps % self.size

        if self.rank == (self.size - 1):
            local_steps = local_steps_constant + remainder
        else:
            local_steps = local_steps_constant

        begIndGlobal = local_steps_constant * self.rank
        endIndGlobal = begIndGlobal + local_steps

        # Loop through all of the local steps, assign the date index, and
        # basin index locally. This will help guide analysis.
        self.dInd = np.empty([local_steps], np.int32)
        self.bInd = np.empty([local_steps], np.int32)

        for stepTmp in range(0,local_steps):
            indCurrent = begIndGlobal + stepTmp
            self.bInd[stepTmp] = math.floor(indCurrent / wbObj.nGlobalSteps)
            self.dInd[stepTmp] = indCurrent % wbObj.nGlobalSteps

        # Initialize local water balance arrays that will be populated.
        wbObj.accEcanLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.accEdirLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.accEtranLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.accPrcpLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.accSneqLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.canIceLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.canLiqLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.sfcRnoffLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.sfcHeadSubRtLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.uGrdRnoffLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.soilMLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.streamVolLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.qLatVolLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.qbdryRtLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.gwInLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.gwOutLocal = np.full([local_steps], -9999.0, np.float64)
        wbObj.zLevLocal = np.full([local_steps], -9999.0, np.float64)
