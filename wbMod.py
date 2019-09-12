from netCDF4 import Dataset
import numpy as np
from skimage.transform import downscale_local_mean
import os

class wbObj:
    """
    Object to contain information on basins being processed, date ranges, etc.
    """
    def __init__(self):
        """
        Initialize the object class that will contain information.
        """
        self.geoPath = None
        self.modelDir = None
        self.outDir = None
        self.fullDomPath = None
        self.rtLinkPath = None
        self.spWtPath = None
        self.nxHydro = None
        self.nyHydro = None
        self.nxLand = None
        self.nyLand = None
        self.basinsGlobal = None
        self.nGlobalBasins = None
        self.nGlobalSteps = None
        self.linksGlobal = None
        self.bDateGlobal = None
        self.eDateGlobal = None
        self.streamVolLocal = None
        self.qLatVolLocal = None
        self.accPrcpLocal = None
        self.accEcanLocal = None
        self.accEtranLocal = None
        self.accEdirLocal = None
        self.accSneqLocal = None
        self.canIceLocal = None
        self.canLiqLocal = None
        self.sfcRnoffLocal = None
        self.uGrdRnoffLocal = None
        self.soilMLocal = None
        self.sfcHeadSubRtLocal = None
        self.qbdryRtLocal = None
        self.qStrmVolRtLocal = None
        self.gwOutLocal = None
        self.zLevLocal = None
        self.geoRes = None
        self.hydRes = None
        self.aggFact = None
        self.upstreamLinks = {}
        self.bsnMskLand = {}
        self.bsnMskHydro = {}
        #self.bsnMskLandArea = {}
        #self.bsnMskHydroArea = {}
        self.gageIDs = {}
        self.linksLocal = {}

    def calcGeoParams(self, MpiConfig):
        """
        Function that calculates all geospatial parameters for each streamflow location
        in order to properly calculate the water budget.
        :return:
        """
        try:
            idGeo = Dataset(self.geoPath,'r')
        except:
            print("Unable to open: " + self.geoPath)
            raise Exception()

        try:
            idFullDom = Dataset(self.fullDomPath,'r')
        except:
            print("Unable to open: " + self.fullDomPath)
            raise Exception()

        try:
            idRt = Dataset(self.rtLinkPath, 'r')
        except:
            print("Unable to open: " + self.rtLinkPath)
            raise Exception()

        try:
            idWt = Dataset(self.spWtPath, 'r')
        except:
            print("Unable to open: " + self.spWtPath)
            raise Exception()

        toVar = idRt.variables['to'][:].data
        linkVar = idRt.variables['link'][:].data
        idMaskVar = idWt.variables['IDmask'][:].data
        regridWeightVar = idWt.variables['regridweight'][:].data
        iVar = idWt.variables['i_index'][:].data
        jVar = idWt.variables['j_index'][:].data

        self.geoRes = idGeo.DX
        # Hard-code for now...
        self.hydRes = 250.0
        self.aggFact = self.geoRes / self.hydRes
        self.nxHydro = idFullDom.variables['x'].shape[0]
        self.nyHydro = idFullDom.variables['y'].shape[0]
        self.nxLand = idGeo.variables['XLAT_M'].shape[2]
        self.nyLand = idGeo.variables['XLAT_M'].shape[1]

        # Process each geospatial info for the gages. Store into a dictionary object.
        for bsnTmp in np.unique(MpiConfig.bInd):
            # Fetch the gage ID and link from the water budget object.
            self.gageIDs[bsnTmp] = self.basinsGlobal[bsnTmp]
            self.linksLocal[bsnTmp] = self.linksGlobal[bsnTmp]

            # First, we will calculate all of the upstream links for this particular gage
            # site.
            upLinks = np.full([1],self.linksGlobal[bsnTmp])
            doneLinks = np.full([0],0)

            while len(upLinks) > 0:
                i = upLinks[0]
                doneLinks = np.append(doneLinks, i)
                newLinks = linkVar[np.where(toVar == i)]
                upLinks = np.append(upLinks, newLinks)
                upLinks = upLinks[np.in1d(upLinks, doneLinks, invert=True)]

            self.upstreamLinks[bsnTmp] = doneLinks

            # Extract which polygons in the weight file are associated with the
            # upstream links.
            polyInd = np.in1d(idMaskVar, doneLinks, invert=False)
            regridSub = regridWeightVar[polyInd]
            iVarSub = iVar[polyInd]
            jVarSub = jVar[polyInd]

            # Create a mask grid on the hydro grid.
            self.bsnMskHydro[bsnTmp] = np.full([self.nyHydro, self.nxHydro], 0.0)

            # Loop over all unique I/J values for this basin. Sum up the weight values
            # and place them into the mask grid.
            for stepTmp in range(0, len(regridSub)):
                iTmp = iVarSub[stepTmp]
                jTmp = jVarSub[stepTmp]
                self.bsnMskHydro[bsnTmp][jTmp,iTmp] = self.bsnMskHydro[bsnTmp][jTmp,iTmp] +\
                                                      regridSub[stepTmp]

            # Use the scikit image processing to resample
            self.bsnMskLand[bsnTmp] = downscale_local_mean(self.bsnMskHydro[bsnTmp], (int(self.aggFact), int(self.aggFact)))

            # Reset temporary arrays for next basin tracing.
            doneLinks = None
            polyInd = None
            regridSub = None
            iVarSub = None
            jVarSub = None

            # Calculate basin area in squared meters.
            #self.bsnMskLandArea[bsnTmp] = self.bsnMskLand[bsnTmp] * (self.geoRes * self.geoRes)
            #self.bsnMskLandArea[bsnTmp] = self.bsnMskLandArea.sum()
            #self.bsnMskHydroArea[bsnTmp] = self.bsnMskHydro[bsnTmp] * (self.hydRes * self.hydRes)
            #self.bsnMskHydroArea[bsnTmp] = self.bsnMskHydroArea.sum()

        # Close the NetCDF files and reset variables for memory purposes
        idFullDom.close()
        idRt.close()
        idWt.close()
        idGeo.close()

        toVar = None
        linkVar = None
        idMaskVar = None
        regridWeightVar = None
        iVar = None
        jVar = None

    def readLdasOut(self, stepCurrent, dCurrent, bCurrent, MpiConfig):
        """
        Function to read in LDASOUT water balance related variables and aggregate to the basin using the
        pre-calculated masks on the land grid.
        :return:
        """
        modelPath = self.modelDir + "/" + dCurrent.strftime('%Y%m%d%H00') + ".LDASOUT_DOMAIN1"

        #if MpiConfig.rank == 0:
        #    print(modelPath)

        # If the file is not present, this may not indicate an issue, but that we are only producing
        # LDASOUT files at an infrequent time period. Simply return to the main calling program and leave
        # values in our local arrays to missing.
        if not os.path.isfile(modelPath):
            return

        # Open the NetCDF file.
        idLdas = Dataset(modelPath, 'r')

        # Read in SWE and aggregate to the basin.
        varTmp = idLdas.variables['SNEQV'][0,:,:]
        indTmp = np.where(varTmp >= 0.0)
        varTmp[indTmp] = varTmp[indTmp] * (self.geoRes * self.geoRes) # Convert to cubic meters
        varTmp = varTmp * self.bsnMskLand[bCurrent]
        self.accSneqLocal[stepCurrent] = varTmp.sum() # Volume of cubic meters.
        varTmp = None

        # Close the NetCDF file.
        idLdas.close()

    def outputWb(self, MpiConfig):
        """
        Output routine that collects water balance variables from various processors, then breaks things up
        by basin into a final NetCDF file.
        :param MpiConfig:
        :return:
        """
        outPath = self.outDir + "/WaterBudget_" + self.bDateGlobal.strftime('%Y%m%d%H') + "_" + \
                  self.eDateGlobal.strftime('%Y%m%d%H') + ".nc"

        if MpiConfig.rank == 0:
            if os.path.isfile(outPath):
                print("Water budget output file: " + outPath + " already exists.")
                raise Exception()

            # Create ouptut file.
            idOut = Dataset(outPath, 'w')

            idOut.createDimension('numBasins', len(self.basinsGlobal))
            idOut.createDimension('numSteps', self.nGlobalSteps)

            idOut.createVariable('SWE_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)

        # Collect arrays
        final = MpiConfig.comm.gather(self.accSneqLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['SWE_Volume'][bTmp,:] = dataOutTmp[bIndTmp:eIndTmp]

        # Close the output netCDF file.
        if MpiConfig.rank == 0:
            idOut.close()
