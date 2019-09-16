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
        self.gwInLocal = None
        self.gwOutLocal = None
        self.zLevLocal = None
        self.geoRes = None
        self.hydRes = None
        self.aggFact = None
        self.soilDepths = [0.1, 0.3, 0.6, 1.0]
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

        if MpiConfig.rank == 0:
            print(modelPath)

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
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.geoRes * self.geoRes) # Convert to cubic meters
        varTmp = varTmp * self.bsnMskLand[bCurrent]
        self.accSneqLocal[stepCurrent] = varTmp.sum() # Volume of cubic meters.
        varTmp = None

        # Read in accumulated precipitation and aggregate to the basin.
        varTmp = idLdas.variables['ACCPRCP'][0,:,:]
        indTmp = np.where(varTmp >= 0.0)
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.geoRes * self.geoRes)  # Convert to cubic meters
        varTmp = varTmp * self.bsnMskLand[bCurrent]
        self.accPrcpLocal[stepCurrent] = varTmp.sum()
        varTmp = None

        # Read in accumulated canopy evaporation and aggregate to the basin.
        varTmp = idLdas.variables['ACCECAN'][0, :, :]
        indTmp = np.where(varTmp >= 0.0)
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.geoRes * self.geoRes)  # Convert to cubic meters
        varTmp = varTmp * self.bsnMskLand[bCurrent]
        self.accEcanLocal[stepCurrent] = varTmp.sum()
        varTmp = None

        # Read in accumulated evapotranspiration and aggregate to the basin.
        varTmp = idLdas.variables['ACCETRAN'][0, :, :]
        indTmp = np.where(varTmp >= 0.0)
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.geoRes * self.geoRes)  # Convert to cubic meters
        varTmp = varTmp * self.bsnMskLand[bCurrent]
        self.accEtranLocal[stepCurrent] = varTmp.sum()
        varTmp = None

        # Read in accumulated ??????? and aggregate to the basin.
        varTmp = idLdas.variables['ACCEDIR'][0, :, :]
        indTmp = np.where(varTmp >= 0.0)
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.geoRes * self.geoRes)  # Convert to cubic meters
        varTmp = varTmp * self.bsnMskLand[bCurrent]
        self.accEdirLocal[stepCurrent] = varTmp.sum()
        varTmp = None

        # Read in canopy ice and aggregate to the basin.
        varTmp = idLdas.variables['CANICE'][0, :, :]
        indTmp = np.where(varTmp >= 0.0)
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.geoRes * self.geoRes)  # Convert to cubic meters
        varTmp = varTmp * self.bsnMskLand[bCurrent]
        self.canIceLocal[stepCurrent] = varTmp.sum()
        varTmp = None

        # Read in canopy liquid and aggregate to the basin.
        varTmp = idLdas.variables['CANLIQ'][0, :, :]
        indTmp = np.where(varTmp >= 0.0)
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.geoRes * self.geoRes)  # Convert to cubic meters
        varTmp = varTmp * self.bsnMskLand[bCurrent]
        self.canLiqLocal[stepCurrent] = varTmp.sum()
        varTmp = None

        # Read in surface runoff and aggregate to the basin.
        varTmp = idLdas.variables['SFCRNOFF'][0, :, :]
        indTmp = np.where(varTmp >= 0.0)
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.geoRes * self.geoRes)  # Convert to cubic meters
        varTmp = varTmp * self.bsnMskLand[bCurrent]
        self.sfcRnoffLocal[stepCurrent] = varTmp.sum()
        varTmp = None

        # Read in underground runoff runoff and aggregate to the basin.
        varTmp = idLdas.variables['UGDRNOFF'][0, :, :]
        indTmp = np.where(varTmp >= 0.0)
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.geoRes * self.geoRes)  # Convert to cubic meters
        varTmp = varTmp * self.bsnMskLand[bCurrent]
        self.uGrdRnoffLocal[stepCurrent] = varTmp.sum()
        varTmp = None

        # Read in and aggregate soil moisture from the various layers.
        #varTmp = idLdas.variables['SOIL_M'][:, 0, :, :]
        self.soilMLocal[stepCurrent] = 0.0
        for lyrTmp in range(4):
            varTmp = idLdas.variables['SOIL_M'][:, 0, lyrTmp, :]
            indTmp = np.where(varTmp >= 0.0)
            varTmp[indTmp] = varTmp[indTmp] * self.soilDepths[lyrTmp]
            varTmp[indTmp] = varTmp[indTmp] * (self.geoRes * self.geoRes)  # Convert to cubic meters
            self.soilMLocal[stepCurrent] = self.soilMLocal[stepCurrent] + varTmp.sum()
            varTmp = None


        # Close the NetCDF file.
        idLdas.close()

    def readRtOut(self, stepCurrent, dCurrent, bCurrent, MpiConfig):
        """
        Function to read in RTOUT files that contain distributed surface head, and accumulated
        spatial variables on the routing grid.
        :param stepCurrent:
        :param dCurrent:
        :param bCurrent:
        :param MpiConfig:
        :return:
        """
        modelPath = self.modelDir + "/" + dCurrent.strftime('%Y%m%d%H00') + ".RTOUT_DOMAIN1"

        # if MpiConfig.rank == 0:
        #    print(modelPath)

        # If the file is not present, this may not indicate an issue, but that we are only producing
        # RTOUT files at an infrequent time period. Simply return to the main calling program and leave
        # values in our local arrays to missing.
        if not os.path.isfile(modelPath):
            return

        # Open the NetCDF file.
        idRt = Dataset(modelPath, 'r')

        # Read in surface head and aggregate to the basin.
        varTmp = idRt.variables['sfcheadsubrt'][0, :, :]
        indTmp = np.where(varTmp >= 0.0)
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.hydRes * self.hydRes)  # Convert to cubic meters
        varTmp = varTmp * self.bsnMskHydro[bCurrent]
        self.sfcHeadSubRtLocal[stepCurrent] = varTmp.sum()  # Volume of cubic meters.
        varTmp = None

        # Read in QBDRYRT and aggregate to the basin.
        varTmp = idRt.variables['QBDRYRT'][0, :, :]
        indTmp = np.where(varTmp != varTmp.fill_value)
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.hydRes * self.hydRes)  # Convert to cubic meters
        varTmp = varTmp * self.bsnMskHydro[bCurrent]
        self.qbdryRtLocal[stepCurrent] = varTmp.sum()  # Volume of cubic meters.
        varTmp = None

        # Read in QSTRMVOLRT and aggregate to the basin.
        varTmp = idRt.variables['QSTRMVOLRT'][0, :, :]
        indTmp = np.where(varTmp >= 0.0)
        varTmp[indTmp] = varTmp[indTmp] / 1000.0  # Convert from mm to meters
        varTmp[indTmp] = varTmp[indTmp] * (self.hydRes * self.hydRes)  # Convert to cubic meters
        varTmp = varTmp * self.bsnMskHydro[bCurrent]
        self.qStrmVolRtLocal[stepCurrent] = varTmp.sum()  # Volume of cubic meters.
        varTmp = None

        # Close the NetCDF file
        idRt.close()

    def readGwOut(self, stepCurrent, dCurrent, bCurrent, MpiConfig):
        """
        Function to read inflow and outflow variables from groundwater output files, aggregated
        to the basin using the uplinks calculated earlier in the program.
        :param MpiConfig:
        :return:
        """
        modelPath = self.modelDir + "/" + dCurrent.strftime('%Y%m%d%H00') + ".GWOUT_DOMAIN1"

        # if MpiConfig.rank == 0:
        #    print(modelPath)

        # If the file is not present, this may not indicate an issue, but that we are only producing
        # GWOUT files at an infrequent time period. Simply return to the main calling program and leave
        # values in our local arrays to missing.
        if not os.path.isfile(modelPath):
            return

        # Open the NetCDF file.
        idGw = Dataset(modelPath, 'r')

        # Read in the inflow/outflow variables, along with the feature_id variable.
        # Using the pre-calculated uplinks, we will subset the inflow/outflow volumes,
        # then sum up total volumes for the entire basin.
        varTmp = idGw.variables['feature_id'][:]
        gwInd = np.in1d(varTmp, self.upstreamLinks[bCurrent], invert=True)
        varTmp = None

        # Read in GW inflow and aggregate to the basin.
        varTmp = idGw.variables['inflow'][:]
        varTmp = varTmp[gwInd]
        indTmp = np.where(varTmp != varTmp.fill_value)
        varTmp = varTmp[indTmp]
        varTmp = varTmp * 3600.0 # Converting m^3/s to cubic meters.
        self.gwInLocal[stepCurrent] = varTmp.sum()  # Volume of cubic meters.
        varTmp = None

        # Read in GW outflow and aggregate to the basin.
        varTmp = idGw.variables['outflow'][:]
        varTmp = varTmp[gwInd]
        indTmp = np.where(varTmp != varTmp.fill_value)
        varTmp = varTmp[indTmp]
        varTmp = varTmp * 3600.0  # Converting m^3/s to cubic meters.
        self.gwOutLocal[stepCurrent] = varTmp.sum()  # Volume of cubic meters.
        varTmp = None

        # Close the NetCDF file.
        idGw.close()

    def readChrtout(self, stepCurrent, dCurrent, bCurrent, MpiConfig):
        """
        Function that will read in hourly streamflow for a given basin and convert to a water volume.
        Assuming hourly outputs here.
        :param stepCurrent:
        :param dCurrent:
        :param bCurrent:
        :param MpiConfig:
        :return:
        """
        modelPath = self.modelDir + "/" + dCurrent.strftime('%Y%m%d%H00') + ".CHRTOUT_DOMAIN1"

        # if MpiConfig.rank == 0:
        #    print(modelPath)

        # If the file is not present, this may not indicate an issue, but that we are only producing
        # GWOUT files at an infrequent time period. Simply return to the main calling program and leave
        # values in our local arrays to missing.
        if not os.path.isfile(modelPath):
            return

        # Open the NetCDF file.
        idCh = Dataset(modelPath, 'r')

        # Read in the streamflow, along with the feature_id variable.
        # Using the pre-calculated uplinks, we will subset the inflow/outflow volumes,
        # then sum up total volumes for the entire basin.
        varTmp = idCh.variables['feature_id'][:]
        varInd = np.where(varTmp == self.linksGlobal[bCurrent])
        varTmp = None

        # Read in GW inflow and aggregate to the basin.
        varTmp = idCh.variables['streamflow'][varInd]
        self.streamVolLocal[stepCurrent] = varTmp * 3600.0  # Volume of cubic meters.
        varTmp = None

        # Close the NetCDF file.
        idCh.close()

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
            idOut.createVariable('PRCP_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('ECAN_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('EDIR_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('ETRAN_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('CANICE_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('CANLIQ_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('SFCRNOFF_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('UGRDRNOFF_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('SOIL_M_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('SFCHEAD_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('QBDRYRT_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('QSTRMVOLRT_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('GWIN_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('GWOUT_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('Stream_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)
            idOut.createVariable('Accumulated_Streamflow_Volume', np.float64, ('numBasins', 'numSteps'), fill_value=-9999.0)

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

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.accPrcpLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['PRCP_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.accEcanLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['ECAN_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.accEdirLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['EDIR_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.accEtranLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['ETRAN_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.canIceLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['CANICE_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.canLiqLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['CANLIQ_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.sfcRnoffLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['SFCRNOFF_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.uGrdRnoffLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['UGRDRNOFF_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.soilMLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['SOIL_M_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.sfcHeadSubRtLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['SFCHEAD_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.qbdryRtLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['QBDRYRT_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.qStrmVolRtLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['QSTRMVOLRT_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.gwInLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['GWIN_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.gwOutLocal, root=0)

        MpiConfig.comm.barrier()

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['GWOUT_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

        MpiConfig.comm.barrier()

        final = MpiConfig.comm.gather(self.streamVolLocal, root=0)

        if MpiConfig.rank == 0:
            dataOutTmp = np.concatenate([final[i] for i in range(MpiConfig.size)], axis=0)

            # Loop through each basin and place final output variables.
            for bTmp in range(len(self.basinsGlobal)):
                bIndTmp = bTmp * self.nGlobalSteps
                eIndTmp = (bTmp + 1) * self.nGlobalSteps

                idOut.variables['Stream_Volume'][bTmp, :] = dataOutTmp[bIndTmp:eIndTmp]

                sumTmp = dataOutTmp[bIndTmp:eIndTmp]
                sumTmp = np.cumsum(sumTmp)
                idOut.variables['Accumulated_Streamflow_Volume'][bTmp, :] = sumTmp
                sumTmp = None


        MpiConfig.comm.barrier()


        # Close the output netCDF file.
        if MpiConfig.rank == 0:
            idOut.close()
