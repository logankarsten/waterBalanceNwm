from netCDF4 import Dataset
import numpy as np
from skimage.transform import downscale_local_mean

class wbObj:
    """
    Object to contain information on basins being processed, date ranges, etc.
    """
    def __init__(self):
        """
        Initialize the object class that will contain information.
        """
        self.geoPath = None
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
                self.bsnMskHydro[bsnTmp][iTmp,jTmp] = self.bsnMskHydro[bsnTmp][iTmp,jTmp] +\
                                                      regridSub[stepTmp]

            # Use the scikit image processing to resample
            self.bsnMskLand[bsnTmp] = downscale_local_mean(self.bsnMskHydro[bsnTmp], (int(self.aggFact), int(self.aggFact)))

            # Reset temporary arrays for next basin tracing.
            doneLinks = None
            polyInd = None
            regridSub = None
            iVarSub = None
            jVarSub = None

            outPathTmp = "BSN_" + str(self.gageIDs[bsnTmp]) + "_RANK_" + str(MpiConfig.rank) + ".nc"
            idTmp = Dataset(outPathTmp, 'w')
            idTmp.createDimension('xHydro', self.nxHydro)
            idTmp.createDimension('yHydro', self.nyHydro)
            idTmp.createDimension('xLand', self.nxLand)
            idTmp.createDimension('yLand', self.nyLand)
            idTmp.createVariable('mask_hydro', np.float32, ('yHydro', 'xHydro'))
            idTmp.variables['mask_hydro'][:, :] = self.bsnMskHydro[bsnTmp]
            idTmp.createVariable('mask_land', np.float32, ('yLand', 'xLand'))
            idTmp.variables['mask_land'][:, :] = self.bsnMskLand[bsnTmp]
            idTmp.close()

        # Close the NetCDF files and reset variables for memory purposes
        idFullDom.close()
        idRt.close()

        toVar = None
        linkVar = None
        idMaskVar = None
        regridWeightVar = None
        iVar = None
        jVar = None

