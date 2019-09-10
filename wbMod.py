class wbObj:
    """
    Object to contain information on basins being processed, date ranges, etc.
    """
    def __init__(self):
        """
        Initialize the object class that will contain information.
        """
        self.basinsGlobal = None
        self.nGlobalBasins = None
        self.nGlobalSteps = None
        self.linksGlobal = None
        self.bDateGlobal = None
        self.eDateGlobal = None
        #self.basinsLocal = None
        #self.bDateLocal = None
        #self.eDateLocal = None