import os
import logging


class UtilParamCheck:
    # Creates the logger and the required util parameters map
    def __init__(self):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")

        # Map with all required parameters per VaSe Util
        self.required_util_params = {'acceptorcheck': ['varcon', 'vasefq1', 'vasefq2'],
                                     'acceptorreadinfo': ['varcon', 'acceptorbam'],
                                     'checkdonorfiles': ['donorfiles'],
                                     'checkfastq': ['varcon', 'templatefq1', 'templatefq2', 'vasefq1', 'vasefq2'],
                                     'compareacceptor': ['infile1', 'infile2'],
                                     'comparedonor': ['infile1', 'infile2'],
                                     'comparefastq': ['vasefq1', 'vasefq2'],
                                     'comparevarcon': ['varcon', 'varcon2'],
                                     'donorcheck': ['varcon', 'vasefq1', 'vasefq2'],
                                     'donorreadinfo': ['donorfiles', 'varcon'],
                                     'loginfo': ['vaselog'],
                                     'subsetvarcon': ['varcon', 'outfile'],
                                     'subvcfvarcon': ['varconin', 'vcflist'],
                                     'subvcfvarlist': ['variantlist', 'vcflist'],
                                     'unmappedinfo': ['acceptorbam', 'donorfiles', 'unmappedmates', 'unmappedmates2'],
                                     'varcondata': ['donorfiles', 'varcon'],
                                     'vcfinvarcon': ['varcon', 'infile1']
                                     }

        # Map with all optional parameters per VaSe Util
        self.optional_util_params = {'acceptorreadinfo': ['samplefilter', 'varconfilter', 'readidfilter'],
                                     'checkdonorfiles': ['samplefilter'],
                                     'compareacceptor': ['samplefilter', 'varconfilter', 'chromfilter'],
                                     'comparedonor': ['samplefilter', 'varconfilter', 'chromfilter'],
                                     'donorreadinfo': ['samplefilter', 'varconfilter', 'readidfilter'],
                                     'loginfo': ['logfilter'],
                                     'subsetvarcon': ['samplefilter', 'varconfilter', 'chromfilter'],
                                     'unmappedinfo': ['samplefilter', 'varconfilter', 'readidfilter'],
                                     'varcondata': ['samplefilter', 'varconfilter', 'chromfilter'],
                                     'vcfinvarcon': ['samplefilter', 'varconfilter', 'chromfilter']
                                     }

    def required_params_set(self, utiltorun, paramlist):
        """Checks and returns whether the required parameters for the specified VaSeUtil are set.

        Parameters
        ----------
        utiltorun : str
            Namee of the VaSeUtil to run
        paramlist : dict
            Set parameter names and values

        Returns
        -------
        bool
            True if all required parameters are set, False if not
        """
        if utiltorun in self.required_util_params:
            for reqparam in self.required_util_params[utiltorun]:
                if paramlist[reqparam] is not None:
                    if reqparam != "outfile":
                        if not os.path.isfile(paramlist[reqparam]):
                            return False
                else:
                    return False
            return True
        return False

    def get_required_parameters(self, vaseutil):
        """Returns the required parameters for the specified VaseUtil.

        If the specified util name is invalid None will be returned.

        Parameters
        ----------
        vaseutil : str
            VaSeUtil name

        Returns
        -------
        list of str or None
            Required parameter names if VaSeUil name if valid, None otherwise
        """
        if vaseutil in self.required_util_params:
            return self.required_util_params[vaseutil]
        return None

    def get_not_set_parameters(self, utiltorun, paramlist):
        """Returns the names of required parameters that were not set for the specified VaSeUtil.

        Parameters
        ----------
        utiltorun : str
            Name of the specified VaSeUtil
        paramlist : dict

        Returns
        -------
        list of str

        """
        params_not_set = []
        if utiltorun in self.required_util_params:
            for utilparam in self.required_util_params[utiltorun]:
                if paramlist[utilparam] is None:
                    params_not_set.append(utilparam)
            return params_not_set
        return ['util']

    def get_set_optional_parameters(self, utiltorun, paramlist):
        """Returns names of set optional parameters.

        Parameters
        ----------
        utiltorun : str
            Name of VaSeutil to run
        paramlist : dict
            Parameters names and values

        Returns
        -------
        list of str
            Optional parameter names that have been set
        """
        if utiltorun in self.optional_util_params:
            setparams = list(paramlist.keys())
            return list(set(self.optional_util_params[utiltorun]) & set(setparams))

    def get_unused_optional_parameters(self, utiltorun, paramlist):
        """Returns the list of unused optional parameters.

        Parameters
        ----------
        utiltorun : str
            Name of the VaSeUtil to run
        paramlist : dict
            Set parameter names and values

        Returns
        -------
        list of str
            Optional parameter names not set
        """
        if utiltorun in self.optional_util_params:
            setparams = list(paramlist.keys())
            return list(set(self.optional_util_params[utiltorun]) - set(setparams))

    def get_valid_vaseutils(self):
        """Returns the list of valid VaSeUtil names

        Returns
        -------
        list of str
            Valid VaSeUtil names
        """
        return self.required_util_params.keys()

    def optional_parameter_is_set(self, paramname, paramvals):
        """Checks and returns whether a specified parameter is set.

        Parameters
        ----------
        paramname : str
            Name of the parameter to check
        paramvals : dict
            Parameter names and values

        Returns
        -------
        bool

        """
        if paramname in paramvals:
            if paramvals[paramname] is None or paramvals[paramname] == "":
                return False
            return True
        return False
