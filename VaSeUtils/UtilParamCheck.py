import os
import logging

class UtilParamCheck:
    # Creates the logger and the required util parameters map
    def __init__(self):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")

        # Map with all required parameters per VaSe Util
        self.required_util_params = {'acceptorcheck': ['varcon', 'vasefq1', 'vasefq2'],
                                     'acceptorreadinfo': ['varcon', 'acceptorbam'],
                                     'checkfastq': ['varcon', 'templatefq1', 'templatefq2', 'vasefq1', 'vasefq2'],
                                     'compareacceptor': ['varcon', 'varcon2'], 'comparedonor': ['varcon', 'varcon2'],
                                     'comparefastq': ['vasefq1', 'vasefq2'],
                                     'comparevarcon': ['varcon', 'varcon2'],
                                     'donorcheck': ['varcon', 'vasefq1', 'vasefq2'],
                                     'donorreadinfo': ['donorfiles', 'varcon'],
                                     'loginfo': ['vaselog'],
                                     'unmappedinfo': ['acceptorbam', 'donorfiles', 'unmappedmates', 'unmappedmates2'],
                                     'varcondata': ['donorfiles', 'varcon']
                                     }

        # Map with all optional parameters per VaSe Util
        self.optional_util_params = {'acceptorreadinfo' : ['samplefilter', 'varconfilter', 'readidfilter'],
                                     'donorreadinfo' : ['samplefilter', 'varconfilter', 'readidfilter'],
                                     'loginfo' : ['logfilter'],
                                     'unmappedinfo' : ['samplefilter', 'varconfilter', 'readidfilter'],
                                     'varcondata' : ['samplefilter', 'varconfilter', 'chromfilter']
                                     }

    # Check that all the required parameters for a util are set
    def required_params_set(self, utiltorun, paramlist):
        if utiltorun in self.required_util_params:
            for reqparam in self.required_util_params[utiltorun]:
                if paramlist[reqparam] is not None:
                    if not os.path.isfile(paramlist[reqparam]):
                        return False
                else:
                    return False
            return True
        return False

    # Returns the names of the required parameters for a specified VaSe util
    def get_required_parameters(self, vaseutil):
        if vaseutil in self.required_util_params:
            return self.required_util_params[vaseutil]
        return None

    # Returns the names of the parameters not set correctly.
    def get_not_set_parameters(self, utiltorun, paramlist):
        params_not_set = []
        if utiltorun in self.required_util_params:
            for utilparam in self.required_util_params[utiltorun]:
                if paramlist[utilparam] is None:
                    params_not_set.append(utilparam)
            return params_not_set
        return ['util']

    # Returns which optional parameters for a certain utility are set.
    def get_set_optional_parameters(self, utiltorun, paramlist):
        if utiltorun in self.optional_util_params:
            setparams = list(paramlist.keys())
            return list(set(self.optional_util_params[utiltorun]) & set(setparams))

    # Returns the list of unused optional parameters if the correct util is set
    def get_unused_optional_parameters(self, utiltorun, paramlist):
        if utiltorun in self.optional_util_params:
            setparams = list(paramlist.keys())
            return list(set(self.optional_util_params[utiltorun]) - set(setparams))
