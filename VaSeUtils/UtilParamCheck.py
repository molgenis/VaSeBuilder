import os
import logging

class UtilParamCheck:
	# Creates the logger and the required util parameters map
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.requiredUtilParams = {
			'acceptorcheck' : ['varcon', 'vasefq1', 'vasefq2'],
			'acceptorreadinfo' : ['varcon', 'acceptorbam'],
			'checkfastq' : ['varcon', 'templatefq1', 'templatefq2', 'vasefq1', 'vasefq2'],
			'compareacceptor' : ['varcon', 'varcon2'],
			'comparedonor' : ['varcon', 'varcon2'],
			'comparefastq' : ['vasefq1', 'vasefq2'],
			'comparevarcon' : ['varcon', 'varcon2'],
			'donorcheck' : ['varcon', 'vasefq1', 'vasefq2'],
			'donorreadinfo' : ['donorfiles', 'varcon'],
			'loginfo' : ['vaselog'],
			'unmappedinfo' : ['acceptorbam', 'donorfiles', 'unmappedmates', 'unmappedmates2'],
			'varcondata' : ['donorfiles', 'varcon']
			}	# Map with all required parameters for each VaSe Util
		
		self.optionalUtilParams = {
			'acceptorreadinfo' : ['samplefilter', 'varconfilter', 'readidfilter'],
			'donorreadinfo' : ['samplefilter', 'varconfilter', 'readidfilter'],
			'loginfo' : ['logfilter'],
			'unmappedinfo' : ['samplefilter', 'varconfilter', 'readidfilter'],
			'varcondata' : ['samplefilter', 'varconfilter', 'chromfilter']
			}	# Map with all optional parameters for each VaSe Util
	
	
	# Check that all the required parameters for a util are set
	def requiredParamsSet(self, utilToRun, paramList):
		if(utilToRun in self.requiredUtilParams):
			for reqparam in self.requiredUtilParams[utilToRun]:
				if(paramList[reqparam] is not None):
					if(not os.path.isfile(paramList[reqparam])):
						return False
				else:
					return False
			return True
		return False
	
	
	# Returns the names of the required parameters for a specified VaSe util
	def getRequiredParameters(self, vaseutil):
		if(vaseutil in self.requiredUtilParams):
			return self.requiredUtilParams[vaseutil]
		return None
	
	
	# Returns the names of the parameters not set correctly.
	def getNotSetParameters(self, utilToRun, paramList):
		paramsNotSet = []
		if(utilToRun in self.requiredUtilParams):
			for utilparam in self.requiredUtilParams[utilToRun]:
				if(paramList[utilparam] is None):
					paramsNotSet.append(utilparam)
			return paramsNotSet
		return ['util']
	
	
	# Returns which optional parameters for a certain utility are set.
	def getSetOptionalParameters(self, utilToRun, paramList):
		if(utilToRun in self.optionalParams):
			setparams = list(paramList.keys())
			return list(set(self.optionalParams[utilToRun]) & set(setparams))
	
	
	# Returns the list of unused optional parameters if the correct util is set
	def getUnusedOptionalParameters(self, utilToRun, paramList):
		if(utilToRun in self.optionalUtilParams):
			setparams = list(paramList.keys())
			return list(set(self.optionalUtilParams[utilToRun]) - set(setparams))
