import os
import logging

class UtilParamCheck:
	# Creates the logger and the required util parameters map
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.requiredUtilParams = {
			'acceptorcheck' : ['acceptorreads', 'vasefastq1', 'vasefastq2'],
			'acceptorreadinfo' : ['acceptorreads', 'acceptorbam'],
			'checkfastq' : ['acceptorreads', 'donorreads', 'templatefastq1', 'templatefastq2', 'vasefastq1', 'vasefastq2'],
			'comparefastq' : [],
			'comparevarcon' : ['varcon', 'varcon2'],
			'donorcheck' : ['donorreads', 'vasefastq1', 'vasefastq2'],
			'donorreadinfo' : ['donorfiles', 'donorreads'],
			'loginfo' : ['vaselog', 'logfilter'],
			'unmappedinfo' : ['acceptorbam', 'donorfiles', 'unmappedmates', 'unmappedmates2'],
			'varcondata' : ['donorfiles', 'varcon']
			}
		
		self.optionalUtilParams = {
			'acceptorreadinfo' : ['samplefilter', 'varconfilter', 'readidfilter'],
			'donorreadinfo' : ['samplefilter', 'varconfilter', 'readidfilter'],
			'loginfo' : ['logfilter'],
			'varcondata' : ['samplefilter', 'varconfilter', 'chromfilter']
			}
	
	
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
	#def getSetOptionalParameters(self, utilToRun, paramList):
	
	
	# Returns if unused optional
	def getUnusedOptionalParameters(self, utilToRun, paramList):
		if(utilToRun in optionalUtilParams):
			
