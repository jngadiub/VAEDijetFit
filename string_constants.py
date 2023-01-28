# *********************************************************** #
#                  EVENT SAMPLE NUMBERS                       #
# *********************************************************** #

def get_gen_events(sample=""):

	gen_events_number_dir = {
	
	'qcdSide' : 134366091,
	'qcdSig' : 199435365,
	'qcdSideExt' : 90490645,
	'qcdSigExt' : 134264102,
	'GtoWW15na' : 979422,
	'GtoWW15br' : 978558,
	'GtoWW25na' : 983609,
	'GtoWW25br' : 982588,
	'GtoWW35na' : 982038,
	'GtoWW35br' : 972050,
	'GtoWW45na' : 976979,
	'GtoWW45br' : 969741,
	'AtoHZ15' : 99935,
	'AtoHZ25' : 99977,
	'AtoHZ35' : 99979,
	'AtoHZ45' : 99954,
	}

	if sample=='qcdAll':
		return gen_events_number_dir['qcdSide']+gen_events_number_dir['qcdSig']+gen_events_number_dir['qcdSideExt']+gen_events_number_dir['qcdSigExt']
	else:
		return gen_events_number_dir[sample]
   
