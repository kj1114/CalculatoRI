import json
import cherrypy
import pandas as pd

from subprocess 			import call
from ErlangenWater 			import realPartRefractiveIndex

def addCorsToResponse():
	cherrypy.response.headers['Access-Control-Allow-Origin'] = '*'
	cherrypy.response.headers['Access-Control-Allow-Methods'] = 'GET, POST, PUT, OPTIONS'
	cherrypy.response.headers['Access-Control-Allow-Headers'] = 'Origin, Accept, Content-Type, X-Requested-With, Content-Length'


class KrstinasServer( object ):
	@cherrypy.expose
	def callCPPFunc( self, lambdaReference=None, concentration=None, Nreference=None, radio1=False, file1=None ):
		
		if radio1 == 'true':
			rad1 = True
		else:
			rad1=False

		if isinstance(file1, unicode):
			krs = None
		else:
			krs = file1.file

		if lambdaReference == '':
			lR = 800
		else:
			lR = lambdaReference
		if concentration =='':
			con = 0
		else:
			con=concentration
		if Nreference =='':
			Nref = 1.4005525054
		else:
			Nref=Nreference

		result = realPartRefractiveIndex(lR, rad1, con, Nref, krs)
		result = json.dumps(result)
		return result

def run ():
	conf = {
		'/': {
			'tools.CORS.on': True
		}
	}

	cherrypy.tools.CORS = cherrypy.Tool('before_handler', addCorsToResponse)
	cherrypy.quickstart( KrstinasServer(), '/', conf )

run()