-----------------------------IND_Inventory-----------------------------

This repository contains all functions and data required to create your
industrial emissions inventory from excel file to netCDF. 

If you want to run IND_Inventory you just need to run the script 
IND2netCDF.py. You can adjust this scrip to create your database. 

This repository also contains:

FOLDERS:
	Input folder: folder with excel files containing the industrial
			emissions in Brazil.


FUNCTIONS:
	gridding.py : function used to create the gridded domain and
			populate this grid with emission data. 

	netCDFcreator.py : creates the netCDF file

	IND2netCDF.py : main function.  


OUTPUTS:

	The IND2netCDF.py will create the Outputs folder in your 
	computer with outputs in netCDF. 


REQUIRED DATA:

	Before you run the IND2netCDF.py script you should install
	the requires python packages in requirements.txt. You can do this by
	pip install -r requirements.txt





