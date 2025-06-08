### get started with rgee


install.packages("rgee")
library(rgee)
ee_install(py_env = "rgee")
ee_check()

ee_clean_user_credentials()
ee_Initialize(user = 'jonas.trepel@bio.au.dk') 
ee_get_earthengine_path()

ee_Initialize(user = 'jonas.trepel@bio.au.dk', drive = TRUE)

ee_Authenticate()

library(rgee)
ee_Initialize()
srtm <- ee$Image("USGS/SRTMGL1_003")


library(reticulate)
py_config() # see the name of your conda (python) environment, in my case "r-reticulate" 
reticulate::py_install('earthengine-api==0.1.370', envname='r-reticulate') 

# Check the installation of "earthengine-api" with 
py_list_packages() 
pyl <- py_list_packages()
pyl[pyl$package == "earthengine-api", ]

# check python version with
py_run_string("import sys; print(sys.version)")
