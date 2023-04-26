from pySmartDL import SmartDL
import time
URL = "https://zinc20.docking.org/substances.smi?count=all"

obj = SmartDL(URL, "./ZINC20_ALL_test.smi",threads=2000) #by default thread = 5
obj.start()
